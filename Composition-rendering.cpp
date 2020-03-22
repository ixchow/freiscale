#include "Composition.hpp"

#include <mutex>
#include <condition_variable>
#include <thread>
#include <iostream>
#include <algorithm>
#include <array>

#include "otfft/otfft.h"
constexpr uint32_t FFTSize = (1<<13);

static std::vector< std::shared_ptr< Composition::RenderBlock > > finished;
static std::vector< std::shared_ptr< Composition::RenderBlock > > pending;
static std::mutex mutex;
static std::condition_variable cv;

static std::list< std::thread > render_threads;
static bool quit_flag = false;

constexpr uint32_t BlockPadding = FFTSize / 2;

void Composition::update_rendered(Time focus) {

	while (render_threads.size() < 4) {
		render_threads.emplace_back([](){
			std::unique_lock< std::mutex > lock(mutex);
			while (!quit_flag) {
				if (pending.empty()) {
					cv.wait(lock);
					continue;
				}
				std::shared_ptr< RenderBlock > block = pending.back();
				pending.pop_back();
				//std::cout << "Doing block " << block->DEBUG_id << std::endl;
				lock.unlock();

				block->render();

				lock.lock();
				//std::cout << "       done " << block->DEBUG_id << std::endl;
				finished.emplace_back(block);
			}
		});
	}

	{ //bring in any finished blocks from render queue / mark ready for use:
		std::unique_lock< std::mutex > lock(mutex);
		for (auto &b : finished) {
			b->dirty = false;
			assert(b->tex == 0); //only way to get into render list is to be a new block, new blocks have tex == 0
			//do texture upload:
			if (!unused_tex.empty()) {
				b->tex = *unused_tex.begin();
				unused_tex.erase(unused_tex.begin());
			} else {
				glGenTextures(1, &(b->tex));
			}
			glBindTexture(GL_TEXTURE_2D, b->tex);
				glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, SpectrumBins, b->spectrums.size(), 0, GL_RED, GL_FLOAT, b->spectrums.data());
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

				glGenerateMipmap(GL_TEXTURE_2D);

			glBindTexture(GL_TEXTURE_2D, 0);
		}
		finished.clear();
	}

	//Sort triggers to blocks.
	std::map< int32_t, std::shared_ptr< RenderBlock > > new_blocks;

	//render a bit past song ends:
	int32_t first_block = int32_t(std::floor( ((std::min(begin, loop_begin) - 4.0f) * SampleRate) / float(BlockSize) ));
	int32_t last_block  = int32_t(std::floor( ((std::max(end, loop_end)     + 4.0f) * SampleRate) / float(BlockSize) ));

	for (auto const &trigger : triggers) {
		int32_t min_block = int32_t(std::floor((trigger->begin_sample() - int32_t(BlockPadding)) / float(BlockSize)));
		int32_t max_block = int32_t(std::floor((trigger->end_sample() + int32_t(BlockPadding)) / float(BlockSize)));

		min_block = std::max(min_block, first_block);
		max_block = std::min(max_block, last_block);

		for (int32_t b = min_block; b <= max_block; ++b) {
			auto &ptr = new_blocks[b];
			if (!ptr) {
				ptr = std::make_shared< RenderBlock >();
				ptr->start_sample = b * int32_t(BlockSize);
			}
			ptr->triggers.emplace_back(trigger);
		}
	}

	//For each block, check if triggers list matches, if not set dirty and queue new block for render:

	auto free_tex = [this](RenderBlock &rb) {
		if (rb.tex != 0) {
			unused_tex.insert(rb.tex);
			rb.tex = 0;
		}
	};

	//trim to the current render range:

	/* These aren't needed, as far as I can tell?
	while (!rendered.empty() && rendered.begin()->first < first_block) {
		free_tex(*rendered.begin()->second);
		rendered.erase(rendered.begin());
	}
	while (!rendered.empty() && rendered.rbegin()->first > last_block) {
		auto end = rendered.end();
		--end;
		free_tex(*end->second);
		rendered.erase(end);
	}
	*/

	auto new_block = new_blocks.begin();
	auto old_block = rendered.begin();
	while (old_block != rendered.end() && new_block != new_blocks.end()) {
		if (old_block->first < new_block->first) {
			//old_block doesn't appear in new list? evict!
			free_tex(*old_block->second);
			old_block = rendered.erase(old_block);
			continue;
		} else if (new_block->first < old_block->first) {
			//new_block doesn't appear in old list? add!
			rendered.emplace_hint(old_block, *new_block);
			++new_block;
		} else {
			assert(new_block->first == old_block->first);
			if (new_block->second->triggers != old_block->second->triggers) {
				free_tex(*old_block->second);
				old_block->second = new_block->second;
			}
			++old_block;
			++new_block;
		}
	}

	while (old_block != rendered.end()) {
		free_tex(*old_block->second);
		old_block = rendered.erase(old_block);
	}

	while (new_block != new_blocks.end()) {
		rendered.emplace_hint(rendered.end(), *new_block);
		++new_block;
	}

	{ //re-make pending list:
		std::vector< std::shared_ptr< Composition::RenderBlock > > new_pending;
		for (auto &[ block, ptr ] : rendered) {
			if (ptr->dirty) {
				//static uint32_t fresh_id = 1;
				//ptr->DEBUG_id = fresh_id++;
				new_pending.emplace_back(ptr);
			}
		}

		std::unique_lock< std::mutex > lock(mutex);
		pending = std::move(new_pending);

		if (!pending.empty()) {
			cv.notify_all();
		}
	}

}



void Composition::RenderBlock::render() {
	std::vector< Sample > padded(BlockPadding + BlockSize + BlockPadding, 0.0f);

	int32_t padded_begin = start_sample - int32_t(BlockPadding);
	int32_t padded_end = padded_begin + int32_t(padded.size());

	for (auto const &t : triggers) {
		int32_t t_begin = t->begin_sample();
		int32_t t_end = t_begin + int32_t(t->sources.size());

		int32_t begin = std::max(padded_begin, t_begin);
		int32_t end = std::min(padded_end, t_end);
		for (int32_t s = begin; s < end; ++s) {
			int32_t from = int32_t(std::round(t->sources[s-t_begin]));
			if (from > 0 && from < int32_t(t->sound->size())) {
				padded[s-padded_begin] += (*t->sound)[from];
			}
		}
	}

	{ //spectrum analysis!
		spectrums.assign( BlockSize / SpectrumStep, std::array< float, SpectrumBins >() );

		static thread_local OTFFT::simd_array< double > x(FFTSize);
		static thread_local OTFFT::simd_array< OTFFT::complex_t > y(FFTSize);

		//FFT processor:
		static thread_local OTFFT::RFFT rfft(FFTSize);

		const uint32_t spectrum_count = BlockSize / SpectrumStep;
		for (uint32_t spectrum_index = 0; spectrum_index < spectrum_count; ++spectrum_index) {
			uint32_t offset = spectrum_index * SpectrumStep;

			assert(offset + FFTSize <= padded.size());
			//NOTE: padding has space before and after block, so this results in FFTs centered on the desired position
			for (uint32_t s = 0; s < FFTSize; ++s) {
				x[s] = padded[offset + s];
			}

			//transform:
			rfft.fwd(x.p, y.p);

			auto &bins = spectrums[spectrum_index];
			for (auto &v : bins) {
				v = 0.0f;
			}

			//convert to power and place in bins:
			// only do first half, since second half is symmetric negative-frequency stuff:
			// also skip DC
			for (uint32_t s = 1; s < FFTSize / 2; ++s) {
				float power = (y[s].Re * y[s].Re + y[s].Im * y[s].Im);

				power = (std::log10(power) + 10.0f) / 10.0f; //some ad-hoc transform

				float min_freq = (s - 0.5f) * (float(SampleRate) / float(FFTSize));
				float max_freq = (s + 0.5f) * (float(SampleRate) / float(FFTSize));

				float min_bin = (std::log2(min_freq) - std::log2(SpectrumMinHz)) / (std::log2(SpectrumMaxHz) - std::log2(SpectrumMinHz)) * SpectrumBins;
				float max_bin = (std::log2(max_freq) - std::log2(SpectrumMinHz)) / (std::log2(SpectrumMaxHz) - std::log2(SpectrumMinHz)) * SpectrumBins;

				int32_t min_bin_int = int32_t(std::floor(min_bin));
				int32_t max_bin_int = int32_t(std::floor(max_bin));

				if (min_bin_int == max_bin_int) {
					if (min_bin_int >= 0 && min_bin_int < int32_t(bins.size())) {
						bins[min_bin_int] += (max_bin - min_bin) * power;
					}
				} else {
					assert(min_bin_int < max_bin_int);
					if (min_bin_int >= 0 && min_bin_int < int32_t(bins.size())) {
						bins[min_bin_int] += (min_bin_int+1 - min_bin) * power;
					}
					for (int32_t b = min_bin_int + 1; b < max_bin_int; ++b) {
						if (b >= 0 && b < int32_t(bins.size())) {
							bins[b] = power;
						}
					}
					if (max_bin_int >= 0 && max_bin_int < int32_t(bins.size())) {
						bins[max_bin_int] += (max_bin - max_bin_int) * power;
					}
				}

			}

		}

	}

	//copy to samples:
	samples.assign(padded.begin() + BlockPadding, padded.end() - BlockPadding);
}
