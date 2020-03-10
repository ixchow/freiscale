#include "Composition.hpp"

#include <mutex>
#include <condition_variable>
#include <thread>
#include <iostream>

static std::vector< std::shared_ptr< Composition::RenderBlock > > finished;
static std::vector< std::shared_ptr< Composition::RenderBlock > > pending;
static std::mutex mutex;
static std::condition_variable cv;

static std::list< std::thread > render_threads;
static bool quit_flag = false;

constexpr uint32_t BlockPadding = SpectrumSize;

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
				lock.unlock();

				block->render();

				lock.lock();
				finished.emplace_back(block);
			}
		});
	}

	{ //bring in any finished blocks from render queue / mark ready for use:
		std::unique_lock< std::mutex > lock(mutex);
		for (auto &b : finished) {
			b->dirty = false;
			//TODO: maybe trigger texture upload too
		}
		finished.clear();
	}

	//Sort triggers to blocks.
	std::map< int32_t, std::shared_ptr< RenderBlock > > new_blocks;

	//render a bit past song ends:
	int32_t first_block = int32_t(std::floor( ((std::min(begin, loop_begin) - 4.0f) * SampleRate) / float(BlockSize) ));
	int32_t last_block  = int32_t(std::floor( ((std::max(end, loop_end)     + 4.0f) * SampleRate) / float(BlockSize) ));

	for (auto const &trigger : triggers) {
		int32_t min_block = int32_t(std::floor((trigger->begin_sample() + int32_t(BlockPadding)) / float(BlockSize)));
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

	//trim to the current render range:
	while (!rendered.empty() && rendered.begin()->first < first_block) {
		rendered.erase(rendered.begin());
	}
	while (!rendered.empty() && rendered.rbegin()->first > last_block) {
		rendered.erase(--rendered.end());
	}

	auto new_block = new_blocks.begin();
	auto old_block = rendered.begin();
	while (old_block != rendered.end() && new_block != new_blocks.end()) {
		if (old_block->first < new_block->first) {
			//old_block doesn't appear in new list? evict!
			old_block = rendered.erase(old_block);
			continue;
		} else if (new_block->first < old_block->first) {
			//new_block doesn't appear in old list? add!
			rendered.emplace_hint(old_block, *new_block);
			++new_block;
		} else {
			assert(new_block->first == old_block->first);
			if (new_block->second->triggers != old_block->second->triggers) {
				old_block->second = new_block->second;
			}
			++old_block;
			++new_block;
		}
	}

	while (old_block != rendered.end()) {
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

	//TODO: spectrum analysis!
	spectrums.assign( BlockSize / SpectrumStep, std::array< float, SpectrumSize >() );

	//copy to samples:
	samples.assign(padded.begin() + BlockPadding, padded.end() - BlockPadding);
}
