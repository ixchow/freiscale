#include "Composition.hpp"

#include <kit/read_chunk.hpp>

#include <cassert>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <algorithm>

std::shared_ptr< Sound const > Composition::add_sound(Sound const &sound) {
	//NOTE: could accelerate with some sort of hash of sounds if content-compare ends up too slow (as might happen with many samples of exactly the same length and initial content)
	for (auto const &have : sounds) {
		if (static_cast< std::vector< Sample > const & >(*have) == static_cast< std::vector< Sample > const & >(sound)) {
			return have;
		}
	}
	sounds.emplace_back(std::make_shared< Sound >(sound));
	if (sounds.back()->peaks.empty()) {
		sounds.back()->compute_viz();
	}
	return sounds.back();
}

void Composition::render(int32_t begin_sample, int32_t end_sample, std::vector< Sample > *buffer_) {
	assert(begin_sample <= end_sample);
	assert(buffer_);
	auto &buffer = *buffer_;
	buffer.assign(end_sample - begin_sample, Sample(0));

	for (auto const &[ idx, block ] : rendered) {
		if (block->dirty) continue;
		assert(block->samples.size() == BlockSize);

		int32_t begin = std::max(begin_sample, block->start_sample);
		int32_t end = std::min< int32_t >(end_sample, block->start_sample + BlockSize);
		for (int32_t s = begin; s < end; ++s) {
			buffer[s-begin_sample] = block->samples[s-block->start_sample];
		}
	}
}

struct FS2_snd0 {
	uint32_t begin;
	uint32_t end;
};

struct FS2_tps0 {
	float t;
	float p;
};

struct FS2_trg0 {
	uint32_t snd;
	uint32_t begin;
	uint32_t end;
};

struct FS2_ifo0 {
	float begin;
	float end;
	float loop_begin;
	float loop_end;
};

Composition Composition::load(std::string const &path) {
	std::ifstream from(path, std::ios::binary);

	std::vector< Sample > smp0;
	std::vector< FS2_snd0 > snd0;
	std::vector< FS2_tps0 > tps0;
	std::vector< FS2_trg0 > trg0;
	FS2_ifo0 ifo0;

	read_chunk(from, "smp0", &smp0);
	read_chunk(from, "snd0", &snd0);
	read_chunk(from, "tps0", &tps0);
	read_chunk(from, "trg0", &trg0);
	read_struct(from, "ifo0", &ifo0);

	Composition ret;

	std::vector< std::shared_ptr< Sound const > > sounds;
	for (auto &snd : snd0) {
		if (snd.begin > snd.end || snd.end > smp0.size()) {
			throw std::runtime_error("snd0 with out-of-range index");
		}
		sounds.emplace_back(ret.add_sound(Sound::from_samples(smp0.data() + snd.begin, smp0.data() + snd.end)));
	}
	std::cout << "Loaded " << sounds.size() << " sounds." << std::endl;

	for (auto &trg : trg0) {
		//NOTE: triggers *must* have at least two points:
		if (!(trg.begin + 1 < trg.end && trg.end <= tps0.size())) {
			throw std::runtime_error("trg0 with out-of-range index");
		}
		if (trg.snd >= sounds.size()) {
			throw std::runtime_error("trg0 with out-of-range sound");
		}
		std::shared_ptr< Trigger > trigger = std::make_shared< Trigger >(sounds[trg.snd]);
		trigger->steps.clear();
		for (uint32_t i = trg.begin; i < trg.end; ++i) {
			trigger->steps.emplace_back(tps0[i].t, tps0[i].p);
		}
		trigger->fix_steps();
		ret.triggers.emplace_back(trigger);
		trigger->compute_sources();
	}

	std::cout << "Loaded " << ret.triggers.size() << " triggers." << std::endl;

	ret.begin = ifo0.begin;
	ret.end = ifo0.end;
	ret.loop_begin = ifo0.loop_begin;
	ret.loop_end = ifo0.loop_end;

	return ret;
}

void Composition::save(std::string const &path) const {
	std::ofstream to(path, std::ios::binary);

	std::vector< Sample > smp0;
	std::vector< FS2_snd0 > snd0;

	std::unordered_map< Sound const *, uint32_t > sound_index;
	for (auto const &sound : sounds) {
		FS2_snd0 snd;
		snd.begin = smp0.size();
		smp0.insert(smp0.end(), sound->begin(), sound->end());
		snd.end = smp0.size();
		snd0.emplace_back(snd);

		sound_index.emplace(sound.get(), sound_index.size());
	}
	write_chunk("smp0", smp0, &to);
	write_chunk("snd0", snd0, &to);

	std::vector< FS2_tps0 > tps0;
	std::vector< FS2_trg0 > trg0;
	for (auto const &trigger : triggers) {
		FS2_trg0 trg;
		trg.snd = sound_index.at(trigger->sound.get());
		trg.begin = tps0.size();
		for (auto const &step : trigger->steps) {
			tps0.emplace_back(FS2_tps0{step.t, step.p});
		}
		trg.end = tps0.size();
		trg0.emplace_back(trg);
	}

	write_chunk("tps0", tps0, &to);
	write_chunk("trg0", trg0, &to);

	FS2_ifo0 ifo0;
	ifo0.begin = begin;
	ifo0.end = end;
	ifo0.loop_begin = loop_begin;
	ifo0.loop_end = loop_end;

	write_struct("ifo0", ifo0, &to);
}
