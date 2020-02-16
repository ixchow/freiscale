#include "Composition.hpp"

#include <kit/read_chunk.hpp>

#include <cassert>
#include <iostream>
#include <fstream>
#include <unordered_map>

Sound const *Composition::add_sound(Sound const &sound) {
	//NOTE: could accelerate with some sort of hash of sounds if content-compare ends up too slow (as might happen with many samples of exactly the same length and initial content)
	for (auto const &have : sounds) {
		if (have == sound) {
			return &have;
		}
	}
	sounds.emplace_back(sound);
	return &sounds.back();
}

void Composition::render(int32_t begin_sample, int32_t end_sample, std::vector< Sample > *buffer_) {
	assert(begin_sample <= end_sample);
	assert(buffer_);
	auto &buffer = *buffer_;
	buffer.assign(end_sample - begin_sample, Sample(0));

	Time begin_time = Time(begin_sample) / Time(SampleRate);
	Time end_time = Time(end_sample) / Time(SampleRate);

	//DEBUG: clear sources info for all triggers:
	for (auto &t : triggers) {
		t.sources.clear();
	}
	//end DEBUG

	for (auto &t : triggers) {
		if (!t.sound) continue;
		if (t.start.t >= end_time) continue;
		if (t.start.t + t.length() < begin_time) continue;

		//samples are located at N / Time(SampleRate) points. So first sample inside trigger is:
		int32_t t_begin_sample = int32_t(std::ceil(t.start.t * SampleRate));

		if (t_begin_sample >= end_sample) continue;

		//So, relative to start of trigger, first sample is...
		float first_sample_time = t_begin_sample / Time(SampleRate) - t.start.t;

		//So write down where (in source sound) time_ofs, time_ofs + 1 / SampleRate, ... come from:
		std::vector< float > samples;
		float length = t.sound->size() / Time(SampleRate);

		float t0 = 0.0f;
		float pos = 0.0f;
		for (uint32_t s = 0; s < t.steps.size(); ++s) {
			if (t.steps[s].t <= 0.0f) continue;
			float p0 = (s == 0 ? t.start.p : t.steps[s-1].p);
			float p1 = t.steps[s].p;
			float t1 = t0 + t.steps[s].t;

			//ax + b == (p1 - p0) / (t1 - t0) * x + p0
			float a = (p1 - p0) / (t1 - t0);
			float b = p0;
			//for every sample that would be taken from this interval...
			while (first_sample_time + (samples.size() / Time(SampleRate)) < t1) {
				//speed = 2^(ax+b) = 2^b * (2^a)^x = 2^b * e^(log(2)ax)
				//position = \int_0^t speed = 2^b / (log(2)a) * ( e^(log(2)ax) - e^(0) )
				//position = 2^b / (log(2)a) * ( e^(log(2)ax) - 1 )
				//compute delta-position using integral:
				float x = first_sample_time + (samples.size() / Time(SampleRate));
				samples.emplace_back(
					pos + std::exp2( b ) / (std::log( 2.0f ) * a ) * ( std::exp2( a * (x - t0) ) - 1.0f )
				);
				if (samples.back() > length) break;
			}
			pos += std::exp2( b ) / (std::log( 2.0f ) * a ) * ( std::exp2( a * (t1 - t0) ) - 1.0f );
			if (pos > length) break;

			t0 = t1;
		}

		if (samples.empty() || samples.back() < length) { //deal with tail:
			float p = (t.steps.empty() ? t.start.p : t.steps.back().p);
			do {
				float x = first_sample_time + (samples.size() / Time(SampleRate));
				samples.emplace_back(
					pos + std::exp2( p ) * (x - t0)
				);
			} while (samples.back() < length);
		}

		//okay, sample sources have been computed!

		t.sources = samples; //DEBUG: store sources info for visualization

		int32_t mix_begin_sample = std::max(t_begin_sample, begin_sample);
		int32_t mix_end_sample = std::min< int32_t >(t_begin_sample + samples.size(), end_sample);

		for (int32_t s = mix_begin_sample; s < mix_end_sample; ++s) {
			float source = samples[s - t_begin_sample] * SampleRate;
			int32_t isource = int32_t(std::round(source));
			if (isource >= 0 && uint32_t(isource) < t.sound->size()) {
				buffer[s - begin_sample] += (*t.sound)[isource];
			}
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

	std::vector< Sound const * > sounds;
	for (auto &snd : snd0) {
		if (snd.begin > snd.end || snd.end > smp0.size()) {
			throw std::runtime_error("snd0 with out-of-range index");
		}
		sounds.emplace_back(ret.add_sound(Sound::from_samples(smp0.data() + snd.begin, smp0.data() + snd.end)));
	}
	std::cout << "Loaded " << sounds.size() << " sounds." << std::endl;

	for (auto &trg : trg0) {
		//NOTE: triggers *must* have at least one point, thus >= in begin/end compare:
		if (trg.begin >= trg.end || trg.end > tps0.size()) {
			throw std::runtime_error("trg0 with out-of-range index");
		}
		if (trg.snd >= sounds.size()) {
			throw std::runtime_error("trg0 with out-of-range sound");
		}
		Trigger trigger;
		trigger.sound = sounds[trg.snd];
		trigger.start = TimeLog2Hz(tps0[trg.begin].t, tps0[trg.begin].p);
		for (uint32_t i = trg.begin + 1; i < trg.end; ++i) {
			trigger.steps.emplace_back(tps0[i].t, tps0[i].p);
		}
		ret.triggers.emplace_back(trigger);
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
		smp0.insert(smp0.end(), sound.begin(), sound.end());
		snd.end = smp0.size();
		snd0.emplace_back(snd);

		sound_index.emplace(&sound, sound_index.size());
	}
	write_chunk("smp0", smp0, &to);
	write_chunk("snd0", snd0, &to);

	std::vector< FS2_tps0 > tps0;
	std::vector< FS2_trg0 > trg0;
	for (auto const &trigger : triggers) {
		FS2_trg0 trg;
		trg.snd = sound_index.at(trigger.sound);
		trg.begin = tps0.size();
		tps0.emplace_back(FS2_tps0{trigger.start.t, trigger.start.p});
		for (auto const &step : trigger.steps) {
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


