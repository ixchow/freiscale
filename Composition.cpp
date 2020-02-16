#include "Composition.hpp"

#include <cassert>
#include <iostream>

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
