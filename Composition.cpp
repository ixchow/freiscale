#include "Composition.hpp"

#include <cassert>

void Composition::render(int32_t begin_sample, int32_t end_sample, std::vector< Sample > *buffer_) {
	assert(begin_sample <= end_sample);
	assert(buffer_);
	auto &buffer = *buffer_;
	buffer.assign(end_sample - begin_sample, Sample(0));

	Time begin_time = Time(begin_sample) / Time(SampleRate);
	Time end_time = Time(end_sample) / Time(SampleRate);

	for (auto &t : triggers) {
		if (t.start.t >= end_time) continue;
		if (t.start.t + t.length() < begin_time) continue;

		//TODO: actual sample-warping for rendering
	}
}
