#include "Composition.hpp"

void Trigger::compute_sources() {
	sources_dirty = false;
	fix_steps();

	int32_t begin_sample = Trigger::begin_sample();
	int32_t end_sample = Trigger::end_sample();

	//So write down where (in source sound) time_ofs, time_ofs + 1 / SampleRate, ... come from:
	sources.clear();
	sources.reserve(end_sample - begin_sample + 1);

	int32_t s0 = int32_t(std::round(steps[0].t * SampleRate));
	float p0 = steps[0].p;
	sources.emplace_back(0.0f);

	for (uint32_t i = 1; i < steps.size(); ++i) {
		int32_t s1 = int32_t(std::round(steps[i].t * SampleRate));
		float p1 = steps[i].p;
		if (s1 <= s0) {
			p0 = p1;
			continue;
		}

		if (p1 == p0) {
			p1 = p0 + 1e-4f;
		}

		float a = (p1 - p0) / float(s1 - s0);
		float b = p0;

		float base = sources.back();
		for (int32_t s = s0 + 1; s <= s1; ++s) {
			//speed = 2^(ax+b) = 2^b * (2^a)^x = 2^b * e^(log(2)ax)
			//position = \int_0^t speed = 2^b / (log(2)a) * ( e^(log(2)ax) - e^(0) )
			//position = 2^b / (log(2)a) * ( e^(log(2)ax) - 1 )
			//compute delta-position using integral:
			sources.emplace_back(
				base + std::exp2( b ) / (std::log( 2.0f ) * a ) * ( std::exp2( a * (s - s0) ) - 1.0f )
			);
		}

		p0 = p1;
		s0 = s1;
	}

	assert(sources.size() == size_t(end_sample - begin_sample + 1));
}
