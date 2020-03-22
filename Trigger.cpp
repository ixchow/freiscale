#include "Composition.hpp"

void Trigger::compute_sources() {
	sources_dirty = false;
	fix_steps();

	int32_t begin_sample = Trigger::begin_sample();
	int32_t end_sample = Trigger::end_sample();

	//So write down where (in source sound) time_ofs, time_ofs + 1 / SampleRate, ... come from:
	sources.clear();
	sources.reserve(end_sample - begin_sample + 1);
	source_speeds.clear();
	source_speeds.reserve(end_sample - begin_sample + 1);

	int32_t s0 = int32_t(std::round(steps[0].t * SampleRate));
	float p0 = steps[0].p;
	sources.emplace_back(0.0f);
	source_speeds.emplace_back(p0);

	for (uint32_t i = 1; i < steps.size(); ++i) {
		int32_t s1 = int32_t(std::round(steps[i].t * SampleRate));
		float p1 = steps[i].p;
		if (s1 <= s0) {
			p0 = p1;
			continue;
		}

		/*
		if (std::abs(p1 - p0) < 1e-4f) {
			p1 = p0 + 1e-4f;
		}
		*/

		float a = (p1 - p0) / float(s1 - s0);
		float b = p0;

		float base = sources.back();
		for (int32_t s = s0 + 1; s <= s1; ++s) {
			//speed = 2^(ax+b) = 2^b * (2^a)^x = 2^b * e^(log(2)ax)
			//position = \int_0^t speed = 2^b / (log(2)a) * ( e^(log(2)ax) - e^(0) )
			//position = 2^b / (log(2)a) * ( e^(log(2)ax) - 1 )
			//compute delta-position using integral:
			/*
			sources.emplace_back(
				base + std::exp2( b ) / (std::log( 2.0f ) * a ) * ( std::exp2( a * (s - s0) ) - 1.0f )
			);
			*/
			//version with more stable behavior around a == 0:
			// (thanks to gp/pari for wrangling taylor series expansion of (2^(a*x)-1)/(log(2)*a) )
			float x = (s - s0);
			sources.emplace_back(
				base + std::exp2( b ) * (
					x * (1.0000000000000000000000000000000000000 +
					a*x * (0.34657359027997265470861606072908828404 +
					a*x * (0.080075502319700237444517087721110828622 +
					a*x * (0.013876027166205394988285565942155439340 +
					a*x * (0.0019236258215256954323958143147317730960 +
					a*x * (0.00022222596910714072372353703313326957885 +
					a*x * (2.2005043419116585649195853332489068516e-5 +
					a*x * (1.9065917255074800350031798765012047712e-6 +
					a*x * (1.4683874211271454987115286920320400840e-7 +
					a*x * (1.0178086009239699727490007597744629255e-8 +
					a*x * (6.4135560189101121180685383468643075080e-10 +
					a*x * (3.7046152265590095813303404657400933775e-11 +
					a*x * (1.9752643071914003955380617224485303465e-12 +
					a*x * (9.7796348956458063434942528538107647042e-14 +
					a*x * (4.5191509032150304222994028759897642408e-15 +
					a*x * (1.9577729419302678885218402756081727125e-16
					/* + ... */
					))))))))))))))))
				)
			);

			source_speeds.emplace_back( a * x + b );
		}

		p0 = p1;
		s0 = s1;
	}

	assert(sources.size() == size_t(end_sample - begin_sample + 1));
	assert(source_speeds.size() == size_t(end_sample - begin_sample + 1));
}
