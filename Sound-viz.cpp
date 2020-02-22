#include "Composition.hpp"

#include <cassert>
#include <iostream>
#include <complex>
#include <algorithm>


#include "otfft/otfft.h"
constexpr uint32_t FFTSize = (1<<13);

void Sound::compute_viz() {
	//because writing '*this' is awkward:
	std::vector< Sample > const &samples = *this;

	uint32_t peaks_count = (samples.size() + PeaksStep-1) / PeaksStep;
	peaks.assign(PeaksSlots * peaks_count, std::make_pair(0.0f, 0.0f));

	//(approximately based on example in OTFFT docs)

	//processing buffers:
	OTFFT::simd_array< double > x(FFTSize);
	OTFFT::simd_array< OTFFT::complex_t > y(FFTSize);

	//std::cout << "X alignment % 64 " << ((reinterpret_cast< uint8_t * >(x.p) - reinterpret_cast< uint8_t * >(NULL)) % 64) << std::endl;
	//std::cout << "Y alignment % 64 " << ((reinterpret_cast< uint8_t * >(y.p) - reinterpret_cast< uint8_t * >(NULL)) % 64) << std::endl;

	//FFT processor:
	OTFFT::RFFT rfft(FFTSize);


	//window:
	std::vector< float > window;

	//...gaussian window:
	window.resize(FFTSize);
	{ //smth like this I guess:
		float sigma = 0.1f * SampleRate; //2.0f * 1.0f / float(PeaksRate) * float(SampleRate);
		for (uint32_t i = 0; i < window.size(); ++i) {
			window[i] = 1.0f / (std::sqrt(2.0f * M_PI) * sigma) * std::exp(-float(i)*float(i) / (2.0f * sigma*sigma));
		}
	}

	for (uint32_t peaks_index = 0; peaks_index < peaks_count; ++peaks_index) {
		uint32_t offset = peaks_index * PeaksStep;

		//copy to input buffer:
		assert(offset < samples.size());
		uint32_t len = std::min< uint32_t >(samples.size() - offset, FFTSize);
		for (uint32_t s = 0; s < len; ++s) {
			x[s] = window[s] * samples[offset+s];
		}
		for (uint32_t s = len; s < FFTSize; ++s) {
			x[s] = 0.0;
		}

		//transform:
		rfft.fwd(x.p, y.p);

		//convert to power:
		for (uint32_t s = 0; s < FFTSize; ++s) {
			x[s] = (y[s].Re * y[s].Re + y[s].Im * y[s].Im);
		}

		//Pull out all local peaks:
		std::vector< std::pair< float, float > > found; //(power, freq)

		for (uint32_t s = 2; s + 2 <= FFTSize/2; ++s) {
			// want to see:
			// s0 < s1 > s2
			double sL = x[s-2];
			double s0 = x[s-1];
			double s1 = x[s];
			double s2 = x[s+1];
			double sH = x[s+2];
			if (!(s1 > s0 && s1 > s2)) continue;
			if (!(s1 * 0.1 > std::min(sL, sH))) continue;

			//Basic estimation:
			//float freq = SampleRate / float(FFTSize) * s;
			//float power = s1;

			//quadratic interpolation peaks, as per:
			// https://ccrma.stanford.edu/~jos/parshl/Peak_Detection_Steps_3.html
			float alpha = std::log10(s0);
			float beta = std::log10(s1);
			float gamma = std::log10(s2);

			float p = 0.5f * (alpha-gamma) / (alpha - 2*beta + gamma);

			float freq = SampleRate / float(FFTSize) * (s + p);

			float h = beta - 0.25f * (alpha - gamma) * p;

			float power = std::pow(10.0f, h);

			found.emplace_back(power, freq);
		}

		std::sort(found.begin(), found.end(), [](std::pair< float, float > const &a, std::pair< float, float > const &b) -> bool {
			return a.first > b.first;
		});

		std::pair< float, float > *out = &(peaks[peaks_index * PeaksSlots]);
		std::cout << peaks_index << ": ";
		for (uint32_t i = 0; i < found.size() && i < PeaksSlots; ++i) {
			out[i] = std::make_pair(found[i].second, found[i].first);
			std::cout << out[i].second << "@" << out[i].first << "Hz ";
		}
		std::cout << std::endl;
	}

	float min = std::numeric_limits< float >::infinity();
	float max = -std::numeric_limits< float >::infinity();
	for (auto &[ freq, power ] : peaks) {
		min = std::min(min, power);
		max = std::max(max, power);
	}
	std::cout << "Peak powers in [" << min << ", " << max << "]." << std::endl;

	for (auto &[ freq, power ] : peaks) {
		power /= max;
	}

}
