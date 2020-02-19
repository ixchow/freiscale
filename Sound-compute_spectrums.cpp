#include "Composition.hpp"

#include "otfft/otfft.h"

#include <cassert>
#include <iostream>

constexpr uint32_t FFTSize = 2048;

void Sound::compute_spectrums() {

	uint32_t spectrum_count = (this->size() + SpectrumStep-1) / SpectrumStep;
	spectrums.assign(SpectrumBins * spectrum_count, 0.0f);

	//Based on example in OTFFT docs!

	//processing buffers:
	OTFFT::simd_array< double > x(FFTSize);
	OTFFT::simd_array< OTFFT::complex_t > y(FFTSize);

	std::cout << "X alignment % 64 " << ((reinterpret_cast< uint8_t * >(x.p) - reinterpret_cast< uint8_t * >(NULL)) % 64) << std::endl;
	std::cout << "Y alignment % 64 " << ((reinterpret_cast< uint8_t * >(y.p) - reinterpret_cast< uint8_t * >(NULL)) % 64) << std::endl;

	//FFT processor:
	OTFFT::RFFT rfft(FFTSize);

	//window:
	std::vector< float > window;

	//...rectangular window:
	//  wikipedia says a lot about this topic :-/
	window.assign(FFTSize, 1.0f);

	auto do_rate = [&](uint32_t sample_rate, std::vector< Sample > const &samples) {
		//compute what ranges of spectrum will pull from this fft:

		uint32_t max_frequency_index = FFTSize / 2;
		float max_frequency = sample_rate / 2.0f;
		uint32_t min_frequency_index = 0;
		float min_frequency = sample_rate / FFTSize;

		std::vector< float > bin_starts;
		bin_starts.reserve(SpectrumBins+1);
		for (uint32_t bin = 0; bin <= SpectrumBins; ++bin) {
			float bin_start_hz = std::exp2( (std::log2(SpectrumMaxFreq) - std::log2(SpectrumMinFreq)) * bin / float(SpectrumBins) + std::log2(SpectrumMinFreq) );
			if (bin_start_hz >= min_frequency && bin_start_hz <= max_frequency) {
				float bin_start_index = (bin_start_hz - min_frequency) / (max_frequency - min_frequency) * (max_frequency_index - min_frequency_index) + min_frequency_index;
				bin_starts.emplace_back(bin_start_index);
			} else {
				bin_starts.emplace_back(std::numeric_limits< float >::quiet_NaN());
			}
		}
		assert(bin_starts.size() == SpectrumBins+1);

		//DEBUG:
		uint32_t first_real = 0;
		while (first_real < bin_starts.size() && !(bin_starts[first_real] == bin_starts[first_real])) ++first_real;
		uint32_t last_real = bin_starts.size() - 1;
		while (last_real < bin_starts.size() && !(bin_starts[last_real] == bin_starts[last_real])) --last_real;

		std::cout << "At " << sample_rate << " and FFT size " << FFTSize << " bins starting from " << first_real << " to " << last_real << " of " << SpectrumBins << " map to computed range." << std::endl;

		//PARANOIA: it's a solid range of valid bins, right?
		for (uint32_t i = first_real; i <= last_real; ++i) {
			assert(bin_starts[i] == bin_starts[i]);
			//std::cout << i << ": " << bin_starts[i] << "\n";
			if (i < last_real) assert(bin_starts[i] < bin_starts[i+1]);
		}
		std::cout.flush();

		assert(sample_rate % SpectrumRate == 0 && "SpectrumRate must respect decimated samples as well.");
		uint32_t rate_step = sample_rate / SpectrumRate;

		for (uint32_t spectrum_index = 0; spectrum_index < spectrum_count; ++spectrum_index) {
			uint32_t offset = spectrum_index * rate_step;

			//copy to input buffer:
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

			//copy to bins:
			float *bins = &(spectrums[spectrum_index * SpectrumBins]);
			for (uint32_t b = first_real; b + 1 <= last_real; ++b) {
				float begin = bin_starts[b];
				float end = bin_starts[b+1];
				//assert(end > begin);
				uint32_t ibegin = uint32_t(std::floor(begin));
				uint32_t iend = uint32_t(std::floor(end));

				//Hmm, should we add or average? I'm going to go with average for now.
				if (ibegin == iend) {
					bins[b] = x[ibegin];
				} else {
					float fbegin = begin - ibegin;
					float fend = end - iend;
					float total = 0.0f;
					total += (1.0f - fbegin) * x[ibegin];
					for (uint32_t i = ibegin + 1; i < iend; ++i) {
						total += x[i];
					}
					total += fend * x[iend];
					bins[b] = total / (end - begin);
				}
			}
		}

	};

	do_rate(SampleRate, *this);

	auto downsample2 = [](std::vector< Sample > const &src) -> std::vector< Sample > {
		std::vector< Sample > dst;
		dst.reserve((src.size()+1)/2);
		//downsample with 1 -> 0.25 0.5 0.25 kernel
		for (uint32_t d = 0; d < (src.size()+1)/2; ++d) {
			uint32_t s = 2*d;
			float sum = 0.0f;
			if (s - 1U < src.size()) sum += 0.25f * src[s-1U];
			if (s < src.size()) sum += 0.5f * src[s];
			if (s + 1U < src.size()) sum += 0.25f * src[s+1U];
			dst.emplace_back(sum);
		}
		return dst;
	};

	std::vector< Sample > div2 = downsample2(*this);
	do_rate(SampleRate/2, div2);

	std::vector< Sample > div4 = downsample2(div2);
	do_rate(SampleRate/4, div4);

	std::vector< Sample > div8 = downsample2(div4);
	do_rate(SampleRate/8, div8);

	std::vector< Sample > div16 = downsample2(div8);
	do_rate(SampleRate/16, div16);

	float min = std::numeric_limits< float >::infinity();
	float max = -std::numeric_limits< float >::infinity();
	for (auto &s : spectrums) {
		min = std::min(min, s);
		max = std::max(max, s);
	}
	std::cout << "Spectrum values in [" << min << ", " << max << "]." << std::endl;

}
