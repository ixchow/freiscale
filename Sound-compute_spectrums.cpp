#include "Composition.hpp"

#include <cassert>
#include <iostream>
#include <complex>
#include <algorithm>


#ifdef USE_CONVOLVE
void Sound::compute_spectrums() {
	//Try using convolution with a "generic impulse" at each target frequency.
	
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

	//build a nicely smoothed pyramid, yo:
	//n.b. could also use gaussians here. would prolly be better. oh well.
	std::vector< std::vector< Sample > > pyramid;
	pyramid.reserve(20);
	pyramid.emplace_back(*this);
	while (pyramid.size() < 20) {
		pyramid.emplace_back(downsample2(pyramid.back()));
	}

	std::vector< std::vector< Sample > > filtered;
	std::vector< uint32_t > levels;

	filtered.reserve(SpectrumBins);

	//float prefilter_sigma = SampleRate / (2.0f * M_PI * (SampleRate / 2.0f) / std::sqrt(2.0f * std::log( std::sqrt( 2.0f ) ) ));

	for (uint32_t bin = 0; bin < SpectrumBins; ++bin) {
		float min_freq = std::exp2(bin / float(SpectrumBins) * (std::log2(SpectrumMaxFreq) - std::log2(SpectrumMinFreq)) + std::log2(SpectrumMinFreq));
		float max_freq = std::exp2((bin + 1) / float(SpectrumBins) * (std::log2(SpectrumMaxFreq) - std::log2(SpectrumMinFreq)) + std::log2(SpectrumMinFreq));

		float freq = 0.5f * (max_freq + min_freq);

		float wavelength, sigma, radius;
		uint32_t l = -1U;
		do {
			++l;
			wavelength = (float(SampleRate) / float(1 << l)) / freq;
			//sigma set relative to wavelength (sure, why not?)
			sigma = 0.5f * (8 * wavelength);
			//radius in terms of threshold for smallest kernel value:
			radius = 6.0f * wavelength; //std::sqrt( -(2.0f * sigma*sigma) * std::log(1e-6 * (std::sqrt(2.0f * M_PI) * sigma)) );
		} while (wavelength > 16.0f && (l + 1 < pyramid.size()));

		if (bin % 10 == 0) std::cout << "Freqency: " << freq << " -> wavelength: " << wavelength << " -> sigma: " << sigma << " [" << std::ceil(radius) << "] @ " << (l+1) << "/" << pyramid.size() << std::endl;

		std::vector< std::complex< float > > kernel(std::ceil(radius));
		for (uint32_t i = 0; i < kernel.size(); ++i) {
			float g = 1.0f / (std::sqrt(2.0f * M_PI) * sigma) * std::exp(-float(i)*float(i) / (2.0f * sigma*sigma));
			g = std::exp2(-(i / wavelength));
			float a = std::sin(float(2.0f * M_PI) / wavelength) * i;
			kernel[i].real(g * std::sin(a));
			kernel[i].imag(g * std::cos(a));
		}

		std::vector< Sample > const &samples = pyramid[l];

		std::vector< Sample > out;
		out.reserve(samples.size());

		for (int32_t s = 0; s < int32_t(samples.size()); ++s) {
			int32_t b = s;
			int32_t e = b + kernel.size();

			int32_t sb = std::max(0, b);
			int32_t se = std::min(e, int32_t(samples.size()));

			std::complex< double > sum = 0.0;
			for (int32_t i = sb; i < se; ++i) {
				sum += kernel[i - b] * samples[i];
			}
			out.emplace_back(std::norm(sum));
		}

		filtered.emplace_back(out);
		levels.emplace_back(l);
	}

	assert(filtered.size() == SpectrumBins);
	assert(levels.size() == filtered.size());

	uint32_t spectrum_count = (this->size() + SpectrumStep-1) / SpectrumStep;
	spectrums.assign(SpectrumBins * spectrum_count, 0.0f);

	for (uint32_t bin = 0; bin < SpectrumBins; ++bin) {
		std::vector< Sample > const &result = filtered[bin];

		for (uint32_t s = 0; s < result.size(); ++s) {
			//NOTE: does using squared amplitude interact badly with downsampling?
			float amt = result[s];
			amt *= (1 << levels[bin]);
			uint32_t spectrum = (s << levels[bin]) / SpectrumStep;
			spectrums[spectrum * SpectrumBins + bin] += amt; //<<--- if this approach is better, should transpose
		}

	}

	float min = std::numeric_limits< float >::infinity();
	float max = -std::numeric_limits< float >::infinity();
	for (auto &s : spectrums) {
		min = std::min(min, s);
		max = std::max(max, s);
	}
	std::cout << "Spectrum values in [" << min << ", " << max << "]." << std::endl;
	for (auto &s : spectrums) {
		s = (s - min) / (max - min);
	}

}

#endif //USE_CONVOLVE



#ifdef USE_FILTERBANK
void Sound::compute_spectrums() {
	//Try using a sort of "gaussian filter bank" to compute spectrum bins.
	//idea: compute gaussian (lowpass) with half-power at divisions in frequency, subtract to bandpass, locally sum for energy.
	
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

	//build a nicely smoothed pyramid, yo:
	//n.b. could also use gaussians here. would prolly be better. oh well.
	std::vector< std::vector< Sample > > pyramid;
	pyramid.reserve(10);
	pyramid.emplace_back(*this);
	while (pyramid.size() < 10) {
		pyramid.emplace_back(downsample2(pyramid.back()));
	}


	std::vector< std::vector< Sample > > filtered;
	std::vector< uint32_t > levels;

	filtered.reserve(SpectrumBins+1);

	//float prefilter_sigma = SampleRate / (2.0f * M_PI * (SampleRate / 2.0f) / std::sqrt(2.0f * std::log( std::sqrt( 2.0f ) ) ));

	for (uint32_t bin = 0; bin <= SpectrumBins; ++bin) {
		float freq = std::exp2(bin / float(SpectrumBins) * (std::log2(SpectrumMaxFreq) - std::log2(SpectrumMinFreq)) + std::log2(SpectrumMinFreq));

		#define USE_GAUSSIAN
		#ifdef USE_GAUSSIAN
		//Wikipedia suggests this is going to be the half-power point 
		float sigma_f = freq / std::sqrt(2.0f * std::log( std::sqrt( 2.0f ) ) );

		//std::cout << "Freqency: " << freq << " -> sigma: " << sigma << " + " << prefilter_sigma << " = ";
		//sigma = std::sqrt(sigma*sigma + prefilter_sigma*prefilter_sigma);
		//std::cout << sigma << std::endl;

		//TODO: kernel width reduction / downsampling!

		float sigma, radius;
		uint32_t l = -1U;
		do {
			++l;
			sigma = (float(SampleRate) / float(1 << l)) / (2.0f * M_PI * sigma_f);
			//radius in terms of threshold for smallest kernel value:
			radius = std::sqrt( -(2.0f * sigma*sigma) * std::log(1e-6 * (std::sqrt(2.0f * M_PI) * sigma)) );
		} while (radius > 5.0f && (l + 1 < pyramid.size()));

		if (bin % 10 == 0) std::cout << "Freqency: " << freq << " -> sigma: " << sigma << " [" << std::ceil(radius) << "] @ " << (l+1) << "/" << pyramid.size() << std::endl;

		std::vector< float > half_kernel(std::ceil(radius));
		for (uint32_t i = 0; i < half_kernel.size(); ++i) {
			half_kernel[i] = 1.0f / (std::sqrt(2.0f * M_PI) * sigma) * std::exp(-float(i)*float(i) / (2.0f * sigma*sigma));
		}

		#endif

		#ifdef USE_SINC

		uint32_t l = -1U;
		do {
			++l;
		} while (SampleRate / float(1 << l) > 5.0f * freq && l + 1 < pyramid.size());
		float radius = 50.0f; //some ad-hoc radius
		std::vector< float > half_kernel(std::ceil(radius));

		if (bin % 10 == 0) std::cout << "Freqency: " << freq << " [" << std::ceil(radius) << "] @ " << (l+1) << "/" << pyramid.size() << std::endl;

		for (uint32_t i = 0; i < half_kernel.size(); ++i) {
			float t = float(i) / (float(SampleRate) / float(1 << l));
			half_kernel[i] = std::sin(2.0f * M_PI * freq * t) / float(M_PI * t);
		}

		#endif

		std::vector< float > kernel;
		kernel.reserve(2*half_kernel.size()-1);
		kernel.insert(kernel.end(), half_kernel.rbegin(), half_kernel.rend());
		kernel.insert(kernel.end(), half_kernel.begin()+1, half_kernel.end());


		std::vector< Sample > const &samples = pyramid[l];

		std::vector< Sample > out;
		out.reserve(samples.size());

		int32_t offset = -int32_t((kernel.size()-1)/2);

		for (int32_t s = 0; s < int32_t(samples.size()); ++s) {
			int32_t b = s + offset;
			int32_t e = b + kernel.size();

			int32_t sb = std::max(0, b);
			int32_t se = std::min(e, int32_t(samples.size()));

			double sum = 0.0;
			for (int32_t i = sb; i < se; ++i) {
				sum += kernel[i - b] * samples[i];
			}
			out.emplace_back(sum);
		}

		filtered.emplace_back(out);
		levels.emplace_back(l);
	}

	assert(filtered.size() == SpectrumBins+1);
	assert(levels.size() == filtered.size());

	uint32_t spectrum_count = (this->size() + SpectrumStep-1) / SpectrumStep;
	spectrums.assign(SpectrumBins * spectrum_count, 0.0f);

	for (uint32_t bin = 0; bin <= SpectrumBins; ++bin) {
		std::vector< Sample > const &above = filtered[bin+1];
		std::vector< Sample > const &below = filtered[bin];

		if (above.size() == below.size()) {
			for (uint32_t s = 0; s < above.size(); ++s) {
				//NOTE: does using squared amplitude interact badly with downsampling?
				float amt = (above[s] - below[s]);
				amt = amt*amt;
				uint32_t spectrum = (s << levels[bin]) / SpectrumStep;
				spectrums[spectrum * SpectrumBins + bin] += amt * (1 << levels[bin]); //<<--- if this approach is better, should transpose
			}
		} else {
			//TODO: resamplin'
		}

	}

	float min = std::numeric_limits< float >::infinity();
	float max = -std::numeric_limits< float >::infinity();
	for (auto &s : spectrums) {
		min = std::min(min, s);
		max = std::max(max, s);
	}
	std::cout << "Spectrum values in [" << min << ", " << max << "]." << std::endl;

}

#endif //USE_FILTERBANK


#define USE_FFT
#ifdef USE_FFT

#include "otfft/otfft.h"
constexpr uint32_t FFTSize = (1<<13);

void Sound::compute_spectrums() {

	uint32_t spectrum_count = (this->size() + SpectrumStep-1) / SpectrumStep;
	spectrums.assign(SpectrumBins * spectrum_count, 0.0f);

	//Based on example in OTFFT docs!

	//processing buffers:
	OTFFT::simd_array< double > x(FFTSize);
	OTFFT::simd_array< OTFFT::complex_t > y(FFTSize);

	//std::cout << "X alignment % 64 " << ((reinterpret_cast< uint8_t * >(x.p) - reinterpret_cast< uint8_t * >(NULL)) % 64) << std::endl;
	//std::cout << "Y alignment % 64 " << ((reinterpret_cast< uint8_t * >(y.p) - reinterpret_cast< uint8_t * >(NULL)) % 64) << std::endl;

	//FFT processor:
	OTFFT::RFFT rfft(FFTSize);

	//window:
	std::vector< float > window;

	//...rectangular window:
	//  wikipedia says a lot about this topic :-/
	//window.assign(FFTSize, 1.0f);


	auto do_rate = [&](uint32_t sample_rate, std::vector< Sample > const &samples) {
		//compute what ranges of spectrum will pull from this fft:

		//...gaussian window:
		window.resize(FFTSize);
		{ //smth like this I guess:
			float sigma = 0.05f * sample_rate;
			for (uint32_t i = 0; i < window.size(); ++i) {
				window[i] = 1.0f / (std::sqrt(2.0f * M_PI) * sigma) * std::exp(-float(i)*float(i) / (2.0f * sigma*sigma));
			}
		}


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

			#define USE_PEAKS
			#ifdef USE_PEAKS
			/* rejecting peaks seems to filter high frequencies too much
			//reject some peaks:
			double max = 0.0;
			for (uint32_t s = 0; s < FFTSize; ++s) {
				max = std::max(max, x[s]);
			}
			double threshold = 0.1 * max;
			for (uint32_t s = 0; s < FFTSize; ++s) {
				if (x[s] < threshold) x[s] = 0.0;
			}
			*/
			//Pull out all local peaks:
			std::vector< std::pair< float, float > > peaks;

			float max_power = 0.0f;
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
				//float freq = sample_rate / float(FFTSize) * s;
				//float power = s1;

				//quadratic interpolation peaks, as per:
				// https://ccrma.stanford.edu/~jos/parshl/Peak_Detection_Steps_3.html
				float alpha = std::log(s0);
				float beta = std::log(s1);
				float gamma = std::log(s2);

				float p = 0.5f * (alpha-gamma) / (alpha - 2*beta + gamma);

				float freq = sample_rate / float(FFTSize) * (s + p);

				float h = beta - 0.25f * (alpha - gamma) * p;

				float power = std::exp(h);


				if (s1 > s0 && s1 > s2) {
					peaks.emplace_back(power, freq);
					max_power = std::max(max_power, power);
				}
			}

			std::sort(peaks.begin(), peaks.end(), std::greater< std::pair< float, float > >());
			if (peaks.size() > 20) {
				peaks.erase(peaks.begin() + 20, peaks.end());
			}

			float *bins = &(spectrums[spectrum_index * SpectrumBins]);
			for (auto const &[ power , freq ] : peaks) {
				if (freq < SpectrumMinFreq || freq > SpectrumMaxFreq) continue;
				//recall: std::log2(freq) = std::exp2( b * (log(max) - log(min)) + log(min) )
				float b = SpectrumBins * (std::log2(freq) - std::log2(SpectrumMinFreq)) / (std::log2(SpectrumMaxFreq) - std::log2(SpectrumMinFreq));
				int32_t bi = int32_t(std::floor(b));
				bi = std::max< int32_t >(0, std::min< int32_t >(SpectrumBins-1, bi));
				float bf = b - bi;
				bins[bi] += (1.0f - bf) * power;
				bins[bi+1] += bf * power;
			}
			/*
			for (uint32_t b = 0; b < FFTSize && b < SpectrumBins; ++b) {
				bins[b] = x[b];
			}
			*/

			#endif //USE_PEAKS

			#ifdef USE_COPY
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
			#endif //USE_COPY
		}

	};

	do_rate(SampleRate, *this);
	/*

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
	*/
/*
	double avg_power = 0.0;
	for (auto &s : *this) {
		avg_power += s*s;
	}
	avg_power /= this->size();

	std::cout << "avg_power: " << avg_power << " * FFTSize: " << avg_power * FFTSize << std::endl;

	for (auto &s : spectrums) {
		s = std::log10(s) - std::log10(1e-6);
		s = std::max(0.0f, s);
	}
*/

	float min = std::numeric_limits< float >::infinity();
	float max = -std::numeric_limits< float >::infinity();
	for (auto &s : spectrums) {
		min = std::min(min, s);
		max = std::max(max, s);
	}
	std::cout << "Spectrum values in [" << min << ", " << max << "]." << std::endl;

	for (auto &s : spectrums) {
		s = (s - min) / (max - min);
	}

}

#endif //USE_FFT


#ifdef USE_AUTOCORRELATION
//Autocorrelation gets first peak okay, but not so good at getting harmonics (maybe?)

#include "otfft/otfft.h"
constexpr uint32_t FFTSize = 4096;

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

			//TEST: let's just remove a bunch of frequency content, eh?
			for (uint32_t s = FFTSize/2; s < FFTSize; ++s) {
				x[s] = 0.0;
			}

			//transform:
			rfft.fwd(x.p, y.p);

			//convolve with self:
			for (uint32_t s = 0; s < FFTSize; ++s) {
				y[s] = y[s] * y[s];
			}

			//transform back:
			rfft.inv(y.p, x.p);

			//copy to bins:
			float *bins = &(spectrums[spectrum_index * SpectrumBins]);
			for (uint32_t lag = 0; lag < FFTSize/2 && lag < SpectrumBins; ++lag) {
				bins[lag] = std::max< float >(0.0f, x[lag]);
			}

#if 0
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
#endif
		}

	};

	do_rate(SampleRate, *this);

/*
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
*/
	float min = std::numeric_limits< float >::infinity();
	float max = -std::numeric_limits< float >::infinity();
	for (auto &s : spectrums) {
		min = std::min(min, s);
		max = std::max(max, s);
	}
	std::cout << "Spectrum values in [" << min << ", " << max << "]." << std::endl;

	for (auto &s : spectrums) {
		s = (s - min) / (max - min);
	}


}

#endif //USE_AUTOCORRELATION
