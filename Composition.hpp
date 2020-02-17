#pragma once

#include <vector>
#include <list>
#include <cstdint>
#include <cmath>
#include <string>
#include <memory>

//sample rate:
constexpr uint32_t SampleRate = 48000;

constexpr uint32_t SpectrumRate = 200; //spectrums per second
constexpr uint32_t SpectrumStep = SampleRate / SpectrumRate; //sample offset between subsequent spectrums
constexpr uint32_t SpectrumSize = (1 << 11); //spectrum is computed from this many samples

static_assert(SampleRate % SpectrumRate == 0, "Spectrums start on sample boundaries.");

//individual sample:
typedef float Sample;

//positioning for samples:
typedef float Log2Hz;
typedef float Time;
struct TimeLog2Hz {
	TimeLog2Hz(Time const &t_, Log2Hz const &p_) : t(t_), p(p_) { }
	Time t;
	Log2Hz p;
};

//Sounds are lists of Samples @ (by default) SampleRate:
struct Sound : std::vector< Sample > {
	//TODO: autocorrelation info (Log2Hz vs time curves for vis/snapping)

	std::vector< float > spectrums;
	void compute_spectrums();

	static Sound load(std::string const &path); //throws on error
	static Sound from_samples(Sample const *begin, Sample const *end);
};

//Trigger a Sound at a given time/pitch:
struct Trigger {
	Sound const *sound = nullptr;
	TimeLog2Hz start = TimeLog2Hz(0.0f, 0.0f);
	std::vector< TimeLog2Hz > steps; //pitch/time warping
	//NOTE: steps are *delta* time, *absolute* pitch
	//TODO: panning? loudness?

	Time length() const {
		Time remain = (sound ? Time(sound->size()) : Time(0)) / Time(SampleRate);
		Time used = 0.0f;
		for (uint32_t s = 0; s < steps.size(); ++s) {
			if (steps[s].t <= 0.0f) continue;
			//speed = 2^(ax+b) = 2^b * (2^a)^x = 2^b * e^(log(2)ax)
			//position = \int_0^t speed = 2^b / (log(2)a) * ( e^(log(2)ax) - e^(0) )
			//position = 2^b / (log(2)a) * ( e^(log(2)ax) - 1 )
			float p0 = (s == 0 ? start.p : steps[s-1].p);
			float p1 = steps[s].p;
			float dt = steps[s].t;

			float a = (p1 - p0) / dt;
			float b = p0;

			float len = std::exp2( b ) / (std::log(2.0f) * a) * (std::exp2( a * dt ) - 1.0f);

			remain -= len;
			used += dt;

			if (remain <= 0.0f) return used;
		}

		float p = (steps.empty() ? start.p : steps.back().p);
		return used + remain / std::exp2(p);
	}

	//DEBUG:
	std::vector< float > sources; //time points in original sound that map to samples starting at std::ceil(start.t * SampleRate)
};

/*
//Blocks group (and tints?) triggers:
struct Block {
	std::list< Trigger > triggers;
};

//BlockTriggers position (and time-warp?) blocks:
struct BlockTrigger {
	Block const *block = nullptr;
	TimeLog2Hz start = TimeLog2Hz(0.0f, 0.0f);
	std::vector< TimeLog2Hz > steps; //pitch/time warping
};

//Layer includes blocks at different offsets:
struct Layer {
	std::vector< std::pair< Time, Block * > > blocks;
	//TODO: reverb send
	//TODO: high/low shelf curves
};
*/

//Composition holds everything:
struct Composition {
	std::list< Sound > sounds;
	std::list< Trigger > triggers;
	/*
	std::list< Block > blocks; //reference sounds via triggers
	std::list< Layer > layers; //reference blocks via
	*/
	//Markers:
	Time begin = 0.0f;
	Time end = 8.0f;
	Time loop_begin = 0.0f;
	Time loop_end = 8.0f;

	//Helper: add (/dedup) sound:
	Sound const *add_sound(Sound const &sound);

	//Computationing:
	//notice that samples can be negative(!)
	//Required:
	//  begin_sample <= end_sample
	//  buffer != nullptr
	//Renders [begin_sample, end_sample) to buffer
	void render(int32_t begin_sample, int32_t end_sample, std::vector< Sample > *buffer);

	//load/save:
	static Composition load(std::string const &path);
	void save(std::string const &path) const;

};
