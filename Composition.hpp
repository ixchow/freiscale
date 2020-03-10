#pragma once

#include <vector>
#include <list>
#include <cstdint>
#include <cmath>
#include <string>
#include <memory>
#include <map>
#include <cassert>

//sample rate:
constexpr uint32_t SampleRate = 48000;

constexpr uint32_t SpectrumRate = 100; //spectrums per second
static_assert(SampleRate % SpectrumRate == 0, "Spectrums start on sample boundaries.");
constexpr uint32_t SpectrumStep = SampleRate / SpectrumRate; //sample offset between subsequent spectrums
constexpr uint32_t SpectrumSize = 2048;

constexpr uint32_t PeaksRate = 100; //peaks per second
constexpr uint32_t PeaksSlots = 20; //peaks to retain per sample
static_assert(SampleRate % PeaksRate == 0, "Peaks frames should start on sample boundaries.");
constexpr uint32_t PeaksStep = SampleRate / PeaksRate; //sample offset between subsequent peaks computations


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

	//visualization/UI stuff:
	//sound handles are shown at this frequency:
	float fundamental = 440.0f;
	//sound is drawn using points at peaks rectangular array of PeaksSlots x (length + PeaksStep-1) / PeaksStep values
	std::vector< std::pair< float, float > > peaks; //(frequency, power)

	void compute_viz();

	static Sound load(std::string const &path); //throws on error
	static Sound from_samples(Sample const *begin, Sample const *end);
};

//Trigger a Sound at a given time/pitch:
struct Trigger {
	Trigger(std::shared_ptr< Sound const > sound_) : sound(sound_) {
		assert(sound);
	}
	std::shared_ptr< Sound const > sound;
	//pitch/time position points:
	std::vector< TimeLog2Hz > steps = { TimeLog2Hz(0.0f, 0.0f), TimeLog2Hz(1.0f, 0.0f) };
	//NOTE: steps are (sorted) *absolute* time, *absolute* speed change
	//*must* have at least two steps.

	bool fix_steps() {
		bool fixed = false;
		if (steps.empty()) {
			steps = { TimeLog2Hz(0.0f, 0.f), TimeLog2Hz(1.0f, 0.0f) };
			fixed = true;
		} else if (steps.size() == 1) {
			steps.emplace_back(steps[0].t + 1.0f, steps[0].p);
			fixed = true;
		} else {
			for (uint32_t i = 1; i < steps.size(); ++i) {
				if (steps[i].t < steps[i-1].t) {
					steps[i].t =  steps[i-1].t;
					fixed = true;
				}
			}
		}
		return fixed;
	}

	const Time begin() { return steps[0].t; }
	const int32_t begin_sample() { return int32_t(std::round(begin() * SampleRate)); }
	const Time end() { return steps.back().t; }
	const int32_t end_sample() { return int32_t(std::round(end() * SampleRate)); }

	//for every output sample in the range [start.t, end.t], use this input sample:
	//(computed from steps)
	std::vector< Time > sources;
	bool sources_dirty = true;
	void compute_sources();
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
	std::vector< std::shared_ptr< Sound > > sounds;
	std::vector< std::shared_ptr< Trigger > > triggers;
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
	std::shared_ptr< Sound const > add_sound(Sound const &sound);

	//Computationing:
	//notice that samples can be negative(!)
	//Required:
	//  begin_sample <= end_sample
	//  buffer != nullptr
	//Renders [begin_sample, end_sample) to buffer (copies from render blocks, actually)
	void render(int32_t begin_sample, int32_t end_sample, std::vector< Sample > *buffer);

	//load/save:
	static Composition load(std::string const &path);
	void save(std::string const &path) const;

	//rendered composition:
	struct RenderBlock {
		int32_t start_sample = 0;
		std::vector< std::shared_ptr< Trigger > > triggers; //triggers assigned to this block for rendering (used to check/set dirty)

		//if true, don't use samples/spectrum (might be in use by render thread):
		bool dirty = true;

		//empty unless block is freshly rendered:
		std::vector< Sample > samples;
		std::vector< std::array< float, SpectrumSize > > spectrums;
		void render(); //compute samples + spectrums

	};
	static constexpr int32_t BlockSize = SampleRate / 2;
	std::map< int32_t, std::shared_ptr< RenderBlock > > rendered;

	//update rendering priorities based on focus:
	void update_rendered(Time focus);
};
