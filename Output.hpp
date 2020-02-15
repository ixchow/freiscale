#pragma once

#include "Composition.hpp"

#include <glm/glm.hpp>

#include <memory>
#include <vector>
#include <string>
#include <cmath>

//Game audio system. Simplified from f18-base3.
//Uses 48kHz sampling rate.

namespace Output {

// ------- global functions -------

const uint32_t SampleRate = 48000;
static_assert(SampleRate == ::SampleRate, "Sound::SampleRate should match composition.");

struct Sample {
	float l, r;
};

void init(); //call Sound::init() from main.cpp before using any member functions

void shutdown(); //call Sound::shutdown() from main.cpp to gracefully(-ish) exit

//the audio callback doesn't run between Sound::lock() and Sound::unlock()
// the set_*/stop/play/... functions already use these helpers, so you shouldn't need
// to call them unless your code is modifying values directly:
void lock();
void unlock();

//sound just plays this buffer:
extern std::vector< Output::Sample > playing_data;
extern uint32_t playing_position;

//at this volume:
extern float volume;

} //namespace Sound
