#include "Output.hpp"

#include <SDL.h>

#include <list>
#include <cassert>
#include <exception>
#include <iostream>
#include <algorithm>

//local (to this file) data used by the audio system:
namespace {

	//handy constants:
	constexpr uint32_t const MIX_SAMPLES = 1024; //number of samples to mix per call of mix_audio callback; n.b. SDL requires this to be a power of two

	//The audio device:
	SDL_AudioDeviceID device = 0;

}

//public-facing data:

std::vector< Output::Sample > Output::playing_data;
uint32_t Output::playing_position = 0;
float Output::volume = 1.0f;


//This audio-mixing callback is defined below:
void mix_audio(void *, Uint8 *buffer_, int len);

//------------------------ public-facing --------------------------------

void Output::init() {
	if (SDL_InitSubSystem(SDL_INIT_AUDIO) != 0) {
		std::cerr << "Failed to initialize SDL audio subsytem:\n" << SDL_GetError() << std::endl;
		std::cerr << "  (Will continue without audio.)\n" << std::endl;
		return;
	}

	//Based on the example on https://wiki.libsdl.org/SDL_OpenAudioDevice
	SDL_AudioSpec want, have;
	SDL_zero(want);
	want.freq = SampleRate;
	want.format = AUDIO_F32SYS;
	want.channels = 2;
	want.samples = MIX_SAMPLES;
	want.callback = mix_audio;

	device = SDL_OpenAudioDevice(nullptr, 0, &want, &have, 0);
	if (device == 0) {
		std::cerr << "Failed to open audio device:\n" << SDL_GetError() << std::endl;
		std::cerr << "  (Will continue without audio.)\n" << std::endl;
	} else {
		//start audio playback:
		SDL_PauseAudioDevice(device, 0);
		std::cout << "Audio output initialized." << std::endl;
	}
}


void Output::shutdown() {
	if (device != 0) {
		//stop audio playback:
		SDL_PauseAudioDevice(device, 1);
		SDL_CloseAudioDevice(device);
		device = 0;
	}
}


void Output::lock() {
	if (device) SDL_LockAudioDevice(device);
}

void Output::unlock() {
	if (device) SDL_UnlockAudioDevice(device);
}

//------------------------ internals --------------------------------

//The audio callback -- invoked by SDL when it needs more sound to play:
void mix_audio(void *, Uint8 *buffer_, int len) {
	assert(buffer_); //should always have some audio buffer

	static_assert(sizeof(Output::Sample) == 8, "Output::Sample is packed");
	assert(len == MIX_SAMPLES * sizeof(Output::Sample)); //should always have the expected number of samples
	Output::Sample *buffer = reinterpret_cast< Output::Sample * >(buffer_);

	uint32_t to_copy = std::min< uint32_t >(Output::playing_data.size() - Output::playing_position, MIX_SAMPLES);
	if (to_copy > 0) {
		memcpy(buffer, &Output::playing_data[Output::playing_position], MIX_SAMPLES * sizeof(Output::Sample));
		Output::playing_position += to_copy;
	}
	if (to_copy < MIX_SAMPLES) {
		memset(buffer + to_copy, 0, (MIX_SAMPLES - to_copy) * sizeof(Output::Sample));
	}
}


