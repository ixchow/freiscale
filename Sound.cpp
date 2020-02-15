#include "Composition.hpp"

#include <SDL.h>

#include <stdexcept>
#include <iostream>
#include <cassert>

void load_wav(std::string const &filename, std::vector< float > *data_) {
	assert(data_);
	auto &data = *data_;

	SDL_AudioSpec audio_spec;
	Uint8 *audio_buf = nullptr;
	Uint32 audio_len = 0;

	SDL_AudioSpec *have = SDL_LoadWAV(filename.c_str(), &audio_spec, &audio_buf, &audio_len);
	if (!have) {
		throw std::runtime_error("Failed to load WAV file '" + filename + "'; SDL says \"" + std::string(SDL_GetError()) + "\"");
	}

	//based on the SDL_AudioCVT example in the docs: https://wiki.libsdl.org/SDL_AudioCVT
	SDL_AudioCVT cvt;
	SDL_BuildAudioCVT(&cvt, have->format, have->channels, have->freq, AUDIO_F32SYS, 1, SampleRate);
	if (cvt.needed) {
		//std::cout << "'" + filename + "' -> " + std::to_string(SampleRate) + " Hz, float32, mono." << std::endl;
		cvt.len = audio_len;
		cvt.buf = (Uint8 *)SDL_malloc(cvt.len * cvt.len_mult);
		SDL_memcpy(cvt.buf, audio_buf, audio_len);
		SDL_ConvertAudio(&cvt);
		int final_size = cvt.len_cvt;
		assert(final_size >= 0 && final_size <= cvt.len * cvt.len_mult && "Converted audio should fit in buffer.");
		assert(final_size % 4 == 0 && "Converted audio should consist of 4-byte elements.");
		data.assign(reinterpret_cast< float * >(cvt.buf), reinterpret_cast< float * >(cvt.buf + final_size));
		SDL_free(cvt.buf);
	} else {
		data.assign(reinterpret_cast< float * >(audio_buf), reinterpret_cast< float * >(audio_buf + audio_len));
	}
	SDL_FreeWAV(audio_buf);
}

Sound Sound::load(std::string const &path) {
	std::string failures = "";
	//attempt to load as "WAV":
	try {
		std::vector< Sample > data;
		load_wav(path, &data);

		return Sound::from_samples(data.data(), data.data() + data.size());
	} catch (std::runtime_error &e) {
		if (failures != "") failures += "\n";
		failures += "Failed to load as WAV: ";
		failures += e.what();
	}

	throw std::runtime_error("Failed to load '" + path + "':\n" + failures);
}

Sound Sound::from_samples(Sample const *begin, Sample const *end) {
	Sound sound;
	sound.assign(begin, end);
	return sound;
}
