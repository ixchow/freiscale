#include "FreiScale.hpp"
#include "Output.hpp"

#include <kit/kit.hpp>
#include <kit/path.hpp>
#include <kit/Load.hpp>
#include <kit/TagValueArgs.hpp>
#include <kit/load_save_png.hpp>

#include <iostream>

kit::Config kit_config() {
	kit::Config config;
	config.size = glm::uvec2(1280, 720);
	config.title = "FreiScale SemiComposer 2";

	return config;
}

std::shared_ptr< kit::Mode > kit_mode() {
	Output::init();

	kit::call_load_functions(); //only needed if using kit::Load< > for resource loading

	std::string library_path = kit::data_path("sounds");
	std::string load_file = "";
	bool do_fss = false;

	TagValueArgs args;
	args.emplace_back(TagValueArg::simple("library", &library_path, "path to sample library"));
	args.emplace_back(TagValueArg::simple("load", &load_file, "load composition from file"));
	args.emplace_back(TagValueArg::simple("volume", &Output::volume, "set volume"));
	args.emplace_back(TagValueArg::simple("do-fss", &do_fss, "make spectrums for samples then quit"));

	std::string errs;
	if (!args.parse(kit::args.begin() + 1, kit::args.end(), &errs)) {
		std::cerr << "ERROR parsing args:\n" << errs << std::endl;
		std::cerr << args.usage("fs2") << std::endl;
		return nullptr;
	}

	auto fs = std::make_shared< FreiScale >(library_path);

	if (do_fss) {
		//std::cout << "I guess we don't do this any more :-(" << std::endl;
		//build spectrums:
		fs->library.foreach_sound([](std::string const &path, Sound &sound) {
			std::cout << path << std::endl;
			sound.compute_viz();
			/*
			auto dot = path.rfind('.');
			std::string out_path = path.substr(0,dot) + ".fss";
			std::cout << path << " -> " << out_path << std::endl;
			sound.compute_spectrums();
			std::vector< glm::u8vec4 > data;
			data.reserve(sound.spectrums.size());
			for (auto const &s : sound.spectrums) {
				float r = s * 1.0f;
				float g = s * 10.0f;
				float b = s * 100.0f;
				data.emplace_back(
					std::max< int32_t >(0, std::min< int32_t >(255, 255.0f * r)),
					std::max< int32_t >(0, std::min< int32_t >(255, 255.0f * g)),
					std::max< int32_t >(0, std::min< int32_t >(255, 255.0f * b)),
					255
				);
			}
			assert(data.size() == sound.spectrums.size());
			save_png(out_path + ".png", SpectrumBins, data.size() / SpectrumBins, reinterpret_cast< uint32_t const * >(data.data()), UpperLeftOrigin);
			*/

		});
		return nullptr;
	}

	if (load_file != "") {
		std::cerr << "Loading from '" << load_file << "'" << std::endl;
		try {
			*fs->composition = Composition::load(load_file);
		} catch (std::exception &e) {
			std::cerr << "ERROR loading from '" << load_file << "':" << e.what() << std::endl;
			return nullptr;
		}
	}

	return fs;
}

