#include "FreiScale.hpp"
#include "Output.hpp"

#include <kit/kit.hpp>
#include <kit/path.hpp>
#include <kit/Load.hpp>
#include <kit/TagValueArgs.hpp>

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

	TagValueArgs args;
	args.emplace_back(TagValueArg::simple("library", &library_path, "path to sample library"));

	std::string errs;
	if (!args.parse(kit::args.begin() + 1, kit::args.end(), &errs)) {
		std::cerr << "ERROR parsing args:\n" << errs << std::endl;
		std::cerr << args.usage("fs2") << std::endl;
	}

	return std::make_shared< FreiScale >(library_path);
}

