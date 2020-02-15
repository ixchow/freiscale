#pragma once

#include "Composition.hpp"

#include <map>
#include <string>

struct Folder {
	std::map< std::string, Folder > folders;
	std::map< std::string, Sound > sounds;
	enum {
		Collapsed,
		Expanded,
		Error
	} state = Collapsed;
};

struct Library {
	Library(std::string const &path);
	std::string path;
	Folder root;
};
