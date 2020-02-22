#pragma once

#include "Composition.hpp"

#include <map>
#include <string>
#include <functional>

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
	void foreach_sound(std::function< void(std::string const &path, Sound &sound) > const &fn) {
		std::function< void(std::string, Folder &) > rec;
		rec = [&fn,&rec](std::string const &path, Folder &folder) {
			for (auto &[ file, sound ] : folder.sounds) {
				fn(path + "/" + file, sound);
			}
			for (auto &[ file, subfolder ] : folder.folders) {
				rec(path + "/" + file, subfolder);
			}
		};
		rec(path, root);
	}
};
