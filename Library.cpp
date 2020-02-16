
#include "Library.hpp"

#include <algorithm>
#include <iostream>

//will use this helper:
bool list_dir(std::string dirname, std::vector< std::string > *dirnames, std::vector< std::string > *filenames);

void expand_folder(std::string const &path, Folder &folder) {
	folder.state = Folder::Collapsed;
	folder.folders.clear();
	folder.sounds.clear();

	std::vector< std::string > dirnames, filenames;
	if (!list_dir(path, &dirnames, &filenames)) {
		std::cerr << "Failed to list '" << path << "'" << std::endl;
		folder.state = Folder::Error;
		return;
	}
	folder.state = Folder::Expanded;
	for (auto const &dirname : dirnames) {
		folder.folders.emplace(dirname, Folder());
	}
	for (auto const &filename : filenames) {
		try {
			Sound sound = Sound::load(path + "/" + filename);
			folder.sounds.emplace(filename, sound);
		} catch (std::runtime_error &e) {
			std::cerr << "Failed to load '" << (path + "/" + filename) << "' as a sound: " << e.what() << std::endl;
		}
	}
	std::cout << path << " > " << folder.sounds.size() << " sounds." << std::endl; //DEBUG
}

Library::Library(std::string const &path_) : path(path_) {
	expand_folder(path, root);

	for (auto & [ name, folder ] : root.folders) {
		std::cout << path << " / " << name << '\n';
		expand_folder(path + "/" + name, folder);
	}
	std::cout.flush();
}



#ifdef _WIN32
//windows-y directory listing:
#include <io.h>
#include <sys/stat.h>
#include <sys/types.h>


//TODO: test on windows(!)
bool list_dir(std::string dirname, std::vector< std::string > *dirnames, std::vector< std::string > *filenames) {
	struct _finddata_t fileinfo;
	intptr_t handle = _findfirst((dirname + "\\*").c_str(), &fileinfo);
	if (handle == -1) {
		return false;
	}
	do {
		std::string name = std::string(fileinfo.name);
		if (fileinfo.attrib & _A_SUBDIR) {
			if (dirnames) {
				dirnames->emplace_back(name);
			}
		} else {
			if (filenames) {
				filenames->emplace_back(name);
			}
		}
	} while (0 == _findnext(handle, &fileinfo));
	_findclose(handle);

	if (dirnames) {
		std::sort(dirnames->begin(), dirnames->end());
	}
	if (filenames) {
		std::sort(filenames->begin(), filenames->end());
	}
	return true;
}


#else //WINDOWS
//linux-y directory listing:
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>


bool list_dir(std::string dirname, std::vector< std::string > *dirnames, std::vector< std::string > *filenames) {
	DIR *dir = opendir(dirname.c_str());
	if (dir == NULL) {
		return false;
	}
	struct dirent *ent = NULL;
	while ((ent = readdir(dir))) {
		std::string name = std::string(ent->d_name);

		if (ent->d_type == DT_DIR) {
			if (dirnames) {
				if (!(name == "." || name == "..")) {
					dirnames->emplace_back(name);
				}
			}
		} else {
			if (filenames) {
				filenames->emplace_back(name);
			}
		}
	}
	closedir(dir);

	if (dirnames) {
		std::sort(dirnames->begin(), dirnames->end());
	}
	if (filenames) {
		std::sort(filenames->begin(), filenames->end());
	}

	return true;
}


#endif //WINDOWS else

