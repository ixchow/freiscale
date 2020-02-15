#pragma once

#include "Composition.hpp"
#include "Library.hpp"
#include "UI.hpp"

#include <kit.hpp>

struct Action {
	virtual void handle_event(SDL_Event const &) = 0;
	virtual void draw() = 0;
};

struct FreiScale : kit::Mode {
	FreiScale(std::string const &library_path);
	virtual ~FreiScale();
	virtual void resized() override;
	virtual void update(float elapsed) override;
	virtual void draw() override;
	virtual void handle_event(SDL_Event const &) override;

	Library library;
	std::unique_ptr< Composition > composition;

	//UI state:

	std::unique_ptr< Action > action;

	//NOTE: ui handled in layout pixels
	UI ui;
	UIBox library_box;
	float library_top = 0.0f; //top, in terms of item heights

	//set during drawing:
	std::vector< std::pair< UIBox, Folder const * > > library_folder_boxes;
	std::vector< std::pair< UIBox, Sound const * > > library_sound_boxes;

	UIBox song_box;
	TimeLog2Hz song_center = TimeLog2Hz(0.0f, 0.0f);
	TimeLog2Hz song_radius = TimeLog2Hz(8.0f, 2.5f);

	glm::vec2 mouse;
	struct {
		Folder const *library_folder = nullptr;
		Sound const *library_sound = nullptr;
		void clear() {
			library_folder = nullptr;
			library_sound = nullptr;
		}
	} hovered;
	void update_hovered();
};
