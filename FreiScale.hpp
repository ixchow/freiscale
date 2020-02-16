#pragma once

#include "Composition.hpp"
#include "Library.hpp"
#include "UI.hpp"
#include "Output.hpp"

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
	UIBox library_box;
	float library_top = 0.0f; //top, in terms of item heights

	//set during drawing:
	std::vector< std::pair< UIBox, Folder const * > > library_folder_boxes;
	std::vector< std::pair< UIBox, Sound const * > > library_sound_boxes;

	UIBox song_box;
	TimeLog2Hz song_center = TimeLog2Hz(8.0f, 0.0f);
	TimeLog2Hz song_radius = TimeLog2Hz(8.5f, 2.5f);

	glm::vec2 mouse; //in layout pixels
	TimeLog2Hz get_song_position(glm::vec2 const &px) const {
		return TimeLog2Hz(
			((px.x - song_box.min.x) / (song_box.max.x - song_box.min.x) * 2.0f - 1.0f) * song_radius.t + song_center.t,
			((px.y - song_box.min.y) / (song_box.max.y - song_box.min.y) * 2.0f - 1.0f) * song_radius.p + song_center.p
		);
	}
	glm::vec2 get_screen_position(TimeLog2Hz const &tp) const {
		return glm::vec2(
			((tp.t - song_center.t) / song_radius.t * 0.5f + 0.5f) * (song_box.max.x - song_box.min.x) + song_box.min.x,
			((tp.p - song_center.p) / song_radius.p * 0.5f + 0.5f) * (song_box.max.y - song_box.min.y) + song_box.min.y
		);
	}


	struct {
		Folder const *library_folder = nullptr;
		Sound const *library_sound = nullptr;
		std::pair< Trigger *, uint32_t > song_trigger_handle = std::make_pair(nullptr, 0U); //NOTE: 0 is the 'start' handle, 1 - N the 'step' handles
		std::pair< Trigger *, uint32_t > song_trigger_segment = std::make_pair(nullptr, 0U); //NOTE: 0 is the segment from the start to the first step
		void clear() {
			library_folder = nullptr;
			library_sound = nullptr;
			song_trigger_handle = std::make_pair(nullptr, 0U);
			song_trigger_segment = std::make_pair(nullptr, 0U);
		}
	} hovered;
	void update_hovered();

	Sound const *current_library_sound = nullptr;

	//DEBUG:

	std::vector< Output::Sample > show_sound;

};
