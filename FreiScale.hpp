#pragma once

#include "Composition.hpp"
#include "UI.hpp"

#include <kit.hpp>

struct FreiScale : kit::Mode {
	FreiScale();
	virtual ~FreiScale();
	virtual void resized() override;
	virtual void update(float elapsed) override;
	virtual void draw() override;
	virtual void handle_event(SDL_Event const &) override;

	std::list< Sample > library;
	std::unique_ptr< Composition > composition;

	//UI state:
	UI ui;
	UIBox library_box;
	float library_top = 0.0f;

	UIBox song_box;
	TimeLog2Hz song_center = TimeLog2Hz(0.0f, 0.0f);
	TimeLog2Hz song_radius = TimeLog2Hz(8.0f, 2.5f);
};
