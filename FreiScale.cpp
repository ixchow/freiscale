#include "FreiScale.hpp"

#include "DrawLines.hpp"
#include "Output.hpp"

#include <SDL_events.h>

#include <iostream>

FreiScale::FreiScale(std::string const &library_path) : library(library_path) {
	ui.add_box(&library_box);
	ui.add_box(&song_box);

	resized();
}

FreiScale::~FreiScale() {
}

void FreiScale::resized() {
	library_box.min = glm::vec2(0.0f, 0.0f);
	library_box.max = glm::vec2(std::min(std::round(0.4f * kit::display.window_size.x), 200.0f), kit::display.window_size.y);
	song_box.min = glm::vec2(library_box.max.x, 0.0f);
	song_box.max = glm::vec2(kit::display.window_size.x, kit::display.window_size.y);
}

void FreiScale::update(float elapsed) {
	update_hovered();
}

void FreiScale::draw() {

	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);

	{ //just something:
		DrawLines draw_lines(glm::mat4(
			1.0f, 0.0f, 0.0f, 0.0f,
			0.0f, 1.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
		));

		glm::vec2 px = glm::vec2(
			2.0f / kit::display.window_size.x,
			2.0f / kit::display.window_size.y
		);

		draw_lines.draw(glm::vec3(-1.0f,-1.0f, 0.0f), glm::vec3( 1.0f, 1.0f, 0.0f), glm::u8vec4(0xff, 0xff, 0x00, 0xff));
		draw_lines.draw(glm::vec3(-1.0f, 1.0f, 0.0f), glm::vec3( 1.0f,-1.0f, 0.0f), glm::u8vec4(0xff, 0x00, 0xff, 0xff));

		draw_lines.draw_text("FreiScale Begin",
			glm::vec3(0.0f, 0.0f, 0.0f),
			10.0f * glm::vec3(px.x, 0.0f, 0.0f),
			10.0f * glm::vec3(0.0f, px.y, 0.0f));

	}

	auto draw_panel = [&](UIBox const &box, glm::u8vec4 const &color) {
		DrawLines draw(glm::mat4(
			2.0f / kit::display.window_size.x, 0.0f, 0.0f, 0.0f,
			0.0f, 2.0f / kit::display.window_size.y, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			-1.0f, -1.0f, 0.0f, 1.0f
		));

		std::vector< glm::vec2 > outline{
			glm::vec2(box.min.x + 2.5f, box.min.y + 0.5f),
			glm::vec2(box.min.x + 0.5f, box.min.y + 2.5f),

			glm::vec2(box.min.x + 0.5f, box.max.y - 2.5f),
			glm::vec2(box.min.x + 2.5f, box.max.y - 0.5f),

			glm::vec2(box.max.x - 2.5f, box.max.y - 0.5f),
			glm::vec2(box.max.x - 0.5f, box.max.y - 2.5f),

			glm::vec2(box.max.x - 0.5f, box.min.y + 2.5f),
			glm::vec2(box.max.x - 2.5f, box.min.y + 0.5f),
		};

		for (uint32_t i = 0; i < outline.size(); ++i) {
			draw.draw(outline[i], outline[(i+1)%outline.size()], color);
		}
	};

	draw_panel(library_box, glm::u8vec4(0x88, 0x00, 0x88, 0xff));
	draw_panel(song_box, glm::u8vec4(0xff, 0xff, 0x00, 0xff));

	{ //song box grid:

		DrawLines draw(glm::mat4(
			2.0f / kit::display.window_size.x, 0.0f, 0.0f, 0.0f,
			0.0f, 2.0f / kit::display.window_size.y, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			-1.0f, -1.0f, 0.0f, 1.0f
		));
		
		{ //pitch scale grid:
			int32_t major = 10;
			int32_t min = int32_t(std::ceil((song_center.p - song_radius.p) * major));
			int32_t max = int32_t(std::floor((song_center.p + song_radius.p) * major));

			for (int32_t p = min; p <= max; ++p) {
				float y = ((p / float(major)) - song_center.p + song_radius.p) / (2.0f * song_radius.p) * (song_box.max.y - song_box.min.y) + song_box.min.y;
				draw.draw(glm::vec2(song_box.min.x, y), glm::vec2(song_box.max.x, y), 
					(p % major == 0 ? glm::u8vec4(0x88, 0x88, 0x88, 0xff) : glm::u8vec4(0x44, 0x44, 0x44, 0xff)) );
			}
		}

		{ //time grid:
			int32_t major = 4;
			int32_t min = int32_t(std::ceil((song_center.t - song_radius.t) * major));
			int32_t max = int32_t(std::floor((song_center.t + song_radius.t) * major));

			for (int32_t t = min; t <= max; ++t) {
				float x = ((t / float(major)) - song_center.t + song_radius.t) / (2.0f * song_radius.t) * (song_box.max.x - song_box.min.x) + song_box.min.x;
				draw.draw(glm::vec2(x, song_box.min.y), glm::vec2(x, song_box.max.y), 
					(t % major == 0 ? glm::u8vec4(0x88, 0x00, 0x88, 0xff) : glm::u8vec4(0x44, 0x00, 0x44, 0xff)) );
			}
		}
	}

	{ //library list:
		DrawLines draw(glm::mat4(
			2.0f / kit::display.window_size.x, 0.0f, 0.0f, 0.0f,
			0.0f, 2.0f / kit::display.window_size.y, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			-1.0f, -1.0f, 0.0f, 1.0f
		));

		uint32_t item = 0;

		library_folder_boxes.clear();
		library_sound_boxes.clear();

		auto draw_item = [&draw,this](uint32_t item, uint32_t depth, std::string const &text, bool highlight, UIBox *box) {
			const float Height = 16.0f;
			const float Margin = 2.0f;
			const float Indent = 5.0f;

			float y = library_box.max.y - (item + 1 - library_top) * (Height + Margin);
			float x = library_box.min.x + depth * Indent;

			draw.draw_text(text, glm::vec2(x,y), Height * glm::vec2(1.0f, 0.0f), Height * glm::vec2(0.0f, 1.0f),
				(highlight ? glm::u8vec4(0xff, 0xff, 0xff, 0xff) : glm::u8vec4(0x88, 0x88, 0x77, 0xff) )
				);

			if (box) {
				box->min.x = library_box.min.x;
				box->max.x = library_box.max.x;
				box->min.y = y - 0.5f * Margin;
				box->max.y = y + Height + 0.5f * Margin;
			}
		};

		std::function< void(uint32_t, Folder const &) > draw_folder;
		draw_folder = [&](uint32_t depth, Folder const &folder) {
			for (auto &[name, sub] : folder.folders) {
				//draw name...
				std::string draw_name;
				if (sub.state == Folder::Collapsed) {
					draw_name = "> " + name;
				} else if (sub.state == Folder::Expanded) {
					draw_name = "v " + name;
				} else {
					draw_name = "! " + name;
				}

				UIBox box;
				draw_item(item, depth, draw_name, (&sub == hovered.library_folder), &box);
				library_folder_boxes.emplace_back(box, &sub);

				++item;
				if (sub.state == Folder::Expanded) {
					draw_folder(depth + 1, sub);	
				}
			}
			for (auto &[name, sound] : folder.sounds) {
				UIBox box;
				draw_item(item, depth, name, (&sound == hovered.library_sound), &box);
				library_sound_boxes.emplace_back(box, &sound);
				++item;
			}
		};

		draw_folder(0, library.root);

	}
}

void FreiScale::handle_event(SDL_Event const &evt) {
	if (action) {
		action->handle_event(evt);
		return;
	}
	if (evt.type == SDL_MOUSEMOTION) {
		mouse = glm::vec2( evt.motion.x + 0.5f, (kit::display.window_size.y - 1 - evt.motion.y) + 0.5f );
	}
	if (evt.type == SDL_MOUSEBUTTONDOWN) {
		mouse = glm::vec2( evt.button.x + 0.5f, (kit::display.window_size.y - 1 - evt.motion.y) + 0.5f );
		update_hovered();
		if (hovered.library_sound) {
			if (evt.button.button == SDL_BUTTON_LEFT) {
				std::cout << "~~~Select Sound!~~~" << std::endl;
			} else if (evt.button.button == SDL_BUTTON_RIGHT) {
				std::cout << "~~~Play Sound!~~~" << std::endl;
				std::vector< Output::Sample > preview;
				preview.reserve(hovered.library_sound->size());
				for (auto s : *hovered.library_sound) {
					Output::Sample o;
					o.l = s;
					o.r = s;
					preview.emplace_back(o);
				}
				Output::lock();
				Output::playing_data = preview;
				Output::playing_position = 0;
				Output::unlock();
			}
		}
	}
	if (evt.type == SDL_MOUSEWHEEL) {
		if (library_box.contains(mouse)) {
			library_top -= evt.wheel.y;
			library_top = std::max(0.0f, library_top);
		}
	}
	if (evt.type == SDL_KEYDOWN) {
		if (evt.key.keysym.scancode == SDL_SCANCODE_SPACE) {
			Output::lock();
			if (Output::playing_position < Output::playing_data.size()) {
				std::cout << "Halting Playback." << std::endl;
				Output::playing_position = Output::playing_data.size();
			}
			Output::unlock();
		}
	}
}


void FreiScale::update_hovered() {
	hovered.clear();

	if (library_box.contains(mouse)) {
		for (auto const &[box, folder] : library_folder_boxes) {
			if (box.contains(mouse)) {
				hovered.library_folder = folder;
			}
		}
		for (auto const &[box, sound] : library_sound_boxes) {
			if (box.contains(mouse)) {
				hovered.library_sound = sound;
			}
		}
	}

}
