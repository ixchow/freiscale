#include "FreiScale.hpp"

#include "DrawLines.hpp"

#include <SDL_events.h>

FreiScale::FreiScale() {
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
}

void FreiScale::handle_event(SDL_Event const &evt) {
	if (evt.type == SDL_MOUSEMOTION) {
	}
	if (evt.type == SDL_MOUSEBUTTONDOWN) {
	}
}

