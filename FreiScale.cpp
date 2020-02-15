#include "FreiScale.hpp"

#include "DrawLines.hpp"

#include <SDL_events.h>

FreiScale::FreiScale() {
}

FreiScale::~FreiScale() {
}

void FreiScale::resized() {
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
}

void FreiScale::handle_event(SDL_Event const &evt) {
	if (evt.type == SDL_MOUSEMOTION) {
	}
	if (evt.type == SDL_MOUSEBUTTONDOWN) {
	}
}

