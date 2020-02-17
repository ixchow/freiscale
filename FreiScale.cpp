#include "FreiScale.hpp"

#include "DrawLines.hpp"
#include "Output.hpp"

#include <SDL_events.h>

#include <iostream>
#include <algorithm>

FreiScale::FreiScale(std::string const &library_path) : library(library_path) {
	resized();
	composition = std::make_unique< Composition >();
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

	{ //song box:
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

			min = std::max(0 * major, min);
			max = std::min(15 * major, max);

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

		{ //time markers:

			auto draw_marker = [&,this](float t, glm::u8vec4 const &color) {
				if (t < song_center.t - song_radius.t) return;
				if (t > song_center.t + song_radius.t) return;

				float x = ((t - song_center.t) / song_radius.t * 0.5f + 0.5f) * (song_box.max.x - song_box.min.x) + song_box.min.x;

				draw.draw( glm::vec2(x, song_box.min.y), glm::vec2(x, song_box.max.y), color );
			};
			

			draw_marker( composition->loop_begin, glm::u8vec4( 0xff, 0x00, 0x88, 0xff ) );
			draw_marker( composition->loop_end, glm::u8vec4( 0x88, 0x00, 0x44, 0xff ) );
		}

		//Triggers:
		for (auto const &t : composition->triggers) {

			{ //show sample line:
				glm::u8vec4 color = glm::u8vec4(0xaa, 0xcc, 0xbb, 0xff);
				glm::u8vec4 color_dim = glm::u8vec4(0x55, 0x77, 0x44, 0xff);
				const float Height = 10;

				TimeLog2Hz at = t.start;
				glm::vec2 px = get_screen_position(at);
				draw.draw( px + glm::vec2(0.0f, -0.5f * Height), px + glm::vec2(0.0f, 0.5f * Height), color );

				float remain = Time(t.sound ? t.sound->size() : 0) / Time(SampleRate);


				for (uint32_t s = 0; s < t.steps.size(); ++s) {
					float p0 = (s == 0 ? t.start.p : t.steps[s-1].p);
					float p1 = t.steps[s].p;
					float dt = t.steps[s].t;

					float a = (p1 - p0) / dt;
					float b = p0;

					float len = std::exp2( b ) / (std::log(2.0f) * a) * (std::exp2( a * dt ) - 1.0f);


					//TODO: figure out if sample lasts for entire segment...
					TimeLog2Hz next = at;
					next.t += dt;
					next.p = p1;
					glm::vec2 next_px = get_screen_position(next);

					if (len < remain) {
						remain -= len;

						draw.draw( px, next_px, color );
					} else if (remain > 0.0f) {
						//remain terminates somewhere in here!
						//TODO: do a two-color line to indicate where sample ends

					} else {
						//remain already terminated
						draw.draw( px, next_px, color_dim );
					}

					at = next;
					px = next_px;
				}

				float tail_px = 0.0f;
				if (remain > 0.0f) {
					float p = (t.steps.empty() ? t.start.p : t.steps.back().p);
					TimeLog2Hz next = at;
					next.t += remain / std::exp2( p );
					glm::vec2 next_px = get_screen_position(next);

					tail_px = next_px.x - px.x;

					draw.draw( px, next_px, color );
					draw.draw( next_px + glm::vec2(0.0f, -0.5f * Height), next_px + glm::vec2(0.0f, 0.5f * Height), color );
					
					px = next_px;
				}

				if (tail_px < 20.0f) {
					glm::vec2 next_px = px;
					next_px.x += (20.0f - tail_px);
					draw.draw( px, next_px, color_dim );
				}

				//TODO: draw enough of a tail to grab
			}

			{ //outline handles:
				glm::u8vec4 color = glm::u8vec4(0x76, 0x9a, 0x4b, 0xff);

				std::vector< glm::vec2 > box{
					glm::vec2(-2.0f, -2.0f),
					glm::vec2( 2.0f, -2.0f),
					glm::vec2( 2.0f,  2.0f),
					glm::vec2(-2.0f,  2.0f)
				};

				TimeLog2Hz at = t.start;

				glm::vec2 px = get_screen_position(at);
				for (uint32_t i = 0; i < box.size(); ++i) {
					draw.draw( px + box[i], px + box[(i+1)%box.size()], color );
				}
				for (uint32_t s = 0; s < t.steps.size(); ++s) {
					at.t += t.steps[s].t;
					at.p = t.steps[s].p;
					px = get_screen_position(at);
					for (uint32_t i = 0; i < box.size(); ++i) {
						draw.draw( px + box[i], px + box[(i+1)%box.size()], color );
					}
				}
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

	if (!show_sound.empty()) { //DEBUG: playback waveform
		DrawLines draw(glm::mat4(
			2.0f / float(show_sound.size()), 0.0f, 0.0f, 0.0f,
			0.0f, 0.5f, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			-1.0f, 0.0f, 0.0f, 1.0f
		));
		for (uint32_t i = 1; i < show_sound.size(); ++i) {
			draw.draw( glm::vec2( i-1, show_sound[i-1].l ), glm::vec2( i, show_sound[i].l ), glm::u8vec4( 0x88, 0x88, 0x88, 0xff ) );
		}
	}
}

struct MoveTriggerAction : public Action {
	MoveTriggerAction(FreiScale &fs_, Trigger &t_) : fs(fs_), t(t_), reference(fs.get_song_position(fs.mouse)) {
		original.emplace_back(t.start);
		for (auto const &s : t.steps) {
			original.emplace_back(s);
		}
	}
	virtual ~MoveTriggerAction() {
	}
	virtual void handle_event(SDL_Event const &evt) override {
		if (evt.type == SDL_MOUSEBUTTONUP) {
			if (evt.button.button == SDL_BUTTON_LEFT) {
				fs.action.reset();
				return;
			}
		}
		if (evt.type == SDL_KEYDOWN) {
			if (evt.key.keysym.scancode == SDL_SCANCODE_SPACE) {
				//TODO: preview playback trigger/cancel
			} else if (evt.key.keysym.scancode == SDL_SCANCODE_ESCAPE) {
				//TODO: cancel action
			}
		}
		if (evt.type == SDL_MOUSEMOTION) {
			fs.mouse = glm::vec2( evt.motion.x + 0.5f, (kit::display.window_size.y - 1 - evt.motion.y) + 0.5f );
			TimeLog2Hz current = fs.get_song_position(fs.mouse);

			TimeLog2Hz offset(current.t - reference.t, current.p - reference.p);

			t.start = original[0];
			t.start.t += offset.t;
			t.start.p += offset.p;
			for (uint32_t s = 0; s < t.steps.size(); ++s) {
				t.steps[s] = original[s+1];
				//t.steps[s].t += offset.t; <-- already relative
				t.steps[s].p += offset.p;
			}

		}
	}
	virtual void draw() {
	}
	FreiScale &fs;
	Trigger &t;
	std::vector< TimeLog2Hz > original;
	TimeLog2Hz reference;
};

struct MoveStepAction : public Action {
	MoveStepAction(FreiScale &fs_, Trigger &t_, uint32_t idx_) : fs(fs_), t(t_), idx(idx_), reference(fs.get_song_position(fs.mouse)) {
		original = t;
	}
	virtual ~MoveStepAction() {
	}
	virtual void handle_event(SDL_Event const &evt) override {
		if (evt.type == SDL_MOUSEBUTTONUP) {
			if (evt.button.button == SDL_BUTTON_LEFT) {
				fs.action.reset();
				return;
			}
		}
		if (evt.type == SDL_KEYDOWN) {
			if (evt.key.keysym.scancode == SDL_SCANCODE_SPACE) {
				//TODO: preview playback trigger/cancel
			} else if (evt.key.keysym.scancode == SDL_SCANCODE_ESCAPE) {
				//cancel action
				t = original;
				fs.action.reset();
			}
		}
		if (evt.type == SDL_MOUSEMOTION) {
			fs.mouse = glm::vec2( evt.motion.x + 0.5f, (kit::display.window_size.y - 1 - evt.motion.y) + 0.5f );
			TimeLog2Hz current = fs.get_song_position(fs.mouse);

			TimeLog2Hz offset(current.t - reference.t, current.p - reference.p);

			t = original;

			if (idx == 0) {
				//moving first step == 'start'
			} else if (0 < idx && idx <= t.steps.size()) {
				float old = t.steps[idx-1].t;
				t.steps[idx-1].t = std::max(0.0f, t.steps[idx-1].t + offset.t);
				t.steps[idx-1].p += offset.p;
				if (idx < t.steps.size()) {
					t.steps[idx].t = std::max(0.0f, t.steps[idx].t - (t.steps[idx-1].t - old));
				}
			} else {
				//??? bad index
			}
		}
	}
	virtual void draw() {
	}
	FreiScale &fs;
	Trigger &t;
	uint32_t idx; //0 -> start, 1 .. N -> steps

	Trigger original;
	TimeLog2Hz reference;
};


void FreiScale::handle_event(SDL_Event const &evt) {
	if (action) {
		action->handle_event(evt);
		return;
	}
	if (evt.type == SDL_MOUSEMOTION) {
		mouse = glm::vec2( evt.motion.x + 0.5f, (kit::display.window_size.y - 1 - evt.motion.y) + 0.5f );
		return;
	}
	if (evt.type == SDL_MOUSEBUTTONDOWN) {
		mouse = glm::vec2( evt.button.x + 0.5f, (kit::display.window_size.y - 1 - evt.motion.y) + 0.5f );
		update_hovered();
		if (song_box.contains(mouse)) {
			//------ song interactions -----
			if (evt.button.button == SDL_BUTTON_LEFT) {
				auto mod_state = SDL_GetModState();
				if (mod_state & KMOD_CTRL) {
					composition->loop_begin = get_song_position(mouse).t;
				} else if (mod_state & KMOD_ALT) {
					composition->loop_end = get_song_position(mouse).t;
				} else {
					Trigger *song_trigger = nullptr;
					if (hovered.song_trigger_segment.first) song_trigger = hovered.song_trigger_segment.first;
					if (hovered.song_trigger_handle.first) song_trigger = hovered.song_trigger_handle.first;
					if (song_trigger) {
						action.reset(new MoveTriggerAction(*this, *song_trigger));
						return;
					}
				}
			} else if (evt.button.button == SDL_BUTTON_MIDDLE) {
				//TODO: pan view
			}
		} else if (library_box.contains(mouse)) {
			//------ library interactions -----
			if (hovered.library_sound) {
				if (evt.button.button == SDL_BUTTON_LEFT) {
					std::cout << "~~~Select Sound!~~~" << std::endl;
					if (hovered.library_sound) {
						current_library_sound = hovered.library_sound;
					}
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
					show_sound = preview;
					Output::lock();
					Output::playing_data = preview;
					Output::playing_position = 0;
					Output::unlock();
				}
			}
		}
	}
	if (evt.type == SDL_MOUSEWHEEL) {
		if (library_box.contains(mouse)) {
			library_top -= evt.wheel.y;
			library_top = std::max(0.0f, library_top);
		} else if (song_box.contains(mouse)) {
			TimeLog2Hz focus = get_song_position(mouse);

			if (SDL_GetModState() & KMOD_SHIFT) {
				song_radius.p = std::min(16.0f, std::max(0.1f, song_radius.p * std::exp2(0.25f * evt.wheel.y)));
			} else {
				song_radius.t = std::min(100.0f, std::max(0.1f, song_radius.t * std::exp2(0.25f * evt.wheel.y)));
			}

			TimeLog2Hz new_focus = get_song_position(mouse);

			song_center.t += focus.t - new_focus.t;
			song_center.p += focus.p - new_focus.p;
		}
	}
	if (evt.type == SDL_KEYDOWN) {
		if (evt.key.keysym.sym == SDLK_F4) {
			//save:
			std::string filename = "temp.fs";
			composition->save(filename);
			std::cout << "Wrote '" << filename << "'." << std::endl;
		} else if (evt.key.keysym.sym == SDLK_c) {
			//"create" new trigger
			if (song_box.contains(mouse) && current_library_sound) {
				std::cout << "Make Trigger!" << std::endl;
				Trigger trigger;
				trigger.start = get_song_position(mouse);
				trigger.sound = composition->add_sound(*current_library_sound);
				composition->triggers.emplace_back(trigger);

				action.reset(new MoveTriggerAction(*this, composition->triggers.back()));
				return;
			}
		} else if (evt.key.keysym.sym == SDLK_s) {
			//add new step to trigger
			if (song_box.contains(mouse) && hovered.song_trigger_segment.first) {
				auto &trigger = *hovered.song_trigger_segment.first;
				TimeLog2Hz pos = get_song_position(mouse);
				float t = pos.t - trigger.start.t;
				uint32_t before = 0;
				while (before < trigger.steps.size() && t > trigger.steps[before].t) {
					t -= trigger.steps[before].t;
					++before;
				}
				trigger.steps.insert(trigger.steps.begin() + before, TimeLog2Hz(0, pos.p));

				action.reset(new MoveStepAction(*this, *hovered.song_trigger_segment.first, before));
				return;
			}
		} else if (evt.key.keysym.sym == SDLK_g) {
			//grab step in trigger
			if (song_box.contains(mouse) && hovered.song_trigger_handle.first) {
				action.reset(new MoveStepAction(*this, *hovered.song_trigger_handle.first, hovered.song_trigger_handle.second));
				return;
			}
		} else if (evt.key.keysym.scancode == SDL_SCANCODE_SPACE) {
			bool start_playback = true;
			Output::lock();
			if (Output::playing_position < Output::playing_data.size()) {
				std::cout << "Halting Playback." << std::endl;
				Output::playing_position = Output::playing_data.size();
				start_playback = false;
			}
			Output::unlock();

			if (start_playback) {
				if (song_box.contains(mouse)) {
					std::cout << "Trying to play loop..." << std::endl;
					int32_t loop_begin = composition->loop_begin * SampleRate;
					int32_t loop_end = composition->loop_end * SampleRate;
					if (loop_begin < loop_end) {
						//HACK: don't render more than 60 seconds for now...
						loop_end = std::min(loop_end, loop_begin + 60 * int32_t(SampleRate));

						std::cout << "  rendering..." << std::endl;
						std::vector< Sample > buffer;
						composition->render(loop_begin, loop_end, &buffer);

						std::cout << "  analyzing..." << std::endl;
						Sound temp = Sound::from_samples(buffer.data(), buffer.data() + buffer.size());
						temp.compute_spectrums();

						std::cout << "  playing..." << std::endl;
						std::vector< Output::Sample > preview;
						preview.reserve(buffer.size());
						for (auto s : buffer) {
							Output::Sample o;
							o.l = s;
							o.r = s;
							preview.emplace_back(o);
						}
						show_sound = preview;
						Output::lock();
						Output::playing_data = preview;
						Output::playing_position = 0;
						Output::unlock();
					}
				}
			} else {
				show_sound.clear();
			}
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
	} else if (song_box.contains(mouse)) {
		float close_handle = 10.0f;
		float close_segment = 10.0f;

		auto check_handle = [&](glm::vec2 const &px, Trigger &t, uint32_t idx) {
			float d = std::max(std::abs(px.x - mouse.x), std::abs(px.y - mouse.y));
			if (d < close_handle) {
				close_handle = d;
				hovered.song_trigger_handle = std::make_pair(&t, idx);
			}
		};

		auto check_segment = [&](glm::vec2 const &a, glm::vec2 const &b, Trigger &t, uint32_t idx) {
			float along = glm::dot(mouse - a, b - a);
			along = std::max(along, 0.0f);
			along = std::min(along, glm::dot(b-a, b-a));
			glm::vec2 px = (along / glm::dot(b-a, b-a)) * (b-a) + a;
			float d = std::max(std::abs(px.x - mouse.x), std::abs(px.y - mouse.y));
			if (d < close_segment) {
				close_segment = d;
				hovered.song_trigger_segment = std::make_pair(&t, idx);
			}
		};
		for (auto &t : composition->triggers) {
			check_handle(get_screen_position(t.start), t, 0);
			TimeLog2Hz at = t.start;
			TimeLog2Hz prev = at;
			for (uint32_t s = 0; s < t.steps.size(); ++s) {
				at.t += t.steps[s].t;
				at.p = t.steps[s].p;
				check_handle(get_screen_position(at), t, 1 + s);
				check_segment(get_screen_position(prev), get_screen_position(at), t, s);
				prev = at;
			}
			//TODO: correct tail length(!!)
			check_segment(get_screen_position(prev), get_screen_position(prev) + glm::vec2(20.0f, 0.0f), t, t.steps.size());
		}
	}

}
