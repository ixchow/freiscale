#include "FreiScale.hpp"

#include "DrawLines.hpp"
#include "DrawStuff.hpp"
#include "Output.hpp"
#include "SpectrumProgram.hpp"

#include <kit/gl_errors.hpp>

#include <SDL_events.h>

#include <iostream>
#include <algorithm>


static GLuint vertex_buffer = 0;
static GLuint vertex_buffer_for_spectrum_program = 0;

struct SpectrumVertex {
	glm::vec2 Position;
	glm::vec2 TexCoord;
};

static kit::Load< void > setup_buffers(kit::LoadTagDefault, [](){
	//you may recognize this init code from DrawLines.cpp:

	glGenBuffers(1, &vertex_buffer);

	glGenVertexArrays(1, &vertex_buffer_for_spectrum_program);
	glBindVertexArray(vertex_buffer_for_spectrum_program);

	glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);

	glVertexAttribPointer(
		spectrum_program->Position_vec4, //attribute
		2, //size
		GL_FLOAT, //type
		GL_FALSE, //normalized
		sizeof(DrawLines::Vertex), //stride
		(GLbyte *)0 + offsetof(SpectrumVertex, Position) //offset
	);
	glEnableVertexAttribArray(spectrum_program->Position_vec4);

	glVertexAttribPointer(
		spectrum_program->TexCoord_vec2, //attribute
		2, //size
		GL_FLOAT, //type
		GL_FALSE, //normalized
		sizeof(DrawLines::Vertex), //stride
		(GLbyte *)0 + offsetof(SpectrumVertex, TexCoord) //offset
	);
	glEnableVertexAttribArray(spectrum_program->TexCoord_vec2);


	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindVertexArray(0);

	GL_ERRORS(); //PARANOIA: make sure nothing strange happened during setup
});



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
	for (auto const &t : composition->triggers) {
		if (t->sources_dirty) {
			t->compute_sources();
		}
	}
	composition->update_rendered(composition->loop_begin);
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

	glm::mat4 px_to_clip = glm::mat4(
		2.0f / kit::display.window_size.x, 0.0f, 0.0f, 0.0f,
		0.0f, 2.0f / kit::display.window_size.y, 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		-1.0f, -1.0f, 0.0f, 1.0f
	);

	//DEBUG: draw this under the grid:
	/*
	{
		if (spectrum_tex) {
			glm::vec2 min = get_screen_position(
				TimeLog2Hz( spectrum_tex_t0, std::log2( SpectrumMinFreq ) )
			);
			glm::vec2 max = get_screen_position(
				TimeLog2Hz( spectrum_tex_t1, std::log2( SpectrumMaxFreq ) )
			);

			//DEBUG:
			//min = song_box.min;
			//max = song_box.max;

			std::vector< SpectrumVertex > attribs{
				{glm::vec2(min.x, min.y), glm::vec2(0.0f, 0.0f ) },
				{glm::vec2(min.x, max.y), glm::vec2(1.0f, 0.0f ) },
				{glm::vec2(max.x, min.y), glm::vec2(0.0f, 1.0f ) },
				{glm::vec2(max.x, max.y), glm::vec2(1.0f, 1.0f ) }
			};
			glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
			glBufferData(GL_ARRAY_BUFFER, attribs.size() * sizeof(attribs[0]), attribs.data(), GL_STREAM_DRAW);
			glBindBuffer(GL_ARRAY_BUFFER, 0);

			glUseProgram(spectrum_program->program);

			glBindVertexArray(vertex_buffer_for_spectrum_program);

			glBindTexture(GL_TEXTURE_2D, spectrum_tex);

			glUniformMatrix4fv(spectrum_program->OBJECT_TO_CLIP_mat4, 1, GL_FALSE, glm::value_ptr(px_to_clip));

			glDrawArrays(GL_TRIANGLE_STRIP, 0, GLsizei(attribs.size()));

			glBindTexture(GL_TEXTURE_2D, 0);

			glBindVertexArray(0);

			glUseProgram(0);

			GL_ERRORS();
		}
	}*/

	{ //song box ready:
		std::vector< DrawStuff::Pos2f_Col4ub > tristrip;

		for (auto const &[ idx, ptr ] : composition->rendered) {
			glm::vec2 min = get_screen_position(TimeLog2Hz( ptr->start_sample / float(SampleRate), song_center.p - song_radius.p));
			glm::vec2 max = get_screen_position(TimeLog2Hz( (ptr->start_sample + Composition::BlockSize) / float(SampleRate), song_center.p + song_radius.p));

			glm::u8vec4 colorB, colorT;
			if (ptr->dirty) {
				colorB = glm::u8vec4(0x55, 0x11, 0x11, 0xff);
				colorT = glm::u8vec4(0x44, 0x44, 0x44, 0xff);
			} else {
				colorB = glm::u8vec4(0x22, 0x22, 0x22, 0xff);
				colorT = glm::u8vec4(0x44, 0x44, 0x66, 0xff);
			}

			if (!tristrip.empty()) tristrip.emplace_back(tristrip.back());
			tristrip.emplace_back(glm::vec2(min.x, min.y), colorB);
			if (tristrip.size() != 1) tristrip.emplace_back(tristrip.back());
			tristrip.emplace_back(glm::vec2(max.x, min.y), colorB);
			tristrip.emplace_back(glm::vec2(min.x, max.y), colorT);
			tristrip.emplace_back(glm::vec2(max.x, max.y), colorT);
		}

		DrawStuff::draw(px_to_clip, GL_TRIANGLE_STRIP, tristrip);
	}

	{ //song box background:
		DrawLines draw(px_to_clip);
			
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
	}



	{ //song box triggers:
		std::vector< DrawStuff::Pos2f_Col4ub > lines;
		std::vector< DrawStuff::Pos2f_Col4ub > points;
	
		//Triggers:
		for (auto const &t : composition->triggers) {
			//float sample_length = Time(t->sound ? t->sound->size() : 0) / Time(SampleRate);
			//uint32_t peaks_count = (t->sound ? t->sound->peaks.size() : 0) / PeaksSlots;

			float offset_p = (t->sound ? std::log2( t->sound->fundamental ) : 0.0f );

			{ //show sample line:
				glm::u8vec4 color = glm::u8vec4(0xaa, 0xcc, 0xbb, 0xff);
				glm::u8vec4 color_dim = glm::u8vec4(0x55, 0x77, 0x44, 0xff);
				const float Height = 10;

				TimeLog2Hz at = t->steps[0];
				at.p += offset_p;
				glm::vec2 at_px = get_screen_position(at);
				lines.emplace_back( at_px + glm::vec2(0.0f, -0.5f * Height), color);
				lines.emplace_back( at_px + glm::vec2(0.0f, 0.5f * Height), color );

				for (uint32_t s = 1; s < t->steps.size(); ++s) {
					TimeLog2Hz next = t->steps[s];
					next.p += offset_p;
					glm::vec2 next_px = get_screen_position(next);

					/* TODO: bring back sample end indication somehow.
					int32_t s0 = int32_t(std::round(at.t * SampleRate));
					int32_t s1 = int32_t(std::round(at.t * SampleRate));
					s1 = std::max(s0, s1);
					*/
					lines.emplace_back( at_px, color_dim );
					lines.emplace_back( next_px, color_dim );

					at = next;
					at_px = next_px;
				}


				/*

				std::vector< TimeLog2Hz > warped_peaks;
				warped_peaks.reserve(peaks_count);

				auto warp_peaks = [&warped_peaks,&peaks_count,&offset_p](TimeLog2Hz tp0, TimeLog2Hz tp1, float s0) {
					//draw peaks from tp0 to tp1 given sample time s0:

					float p0 = tp0.p;
					float p1 = tp1.p;
					if (p1 == p0) p1 += 1e-3;
					float t0 = tp0.t;
					float t1 = tp1.t;

					float a = (p1 - p0) / (t1 - t0);
					float b = p0 + offset_p;

					
					//float len = std::exp2( b ) / (std::log(2.0f) * a) * (std::exp2( a * (t1-t0) ) - 1.0f);

					//std::cout << "[" << t0 << " - " << t1 << "] => " << p0 << " - " << p1 << " ==> " << len << std::endl; //DEBUG

					while (warped_peaks.size() < peaks_count) {
						float ds = warped_peaks.size() / float(PeaksRate) - s0;
						float t = std::log2( ds * ( (std::log(2.0f) * a) / std::exp2( b ) ) + 1.0f ) / a;
						//std::cout << "ds: " << ds << " -> t " << t << " of " << len << std::endl;
						if (t > t1 - t0) break;
						warped_peaks.emplace_back(t0 + t, a * t + b);
					}
				};

				float play_t = 0.0f;

				for (uint32_t s = 0; s < t.steps.size(); ++s) {
					TimeLog2Hz next = t.steps[s];
					next.t += at.t;
					glm::vec2 next_px = get_screen_position(next);

					if (next.t == at.t) {
						at = next;
						px = next_px;
						continue;
					}

					warp_peaks(at, next, play_t);

					//draw line:

					float p0 = at.p;
					float p1 = next.p;
					float dt = t.steps[s].t;

					float a = (p1 - p0) / dt;
					float b = p0 + offset_p;

					float len = std::exp2( b ) / (std::log(2.0f) * a) * (std::exp2( a * dt ) - 1.0f);

					if (play_t > sample_length) {
						//already done
						lines.emplace_back( px, color_dim );
						lines.emplace_back( next_px, color_dim );
					} else if (play_t + len > sample_length) {
						//going to finish in here
						//TODO: fancy hard edge
						lines.emplace_back( px, color );
						lines.emplace_back( next_px, color_dim );
					} else {
						//not going to finish in this segment
						lines.emplace_back( px, color );
						lines.emplace_back( next_px, color );
					}

					//update playback position etc:

					play_t += len;
					at = next;
					px = next_px;
				}

				float tail_px = 0.0f;
				if (play_t < sample_length) {

					float p = (t.steps.empty() ? t.start.p : t.steps.back().p);
					TimeLog2Hz next = at;
					next.t += (sample_length - play_t) / std::exp2( p + offset_p );
					glm::vec2 next_px = get_screen_position(next);

					warp_peaks(at, next, play_t);

					tail_px = next_px.x - px.x;

					lines.emplace_back( px, color );
					lines.emplace_back( next_px, color );
					lines.emplace_back( next_px + glm::vec2(0.0f, -0.5f * Height), color );
					lines.emplace_back( next_px + glm::vec2(0.0f, 0.5f * Height), color );
					
					px = next_px;
				}

				if (tail_px < 20.0f) {
					glm::vec2 next_px = px;
					next_px.x += (20.0f - tail_px);
					lines.emplace_back( px, color_dim );
					lines.emplace_back( next_px, color_dim );
				}

				//peaks:
				for (uint32_t peak = 0; peak < warped_peaks.size(); ++peak) {
					std::pair< float, float > const *slots = &(t.sound->peaks[PeaksSlots * peak]);
					for (uint32_t s = 0; s < PeaksSlots; ++s) {
						if (slots[s].first <= 0.0f) continue;
						//std::cout << slots[s].first << "/" << warped_peaks[peak].t << "/" << warped_peaks[peak].p << std::endl;
						int32_t amt = 255.0f * slots[s].second;
						amt = std::max(0, std::min(255, amt));
						glm::u8vec4 color = glm::u8vec4(0x0, amt, 0x0, 0xff);

						if(&t == hovered.song_trigger_handle.first || &t == hovered.song_trigger_segment.first) {
							color.r = color.g;
						}

						TimeLog2Hz from(warped_peaks[peak].t, glm::log2(slots[s].first) + warped_peaks[peak].p);
						TimeLog2Hz to = from;
						if (peak + 1 < warped_peaks.size()) {
							to.t = warped_peaks[peak+1].t;
							to.p = std::log2(slots[s].first) + warped_peaks[peak+1].p;
						}

						lines.emplace_back( get_screen_position(from), color );
						lines.emplace_back( get_screen_position(to), color );
					}
				}
				*/
			}

			{ //outline handles:
				glm::u8vec4 color = glm::u8vec4(0x76, 0x9a, 0x4b, 0xff);

				std::vector< glm::vec2 > box{
					glm::vec2(-2.0f, -2.0f),
					glm::vec2( 2.0f, -2.0f),
					glm::vec2( 2.0f,  2.0f),
					glm::vec2(-2.0f,  2.0f)
				};

				for (uint32_t s = 0; s < t->steps.size(); ++s) {
					TimeLog2Hz at = t->steps[s];
					at.p += offset_p;
					glm::vec2 px = get_screen_position(at);
					for (uint32_t i = 0; i < box.size(); ++i) {
						lines.emplace_back( px + box[i], color );
						lines.emplace_back( px + box[(i+1)%box.size()], color );
					}
				}
			}

			DrawStuff::draw(px_to_clip, GL_POINTS, points);
			DrawStuff::draw(px_to_clip, GL_LINES, lines);

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

		std::function< void(uint32_t, Folder &) > draw_folder;
		draw_folder = [&](uint32_t depth, Folder &folder) {
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

	/*
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
	*/
}

struct PanViewAction : public Action {
	PanViewAction(FreiScale &fs_) : fs(fs_), reference(fs.get_song_position(fs.mouse)) {
	}
	virtual ~PanViewAction() {
	}
	virtual void handle_event(SDL_Event const &evt) override {
		if (evt.type == SDL_MOUSEBUTTONUP) {
			if (evt.button.button == SDL_BUTTON_MIDDLE) {
				fs.action.reset();
				return;
			}
		}
		if (evt.type == SDL_MOUSEMOTION) {
			fs.mouse = glm::vec2( evt.motion.x + 0.5f, (kit::display.window_size.y - 1 - evt.motion.y) + 0.5f );
			TimeLog2Hz current = fs.get_song_position(fs.mouse);

			TimeLog2Hz offset(current.t - reference.t, current.p - reference.p);

			fs.song_center.t -= offset.t;
			fs.song_center.p -= offset.p;

		}
	}
	virtual void draw() override {
	}
	FreiScale &fs;
	TimeLog2Hz reference;
};


struct MoveTriggerAction : public Action {
	MoveTriggerAction(FreiScale &fs_, Trigger &t_) : fs(fs_), t(t_), reference(fs.get_song_position(fs.mouse)) {
		original = t.steps;
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
				t.steps = original;
				t.sources_dirty = true;
				fs.action.reset();
			}
		}
		if (evt.type == SDL_MOUSEMOTION) {
			fs.mouse = glm::vec2( evt.motion.x + 0.5f, (kit::display.window_size.y - 1 - evt.motion.y) + 0.5f );
			TimeLog2Hz current = fs.get_song_position(fs.mouse);

			TimeLog2Hz offset(current.t - reference.t, current.p - reference.p);

			t.steps = original;
			for (uint32_t s = 0; s < t.steps.size(); ++s) {
				t.steps[s].t += offset.t;
				t.steps[s].p += offset.p;
			}
			t.sources_dirty = true;
		}
	}
	virtual void draw() override {
	}
	FreiScale &fs;
	Trigger &t;
	std::vector< TimeLog2Hz > original;
	TimeLog2Hz reference;
};

struct MoveStepAction : public Action {
	MoveStepAction(FreiScale &fs_, Trigger &t_, uint32_t idx_) : fs(fs_), t(t_), idx(idx_), reference(fs.get_song_position(fs.mouse)) {
		original = t.steps;
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
				t.steps = original;
				t.sources_dirty = true;
				fs.action.reset();
			}
		}
		if (evt.type == SDL_MOUSEMOTION) {
			fs.mouse = glm::vec2( evt.motion.x + 0.5f, (kit::display.window_size.y - 1 - evt.motion.y) + 0.5f );
			TimeLog2Hz current = fs.get_song_position(fs.mouse);

			TimeLog2Hz offset(current.t - reference.t, current.p - reference.p);

			t.steps = original;
			if (idx < t.steps.size()) {
				t.steps[idx].t = original[idx].t + offset.t;
				t.steps[idx].p = original[idx].p + offset.p;

				if (idx > 0 && t.steps[idx].t < t.steps[idx-1].t) {
					float delta = t.steps[idx].t - t.steps[idx-1].t;
					for (uint32_t i = 0; i < idx; ++i) {
						t.steps[i].t += delta;
					}
				}

				if (idx + 1 < t.steps.size() && t.steps[idx].t > t.steps[idx+1].t) {
					float delta = t.steps[idx].t - t.steps[idx+1].t;
					for (uint32_t i = idx + 1; i < t.steps.size(); ++i) {
						t.steps[i].t += delta;
					}
				}

				t.sources_dirty = true;
			} else {
				//??? bad index
			}
		}
	}
	virtual void draw() override {
	}
	FreiScale &fs;
	Trigger &t;
	uint32_t idx; //0 -> start, 1 .. N -> steps

	std::vector< TimeLog2Hz > original;
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
				action.reset(new PanViewAction(*this));
			}
		} else if (library_box.contains(mouse)) {
			//------ library interactions -----
			if (hovered.library_sound) {
				if (hovered.library_sound->peaks.empty()) {
					hovered.library_sound->compute_viz();
				}
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
				song_radius.t = std::min(100.0f, std::max(0.1f, song_radius.t * std::exp2(0.25f * evt.wheel.y)));
			} else {
				song_radius.p = std::min(16.0f, std::max(0.1f, song_radius.p * std::exp2(0.25f * evt.wheel.y)));
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
				std::shared_ptr< Trigger > trigger = std::make_shared< Trigger >(composition->add_sound(*current_library_sound));
				trigger->steps[0] = get_song_position(mouse);
				trigger->steps[0].p -= std::log2(trigger->sound->fundamental);
				trigger->steps[1] = trigger->steps[0];
				trigger->steps[1].t += 1.0f;
				composition->triggers.emplace_back(trigger);

				action.reset(new MoveTriggerAction(*this, *composition->triggers.back()));
				return;
			}
		} else if (evt.key.keysym.sym == SDLK_s) {
			//add new step to trigger
			if (song_box.contains(mouse) && hovered.song_trigger_segment.first) {
				auto &trigger = *hovered.song_trigger_segment.first;
				TimeLog2Hz pos = get_song_position(mouse);
				uint32_t before = 0;
				while (before < trigger.steps.size() && pos.t > trigger.steps[before].t) {
					++before;
				}
				trigger.steps.insert(trigger.steps.begin() + before, TimeLog2Hz(pos.t, pos.p - std::log2(trigger.sound->fundamental)));

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

/*
						std::cout << "  analyzing..." << std::endl;
						Sound temp = Sound::from_samples(buffer.data(), buffer.data() + buffer.size());
						temp.compute_spectrums();

						if (spectrum_tex == 0) glGenTextures(1, &spectrum_tex);
						glBindTexture(GL_TEXTURE_2D, spectrum_tex);
						glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, SpectrumBins, temp.spectrums.size() / SpectrumBins, 0, GL_RED, GL_FLOAT, temp.spectrums.data());

						glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
						glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

						glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
						glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

						glGenerateMipmap(GL_TEXTURE_2D);

						glBindTexture(GL_TEXTURE_2D, 0);

						spectrum_tex_t0 = composition->loop_begin;
						spectrum_tex_t1 = composition->loop_end;
*/

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
			float offset_p = std::log2(t->sound->fundamental);
			glm::vec2 prev_px = glm::vec2(0.0f);
			for (uint32_t s = 0; s < t->steps.size(); ++s) {
				glm::vec2 px = get_screen_position(TimeLog2Hz(t->steps[s].t, t->steps[s].p + offset_p));
				check_handle(px, *t, s);
				if (s > 0) {
					check_segment(prev_px, px, *t, s);
				}
				prev_px = px;
			}
		}
	}

}
