#include "UI.hpp"

#include <iostream>

UIBox::UIBox(UI *_ui, glm::vec2 const &_min, glm::vec2 const &_max,
	std::function< void() > _touch_up_inside) :
		min(_min), max(_max),
		touch_down( [](FingerID, glm::vec2){} ),
		touch_moved( [](FingerID, glm::vec2, glm::vec2){} ),
		touch_up( [](FingerID){} ),
		touch_up_inside(_touch_up_inside),
		owned_touches(0), inside_touches(0),
		ui(NULL), prev_box(NULL), next_box(NULL) {
	assert(min.x <= max.x);
	assert(min.y <= max.y);
	if (_ui) {
		_ui->add_box(this);
		assert(ui);
	}
}
UIBox::~UIBox() {
	remove_if_added();

	assert(!ui);
	assert(!next_box);
	assert(!prev_box);
}

void UIBox::remove_if_added() {
	if (ui) {
		ui->remove_box(this);
	}
}

UI::UI() : first_box(NULL), last_box(NULL) {
}

UI::~UI() {
	while (last_box) {
		remove_box(last_box);
	}
	assert(!first_box);
}

void UI::finger_action(FingerID finger, FingerAction action, glm::vec2 at) {
	if (action == FingerDown) {
		UIFingerInfo info(at);
		for (UIBox *box = first_box; box != NULL; box = box->next_box) {
			if (at.x >= box->min.x && at.x <= box->max.x
			 && at.y >= box->min.y && at.y <= box->max.y) {
				info.owner = box;
				info.inside = true;
				break;
			}
		}
		if (info.owner) {
			info.owner->owned_touches += 1;
			info.owner->inside_touches += 1;
			info.owner->touch_down(finger, at);
		}
		fingers.insert(std::make_pair(finger, info));
	} else {
		auto f = fingers.find(finger);
		if (f == fingers.end()) {
			std::cerr << "Odd, missing finger." << std::endl;
			return;
		}
		UIFingerInfo &info = f->second;
		if (action == FingerMove) {
			if (info.owner) {
				bool was_inside = info.inside;
				if (at.x >= info.owner->min.x && at.x <= info.owner->max.x
				 && at.y >= info.owner->min.y && at.y <= info.owner->max.y) {
				 	if (!was_inside) {
						info.inside = true;
						info.owner->inside_touches += 1;
						assert(info.owner->inside_touches <= info.owner->owned_touches);
					}
				} else {
					if (was_inside) {
						info.inside = false;
						assert(info.owner->inside_touches > 0);
						info.owner->inside_touches -= 1;
					}
				}
				info.owner->touch_moved(finger, info.at, at);
				info.at = at;
			}
		} else if (action == FingerUp || action == FingerCancel) {
			if (info.owner) {
				if (info.inside) {
					assert(info.owner->inside_touches > 0);
					info.owner->inside_touches -= 1;
				}
				assert(info.owner->owned_touches > 0);
				info.owner->owned_touches -= 1;
				assert(info.owner->inside_touches <= info.owner->owned_touches);
			}
			UIBox *owner = info.owner;
			bool call_touch_up_inside = (action == FingerUp && info.inside);
			fingers.erase(f);
			if (owner) {
				owner->touch_up(finger);
				if (call_touch_up_inside) {
					owner->touch_up_inside();
				}
			}
		} else {
			assert(0 && "unknown finger action");
		}
	}

}

void UI::add_box(UIBox *box, AddWhere where) {
	assert(box);

	if (box->ui) {
		box->ui->remove_box(box);
	}
	assert(!box->ui);
	assert(!box->prev_box);
	assert(!box->next_box);
	assert(box->owned_touches == 0);
	assert(box->inside_touches == 0);

	box->ui = this;

	if (!first_box) {
		assert(!last_box);
		first_box = last_box = box;
		return;
	}
	if (where == AtFront) {
		assert(first_box);
		assert(!first_box->prev_box);
		first_box->prev_box = box;
		box->next_box = first_box;
		first_box = box;
	} else {
		assert(last_box);
		assert(!last_box->next_box);
		last_box->next_box = box;
		box->prev_box = last_box;
		last_box = box;
	}
}

void UI::remove_box(UIBox *box) {
	assert(box);
	assert(box->ui == this);

	box->owned_touches = 0;
	box->inside_touches = 0;
	//slow, yeah, but this probably doesn't happen much:
	for (auto f = fingers.begin(); f != fingers.end(); ++f) {
		if (f->second.owner == box) {
			f->second.owner = NULL;
		}
	}

	if (box->prev_box == NULL) {
		assert(first_box == box);
		first_box = box->next_box;
	} else {
		box->prev_box->next_box = box->next_box;
	}
	if (box->next_box == NULL) {
		assert(last_box == box);
		last_box = box->prev_box;
	} else {
		box->next_box->prev_box = box->prev_box;
	}
	box->ui = NULL;
	box->next_box = NULL;
	box->prev_box = NULL;
}
