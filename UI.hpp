#pragma once

#include "Touch.hpp"

#include <glm/glm.hpp>

#include <map>
#include <functional>

class UI;

class UIBox {
public:
	UIBox(UI *ui = NULL, glm::vec2 const &min = glm::vec2(0.0, 0.0), glm::vec2 const &max = glm::vec2(1.0, 1.0),
		std::function< void() > touch_up_inside = [](){});
	~UIBox();

	glm::vec2 min;
	glm::vec2 max;
	glm::vec2 size() const {
		return max - min;
	}
	glm::vec2 center() const {
		return 0.5f * (max + min);
	}
	bool contains(glm::vec2 const &pt) const {
		return (min.x <= pt.x && pt.x < max.x && min.y <= pt.y && pt.y < max.y);
	}

	std::function< void(FingerID, glm::vec2) > touch_down;
	std::function< void(FingerID, glm::vec2, glm::vec2) > touch_moved;
	std::function< void(FingerID) > touch_up;
	std::function< void() > touch_up_inside;

	uint32_t owned_touches;
	uint32_t inside_touches;

	UI *ui;
	UIBox *prev_box;
	UIBox *next_box;

	void remove_if_added();
};

class UIFingerInfo {
public:
	UIFingerInfo(glm::vec2 _at, UIBox *_owner = NULL, bool _inside = false) : at(_at), owner(_owner), inside(_inside) {
	}
	glm::vec2 at;
	UIBox *owner;
	bool inside;
};

class UI {
public:
	UI();
	~UI();

	void finger_action(FingerID finger, FingerAction action, glm::vec2 at);

	std::map< FingerID, UIFingerInfo > fingers;

	enum AddWhere {
		AtFront,
		AtBack,
	};

	void add_box(UIBox *, AddWhere where = AtFront);
	void remove_box(UIBox *);

	//Boxes get to handle touches in first-to-last order.

	UIBox *first_box;
	UIBox *last_box;
};
