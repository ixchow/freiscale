#pragma once

#include <stdint.h>

typedef uint64_t FingerID;
const FingerID NoFinger = 0;
const FingerID MouseFinger = 1;

enum FingerAction {
	FingerDown,
	FingerMove,
	FingerUp,
	FingerCancel, //used when Mode becomes no longer current.
};
