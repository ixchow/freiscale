#pragma once

/*
 * PathFont -- line based font used by DrawLines.
 *
 * Based on code from Chesskoban (c) 2017-2019 Jim McCann;
 * this adapted-for-15-466 code is released into the public domain.
 *
 */

#include <glm/glm.hpp>

#include <string>
#include <vector>
#include <map>

struct PathFont {
	//meant to be intitialized with some pointers to constant data:
	PathFont(uint32_t glyphs,
		const float *glyph_widths,
		const uint32_t *glyph_char_starts, const uint8_t *chars,
		const uint32_t *glyph_coord_starts, const float *coords
		);
	const uint32_t glyphs = 0;
	const float *glyph_widths = nullptr;

	const uint32_t *glyph_char_starts = nullptr; //indices into 'chars' table
	const uint8_t *chars = nullptr;

	const uint32_t *glyph_coord_starts = nullptr; //indices into 'coords' table
	const float *coords = nullptr;

	//computed in constructor:
	std::map< std::string, uint32_t > glyph_map;

	//the default font:
	static PathFont font;
};

