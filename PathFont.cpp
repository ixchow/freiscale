// Based on code from Chesskoban (c) 2017-2019 Jim McCann;
// this adapted-for-15-466 code is released into the public domain.

#include "PathFont.hpp"

#include <iostream>

PathFont::PathFont(uint32_t glyphs_,
	const float *glyph_widths_,
	const uint32_t *glyph_char_starts_, const uint8_t *chars_,
	const uint32_t *glyph_coord_starts_, const float *coords_
	) : glyphs(glyphs_),
		glyph_widths(glyph_widths_),
		glyph_char_starts(glyph_char_starts_), chars(chars_),
		glyph_coord_starts(glyph_coord_starts_), coords(coords_) {

	for (uint32_t i = 0; i < glyphs; ++i) {
		std::string str(reinterpret_cast< const char * >(chars + glyph_char_starts[i]), reinterpret_cast< const char * >(chars + glyph_char_starts[i+1]));
		auto res = glyph_map.insert(std::make_pair(str, i));
		if (!res.second) {
			std::cerr << "WARNING: ignoring duplicate glyph for '" << str << "'." << std::endl;
		}
	}
}
