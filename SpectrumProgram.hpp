#pragma once

#include <kit/Load.hpp>
#include <kit/GLProgram.hpp>

//Shader program that draws transformed, colored vertices:
struct SpectrumProgram : GLProgram {
	SpectrumProgram();

	//Attribute (per-vertex variable) locations:
	GLuint Position_vec4 = -1U;
	GLuint TexCoord_vec2 = -1U;

	//Uniform (per-invocation variable) locations:
	GLuint OBJECT_TO_CLIP_mat4 = -1U;

	//Textures:
	// 0 -- spectrum
};

extern kit::Load< SpectrumProgram > spectrum_program;
