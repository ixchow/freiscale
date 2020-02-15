#pragma once

#include <kit/Load.hpp>
#include <kit/GLProgram.hpp>

//Shader program that draws transformed, colored vertices:
struct ColorProgram : GLProgram {
	ColorProgram();

	//Attribute (per-vertex variable) locations:
	GLuint Position_vec4 = -1U;
	GLuint Color_vec4 = -1U;

	//Uniform (per-invocation variable) locations:
	GLuint OBJECT_TO_CLIP_mat4 = -1U;

	//Textures:
	// none
};

extern kit::Load< ColorProgram > color_program;
