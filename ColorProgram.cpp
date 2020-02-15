#include "ColorProgram.hpp"

kit::Load< ColorProgram > color_program(kit::LoadTagInit);

ColorProgram::ColorProgram() : GLProgram(
		//vertex shader:
		"#version 330\n"
		"uniform mat4 OBJECT_TO_CLIP;\n"
		"in vec4 Position;\n"
		"in vec4 Color;\n"
		"out vec4 color;\n"
		"void main() {\n"
		"	gl_Position = OBJECT_TO_CLIP * Position;\n"
		"	color = Color;\n"
		"}\n"
	,
		//fragment shader:
		"#version 330\n"
		"in vec4 color;\n"
		"out vec4 fragColor;\n"
		"void main() {\n"
		"	fragColor = color;\n"
		"}\n"
	) {

	//look up the locations of vertex attributes:
	Position_vec4 = getAttribLocation("Position");
	Color_vec4 = getAttribLocation("Color");

	//look up the locations of uniforms:
	OBJECT_TO_CLIP_mat4 = getUniformLocation("OBJECT_TO_CLIP");
}
