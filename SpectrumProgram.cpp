#include "SpectrumProgram.hpp"

kit::Load< SpectrumProgram > spectrum_program(kit::LoadTagInit);

SpectrumProgram::SpectrumProgram() : GLProgram(
		//vertex shader:
		"#version 330\n"
		"uniform mat4 OBJECT_TO_CLIP;\n"
		"in vec4 Position;\n"
		"in vec2 TexCoord;\n"
		"out vec2 texCoord;\n"
		"void main() {\n"
		"	gl_Position = OBJECT_TO_CLIP * Position;\n"
		"	texCoord = TexCoord;\n"
		"}\n"
	,
		//fragment shader:
		"#version 330\n"
		"uniform sampler2D tex;\n"
		"in vec2 texCoord;\n"
		"out vec4 fragColor;\n"
		"void main() {\n"
		"	float p = texture(tex, texCoord).r;\n"
		"	fragColor = vec4(p * 0.25, p * 0.5, p * 1.0, 1.0);\n"
		"}\n"
	) {

	//look up the locations of vertex attributes:
	Position_vec4 = getAttribLocation("Position");
	TexCoord_vec2 = getAttribLocation("TexCoord");

	//look up the locations of uniforms:
	OBJECT_TO_CLIP_mat4 = getUniformLocation("OBJECT_TO_CLIP");

	//set 'tex' to texture 0:
	glUseProgram(program);
	glUniform1i(getUniformLocation("tex"), 0);
	glUseProgram(0);
}
