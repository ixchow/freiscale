#include "DrawStuff.hpp"
#include "ColorProgram.hpp"

#include <kit/GLVertexArray.hpp>
#include <kit/gl_errors.hpp>

static GLAttribBuffer< glm::vec3, glm::u8vec4 > *pos3_col4ub_buffer = nullptr;
static kit::Load< GLVertexArray > pos3_col4ub_buffer_for_color_program(kit::LoadTagDefault, []() -> GLVertexArray * {
	pos3_col4ub_buffer = new GLAttribBuffer< glm::vec3, glm::u8vec4 >();
	return new GLVertexArray(GLVertexArray::make_binding(color_program->program, {
		{color_program->getAttribLocation("Position"), (*pos3_col4ub_buffer)[0]},
		{color_program->getAttribLocation("Color"), (*pos3_col4ub_buffer)[1]}
	}));
});

static GLAttribBuffer< glm::vec2, glm::u8vec4 > *pos2_col4ub_buffer = nullptr;
static kit::Load< GLVertexArray > pos2_col4ub_buffer_for_color_program(kit::LoadTagDefault, []() -> GLVertexArray * {
	pos2_col4ub_buffer = new GLAttribBuffer< glm::vec2, glm::u8vec4 >();
	return new GLVertexArray(GLVertexArray::make_binding(color_program->program, {
		{color_program->getAttribLocation("Position"), (*pos2_col4ub_buffer)[0]},
		{color_program->getAttribLocation("Color"), (*pos2_col4ub_buffer)[1]}
	}));
});

void DrawStuff::draw(glm::mat4 const &world_to_clip, GLenum mode, std::vector< DrawStuff::Pos3f_Col4ub > const &verts) {
	if (verts.empty()) return;

	pos3_col4ub_buffer->set(verts, GL_STREAM_DRAW);

	glUseProgram(color_program->program);
	glUniformMatrix4fv(color_program->OBJECT_TO_CLIP_mat4, 1, GL_FALSE, glm::value_ptr(world_to_clip));
	glBindVertexArray(pos3_col4ub_buffer_for_color_program->array);
	glDrawArrays(mode, 0, GLsizei(verts.size()));
	glBindVertexArray(0);
	glUseProgram(0);
}


void DrawStuff::draw(glm::mat4 const &world_to_clip, GLenum mode, std::vector< DrawStuff::Pos2f_Col4ub > const &verts) {
	if (verts.empty()) return;

	pos2_col4ub_buffer->set(verts, GL_STREAM_DRAW);

	glUseProgram(color_program->program);
	glUniformMatrix4fv(color_program->OBJECT_TO_CLIP_mat4, 1, GL_FALSE, glm::value_ptr(world_to_clip));
	glBindVertexArray(pos2_col4ub_buffer_for_color_program->array);
	glDrawArrays(mode, 0, GLsizei(verts.size()));
	glBindVertexArray(0);
	glUseProgram(0);
}
