#pragma once

/*
 * Helper to throw geometry at OpenGL.
 * not quite immediate mode, but close to it.
 *
 */

#include <kit/GLBuffer.hpp>
#include <glm/glm.hpp>

#include <string>
#include <vector>

struct DrawStuff {
	typedef GLAttribBuffer< glm::vec3, glm::u8vec4 >::Vertex Pos3f_Col4ub;
	typedef GLAttribBuffer< glm::vec2, glm::u8vec4 >::Vertex Pos2f_Col4ub;

	static void draw(glm::mat4 const &world_to_clip, GLenum mode, std::vector< Pos3f_Col4ub > const &verts);
	static void draw(glm::mat4 const &world_to_clip, GLenum mode, std::vector< Pos2f_Col4ub > const &verts);
};
