#pragma once
#include "input.h"
#include <vector>
#include <string>
#include "glm/glm.hpp"
#include <windows.h>
#include <glm/gtx/intersect.hpp>
using std::vector;
using std::pair;

struct Ray {
	glm::vec3 ray_position;
	glm::vec3 ray_direction;
};

class rayTrace
{
public:
	rayTrace(char* file);
	unsigned char* getPicture(int DISPLAY_WIDTH, int DISPLAY_HEIGHT);
private:
	uint32_t rayTrace::perPixel(glm::vec3 ray_direction);
	glm::vec3 getHitCoords(Ray ray, int& hit_index, bool& plane, bool& dark);
	glm::vec4 getFinalColor(Ray ray, glm::vec3 hit_coordinate, int hit_index, bool is_plane, bool dark, int depth);
	glm::vec3 rayTrace::getDiffuseColor(glm::vec3 hit, glm::vec3 direction, glm::vec3 position, glm::vec3 color_of_light, glm::vec3 color, bool is_directional, int hit_index, bool is_plane, bool dark, float cos_value);
	glm::vec3 rayTrace::handle_specular(glm::vec3 hit, glm::vec3 direction, glm::vec3 position, glm::vec3 light_color, glm::vec3 color, bool is_directional, int hit_index, bool is_plane, bool dark, float cos_value, float shiny_value);
	bool rayTrace::hasShadow(glm::vec3 hit, glm::vec3 direction, glm::vec3 position, glm::vec3 color_of_light, glm::vec3 color, bool is_direction, int hit_index, bool is_plane, bool dark, float cos_value, float shiny_value);
	glm::vec3 calculatePlane(Ray ray, glm::vec4 plane);
	glm::vec3 calculateSphere(glm::vec4 sphere, Ray ray);
	glm::vec3 rayTrace::handle_object(glm::vec3 hit_coordinate, glm::vec3 phong_model, glm::vec3 color, float shiny_value, int hit_index, bool is_plane, bool dark);
	glm::vec3 rayTrace::handle_reflective(glm::vec3 hit_coordinate, glm::vec3 ray_direction, glm::vec3 phongModel, glm::vec3 color, float shiny_value, int hit_index, bool is_plane, bool dark, int depth);
	glm::vec3 rayTrace::handle_transparent(glm::vec3 hit_coordinate, glm::vec3 ray_direction, glm::vec3 phong_model, glm::vec3 color, float shiny_value, int hit_index, bool is_plane, bool dark, int depth);
	glm::vec3 rayTrace::getVecNormal(bool is_plane, glm::vec3 hit_coardinates, int hit_index);

	input in;
};

