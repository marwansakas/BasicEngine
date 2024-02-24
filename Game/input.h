#pragma once
#include <vector>
#include "glm/glm.hpp"
using std::vector;
using std::pair;
using namespace glm;

struct eye {
	float x;
	float y;
	float z;
	vec4 coords;
};

struct sphere {
	float x;
	float y;
	float z;
	float r;
	char type;
	float color_r;
	float color_g;
	float color_b;
	vec4 coords;
	vec4 colors;
};

struct plane {
	float x;
	float y;
	float z;
	float r;
	char type;
	float color_r;
	float color_g;
	float color_b;
	vec4 coords;
	vec4 colors;
};

struct directional_light {
	float color_r;
	float color_g;
	float color_b;
	float color_a;
	float x_direction;
	float y_direction;
	float z_direction;
	vec4 direction;
	vec4 color;
};

struct spot_light {
	float y_position;
	float z_position;
	float w_position;
	float x_direction;
	float y_direction;
	float z_direction;
	float x_position;
	float color_r;
	float color_g;
	float color_b;
	float color_a;
	vec4 direction;
	vec4 color;
	vec4 position;
};

struct ambient_light {
	float r;
	float g;
	float b;
	float a;
};

class input {
public:
	input(char* textFile);
	eye eye;
	ambient_light amb_light;
	vector<spot_light> spots;
	vector<directional_light> light_directions;
	vector<plane> planes;
	vector<sphere> balls;
	
};