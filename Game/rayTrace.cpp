#include "rayTrace.h"
rayTrace::rayTrace(char* file): in(input(file)) {}


/*									Help Functions									*/
float distance(glm::vec3 vector1, glm::vec4 vector2) {
	float x = (vector1.x - vector2.x);
	float y = (vector1.y - vector2.y);
	float z = (vector1.z - vector2.z);
	
	return sqrt(x * x + y * y + z * z);
}



glm::vec3 handle_spotlight(glm::vec3 hit, glm::vec3 position, glm::vec3 direction_ray_normal, float cos_value) {
	
	glm::vec3 ray_spotlight_normal = glm::normalize(hit - position);
	
	float light_cos_vlue = glm::dot(ray_spotlight_normal, direction_ray_normal);
	
	if (light_cos_vlue < cos_value) { 
		return glm::vec3(0.0f); 
	} else {
		return ray_spotlight_normal;
	}
}



glm::vec3 rayTrace::getVecNormal(bool is_plane, glm::vec3 hit_coardinates, int hit_index) {
	if (!is_plane)
		return glm::normalize(hit_coardinates - glm::vec3(in.balls[hit_index].coords));
	if (in.planes[hit_index].coords != glm::vec4(0.0f))
		return glm::normalize(glm::vec3(in.planes[hit_index].coords));
	else return glm::vec3(in.planes[hit_index].coords);
}



/*								   Render picture									*/
unsigned char* rayTrace::getPicture(int DISPLAY_WIDTH, int DISPLAY_HEIGHT)
{
	uint32_t* data = new uint32_t[DISPLAY_WIDTH * DISPLAY_HEIGHT];
	for (int y = 0; y < DISPLAY_HEIGHT; y++) {
		for (int x = 0; x < DISPLAY_WIDTH; x++) {
			glm::vec2 coardinate = { ((float)x) / ((float) DISPLAY_WIDTH), ((float) y) / ((float) DISPLAY_HEIGHT) };
			coardinate = coardinate * 2.0f - 1.0f;
			data[x + y * DISPLAY_WIDTH] = perPixel(glm::vec3(coardinate.x, -coardinate.y, 0.0f));
		}
	}
	return (unsigned char*) data;
}



uint32_t rayTrace::perPixel(glm::vec3 ray_direction) {
	int hit_index = -1;
	bool is_plane = false;
	bool dark = false;

	Ray ray1 = { glm::vec3(in.eye.x, in.eye.y, in.eye.z) , ray_direction - glm::vec3(in.eye.x, in.eye.y, in.eye.z) };
	Ray ray2 = { glm::vec3(in.eye.x, in.eye.y, in.eye.z) , ray_direction};

	glm::vec3 ray_position = glm::vec3(in.eye.x, in.eye.y, in.eye.z);
	glm::vec3 hit_coordinate = getHitCoords(ray1, hit_index, is_plane, dark);
	glm::vec4 color_of_pixel = getFinalColor(ray1, hit_coordinate, hit_index, is_plane, dark, 0);

	unsigned int r = color_of_pixel.x * 255.0;
	unsigned int g = color_of_pixel.y * 255.0;
	unsigned int b = color_of_pixel.z * 255.0;
	unsigned int a = color_of_pixel.w * 255.0;

	g = g << 8;
	b = b << 16;
	a = a << 24;

	return r + g + b + a;
}



glm::vec3 rayTrace::getHitCoords(Ray ray, int& hit_index, bool& plane, bool& dark) {
	int hit_prev = hit_index;
	bool prev_plane = plane;
	float hit_distance = -1;
	glm::vec3 hit_coordinate(0.0f);

	if (ray.ray_direction != glm::vec3(0.0f))
		ray.ray_direction = glm::normalize(ray.ray_direction);

	for (unsigned int i = 0; i < in.planes.size(); i++) {
		if (!(hit_prev != -1 && in.planes[i].type == 't')) {
			if (!(prev_plane && hit_prev == i)) {
				glm::vec3 val = calculatePlane(ray, in.planes[i].coords);
				if (val != glm::vec3(0.0f))
					if (hit_distance == -1 || distance(val, in.eye.coords) < hit_distance) {
						hit_coordinate = val;
						hit_distance = distance(val, in.eye.coords);
						plane = true;
						hit_index = i;
						if (in.planes[i].type == 'o')
							if (abs(int(2.0 * (val.x + 10))) % 2 != abs(int(2.0 * (val.y + 10))) % 2) dark = true;
							else dark = false;
					}
			}
		}
	}
	for (unsigned int i = 0; i < in.balls.size(); i++) {
		if (!(!prev_plane && hit_prev == i)) {
			if (!(hit_prev != -1 && in.balls[i].type == 't')) {
				glm::vec3 val = calculateSphere(in.balls[i].coords, ray);
				if (val != glm::vec3(0.0f))
					if (hit_distance == -1 || distance(val, in.eye.coords) < hit_distance) {
						hit_coordinate = val;
						hit_distance = distance(val, in.eye.coords);
						plane = false;
						dark = false;
						hit_index = i;
					}
			}
		}
	}
	return hit_coordinate;
}

glm::vec4 rayTrace::getFinalColor(Ray ray, glm::vec3 hit_coordinate, int hit_index, bool is_plane, bool dark, int depth) {
	bool is_plane2 = is_plane, isDark2 = dark;
	int hitIndex2 = hit_index;
	glm::vec3 phong_model = glm::vec3(0, 0, 0);
	if (is_plane) {
		plane plane_value = in.planes[hit_index];
		if (plane_value.type == 'o') {
			phong_model = handle_object(hit_coordinate, phong_model, glm::vec3(plane_value.colors), plane_value.colors.w, hitIndex2, is_plane2, isDark2);
		}
		else if (plane_value.type == 'r') {
			phong_model = handle_reflective(hit_coordinate, ray.ray_direction, phong_model, glm::vec3(plane_value.colors), plane_value.colors.w, hitIndex2, is_plane2, isDark2, depth);
		}
		else if (plane_value.type == 't') {
			phong_model = handle_transparent(hit_coordinate, ray.ray_direction, phong_model, glm::vec3(plane_value.colors), plane_value.colors.w, hitIndex2, is_plane2, isDark2, depth);
		}
		if (dark)
			return glm::vec4(phong_model.r, phong_model.g, phong_model.b, 0.0) * 0.5f;
		else return glm::vec4(phong_model.r, phong_model.g, phong_model.b, 0.0);
	}
	else if (hit_index != -1) {
		sphere sphere = in.balls[hit_index];
		if (sphere.type == 'o') {
			phong_model = handle_object(hit_coordinate, phong_model, glm::vec3(sphere.colors), sphere.colors.w, hitIndex2, is_plane2, isDark2);
		}
		else if (sphere.type == 'r') {
			phong_model = handle_reflective(hit_coordinate, ray.ray_direction, phong_model, glm::vec3(sphere.colors), sphere.colors.w, hitIndex2, is_plane2, isDark2, depth);
		}
		else if (sphere.type == 't') {
			phong_model = handle_transparent(hit_coordinate, ray.ray_direction, phong_model, glm::vec3(sphere.colors), sphere.colors.w, hitIndex2, is_plane2, isDark2, depth);
		}
		if (dark)
			return glm::vec4(phong_model.r, phong_model.g, phong_model.b, 0.0) * 0.5f;
		else 
			return glm::vec4(phong_model.r, phong_model.g, phong_model.b, 0.0);
	}
	return glm::vec4(0.0f);
}



/*									Intersections									*/
glm::vec3 rayTrace::calculatePlane(Ray ray, glm::vec4 plane) {
	glm::vec3 plane_normal(plane.x, plane.y, plane.z);
	float d = plane.w;
	if (glm::dot(ray.ray_direction, plane_normal) != 0) {
		float ans = -(glm::dot(plane_normal, ray.ray_position) + d) / glm::dot(plane_normal, ray.ray_direction);
		if (ans >= 0) {
			glm::vec3 res = ray.ray_position+ ray.ray_direction * ans;
			return res;
		}
	}
	return glm::vec3(0.0f);
}



glm::vec3 rayTrace::calculateSphere(glm::vec4 sphere, Ray ray) {
	float discriminant;
	glm::vec3 center = ray.ray_position - glm::vec3(sphere.x, sphere.y, sphere.z);

	float a = 1;
	float b = 2.0f * glm::dot(ray.ray_direction, center);
	float c = glm::dot(center, center) - sphere.w * sphere.w;

	discriminant = b * b - 4.0f * a * c;
	
	float delta;
	float ans1;
	float ans2;
	if (discriminant >= 0.0f) {
		delta = sqrt(discriminant);
		ans1 = (-b + delta) / (2 * a);
		ans2 = (-b - delta) / (2 * a);
		discriminant = glm::min(ans1, ans2);
		return ray.ray_position+ ray.ray_direction * discriminant;
	}
	return glm::vec3(0.0f);
}




/*									Handle Sphere/Plane									*/
glm::vec3 rayTrace::handle_object(glm::vec3 hit_coordinate, glm::vec3 phong_model, glm::vec3 color, float shiny_value, int hit_index, bool is_plane, bool dark) {
	phong_model = color * glm::vec3(in.amb_light.r, in.amb_light.g, in.amb_light.b);
	for (int i = 0; i < in.spots.size(); i++) {
		glm::vec3 color_for_diffuse = getDiffuseColor(hit_coordinate, glm::vec3(in.spots[i].direction), glm::vec3(in.spots[i].position), glm::vec3(in.spots[i].color), color, false, hit_index, is_plane, dark, in.spots[i].position.w);
		color_for_diffuse = glm::max(color_for_diffuse, glm::vec3(0, 0, 0));
		glm::vec3 color_for_specular = handle_specular(hit_coordinate, glm::vec3(in.spots[i].direction), glm::vec3(in.spots[i].position), glm::vec3(in.spots[i].color), color, false, hit_index, is_plane, dark, in.spots[i].position.w, shiny_value);
		color_for_specular = glm::max(color_for_specular, glm::vec3(0, 0, 0));
		if (hasShadow(hit_coordinate, glm::vec3(in.spots[i].direction), glm::vec3(in.spots[i].position), glm::vec3(in.spots[i].color), color, false, hit_index, is_plane, dark, in.spots[i].position.w, shiny_value))
			phong_model += color_for_diffuse + color_for_specular;
	}
	for (int i = 0; i < in.light_directions.size(); i++) {
		
		glm::vec3 color_for_diffuse = getDiffuseColor(hit_coordinate, glm::vec3(in.light_directions[i].direction), glm::vec3(glm::vec4(0.0f)), glm::vec3(in.light_directions[i].color), color, true, hit_index, is_plane, dark, glm::vec4(0.0f).w);
		color_for_diffuse = glm::max(color_for_diffuse, glm::vec3(0, 0, 0));
		
		glm::vec3 color_for_specular = handle_specular(hit_coordinate, glm::vec3(in.light_directions[i].direction), glm::vec3(glm::vec4(0.0f)), glm::vec3(in.light_directions[i].color), color, true, hit_index, is_plane, dark, glm::vec4(0.0f).w, shiny_value);
		color_for_specular = glm::max(color_for_specular, glm::vec3(0, 0, 0));
		
		if (hasShadow(hit_coordinate, glm::vec3(in.light_directions[i].direction), glm::vec3(glm::vec4(0.0f)), glm::vec3(in.light_directions[i].color), color, true, hit_index, is_plane, dark, glm::vec4(0.0f).w, shiny_value))
			phong_model += color_for_diffuse + color_for_specular;
	}
	return phong_model = glm::min(phong_model, glm::vec3(1.0f, 1.0f, 1.0f));
}





glm::vec3 rayTrace::handle_reflective(glm::vec3 hit_coordinate, glm::vec3 ray_direction, glm::vec3 phongModel, glm::vec3 color, float shiny_value, int hit_index, bool is_plane, bool dark, int depth) {
	if (depth == 5) {
		return glm::vec3(0.f, 0.f, 0.f);
	}
	Ray ray = {hit_coordinate, ray_direction};
	
	glm::vec3 normal_of_object = getVecNormal(is_plane, hit_coordinate, hit_index);
	glm::vec3 reflection_direction = ray_direction - 2.0f * normal_of_object * (glm::dot(ray_direction, normal_of_object));
	
	Ray ray1 = { hit_coordinate, reflection_direction};
	
	glm::vec3 hit_coordinate2 = getHitCoords(ray1, hit_index, is_plane, dark);
	
	if (hit_index == -1) {
		return glm::vec3(0.f, 0.f, 0.f);
	}
	
	glm::vec4 color_of_reflection = getFinalColor(ray, hit_coordinate2, hit_index, is_plane, dark, depth + 1);
	phongModel = glm::vec3(color_of_reflection.r, color_of_reflection.g, color_of_reflection.b);
	
	return glm::min(phongModel, glm::vec3(1.0, 1.0, 1.0));
}





glm::vec3 rayTrace::handle_transparent(glm::vec3 hit_coordinate, glm::vec3 ray_direction, glm::vec3 phong_model, glm::vec3 color, float shiny_value, int hit_index, bool is_plane, bool dark, int depth) {
	if (depth == 5) {
		return glm::vec3(0.f, 0.f, 0.f);
	}
	
	glm::vec4 color_of_transparency = glm::vec4(0.f, 0.f, 0.f, 0.f);
	Ray ray = {hit_coordinate, ray_direction};
	
	if (is_plane) {
		glm::vec3 hit_coordinate2 = getHitCoords(ray, hit_index, is_plane, dark);
		if (hit_index == -1) {
			return glm::vec3(0.f, 0.f, 0.f);
		}
		color_of_transparency = getFinalColor(ray, hit_coordinate2, hit_index, is_plane, dark, depth + 1);
	} else {
		glm::vec3 normalForObj = getVecNormal(is_plane, hit_coordinate, hit_index);
		
		float pi = 3.14159265f;
		float cos1 = glm::dot(normalForObj, -ray_direction);
		float theta1 = acos(cos1) * (180.0f / pi);
		float sin1 = sin(theta1 * (pi / 180.0f));
		
		float snellFraction = (1.0f / 1.5f);
		
		float sin2 = snellFraction * sin1;
		float theta2 = asin(sin2) * (180.0f / pi);
		float cos2 = cos(theta2 * (pi / 180.0f));

		glm::vec3 ray_directional_to_in = (snellFraction * cos1 - cos2) * normalForObj - snellFraction * (-ray_direction);
		ray_directional_to_in = glm::normalize(ray_directional_to_in);
		glm::vec3 hit_coordinate2 = getHitCoords(ray, hit_index, is_plane, dark);
		
		Ray ray2 = { hit_coordinate, ray_directional_to_in };

		if (hit_index != -1) {
			color_of_transparency = getFinalColor(ray2, hit_coordinate2, hit_index, is_plane, dark, depth + 1);
		} else {
			glm::vec3 theSecondHitPoint = calculateSphere(in.balls[hit_index].coords, ray2);

			cos1 = glm::dot(-glm::normalize(theSecondHitPoint), -ray_directional_to_in);
			theta1 = acos(cos1) * (180.0f / pi);
			snellFraction = (1.5f / 1.0f);
			sin1 = sin(theta1 * (pi / 180.0f));
			sin2 = snellFraction * sin1;
			theta2 = asin(sin2) * (180.0f / pi);
			cos2 = cos(theta2 * (pi / 180.0f));

			glm::vec3 rayDirectionToOut = (snellFraction * cos1 - cos2) * -normalForObj - snellFraction * (-ray_directional_to_in);
			rayDirectionToOut = glm::normalize(rayDirectionToOut);
			
			Ray ray3 = { theSecondHitPoint, rayDirectionToOut };
		
			glm::vec3 hitCoord2 = getHitCoords(ray3, hit_index, is_plane, dark);

			if (hit_index == -1) {
				return glm::vec3(0.f, 0.f, 0.f);
			}
			
			color_of_transparency = getFinalColor(ray3, hitCoord2, hit_index, is_plane, dark, depth + 1);
		}
	}

	phong_model = glm::vec3(color_of_transparency.r, color_of_transparency.g, color_of_transparency.b);
	phong_model = glm::min(phong_model, glm::vec3(1.0, 1.0, 1.0));
	
	return phong_model;
}




/*									Handle Lights									*/
glm::vec3 rayTrace::getDiffuseColor(glm::vec3 hit, glm::vec3 direction, glm::vec3 position, glm::vec3 color_of_light, glm::vec3 color, bool is_directional, int hit_index, bool is_plane, bool dark, float cos_value) {
	float object_factor;
	if (is_plane) {
		object_factor = -1.0f;
	}
	else {
		object_factor = 1.0f;
	}
	glm::vec3 ray_direction_normal = object_factor * glm::normalize(direction);
	if (!is_directional) {
		ray_direction_normal = handle_spotlight(hit, position, object_factor * ray_direction_normal, cos_value) * object_factor;
		if (ray_direction_normal == glm::vec3(0.0f))
			return ray_direction_normal;
	}

	glm::vec3 normal_of_object = getVecNormal(is_plane, hit, hit_index);
	float cosvalue_of_hit = glm::dot(normal_of_object, -ray_direction_normal);

	glm::vec3 color_of_diffuse = color * cosvalue_of_hit * color_of_light;
	return color_of_diffuse;
}






glm::vec3 rayTrace::handle_specular(glm::vec3 hit, glm::vec3 direction, glm::vec3 position, glm::vec3 light_color, glm::vec3 color, bool is_directional, int hit_index, bool is_plane, bool dark, float cos_value, float shiny_value) {
	glm::vec3 ray_direction_normal = direction;
	if (direction != glm::vec3(0.0f))
		ray_direction_normal = glm::normalize(direction);

	if (!is_directional) {
		ray_direction_normal = handle_spotlight(hit, position, ray_direction_normal, cos_value);
		if (ray_direction_normal == glm::vec3(0.0f))
			return ray_direction_normal;
	}

	glm::vec3 normalized_Object = getVecNormal(is_plane, hit, hit_index);
	glm::vec3 reflected_lightRay = ray_direction_normal - 2.0f * normalized_Object * glm::dot(ray_direction_normal, normalized_Object);
	glm::vec3 ray_to_veye = glm::normalize(position - hit);

	float cosValue_of_hit = glm::dot(ray_to_veye, reflected_lightRay);
	
	cosValue_of_hit = glm::max(0.0f, cosValue_of_hit);
	cosValue_of_hit = pow(cosValue_of_hit, shiny_value);

	glm::vec3 color_rizz = 0.7f * cosValue_of_hit * light_color;
	
	return color_rizz;
}




bool rayTrace::hasShadow(glm::vec3 hit, glm::vec3 direction, glm::vec3 position, glm::vec3 lightColor, glm::vec3 color, bool is_direction, int hit_index, bool is_plane, bool dark, float cos_value, float shiny_value) {
	glm::vec3 ray_direction_normal = direction;
	
	if (direction != glm::vec3(0.0f))
		ray_direction_normal = glm::normalize(direction);

	float min = INFINITY;

	if (!is_direction) {
		ray_direction_normal = handle_spotlight(hit, position, ray_direction_normal, cos_value);
		if (ray_direction_normal == glm::vec3(0.0f)) {
			return false;
		}
		else 
			min = -(glm::dot(hit, position)) / std::abs(glm::dot(-ray_direction_normal, position));
	}

	int hit_index2 = hit_index;
	bool is_plane2 = is_plane;
	Ray ray = { hit, -ray_direction_normal };
	
	glm::vec3 hit_coardenate = getHitCoords(ray, hit_index2, is_plane2, dark);
	
	if (hit_coardenate != glm::vec3(0.0f) && hit_index2 != -1 && distance(hit_coardenate, glm::vec4(hit.x, hit.y, hit.z, 0.0f)) < min) {
		return false;
	}
	return true;
}