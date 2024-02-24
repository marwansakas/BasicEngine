#include "input.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

/*									Parse file									*/
input::input(char* textFile) :amb_light({ 0,0,0,0 })
{
	std::ifstream infile(textFile);
	std::string line;
	vector<pair<char, glm::vec4>> objects;
	vector<glm::vec4> light_temperatues;
	vector<glm::vec4> color_temperatures;
	vector<glm::vec4> colors_intensities;
	vector<glm::vec4> positions;
	while (std::getline(infile, line))
	{
		std::istringstream in_s(line);
		char mode;
		float x, y, z, w;
		if (!(in_s >> mode >> x >> y >> z >> w)) { break; }

		glm::vec4 vec(x, y, z, w);
		if (mode == 'e') {
			eye = { x, y, z, vec };
		}
		else if (mode == 'a') {
			amb_light = { x, y, z, w };
		}
		else if (mode == 'd') {
			light_temperatues.push_back(vec);
		}
		else if (mode == 'p') {
			positions.push_back(vec);
		}
		else if (mode == 'i') {
			colors_intensities.push_back(vec);
		}
		else if (mode == 'o' || mode == 'r' || mode == 't') {
			objects.push_back({ mode,vec });
		}
		else if (mode == 'c') {
			color_temperatures.push_back(vec);
		}
		
	}

	int count = 0;
	int temp_size = light_temperatues.size();
	int num_objects = objects.size();
	for (unsigned int i = 0; i < temp_size; i++) {
		if (light_temperatues[i].w == 0) {
			light_directions.push_back(
				{ light_temperatues[i].x, light_temperatues[i].y, light_temperatues[i].z, 
				colors_intensities[i].x, colors_intensities[i].y, colors_intensities[i].z ,colors_intensities[i].w,
				light_temperatues[i], colors_intensities[i]});
		}
		else {
			spots.push_back(
				{ light_temperatues[i].x, light_temperatues[i].y, light_temperatues[i].z, 
				positions[count].x, positions[count].y, positions[count].z,positions[count].w,
				colors_intensities[i].x, colors_intensities[i].y, colors_intensities[i].z,colors_intensities[i].w,
				light_temperatues[i], colors_intensities[i], positions[count]});
			count++;
		}
	}
	for (unsigned int i = 0; i < num_objects; i++) {
		if (objects[i].second.w > 0) 
			balls.push_back({ objects[i].second.x, objects[i].second.y,objects[i].second.z,objects[i].second.w
				,objects[i].first, color_temperatures[i].x, color_temperatures[i].y,color_temperatures[i].z, objects[i].second, color_temperatures[i]});
		else 
			planes.push_back({ objects[i].second.x, objects[i].second.y,objects[i].second.z,objects[i].second.w
				,objects[i].first, color_temperatures[i].x, color_temperatures[i].y,color_temperatures[i].z, objects[i].second, color_temperatures[i] });
	}
}



