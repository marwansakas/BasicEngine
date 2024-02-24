#version 330 core

in vec3 position;
in vec3 color;
in vec3 normal;
in vec2 texCoords;

uniform mat4 MVP;
uniform mat4 Normal;

//out vec3 color0;

void main()
{
	//color0 = color;
	gl_Position = MVP *Normal* vec4(position, 1.0);
}
