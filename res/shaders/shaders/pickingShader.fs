#version 330

in vec3 color0;
in vec3 normal0;

uniform sampler2D sampler;
uniform vec4 lightDirection;
uniform vec4 lightColor;

out vec4 fragColor;

void main()
{
	fragColor = lightColor;
}
