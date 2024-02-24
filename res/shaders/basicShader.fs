#version 330

in vec2 texCoord0;
in vec3 normal0;
in vec3 color0;
in vec3 position0;

uniform vec4 lightColor;
uniform sampler2D sampler;
uniform vec4 lightDirection;

out vec4 fragColor;
void main()
{
	fragColor = texture(sampler, texCoord0)* vec4(color0,1.0); //you must have gl_FragColor
}