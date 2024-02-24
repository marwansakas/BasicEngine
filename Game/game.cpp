#include "game.h"
#include <iostream>
#include <glm/gtc/matrix_transform.hpp>
#include "../Engine3D/stb_image.h"
#include <stdio.h>
#include <iostream>

static void printMat(const glm::mat4 mat)
{
	std::cout << " matrix:" << std::endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			std::cout << mat[j][i] << " ";
		std::cout << std::endl;
	}
}

Game::Game() : Scene()
{
}

Game::Game(float angle, float relationWH, float near1, float far1) : Scene(angle, relationWH, near1, far1)
{
}

void Game::Init()
{
	AddShader("../res/shaders/pickingShader");
	AddShader("../res/shaders/basicShader");

	AddShape(Plane, -1, TRIANGLES);

	pickedShape = 0;

	SetShapeTex(0, 0);
	MoveCamera(0, zTranslate, 10);
	pickedShape = -1;

	//ReadPixel(); //uncomment when you are reading from the z-buffer
}

void Game::Update(const glm::mat4& MVP, const glm::mat4& Model, const int  shaderIndx)
{
	Shader* s = shaders[shaderIndx];
	s->Bind();
	s->SetUniformMat4f("MVP", MVP);
	s->SetUniformMat4f("Normal", Model);
	s->Unbind();
}

void Game::WhenRotate()
{
}

void Game::WhenTranslate()
{
}

void Game::Motion()
{
	if (isActive)
	{

	}
}

Game::~Game(void)
{
}
