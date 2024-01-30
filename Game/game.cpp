#include "game.h"
#include <iostream>
#include <glm/gtc/matrix_transform.hpp>
#include "stb_image.h"
#include <vector>

#include <cstdio>
// do both exports in one function
void export_output(const char* fileName, const std::vector<std::vector<unsigned char>>& data) {
	FILE* file = fopen(fileName, "w");

	if (file != nullptr) {
		for (int i = 0; i < data.size(); i++) {
			for (int j = 0; j < data[i].size(); j++) {
				if (i < 255 || j < 255) {
					fprintf(file, "%d,", data[i][j] == 0 ? 0 : 1);
				}
				else {
					fprintf(file, "%d", data[i][j] == 0 ? 0 : 1);
				}
			}
		}

		fclose(file);
	}
}

void export_floyd(const char* fileName, const std::vector<std::vector<unsigned char>>& data) {
	FILE* file = fopen(fileName, "w");

	if (file != nullptr) {
		for (int i = 0; i < data.size(); i++) {
			for (int j = 0; j < data[i].size(); j++) {
				fprintf(file, "%d,", data[i][j] / 16);
			}
		}

		fclose(file);
	}
}


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

const int WHITE_PIXEL = 255;
const int BLACK_PIXEL = 0;

std::vector<std::vector<unsigned char>> pixelsAverage(const unsigned char* data, int square) {
	std::vector<std::vector<unsigned char>> mat(square, std::vector<unsigned char>(square));

	unsigned char pixel1;
	unsigned char pixel2;
	unsigned char pixel3;
	for (int i = 0; i < square; i++) {
		for (int j = 0; j < square; j++) {
			pixel1 = data[4 * (i * square + j)];
			pixel2 = data[4 * (i * square + j) + 1];
			pixel3 = data[4 * (i * square + j) + 2];
			unsigned char avg = (pixel1 + pixel2 + pixel3) / 3;
			mat[i][j] = avg;
		}
	}

	return mat;
}


#include <vector>

const unsigned char EDGE_VALUE = 255;
const unsigned char ORIGINAL_VALUE = 75;


//checks if the pixel is an edge pixel
bool isEdgePixel(const std::vector<std::vector<unsigned char>>& data, size_t i, size_t j) {
	if (i == 0 || j == 0 || i == data.size() - 1 || j == data[i].size() - 1)
		return true;

	const auto& rowAbove = data[i - 1];
	const auto& currentRow = data[i];
	const auto& rowBelow = data[i + 1];

	// the condition for it being an edge pixel
	return (
		rowAbove[j - 1] == EDGE_VALUE || rowAbove[j] == EDGE_VALUE || rowAbove[j + 1] == EDGE_VALUE ||
		currentRow[j - 1] == EDGE_VALUE || currentRow[j + 1] == EDGE_VALUE ||
		rowBelow[j - 1] == EDGE_VALUE || rowBelow[j] == EDGE_VALUE || rowBelow[j + 1] == EDGE_VALUE
		);
}

//the function that does the hysteresis algorithm
std::vector<std::vector<unsigned char>> hysteresisAlgo(const std::vector<std::vector<unsigned char>>& image) {
	std::vector<std::vector<unsigned char>> mat(image.size(), std::vector<unsigned char>(image[0].size()));

	for (size_t i = 0; i < image.size(); ++i) {
		for (size_t j = 0; j < image[i].size(); ++j) {
			if (image[i][j] == ORIGINAL_VALUE) {
				if (isEdgePixel(image, i, j))
					mat[i][j] = EDGE_VALUE;
				else
					mat[i][j] = 0;
			}
			else {
				mat[i][j] = image[i][j];
			}
		}
	}

	return mat;
}



//gets the angle
unsigned char getAngle(unsigned char byte) {
	//calculate the degree of the angle
	unsigned char val = static_cast<unsigned char>((byte * 180) / 3.14);
	//if the angle was negative make it between 0 and 180
	//tranformation of the angle
	if (val < 0) {
		return val + 180;
	}
	return val;
}


std::vector<std::vector<unsigned char>> nonMaxSuppression(const std::vector<std::vector<unsigned char>>& data) {
	const size_t size = data.size();
	std::vector<std::vector<unsigned char>> nonmax(size, std::vector<unsigned char>(data[0].size(), 0));

	for (size_t i = 1; i < size - 1; ++i) {
		for (size_t j = 1; j < data[i].size() - 1; ++j) {
			unsigned char val2 = 0;

			unsigned char val = getAngle(data[i][j]);

			int q = 255;
			int r = 255;

			if ((0 <= val && val < 22.5) || (157.5 <= val && val <= 180)) {
				q = data[i][j + 1];
				r = data[i][j - 1];
			}
			else if (22.5 <= val && val < 67.5) {
				q = data[i + 1][j - 1];
				r = data[i - 1][j + 1];
			}
			else if (67.5 <= val && val < 112.5) {
				q = data[i + 1][j];
				r = data[i - 1][j];
			}
			else if (112.5 <= val && val < 157.5) {
				q = data[i - 1][j - 1];
				r = data[i + 1][j + 1];
			}

			if ((data[i][j] >= q) && (data[i][j] >= r)) {
				val2 = data[i][j];
			}

			// Threshold
			val2 = (val2 >= 180) ? 255 : ((val2 <= 25) ? 0 : 75);

			nonmax[i][j] = val2;
		}
	}

	return hysteresisAlgo(nonmax);
}

//edge algorithm is the name
//change 
#include <vector>
#include <cmath>
#include <algorithm>

std::vector<std::vector<unsigned char>> edgeAlgorithm(const std::vector<std::vector<unsigned char>>& data) {
	const size_t rows = data.size();
	const size_t cols = data[0].size();

	std::vector<std::vector<unsigned char>> edge(rows, std::vector<unsigned char>(cols, 0));

	for (size_t i = 1; i < rows - 1; ++i) {
		for (size_t j = 1; j < cols - 1; ++j) {
			char gx = data[i - 1][j - 1] + 2 * data[i - 1][j] + data[i - 1][j + 1]
				- data[i + 1][j - 1] - 2 * data[i + 1][j] - data[i + 1][j + 1];
			char gy = -data[i + 1][j - 1] - 2 * data[i][j - 1] - data[i - 1][j - 1]
				+ data[i + 1][j + 1] + 2 * data[i][j + 1] + data[i - 1][j + 1];

			unsigned char val = static_cast<unsigned char>(std::abs(gx) + std::abs(gy));

			edge[i][j] = val;
		}
	}

	return nonMaxSuppression(edge);
}



std::vector<std::vector<unsigned char>> halftoneAlgorithm(const std::vector<std::vector<unsigned char>>& averageData, const unsigned char* original) {
    std::vector<std::vector<unsigned char>> halftone;

    for (int i = 0; i < 256; i++) {
        std::vector<unsigned char> line1;
        std::vector<unsigned char> line2;

        for (int j = 0; j < 256; j++) {
            int steps = averageData[i][j] / 51;

            switch (steps) {
                case 1: {
                    line1.push_back(BLACK_PIXEL);
                    line1.push_back(BLACK_PIXEL);
                    line2.push_back(WHITE_PIXEL);
                    line2.push_back(BLACK_PIXEL);
                    break;
                }
                case 2: {
                    line1.push_back(BLACK_PIXEL);
                    line1.push_back(WHITE_PIXEL);
                    line2.push_back(WHITE_PIXEL);
                    line2.push_back(BLACK_PIXEL);
                    break;
                }
                case 3: {
                    line1.push_back(BLACK_PIXEL);
                    line1.push_back(WHITE_PIXEL);
                    line2.push_back(WHITE_PIXEL);
                    line2.push_back(WHITE_PIXEL);
                    break;
                }
                case 4: {
                    line1.push_back(WHITE_PIXEL);
                    line1.push_back(WHITE_PIXEL);
                    line2.push_back(WHITE_PIXEL);
                    line2.push_back(WHITE_PIXEL);
                    break;
                }
                default:
                    line1.push_back(BLACK_PIXEL);
                    line1.push_back(BLACK_PIXEL);
                    line2.push_back(BLACK_PIXEL);
                    line2.push_back(BLACK_PIXEL);
                    break;
            }
        }

        halftone.push_back(line1);
        halftone.push_back(line2);
    }

    return halftone;
}


std::vector<std::vector<unsigned char>> floydAlgorithm(std::vector<std::vector<unsigned char>> data, unsigned char* original) {
	std::vector<std::vector<unsigned char>> floydData = pixelsAverage(original, 256);
	for (int x = 1; x < 255; x++) {
		for (int y = 1; y < 255; y++) {
			char oldpixel = data[x][y];
			char newpixel = oldpixel / 255;
			floydData[x][y] = newpixel;
			char quant_error = oldpixel - newpixel;
			floydData[x + 1][y] = data[x + 1][y] + quant_error * 7 / 16;
			floydData[x - 1][y + 1] = data[x - 1][y + 1] + quant_error * 1 / 16;
			floydData[x][y + 1] = data[x][y + 1] + quant_error * 5 / 16;
			floydData[x + 1][y + 1] = data[x + 1][y + 1] + quant_error * 3 / 16;
		}
	}
	export_floyd("../img6.txt", floydData);
	return floydData;
}

unsigned char* averageToRGB(const std::vector<std::vector<unsigned char>>& average, const unsigned char* data, int multiplier) {
	unsigned char* normal = static_cast<unsigned char*>(std::calloc(512 * 512, multiplier * multiplier));
	if (normal != nullptr) {
		int counter = 0;

		for (const auto& row : average) {
			for (size_t j = 0; j < row.size(); ++j) {
				normal[counter] = row[j];
				normal[counter + 1] = row[j];
				normal[counter + 2] = row[j];
				normal[counter + 3] = data[counter / 4 * multiplier + 3];
				counter += 4;
			}
		}
	}

	return normal;
}

void Game::Init()
{
	int width, height, numComponents;
	unsigned char* data = stbi_load("../lena256.jpg", &width, &height, &numComponents, 4);

	AddShader("../res/shaders/pickingShader");
	AddShader("../res/shaders/basicShader");
	std::vector<std::vector<unsigned char>> grayscale = pixelsAverage(data, 256);
	std::vector<std::vector<unsigned char>> halftone = halftoneAlgorithm(pixelsAverage(data, 256), data);
	std::vector<std::vector<unsigned char>> floyd = floydAlgorithm(pixelsAverage(data, 256), data);
	std::vector<std::vector<unsigned char>> edge = edgeAlgorithm(pixelsAverage(data, 256));
	export_output("../img4.txt", grayscale);
	export_output("../img5.txt", halftone);
	export_floyd("../img6.txt", floyd);
	export_output("../img7.txt", edge);
	AddTexture(256, 256, averageToRGB(grayscale, data, 1));
	AddTexture(512, 512, averageToRGB(halftone, data, 2));
	AddTexture(256, 256, averageToRGB(floyd, data, 1));
	AddTexture(256, 256, averageToRGB(edge, data, 1));

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
	int r = ((pickedShape + 1) & 0x000000FF) >> 0;
	int g = ((pickedShape + 1) & 0x0000FF00) >> 8;
	int b = ((pickedShape + 1) & 0x00FF0000) >> 16;
	s->Bind();
	s->SetUniformMat4f("MVP", MVP);
	s->SetUniformMat4f("Normal", Model);
	s->SetUniform4f("lightDirection", 0.0f, 0.0f, -1.0f, 0.0f);
	if (shaderIndx == 0)
		s->SetUniform4f("lightColor", r / 255.0f, g / 255.0f, b / 255.0f, 1.0f);
	else
		s->SetUniform4f("lightColor", 0.7f, 0.8f, 0.1f, 1.0f);
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
