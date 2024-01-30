#include <glfw/deps/glad/gl.h>
#include "InputManager.h"
// #include "../DisplayGLFW/display.h"
#include "game.h"
#include "../res/includes/glm/glm.hpp"
#include "stb_image.h";



int main(int argc, char* argv[])
{
	fprintf(stderr, "test\n");
	const int DISPLAY_WIDTH = 512;
	const int DISPLAY_HEIGHT = 512;
	const float CAMERA_ANGLE = 0.0f;
	const float NEAR = 1.0f;
	const float FAR = 100.0f;

	Game* scn = new Game(CAMERA_ANGLE, (float)DISPLAY_WIDTH / DISPLAY_HEIGHT, NEAR, FAR);

	Display display(DISPLAY_WIDTH, DISPLAY_HEIGHT, "OpenGL");

	Init(display);

	scn->Init();

	display.SetScene(scn);

	while (!display.CloseWindow())
	{

		// GrayScale
		scn->SetShapeTex(0, 0);
		glViewport(0, DISPLAY_HEIGHT / 2, DISPLAY_WIDTH / 2, DISPLAY_HEIGHT / 2);
		scn->Draw(1, 0, scn->BACK, true, false);

		// half tone
		scn->SetShapeTex(0, 1);
		glViewport(0, 0, DISPLAY_WIDTH / 2, DISPLAY_HEIGHT / 2);
		scn->Draw(1, 0, scn->BACK, false, false);

		// floyd
		scn->SetShapeTex(0, 2);
		glViewport(DISPLAY_WIDTH / 2, 0, DISPLAY_WIDTH / 2, DISPLAY_HEIGHT / 2);
		scn->Draw(1, 0, scn->BACK, false, false);

		// edges
		scn->SetShapeTex(0, 3);
		glViewport(DISPLAY_WIDTH / 2, DISPLAY_HEIGHT / 2, DISPLAY_WIDTH / 2, DISPLAY_HEIGHT / 2);
		scn->Draw(1, 0, scn->BACK, false, false);

		scn->Motion();
		display.SwapBuffers();
		display.PollEvents();

	}
	delete scn;
	return 0;
}
