//#define GLEW_STATIC
#include <iostream>
#include "../res/includes/glad/include/glad/glad.h"
#include "display.h"

Display::Display(int width, int height, const std::string& title)
{
	/* Initialize the library */
    if (!glfwInit())
        error =  -1;

    glfwWindowHint(GLFW_SAMPLES, 8);
   
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);

	#ifdef __APPLE__
	    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	#endif

	m_window = glfwCreateWindow(width, height, title.c_str(), NULL, NULL);
	if(!m_window)
	{
		glfwTerminate();
		error = -1;
	}
	glfwMakeContextCurrent(m_window);

	// Load OpenGL and its extensions
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		printf("Failed to load OpenGL and its extensions\n");
		exit(EXIT_FAILURE);
	}
	 glfwSetInputMode(m_window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
	/*GLenum res = glewInit();
    if(res != GLEW_OK)
   {
		std::cerr << "Glew failed to initialize!" << std::endl;
    }
	isFullScreen = false;*/
	glLineWidth(5);

}

void Display::AddKeyCallBack(void(*keyCallback)(GLFWwindow *,int,int,int,int))
{
	glfwSetKeyCallback(m_window,(void(*)(GLFWwindow *,int,int,int,int))keyCallback);//{

}

void Display::AddMouseCallBacks(void (* mousebuttonfun)(GLFWwindow*,int,int,int),void (* scrollfun)(GLFWwindow*,double,double),void (* cursorposfun)(GLFWwindow*,double,double))
{
	glfwSetMouseButtonCallback(m_window,mousebuttonfun);
	glfwSetScrollCallback(m_window, scrollfun);
	glfwSetCursorPosCallback(m_window , cursorposfun);
}

void Display::AddResizeCallBack(void (*windowsizefun)(GLFWwindow*,int,int))
{
	glfwSetWindowSizeCallback(m_window,windowsizefun);
}

void Display::SetScene(void *userPointer)
{
	
	glfwSetWindowUserPointer(m_window,userPointer);
}

void* Display::GetScene()
{
	return glfwGetWindowUserPointer(m_window);
}

bool Display::CloseWindow()
{
	return glfwWindowShouldClose(m_window) !=0;
}

Display::~Display()
{
	glfwDestroyWindow(m_window);
	glfwTerminate();
}

bool Display::IsFullscreen( void )
{
	if((glfwGetWindowMonitor( m_window ) == nullptr) != isFullScreen)
	{
		isFullScreen = (glfwGetWindowMonitor( m_window ) == nullptr);
		return true;
	}
	else 
		return false;
} 

void Display::SwapBuffers()
{
	glfwSwapBuffers(m_window);
	// #ifdef __APPLE__
    //     static bool first_time_hack = true;
    //     if (first_time_hack) {
    //         glfwHideWindow(m_window);
    //         glfwShowWindow(m_window);
    //         first_time_hack = false;
    //     }
	// #endif
}

void Display::PollEvents()
{
	glfwPollEvents();

}
