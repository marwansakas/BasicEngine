# Architecture Notes

BasicEngine is organized as a small C++ graphics engine plus a sample game/application layer.

## Engine Layer

`Engine3D/` contains reusable rendering abstractions such as shapes, meshes, shaders, cameras, scene logic, and transformation handling.

## Display Layer

`DisplayGLFW/` owns the GLFW windowing integration and input callbacks. It is the boundary between the engine and the operating-system/windowing APIs.

## Application Layer

`Game/` contains the executable target. It configures camera/display values, creates the scene, and runs the render loop.

## Build

`CMakeLists.txt` builds `Display`, `Engine`, and the `Game` executable, then links GLFW and GLAD from `res/includes`.

