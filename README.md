# BasicEngine

BasicEngine is a C++ OpenGL 3D graphics engine fork used for a computer graphics course. It helps students experiment with scenes, meshes, shaders, camera setup, input handling, and GLFW-based rendering.

The project is based on an educational OpenGL engine structure and keeps the original GPL-3.0 license from the fork.

## Tech Stack

- C++11
- CMake
- OpenGL
- GLFW
- GLAD

## Quick Start

Install CMake and a C++ compiler, then configure and build:

```bash
cmake -S . -B build
cmake --build build
```

Run the generated `Game` executable from the build output.

## Environment Notes

No `.env` file or environment variables are required. The important local configuration is the graphics toolchain: CMake, a C++11 compiler, OpenGL-capable drivers, and the bundled GLFW/GLAD dependencies under `res/includes`.

## Usage Notes

The `Game/` folder is the application layer. It creates the display, attaches a scene, and drives the render loop. `Engine3D/` contains reusable engine abstractions, and `DisplayGLFW/` contains the GLFW display/input integration.

## Testing

There is no automated rendering test suite yet. The CI workflow configures and builds the CMake project so compile regressions are visible.

## Demo / Walkthrough

This is a desktop OpenGL project rather than a hosted demo. A reviewer can build the project locally, run the `Game` target, and inspect the scene implementation under `Game/`.

## Project Structure

- `Game/`: sample application and scene code.
- `Engine3D/`: engine, mesh, shader, shape, camera, and scene abstractions.
- `DisplayGLFW/`: display, input, and GLFW integration.
- `res/`: bundled third-party headers/libraries and runtime resources.
- `readme.txt`: original course notes.
- `docs/`: architecture notes added for GitHub review.
