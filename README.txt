


On Windows, running "VisualStudio\Build\x64\Release\BasisFluid.exe" directly might work.


--------
DEPENDENCIES
--------
- glm: headers only, includes in Libs/, should work as is for any platform.
- GLFW and GLEW: Windows x64 binaries included in Libs/. You might need to download the binaries for you platform or compile the sources yourself. The Visual Studio project links to the binaries directly, so change its dependencies to your own if needed.
