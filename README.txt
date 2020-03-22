
2D implementation of "Local Bases for Model-reduced Smoke Simulations" with dynamic obstacles.

Implementation by Olivier Mercier <oli.mercier@gmail.com>.

On Windows, running "VisualStudio\Build\x64\Release\BasisFluid.exe" directly should work. The application will create a Data folder in the folder from which the application is started, so we suggest to run the application from the project's root folder (same folder as this readme file).

All useful parameters are near the top of Application.h .

--------
CONTROLS
--------

Space: toggle stepping of simulation
S: toggle particle seeding
V: toggle velocity grid display
F: toggle particle buoyancy
P: toggle draw particles
M: toggle move obstacles
Mouse scroll: change velocity arrow display size

--------
DEPENDENCIES
--------

glm (glm.g-truc.net)
headers only, included in Libs/, should work as is on any platform.

GLFW (glfw.org) and GLEW (glew.sourceforge.net)
Windows x64 binaries are included in Libs/. You might need to download the binaries for you platform or compile their sources yourself. The Visual Studio project links to the binaries's .lib and copies the .dll to the executable folder. Change these dependencies to match your own if needed.

stb (https://github.com/nothings/stb)
header only, should work as is on any platform. Only stb_image_write.h is needed and is included in Libs/ .
