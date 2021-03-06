
C++17 implementation of "Local Bases for Model-reduced Smoke Simulations" in 2D with dynamic obstacles. See www.olivier-mercier.com/res/publications/basisFluid/ for more details.

Implementation by Olivier Mercier <oli.mercier@gmail.com>.

All useful parameters are near the top of BasisFluid/Source/Application.h .

On Windows, running "VisualStudio\Build\x64\Release\BasisFluid.exe" directly should work. The application can create a Data and Output folder in the folder from which the application is started. Although not required, we suggest to run the application from the project's root folder (same folder as this readme file).

We also include a Visual Studio 2017 solution (VisualStudio/BasisFluid.sln).

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
Headers only, included in Libs/, should work as is on any platform.

GLFW (glfw.org) and GLEW (glew.sourceforge.net)
Windows x64 binaries are included in Libs/. You might need to download the binaries for you platform or compile their sources yourself. The Visual Studio project links to the binaries's .lib and copies the .dll to the executable folder. Change these dependencies to match your own if needed.

stb (https://github.com/nothings/stb)
Headers only, should work as is on any platform. Only stb_image_write.h is needed and included in Libs/ .
