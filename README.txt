
2D implementation of "Local Bases for Model-reduced Smoke Simulations" with dynamic obstacles.
Implementation by Olivier Mercier.

On Windows, running "VisualStudio\Build\x64\Release\BasisFluid.exe" directly should work.

All useful paramters are near the top of Application.h .

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

- glm: headers only, included in Libs/, should work as is on any platform.
- GLFW and GLEW: Windows x64 binaries are included in Libs/. You might need to download the binaries for you platform or compile their sources yourself. The Visual Studio project links to the binaries's .lib and copies the .dll to the executable folder. Change these dependencies to match your own if needed.
