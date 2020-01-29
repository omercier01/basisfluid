
#include "Application.h"

bool Application::Init_ComputeShaders() {
    
    Init_CSIntegrateAvgBasis();
    Init_CSIntegrateBasisBasis();
    Init_CSIntegrateBasisGradBasis();
    Init_CSIntegrateBasisGrid_onePerDispatch();
    Init_CSIntegrateBasisGrid_onePerInvocation();

    return true;
}

bool Application::Init_CSIntegrateAvgBasis() {
    return true;
}

bool Application::Init_CSIntegrateBasisBasis() {
    return true;
}

bool Application::Init_CSIntegrateBasisGradBasis() {
    return true;
}

bool Application::Init_CSIntegrateBasisGrid_onePerDispatch() {
    return true;
}

bool Application::Init_CSIntegrateBasisGrid_onePerInvocation() {
    return true;
}