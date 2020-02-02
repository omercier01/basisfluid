
#include "Application.h"

void Application::Draw()
{

    if (_drawObstacles) {
        // update obstacle display

        if (obstacleDisplayNeedsUpdating) {

            obstacleLines->resize(0);

            for (Obstacle* obs : obstacles) {

                int nbOffsets = 9;
                float thickness = 0.005;

                for (int offset = -nbOffsets; offset <= nbOffsets; offset++) {

                    // obstacle boundary with marching square
                    float dx = (domainRight - domainLeft) / obstacleDisplayRes;
                    float dy = (domainTop - domainBottom) / obstacleDisplayRes;
                    for (int i = -2; i <= int(obstacleDisplayRes) + 1; i++) {
                        for (int j = -2; j <= int(obstacleDisplayRes) + 1; j++) {

                            // bl = bottom left, tr = top right.
                            vec2 blPos = vec2(domainLeft + i * dx, domainBottom + j * dy);
                            float bl = thickness * float(offset) / float(nbOffsets) + obs->phi(blPos);
                            float br = thickness * float(offset) / float(nbOffsets) + obs->phi(blPos + vec2(dx, 0));
                            float tl = thickness * float(offset) / float(nbOffsets) + obs->phi(blPos + vec2(0, dy));
                            float tr = thickness * float(offset) / float(nbOffsets) + obs->phi(blPos + vec2(dx, dy));

                            short flags = 8 * (bl >= 0 ? 1 : 0) + 4 * (br >= 0 ? 1 : 0) + 2 * (tl >= 0 ? 1 : 0) + 1 * (tr >= 0 ? 1 : 0);
                            switch (flags) {
                            case 0:
                            case 15:
                                break;
                            case 8:
                            case 7:
                                obstacleLines->appendCpu(blPos + vec2(0, dy*ratioZero(bl, tl)));
                                obstacleLines->appendCpu(blPos + vec2(dx*ratioZero(bl, br), 0));
                                break;
                            case 4:
                            case 11:
                                obstacleLines->appendCpu(blPos + vec2(dx*ratioZero(bl, br), 0));
                                obstacleLines->appendCpu(blPos + vec2(dx, dy*ratioZero(br, tr)));
                                break;
                            case 1:
                            case 14:
                                obstacleLines->appendCpu(blPos + vec2(dx, dy*ratioZero(br, tr)));
                                obstacleLines->appendCpu(blPos + vec2(dx*ratioZero(tl, tr), dy));
                                break;
                            case 2:
                            case 13:
                                obstacleLines->appendCpu(blPos + vec2(dx*ratioZero(tl, tr), dy));
                                obstacleLines->appendCpu(blPos + vec2(0, dy*ratioZero(bl, tl)));
                                break;
                            case 12:
                            case 3:
                                obstacleLines->appendCpu(blPos + vec2(0, dy*ratioZero(bl, tl)));
                                obstacleLines->appendCpu(blPos + vec2(dx, dy*ratioZero(br, tr)));
                                break;
                            case 10:
                            case 5:
                                obstacleLines->appendCpu(blPos + vec2(dx*ratioZero(bl, br), 0));
                                obstacleLines->appendCpu(blPos + vec2(dx*ratioZero(tl, tr), dy));
                                break;
                            case 9:
                            case 6: // crossed case, arbitrary choice is made for line orientation.
                                obstacleLines->appendCpu(blPos + vec2(0, dy*ratioZero(bl, tl)));
                                obstacleLines->appendCpu(blPos + vec2(dx*ratioZero(bl, br), 0));
                                obstacleLines->appendCpu(blPos + vec2(dx, dy*ratioZero(br, tr)));
                                obstacleLines->appendCpu(blPos + vec2(dx*ratioZero(tl, tr), dy));
                                break;
                            }
                        }
                    }
                }
            }

            obstacleDisplayNeedsUpdating = false;
        }

        pipelineObstacle->goglu_nbPrimitives.set(obstacleLines->metadataBuffer.nbElements / 2);
    }


    // display all simulation elements
    glClearColor(1, 1, 1, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    glClear(GL_DEPTH_BUFFER_BIT);
    if (tw->drawGrid.get()) pipelineGrid->execute();
    pipelineParticles->goglu_nbPrimitives.set(particles->metadataBuffer.nbElements);
    if (tw->drawBasisContours.get()) pipelineBasisContour->execute();
    //    pipelineBasisContour->execute();

    if (tw->drawObstacle.get() && obstacleLines->metadataBuffer.nbElements > 0) pipelineObstacle->execute();


    if (tw->displayVelocityGrid.get()) {
        pipelineArrows->execute();
    }
    if (tw->drawParticles.get()) { pipelineParticles->execute(); }

    // tweak bar
    controlWindow->makeContextCurrent();
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    tw->draw();
    mainWindow->makeContextCurrent();



    glfwSwapBuffers(_glfwWindow);

}


// used for display only
void Application::ComputeVelocityGridForDisplay()
{


    // clear velocity grid
    // TODO: speed this up
    _velocityField->populateWithFunction([](float /*x*/, float /*y*/) {return vec2(0); });

    bool displayVelocityStretched = tw->displayVelocityStretched.get();




    for (unsigned int i = 0; i < _basisFlowParams->_nbElements; ++i) {
        BasisFlow b = _basisFlowParams->getCpuData_noRefresh(i);


        //TODO: Use the decompressed intersection per grid node instead, should be much faster.
        _velocityField->addFunction(
            [&](float x, float y) {
            vec2 p = vec2(x, y);

            vec2 vec(0);
            if (AllBitsSet(b.bitFlags, INTERIOR)) {
                vec += vecObstacle_stretch(p, b);
            }
            vec += b.coeffBoundary * (vec2)(
                TranslatedBasisEval(
                    p, b.freqLvl, b.center
                )
                )
                ;
            return vec;
        }
        );
    }
}