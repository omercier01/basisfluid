
#include "Application.h"
#include "Obstacles.h"

using namespace std;
using namespace glm;

bool Application::Init_Obstacles() {
    
    // domain boundaries
    _obstacles.push_back(new Obstacle(
        [=](vec2 p) {return p.x - _domainLeft; },
        [=](vec2 /*p*/) {return vec2(1, 0); }
    ));
    _obstacles.push_back(new Obstacle(
        [=](vec2 p) {return -(p.x - _domainRight); },
        [=](vec2 /*p*/) {return vec2(-1, 0); }
    ));
    _obstacles.push_back(new Obstacle(
        [=](vec2 p) {return p.y - _domainBottom; },
        [=](vec2 /*p*/) {return vec2(0, 1); }
    ));
    _obstacles.push_back(new Obstacle(
        [=](vec2 p) {return -(p.y - _domainTop); },
        [=](vec2 /*p*/) {return vec2(0, -1); }
    ));

    // rotating circle
    _obstacles.push_back(new ObstacleCircle());


    for (Obstacle* obs : _obstacles)
    {
        if (obs->dynamic)
        {
            obs->updatePhi();
        }
    }

    return true;
}