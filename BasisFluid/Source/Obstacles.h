#include "Application.h"

class Obstacle {
public:
    std::function<float(glm::vec2)> phi; // assumed eikonal, to be updated by updatePhi().
    std::function<glm::vec2(glm::vec2)> gradPhi; // assumed unit norm, to be updated by updatePhi().

    bool dynamic; //  can move during the simulation
    bool visible;
    std::function<float(glm::vec2)> prevPhi; // only used for dynamic obstacles, copied from phi;
    std::function<glm::vec2(glm::vec2)> prevGradPhi; // only used for dynamic obstacles, copied from gradPhi;


    Obstacle(
    ) {
        visible = true;
        dynamic = true;
    }

    Obstacle(
        std::function<float(glm::vec2)> inPhi,
        std::function<glm::vec2(glm::vec2)> inGradPhi,
        bool inVisible = true
    ) {
        phi = inPhi;
        gradPhi = inGradPhi;

        prevPhi = phi;
        prevGradPhi = gradPhi;

        dynamic = false;
        visible = inVisible;
    }

    virtual void updatePhi() {} // to update dynamic obstacles
};

class ObstacleCircle : public Obstacle {
public:

    ObstacleCircle()
    {
        dynamic = true;
    }

    void updatePhi() override {

        _frameCount++;

        glm::vec2 center = 0.25f*glm::vec2(sin(app->_obstacleSpeed*app->_dt*_frameCount),
            cos(app->_obstacleSpeed*app->_dt*_frameCount));
        phi = [=](glm::vec2 p) {
            float r = app->_obstacleRadius;
            glm::vec2 diff = p - center;
            return -r + glm::length(diff);
        };
        gradPhi = [=](vec2 p) {
            glm::vec2 diff = p - center;
            return glm::length(diff) < 0.0001 ? glm::vec2(0) : glm::normalize(diff);
        };
    }

    unsigned int _frameCount = 0;

};
