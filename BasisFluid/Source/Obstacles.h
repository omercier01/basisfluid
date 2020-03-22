#include "Application.h"

class Obstacle {
public:
    std::function<float(glm::vec2)> phi; // assumed eikonal, to be updated by updatePhi().
    std::function<glm::vec2(glm::vec2)> gradPhi; // assumed unit norm, to be updated by updatePhi().

    bool dynamic; // true is obstacle moves or deforms during simulation
    std::function<float(glm::vec2)> prevPhi; // only used for dynamic obstacles, copied from phi;
    std::function<glm::vec2(glm::vec2)> prevGradPhi; // only used for dynamic obstacles, copied from gradPhi;

    Obstacle() {
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
    }

    virtual void updatePhi() {} // to update dynamic obstacles
};

class ObstacleCircle : public Obstacle {
public:

    ObstacleCircle() {
        dynamic = true;
    }

    void updatePhi() override {
        _frameCount++;

        glm::vec2 center = app->_obstacleCircleMotionRadius*glm::vec2(sin(app->_obstacleCircleMotionSpeed*app->_dt*_frameCount),
            cos(app->_obstacleCircleMotionSpeed*app->_dt*_frameCount));
        phi = [=](glm::vec2 p) {
            float r = app->_obstacleCircleRadius;
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

struct ObstacleBar : Obstacle {

    ObstacleBar() {
        dynamic = true;
    }

    void updatePhi() override {

        _frameCount++;

        float theta = app->_obstacleBarRotationSpeed * app->_dt * _frameCount;
        float widthX = app->_obstacleBarWidth;
        float widthY = app->_obstacleBarHeight;

        float t2 = app->_obstacleBarMotionSpeed * app->_dt * _frameCount;
        vec2 c = vec2(app->_obstacleBarMotionAmplitude * sin(t2), 0);

        phi = [=](vec2 p) {

            vec2 axisX = vec2(-sin(theta), cos(theta));
            vec2 axisY = vec2(-axisX.y, axisX.x);

            float x = glm::dot(p - c, axisX);
            float y = glm::dot(p - c, axisY);

            if (x > widthX && y > widthY) {
                return std::sqrt(Sqr(x - widthX) + Sqr(y - widthY));
            }
            else if (x > widthX && y < -widthY) {
                return std::sqrt(Sqr(x - widthX) + Sqr(y + widthY));
            }
            else if (x < -widthX && y > widthY) {
                return std::sqrt(Sqr(x + widthX) + Sqr(y - widthY));
            }
            else if (x < -widthX && y < -widthY) {
                return std::sqrt(Sqr(x + widthX) + Sqr(y + widthY));
            }
            else {
                return glm::max<float> (abs(x) - widthX, abs(y) - widthY);
            }
        };

        gradPhi = [=](vec2 p) {
            float eps = 1e-4f;
            vec2 diff = vec2(
                phi(p + vec2(eps, 0)) - phi(p + vec2(-eps, 0)),
                phi(p + vec2(0, eps)) - phi(p + vec2(0, -eps))
            );
            if (glm::length(diff) > 1e-4) { diff = glm::normalize(diff); }
            return diff;
        };
    }

    unsigned int _frameCount = 0;
};
