#ifndef UTILS_H
#define UTILS_H

#include <time.h>
#include <chrono>
#include <iostream>
#include <sstream>

#define M_PI 3.141592653589793238462643

inline bool AllBitsSet(unsigned int bitfield, unsigned int mask) {
    return !(~bitfield & mask);
}

inline bool AtLeastOneBitNotSet(unsigned int bitfield, unsigned int mask) {
    return (~bitfield & mask) != 0;
}

inline unsigned int SetBits(unsigned int bitfield, unsigned int mask) {
    return bitfield | mask;
}

inline unsigned int UnsetBits(unsigned int bitfield, unsigned int mask) {
    return bitfield & ~mask;
}

inline float RoundToMultiple(float val, float step) {
    return round(val / step)*step;
}

template <typename T>
inline T Sqr(T a) { return a * a; }

inline unsigned int SizeOfEnumType(GLenum enumType) {
    switch (enumType) {
    case GL_BYTE:
    case GL_UNSIGNED_BYTE:
        return 1;
        break;
    case GL_SHORT:
    case GL_UNSIGNED_SHORT:
    case GL_HALF_FLOAT:
        return 2;
        break;
    case GL_INT:
    case GL_UNSIGNED_INT:
    case GL_FLOAT:
    case GL_FIXED:
    case GL_INT_2_10_10_10_REV:
    case GL_UNSIGNED_INT_2_10_10_10_REV:
    case GL_UNSIGNED_INT_10F_11F_11F_REV:
        return 4;
        break;
    case GL_DOUBLE:
        return 8;
        break;
    default:
        return 0;
        break;
    }
}

template <typename T>
inline bool IsInClosedInterval(T x, T min, T max) {
    return min <= x && x <= max;
}


// if a linear function has value a at 0 and b at 1, gives the location of the zero. Gives 0.5 if a==b.
inline float RatioZero(float a, float b) {
    if (abs(a - b) < 1e-5) return 0.5;
    else return a / (a - b);
}

inline float VecNorm(glm::vec2 vec) {
    return std::sqrt(vec.x*vec.x + vec.y*vec.y);
}


inline float VecNorm(glm::vec3 vec) {
    return std::sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}

inline void PrintTime() {
    time_t rawtime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    tm timeinfo;
    localtime_s(&timeinfo, &rawtime);
    char buffer[80];
    strftime(buffer, sizeof(buffer), "%d-%m-%Y %H:%M:%S", &timeinfo);
    std::cout << std::string(buffer) << std::endl;
}


#endif // UTILS_H
