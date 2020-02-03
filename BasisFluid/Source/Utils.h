#ifndef UTILS_H
#define UTILS_H

#define M_PI 3.141592653589793238462643

inline bool AllBitsSet(unsigned int bitfield, unsigned int mask) {
    return !(~bitfield & mask);
}

inline bool AtLeastOneBitNotSet(unsigned int bitfield, unsigned int mask) {
    return (~bitfield & mask) != 0;
}

unsigned int SetBits(unsigned int bitfield, unsigned int mask) {
    return bitfield | mask;
}

unsigned int UnsetBits(unsigned int bitfield, unsigned int mask) {
    return bitfield & ~mask;
}

inline float RoundToMultiple(float val, float step) {
    return round(val/step)*step;
}

template <typename T>
inline T Sqr(T a) {return a*a;}

inline unsigned int SizeOfEnumType(GLenum enumType) {
    switch(enumType) {
    case GL_BYTE:
    case GL_UNSIGNED_BYTE:
//        std::cout << "RETURNING 1" << std::endl;
        return 1;
        break;
    case GL_SHORT:
    case GL_UNSIGNED_SHORT:
    case GL_HALF_FLOAT:
//        std::cout << "RETURNING 2" << std::endl;
        return 2;
        break;
    case GL_INT:
    case GL_UNSIGNED_INT:
    case GL_FLOAT:
    case GL_FIXED:
    case GL_INT_2_10_10_10_REV:
    case GL_UNSIGNED_INT_2_10_10_10_REV:
    case GL_UNSIGNED_INT_10F_11F_11F_REV:
//        std::cout << "RETURNING 4" << std::endl;
        return 4;
        break;
    case GL_DOUBLE:
//        std::cout << "RETURNING 8" << std::endl;
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
float RatioZero(float a, float b) {
    if(abs(a-b) < 1e-5) return 0.5;
    else return a/(a-b);
}



#endif // UTILS_H
