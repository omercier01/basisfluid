
// Generic 2D buffer. The order for [i,j] is the same as [x,y] on a cartesian plane, x is the
// fast moving variable in linear order.

#ifndef DATABUFFER2D_H
#define DATABUFFER2D_H

#include "GL/glew.h"
#include "glm/glm.hpp"

struct Metadata2DCpu {
    void* dataPointer;
};


struct Metadata2DTexture2D {
    GLuint textureId;
    unsigned int nbElementsX;
    unsigned int nbElementsY;
    GLenum internalFormat; // internal texel format. (http://docs.gl/gl4/glTexStorage1D)
    GLenum sizedInternalFormat;
    GLenum externalFormat; // external texel format, used to communicate with user, e.g. set or get the texture data. (http://docs.gl/gl4/glGetTexImage)
    GLenum sizedExternalFormat;
};


struct Metadata2DImage2D {
    GLuint textureId;
    unsigned int nbElementsX;
    unsigned int nbElementsY;
    GLenum internalFormat; // internal texel format. (http://docs.gl/gl4/glTexStorage1D)
    GLenum sizedInternalFormat;
    GLenum externalFormat; // external texel format, used to communicate with user, e.g. set or get the texture data. (http://docs.gl/gl4/glGetTexImage)
    GLenum sizedExternalFormat;
    unsigned int level;
};


template<class T>
class DataBuffer2D
{
public:
    Metadata2DTexture2D _metadataTexture2D;
    Metadata2DImage2D _metadataImage2D;

    typedef typename T _dataType;
    unsigned int _nbElementsX, _nbElementsY;
    GLuint _glidTexture2D;
    bool _hasCpuStorage;
    bool _hasTexture2DStorage;
    GLenum _texture2DInternalFormat; // internal texel format. (http://docs.gl/gl4/glTexStorage1D)
    GLenum _texture2DSizedInternalFormat;
    GLenum _texture2DExternalFormat; // external texel format, used to communicate with user, e.g. set or get the texture data. (http://docs.gl/gl4/glGetTexImage)
    GLenum _texture2DSizedExternalFormat;

public:
    T* _dataCpu;
    GLuint _dataTexture2D;

public:

    DataBuffer2D(unsigned int sizeX, unsigned int sizeY);
    void createCpuStorage();
    void deleteCpuStorage();
    void createTexture2DStorage(GLenum internalFormat, GLenum sizedInternalFormat,
                                GLenum externalFormat, GLenum sizedExternalFormat,
                                unsigned int nbMipmapLevels);
    void createTexture2DStorage(GLenum internalFormat, GLenum sizedInternalFormat,
                                GLenum externalFormat, GLenum sizedExternalFormat);
    void deleteTexture2DStorage();

    void resize(unsigned int newSizeX, unsigned int newSizeY);

    T getCpuData(unsigned int i, unsigned int j);
    T* getCpuDataPointer();
    void setCpuData(unsigned int i, unsigned int j, T data);
    void addCpuData(unsigned int i, unsigned int j, T data);

    void setFromCpuData(T* data);
    void setFromImage(Metadata2DImage2D metadataImage2D, unsigned int level);

    T interp(glm::vec2 pos);

    unsigned int dataSizeInBytes() { return _nbElementsX * _nbElementsY * sizeof(T);}
};

// include definitions because the class is templated.
#include "DataBuffer2D.tpp"

#endif // DATABUFFER2D_H
