
// Generic 1D buffer. Can be stored on CPU or GPU, and transfer between the two.

#ifndef DATABUFFER1D_H
#define DATABUFFER1D_H

#include "GL/glew.h"


struct Metadata1DCpu {
    void* dataPointer;
};


struct Metadata1DBuffer {
    GLuint bufferId;
    unsigned int nbElements;
    GLenum dataType; // e.g. GL_FLOAT
    unsigned int nbElementsPerComponent; // e.g. 3 for vec3
};


template<class T>
class DataBuffer1D
{
public:
    Metadata1DBuffer _metadataBuffer;
    Metadata1DCpu _metadataCpu;

    unsigned int _capacity;
    unsigned int _nbElements;
    GLuint _glidBuffer;
    GLuint _glidTexture1D;
    GLuint _glidTextureBuffer;
    bool _hasCpuStorage;
    bool _hasBufferStorage;
    bool _hasTexture1DStorage;
    bool _hasTextureBufferStorage;
    GLenum _texture1DInternalFormat; // internal texel format. (http://docs.gl/gl4/glTexStorage1D)
    GLenum _texture1DSizedInternalFormat;
    GLenum _texture1DExternalFormat; // external texel format, used to communicate with user, e.g. set or get the texture data. (http://docs.gl/gl4/glGetTexImage)
    GLenum _texture1DSizedExternalFormat;
    GLenum _textureBufferSizedInternalFormat; // internal texture buffer format. (http://docs.gl/gl4/glTexBuffer)

public:
    T* _dataCpu;
    GLuint dataBuffer;

public:

    DataBuffer1D(unsigned int size);
    void createCpuStorage();
    void deleteCpuStorage();
    void createBufferStorage(GLenum dataType, unsigned int nbElementsPerComponent, GLbitfield flags = GL_DYNAMIC_STORAGE_BIT | GL_MAP_READ_BIT | GL_MAP_WRITE_BIT);
    void deleteBufferStorage();
    void createTexture1DStorage(GLenum internalFormat, GLenum sizedInternalFormat,
                                GLenum externalFormat, GLenum sizedExternalFormat);
    void deleteTexture1DStorage();
    void createTextureBufferStorage(GLenum sizedFormat);
    void deleteTextureBufferStorage();

    void resize(unsigned int size);
    void appendCpu(T elem);

    T getCpuData(unsigned int i);
    T* getCpuDataPointer();
    void setCpuData(unsigned int i, T data);
    void TransferDataCpuToBuffer();
    
    unsigned int dataCpuSizeInBytes() {
        return _nbElements * sizeof(T);
    }
};

// include definitions because the class is templated.
#include "DataBuffer1D.tpp"

#endif // DATABUFFER1D_H
