//the order for [i,j] is the same as [x,y] on a cartesian plane, and we store
//consecutive rows (fixed y, moving x).

// TODO: can only have one PULL connection, enforce this at connection time.

#ifndef DATABUFFER2D_H
#define DATABUFFER2D_H

#include "GL/glew.h"
#include "glm/glm.hpp"
//#include "glm/core/setup.hpp"
//#include "glm/core/type.hpp"

struct Metadata2DCpu {
    void* dataPointer;
//    unsigned int nbElements;
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
//private:

    Metadata2DTexture2D _metadataTexture2D;
    Metadata2DImage2D _metadataImage2D;

    enum class StorageType {CPU, TEXTURE2D};
    StorageType _sourceStorageType;

    typedef typename T _dataType;
    unsigned int _nbElementsX, _nbElementsY;
//    GLuint mGlidBuffer; // if you want to transfer to buffer, transfer through a DataBuffer1D.
    GLuint _glidTexture2D;
    bool _hasCpuStorage;
//    bool mbHasBufferStorage;
    bool _hasTexture2DStorage;
    GLenum _texture2DInternalFormat; // internal texel format. (http://docs.gl/gl4/glTexStorage1D)
    GLenum _texture2DSizedInternalFormat;
    GLenum _texture2DExternalFormat; // external texel format, used to communicate with user, e.g. set or get the texture data. (http://docs.gl/gl4/glGetTexImage)
    GLenum _texture2DSizedExternalFormat;

public:
    T* _dataCpu;
    GLuint _dataTexture2D;
//    DIRTYABLE_DATA(dataBuffer, GLuint, DataBuffer2D<T>);

public:

    // TODO: when deleting a storage, update all other storages (e.g. so that it gets to the GPU before we delete it from the CPU)

    DataBuffer2D(unsigned int sizeX, unsigned int sizeY);
//    T& operator()(unsigned int x, unsigned int y);
    void createCpuStorage();
    void deleteCpuStorage();
    //void createBufferStorage();
    //void deleteBufferStorage();
    void createTexture2DStorage(GLenum internalFormat, GLenum sizedInternalFormat,
                                GLenum externalFormat, GLenum sizedExternalFormat,
                                unsigned int nbMipmapLevels);
    void createTexture2DStorage(GLenum internalFormat, GLenum sizedInternalFormat,
                                GLenum externalFormat, GLenum sizedExternalFormat);
    void deleteTexture2DStorage();

//    void transferCpuDataToGpu();
//    void transferBufferDataToCpu();
//    void transferTexture2DDataToCpu();

    void resize(unsigned int newSizeX, unsigned int newSizeY);

    void bindBufferToTarget(GLenum target);

    T getCpuData(unsigned int i, unsigned int j);
    T getCpuData_noRefresh(unsigned int i, unsigned int j) const;
    T* getCpuDataPointer();
    void refreshCpuData();
    void setCpuData(unsigned int i, unsigned int j, T data);
    void addCpuData(unsigned int i, unsigned int j, T data);
    void dirtyData();
    // TODO: getter and setters for Texture and TextureBuffer similar to set/getCpu().

    void setFromCpuData(T* data);
    void setFromImage(Metadata2DImage2D metadataImage2D, unsigned int level);
    // TODO:
    //void copyFrom(DataBuffer* src);

    T interp(glm::vec2 pos);

    unsigned int dataSizeInBytes() { return _nbElementsX * _nbElementsY * sizeof(T);}
};

// include definitions because the class is templated.
#include "DataBuffer2D.tpp"

#endif // DATABUFFER2D_H
