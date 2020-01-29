//the order for [i,j] is the same as [x,y] on a cartesian plane, and we store
//consecutive rows (fixed y, moving x).

// TODO: can only have one PULL connection, enforce this at connection time.

#ifndef DATABUFFER2D_H
#define DATABUFFER2D_H

#ifndef GLEW_MX
#define GLEW_MX
#endif
#include "GL/glew.h"
GLEWContext* glewGetContext();

#include "glm/glm.hpp"
//#include "glm/core/setup.hpp"
//#include "glm/core/type.hpp"

#include "../connections/InputPlug.h"
#include "../connections/OutputPlug.h"
#include "../connections/Dirtyable.h"

namespace goglu{

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

    Metadata2DTexture2D metadataTexture2D;
    Metadata2DImage2D metadataImage2D;

    enum class StorageType {CPU, TEXTURE2D};
    StorageType sourceStorageType;

    typedef typename T dataType;
    unsigned int mNbElementsX, mNbElementsY;
//    GLuint mGlidBuffer; // if you want to transfer to buffer, transfer through a DataBuffer1D.
    GLuint mGlidTexture2D;
    bool mbHasCpuStorage;
//    bool mbHasBufferStorage;
    bool mbHasTexture2DStorage;
    GLenum mTexture2DInternalFormat; // internal texel format. (http://docs.gl/gl4/glTexStorage1D)
    GLenum mTexture2DSizedInternalFormat;
    GLenum mTexture2DExternalFormat; // external texel format, used to communicate with user, e.g. set or get the texture data. (http://docs.gl/gl4/glGetTexImage)
    GLenum mTexture2DSizedExternalFormat;

public:
    // TODO: also allow tranfers through GPU data directly.
    INPUTPLUG_EXT_PULLPUSH(metadataCpu, Metadata2DCpu, DataBuffer2D<T>);
//    INPUTPLUG_PULLPUSH(DataTexture, GLuint, DataBuffer2D<T>);
//    INPUTPLUG_PULLPUSH(DataTextureBuffer, GLuint, DataBuffer2D<T>);
    // TODO: turn into a PULL_QUERY and query for only certain parts of the data. QueryType would be std::tuple<x, y, sizeX, sizeY>.
//    OUTPUTPLUG_PULL(DataCpu, T*, DataBuffer2D<T>);
//    OUTPUTPLUG_PULLPUSH(metadataCpu, std::tuple<T*, unsigned int, unsigned int>, DataBuffer2D<T>); // <dataPointer, nbElementsX, nbElementsY>
//    OUTPUTPLUG_PULLPUSH(DataTexture2D, GLuint, DataBuffer2D<T>);
    OUTPUTPLUG_EXT_PULLPUSH(metadataTexture2D, Metadata2DTexture2D, DataBuffer2D<T>);
//    OUTPUTPLUG_PULLPUSH(DataBuffer, GLuint, DataBuffer2D<T>);
    INPUTPLUG_EXT_PULLPUSH(metadataImage2D, Metadata2DImage2D, DataBuffer2D<T>);
    OUTPUTPLUG_EXT_PULL_QUERY(metadataImage2D, Metadata2DImage2D, unsigned int, DataBuffer2D<T>);

public:
    DIRTYABLE_DATA(dataCpu, T*, DataBuffer2D<T>);
    DIRTYABLE_DATA(dataTexture2D, GLuint, DataBuffer2D<T>);
//    DIRTYABLE_DATA(dataBuffer, GLuint, DataBuffer2D<T>);

public:

    // TODO: when deleting a storage, update all other storages (e.g. so that it gets to the GPU before we delete it from the CPU)

    DataBuffer2D(unsigned int sizeX, unsigned int sizeY);
//    T& operator()(unsigned int x, unsigned int y);
    void createCpuStorage();
    void deleteCpuStorage();
    void createBufferStorage();
    void deleteBufferStorage();
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

    void loadFromFile(std::string filename);

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

    unsigned int dataSizeInBytes() { return mNbElementsX * mNbElementsY * sizeof(T);}
};

} // namespace goglu.

// include definitions because the class is templated.
#include "DataBuffer2D.cpp"

#endif // DATABUFFER2D_H
