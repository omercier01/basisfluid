// 1D array of data. We can decide if the data lives on the CPU or GPU, what is
// format is (templates?). Can be reinterpreted as 2D, 3D, or any of the
// formats in OpenGL (maybe even more).

// The order of texture data e.g. 3D seems to be W.H.z + W.y + x, but I have to
// verify that.

//store on CPU and send to GPU if needed. Then we can map the GPU data to the
//CPU to edit it. If we edit it and the mapped region is not mapped, we
//generate an error. If data is only stores on the GPU, the mapped region is 0
//for X Y and Z. THe DataBuffer can of course be used as a CPU data storage if
//the data is not sent to the GPU.

//For, either store all on CPU, all on GPU, or both and have a system to
//transfer the info. I can then optimize and specialize this later. FOr
//instance, for now the transferCpuDataToGpu and transferGpuDataToCpu functions
//wil update the whole data, but later it will be possible to only transfer a
//certain range of this information. This is trickier for non-1D textures, e.g.
//since a 2D range has to be mapped to multiple 1D slices on the GPU. Don't
//know yet how to implement it, either with bufferSubData stuff o with
//glMapBuffer stuff.

//we use immutable textures, but we can still resize the CPU texture and create
//a new GPU texture with the new size.

//we can store the data on the GPU in a buffer or in a texture, or both.
//Textures are useful for rendering the data, and buffers are useful for
//writting to the data from the GPU, or using the data as e.g. a vertex buffer.
//1D textures can also be stored as buffer textures. 2D textures can have
//texture storage or 1D_array storage. etc for 3D.

// TODO: use data format conversion with lambdas. i.e. say something like
// copyDataTo(srcBuffer, dstBuffer, <optional lambda function for conversion, or else
// default casting is used>)

// TODO: Use texture views to reinterpret the data?

// TODO: "const" some of this stuff?

// TODO: add flags specification for glBufferStorage

// TODO: make sure to delete all CPU and GPU storage when the DataBuffer object is deleted.

//TODO: do I really need to specify the number of mipmap levels in
//glTextureStorage1D, or does glGenerateTextureMipmap create them anyway?

//TODO: utiliser des samplers pour controler les stats de la texture, plutot que
//d'utiliser le sampler par defaut de la texture?

//TODO: some warning messages are missing for transferCPUDataToGpu. Same for Texture2D.

//TODO: rename the size variables (e.g. in resize) like in the 2D definitions,
//it is much clearer than always using the nbDataElements.

//TODO: check size compatibility at plug time.

//TODO: problems if the buffer has zero elements?


#ifndef DATABUFFER1D_H
#define DATABUFFER1D_H

#ifndef GLEW_MX
#define GLEW_MX
#endif
#include "GL/glew.h"
GLEWContext* glewGetContext();

#include "../connections/InputPlug.h"
#include "../connections/OutputPlug.h"
#include "../connections/Dirtyable.h"

namespace goglu{

struct Metadata1DCpu {
    void* dataPointer;
//    unsigned int nbElements;
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
//private:

    Metadata1DBuffer metadataBuffer;
    Metadata1DCpu metadataCpu;

    enum class StorageType {UNKNOWN, CPU, BUFFER, TEXTURE1D, TEXTUREBUFFER};
    StorageType sourceStorageType;

    unsigned int capacity;
    unsigned int mNbElements;
    GLuint mGlidBuffer;
    GLuint mGlidTexture1D;
    GLuint mGlidTextureBuffer;
    bool mbHasCpuStorage;
    bool mbHasBufferStorage;
    bool mbHasTexture1DStorage;
    bool mbHasTextureBufferStorage;
    GLenum mTexture1DInternalFormat; // internal texel format. (http://docs.gl/gl4/glTexStorage1D)
    GLenum mTexture1DSizedInternalFormat;
    GLenum mTexture1DExternalFormat; // external texel format, used to communicate with user, e.g. set or get the texture data. (http://docs.gl/gl4/glGetTexImage)
    GLenum mTexture1DSizedExternalFormat;
    GLenum mTextureBufferSizedInternalFormat; // internal texture buffer format. (http://docs.gl/gl4/glTexBuffer)

public:
    // TODO: also allow tranfers through GPU data directly.
    // TODO: only have outputs for the metadata, which includes the data pointer, sizes, format, etc, and no outputs for the data itself?
    INPUTPLUG_EXT_PULLPUSH(metadataCpu, Metadata1DCpu, DataBuffer1D<T>);
//    INPUTPLUG_PULLPUSH(DataTexture, GLuint, DataBuffer1D<T>);
//    INPUTPLUG_PULLPUSH(DataTextureBuffer, GLuint, DataBuffer1D<T>);
//    OUTPUTPLUG_PULLPUSH(DataCpu, T*, DataBuffer1D<T>);
//    OUTPUTPLUG_PULLPUSH(metadataCpu, std::tuple<T*, unsigned int>, DataBuffer1D<T>); // <dataPointer, nbElements>, format is inferred from the *dataPointer
//    OUTPUTPLUG_PULLPUSH(DataTexture, GLuint, DataBuffer1D<T>);
//    OUTPUTPLUG_PULLPUSH(DataBuffer, GLuint, DataBuffer1D<T>);
    INPUTPLUG_EXT_PULLPUSH(metadataBuffer, Metadata1DBuffer, DataBuffer1D<T>);
    OUTPUTPLUG_EXT_PULLPUSH(metadataBuffer, Metadata1DBuffer, DataBuffer1D<T>);
//    OUTPUTPLUG_PULLPUSH(DataTextureBuffer, GLuint, DataBuffer1D<T>);
    OUTPUTPLUG_EXT_PULL(nbParticles, unsigned int, DataBuffer1D<T>);

public:
    DIRTYABLE_DATA(dataCpu, T*, DataBuffer1D<T>);
    DIRTYABLE_DATA(dataBuffer, GLuint, DataBuffer1D<T>);

public:

    // TODO: when deleting a storage, update all other storages (e.g. so that it gets to the GPU before we delete it from the CPU)

    DataBuffer1D(unsigned int size);
//    T& operator()(unsigned int x);
    void createCpuStorage();
    void deleteCpuStorage();
    void createBufferStorage(GLenum dataType, unsigned int nbElementsPerComponent, GLbitfield flags = GL_DYNAMIC_STORAGE_BIT | GL_MAP_READ_BIT | GL_MAP_WRITE_BIT);
    void deleteBufferStorage();
    void createTexture1DStorage(GLenum internalFormat, GLenum sizedInternalFormat,
                                GLenum externalFormat, GLenum sizedExternalFormat);
    void deleteTexture1DStorage();
    void createTextureBufferStorage(GLenum sizedFormat);
    void deleteTextureBufferStorage();

//    void transferCpuDataToGpu();
//    void transferBufferDataToCpu();
//    void transferTexture1DDataToCpu();

    void resize(unsigned int size);
    void appendCpu(T elem);

//    void bindBufferToTarget(GLenum target);

    T getCpuData(unsigned int i);
    T getCpuData_noRefresh(unsigned int i) const;
    T* getCpuDataPointer();
    void refreshCpuData();
    void setCpuData(unsigned int i, T data);
    void setCpuData_noDirty(unsigned int i, T data);
    void dirtyData();
    
    // TODO:
    //void copyFrom(DataBuffer* src);

    unsigned int dataCpuSizeInBytes() {
        return mNbElements * sizeof(T);
        //return metadataCpu.mNbElements * metadataCpu.nbElementsPerComponent * sizeof(T)
    }
};

} // namespace goglu.

// include definitions because the class is templated.
#include "DataBuffer1D.cpp"

#endif // DATABUFFER1D_H
