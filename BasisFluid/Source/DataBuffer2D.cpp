// TODO: load different formats, right now it can only load BMP (maybe other
// formats too, but not PNG.)
//

#include "DataBuffer2D.h"
#include "FreeImage.h"
#include "../Application.h"

using namespace goglu;
using namespace std;
using namespace glm;

template<class T>
DataBuffer2D<T>::DataBuffer2D(unsigned int sizeX, unsigned int sizeY) :
    in_metadataCpu(*this),
//    outDataCpu(*this),
    out_metadataTexture2D(*this),
    in_metadataImage2D(*this),
    out_metadataImage2D(*this),
    dataCpu(*this),
    dataTexture2D(*this)
//    dataBuffer(*this)
{
    mNbElementsX = sizeX;
    mNbElementsY = sizeY;
//    mNbDataElements = mNbElementsX*mNbElementsY;
    mbHasCpuStorage = false;
//    mbHasBufferStorage = false;
    mbHasTexture2DStorage = false;

}

//template<class T>
//T& goglu::DataBuffer2D<T>::operator()(unsigned int x, unsigned int y)
//{
//    if(mbHasCpuStorage) {
//        //return data.get()[mSizeX*y + x];
//        return get(x, y);
//    } else {
//        Application::sInsertDebugMessage(
//                "goglu::DataBuffer2D::operator() : DataBuffer has no "
//                "CPU storage, accessing it returns a default value.",
//                GL_DEBUG_SEVERITY_HIGH);
//        T temp;
//        return temp;
//    }

//}

// based of TextureManager::LoadTexture that comes with the FreeImage library.
template<class T>
void goglu::DataBuffer2D<T>::loadFromFile(std::string filename)
{
    //image format
    FREE_IMAGE_FORMAT fif = FIF_UNKNOWN;
    //pointer to the image, once loaded
    FIBITMAP *dib(0);
    //pointer to the image data
    BYTE* bits(0);
    //image width and height
    unsigned int fiWidth(0), fiHeight(0);

    //check the file signature and deduce its format
    fif = FreeImage_GetFileType(filename.c_str(), 0);
    //if still unknown, try to guess the file format from the file extension
    if(fif == FIF_UNKNOWN) {
        fif = FreeImage_GetFIFFromFilename(filename.c_str());
    }
    //if still unkown, return failure
    if(fif == FIF_UNKNOWN) {
        Application::sInsertDebugMessage(
                "goglu::DataBuffer2D::loadFromFile : file type unknown, file "
                + filename + " was not loaded.");
        return;
    }

    if(FreeImage_FIFSupportsReading(fif))
        dib = FreeImage_Load(fif, filename.c_str());
    //if the image failed to load, return failure
    if(!dib) {
        Application::sInsertDebugMessage(
                "goglu::DataBuffer2D::loadFromFile : cannot load file, file "
                + filename + " was not loaded.");
        return;
    }

    bits = FreeImage_GetBits(dib);
    unsigned int width = FreeImage_GetWidth(dib);
    unsigned int height = FreeImage_GetHeight(dib);

    resize(width, height);

    if(!mbHasCpuStorage) {
        createCpuStorage();
    }
    if(!mbHasTexture2DStorage) {
        createTexture2DStorage(GL_RGB, GL_RGB32F, GL_BGR, GL_FLOAT);
    }

    // TODO: do this the clean way.
    for(int i=0; i<width*height; i++) {
        dataCpu.get()[i] = vec3(
                        bits[3*i+0]*1.f/255.f,
                        bits[3*i+1]*1.f/255.f,
                        bits[3*i+2]*1.f/255.f);
    }
    dataCpu.dirtyItselfAndDependents();
}

template<class T>
void DataBuffer2D<T>::createCpuStorage()
{
    if(mbHasCpuStorage) {
        Application::sInsertDebugMessage(
                "goglu::DataBuffer2D::createCpuStorage : DataBuffer already has "
                "CPU storage, creating the storage again has no effect.",
                GL_DEBUG_SEVERITY_LOW);
    } else {
        dataCpu.set(new T[mNbElementsX*mNbElementsY]);
        mbHasCpuStorage = true;
    }
}

template<class T>
void DataBuffer2D<T>::deleteCpuStorage()
{
    if(mbHasCpuStorage) {
        delete[] mpData;
        mbHasCpuStorage = false;
    } else {
        Application::sInsertDebugMessage(
                "goglu::DataBuffer2D::deleteCpuStorage : DataBuffer has no CPU "
                "storage, deleting the storage has no effect.",
                GL_DEBUG_SEVERITY_LOW);
    }
}

template<class T>
void DataBuffer2D<T>::createBufferStorage()
{
    if(mbHasBufferStorage) {
        Application::sInsertDebugMessage(
                "goglu::DataBuffer2D::createBufferStorage : DataBuffer already "
                "has BUFFER storage, creating the storage again has no effect.",
                GL_DEBUG_SEVERITY_LOW);
    } else {
        glCreateBuffers(1, &mGlidBuffer);
        glNamedBufferStorage(mGlidBuffer, dataSizeInBytes(), NULL,
                             GL_DYNAMIC_STORAGE_BIT | GL_MAP_READ_BIT |
                             GL_MAP_WRITE_BIT);
        mbHasBufferStorage = true;
    }
}

template<class T>
void DataBuffer2D<T>::deleteBufferStorage()
{
    if(mbHasBufferStorage) {
        glDeleteBuffers(1, &mGlidBuffer);
        mbHasBufferStorage = false;
    } else {
        Application::sInsertDebugMessage(
                "goglu::DataBuffer2D::deleteBufferStorage : DataBuffer has no "
                "BUFFER storage, deleting the storage has no effect.",
                GL_DEBUG_SEVERITY_LOW);
    }
}

template<class T>
void DataBuffer2D<T>::createTexture2DStorage(
        GLenum internalFormat, // e.g. GL_RGB
        GLenum sizedInternalFormat, // e.g. GL_RGB32F
        GLenum externalFormat, // e.g. GL_RGB
        GLenum sizedExternalFormat) // e.g. GL_FLOAT
{
    // this computes nbLevels = log_2(max(mSizeX,mSizeY)), which is the maximum number of
    // mipmap level we can have for this size.
    unsigned int nbLevels = 1;
    unsigned int tempSize = glm::max(mNbElementsX, mNbElementsY);
    while (tempSize >>= 1) ++nbLevels;

    createTexture2DStorage(internalFormat, sizedInternalFormat,
                           externalFormat, sizedExternalFormat,
                           nbLevels);
}


//nbMipmapLevels=0 will put the maximum number of mipmaps. nbMipmapLevels=1
//will only put the top layer without any mipmaps below it.
template<class T>
void DataBuffer2D<T>::createTexture2DStorage(
        GLenum internalFormat, // e.g. GL_RGB
        GLenum sizedInternalFormat, // e.g. GL_RGB32F
        GLenum externalFormat, // e.g. GL_RGB
        GLenum sizedExternalFormat, // e.g. GL_FLOAT
        unsigned int nbMipmapLevels)
{
    if(mbHasTexture2DStorage) {
        Application::sInsertDebugMessage(
                "goglu::DataBuffer2D::createTexture1DStorage : DataBuffer "
                "already has TEXTURE_2D storage, creating the storage again "
                "has no effect.",
                GL_DEBUG_SEVERITY_LOW);
    } else {
        mTexture2DInternalFormat = internalFormat;
        mTexture2DSizedInternalFormat = sizedInternalFormat;
        mTexture2DExternalFormat = externalFormat;
        mTexture2DSizedExternalFormat = sizedExternalFormat;
        glCreateTextures(GL_TEXTURE_2D, 1, &mGlidTexture2D);
//        unsigned int nbLevels;
//        // this computes nbLevels = log_2(max(mSizeX,mSizeY)), which is the maximum number of
//        // mipmap level we can have for this size.
//        nbLevels = 1;
//        unsigned int tempSize = max(mNbElementsX, mNbElementsY);
//        while (tempSize >>= 1) ++nbLevels;

        glTextureStorage2D(mGlidTexture2D, nbMipmapLevels,
                           mTexture2DSizedInternalFormat, mNbElementsX, mNbElementsY);
        glClearTexImage(mGlidTexture2D, 0, mTexture2DExternalFormat,
                        mTexture2DSizedExternalFormat, NULL);
        glGenerateTextureMipmap(mGlidTexture2D);
        glTextureParameteri(mGlidTexture2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTextureParameteri(mGlidTexture2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        // clamp to zero by default, override if different boundary conditions desired.
        glTextureParameteri(mGlidTexture2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
        glTextureParameteri(mGlidTexture2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
        float borderColor[] = { 0.0f, 0.0f, 0.0f, 1.0f };
        glTextureParameterfv(mGlidTexture2D, GL_TEXTURE_BORDER_COLOR, borderColor);
        mbHasTexture2DStorage = true;

        // metadata
        metadataTexture2D.textureId = mGlidTexture2D;
        metadataTexture2D.nbElementsX = mNbElementsX;
        metadataTexture2D.nbElementsY = mNbElementsY;
        metadataTexture2D.internalFormat = mTexture2DInternalFormat;
        metadataTexture2D.sizedInternalFormat = mTexture2DSizedInternalFormat;
        metadataTexture2D.externalFormat = mTexture2DExternalFormat;
        metadataTexture2D.sizedExternalFormat = mTexture2DSizedExternalFormat;

    }
}

template<class T>
void DataBuffer2D<T>::deleteTexture2DStorage()
{
    if(mbHasTexture2DStorage) {
        glDeleteTextures(1, &mGlidTexture2D);
        mbHasTexture1DStorage = false;
    } else {
        Application::sInsertDebugMessage(
                "goglu::DataBuffer2D::deleteTexture1DSorage : DataBuffer has "
                "no TEXTURE_2D storage, deleting the storage has no effect.",
                GL_DEBUG_SEVERITY_LOW);
    }
}

//template<class T>
//void DataBuffer2D<T>::transferCpuDataToGpu()
//{
//    if(mbHasTexture2DStorage) {
//        glTextureSubImage2D(mGlidTexture2D, 0, 0, 0, mSizeX, mSizeY,
//                            mTexture2DExternalFormat,
//                            mTexture2DSizedExternalFormat, mpData);
//        glGenerateTextureMipmap(mGlidTexture2D);
//    }
//    if(mbHasBufferStorage) {
//        glNamedBufferSubData(mGlidBuffer, 0, dataSizeInBytes(), mpData);
//    }
//}

//template<class T>
//void DataBuffer2D<T>::transferBufferDataToCpu()
//{
//    if(mbHasCpuStorage && mbHasBufferStorage) {
//        glGetNamedBufferSubData( mGlidBuffer, 0, dataSizeInBytes(), mpData);
//    } else {
//        if(!mbHasCpuStorage) {
//            Application::sInsertDebugMessage(
//                    "goglu::DataBuffer2D::transferBufferDataToCpu : DataBuffer "
//                    "has no CPU storage, cannot transfer BUFFER data to it.",
//                    GL_DEBUG_SEVERITY_HIGH);
//        } else {
//            Application::sInsertDebugMessage(
//                    "goglu::DataBuffer2D::transferBufferDataToCpu : DataBuffer "
//                    "has no BUFFER storage, cannot transfer it to CPU.",
//                    GL_DEBUG_SEVERITY_HIGH);
//        }
//    }
//}

//template<class T>
//void DataBuffer2D<T>::transferTexture2DDataToCpu()
//{
//    if(mbHasCpuStorage && mbHasTexture2DStorage) {
//        glGetTextureImage(mGlidTexture2D, 0, mTexture2DExternalFormat,
//                          mTexture2DSizedExternalFormat, dataSizeInBytes(),
//                          mpData);
//    } else {
//        if(!mbHasCpuStorage) {
//            Application::sInsertDebugMessage(
//                    "goglu::DataBuffer2D::transferTexture2DToCpu : DataBuffer "
//                    "has no CPU storage, cannot transfer TEXTURE_2D data to it.",
//                    GL_DEBUG_SEVERITY_HIGH);
//        } else {
//            Application::sInsertDebugMessage(
//                    "goglu::DataBuffer2D::transferTexture2DDataToCpu : "
//                    "DataBuffer has no TEXTURE_2D storage, cannot transfer it "
//                    "to CPU.",
//                    GL_DEBUG_SEVERITY_HIGH);
//        }
//    }
//}


// assumes the size of data is the same as the DataBuffer's size.
template<class T>
void DataBuffer2D<T>::setFromCpuData(T* inData) {

    if(mbHasCpuStorage) {
        for(unsigned int i=0; i<mNbElementsX; i++) {
            for(unsigned int j=0; j<mNbElementsY; j++) {
                setCpuData(i, j, inData[mNbElementsY*j + i]);
            }
        }
    }

    if(mbHasTexture2DStorage) {
        glTextureSubImage2D(mGlidTexture2D, 0, 0, 0, mNbElementsX, mNbElementsY,
                            mTexture2DExternalFormat,
                            mTexture2DSizedExternalFormat, inData);
//        glGenerateTextureMipmap(mGlidTexture2D);
    }

//    if(mbHasBufferStorage) {
//        glNamedBufferSubData(mGlidBuffer, 0, dataSizeInBytes(), inData);
//    }

}

// assumes the size of data is the same as the DataBuffer's size.
template<class T>
void DataBuffer2D<T>::setFromImage(Metadata2DImage2D srcImage, unsigned int level) {

    if(mbHasTexture2DStorage) {
        glCopyImageSubData(
                srcImage.textureId, GL_TEXTURE_2D, srcImage.level,
                0, 0, 0,
                metadataTexture2D.textureId, GL_TEXTURE_2D, level,
                0, 0, 0,
                metadataTexture2D.nbElementsX, metadataTexture2D.nbElementsY, 1);
        sourceStorageType = StorageType::TEXTURE2D;
        dataCpu.dirtyItselfAndDependents();
    } else {
        // TODO: error, texture storage required.
    }


}

template <class T>
T DataBuffer2D<T>::interp(vec2 pos)
{
    float x = pos.x;
    float y = pos.y;

    int indexXLeft, indexXRight;
    int indexYLeft, indexYRight;
    float weightX;
    float weightY;

    // TODO: this is not clean...
    float mBoundXMin = 0;
    float mBoundYMin = 0;
    float mBoundXMax = 1;
    float mBoundYMax = 1;
    int mNbCellsX = mNbElementsX;
    int mNbCellsY = mNbElementsY;

    float normalizedX = (x - mBoundXMin)/(mBoundXMax - mBoundXMin)*mNbCellsX;
    float normalizedY = (y - mBoundYMin)/(mBoundYMax - mBoundYMin)*mNbCellsY;



    if(normalizedX < 0.5) {
        return T();
    } else if(normalizedX >= mNbCellsX - 0.5f) {
        return T();
    } else {
        indexXLeft = clamp<int>(int(floor(normalizedX - 0.5f)), 0, mNbCellsX-2);
        indexXRight = indexXLeft + 1;
        weightX = normalizedX - (indexXLeft + 0.5f);
    }
    if(normalizedY < 0.5) {
        return T();
    } else if(normalizedY >= mNbCellsY - 0.5f) {
        return T();
    } else {
        indexYLeft = clamp<int>(int(floor(normalizedY - 0.5f)), 0, mNbCellsY-2);
        indexYRight = indexYLeft + 1;
        weightY = normalizedY - (indexYLeft + 0.5f);
    }



    //TODO: indices might still be out of range because of numerical errors,
    //check for that and if it happens, the out-of-bound indices will have to
    //be dealt with in terms of indices, not positions.

//    return (1-weightX)*(1-weightY)*mVectors(indexXLeft , indexYLeft ) +
//           (1-weightX)*(  weightY)*mVectors(indexXLeft , indexYRight) +
//           (  weightX)*(1-weightY)*mVectors(indexXRight, indexYLeft ) +
//           (  weightX)*(  weightY)*mVectors(indexXRight, indexYRight);

    return (1-weightX)*(1-weightY)*getCpuData(indexXLeft , indexYLeft ) +
           (1-weightX)*(  weightY)*getCpuData(indexXLeft , indexYRight) +
           (  weightX)*(1-weightY)*getCpuData(indexXRight, indexYLeft ) +
            (  weightX)*(  weightY)*getCpuData(indexXRight, indexYRight);
}

// TODO: instead do something similar to DataBuffer1D's resize with the *2 scheme to amortize costs.
template<class T>
void DataBuffer2D<T>::resize(unsigned int newSizeX, unsigned int newSizeY)
{
    int nbElementsToCopyX = glm::min(newSizeX, mNbElementsX);
    int nbElementsToCopyY = glm::min(newSizeY, mNbElementsY);
//    mNbDataElements = newSizeX*newSizeY;

    if(mbHasCpuStorage) {
        T* newData = new T[newSizeX*newSizeY];
        for(int i=0; i<nbElementsToCopyX; i++) {
            for(int j=0; j<nbElementsToCopyY; j++) {
                //newData[newSizeY*j + i] = data.get()[mSizeY*j + i];
                newData[newSizeY*j + i] = getCpuData(i, j);
            }
//            for(int j=nbElementsToCopyY; j<newSizeY; j++) {
//                newData[newSizeY*j + i] = ;
//            }
        }
//        for(int i=nbElementsToCopyX; i<newSizeX; i++)
//        for(int j=0; j<newSizeY; j++) {
//            newData[newSizeY*j + i] = 0;
//        }
        delete[] dataCpu.get();
        dataCpu.set(newData);
    }

//    if(mbHasBufferStorage) {
//        GLuint newBuffer;
//        glCreateBuffers(1, &newBuffer);
//        glNamedBufferStorage(newBuffer, dataSizeInBytes(), NULL,
//                             GL_DYNAMIC_STORAGE_BIT | GL_MAP_READ_BIT |
//                             GL_MAP_WRITE_BIT);
//        glCopyNamedBufferSubData(mGlidBuffer, newBuffer, 0, 0,
//                                 sizeof(T)*nbElementsToCopyX*nbElementsToCopyY);
//        glDeleteBuffers(1, &mGlidBuffer);
//        mGlidBuffer = newBuffer;
//    }

    if(mbHasTexture2DStorage) {
        GLuint newTexture;
        glCreateTextures(GL_TEXTURE_2D, 1, &newTexture);
        unsigned int nbLevels;
        // this computes nbLevels = log_2(mSize), which is the maximum number of
        // mipmap level we can have for this size.
        nbLevels = 1;
        unsigned int tempSize = glm::max(newSizeX, newSizeY);
        while (tempSize >>= 1) ++nbLevels;
        glTextureStorage2D(newTexture, nbLevels, mTexture2DSizedInternalFormat,
                           newSizeX, newSizeY);

        glClearTexImage(newTexture, 0, mTexture2DExternalFormat,
                        mTexture2DSizedExternalFormat, NULL);
        glCopyImageSubData(mGlidTexture2D, GL_TEXTURE_2D, 0, 0, 0, 0,
                           newTexture, GL_TEXTURE_2D, 0, 0, 0, 0,
                           nbElementsToCopyX, nbElementsToCopyY, 1);
        glDeleteTextures(1, &mGlidTexture2D);
        mGlidTexture2D = newTexture;
//        glGenerateTextureMipmap(mGlidTexture2D);
    }

    mNbElementsX = newSizeX;
    mNbElementsY = newSizeY;
}

template<class T>
void DataBuffer2D<T>::bindBufferToTarget(GLenum target)
{
    if(mbHasBufferStorage) {
        glBindBuffer(target, mGlidBuffer);
    } else {
        Application::sInsertDebugMessage(
                "goglu::DataBuffer2D::bindToTarget : DataBuffer has no "
                "BUFFER storage, cannot bind it to a target.",
                GL_DEBUG_SEVERITY_HIGH);
    }
}


template <class T>
T goglu::DataBuffer2D<T>::getCpuData(unsigned int i, unsigned int j)
{
    dataCpu.refresh();
    return dataCpu.get()[mNbElementsX*j + i];
}


template <class T>
T goglu::DataBuffer2D<T>::getCpuData_noRefresh(unsigned int i, unsigned int j) const
{
    return dataCpu.get_noRefresh()[mNbElementsX*j + i];
}

template <class T>
T* goglu::DataBuffer2D<T>::getCpuDataPointer()
{
    return dataCpu.get();
}


template <class T>
void goglu::DataBuffer2D<T>::refreshCpuData()
{
    dataCpu.refresh();
}


template <class T>
void goglu::DataBuffer2D<T>::setCpuData(unsigned int i, unsigned int j, T inData)
{
    dataCpu.refresh();
    sourceStorageType = StorageType::CPU;
    dataCpu.get()[mNbElementsX*j + i] = inData;
    if(!dataCpu.getIsDirty()) { // WARNING: is this correct or should I dirty all the time?
        dataCpu.dirtyItselfAndDependents();
    }
}

template <class T>
void DataBuffer2D<T>::addCpuData(unsigned int i, unsigned int j, T inData)
{
    dataCpu.refresh();
    sourceStorageType = StorageType::CPU;
    dataCpu.get()[mNbElementsX*j + i] += inData;
    if(!dataCpu.getIsDirty()) { // WARNING: is this correct or should I dirty all the time?
        dataCpu.dirtyItselfAndDependents();
    }
}


template <class T>
void DataBuffer2D<T>::dirtyData()
{
    dataCpu.dirtyItselfAndDependents();
    //dataTexture2D.dirtyItselfAndDependents();
}



//==============================================================================
// PLUGS ACTIONS
//==============================================================================

//template <class T>
//void DataBuffer2D<T>::InDataCpu_dirtyDependents() {
//    if(mbHasCpuStorage) {
////        dataCpu.dirtyItselfAndDependents();
//    }
//    if(mbHasTexture2DStorage) {
////        dataTexture2D.dirtyItselfAndDependents();
//    }
////    if(mbHasBufferStorage) {
////        dataBuffer.dirtyItselfAndDependents();
////    }
//}

//template <class T>
//void DataBuffer2D<T>::InDataCpu_receive(T *data) {
//    setFromCpuData(data);
//}

//template <class T>
//T* DataBuffer2D<T>::OutDataCpu_getData() {
//    return dataCpu.get();
//}

template <class T>
Metadata2DTexture2D DataBuffer2D<T>::out_metadataTexture2D_getData() {
    dataTexture2D.refresh();
    return metadataTexture2D;
}
template <class T>
void DataBuffer2D<T>::out_metadataTexture2D_receivePushConnectionUpdate(ConnectionAbstract* /*connection*/) {}


template <class T>
void DataBuffer2D<T>::in_metadataCpu_dirtyDependents()
{
    dataCpu.dirtyItselfAndDependents();
    dataTexture2D.dirtyItselfAndDependents();
}
template <class T>
void DataBuffer2D<T>::in_metadataCpu_receive(Metadata2DCpu /*data*/)
{
    // TODO
}
template <class T>
void DataBuffer2D<T>::in_metadataCpu_receivePullConnectionUpdate(ConnectionAbstract* /*connection*/) {}


template <class T>
void DataBuffer2D<T>::in_metadataImage2D_dirtyDependents()
{
    dataTexture2D.dirtyItselfAndDependents();
}
template <class T>
void DataBuffer2D<T>::in_metadataImage2D_receive(Metadata2DImage2D data)
{
    sourceStorageType = StorageType::TEXTURE2D;
    dataCpu.dirtyItselfAndDependents();
    metadataImage2D = data;
}
template <class T>
void DataBuffer2D<T>::in_metadataImage2D_receivePullConnectionUpdate(ConnectionAbstract* /*connection*/) {}


template <class T>
Metadata2DImage2D DataBuffer2D<T>::out_metadataImage2D_getData(unsigned int level) {
//    cout << "OutmetadataTexture2D_getData called." << endl;
    dataTexture2D.refresh();
    Metadata2DImage2D metadata;
    metadata.externalFormat = metadataTexture2D.externalFormat;
    metadata.internalFormat = metadataTexture2D.internalFormat;
    metadata.nbElementsX = metadataTexture2D.nbElementsX;
    metadata.nbElementsY = metadataTexture2D.nbElementsY;
    metadata.sizedExternalFormat = metadataTexture2D.sizedExternalFormat;
    metadata.sizedInternalFormat = metadataTexture2D.sizedInternalFormat;
    metadata.textureId = metadataTexture2D.textureId;
    metadata.level = level;
    return metadata;
}




//==============================================================================
// Dirtyable actions
//==============================================================================

template <class T>
void DataBuffer2D<T>::dataCpu_dirtyDependents() {
    dataTexture2D.dirtyItselfAndDependents();
//    out_metadataTexture2D.dirtyConnections();
}

template <class T>
void DataBuffer2D<T>::dataCpu_updateData() {

    switch(sourceStorageType) {
    case StorageType::CPU:
        if(in_metadataCpu.hasPullConnection()) {
            //metadata2DCpu metadataIn = inmetadataCpu.pullUpdate();
            dataCpu.set((T*)(in_metadataCpu.pullUpdate().dataPointer));
        } else {
            // nothing?
        }
        break;
    case StorageType::TEXTURE2D:
        // copy data from buffer to cpu
//        {
//        T* bufferDataCopy = new T[metadataBuffer.nbElements];
//        for(int i=0; i<metadataBuffer.nbElements; i++) {
//            setCpuData(i, bufferDataCopy[i]);
//            cout << "buffer to CPU : " << metadataBuffer.nbElements << endl;
//        }
//        }

//        cout << "Texture2D -> CPU BEGIN" << endl;

//cout << "texture -> cpu" << endl;

        //glGetNamedBufferSubData(metadataBuffer.bufferId, 0, dataCpuSizeInBytes(), dataCpu.get());
        glGetTextureImage(metadataTexture2D.textureId, 0, metadataTexture2D.externalFormat,
                          metadataTexture2D.sizedExternalFormat, dataSizeInBytes(), dataCpu.get());

//        cout << "textureid : " << metadataTexture2D.textureId << endl;
//        cout << "externalFormat : " << metadataTexture2D.externalFormat << " <> " << GL_RED_INTEGER << " NOT " << GL_FLOAT << endl;
//        cout << "sizedExternalFormat : " << metadataTexture2D.sizedExternalFormat << " <> " << GL_INT << " NOT " << GL_RGB << endl;
//        cout << "dataSizeInBytes : " << dataSizeInBytes() << endl;
//        cout << "data[0] = " << ((int*)dataCpu.get())[0] << endl;

//        cout << "Texture2D -> CPU END" << endl;

//        cout << "data from texture to cpu, " << typeid(dataCpu.get()).name() << endl;

//        sourceStorageType = StorageType::CPU;
//        dataCpu.dirtyItselfAndDependents();
//        dataBuffer.dirtyItselfAndDependents();
        break;
    default:
        cout << "DataBuffer2D<T>::Dirt_dataCpu_updateData() : unknown source storage type." << endl;
    }


}


template <class T>
void DataBuffer2D<T>::dataTexture2D_dirtyDependents() {
    dataCpu.dirtyItselfAndDependents();
    out_metadataTexture2D.dirtyConnections();
    out_metadataImage2D.dirtyConnections();
}

template <class T>
void DataBuffer2D<T>::dataTexture2D_updateData() {
//    dataCpu.refresh();
////    cout << "transferring texture to GPU" << endl;
//    glTextureSubImage2D(mGlidTexture2D, 0, 0, 0, mNbElementsX, mNbElementsY,
//                        mTexture2DExternalFormat,
//                        mTexture2DSizedExternalFormat, dataCpu.get());
//    glGenerateTextureMipmap(mGlidTexture2D);

    switch(sourceStorageType) {
    case StorageType::CPU:
//        cout << "CPU" << endl;
        dataCpu.refresh();
        glTextureSubImage2D(mGlidTexture2D, 0, 0, 0, mNbElementsX, mNbElementsY,
                            mTexture2DExternalFormat,
                            mTexture2DSizedExternalFormat, dataCpu.get());
        glGenerateTextureMipmap(mGlidTexture2D);
        break;
    case StorageType::TEXTURE2D:
//        cout << "BUFFER" << endl;
        // TODO: check if pull connection for the buffer input plug.
        break;
    default:
        cout << "Dirt_dataBuffer_updateData : unknown source storage type." << endl;
    }
}








