// TODO: load different formats, right now it can only load BMP (maybe other
// formats too, but not PNG.)
//

#include "DataBuffer2D.h" // not necessary, but helps IDEs.

//#include "../Application.h"

template<class T>
DataBuffer2D<T>::DataBuffer2D(unsigned int sizeX, unsigned int sizeY)
    //dataCpu(*this),
    //dataTexture2D(*this)
//    dataBuffer(*this)
{
    _nbElementsX = sizeX;
    _nbElementsY = sizeY;
//    mNbDataElements = _nbElementsX*_nbElementsY;
    _hasCpuStorage = false;
//    mbHasBufferStorage = false;
    _hasTexture2DStorage = false;

}

template<class T>
void DataBuffer2D<T>::createCpuStorage()
{
    _dataCpu = new T[_nbElementsX*_nbElementsY];
    _hasCpuStorage = true;
}

template<class T>
void DataBuffer2D<T>::deleteCpuStorage()
{
    delete[] mpData;
    _hasCpuStorage = false;
}

template<class T>
void DataBuffer2D<T>::createBufferStorage()
{
    glCreateBuffers(1, &_glidBuffer);
    glNamedBufferStorage(_glidBuffer, dataSizeInBytes(), NULL,
                            GL_DYNAMIC_STORAGE_BIT | GL_MAP_READ_BIT |
                            GL_MAP_WRITE_BIT);
    _hasBufferStorage = true;
}

template<class T>
void DataBuffer2D<T>::deleteBufferStorage()
{
    glDeleteBuffers(1, &_glidBuffer);
    mbHasBufferStorage = false;
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
    unsigned int tempSize = glm::max(_nbElementsX, _nbElementsY);
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
    _texture2DInternalFormat = internalFormat;
    _texture2DSizedInternalFormat = sizedInternalFormat;
    _texture2DExternalFormat = externalFormat;
    _texture2DSizedExternalFormat = sizedExternalFormat;
    glCreateTextures(GL_TEXTURE_2D, 1, &_glidTexture2D);
//        unsigned int nbLevels;
//        // this computes nbLevels = log_2(max(mSizeX,mSizeY)), which is the maximum number of
//        // mipmap level we can have for this size.
//        nbLevels = 1;
//        unsigned int tempSize = max(_nbElementsX, _nbElementsY);
//        while (tempSize >>= 1) ++nbLevels;

    glTextureStorage2D(_glidTexture2D, nbMipmapLevels,
                        _texture2DSizedInternalFormat, _nbElementsX, _nbElementsY);
    glClearTexImage(_glidTexture2D, 0, _texture2DExternalFormat,
                    _texture2DSizedExternalFormat, NULL);
    glGenerateTextureMipmap(_glidTexture2D);
    glTextureParameteri(_glidTexture2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTextureParameteri(_glidTexture2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    // clamp to zero by default, override if different boundary conditions desired.
    glTextureParameteri(_glidTexture2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTextureParameteri(_glidTexture2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    float borderColor[] = { 0.0f, 0.0f, 0.0f, 1.0f };
    glTextureParameterfv(_glidTexture2D, GL_TEXTURE_BORDER_COLOR, borderColor);
    _hasTexture2DStorage = true;

    // metadata
    _metadataTexture2D.textureId = _glidTexture2D;
    _metadataTexture2D.nbElementsX = _nbElementsX;
    _metadataTexture2D.nbElementsY = _nbElementsY;
    _metadataTexture2D.internalFormat = _texture2DInternalFormat;
    _metadataTexture2D.sizedInternalFormat = _texture2DSizedInternalFormat;
    _metadataTexture2D.externalFormat = _texture2DExternalFormat;
    _metadataTexture2D.sizedExternalFormat = _texture2DSizedExternalFormat;

}

template<class T>
void DataBuffer2D<T>::deleteTexture2DStorage()
{
    glDeleteTextures(1, &_glidTexture2D);
    _hasTexture1DStorage = false;
}

//template<class T>
//void DataBuffer2D<T>::transferCpuDataToGpu()
//{
//    if(_hasTexture2DStorage) {
//        glTextureSubImage2D(_glidTexture2D, 0, 0, 0, mSizeX, mSizeY,
//                            _texture2DExternalFormat,
//                            _texture2DSizedExternalFormat, mpData);
//        glGenerateTextureMipmap(_glidTexture2D);
//    }
//    if(mbHasBufferStorage) {
//        glNamedBufferSubData(_glidBuffer, 0, dataSizeInBytes(), mpData);
//    }
//}

//template<class T>
//void DataBuffer2D<T>::transferBufferDataToCpu()
//{
//    if(_hasCpuStorage && mbHasBufferStorage) {
//        glGetNamedBufferSubData( _glidBuffer, 0, dataSizeInBytes(), mpData);
//    } else {
//        if(!_hasCpuStorage) {
//            Application::sInsertDebugMessage(
//                    "DataBuffer2D::transferBufferDataToCpu : DataBuffer "
//                    "has no CPU storage, cannot transfer BUFFER data to it.",
//                    GL_DEBUG_SEVERITY_HIGH);
//        } else {
//            Application::sInsertDebugMessage(
//                    "DataBuffer2D::transferBufferDataToCpu : DataBuffer "
//                    "has no BUFFER storage, cannot transfer it to CPU.",
//                    GL_DEBUG_SEVERITY_HIGH);
//        }
//    }
//}

//template<class T>
//void DataBuffer2D<T>::transferTexture2DDataToCpu()
//{
//    if(_hasCpuStorage && _hasTexture2DStorage) {
//        glGetTextureImage(_glidTexture2D, 0, _texture2DExternalFormat,
//                          _texture2DSizedExternalFormat, dataSizeInBytes(),
//                          mpData);
//    } else {
//        if(!_hasCpuStorage) {
//            Application::sInsertDebugMessage(
//                    "DataBuffer2D::transferTexture2DToCpu : DataBuffer "
//                    "has no CPU storage, cannot transfer TEXTURE_2D data to it.",
//                    GL_DEBUG_SEVERITY_HIGH);
//        } else {
//            Application::sInsertDebugMessage(
//                    "DataBuffer2D::transferTexture2DDataToCpu : "
//                    "DataBuffer has no TEXTURE_2D storage, cannot transfer it "
//                    "to CPU.",
//                    GL_DEBUG_SEVERITY_HIGH);
//        }
//    }
//}


// assumes the size of data is the same as the DataBuffer's size.
template<class T>
void DataBuffer2D<T>::setFromCpuData(T* inData) {

    if(_hasCpuStorage) {
        for(unsigned int i=0; i<_nbElementsX; i++) {
            for(unsigned int j=0; j<_nbElementsY; j++) {
                setCpuData(i, j, inData[_nbElementsY*j + i]);
            }
        }
    }

    if(_hasTexture2DStorage) {
        glTextureSubImage2D(_glidTexture2D, 0, 0, 0, _nbElementsX, _nbElementsY,
                            _texture2DExternalFormat,
                            _texture2DSizedExternalFormat, inData);
//        glGenerateTextureMipmap(_glidTexture2D);
    }

//    if(mbHasBufferStorage) {
//        glNamedBufferSubData(_glidBuffer, 0, dataSizeInBytes(), inData);
//    }

}

// assumes the size of data is the same as the DataBuffer's size.
template<class T>
void DataBuffer2D<T>::setFromImage(Metadata2DImage2D srcImage, unsigned int level) {

    if(_hasTexture2DStorage) {
        glCopyImageSubData(
                srcImage.textureId, GL_TEXTURE_2D, srcImage.level,
                0, 0, 0,
                _metadataTexture2D.textureId, GL_TEXTURE_2D, level,
                0, 0, 0,
                _metadataTexture2D.nbElementsX, _metadataTexture2D.nbElementsY, 1);
        _sourceStorageType = StorageType::TEXTURE2D;
    } else {
        // TODO: error, texture storage required.
    }


}

template <class T>
T DataBuffer2D<T>::interp(glm::vec2 pos)
{
    float x = pos.x;
    float y = pos.y;

    int indexXLeft, indexXRight;
    int indexYLeft, indexYRight;
    float weightX;
    float weightY;

    // TODO: this is not clean...
    float _boundXMin = 0;
    float _boundYMin = 0;
    float _boundXMax = 1;
    float _boundYMax = 1;
    int mNbCellsX = _nbElementsX;
    int mNbCellsY = _nbElementsY;

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
    int nbElementsToCopyX = glm::min(newSizeX, _nbElementsX);
    int nbElementsToCopyY = glm::min(newSizeY, _nbElementsY);
//    mNbDataElements = newSizeX*newSizeY;

    if(_hasCpuStorage) {
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
        delete[] _dataCpu;
        _dataCpu = newData;
    }

//    if(mbHasBufferStorage) {
//        GLuint newBuffer;
//        glCreateBuffers(1, &newBuffer);
//        glNamedBufferStorage(newBuffer, dataSizeInBytes(), NULL,
//                             GL_DYNAMIC_STORAGE_BIT | GL_MAP_READ_BIT |
//                             GL_MAP_WRITE_BIT);
//        glCopyNamedBufferSubData(_glidBuffer, newBuffer, 0, 0,
//                                 sizeof(T)*nbElementsToCopyX*nbElementsToCopyY);
//        glDeleteBuffers(1, &_glidBuffer);
//        _glidBuffer = newBuffer;
//    }

    if(_hasTexture2DStorage) {
        GLuint newTexture;
        glCreateTextures(GL_TEXTURE_2D, 1, &newTexture);
        unsigned int nbLevels;
        // this computes nbLevels = log_2(mSize), which is the maximum number of
        // mipmap level we can have for this size.
        nbLevels = 1;
        unsigned int tempSize = glm::max(newSizeX, newSizeY);
        while (tempSize >>= 1) ++nbLevels;
        glTextureStorage2D(newTexture, nbLevels, _texture2DSizedInternalFormat,
                           newSizeX, newSizeY);

        glClearTexImage(newTexture, 0, _texture2DExternalFormat,
                        _texture2DSizedExternalFormat, NULL);
        glCopyImageSubData(_glidTexture2D, GL_TEXTURE_2D, 0, 0, 0, 0,
                           newTexture, GL_TEXTURE_2D, 0, 0, 0, 0,
                           nbElementsToCopyX, nbElementsToCopyY, 1);
        glDeleteTextures(1, &_glidTexture2D);
        _glidTexture2D = newTexture;
//        glGenerateTextureMipmap(_glidTexture2D);
    }

    _nbElementsX = newSizeX;
    _nbElementsY = newSizeY;
}

template<class T>
void DataBuffer2D<T>::bindBufferToTarget(GLenum target)
{
    glBindBuffer(target, _glidBuffer);
}


template <class T>
T DataBuffer2D<T>::getCpuData(unsigned int i, unsigned int j)
{
    return _dataCpu[_nbElementsX*j + i];
}


template <class T>
T DataBuffer2D<T>::getCpuData_noRefresh(unsigned int i, unsigned int j) const
{
    return _dataCpu[_nbElementsX*j + i];
}

template <class T>
T* DataBuffer2D<T>::getCpuDataPointer()
{
    return _dataCpu;
}


template <class T>
void DataBuffer2D<T>::setCpuData(unsigned int i, unsigned int j, T inData)
{
    _sourceStorageType = StorageType::CPU;
    _dataCpu[_nbElementsX*j + i] = inData;
}

template <class T>
void DataBuffer2D<T>::addCpuData(unsigned int i, unsigned int j, T inData)
{
    _sourceStorageType = StorageType::CPU;
    _dataCpu[_nbElementsX*j + i] += inData;
}



//
//
////==============================================================================
//// PLUGS ACTIONS
////==============================================================================
//
////template <class T>
////void DataBuffer2D<T>::InDataCpu_dirtyDependents() {
////    if(_hasCpuStorage) {
//////        dataCpu.dirtyItselfAndDependents();
////    }
////    if(_hasTexture2DStorage) {
//////        dataTexture2D.dirtyItselfAndDependents();
////    }
//////    if(mbHasBufferStorage) {
//////        dataBuffer.dirtyItselfAndDependents();
//////    }
////}
//
////template <class T>
////void DataBuffer2D<T>::InDataCpu_receive(T *data) {
////    setFromCpuData(data);
////}
//
////template <class T>
////T* DataBuffer2D<T>::OutDataCpu_getData() {
////    return dataCpu.get();
////}
//
//template <class T>
//Metadata2DTexture2D DataBuffer2D<T>::out_metadataTexture2D_getData() {
//    dataTexture2D.refresh();
//    return metadataTexture2D;
//}
//template <class T>
//void DataBuffer2D<T>::out_metadataTexture2D_receivePushConnectionUpdate(ConnectionAbstract* /*connection*/) {}
//
//
//template <class T>
//void DataBuffer2D<T>::in_metadataCpu_dirtyDependents()
//{
//    dataCpu.dirtyItselfAndDependents();
//    dataTexture2D.dirtyItselfAndDependents();
//}
//template <class T>
//void DataBuffer2D<T>::in_metadataCpu_receive(Metadata2DCpu /*data*/)
//{
//    // TODO
//}
//template <class T>
//void DataBuffer2D<T>::in_metadataCpu_receivePullConnectionUpdate(ConnectionAbstract* /*connection*/) {}
//
//
//template <class T>
//void DataBuffer2D<T>::in_metadataImage2D_dirtyDependents()
//{
//    dataTexture2D.dirtyItselfAndDependents();
//}
//template <class T>
//void DataBuffer2D<T>::in_metadataImage2D_receive(Metadata2DImage2D data)
//{
//    sourceStorageType = StorageType::TEXTURE2D;
//    dataCpu.dirtyItselfAndDependents();
//    metadataImage2D = data;
//}
//template <class T>
//void DataBuffer2D<T>::in_metadataImage2D_receivePullConnectionUpdate(ConnectionAbstract* /*connection*/) {}
//
//
//template <class T>
//Metadata2DImage2D DataBuffer2D<T>::out_metadataImage2D_getData(unsigned int level) {
////    cout << "OutmetadataTexture2D_getData called." << endl;
//    dataTexture2D.refresh();
//    Metadata2DImage2D metadata;
//    metadata.externalFormat = metadataTexture2D.externalFormat;
//    metadata.internalFormat = metadataTexture2D.internalFormat;
//    metadata.nbElementsX = metadataTexture2D.nbElementsX;
//    metadata.nbElementsY = metadataTexture2D.nbElementsY;
//    metadata.sizedExternalFormat = metadataTexture2D.sizedExternalFormat;
//    metadata.sizedInternalFormat = metadataTexture2D.sizedInternalFormat;
//    metadata.textureId = metadataTexture2D.textureId;
//    metadata.level = level;
//    return metadata;
//}
//
//
//
//
////==============================================================================
//// Dirtyable actions
////==============================================================================
//
//template <class T>
//void DataBuffer2D<T>::dataCpu_dirtyDependents() {
//    dataTexture2D.dirtyItselfAndDependents();
////    out_metadataTexture2D.dirtyConnections();
//}
//
//template <class T>
//void DataBuffer2D<T>::dataCpu_updateData() {
//
//    switch(sourceStorageType) {
//    case StorageType::CPU:
//        if(in_metadataCpu.hasPullConnection()) {
//            //metadata2DCpu metadataIn = inmetadataCpu.pullUpdate();
//            dataCpu.set((T*)(in_metadataCpu.pullUpdate().dataPointer));
//        } else {
//            // nothing?
//        }
//        break;
//    case StorageType::TEXTURE2D:
//        // copy data from buffer to cpu
////        {
////        T* bufferDataCopy = new T[metadataBuffer.nbElements];
////        for(int i=0; i<metadataBuffer.nbElements; i++) {
////            setCpuData(i, bufferDataCopy[i]);
////            cout << "buffer to CPU : " << metadataBuffer.nbElements << endl;
////        }
////        }
//
////        cout << "Texture2D -> CPU BEGIN" << endl;
//
////cout << "texture -> cpu" << endl;
//
//        //glGetNamedBufferSubData(metadataBuffer.bufferId, 0, dataCpuSizeInBytes(), dataCpu.get());
//        glGetTextureImage(metadataTexture2D.textureId, 0, metadataTexture2D.externalFormat,
//                          metadataTexture2D.sizedExternalFormat, dataSizeInBytes(), dataCpu.get());
//
////        cout << "textureid : " << metadataTexture2D.textureId << endl;
////        cout << "externalFormat : " << metadataTexture2D.externalFormat << " <> " << GL_RED_INTEGER << " NOT " << GL_FLOAT << endl;
////        cout << "sizedExternalFormat : " << metadataTexture2D.sizedExternalFormat << " <> " << GL_INT << " NOT " << GL_RGB << endl;
////        cout << "dataSizeInBytes : " << dataSizeInBytes() << endl;
////        cout << "data[0] = " << ((int*)dataCpu.get())[0] << endl;
//
////        cout << "Texture2D -> CPU END" << endl;
//
////        cout << "data from texture to cpu, " << typeid(dataCpu.get()).name() << endl;
//
////        sourceStorageType = StorageType::CPU;
////        dataCpu.dirtyItselfAndDependents();
////        dataBuffer.dirtyItselfAndDependents();
//        break;
//    default:
//        cout << "DataBuffer2D<T>::Dirt_dataCpu_updateData() : unknown source storage type." << endl;
//    }
//
//
//}
//
//
//template <class T>
//void DataBuffer2D<T>::dataTexture2D_dirtyDependents() {
//    dataCpu.dirtyItselfAndDependents();
//    out_metadataTexture2D.dirtyConnections();
//    out_metadataImage2D.dirtyConnections();
//}
//
//template <class T>
//void DataBuffer2D<T>::dataTexture2D_updateData() {
////    dataCpu.refresh();
//////    cout << "transferring texture to GPU" << endl;
////    glTextureSubImage2D(_glidTexture2D, 0, 0, 0, _nbElementsX, _nbElementsY,
////                        _texture2DExternalFormat,
////                        _texture2DSizedExternalFormat, dataCpu.get());
////    glGenerateTextureMipmap(_glidTexture2D);
//
//    switch(sourceStorageType) {
//    case StorageType::CPU:
////        cout << "CPU" << endl;
//        dataCpu.refresh();
//        glTextureSubImage2D(_glidTexture2D, 0, 0, 0, _nbElementsX, _nbElementsY,
//                            _texture2DExternalFormat,
//                            _texture2DSizedExternalFormat, dataCpu.get());
//        glGenerateTextureMipmap(_glidTexture2D);
//        break;
//    case StorageType::TEXTURE2D:
////        cout << "BUFFER" << endl;
//        // TODO: check if pull connection for the buffer input plug.
//        break;
//    default:
//        cout << "Dirt_dataBuffer_updateData : unknown source storage type." << endl;
//    }
//}
//







