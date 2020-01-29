#include "DataBuffer1D.h"

using namespace goglu;
using namespace std;
using namespace glm;

template<class T>
DataBuffer1D<T>::DataBuffer1D(unsigned int size) :
    in_metadataCpu(*this),
//    outDataCpu(*this),
//    outDataBuffer(*this),
    in_metadataBuffer(*this),
    out_metadataBuffer(*this),
    dataCpu(*this),
    dataBuffer(*this),
    out_nbParticles(*this)
{
//    if(size == 0) {
//        mNbElements = 0;
//        capacity = 1;
//    } else {
        mNbElements = size;
        capacity = size;
//    }
    mbHasCpuStorage = false;
    mbHasBufferStorage = false;
    mbHasTexture1DStorage = false;
    mbHasTextureBufferStorage = false;
    sourceStorageType = StorageType::UNKNOWN;
}

//template<class T>
//T& goglu::DataBuffer1D<T>::operator()(unsigned int x)
//{
//    if(mbHasCpuStorage) {
//        return mpData[x];
//    } else {
//        Application::sInsertDebugMessage(
//                "goglu::DataBuffer1D::operator() : DataBuffer has no "
//                "CPU storage, accessing it returns a default value.",
//                GL_DEBUG_SEVERITY_HIGH);
//        T temp;
//        return temp;
//    }
//}

//template <class T>
//T goglu::DataBuffer1D<T>::getCpuData(unsigned int i)
//{
//    return dataCpu.get()[i];
//}

template <class T>
T goglu::DataBuffer1D<T>::getCpuData(unsigned int i)
{
    return dataCpu.get()[i];
}

template <class T>
T goglu::DataBuffer1D<T>::getCpuData_noRefresh(unsigned int i) const
{
    return dataCpu.get_noRefresh()[i];
}

template <class T>
T* goglu::DataBuffer1D<T>::getCpuDataPointer()
{
    return dataCpu.get();
}


template <class T>
void goglu::DataBuffer1D<T>::refreshCpuData()
{
    dataCpu.refresh();
}


template <class T>
void goglu::DataBuffer1D<T>::setCpuData(unsigned int i, T data)
{
    sourceStorageType = StorageType::CPU;
    dataCpu.get()[i] = data;
    dataCpu.dirtyItselfAndDependents();
    dataBuffer.dirtyItselfAndDependents();
}

template <class T>
void goglu::DataBuffer1D<T>::setCpuData_noDirty(unsigned int i, T data)
{
    sourceStorageType = StorageType::CPU;
    dataCpu.get_noRefresh()[i] = data;
}

template <class T>
void goglu::DataBuffer1D<T>::dirtyData()
{
    dataCpu.dirtyItselfAndDependents();
    dataBuffer.dirtyItselfAndDependents();
}

template<class T>
void DataBuffer1D<T>::createCpuStorage()
{
    if(mbHasCpuStorage) {
        Application::sInsertDebugMessage(
                "goglu::DataBuffer1D::createCpuStorage : DataBuffer already has "
                "CPU storage, creating the storage again has no effect.",
                GL_DEBUG_SEVERITY_LOW);
    } else {
        dataCpu.set(new T[mNbElements]);
        mbHasCpuStorage = true;

        // metadata
        metadataCpu.dataPointer = dataCpu.get();
    }
}

template<class T>
void DataBuffer1D<T>::deleteCpuStorage()
{
    if(mbHasCpuStorage) {
        delete[] dataCpu.get();
        mbHasCpuStorage = false;
    } else {
        Application::sInsertDebugMessage(
                "goglu::DataBuffer1D::deleteCpuStorage : DataBuffer has no CPU "
                "storage, deleting the storage has no effect.",
                GL_DEBUG_SEVERITY_LOW);
    }
}

template<class T>
void DataBuffer1D<T>::createBufferStorage(GLenum dataType, unsigned int nbElementsPerComponent, GLbitfield flags)
{
    if(mbHasBufferStorage) {
        Application::sInsertDebugMessage(
                "goglu::DataBuffer1D::createBufferStorage : DataBuffer already "
                "has BUFFER storage, creating the storage again has no effect.",
                GL_DEBUG_SEVERITY_LOW);
    } else {
        glCreateBuffers(1, &mGlidBuffer);
        glNamedBufferStorage(mGlidBuffer, mNbElements * nbElementsPerComponent * sizeOfEnumType(dataType), NULL, 
                             flags);
        mbHasBufferStorage = true;

        // metadata
        metadataBuffer.bufferId = mGlidBuffer;
        metadataBuffer.dataType = dataType;
        metadataBuffer.nbElements = mNbElements;
        metadataBuffer.nbElementsPerComponent = nbElementsPerComponent;
    }
}

template<class T>
void DataBuffer1D<T>::deleteBufferStorage()
{
    if(mbHasBufferStorage) {
        glDeleteBuffers(1, &mGlidBuffer);
        mbHasBufferStorage = false;
    } else {
        Application::sInsertDebugMessage(
                "goglu::DataBuffer1D::deleteBufferStorage : DataBuffer has no "
                "BUFFER storage, deleting the storage has no effect.",
                GL_DEBUG_SEVERITY_LOW);
    }
}

//nbMipmapLevels=0 will put the maximum number of mipmaps. nbMipmapLevels=1
//will only put the top layer without any mipmaps below it.
template<class T>
void DataBuffer1D<T>::createTexture1DStorage(GLenum internalFormat,
                                             GLenum sizedInternalFormat,
                                             GLenum externalFormat,
                                             GLenum sizedExternalFormat)
{
    if(mbHasTexture1DStorage) {
        Application::sInsertDebugMessage(
                "goglu::DataBuffer1D::createTexture1DStorage : DataBuffer "
                "already has TEXTURE_1D storage, creating the storage again "
                "has no effect.",
                GL_DEBUG_SEVERITY_LOW);
    } else {
        mTexture1DInternalFormat = internalFormat;
        mTexture1DSizedInternalFormat = sizedInternalFormat;
        mTexture1DExternalFormat = externalFormat;
        mTexture1DSizedExternalFormat = sizedExternalFormat;
        glCreateTextures(GL_TEXTURE_1D, 1, &mGlidTexture1D);
        unsigned int nbLevels;
        // this computes nbLevels = log_2(mSize), which is the maximum number of
        // mipmap level we can have for this size.
        nbLevels = 1;
        unsigned int tempSize = mNbElements;
        while (tempSize >>= 1) ++nbLevels;

        glTextureStorage1D(mGlidTexture1D, nbLevels,
                           mTexture1DSizedInternalFormat, mNbElements);
        glClearTexImage(mGlidTexture1D, 0, mTexture1DExternalFormat,
                        mTexture1DSizedExternalFormat, NULL);
        glGenerateTextureMipmap(mGlidTexture1D);
        mbHasTexture1DStorage = true;
    }
}



template<class T>
void DataBuffer1D<T>::deleteTexture1DStorage()
{
    if(mbHasTexture1DStorage) {
        glDeleteTextures(1, &mGlidTexture1D);
        mbHasTexture1DStorage = false;
    } else {
        Application::sInsertDebugMessage(
                "goglu::DataBuffer1D::deleteTexture1DSorage : DataBuffer has "
                "no TEXTURE_1D storage, deleting the storage has no effect.",
                GL_DEBUG_SEVERITY_LOW);
    }
}


template<class T>
void DataBuffer1D<T>::createTextureBufferStorage(GLenum internalSizedFormat)
{
    if(mbHasTextureBufferStorage) {
        Application::sInsertDebugMessage(
                "goglu::DataBuffer1D::createTextureBufferStorage : DataBuffer "
                "already has TEXTURE_BUFFER storage, creating the storage again "
                "has no effect.",
                GL_DEBUG_SEVERITY_LOW);
    } else {
        if(!mbHasBufferStorage) {
            createBufferStorage();
        }
        mTextureBufferSizedInternalFormat = internalSizedFormat;
        glCreateTextures(GL_TEXTURE_BUFFER, 1, &mGlidTextureBuffer);
        glTextureBuffer(mGlidTextureBuffer, mTextureBufferSizedInternalFormat,
                        mGlidBuffer);
        mbHasTextureBufferStorage = true;
    }
}


template<class T>
void DataBuffer1D<T>::deleteTextureBufferStorage()
{
    if(mbHasTextureBufferStorage) {
        glDeleteTextures(1, &mGlidTextureBuffer);
        mbHasTextureBufferStorage = false;
    } else {
        Application::sInsertDebugMessage(
                "goglu::DataBuffer1D::deleteTextureBufferStorage : DataBuffer "
                "has no TEXTURE_BUFFER storage, deleting the storage has no effect.",
                GL_DEBUG_SEVERITY_LOW);
    }
}


//template<class T>
//void DataBuffer1D<T>::transferCpuDataToGpu()
//{
//    if(mbHasTexture1DStorage) {
//        glTextureSubImage1D(mGlidTexture1D, 0, 0, mNbDataElements,
//                            mTexture1DExternalFormat,
//                            mTexture1DSizedExternalFormat, mpData);
//        glGenerateTextureMipmap(mGlidTexture1D);
//    }
//    if(mbHasBufferStorage||mbHasTextureBufferStorage) {
//        glNamedBufferSubData(mGlidBuffer, 0, dataSizeInBytes(), mpData);
//    }
//}


//template<class T>
//void DataBuffer1D<T>::transferBufferDataToCpu()
//{
//    if(mbHasCpuStorage && mbHasBufferStorage) {
//        glGetNamedBufferSubData( mGlidBuffer, 0, dataSizeInBytes(), mpData);
//    } else {
//        if(!mbHasCpuStorage) {
//            Application::sInsertDebugMessage(
//                    "goglu::DataBuffer1D::transferBufferDataToCpu : DataBuffer "
//                    "has no CPU storage, cannot transfer BUFFER data to it.",
//                    GL_DEBUG_SEVERITY_HIGH);
//        } else {
//            Application::sInsertDebugMessage(
//                    "goglu::DataBuffer1D::transferBufferDataToCpu : DataBuffer "
//                    "has no BUFFER storage, cannot transfer it to CPU.",
//                    GL_DEBUG_SEVERITY_HIGH);
//        }
//    }
//}


//template<class T>
//void DataBuffer1D<T>::transferTexture1DDataToCpu()
//{
//    if(mbHasCpuStorage && mbHasTexture1DStorage) {
//        glGetTextureImage(mGlidTexture1D, 0, mTexture1DExternalFormat,
//                          mTexture1DSizedExternalFormat, dataSizeInBytes(),
//                          mpData);
//    } else {
//        if(!mbHasCpuStorage) {
//            Application::sInsertDebugMessage(
//                    "goglu::DataBuffer1D::transferTexture1DToCpu : DataBuffer "
//                    "has no CPU storage, cannot transfer TEXTURE_1D data to it.",
//                    GL_DEBUG_SEVERITY_HIGH);
//        } else {
//            Application::sInsertDebugMessage(
//                    "goglu::DataBuffer1D::transferTexture1DDataToCpu : "
//                    "DataBuffer has no TEXTURE_1D storage, cannot transfer it "
//                    "to CPU.",
//                    GL_DEBUG_SEVERITY_HIGH);
//        }
//    }
//}


template<class T>
void DataBuffer1D<T>::resize(unsigned int size)
{
    int nbElementsToCopy = glm::min(size, mNbElements);
    mNbElements = size;

    if(size <= capacity) {
        if(mbHasBufferStorage) {
            metadataBuffer.nbElements = mNbElements;
        }
    } else {

//        if(mbHasCpuStorage) {
//            T* newData = new T[size];
//            for(int i=0; i<nbElementsToCopy; i++) {
//                newData[i] = getCpuData(i);
//            }
//            for(int i=nbElementsToCopy; i<size; i++) {
//    //            newData[i] = 0;
//                // initialize? With what?
//            }
//            delete[] dataCpu.get();
//            dataCpu.set(newData);
//            metadataCpu.dataPointer = dataCpu.get();
//        }

//        if(mbHasBufferStorage) {
//            GLuint newBuffer;
//            glCreateBuffers(1, &newBuffer);
//            //glNamedBufferStorage(newBuffer, dataCpuSizeInBytes(), NULL,
//            glNamedBufferStorage(newBuffer, mNbElements * metadataBuffer.nbElementsPerComponent * sizeOfEnumType(metadataBuffer.dataType), NULL,
//                                 GL_DYNAMIC_STORAGE_BIT | GL_MAP_READ_BIT |
//                                 GL_MAP_WRITE_BIT);
//            glCopyNamedBufferSubData(mGlidBuffer, newBuffer, 0, 0,
//                                     sizeof(T)*nbElementsToCopy);
//            glDeleteBuffers(1, &mGlidBuffer);
//            mGlidBuffer = newBuffer;

//            metadataBuffer.bufferId = newBuffer;
//    //        metadataBuffer.dataType = dataType;
//            metadataBuffer.nbElements = mNbElements;
//    //        metadataBuffer.nbElementsPerComponent = nbElementsPerComponent;


//        }

//        if(mbHasTexture1DStorage) {
//            GLuint newTexture;
//            glCreateTextures(GL_TEXTURE_1D, 1, &newTexture);
//            unsigned int nbLevels;
//            // this computes nbLevels = log_2(mSize), which is the maximum number of
//            // mipmap level we can have for this size.
//            nbLevels = 1;
//            unsigned int tempSize = mNbElements;
//            while (tempSize >>= 1) ++nbLevels;
//            glTextureStorage1D(newTexture, nbLevels, mTexture1DSizedInternalFormat,
//                               mNbElements);
//            glClearTexImage(newTexture, 0, mTexture1DExternalFormat,
//                            mTexture1DSizedExternalFormat, NULL);
//            glCopyImageSubData(mGlidTexture1D, GL_TEXTURE_1D, 0, 0, 0, 0,
//                               newTexture, GL_TEXTURE_1D, 0, 0, 0, 0,
//                               nbElementsToCopy, 1, 1);
//            glDeleteTextures(1, &mGlidTexture1D);
//            mGlidTexture1D = newTexture;
//            glGenerateTextureMipmap(mGlidTexture1D);
//        }

//        if(mbHasTextureBufferStorage) {
//            glDeleteTextures(1, &mGlidTextureBuffer);
//            glCreateTextures(GL_TEXTURE_BUFFER, 1, &mGlidTextureBuffer);
//            glTextureBuffer(mGlidTextureBuffer, mTextureBufferSizedInternalFormat,
//                            mGlidBuffer);
//        }

//        outmetadataBuffer.dirtyConnections();
//        dataCpu.dirtyItselfAndDependents();
//        dataBuffer.dirtyItselfAndDependents();

        while(capacity < size) {
            if(capacity == 0) {
                capacity = 1;
            } else {
                capacity *= 2;
            }
        }

        if(mbHasCpuStorage) {
            T* newData = new T[capacity];
            for(int i=0; i<nbElementsToCopy; i++) {
                newData[i] = getCpuData(i);
            }
//            for(unsigned int i=nbElementsToCopy; i<size; i++) {
    //            newData[i] = 0;
                // initialize? With what?
//            }
            delete[] dataCpu.get();
            dataCpu.set(newData);
            metadataCpu.dataPointer = dataCpu.get();
        }

        if(mbHasBufferStorage) {
            GLuint newBuffer;
            glCreateBuffers(1, &newBuffer);
            //glNamedBufferStorage(newBuffer, dataCpuSizeInBytes(), NULL,
            glNamedBufferStorage(newBuffer, capacity * metadataBuffer.nbElementsPerComponent * sizeOfEnumType(metadataBuffer.dataType), NULL,
                                 GL_DYNAMIC_STORAGE_BIT | GL_MAP_READ_BIT |
                                 GL_MAP_WRITE_BIT);
            glCopyNamedBufferSubData(mGlidBuffer, newBuffer, 0, 0,
                                     sizeof(T)*nbElementsToCopy);
            glDeleteBuffers(1, &mGlidBuffer);
            mGlidBuffer = newBuffer;

            metadataBuffer.bufferId = newBuffer;
    //        metadataBuffer.dataType = dataType;
            metadataBuffer.nbElements = mNbElements;
    //        metadataBuffer.nbElementsPerComponent = nbElementsPerComponent;


        }

        if(mbHasTexture1DStorage) {
            GLuint newTexture;
            glCreateTextures(GL_TEXTURE_1D, 1, &newTexture);
            unsigned int nbLevels;
            // this computes nbLevels = log_2(mSize), which is the maximum number of
            // mipmap level we can have for this size.
            nbLevels = 1;
            unsigned int tempSize = mNbElements;
            while (tempSize >>= 1) ++nbLevels;
            glTextureStorage1D(newTexture, nbLevels, mTexture1DSizedInternalFormat,
                               mNbElements);
            glClearTexImage(newTexture, 0, mTexture1DExternalFormat,
                            mTexture1DSizedExternalFormat, NULL);
            glCopyImageSubData(mGlidTexture1D, GL_TEXTURE_1D, 0, 0, 0, 0,
                               newTexture, GL_TEXTURE_1D, 0, 0, 0, 0,
                               nbElementsToCopy, 1, 1);
            glDeleteTextures(1, &mGlidTexture1D);
            mGlidTexture1D = newTexture;
            glGenerateTextureMipmap(mGlidTexture1D);
        }

        if(mbHasTextureBufferStorage) {
            glDeleteTextures(1, &mGlidTextureBuffer);
            glCreateTextures(GL_TEXTURE_BUFFER, 1, &mGlidTextureBuffer);
            glTextureBuffer(mGlidTextureBuffer, mTextureBufferSizedInternalFormat,
                            mGlidBuffer);
        }

        out_metadataBuffer.dirtyConnections();
        dataCpu.dirtyItselfAndDependents();
        dataBuffer.dirtyItselfAndDependents();




    }


    out_nbParticles.dirtyConnections();

}

template <class T>
void DataBuffer1D<T>::appendCpu(T elem)
{
    dataCpu.refresh();
    resize(mNbElements+1);
    if(mbHasCpuStorage) {
        setCpuData(mNbElements-1, elem);
    }
//    if(mbHasBufferStorage) {
//        setBufferData(mNbElements-1, elem);
//    }

//    outmetadataBuffer.dirtyConnections();
//    dataCpu.dirtyItselfAndDependents();
//    dataBuffer.dirtyItselfAndDependents();

    out_nbParticles.dirtyConnections();
}

//template<class T>
//void DataBuffer1D<T>::bindBufferToTarget(GLenum target)
//{
//    if(mbHasBufferStorage) {
//        glBindBuffer(target, mGlidBuffer);
//    } else {
//        Application::sInsertDebugMessage(
//                "goglu::DataBuffer1D::bindToTarget : DataBuffer has no "
//                "BUFFER storage, cannot bind it to a target.",
//                GL_DEBUG_SEVERITY_HIGH);
//    }
//}




//==============================================================================
// PLUGS ACTIONS
//==============================================================================

template <class T>
Metadata1DBuffer DataBuffer1D<T>::out_metadataBuffer_getData() {
    dataBuffer.refresh();
    return metadataBuffer;
}
template <class T>
void DataBuffer1D<T>::out_metadataBuffer_receivePushConnectionUpdate(ConnectionAbstract* /*connection*/) {}


template <class T>
void DataBuffer1D<T>::in_metadataCpu_dirtyDependents()
{
    dataCpu.dirtyItselfAndDependents();
}
template <class T>
void DataBuffer1D<T>::in_metadataCpu_receive(Metadata1DCpu data)
{
    sourceStorageType = StorageType::CPU;
    metadataCpu = data;
//    dataBuffer.dirtyItselfAndDependents();
}
template <class T>
void DataBuffer1D<T>::in_metadataCpu_receivePullConnectionUpdate(ConnectionAbstract* /*connection*/) {
    if(sourceStorageType==StorageType::UNKNOWN) {
        sourceStorageType = StorageType::CPU;
    }
}


template <class T>
void DataBuffer1D<T>::in_metadataBuffer_dirtyDependents()
{
    dataBuffer.dirtyItselfAndDependents();
}
template <class T>
void DataBuffer1D<T>::in_metadataBuffer_receive(Metadata1DBuffer data)
{
    sourceStorageType = StorageType::BUFFER;
    dataCpu.dirtyItselfAndDependents();
    metadataBuffer = data;
}
template <class T>
void DataBuffer1D<T>::in_metadataBuffer_receivePullConnectionUpdate(ConnectionAbstract* /*connection*/) {}


template <class T>
unsigned int DataBuffer1D<T>::out_nbParticles_getData() {
    return mNbElements;
}


//==============================================================================
// Dirtyable actions
//==============================================================================

template <class T>
void DataBuffer1D<T>::dataCpu_dirtyDependents() {
    dataBuffer.dirtyItselfAndDependents();
//    outmetadataBuffer.dirtyConnections();
}

template <class T>
void DataBuffer1D<T>::dataCpu_updateData() {
    // TODO: is the is a metadata input for the cpu, check if it has changed, and modify the local data and metadata accordingly.
//    cout << "Dirt_dataCpu_updateData called." << endl;

    switch(sourceStorageType) {
    case StorageType::CPU:
        if(in_metadataCpu.hasPullConnection()) {
            delete dataCpu.get(); // TODO: safe to delete this pointer?
            dataCpu.set((T*)(in_metadataCpu.pullUpdate().dataPointer));
        } else {
            // nothing?
        }
        break;
    case StorageType::BUFFER:
        // copy data from buffer to cpu
//        {
//        T* bufferDataCopy = new T[metadataBuffer.nbElements];
//        for(int i=0; i<metadataBuffer.nbElements; i++) {
//            setCpuData(i, bufferDataCopy[i]);
//            cout << "buffer to CPU : " << metadataBuffer.nbElements << endl;
//        }
//        }

        glGetNamedBufferSubData(metadataBuffer.bufferId, 0, dataCpuSizeInBytes(), dataCpu.get());
//        sourceStorageType = StorageType::CPU;
        dataCpu.dirtyItselfAndDependents();
//        dataBuffer.dirtyItselfAndDependents();
        break;
    default:
        cout << "Dirt_dataCpu_updateData : unknown source storage type." << endl;
    }

}

template <class T>
void DataBuffer1D<T>::dataBuffer_dirtyDependents() {
    dataCpu.dirtyItselfAndDependents();
    out_metadataBuffer.dirtyConnections();
}

template <class T>
void DataBuffer1D<T>::dataBuffer_updateData() {
    switch(sourceStorageType) {
    case StorageType::CPU:
//        cout << "CPU" << endl;
        dataCpu.refresh();
        glNamedBufferSubData(metadataBuffer.bufferId, 0, dataCpuSizeInBytes(), dataCpu.get());
        break;
    case StorageType::BUFFER:
//        cout << "BUFFER" << endl;
        // TODO: check if pull connection for the bffer input plug.
        break;
    default:
        cout << "Dirt_dataBuffer_updateData : unknown source storage type." << endl;
    }
}

