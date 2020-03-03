
//#include "DataBuffer2D.h" // not necessary, but helps IDEs.

using namespace glm;

template<class T>
DataBuffer2D<T>::DataBuffer2D(unsigned int sizeX, unsigned int sizeY)
{
    _nbElementsX = sizeX;
    _nbElementsY = sizeY;
    _hasCpuStorage = false;
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
    delete[] _dataCpu;
    _hasCpuStorage = false;
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
    _hasTexture2DStorage = false;
}


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
    }
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

    float normalizedX = (x - _boundXMin)/(_boundXMax - _boundXMin)*mNbCellsX;
    float normalizedY = (y - _boundYMin)/(_boundYMax - _boundYMin)*mNbCellsX;



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

    if(_hasCpuStorage) {
        T* newData = new T[newSizeX*newSizeY];
        for(int i=0; i<nbElementsToCopyX; i++) {
            for(int j=0; j<nbElementsToCopyY; j++) {
                newData[newSizeY*j + i] = getCpuData(i, j);
            }
        }
        delete[] _dataCpu;
        _dataCpu = newData;
    }

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
    }

    _nbElementsX = newSizeX;
    _nbElementsY = newSizeY;
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

