#include "DataBuffer1D.h"

#include "Utils.h"

template<class T>
DataBuffer1D<T>::DataBuffer1D(unsigned int size) {
    _nbElements = size;
    _capacity = size;
    _hasCpuStorage = false;
    _hasBufferStorage = false;
    _hasTexture1DStorage = false;
    _hasTextureBufferStorage = false;
}

template <class T>
T DataBuffer1D<T>::getCpuData(unsigned int i)
{
    return _dataCpu[i];
}

template <class T>
T* DataBuffer1D<T>::getCpuDataPointer()
{
    return _dataCpu;
}

template <class T>
void DataBuffer1D<T>::setCpuData(unsigned int i, T data)
{
    _dataCpu[i] = data;
}

template<class T>
void DataBuffer1D<T>::createCpuStorage()
{
    _dataCpu = new T[_nbElements];
    _hasCpuStorage = true;

    // metadata
    _metadataCpu.dataPointer = _dataCpu;
}

template<class T>
void DataBuffer1D<T>::deleteCpuStorage()
{
    delete[] _dataCpu;
    _hasCpuStorage = false;
}

template<class T>
void DataBuffer1D<T>::createBufferStorage(GLenum dataType, unsigned int nbElementsPerComponent, GLbitfield flags)
{
    glCreateBuffers(1, &_glidBuffer);
    glNamedBufferStorage(_glidBuffer, _nbElements * nbElementsPerComponent * SizeOfEnumType(dataType), NULL,
        flags);
    _hasBufferStorage = true;

    // metadata
    _metadataBuffer.bufferId = _glidBuffer;
    _metadataBuffer.dataType = dataType;
    _metadataBuffer.nbElements = _nbElements;
    _metadataBuffer.nbElementsPerComponent = nbElementsPerComponent;
}

template<class T>
void DataBuffer1D<T>::deleteBufferStorage()
{
    glDeleteBuffers(1, &_glidBuffer);
    _hasBufferStorage = false;
}

//nbMipmapLevels=0 will put the maximum number of mipmaps. nbMipmapLevels=1
//will only put the top layer without any mipmaps below it.
template<class T>
void DataBuffer1D<T>::createTexture1DStorage(GLenum internalFormat,
    GLenum sizedInternalFormat,
    GLenum externalFormat,
    GLenum sizedExternalFormat)
{
    _texture1DInternalFormat = internalFormat;
    _texture1DSizedInternalFormat = sizedInternalFormat;
    _texture1DExternalFormat = externalFormat;
    _texture1DSizedExternalFormat = sizedExternalFormat;
    glCreateTextures(GL_TEXTURE_1D, 1, &_glidTexture1D);
    unsigned int nbLevels;
    // this computes nbLevels = log_2(mSize), which is the maximum number of
    // mipmap level we can have for this size.
    nbLevels = 1;
    unsigned int tempSize = _nbElements;
    while (tempSize >>= 1) ++nbLevels;

    glTextureStorage1D(_glidTexture1D, nbLevels,
        _texture1DSizedInternalFormat, _nbElements);
    glClearTexImage(_glidTexture1D, 0, _texture1DExternalFormat,
        _texture1DSizedExternalFormat, NULL);
    glGenerateTextureMipmap(_glidTexture1D);
    _hasTexture1DStorage = true;
}


template<class T>
void DataBuffer1D<T>::deleteTexture1DStorage()
{
    glDeleteTextures(1, &_glidTexture1D);
    _hasTexture1DStorage = false;
}


template<class T>
void DataBuffer1D<T>::createTextureBufferStorage(GLenum internalSizedFormat)
{
    if (!_hasBufferStorage) {
        createBufferStorage();
    }
    _textureBufferSizedInternalFormat = internalSizedFormat;
    glCreateTextures(GL_TEXTURE_BUFFER, 1, &_glidTextureBuffer);
    glTextureBuffer(_glidTextureBuffer, _textureBufferSizedInternalFormat,
        _glidBuffer);
    _hasTextureBufferStorage = true;
}


template<class T>
void DataBuffer1D<T>::deleteTextureBufferStorage()
{
    glDeleteTextures(1, &_glidTextureBuffer);
    _hasTextureBufferStorage = false;
}


template<class T>
void DataBuffer1D<T>::resize(unsigned int size)
{
    int nbElementsToCopy = glm::min(size, _nbElements);
    _nbElements = size;

    if (size <= _capacity) {
        if (_hasBufferStorage) {
            _metadataBuffer.nbElements = _nbElements;
        }
    }
    else {

        while (_capacity < size) {
            if (_capacity == 0) {
                _capacity = 1;
            }
            else {
                _capacity *= 2;
            }
        }

        if (_hasCpuStorage) {
            T* newData = new T[_capacity];
            for (int i = 0; i < nbElementsToCopy; i++) {
                newData[i] = getCpuData(i);
            }
            delete[] _dataCpu;
            _dataCpu = newData;
            _metadataCpu.dataPointer = _dataCpu;
        }

        if (_hasBufferStorage) {
            GLuint newBuffer;
            glCreateBuffers(1, &newBuffer);
            glNamedBufferStorage(newBuffer, _capacity * _metadataBuffer.nbElementsPerComponent * SizeOfEnumType(_metadataBuffer.dataType), NULL,
                GL_DYNAMIC_STORAGE_BIT | GL_MAP_READ_BIT |
                GL_MAP_WRITE_BIT);
            glCopyNamedBufferSubData(_glidBuffer, newBuffer, 0, 0,
                sizeof(T)*nbElementsToCopy);
            glDeleteBuffers(1, &_glidBuffer);
            _glidBuffer = newBuffer;

            _metadataBuffer.bufferId = newBuffer;
            _metadataBuffer.nbElements = _nbElements;
        }

        if (_hasTexture1DStorage) {
            GLuint newTexture;
            glCreateTextures(GL_TEXTURE_1D, 1, &newTexture);
            unsigned int nbLevels;
            // this computes nbLevels = log_2(mSize), which is the maximum number of
            // mipmap level we can have for this size.
            nbLevels = 1;
            unsigned int tempSize = _nbElements;
            while (tempSize >>= 1) ++nbLevels;
            glTextureStorage1D(newTexture, nbLevels, _texture1DSizedInternalFormat,
                _nbElements);
            glClearTexImage(newTexture, 0, _texture1DExternalFormat,
                _texture1DSizedExternalFormat, NULL);
            glCopyImageSubData(_glidTexture1D, GL_TEXTURE_1D, 0, 0, 0, 0,
                newTexture, GL_TEXTURE_1D, 0, 0, 0, 0,
                nbElementsToCopy, 1, 1);
            glDeleteTextures(1, &_glidTexture1D);
            _glidTexture1D = newTexture;
            glGenerateTextureMipmap(_glidTexture1D);
        }

        if (_hasTextureBufferStorage) {
            glDeleteTextures(1, &_glidTextureBuffer);
            glCreateTextures(GL_TEXTURE_BUFFER, 1, &_glidTextureBuffer);
            glTextureBuffer(_glidTextureBuffer, _textureBufferSizedInternalFormat,
                _glidBuffer);
        }
    }
}


template <class T>
void DataBuffer1D<T>::appendCpu(T elem)
{
    resize(_nbElements + 1);
    if (_hasCpuStorage) {
        setCpuData(_nbElements - 1, elem);
    }
}


template<class T>
void DataBuffer1D<T>::TransferDataCpuToBuffer()
{
    if(_hasBufferStorage||_hasTextureBufferStorage) {
        glNamedBufferSubData(_glidBuffer, 0, dataCpuSizeInBytes(), _dataCpu);
    }
}