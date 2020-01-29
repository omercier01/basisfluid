#include "VectorField2D.h"

using namespace goglu;
using namespace glm;

VectorField2D::VectorField2D(float parBoundXMin, float parBoundXMax, float parBoundYMin,
                             float parBoundYMax, unsigned int nbCellsX,
                             unsigned int nbCellsY,
                             VectorField2D::GridNodeLocation gridNodeLocation,
                             VectorField2D::BoundaryCondition boundaryCondition) :
    mNbCellsX(nbCellsX),
    mNbCellsY(nbCellsY),
    mVectors(1, 1),
    mGridNodeLocation(gridNodeLocation),
    mBoundaryCondition(boundaryCondition),

    in_vectorsMetadataImage2D(*this),
    out_gridNodeLocationsCpu(*this),
    out_gridNodeLocationsBuffer(*this),
    out_vectorValuesCpu(*this),
    out_divergenceCpu(*this),
    out_vectorsMetadataTexture2D(*this),
    out_vectorsMetadataImage2D(*this),

    in_int_vectorsMetadataImage2D(*this),
    out_int_vectorsMetadataImage2D(*this),

//    <<goglu>> {
//        Beacon VectorField2D_dirtyblesInit;
//    } <<goglu>>;
    
    GOGLU_BEGIN(
        component_constructor VectorField2D
    )GOGLU_END
    #include "../../gogluGeneratedCode/src/dataStructures/VectorField2D.cpp_d/snippet0.goglu"
{
    boundXMin.set(parBoundXMin);
    boundXMax.set(parBoundXMax);
    boundYMin.set(parBoundYMin);
    boundYMax.set(parBoundYMax);

    // resizes, but does not necessarily allocate the data.
    switch(gridNodeLocation) {
    case GridNodeLocation::CENTER:
        mVectors.resize(nbCellsX, nbCellsY);
        break;
    case GridNodeLocation::CORNER:
        mVectors.resize(nbCellsX+1, nbCellsY+1);
        break;
    }

    // internal connections
    createPullConnection(&out_int_vectorsMetadataImage2D,
                         &mVectors.in_metadataImage2D)->activate();
    createPushConnection(&out_int_vectorsMetadataImage2D,
                         &mVectors.in_metadataImage2D)->activate();
    createPullConnection(&mVectors.out_metadataImage2D,
                         &in_int_vectorsMetadataImage2D)->activate();

}

void VectorField2D::populateWithFunction(std::function<vec2 (float x, float y)> function)
{
    if(!mVectors.mbHasCpuStorage && !mVectors.mbHasTexture2DStorage) {
        Application::sInsertDebugMessage(
                "goglu::DataBuffer2D::populateWithFunction : No data storage"
                "created, populating has no effect.",
                GL_DEBUG_SEVERITY_LOW);
        return;
    }

    if(mVectors.mbHasCpuStorage) {
    
        vec2* mVectorsPointer = mVectors.getCpuDataPointer();
        unsigned int nxVec = mVectors.mNbElementsX;
    
        switch(mGridNodeLocation) {
        case GridNodeLocation::CENTER:
            for(unsigned int i=0; i<mNbCellsX; i++)
            for(unsigned int j=0; j<mNbCellsY; j++) {
                float x = boundXMin.get() + (i+0.5f)/mNbCellsX*(boundXMax.get()-boundXMin.get());
                float y = boundYMin.get() + (j+0.5f)/mNbCellsY*(boundYMax.get()-boundYMin.get());
    //            mVectors(i,j) = function(x, y);
                // TODO: setting dtaa one by one like this is probably sloooow!
                //mVectors.setCpuData(i, j, function(x, y));
                mVectorsPointer[nxVec*j + i] = function(x, y);
            }
            break;
        case GridNodeLocation::CORNER:
            for(unsigned int i=0; i<mNbCellsX+1; i++)
            for(unsigned int j=0; j<mNbCellsY+1; j++) {
                float x = boundXMin.get() + float(i)/mNbCellsX*(boundXMax.get()-boundXMin.get());
                float y = boundYMin.get() + float(j)/mNbCellsY*(boundYMax.get()-boundYMin.get());
    //            mVectors(i,j) = function(x, y);
                //mVectors.setCpuData(i, j, function(x, y));
                mVectorsPointer[nxVec*j + i] = function(x, y);
            }
            break;
        }
    }
    if(mVectors.mbHasTexture2DStorage) {
        // TODO: set data on texture too.
    }
    
    mVectors.sourceStorageType = DataBuffer2D<vec2>::StorageType::CPU;
    mVectors.dirtyData();
    
}

void VectorField2D::addFunction(std::function<vec2 (float, float)> function)
{
    if(!mVectors.mbHasCpuStorage && !mVectors.mbHasTexture2DStorage) {
        Application::sInsertDebugMessage(
                "goglu::DataBuffer2D::addFunction : No data storage"
                "created, populating has no effect.",
                GL_DEBUG_SEVERITY_LOW);
        return;
    }

    if(mVectors.mbHasCpuStorage) {
    
        vec2* mVectorsPointer = mVectors.getCpuDataPointer();
        unsigned int nxVec = mVectors.mNbElementsX;
    
        switch(mGridNodeLocation) {
        case GridNodeLocation::CENTER:
            for(unsigned int i=0; i<mNbCellsX; i++)
            for(unsigned int j=0; j<mNbCellsY; j++) {
                float x = boundXMin.get() + (i+0.5f)/mNbCellsX*(boundXMax.get()-boundXMin.get());
                float y = boundYMin.get() + (j+0.5f)/mNbCellsY*(boundYMax.get()-boundYMin.get());
                // TODO: setting data one by one like this is probably sloooow!
                //mVectors.addCpuData(i, j, function(x, y));
                mVectorsPointer[nxVec*j + i] += function(x, y);
            }
            break;
        case GridNodeLocation::CORNER:
            for(unsigned int i=0; i<mNbCellsX+1; i++)
            for(unsigned int j=0; j<mNbCellsY+1; j++) {
                float x = boundXMin.get() + float(i)/mNbCellsX*(boundXMax.get()-boundXMin.get());
                float y = boundYMin.get() + float(j)/mNbCellsY*(boundYMax.get()-boundYMin.get());
                //mVectors.addCpuData(i, j, function(x, y));
                mVectorsPointer[nxVec*j + i] += function(x, y);
            }
            break;
        }
    }
    if(mVectors.mbHasTexture2DStorage) {
        // TODO: set data on texture too.
    }
    
    mVectors.sourceStorageType = DataBuffer2D<vec2>::StorageType::CPU;
    mVectors.dirtyData();
    
}

vec2 VectorField2D::interp(vec2 pos)
{
    float x = pos.x;
    float y = pos.y;

    int indexXLeft, indexXRight;
    int indexYLeft, indexYRight;
    float weightX;
    float weightY;

    float normalizedX = (x - boundXMin.get())/(boundXMax.get() - boundXMin.get())*mNbCellsX;
    float normalizedY = (y - boundYMin.get())/(boundYMax.get() - boundYMin.get())*mNbCellsY;


    switch(mBoundaryCondition) {
    case BoundaryCondition::LINEAR:
        switch(mGridNodeLocation) {
        case GridNodeLocation::CENTER :
            if(normalizedX < 0.5) {
                indexXLeft = 0;
                indexXRight = 1;
                weightX = normalizedX - 0.5f;
            } else if(normalizedX >= mNbCellsX - 0.5f) {
                indexXLeft = mNbCellsX - 2;
                indexXRight = mNbCellsX - 1;
                weightX = normalizedX - (indexXLeft + 0.5f);
            } else {
                indexXLeft = clamp<int>(int(floor(normalizedX - 0.5f)), 0, mNbCellsX-2);
                indexXRight = indexXLeft + 1;
                weightX = normalizedX - (indexXLeft + 0.5f);
            }
            if(normalizedY < 0.5) {
                indexYLeft = 0;
                indexYRight = 1;
                weightY = normalizedY - 0.5f;
            } else if(normalizedY >= mNbCellsY - 0.5f) {
                indexYLeft = mNbCellsY - 2;
                indexYRight = mNbCellsY - 1;
                weightY = normalizedY - (indexYLeft + 0.5f);
            } else {
                indexYLeft = clamp<int>(int(floor(normalizedY - 0.5f)), 0, mNbCellsY-2);
                indexYRight = indexYLeft + 1;
                weightY = normalizedY - (indexYLeft + 0.5f);
            }
            break;
        case GridNodeLocation::CORNER :
            if(normalizedX < 0) {
                indexXLeft = 0;
                indexXRight = 1;
                weightX = normalizedX;
            } else if(normalizedX >= mNbCellsX) {
                indexXLeft = mNbCellsX - 1;
                indexXRight = mNbCellsX;
                weightX = normalizedX - indexXLeft;
            } else {
                indexXLeft = clamp<int>(int(floor(normalizedX)), 0, mNbCellsX-1);
                indexXRight = indexXLeft + 1;
                weightX = normalizedX - indexXLeft;
            }
            if(normalizedY < 0) {
                indexYLeft = 0;
                indexYRight = 1;
                weightY = normalizedY;
            } else if(normalizedY >= mNbCellsY) {
                indexYLeft = mNbCellsY - 1;
                indexYRight = mNbCellsY;
                weightY = normalizedY - indexYLeft;
            } else {
                indexYLeft = clamp<int>(int(floor(normalizedY)), 0, mNbCellsY-1);
                indexYRight = indexYLeft + 1;
                weightY = normalizedY - indexYLeft;
            }
            break;
        }
        break;
    case BoundaryCondition::FLAT:
        switch(mGridNodeLocation) {
        case GridNodeLocation::CENTER :
            if(normalizedX < 0.5) {
                indexXLeft = 0;
                indexXRight = 1;
                weightX = 0;
            } else if(normalizedX >= mNbCellsX - 0.5f) {
                indexXLeft = mNbCellsX - 2;
                indexXRight = mNbCellsX - 1;
                weightX = 1;
            } else {
                indexXLeft = clamp<int>(int(floor(normalizedX - 0.5f)), 0, mNbCellsX-2);
                indexXRight = indexXLeft + 1;
                weightX = normalizedX - (indexXLeft + 0.5f);
            }
            if(normalizedY < 0.5) {
                indexYLeft = 0;
                indexYRight = 1;
                weightY = 0;
            } else if(normalizedY >= mNbCellsY - 0.5f) {
                indexYLeft = mNbCellsY - 2;
                indexYRight = mNbCellsY - 1;
                weightY = 1;
            } else {
                indexYLeft = clamp<int>(int(floor(normalizedY - 0.5f)), 0, mNbCellsY-2);
                indexYRight = indexYLeft + 1;
                weightY = normalizedY - (indexYLeft + 0.5f);
            }
            break;
        case GridNodeLocation::CORNER :
            if(normalizedX < 0) {
                indexXLeft = 0;
                indexXRight = 1;
                weightX = 0;
            } else if(normalizedX >= mNbCellsX) {
                indexXLeft = mNbCellsX - 1;
                indexXRight = mNbCellsX;
                weightX = 1;
            } else {
                indexXLeft = clamp<int>(int(floor(normalizedX)), 0, mNbCellsX-1);
                indexXRight = indexXLeft + 1;
                weightX = normalizedX - indexXLeft;
            }
            if(normalizedY < 0) {
                indexYLeft = 0;
                indexYRight = 1;
                weightY = 0;
            } else if(normalizedY >= mNbCellsY) {
                indexYLeft = mNbCellsY - 1;
                indexYRight = mNbCellsY;
                weightY = 1;
            } else {
                indexYLeft = clamp<int>(int(floor(normalizedY)), 0, mNbCellsY-1);
                indexYRight = indexYLeft + 1;
                weightY = normalizedY - indexYLeft;
            }
            break;
        }
        break;
    case BoundaryCondition::ZERO:
        switch(mGridNodeLocation) {
        case GridNodeLocation::CENTER :
            if(normalizedX < 0.5) {
                return vec2(0, 0);
            } else if(normalizedX >= mNbCellsX - 0.5f) {
                return vec2(0, 0);
            } else {
                indexXLeft = clamp<int>(int(floor(normalizedX - 0.5f)), 0, mNbCellsX-2);
                indexXRight = indexXLeft + 1;
                weightX = normalizedX - (indexXLeft + 0.5f);
            }
            if(normalizedY < 0.5) {
                return vec2(0, 0);
            } else if(normalizedY >= mNbCellsY - 0.5f) {
                return vec2(0, 0);
            } else {
                indexYLeft = clamp<int>(int(floor(normalizedY - 0.5f)), 0, mNbCellsY-2);
                indexYRight = indexYLeft + 1;
                weightY = normalizedY - (indexYLeft + 0.5f);
            }
            break;
        case GridNodeLocation::CORNER :
            if(normalizedX < 0) {
                return vec2(0, 0);
            } else if(normalizedX >= mNbCellsX) {
                return vec2(0, 0);
            } else {
                indexXLeft = clamp<int>(int(floor(normalizedX)), 0, mNbCellsX-1);
                indexXRight = indexXLeft + 1;
                weightX = normalizedX - indexXLeft;
            }
            if(normalizedY < 0) {
                return vec2(0, 0);
            } else if(normalizedY >= mNbCellsY) {
                return vec2(0, 0);
            } else {
                indexYLeft = clamp<int>(int(floor(normalizedY)), 0, mNbCellsY-1);
                indexYRight = indexYLeft + 1;
                weightY = normalizedY - indexYLeft;
            }
            break;
        }
        break;
    case BoundaryCondition::MIRROR:
        switch(mGridNodeLocation) {
        case GridNodeLocation::CENTER :
            if(normalizedX < 0) {
                normalizedX = -normalizedX;
            }
            if(normalizedX >= mNbCellsX - 0.5f) {
                float xInDomain01 = normalizedX/mNbCellsX;
                unsigned int temp = int(floor(xInDomain01)) % 2;
                if(temp == 0) {
                    normalizedX = mNbCellsX * fract(xInDomain01);
                } else {
                    normalizedX = mNbCellsX * (1 - fract(xInDomain01));
                }
            }
            if(normalizedX < 0.5) {
                indexXLeft = 0;
                indexXRight = indexXLeft;
                weightX = 0;
            } else if(normalizedX >= mNbCellsX - 0.5f) {
                indexXLeft = mNbCellsX - 1;
                indexXRight = indexXLeft;
                weightX = 0;
            } else {
                indexXLeft = clamp<int>(int(floor(normalizedX - 0.5f)), 0, mNbCellsX-2);
                indexXRight = indexXLeft + 1;
                weightX = normalizedX - (indexXLeft + 0.5f);
            }

            if(normalizedY < 0) {
                normalizedY = -normalizedY;
            }
            if(normalizedY >= mNbCellsY - 0.5f) {
                float yInDomain01 = normalizedY/mNbCellsY;
                unsigned int temp = int(floor(yInDomain01)) % 2;
                if(temp == 0) {
                    normalizedY = mNbCellsY * fract(yInDomain01);
                } else {
                    normalizedY = mNbCellsY * (1 - fract(yInDomain01));
                }
            }
            if(normalizedY < 0.5) {
                indexYLeft = 0;
                indexYRight = indexYLeft;
                weightY = normalizedY - indexYLeft;
            } else if(normalizedY >= mNbCellsY - 0.5f) {
                indexYLeft = mNbCellsY - 1;
                indexYRight = indexYLeft;
                weightY = normalizedY - indexYLeft;
            } else {
                indexYLeft = clamp<int>(int(floor(normalizedY - 0.5f)), 0, mNbCellsY-2);
                indexYRight = indexYLeft + 1;
                weightY = normalizedY - (indexYLeft + 0.5f);
            }
            break;
        case GridNodeLocation::CORNER :
            if(normalizedX < 0) {
                normalizedX = -normalizedX;
            }
            if(normalizedX >= mNbCellsX) {
                float xInDomain01 = normalizedX/mNbCellsX;
                unsigned int temp = int(floor(xInDomain01)) % 2;
                if(temp == 0) {
                    normalizedX = mNbCellsX * fract(xInDomain01);
                } else {
                    normalizedX = mNbCellsX * (1 - fract(xInDomain01));
                }
            }
            indexXLeft = clamp<int>(int(floor(normalizedX)), 0, mNbCellsX-1);
            indexXRight = indexXLeft + 1;
            weightX = normalizedX - indexXLeft;

            if(normalizedY < 0) {
                normalizedY = -normalizedY;
            }
            if(normalizedY >= mNbCellsY) {
                float yInDomain01 = normalizedY/mNbCellsY;
                unsigned int temp = int(floor(yInDomain01)) % 2;
                if(temp == 0) {
                    normalizedY = mNbCellsY * fract(yInDomain01);
                } else {
                    normalizedY = mNbCellsY * (1 - fract(yInDomain01));
                }
            }
            indexYLeft = clamp<int>(int(floor(normalizedY)), 0, mNbCellsY-1);
            indexYRight = indexYLeft + 1;
            weightY = normalizedY - indexYLeft;
            break;
        }
        break;
    case BoundaryCondition::PERIODIC:
        switch(mGridNodeLocation) {
        case GridNodeLocation::CENTER :
            normalizedX -= mNbCellsX * floor(normalizedX/mNbCellsX);
            indexXLeft = clamp<int>(int(floor(normalizedX - 0.5f)), 0, mNbCellsX-2);
            indexXRight = indexXLeft + 1;
            weightX = normalizedX - (indexXLeft+0.5f);
            normalizedY -= mNbCellsY * floor(normalizedY/mNbCellsY);
            indexYLeft = clamp<int>(int(floor(normalizedY - 0.5f)), 0, mNbCellsY-2);
            indexYRight = indexYLeft + 1;
            weightY = normalizedY - (indexYLeft+0.5f);
            break;
        case GridNodeLocation::CORNER :
            normalizedX -= mNbCellsX * floor(normalizedX/mNbCellsX);
            indexXLeft = clamp<int>(int(floor(normalizedX)), 0, mNbCellsX-1);
            indexXRight = indexXLeft + 1;
            weightX = normalizedX - indexXLeft;
            normalizedY -= mNbCellsY * floor(normalizedY/mNbCellsY);
            indexYLeft = clamp<int>(int(floor(normalizedY)), 0, mNbCellsY-1);
            indexYRight = indexYLeft + 1;
            weightY = normalizedY - indexYLeft;
            break;
        }
        break;
    }

    //TODO: indices might still be out of range because of numerical errors,
    //check for that and if it happens, the out-of-bound indices will have to
    //be dealt with in terms of indices, not positions.

//    return (1-weightX)*(1-weightY)*mVectors(indexXLeft , indexYLeft ) +
//           (1-weightX)*(  weightY)*mVectors(indexXLeft , indexYRight) +
//           (  weightX)*(1-weightY)*mVectors(indexXRight, indexYLeft ) +
//           (  weightX)*(  weightY)*mVectors(indexXRight, indexYRight);

    mVectors.refreshCpuData();
    vec2* mVectorsPointer = mVectors.getCpuDataPointer();
    
    const unsigned int nx = mVectors.mNbElementsX;
    
    return (1-weightX)*(
               (1-weightY)*mVectorsPointer[nx*indexYLeft + indexXLeft] +
               (  weightY)*mVectorsPointer[nx*indexYRight + indexXLeft]
           )
           +
           (weightX)*(
               (1-weightY)*mVectorsPointer[nx*indexYLeft + indexXRight] +
               (  weightY)*mVectorsPointer[nx*indexYRight + indexXRight]
           );
           
}

void VectorField2D::addVectorCpuData(unsigned int i, unsigned int j, vec2 data)
{
    mVectors.addCpuData(i, j, data);
    out_vectorsMetadataImage2D.dirtyConnections();
    out_vectorsMetadataTexture2D.dirtyConnections();
    out_vectorValuesCpu.dirtyConnections();
}

void VectorField2D::setVectorCpuData(unsigned int i, unsigned int j, vec2 data)
{
    mVectors.setCpuData(i, j, data);
    out_vectorsMetadataImage2D.dirtyConnections();
    out_vectorsMetadataTexture2D.dirtyConnections();
    out_vectorValuesCpu.dirtyConnections();
}

glm::vec2 VectorField2D::getVectorCpuData(unsigned int i, unsigned int j)
{
    return mVectors.getCpuData(i, j);
}

void VectorField2D::createVectorCpuStorage()
{
    mVectors.createCpuStorage();
}


void VectorField2D::createVectorTexture2DStorage(
    GLenum internalFormat,
    GLenum sizedInternalFormat,
    GLenum externalFormat,
    GLenum sizedExternalFormat,
    unsigned int nbMipmapLevels)
{
    mVectors.createTexture2DStorage(internalFormat, sizedInternalFormat,
                                    externalFormat, sizedExternalFormat,
                                    nbMipmapLevels);
    vectorTexture2DStorageAdjustBoundaryCondition();
}


void VectorField2D::createVectorTexture2DStorage(
    GLenum internalFormat,
    GLenum sizedInternalFormat,
    GLenum externalFormat,
    GLenum sizedExternalFormat)
{
    mVectors.createTexture2DStorage(internalFormat, sizedInternalFormat,
                                    externalFormat, sizedExternalFormat);
    vectorTexture2DStorageAdjustBoundaryCondition();
}


void VectorField2D::vectorTexture2DStorageAdjustBoundaryCondition()
{
    //LINEAR, FLAT, ZERO, MIRROR, PERIODIC
    float borderColor[] = { 0.0f, 0.0f, 0.0f, 0.0f };
    switch(mBoundaryCondition) {
        case VectorField2D::BoundaryCondition::ZERO:
            glTextureParameteri(mVectors.mGlidTexture2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
            glTextureParameteri(mVectors.mGlidTexture2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
            glTextureParameterfv(mVectors.mGlidTexture2D, GL_TEXTURE_BORDER_COLOR, borderColor);
            break;
        case VectorField2D::BoundaryCondition::LINEAR:
            // not supported by OpgnGL, substiture with FLAT boundary condition
            glTextureParameteri(mVectors.mGlidTexture2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTextureParameteri(mVectors.mGlidTexture2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
            Application::sInsertDebugMessage(
                "goglu::VectorField2D::vectorTexture2DStorageAdjustBoundaryCondition : "
                "LINEAR boundary condition not supported on GPU, FLAT used instead.",
                GL_DEBUG_SEVERITY_LOW);
            break;
        case VectorField2D::BoundaryCondition::FLAT:
            glTextureParameteri(mVectors.mGlidTexture2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTextureParameteri(mVectors.mGlidTexture2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
            break;
        case VectorField2D::BoundaryCondition::MIRROR:
            glTextureParameteri(mVectors.mGlidTexture2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
            glTextureParameteri(mVectors.mGlidTexture2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
            break;
        case VectorField2D::BoundaryCondition::PERIODIC:
            glTextureParameteri(mVectors.mGlidTexture2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
            glTextureParameteri(mVectors.mGlidTexture2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
            break;
        default:
            // TODO: should not happen
        break;
    }
}


void VectorField2D::setVectorsFromImage(Metadata2DImage2D metadataImage2D, unsigned int level)
{
    mVectors.setFromImage(metadataImage2D, level);
}

unsigned int VectorField2D::nbElementsX()
{
    switch(mGridNodeLocation) {
    case GridNodeLocation::CENTER:
        return mNbCellsX;
        break;
    case GridNodeLocation::CORNER:
        return mNbCellsX + 1;
        break;
    default:
        // TODO: error, unknows grid location type.
        return 0;
    }
}

unsigned int VectorField2D::nbElementsY()
{
    switch(mGridNodeLocation) {
    case GridNodeLocation::CENTER:
        return mNbCellsY;
        break;
    case GridNodeLocation::CORNER:
        return mNbCellsY + 1;
        break;
    default:
        // TODO: error, unknows grid location type.
        return 0;
    }
}


void VectorField2D::setBounds(float inBoundXMin, float inBoundXMax, float inBoundYMin, float inBoundYMax)
{
    boundXMin.set(inBoundXMin);
    boundXMax.set(inBoundXMax);
    boundYMin.set(inBoundYMin);
    boundYMax.set(inBoundYMax);
}


uvec2 VectorField2D::pointToClosestIndex(vec2 point)
{
    uvec2 index;

    float normalizedX = (point.x - boundXMin.get())/(boundXMax.get() - boundXMin.get())*mNbCellsX;
    float normalizedY = (point.y - boundYMin.get())/(boundYMax.get() - boundYMin.get())*mNbCellsY;


    switch(mGridNodeLocation) {
    case GridNodeLocation::CENTER :
        if(normalizedX < 0) {
            index.x = 0;
        } else if(normalizedX >= mNbCellsX) {
            index.x = mNbCellsX-1;
        } else {
            index.x = clamp<int>(int(floor(normalizedX)), 0, mNbCellsX-1);
        }
        if(normalizedY < 0) {
            index.y = 0;
        } else if(normalizedY >= mNbCellsY) {
            index.y = mNbCellsY-1;
        } else {
            index.y = clamp<int>(int(floor(normalizedY)), 0, mNbCellsY-1);
        }
        break;
    case GridNodeLocation::CORNER :
        if(normalizedX < 0) {
            index.x = 0;
        } else if(normalizedX >= mNbCellsX) {
            index.x = mNbCellsX;
        } else {
            index.x = clamp<int>(int(floor(normalizedX + 0.5)), 0, mNbCellsX);
        }
        if(normalizedY < 0) {
            index.y = 0;
        } else if(normalizedY >= mNbCellsY) {
            index.y = mNbCellsY;
        } else {
            index.y = clamp<int>(int(floor(normalizedY + 0.5)), 0, mNbCellsY);
        }
        break;
    }

    return index;
}

uvec2 VectorField2D::pointToFlooredIndex(vec2 point)
{
    uvec2 index;

    float normalizedX = (point.x - boundXMin.get())/(boundXMax.get() - boundXMin.get())*mNbCellsX;
    float normalizedY = (point.y - boundYMin.get())/(boundYMax.get() - boundYMin.get())*mNbCellsY;


    switch(mGridNodeLocation) {
    case GridNodeLocation::CENTER :
        if(normalizedX < 0) {
            index.x = 0;
        } else if(normalizedX >= mNbCellsX) {
            index.x = mNbCellsX-1;
        } else {
            index.x = clamp<int>(int(floor(normalizedX-0.5)), 0, mNbCellsX-1);
        }
        if(normalizedY < 0) {
            index.y = 0;
        } else if(normalizedY >= mNbCellsY) {
            index.y = mNbCellsY-1;
        } else {
            index.y = clamp<int>(int(floor(normalizedY-0.5)), 0, mNbCellsY-1);
        }
        break;
    case GridNodeLocation::CORNER :
        if(normalizedX < 0) {
            index.x = 0;
        } else if(normalizedX >= mNbCellsX) {
            index.x = mNbCellsX;
        } else {
            index.x = clamp<int>(int(floor(normalizedX)), 0, mNbCellsX);
        }
        if(normalizedY < 0) {
            index.y = 0;
        } else if(normalizedY >= mNbCellsY) {
            index.y = mNbCellsY;
        } else {
            index.y = clamp<int>(int(floor(normalizedY)), 0, mNbCellsY);
        }
        break;
    }

    return index;
}



vec2 VectorField2D::indexToPosition(uvec2 index)
{
    vec2 pos;

    switch(mGridNodeLocation) {
    case GridNodeLocation::CENTER :
        pos.x = boundXMin.get() + (index.x + 0.5f)*(boundXMax.get()-boundXMin.get())/mNbCellsX;
        pos.y = boundYMin.get() + (index.y + 0.5f)*(boundYMax.get()-boundYMin.get())/mNbCellsY;
        break;
    case GridNodeLocation::CORNER :
        pos.x = boundXMin.get() + (index.x)*(boundXMax.get()-boundXMin.get())/mNbCellsX;
        pos.y = boundYMin.get() + (index.y)*(boundYMax.get()-boundYMin.get())/mNbCellsY;
        break;
    }

    return pos;
}




//void VectorField2D::createGridNodeLocationsCpuStorage() {
//    gridNodeLocations.get()->createCpuStorage();
//    gridNodeLocations.get()->resize(nbElementsX() * nbElementsY());
//    gridNodeLocations.get()->dataCpu.dirtyItselfAndDependents();
//}

//void VectorField2D::createGridNodeLocationsBufferStorage(GLenum dataType) {
//    gridNodeLocations.createBufferStorage(dataType, 2);
//    gridNodeLocations.resize(nbElementsX() * nbElementsY());
//    gridNodeLocations.get()->dataBuffer.dirtyItselfAndDependents();
//}




//==============================================================================
// PLUGS ACTIONS
//==============================================================================

// copies and returns a linear buffer containing the grid positions.
Metadata1DCpu VectorField2D::out_gridNodeLocationsCpu_getData() {

    mVectors.dataCpu.refresh();

    Metadata1DCpu metadata;
    vec2* data;
    float cellWidthX = (boundXMax.get() - boundXMin.get())/mNbCellsX;
    float cellWidthY = (boundYMax.get() - boundYMin.get())/mNbCellsY;

    switch(mGridNodeLocation) {
    case GridNodeLocation::CENTER:
        data = new vec2[mNbCellsX*mNbCellsY];
        for(unsigned int i=0; i<mNbCellsX; i++) {
        for(unsigned int j=0; j<mNbCellsY; j++) {
            data[mNbCellsY*j + i] = vec2(boundXMin.get() + (i+0.5f)*cellWidthX,
                                         boundYMin.get() + (j+0.5f)*cellWidthY);
        }}
        break;
    case GridNodeLocation::CORNER:
        data = new vec2[(mNbCellsX+1)*(mNbCellsY+1)];
        for(unsigned int i=0; i<=mNbCellsX; i++) {
        for(unsigned int j=0; j<=mNbCellsY; j++) {
            data[(mNbCellsY+1)*j + i] = vec2(boundXMin.get() + i*cellWidthX,
                                             boundYMin.get() + j*cellWidthY);
        }}
        break;
    }

    metadata.dataPointer = data;
    return metadata;
}
void VectorField2D::out_gridNodeLocationsCpu_receivePushConnectionUpdate(ConnectionAbstract* /*connection*/) {}


GLuint VectorField2D::out_gridNodeLocationsBuffer_getData() {
    //TODO: generate grid node postions.
    return 0;
}
void VectorField2D::out_gridNodeLocationsBuffer_receivePushConnectionUpdate(ConnectionAbstract* /*connection*/) {}


// creates and return a linear copy of the vector data.
Metadata1DCpu VectorField2D::out_vectorValuesCpu_getData() {

    Metadata1DCpu metadata;
    vec2* data;

    switch(mGridNodeLocation) {
    case GridNodeLocation::CENTER:
        data = new vec2[mNbCellsX*mNbCellsY];
        break;
    case GridNodeLocation::CORNER:
        data = new vec2[(mNbCellsX+1)*(mNbCellsY+1)];
        break;
    }

//    mVectors.sourceStorageType = DataBuffer2D<vec2>::StorageType::TEXTURE2D;
    mVectors.dataCpu.dirtyItselfAndDependents();
    memcpy(data, mVectors.dataCpu.get(), mVectors.dataSizeInBytes());

    metadata.dataPointer = data;
    return metadata;
}
void VectorField2D::out_vectorValuesCpu_receivePushConnectionUpdate(ConnectionAbstract* /*connection*/) {}


Metadata2DCpu VectorField2D::out_divergenceCpu_getData() {
    Metadata2DCpu metadata;
    float* data;

    float hx = (boundXMax.get() - boundXMin.get())/mNbCellsX;
    float hy = (boundYMax.get() - boundYMin.get())/mNbCellsY;
    vec2 vecHx = vec2(hx, 0);
    vec2 vecHy = vec2(0, hy);

    switch(mGridNodeLocation) {
    case GridNodeLocation::CENTER:
        data = new float[mNbCellsX*mNbCellsY];

        for(unsigned int i=0; i<mNbCellsX; ++i) {
        for(unsigned int j=0; j<mNbCellsY; ++j) {
            vec2 c = vec2(boundXMin.get() + (i+0.5)*hx, boundYMin.get() + (j+0.5)*hy);
            data[mNbCellsY*j + i] = 0.5f/hx*(interp(c+vecHx).x - interp(c-vecHx).x) +
                                    0.5f/hy*(interp(c+vecHy).y - interp(c-vecHy).y);
//            data[mNbCellsY*j + i] *= 100;
//            data[mNbCellsY*j + i] = 0.5;
        }}

        break;
    case GridNodeLocation::CORNER:
        data = new float[(mNbCellsX+1)*(mNbCellsY+1)];

        for(unsigned int i=0; i<mNbCellsX+1; ++i) {
        for(unsigned int j=0; j<mNbCellsY+1; ++j) {
            vec2 c = vec2(boundXMin.get() + i*hx, boundYMin.get() + j*hy);
            data[(mNbCellsY+1)*j + i] = 0.5f/hx*(interp(c+vecHx).x - interp(c-vecHx).x) +
                                        0.5f/hy*(interp(c+vecHy).y - interp(c-vecHy).y);
//            data[(mNbCellsY+1)*j + i] = 1;
        }}

        break;
    }

//    memcpy(data, mVectors.dataCpu.get(), mVectors.dataSizeInBytes());


    metadata.dataPointer = data;
    return metadata;
}
void VectorField2D::out_divergenceCpu_receivePushConnectionUpdate(ConnectionAbstract* /*connection*/) {}


Metadata2DTexture2D VectorField2D::out_vectorsMetadataTexture2D_getData()
{
    return mVectors.out_metadataTexture2D.getData();
}
void VectorField2D::out_vectorsMetadataTexture2D_receivePushConnectionUpdate(ConnectionAbstract* /*connection*/) {}


Metadata2DImage2D VectorField2D::out_vectorsMetadataImage2D_getData(unsigned int level)
{
    return mVectors.out_metadataImage2D.getData(level);
}
void VectorField2D::out_vectorsMetadataImage2D_receivePushConnectionUpdate(ConnectionAbstract* /*connection*/) {}


void VectorField2D::in_vectorsMetadataImage2D_dirtyDependents()
{
    out_int_vectorsMetadataImage2D.dirtyConnections();
}
void VectorField2D::in_vectorsMetadataImage2D_receive(Metadata2DImage2D data)
{
    out_int_vectorsMetadataImage2D.pushUpdate(data);
}
void VectorField2D::in_vectorsMetadataImage2D_receivePullConnectionUpdate(ConnectionAbstract* /*connection*/) {}


Metadata2DImage2D VectorField2D::out_int_vectorsMetadataImage2D_getData()
{
    return in_vectorsMetadataImage2D.pullUpdate();
}
void VectorField2D::out_int_vectorsMetadataImage2D_receivePushConnectionUpdate(ConnectionAbstract* /*connection*/) {}


void VectorField2D::in_int_vectorsMetadataImage2D_dirtyDependents()
{
    out_int_vectorsMetadataImage2D.dirtyConnections();
}
void VectorField2D::in_int_vectorsMetadataImage2D_receivePullConnectionUpdate(ConnectionAbstract* /*connection*/) {}


//<<goglu>> {
//    Beacon VectorField2D_dirtyblesDefs;
//} <<goglu>>;

GOGLU_BEGIN(
    component_definitions VectorField2D
)GOGLU_END
#include "../../gogluGeneratedCode/src/dataStructures/VectorField2D.cpp_d/snippet1.goglu"
