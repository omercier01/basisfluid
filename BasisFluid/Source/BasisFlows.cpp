//#include "glm/glm.hpp"

#include "Application.h"
#include "BasisFlows.h"

#include "Utils.h"

#include <set>
#include <iostream>
#include <vector>
#include <fstream>
#include <tuple>
#include <sstream>

using namespace glm;
using namespace std;

template <class T>
T min3(T a, T b, T c) { return glm::min(a, glm::min(b, c)); }

template <class T>
T max3(T a, T b, T c) { return glm::max(a, glm::max(b, c)); }

template <class T>
int minVec(T v) { return glm::min<int>(v.x, v.y); }


bool IntersectionInteriorEmpty(const BasisSupport& sup1, const BasisSupport& sup2)
{
    return (glm::max(sup1.left, sup2.left) >= glm::min(sup1.right, sup2.right))
        || (glm::max(sup1.bottom, sup2.bottom) >= glm::min(sup1.top, sup2.top))
        ;
}


bool BasisFlow::EmptyIntersectionWithBasis(const BasisFlow& b) const {
    return IntersectionInteriorEmpty(getSupport(), b.getSupport());
}


scalar_inversion_inner Application::InverseBBMatrixMain(
    unsigned int iRow, scalar_inversion_storage* vecX, scalar_inversion_storage* vecB,
    BasisFlow* basisDataPointer, unsigned int basisBitMask, float alpha)
{
    scalar_inversion_inner xSquaredDiff = 0;
    scalar_inversion_inner tempX = scalar_inversion_inner(vecB[iRow]);

    vector<CoeffBBDecompressedIntersectionInfo>& intersectionInfos = _coeffsBBDecompressedIntersections[iRow];
    for (const CoeffBBDecompressedIntersectionInfo& inter : intersectionInfos) {
        //            if( atLeastOneBitNotSet(basisDataPointer[inter.j].bitFlags, basisBitMask) ) {
        //                continue;
        //            }
        //            tempX -= inter.coeff * scalar_inversion_inner(vecX[inter.j]);
        if (AllBitsSet(basisDataPointer[inter.j].bitFlags, basisBitMask)) {
            tempX -= inter.coeff * scalar_inversion_inner(vecX[inter.j]);
        }

    }

    // TODO: make all bases normalized in 3D so we can save a fetch + division here
    scalar_inversion_storage newX = scalar_inversion_storage(tempX / scalar_inversion_inner(basisDataPointer[iRow].normSquared));
    scalar_inversion_storage oldX = vecX[iRow];

    xSquaredDiff = Sqr(scalar_inversion_inner(newX - oldX));
    vecX[iRow] = (1 - alpha)*oldX + alpha * newX;
    return scalar_inversion_inner(xSquaredDiff);
}



// Solves BB * x = b using Gauss-Seidel.
// Since bases are ordered by wavenumber, doing only one step of GS is equivalent to projection by frequency layer.
// The bitMask indicates the bits that must be set for the basis to be used in the inversion.
void Application::InverseBBMatrix(
    DataBuffer1D<scalar_inversion_storage>* vecX,
    DataBuffer1D<scalar_inversion_storage>* vecB,
    float tol, unsigned int basisBitMask, unsigned int minFreq)
{



    uint maxNbItMatBBInversion = _maxNbItMatBBInversion;
    uint n = vecX->_nbElements;
    float alpha = _relaxationAlpha;

    // get references
    //vecX->refreshCpuData();
    scalar_inversion_storage* vecXPointer = vecX->getCpuDataPointer();
    //_vecTemp->refreshCpuData();
    scalar_inversion_storage* vecTempPointer = _vecTemp->getCpuDataPointer();
    //vecB->refreshCpuData();
    scalar_inversion_storage* vecBPointer = vecB->getCpuDataPointer();
    //_basisFlowParams->refreshCpuData();
    BasisFlow* basisFlowParamsPointer = _basisFlowParams->getCpuDataPointer();
    //_intersectingBasesSignificantBBIds->refreshCpuData();

    if (_bInversionGaussSeidel)
    {

        // zero x
        //#pragma omp parallel for
        for (int i = 0; i < int(n); i++) {
            vecXPointer[i] = 0.;
        }

        scalar_inversion_inner xSquaredDiff = 999999;
        uint iIt = 0;
        if (maxNbItMatBBInversion == 0)
        {
            for (uint iRow = 0; iRow < n; iRow++) {
                if (
                    AllBitsSet(basisFlowParamsPointer[iRow].bitFlags, basisBitMask) &&
                    uint(basisFlowParamsPointer[iRow].freqLvl.x) >= minFreq
                    && uint(basisFlowParamsPointer[iRow].freqLvl.y) >= minFreq
                    ) {
                    vecXPointer[iRow] = vecBPointer[iRow];
                }
            }
        }
        else {
            while (xSquaredDiff > tol && iIt < maxNbItMatBBInversion) {
                xSquaredDiff = 0;

#if INVERSION_OPENMP

                for (int iBasisGroup = 0; iBasisGroup < _orthogonalBasisGroupIds.size(); iBasisGroup++)
                {
                    vector<unsigned int> ids = _orthogonalBasisGroupIds[iBasisGroup];
                    if (ids.size() >= _minSizeParallelInverse)
                    {
#pragma omp parallel for num_threads(7) reduction(+:xSquaredDiff) shared(vecXPointer,vecBPointer,ids,basisFlowParamsPointer) firstprivate(basisBitMask,alpha,minFreq) default(none)
                        for (int idid = 0; idid < ids.size(); idid++) {
                            int iRow = ids[idid];
                            if (
                                AllBitsSet(basisFlowParamsPointer[iRow].bitFlags, basisBitMask) &&
                                uint(basisFlowParamsPointer[iRow].freqLvl.x) >= minFreq
                                && uint(basisFlowParamsPointer[iRow].freqLvl.y) >= minFreq
                                ) {
                                xSquaredDiff += InverseBBMatrixMain(iRow, vecXPointer, vecBPointer, basisFlowParamsPointer, basisBitMask, alpha);
                            }
                        }
                    }
                    else {
                        for (int idid = 0; idid < ids.size(); idid++) {
                            int iRow = ids[idid];
                            if (
                                AllBitsSet(basisFlowParamsPointer[iRow].bitFlags, basisBitMask) &&
                                uint(basisFlowParamsPointer[iRow].freqLvl.x) >= minFreq
                                && uint(basisFlowParamsPointer[iRow].freqLvl.y) >= minFreq
                                ) {
                                xSquaredDiff += InverseBBMatrixMain(iRow, vecXPointer, vecBPointer, basisFlowParamsPointer, basisBitMask, alpha);
                            }
                        }
                    }
                }

#else


                //                    for(uint iRow=0; iRow < n; iRow++) {
                //                        cout << iRow << " - " << vecXPointer[iRow] << endl;
                //                    }


                for (uint iRow = 0; iRow < n; iRow++) {

                    //                        if(iRow == 6858) {
                    //                            cout << "QWEQWEQWE" << endl;
                    //                        }

                    if (atLeastOneBitNotSet(basisFlowParamsPointer[iRow].bitFlags, basisBitMask)) {
                        continue;
                    }

                    //if( allBitsSet(basisFlowParamsPointer[iRow].bitFlags, basisBitMask) ) {
                    xSquaredDiff += inverseBBMatrixMain(iRow, vecXPointer, vecBPointer, basisFlowParamsPointer, basisBitMask, alpha);
                    //}
                }

#endif

#if INVERSION_INNER_DOUBLE_PRECISION
                xSquaredDiff = std::sqrt(xSquaredDiff / n);
#else
                xSquaredDiff = std::sqrtf(xSquaredDiff / n);
#endif

                iIt++;

            }
        }
    }
    else {

        //Invert with Jacobi
        // TODO: make this paralll and as fast as Gauss-Seidel if we want to compare them fairly.

        // zero x
        for (uint i = 0; i < n; i++) {
            vecXPointer[i] = 0.f;
        }

        scalar_inversion_inner xSquaredDiff = 999999;
        uint iIt = 0;
        while (xSquaredDiff > tol && iIt < maxNbItMatBBInversion) {
            xSquaredDiff = 0;
            for (uint iRow = 0; iRow < n; iRow++) {

                if (AtLeastOneBitNotSet(basisFlowParamsPointer[iRow].bitFlags, basisBitMask)) {
                    continue;
                }
                scalar_inversion_inner tempX = scalar_inversion_inner(vecBPointer[iRow]);

                vector<CoeffBBDecompressedIntersectionInfo>& intersectionInfos = _coeffsBBDecompressedIntersections[iRow];
                for (CoeffBBDecompressedIntersectionInfo& inter : intersectionInfos) {
                    if (AllBitsSet(basisFlowParamsPointer[inter.j].bitFlags, basisBitMask)) {
                        tempX -= inter.coeff * scalar_inversion_inner(vecXPointer[inter.j]);
                    }
                }
                // TODO: remove division once 3D bases have norm 1
                BasisFlow bi = basisFlowParamsPointer[iRow];
                scalar_inversion_storage newX = scalar_inversion_storage(tempX / scalar_inversion_inner(bi.normSquared));
                xSquaredDiff += Sqr(scalar_inversion_inner(newX - vecXPointer[iRow]));
                vecTempPointer[iRow] = newX;

            }
            for (uint i = 0; i < n; i++) {
                vecXPointer[i] = (1 - alpha)*vecXPointer[i] + alpha * vecTempPointer[i];
            }

#if INVERSION_INNER_DOUBLE_PRECISION
            xSquaredDiff = std::sqrt(xSquaredDiff / n);
#else
            xSquaredDiff = std::sqrtf(xSquaredDiff / n);
#endif

            iIt++;
        }

    }




    vecX->_sourceStorageType = DataBuffer1D<scalar_inversion_storage>::StorageType::CPU;
    //vecX->dirtyData();
    _vecTemp->_sourceStorageType = DataBuffer1D<scalar_inversion_storage>::StorageType::CPU;
    //vecTemp->dirtyData();
    vecB->_sourceStorageType = DataBuffer1D<scalar_inversion_storage>::StorageType::CPU;
    //vecB->dirtyData();
    _basisFlowParams->_sourceStorageType = DataBuffer1D<BasisFlow>::StorageType::CPU;
    //basisFlowParams->dirtyData();
}


dvec2 eigenLaplace(dvec2 p, dvec2 k) {
    return dvec2(
        k.y*sin(M_PI*p.x*k.x)*cos(M_PI*p.y*k.y),
        -k.x*cos(M_PI*p.x*k.x)*sin(M_PI*p.y*k.y)
    );
}



// jacobian matrix
dmat2 gradEigenLaplace(dvec2 p, dvec2 k) {
    // glm is column-major
    return dmat2(
        k.x*k.y * M_PI * cos(k.x*M_PI*p.x) * cos(k.y*M_PI*p.y), // xdx
        k.x*k.x * M_PI * sin(k.x*M_PI*p.x) * sin(k.y*M_PI*p.y), // ydx
        -k.y*k.y * M_PI * sin(k.x*M_PI*p.x) * sin(k.y*M_PI*p.y), // xdy
        -k.x*k.y * M_PI * cos(k.x*M_PI*p.x) * cos(k.y*M_PI*p.y)  // ydy
    );
}




// integrates basis(...) dot a vector field defined on a grid
// TODO: pass basis info by reference when the function does not modify the basis? would avoid copying a lot of data...
float Application::IntegrateBasisGrid(BasisFlow& b, VectorField2D* velField)
{

    // compute intersection of supports
    BasisSupport supportBasis = b.getSupport();
    BasisSupport supportField = BasisSupport(_domainLeft, _domainRight, _domainBottom, _domainTop);
    float supLeft = glm::max(supportBasis.left, supportField.left);
    float supRight = glm::min(supportBasis.right, supportField.right);
    float supBottom = glm::max(supportBasis.bottom, supportField.bottom);
    float supTop = glm::min(supportBasis.top, supportField.top);

    if (supLeft >= supRight || supBottom >= supTop
        ) {
        return 0.f;
    }

    // compute integral as discretized sum at grid centers
#if INTEGRATION_SUM_DOUBLE_PRECISION
    double sum = 0;
#else
    float sum = 0;
#endif




#if DEF_COEFF_COMPUTE_GPU

#if INTEGRATE_BASIS_ONE_BASIS_PER_DISPATCH

    // data basis
    auto connBasis = createPullConnection(
        &basisFlowTemplates[abs(b.freqLvl.x - b.freqLvl.y)]->out_vectorsMetadataTexture2D,
        &pipelineIntegrateBasisGrid_onePerDispatch->in_texBasis);
    connBasis->activate();
    int minLvl = glm::min<int>(b.freqLvl.x, b.freqLvl.y);
    pipelineIntegrateBasisGrid_onePerDispatch->in_postScalingBasis.receive(float(1 << minLvl));
    pipelineIntegrateBasisGrid_onePerDispatch->in_centerBasis.receive(b.center);
    pipelineIntegrateBasisGrid_onePerDispatch->in_halfWidthsBasis.receive(0.5f*vec2(supportBasis.right - supportBasis.left, supportBasis.top - supportBasis.bottom));
    if (b.freqLvl.x <= b.freqLvl.y) {
        pipelineIntegrateBasisGrid_onePerDispatch->in_forwardRotBasis.receive(mat2(1, 0, 0, 1));
        pipelineIntegrateBasisGrid_onePerDispatch->in_backwardRotBasis.receive(mat2(1, 0, 0, 1));
    }
    else {
        pipelineIntegrateBasisGrid_onePerDispatch->in_forwardRotBasis.receive(mat2(0, 1, -1, 0));
        pipelineIntegrateBasisGrid_onePerDispatch->in_backwardRotBasis.receive(mat2(0, -1, 1, 0));
    }

    // data field
    auto connField = createPullConnection(
        &velField->out_vectorsMetadataTexture2D,
        &pipelineIntegrateBasisGrid_onePerDispatch->in_texField);
    connField->activate();
    pipelineIntegrateBasisGrid_onePerDispatch->in_centerField.receive(domainCenter);
    pipelineIntegrateBasisGrid_onePerDispatch->in_halfWidthsField.receive(0.5f*vec2(supportField.right - supportField.left, supportField.top - supportField.bottom));

    // other data
    auto connTransfer = createPushConnection(
        &pipelineIntegrateBasisGrid_onePerDispatch->out_transferBuffer,
        &integrationTransferBufferGpu->in_metadataBuffer);
    connTransfer->activate();
    pipelineIntegrateBasisGrid_onePerDispatch->in_integralGridRes.receive(integralGridRes);
    pipelineIntegrateBasisGrid_onePerDispatch->in_supportInterLeftBottom.receive(vec2(supLeft, supBottom));
    pipelineIntegrateBasisGrid_onePerDispatch->in_supportInterRightTop.receive(vec2(supRight, supTop));

    unsigned int nbGroupDiv = 1;
    unsigned int it = 0;
    uvec2 nbGroups;
    do
    {
        nbGroups = uvec2(((integralGridRes + 1) - 1) / nbGroupDiv / INTEGRAL_GPU_GROUP_DIM + 1,
            ((integralGridRes + 1) - 1) / nbGroupDiv / INTEGRAL_GPU_GROUP_DIM + 1);

        pipelineIntegrateBasisGrid_onePerDispatch->in_globalIteration.receive(it++);
        pipelineIntegrateBasisGrid_onePerDispatch->goglu_nbWorkGroups.set(glm::uvec3(nbGroups.x, nbGroups.y, 1));
        pipelineIntegrateBasisGrid_onePerDispatch->execute();
        glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

        nbGroupDiv *= INTEGRAL_GPU_GROUP_DIM;

    } while (nbGroups.x > 1 && nbGroups.y > 1);

    sum = integrationTransferBufferGpu->getCpuData(0).x;

    connBasis->deactivate(); delete connBasis;
    connField->deactivate(); delete connField;
    connTransfer->deactivate(); delete connTransfer;

#elif INTEGRATE_BASIS_ONE_BASIS_PER_INVOCATION
    // SHOULD NEVER REACH THIS
    velField; // to remove warning
    cout << "********" << endl;
    cout << "********" << endl;
    cout << "********" << endl;
    cout << "********" << endl;
    cout << "SHOULD NOT REACH THIS" << endl;
    cout << "********" << endl;
    cout << "********" << endl;
    cout << "********" << endl;
    cout << "********" << endl;
#endif

#else

    for (uint i = 0; i <= _integralGridRes; i++) {
        for (uint j = 0; j <= _integralGridRes; j++) {

            vec2 p = vec2(supLeft + float(i) / _integralGridRes * (supRight - supLeft),
                supBottom + float(j) / _integralGridRes * (supTop - supBottom));

            sum += ((i == 0 || i == _integralGridRes) ? 0.5f : 1.f) * ((j == 0 || j == _integralGridRes) ? 0.5f : 1.f) *
                glm::dot(TranslatedBasisEval(p, b.freqLvl, b.center),
                    velField->interp(p));

        }
    }

#endif

    return float(sum) * (supRight - supLeft)*(supTop - supBottom) / Sqr(_integralGridRes);
}






// integrates basis(...) dot basis(...)
float Application::IntegrateBasisBasis(BasisFlow b1, BasisFlow b2) {

    BasisSupport sup1 = b1.getSupport();
    BasisSupport sup2 = b2.getSupport();
    float supLeft = glm::max(sup1.left, sup2.left);
    float supRight = glm::min(sup1.right, sup2.right);
    float supBottom = glm::max(sup1.bottom, sup2.bottom);
    float supTop = glm::min(sup1.top, sup2.top);

    if (supLeft >= supRight || supBottom >= supTop
        ) {
        return 0.f;
    }


    // compute integral as discretized sum at grid centers
#if INTEGRATION_SUM_DOUBLE_PRECISION
    double sum = 0;
#else
    float sum = 0;
#endif



    for (int i = 0; i <= int(_integralGridRes); i++) {
        for (int j = 0; j <= int(_integralGridRes); j++) {

            vec2 p = vec2(supLeft + float(i) / _integralGridRes * (supRight - supLeft),
                supBottom + float(j) / _integralGridRes * (supTop - supBottom));
            sum += ((i == 0 || i == _integralGridRes) ? 0.5f : 1.f) * ((j == 0 || j == _integralGridRes) ? 0.5f : 1.f) *
                glm::dot(TranslatedBasisEval(p, b1.freqLvl, b1.center),
                    TranslatedBasisEval(p, b2.freqLvl, b2.center));
        }
    }


    return float(sum) * (supRight - supLeft)*(supTop - supBottom) / Sqr(_integralGridRes);


}




#if EXPLICIT_TRANSPORT_ROTATION

// gives the average value of bVec over the support of bSupport
vec2 Application::AverageBasisOnSupport(BasisFlow bVec, BasisFlow bSupport) {

    // compute intersection of supports
    BasisSupport supportVec = bVec.getSupport();
    BasisSupport supportSup = bSupport.getSupport();
    float supLeft = glm::max(supportVec.left, supportSup.left);
    float supRight = glm::min(supportVec.right, supportSup.right);
    float supBottom = glm::max(supportVec.bottom, supportSup.bottom);
    float supTop = glm::min(supportVec.top, supportSup.top);

    if (supLeft >= supRight || supBottom >= supTop
        ) {
        return vec2(0);
    }

    // compute integral as discretized sum at grid centers
#if INTEGRATION_SUM_DOUBLE_PRECISION
    dvec2 sum(0);
#else
    vec2 sum(0);
#endif

#if DEF_COEFF_COMPUTE_GPU

    int minLvl;

    // data basis
    auto connI = createPullConnection(
        &basisFlowTemplates[abs(bVec.freqLvl.x - bVec.freqLvl.y)]->out_vectorsMetadataTexture2D,
        &pipelineIntegrateAvgBasis->in_texBi);
    connI->activate();
    minLvl = glm::min<int>(bVec.freqLvl.x, bVec.freqLvl.y);
    pipelineIntegrateAvgBasis->in_postScalingI.receive(float(1 << minLvl));
    pipelineIntegrateAvgBasis->in_centerI.receive(bVec.center);
    pipelineIntegrateAvgBasis->in_halfWidthsI.receive(0.5f*vec2(supportVec.right - supportVec.left, supportVec.top - supportVec.bottom));
    if (bVec.freqLvl.x <= bVec.freqLvl.y) {
        pipelineIntegrateAvgBasis->in_forwardRotI.receive(mat2(1, 0, 0, 1));
        pipelineIntegrateAvgBasis->in_backwardRotI.receive(mat2(1, 0, 0, 1));
    }
    else {
        pipelineIntegrateAvgBasis->in_forwardRotI.receive(mat2(0, 1, -1, 0));
        pipelineIntegrateAvgBasis->in_backwardRotI.receive(mat2(0, -1, 1, 0));
    }

    // other data
    auto connTransfer = createPushConnection(
        &pipelineIntegrateAvgBasis->out_transferBuffer,
        &integrationTransferBufferGpu->in_metadataBuffer);
    connTransfer->activate();
    pipelineIntegrateAvgBasis->in_integralGridRes.receive(integralGridRes);
    pipelineIntegrateAvgBasis->in_supportInterLeftBottom.receive(vec2(supLeft, supBottom));
    pipelineIntegrateAvgBasis->in_supportInterRightTop.receive(vec2(supRight, supTop));


    unsigned int nbGroupDiv = 1;
    unsigned int it = 0;
    uvec2 nbGroups;
    do
    {
        nbGroups = uvec2(((integralGridRes + 1) - 1) / nbGroupDiv / INTEGRAL_GPU_GROUP_DIM + 1,
            ((integralGridRes + 1) - 1) / nbGroupDiv / INTEGRAL_GPU_GROUP_DIM + 1);

        pipelineIntegrateAvgBasis->in_globalIteration.receive(it++);
        pipelineIntegrateAvgBasis->goglu_nbWorkGroups.set(glm::uvec3(nbGroups.x, nbGroups.y, 1));
        pipelineIntegrateAvgBasis->execute();
        glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

        nbGroupDiv *= INTEGRAL_GPU_GROUP_DIM;

    } while (nbGroups.x > 1 && nbGroups.y > 1);


    sum = vec2(integrationTransferBufferGpu->getCpuData(0)); // vec4 -> vec2

    connI->deactivate(); delete connI;
    connTransfer->deactivate(); delete connTransfer;

#else

    for (uint i = 0; i <= _integralGridRes; i++) {
        for (uint j = 0; j <= _integralGridRes; j++) {

            vec2 p = vec2(supLeft + float(i) / _integralGridRes * (supRight - supLeft),
                supBottom + float(j) / _integralGridRes * (supTop - supBottom));




            sum += ((i == 0 || i == _integralGridRes) ? 0.5f : 1.f) * ((j == 0 || j == _integralGridRes) ? 0.5f : 1.f) *
                TranslatedBasisEval(p, bVec.freqLvl, bVec.center);
        }
    }

#endif

    // divide by domain size to get averaged value
    return vec2(sum) * (supRight - supLeft)*(supTop - supBottom) / float(Sqr(_integralGridRes)) / ((1.f / float(1 << bSupport.freqLvl.x))*(1.f / float(1 << bSupport.freqLvl.y)));

}

#endif



void Application::SaveCoeffsBB(string filename)
{
    if (_newBBCoeffComputed) {
        ofstream file;
        file.open(filename);
        for (auto coeffPairIt = _coeffsBB.begin(); coeffPairIt != _coeffsBB.end(); coeffPairIt++) {
            KeyTypeBB key = coeffPairIt->first;
            float val = coeffPairIt->second;
            file << get<0>(key) << " " <<
                get<1>(key) << " " <<
                get<2>(key) << " " <<
                get<3>(key) << " " <<
                get<4>(key) << " " <<
                get<5>(key) << " " <<
                val << endl;
        }
        file.close();
        std::cout << "saved BB coefficient to " << filename << endl;
    }
}


void Application::LoadCoeffsBB(string filename)
{
    ifstream file;
    file.open(filename);
    string line;
    stringstream ss;
    while (getline(file, line)) {
        int dataUI[4];
        float dataF[3];
        ss.clear();
        ss.str(line);
        ss >> dataUI[0] >>
            dataUI[1] >>
            dataUI[2] >>
            dataUI[3] >>
            dataF[0] >>
            dataF[1] >>
            dataF[2];
        KeyTypeBB key = make_tuple(
            dataUI[0],
            dataUI[1],
            dataUI[2],
            dataUI[3],
            RoundToMultiple(dataF[0], _coeffSnapSize),
            RoundToMultiple(dataF[1], _coeffSnapSize)
        );
        float coeff = dataF[2];
        _coeffsBB.insert(std::pair<KeyTypeBB, float>(key, coeff));
    }
    file.close();
}



void Application::SaveCoeffsT(string filename)
{
    if (_newTCoeffComputed) {
        ofstream file;
        file.open(filename);
        for (auto coeffPairIt = _coeffsT.begin(); coeffPairIt != _coeffsT.end(); coeffPairIt++) {
            KeyTypeT key = coeffPairIt->first;
            vec2 val = coeffPairIt->second;
            file << get<0>(key) << " " <<
                get<1>(key) << " " <<
                get<2>(key) << " " <<
                get<3>(key) << " " <<
                get<4>(key) << " " <<
                get<5>(key) << " " <<
                val.x << " " <<
                val.y << endl;
        }
        file.close();
        std::cout << "saved R coefficient to " << filename << endl;
    }
}


void Application::LoadCoeffsT(string filename)
{
    ifstream file;
    file.open(filename);
    string line;
    stringstream ss;
    while (getline(file, line)) {
        int dataUI[4];
        float dataF[4];
        ss.clear();
        ss.str(line);
        ss >> dataUI[0] >>
            dataUI[1] >>
            dataUI[2] >>
            dataUI[3] >>
            dataF[0] >>
            dataF[1] >>
            dataF[2] >>
            dataF[3];
        KeyTypeT key = make_tuple(
            dataUI[0],
            dataUI[1],
            dataUI[2],
            dataUI[3],
            RoundToMultiple(dataF[0], _coeffSnapSize),
            RoundToMultiple(dataF[1], _coeffSnapSize)
        );
        vec2 coeff(dataF[2], dataF[3]);
        _coeffsT.insert(std::pair<KeyTypeT, vec2>(key, coeff));
    }
    file.close();
}




// TODO: optimize this using const&
vec2 Application::MatTCoeff(int i, int j) {
    BasisFlow bTransported = _basisFlowParams->getCpuData(i);
    BasisFlow bTransporting = _basisFlowParams->getCpuData(j);
    return MatTCoeff(bTransported, bTransporting);
}


// compute the integral and store it for future use
vec2 Application::MatTCoeff(BasisFlow bTransported, BasisFlow bTransporting)
{

    if (IntersectionInteriorEmpty(bTransported.getSupport(), bTransporting.getSupport())) { return vec2(0); }


    int baseLvl;
    baseLvl = glm::min<int>(
        glm::min<int>(bTransported.freqLvl.x, bTransported.freqLvl.y),
        glm::min<int>(bTransporting.freqLvl.x, bTransporting.freqLvl.y)
        );
    float baseFreq = powf(2.f, float(baseLvl));

    // remove common frequency factors
    ivec2 normFreqLvlTransported = bTransported.freqLvl - baseLvl;
    ivec2 normFreqLvlTransporting = bTransporting.freqLvl - baseLvl;

    vec2 relativeOffset = baseFreq * vec2(
        bTransported.center.x - bTransporting.center.x,
        bTransported.center.y - bTransporting.center.y
    );

    // snap offset to very fine grid to avoid float errors
    vec2 snappedRelativeOffset(RoundToMultiple(relativeOffset.x, _coeffSnapSize),
        RoundToMultiple(relativeOffset.y, _coeffSnapSize)
    );

    // return coefficient if already computed, or else compute it and store it.
    KeyTypeT key = make_tuple(normFreqLvlTransported.x, normFreqLvlTransported.y,
        normFreqLvlTransporting.x, normFreqLvlTransporting.y,
        snappedRelativeOffset.x, snappedRelativeOffset.y);
    auto coeffPair = _coeffsT.find(key);
    vec2 result;

    if (coeffPair != _coeffsT.end()) {
        result = coeffPair->second;
    }
    else {
        vec2 coeff;

        // compute coefficient with numerical integration
        BasisFlow bRelativeTransporting(normFreqLvlTransporting, vec2(0));
        BasisFlow bRelativeTransported(normFreqLvlTransported, relativeOffset);

        coeff = AverageBasisOnSupport(bRelativeTransporting, bRelativeTransported);

        _coeffsT.insert(std::pair<KeyTypeT, vec2>(key, coeff));
        result = coeff;

        if (_coeffsT.size() % 1000 == 0) {
            std::cout << "coeffs T : " << _coeffsT.size() << endl;
        }

        _newTCoeffComputed = true;
    }


    // scaled coefficient
    return result * baseFreq;

}








//TODO: optimize this using const&
float Application::MatBBCoeff(int i, int j) {
    BasisFlow b1 = _basisFlowParams->getCpuData(i);
    BasisFlow b2 = _basisFlowParams->getCpuData(j);
    return MatBBCoeff(b1, b2);
}



// compute the integral and store it for future use
float Application::MatBBCoeff(const BasisFlow& b1, const BasisFlow& b2)
{
    if (IntersectionInteriorEmpty(b1.getSupport(), b2.getSupport())) { return 0.0; }

#if ENFORCE_EXACT_ORTHOGONALITY
    if (b1.orthoGroup == b2.orthoGroup && b1.center != b2.center) {
        return 0;
    }
#endif

    int baseLvl;
    baseLvl = glm::min<int>(
        glm::min<int>(b1.freqLvl.x, b1.freqLvl.y),
        glm::min<int>(b2.freqLvl.x, b2.freqLvl.y)
        );

    float baseFreq = powf(2.f, float(baseLvl));

    // remove common frequency factors
    ivec2 normFreqLvl1 = b1.freqLvl - baseLvl;
    ivec2 normFreqLvl2 = b2.freqLvl - baseLvl;

    // all bases are symmetric w.r.t. x and y, so we take the relative offset in the first quadrant.
//    vec2 relativeOffset = baseFreq*vec2( 
//                abs(b2.center.x-b1.center.x),
//                abs(b2.center.y-b1.center.y) 
//                #if SIM3D
//                    ,abs(b2.center.z-b1.center.z)
//                #endif
//                );
    vec2 relativeOffset = baseFreq * vec2(
        b2.center.x - b1.center.x,
        b2.center.y - b1.center.y
    );



    // snap offset to very fine grid to avoid float errors
    vec2 snappedRelativeOffset = vec2(RoundToMultiple(relativeOffset.x, _coeffSnapSize),
        RoundToMultiple(relativeOffset.y, _coeffSnapSize)
    );

    // return coefficient if already computed, or else compute it and store it.
    KeyTypeBB key = make_tuple(normFreqLvl1.x, normFreqLvl1.y,
        normFreqLvl2.x, normFreqLvl2.y,
        snappedRelativeOffset.x, snappedRelativeOffset.y);
    auto coeffPair = _coeffsBB.find(key);
    float result;

    if (coeffPair != _coeffsBB.end()) {
        result = coeffPair->second;
    }
    else {
        float coeff;

        // compute coefficient with numerical integration
        BasisFlow bRelative1(normFreqLvl1, vec2(0));
        BasisFlow bRelative2(normFreqLvl2, relativeOffset);
        coeff = float(IntegrateBasisBasis(bRelative1, bRelative2));

        _coeffsBB.insert(std::pair<KeyTypeBB, float>(key, coeff));
        result = coeff;

        if (_coeffsBB.size() % 1000 == 0) {
            std::cout << "coeffs BB : " << _coeffsBB.size() << endl;
        }

        _newBBCoeffComputed = true;
    }

    // scaled coefficient
    return result;

}





// using harmonics 1-3-5, see frequencyExploration_cleaned.nb. It's the
// \widehat{h} from the paper. We include the normalization parameter in the
// definition, and only compute cases where kx = 1, so ky = 2^log2Aniso.
// Oriented along Y axis
dvec2 flowBasisHat(dvec2 p, int log2Aniso)
{
    // frequencies are 2^level
    int kx = 1;
    int ky = 1 << log2Aniso;

#if SAFETY_ASSERTS
    // out of bound asserts
    if (
        log2Aniso > 2 ||
        !isInClosedInterval(p.x, -0.5 / kx, 0.5 / kx) ||
        !isInClosedInterval(p.y, -0.5 / ky, 0.5 / ky)
        ) {
        Application::spApp->insertDebugMessage("flowBasisHat : out of bounds.");
        return dvec2(0, 0);
    }
#endif

    double coeffs[3][3];
    double norm;

    switch (log2Aniso) {
    case 0:
        coeffs[0][0] = 1.;
        coeffs[0][1] = -0.1107231462697129;
        coeffs[0][2] = -0.1335661122381723;
        coeffs[1][0] = -0.1107231462697125;
        coeffs[1][1] = 0.126276763526633;
        coeffs[1][2] = -0.05362142886203727;
        coeffs[2][0] = -0.1335661122381725;
        coeffs[2][1] = -0.05362142886203719;
        coeffs[2][2] = 0.05888607976485681;
        norm = 1.0221139695997405;
        break;
    case 1:
        coeffs[0][0] = 1.;
        coeffs[0][1] = 0.02773519585551282;
        coeffs[0][2] = -0.2166411175133077;
        coeffs[1][0] = -0.4866818264236261;
        coeffs[1][1] = 0.05437868400363529;
        coeffs[1][2] = 0.06470915488254406;
        coeffs[2][0] = 0.0920090958541756;
        coeffs[2][1] = -0.03817424957328374;
        coeffs[2][2] = 0.00450273057313512;
        norm = 0.7620965477955399;
        break;
    case 2:
        coeffs[0][0] = 1.;
        coeffs[0][1] = 0.03365588438312442;
        coeffs[0][2] = -0.2201935306298747;
        coeffs[1][0] = -0.5578126028758348;
        coeffs[1][1] = -0.00367012134209574;
        coeffs[1][2] = 0.1137645933804244;
        coeffs[2][0] = 0.1346875617255008;
        coeffs[2][1] = -0.004529104071367439;
        coeffs[2][2] = -0.02422004990227971;
        norm = 0.5618900800300474;
        break;
    default:
        std::cout << "unknown basis parameters" << endl;
        return dvec2(0, 0);
        break;
    }


    dvec2 p2(p.x + 0.5 / kx, p.y + 0.5 / ky);

    return norm * (
        coeffs[0][0] * eigenLaplace(p2, dvec2(1 * kx, 1 * ky)) +
        coeffs[0][1] * eigenLaplace(p2, dvec2(1 * kx, 3 * ky)) +
        coeffs[0][2] * eigenLaplace(p2, dvec2(1 * kx, 5 * ky)) +
        coeffs[1][0] * eigenLaplace(p2, dvec2(3 * kx, 1 * ky)) +
        coeffs[1][1] * eigenLaplace(p2, dvec2(3 * kx, 3 * ky)) +
        coeffs[1][2] * eigenLaplace(p2, dvec2(3 * kx, 5 * ky)) +
        coeffs[2][0] * eigenLaplace(p2, dvec2(5 * kx, 1 * ky)) +
        coeffs[2][1] * eigenLaplace(p2, dvec2(5 * kx, 3 * ky)) +
        coeffs[2][2] * eigenLaplace(p2, dvec2(5 * kx, 5 * ky))
        );

}

dmat2 flowBasisHatGrad(dvec2 p, int log2Aniso)
{
    // frequencies are 2^level
    int kx = 1;
    int ky = 1 << log2Aniso;

#if SAFETY_ASSERTS
    // out of bound asserts
    if (
        log2Aniso > 2 ||
        !isInClosedInterval(p.x, -0.5 / kx, 0.5 / kx) ||
        !isInClosedInterval(p.y, -0.5 / ky, 0.5 / ky)
        ) {
        Application::spApp->insertDebugMessage("flowBasisHatGrad : out of bounds.");
        return dmat2(0, 0, 0, 0);
    }
#endif

    double coeffs[3][3];
    double norm; // actually 1/norm? Like,the factor you need to multiply by to get unit norm ?

    switch (log2Aniso) {
    case 0:
        coeffs[0][0] = 1.;
        coeffs[0][1] = -0.1107231462697129;
        coeffs[0][2] = -0.1335661122381723;
        coeffs[1][0] = -0.1107231462697125;
        coeffs[1][1] = 0.126276763526633;
        coeffs[1][2] = -0.05362142886203727;
        coeffs[2][0] = -0.1335661122381725;
        coeffs[2][1] = -0.05362142886203719;
        coeffs[2][2] = 0.05888607976485681;
        norm = 1.0221139695997405;
        break;
    case 1:
        coeffs[0][0] = 1.;
        coeffs[0][1] = 0.02773519585551282;
        coeffs[0][2] = -0.2166411175133077;
        coeffs[1][0] = -0.4866818264236261;
        coeffs[1][1] = 0.05437868400363529;
        coeffs[1][2] = 0.06470915488254406;
        coeffs[2][0] = 0.0920090958541756;
        coeffs[2][1] = -0.03817424957328374;
        coeffs[2][2] = 0.00450273057313512;
        norm = 0.7620965477955399;
        break;
    case 2:
        coeffs[0][0] = 1.;
        coeffs[0][1] = 0.03365588438312442;
        coeffs[0][2] = -0.2201935306298747;
        coeffs[1][0] = -0.5578126028758348;
        coeffs[1][1] = -0.00367012134209574;
        coeffs[1][2] = 0.1137645933804244;
        coeffs[2][0] = 0.1346875617255008;
        coeffs[2][1] = -0.004529104071367439;
        coeffs[2][2] = -0.02422004990227971;
        norm = 0.5618900800300474;
        break;
    default:
        std::cout << "unknown basis parameters" << endl;
        return dmat2(0, 0, 0, 0);
        break;
    }

    dvec2 p2(p.x + 0.5 / kx, p.y + 0.5 / ky);

    return norm * (
        coeffs[0][0] * gradEigenLaplace(p2, dvec2(1 * kx, 1 * ky)) +
        coeffs[0][1] * gradEigenLaplace(p2, dvec2(1 * kx, 3 * ky)) +
        coeffs[0][2] * gradEigenLaplace(p2, dvec2(1 * kx, 5 * ky)) +
        coeffs[1][0] * gradEigenLaplace(p2, dvec2(3 * kx, 1 * ky)) +
        coeffs[1][1] * gradEigenLaplace(p2, dvec2(3 * kx, 3 * ky)) +
        coeffs[1][2] * gradEigenLaplace(p2, dvec2(3 * kx, 5 * ky)) +
        coeffs[2][0] * gradEigenLaplace(p2, dvec2(5 * kx, 1 * ky)) +
        coeffs[2][1] * gradEigenLaplace(p2, dvec2(5 * kx, 3 * ky)) +
        coeffs[2][2] * gradEigenLaplace(p2, dvec2(5 * kx, 5 * ky))
        );
}

//evaluate basis from stored basis templates. Note that in exponent space,
//division is substraction, so lvlY-minLvl is ky/minK
vec2 Application::TranslatedBasisEval(
    const vec2 p,
    const ivec2 freqLvl,
    const vec2 center)
{
#if SAFETY_ASSERTS
    if (abs(int(freqLvl.x) - int(freqLvl.y)) > maxAnisoLvl) { cout << "translateBasisEval: frequency out of bounds." << endl; return vec2(0, 0); }
#endif

    vec2 result;

    int minLvl = glm::min<int>(freqLvl.x, freqLvl.y);
    if (freqLvl.x <= freqLvl.y) {
        result = _basisFlowTemplates[freqLvl.y - freqLvl.x]->interp(
            vec2(float(1 << minLvl)*(p.x - center.x) / _lengthLvl0,
                float(1 << minLvl)*(p.y - center.y) / _lengthLvl0)
        );
    }
    else {
        // reverse coordinates
        result = -_basisFlowTemplates[freqLvl.x - freqLvl.y]->interp(
            vec2((1 << minLvl)*(p.y - center.y) / _lengthLvl0,
            (1 << minLvl)*(p.x - center.x) / _lengthLvl0)
        );
        result = vec2(result.y, result.x);
    }
    return float(1 << minLvl)*result;
}


// evaluates without using stored grid, using double precision.
dvec2 Application::TranslatedBasisEvalPrecise(const dvec2 p, const ivec2 freqLvl, const dvec2 center)
{
    dvec2 result;
    int minLvl = glm::min<int>(freqLvl.x, freqLvl.y);

    if (freqLvl.x <= freqLvl.y) {
        result = flowBasisHat(
            double(1 << minLvl)*(p - center) / double(_lengthLvl0),
            freqLvl.y - minLvl
        );
    }
    else {
        // rotated
        dvec2 r = p - center;
        result = flowBasisHat(
            double(1 << minLvl)*dvec2(r.y, -r.x) / double(_lengthLvl0),
            freqLvl.x - minLvl
        );
        result = dvec2(-result.y, result.x);
    }
    return double(1 << minLvl)*result;
}



mat2 Application::TranslatedBasisGradEval(
    const vec2 p,
    const ivec2 freqLvl,
    const vec2 center)
{
#if USE_PRECISE_BASIS_EVAL
    return FMATD(translatedBasisGradEvalPrecise(dvec2(p), freqLvl, dvec2(center)));
#else
    //vec2 eps = 1.f / (float(_nbCells) * vec2(1 << freqLvl) * 2.0f);
    vec2 eps = 1.f / (float(BASE_GRID_SIZE - 1) * vec2(1 << freqLvl) * 2.0f);

    return mat2(
        (TranslatedBasisEval(p + vec2(eps.x, 0), freqLvl, center) - TranslatedBasisEval(p - vec2(eps.x, 0), freqLvl, center)) / (2.f*eps.x),
        (TranslatedBasisEval(p + vec2(0, eps.y), freqLvl, center) - TranslatedBasisEval(p - vec2(0, eps.y), freqLvl, center)) / (2.f*eps.y)
    );
#endif

}



dmat2 Application::TranslatedBasisGradEvalPrecise(const dvec2 p, const ivec2 freqLvl, const dvec2 center)
{
    dmat2 result;
    int minLvl = glm::min<int>(freqLvl.x, freqLvl.y);
    if (freqLvl.x <= freqLvl.y) {
        result = flowBasisHatGrad(
            double(1 << minLvl)*(p - center) / double(_lengthLvl0),
            freqLvl.y - minLvl
        );
    }
    else {
        // rotated
        dvec2 r = p - center;
        result = -flowBasisHatGrad(
            double(1 << minLvl)*dvec2(r.y, r.x) / double(_lengthLvl0),
            freqLvl.x - minLvl
        );
        result = dmat2(result[1][1], result[1][0], result[0][1], result[0][0]);
        //        result = -flowBasisHatGrad(
        //                    double(1<<minLvl)*(p.y-center.y)/lengthLvl0,
        //                    double(1<<minLvl)*(p.x-center.x)/lengthLvl0,
        //                    freqLvl.x-minLvl
        //                );
    }

    return double(Sqr(1 << minLvl))*result;
}


//#endif



BasisSupport BasisFlow::getSupport() const
{
    vec2 halfSize(
        0.5f / float(1 << freqLvl.x),
        0.5f / float(1 << freqLvl.y)
    );
    return BasisSupport(
        center.x - halfSize.x,
        center.x + halfSize.x,
        center.y - halfSize.y,
        center.y + halfSize.y
    );
}



vec2 BasisFlow::supportHalfSize() const
{
    vec2 freq(1 << freqLvl.x, 1 << freqLvl.y);
    return 0.5f / freq;
}



vec2 BasisFlow::normalizedPositionInSupport(vec2 p) {
    return (p - center) / (2.f*supportHalfSize());
}


bool BasisFlow::pointIsInSupport(vec2 p) {
    vec2 normalizedPos = normalizedPositionInSupport(p);
    return IsInClosedInterval(normalizedPos.x, -0.5f, 0.5f) &&
        IsInClosedInterval(normalizedPos.y, -0.5f, 0.5f);
}


bool IntersectionInteriorEmpty(BasisSupport& sup1, BasisSupport& sup2)
{
    return (glm::max(sup1.left, sup2.left) >= glm::min(sup1.right, sup2.right))
        || (glm::max(sup1.bottom, sup2.bottom) >= glm::min(sup1.top, sup2.top));
}



