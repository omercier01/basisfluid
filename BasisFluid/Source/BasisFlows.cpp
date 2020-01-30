//#include "glm/glm.hpp"

#include "BasisFlows.h"

#include <set>
#include <iostream>

using namespace glm;
using namespace std;

template <class T>
T min3(T a, T b, T c) { return glm::min(a, glm::min(b, c)); }

template <class T>
T max3(T a, T b, T c) { return glm::max(a, glm::max(b, c)); }

template <class T>
int minVec(T v) { return glm::min<int>(v.x, v.y); }


bool intersectionInteriorEmpty(BasisSupport& sup1, BasisSupport& sup2)
{
    return (glm::max(sup1.left, sup2.left) >= glm::min(sup1.right, sup2.right))
        || (glm::max(sup1.bottom, sup2.bottom) >= glm::min(sup1.top, sup2.top))
#if SIM3D
        || (glm::max(sup1.back, sup2.back) >= glm::min(sup1.front, sup2.front))
#endif
        ;
}


bool BasisFlow::emptyIntersectionWithBasis(BasisFlow& b) const
{
    return intersectionInteriorEmpty(getSupport(), b.getSupport());
}


scalar_inversion_inner Application::inverseBBMatrixMain(unsigned int iRow, scalar_inversion_storage* vecX, scalar_inversion_storage* vecB, BasisFlow* basisDataPointer, unsigned int basisBitMask, float alpha)
{
    scalar_inversion_inner xSquaredDiff = 0;
    scalar_inversion_inner tempX = scalar_inversion_inner(vecB[iRow]);

#if USE_DECOMPRESSED_COEFFICIENTS
    vector<CoeffBBDecompressedIntersectionInfo>& intersectionInfos = coeffsBBDecompressedIntersections[iRow];
    for (const CoeffBBDecompressedIntersectionInfo& inter : intersectionInfos) {
        //            if( atLeastOneBitNotSet(basisDataPointer[inter.j].bitFlags, basisBitMask) ) {
        //                continue;
        //            }
        //            tempX -= inter.coeff * scalar_inversion_inner(vecX[inter.j]);
        if (allBitsSet(basisDataPointer[inter.j].bitFlags, basisBitMask)) {
            tempX -= inter.coeff * scalar_inversion_inner(vecX[inter.j]);
        }

    }

    // TODO: make all bases normalized in 3D so we can save a fetch + division here
    scalar_inversion_storage newX = scalar_inversion_storage(tempX / scalar_inversion_inner(basisDataPointer[iRow].normSquared));
    scalar_inversion_storage oldX = vecX[iRow];

    xSquaredDiff = sqr(scalar_inversion_inner(newX - oldX));
    vecX[iRow] = (1 - alpha)*oldX + alpha * newX;
#else
    vector<uint>* localIntersectingBasesIds = intersectingBasesSignificantBBIds->getCpuData_noRefresh(iRow);
    for (auto it = localIntersectingBasesIds->begin(); it != localIntersectingBasesIds->end(); ++it) {
        if (*it == iRow) { continue; }
        if (atLeastOneBitNotSet(basisDataPointer[*it].bitFlags, basisBitMask)) {
            continue;
        }
        tempX -= matBBCoeff(iRow, (*it)) * scalar_inversion_inner(vecX[*it]);
    }
    // TODO: make all bases normalized in 3D so we can save a fetch + division here
    scalar_inversion_storage newX = scalar_inversion_storage(tempX / scalar_inversion_inner(matBBCoeff(iRow, iRow)));
    scalar_inversion_storage oldX = vecX[iRow];

    xSquaredDiff = sqr(scalar_inversion_inner(newX - oldX));
    vecX[iRow] = (1 - alpha)*oldX + alpha * newX;
#endif
    return scalar_inversion_inner(xSquaredDiff);
}



// Solves BB * x = b using Gauss-Seidel.
// Since bases are ordered by wavenumber, doing only one step of GS is equivalent to projection by frequency layer.
// The bitMask indicates the bits that must be set for the basis to be used in the inversion.
void InverseBBMatrix(
    DataBuffer1D<scalar_inversion_storage>* vecX,
    DataBuffer1D<scalar_inversion_storage>* vecB,
    float tol, unsigned int basisBitMask, unsigned int minFreq)
{


    uint maxNbItMatBBInversion = tw->maxNbItMatBBInversion.get();
    uint n = vecX->_nbElements;
    float alpha = tw->relaxationAlpha.get();

    // get references
    vecX->refreshCpuData();
    scalar_inversion_storage* vecXPointer = vecX->getCpuDataPointer();
    vecTemp->refreshCpuData();
    scalar_inversion_storage* vecTempPointer = vecTemp->getCpuDataPointer();
    vecB->refreshCpuData();
    scalar_inversion_storage* vecBPointer = vecB->getCpuDataPointer();
    basisFlowParams->refreshCpuData();
    BasisFlow* basisFlowParamsPointer = basisFlowParams->getCpuDataPointer();
    intersectingBasesSignificantBBIds->refreshCpuData();

    if (tw->bInversionGaussSeidel.get())
    {

        // zero x
        //#pragma omp parallel for
        for (int i = 0; i < int(n); i++) {
            vecXPointer[i] = 0.;
        }

        scalar_inversion_inner xSquaredDiff = 999999;
        uint iIt = 0;
        double time1 = elapsedTime();
        if (maxNbItMatBBInversion == 0)
        {
            for (uint iRow = 0; iRow < n; iRow++) {
                if (
                    allBitsSet(basisFlowParamsPointer[iRow].bitFlags, basisBitMask) &&
                    basisFlowParamsPointer[iRow].freqLvl.x >= minFreq
                    && basisFlowParamsPointer[iRow].freqLvl.y >= minFreq
                    ) {
                    vecXPointer[iRow] = vecBPointer[iRow];
                }
            }
        }
        else {
            while (xSquaredDiff > tol && iIt < maxNbItMatBBInversion) {
                xSquaredDiff = 0;

#if INVERSION_OPENMP

                for (int iBasisGroup = 0; iBasisGroup < orthogonalBasisGroupIds.size(); iBasisGroup++)
                {
                    vector<unsigned int> ids = orthogonalBasisGroupIds[iBasisGroup];
                    if (ids.size() >= minSizeParallelInverse)
                    {
#pragma omp parallel for num_threads(7) reduction(+:xSquaredDiff) shared(vecXPointer,vecBPointer,ids,basisFlowParamsPointer) firstprivate(basisBitMask,alpha,minFreq) default(none)
                        for (int idid = 0; idid < ids.size(); idid++) {
                            int iRow = ids[idid];
                            if (
                                allBitsSet(basisFlowParamsPointer[iRow].bitFlags, basisBitMask) &&
                                basisFlowParamsPointer[iRow].freqLvl.x >= minFreq
                                && basisFlowParamsPointer[iRow].freqLvl.y >= minFreq
                                ) {
                                xSquaredDiff += inverseBBMatrixMain(iRow, vecXPointer, vecBPointer, basisFlowParamsPointer, basisBitMask, alpha);
                            }
                        }
                    }
                    else {
                        for (int idid = 0; idid < ids.size(); idid++) {
                            int iRow = ids[idid];
                            if (
                                allBitsSet(basisFlowParamsPointer[iRow].bitFlags, basisBitMask) &&
                                basisFlowParamsPointer[iRow].freqLvl.x >= minFreq
                                && basisFlowParamsPointer[iRow].freqLvl.y >= minFreq
                                ) {
                                xSquaredDiff += inverseBBMatrixMain(iRow, vecXPointer, vecBPointer, basisFlowParamsPointer, basisBitMask, alpha);
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
        double time2 = elapsedTime();
        if (tw->outputProfileText.get()) {
            cout << "time inverse inner loop : " << time2 - time1 << endl;
        }

        if (tw->outputProfileText.get() && tw->outputProfileText_details.get()) {
            cout << "Gauss-Seidel xSquaredDiff: " << xSquaredDiff << endl;
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

                if (atLeastOneBitNotSet(basisFlowParamsPointer[iRow].bitFlags, basisBitMask)) {
                    continue;
                }
                scalar_inversion_inner tempX = scalar_inversion_inner(vecBPointer[iRow]);

#if USE_DECOMPRESSED_COEFFICIENTS
                vector<CoeffBBDecompressedIntersectionInfo>& intersectionInfos = coeffsBBDecompressedIntersections[iRow];
                for (CoeffBBDecompressedIntersectionInfo& inter : intersectionInfos) {
                    if (allBitsSet(basisFlowParamsPointer[inter.j].bitFlags, basisBitMask)) {
                        tempX -= inter.coeff * scalar_inversion_inner(vecXPointer[inter.j]);
                    }
                }
                // TODO: remove division once 3D bases have norm 1
                BasisFlow bi = basisFlowParamsPointer[iRow];
                scalar_inversion_storage newX = scalar_inversion_storage(tempX / scalar_inversion_inner(bi.normSquared));
                xSquaredDiff += sqr(scalar_inversion_inner(newX - vecXPointer[iRow]));
                vecTempPointer[iRow] = newX;
#else
                vector<uint>* localIntersectingBasesIds = intersectingBasesSignificantBBIds->getCpuData_noRefresh(iRow);
                for (auto it = localIntersectingBasesIds->begin(); it != localIntersectingBasesIds->end(); ++it) {
                    if (*it == iRow) { continue; }
                    if (allBitsSet(basisFlowParamsPointer[*it].bitFlags, basisBitMask)) {
                        tempX -= matBBCoeff(iRow, (*it)) * scalar_inversion_inner(vecXPointer[*it]);
                    }
                }
                // TODO: remove division once 3D bases have norm 1
                scalar_inversion_storage newX = scalar_inversion_storage(tempX / scalar_inversion_inner(matBBCoeff(iRow, iRow)));
                xSquaredDiff += sqr(scalar_inversion_inner(newX - vecXPointer[iRow]));
                vecTempPointer[iRow] = newX;
#endif

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

        if (tw->outputProfileText.get() && tw->outputProfileText_details.get())
        {
            cout << "Jacobi xSquaredDiff: " << xSquaredDiff << endl;
        }

    }


    // compute residual in vecB
    if (tw->outputProfileText.get() && tw->outputProfileText_details.get()) {
        double bSquaredDiff = 0;
        for (uint iRow = 0; iRow < n; iRow++) {
            if (atLeastOneBitNotSet(basisFlowParamsPointer[iRow].bitFlags, basisBitMask)) {
                continue;
            }
            double tempError = vecBPointer[iRow];

#if USE_DECOMPRESSED_COEFFICIENTS
            vector<CoeffBBDecompressedIntersectionInfo>& intersectionInfos = coeffsBBDecompressedIntersections[iRow];
            for (CoeffBBDecompressedIntersectionInfo& inter : intersectionInfos) {
                if (allBitsSet(basisFlowParamsPointer[inter.j].bitFlags, basisBitMask)) {
                    tempError -= inter.coeff * vecXPointer[inter.j];
                }

            }
            // adding diagonal term
            tempError -= basisFlowParamsPointer[iRow].normSquared * vecXPointer[iRow];
            bSquaredDiff += sqr(tempError);
#else
            vector<uint>* localIntersectingBasesIds = intersectingBasesSignificantBBIds->getCpuData_noRefresh(iRow);
            for (auto it = localIntersectingBasesIds->begin(); it != localIntersectingBasesIds->end(); ++it) {
                if (*it == iRow) { continue; }
                if (allBitsSet(basisFlowParamsPointer[*it].bitFlags, basisBitMask)) {
                    tempError -= matBBCoeff(iRow, (*it)) * double(vecXPointer[*it]);
                }
            }
            bSquaredDiff += sqr(tempError);
#endif
        }
        if (tw->bInversionGaussSeidel.get()) {
            cout << "Gauss-Seidel bSquaredDiff: " << bSquaredDiff << endl;
        }
        else {
            cout << "Jacobi bSquaredDiff: " << bSquaredDiff << endl;
        }
    }


    vecX->sourceStorageType = DataBuffer1D<scalar_inversion_storage>::StorageType::CPU;
    vecX->dirtyData();
    vecTemp->sourceStorageType = DataBuffer1D<scalar_inversion_storage>::StorageType::CPU;
    vecTemp->dirtyData();
    vecB->sourceStorageType = DataBuffer1D<scalar_inversion_storage>::StorageType::CPU;
    vecB->dirtyData();
    basisFlowParams->sourceStorageType = DataBuffer1D<BasisFlow>::StorageType::CPU;
    basisFlowParams->dirtyData();

}


dvec2 eigenLaplace(dvec2 p, dvec2 k) {
    return dvec2(
        k.y*sin(PI*p.x*k.x)*cos(PI*p.y*k.y),
        -k.x*cos(PI*p.x*k.x)*sin(PI*p.y*k.y)
    );
}



// jacobian matrix
#if SIM3D
    // Z-aligned
dmat3 gradEigenLaplace(dvec3 p, dvec3 k) {
    // glm is column-major
    return dmat3(
        sqr(k.x)*M_PI * k.z * cos(k.x*M_PI*p.x)*cos(k.y*M_PI*p.y)*cos(k.z*M_PI*p.z),            // xdx
        -k.x*M_PI * k.y * k.z * sin(k.x*M_PI*p.x)*sin(k.y*M_PI*p.y)*cos(k.z*M_PI*p.z),          // ydx
        k.x*M_PI*(sqr(k.x) + sqr(k.y)) * sin(k.x*M_PI*p.x)*cos(k.y*M_PI*p.y)*sin(k.z*M_PI*p.z), // zdx

        -k.y*M_PI*k.x * k.z * sin(k.x*M_PI*p.x)*sin(k.y*M_PI*p.y)*cos(k.z*M_PI*p.z),            // xdy
        k.y*M_PI*k.y * k.z * cos(k.x*M_PI*p.x)*cos(k.y*M_PI*p.y)*cos(k.z*M_PI*p.z),             // ydy
        k.y*M_PI*(sqr(k.x) + sqr(k.y)) * cos(k.x*M_PI*p.x)*sin(k.y*M_PI*p.y)*sin(k.z*M_PI*p.z), // zdy

        -k.z*M_PI*k.x * k.z * sin(k.x*M_PI*p.x)*cos(k.y*M_PI*p.y)*sin(k.z*M_PI*p.z),            // xdz
        -k.z*M_PI*k.y * k.z * cos(k.x*M_PI*p.x)*sin(k.y*M_PI*p.y)*sin(k.z*M_PI*p.z),            // ydz
        -k.z*M_PI*(sqr(k.x) + sqr(k.y)) * cos(k.x*M_PI*p.x)*cos(k.y*M_PI*p.y)*cos(k.z*M_PI*p.z) // zdz
    );
}
#else
dmat2 gradEigenLaplace(dvec2 p, dvec2 k) {
    // glm is column-major
    return dmat2(
        k.x*k.y * M_PI * cos(k.x*M_PI*p.x) * cos(k.y*M_PI*p.y), // xdx
        k.x*k.x * M_PI * sin(k.x*M_PI*p.x) * sin(k.y*M_PI*p.y), // ydx
        -k.y*k.y * M_PI * sin(k.x*M_PI*p.x) * sin(k.y*M_PI*p.y), // xdy
        -k.x*k.y * M_PI * cos(k.x*M_PI*p.x) * cos(k.y*M_PI*p.y)  // ydy
    );
}
#endif


// integrates basis(...) dot a vector field defined on a grid
// TODO: pass basis info by reference when the function does not modify the basis? would avoid copying a lot of data...
float IntegrateBasisGrid(BasisFlow& b, VectorField2D* velField)
{

    // compute intersection of supports
    BasisSupport supportBasis = b.getSupport();
#if SIM3D
    BasisSupport supportField = BasisSupport(domainLeft, domainRight, domainBottom, domainTop, domainBack, domainFront);
#else
    BasisSupport supportField = BasisSupport(domainLeft, domainRight, domainBottom, domainTop);
#endif
    float supLeft = glm::max(supportBasis.left, supportField.left);
    float supRight = glm::min(supportBasis.right, supportField.right);
    float supBottom = glm::max(supportBasis.bottom, supportField.bottom);
    float supTop = glm::min(supportBasis.top, supportField.top);
#if SIM3D
    float supBack = glm::max(supportBasis.back, supportField.back);
    float supFront = glm::min(supportBasis.front, supportField.front);
#endif

    if (supLeft >= supRight || supBottom >= supTop
#if SIM3D
        || supBack >= supFront
#endif
        ) {
        return 0.f;
    }

    // compute integral as discretized sum at grid centers
#if INTEGRATION_SUM_DOUBLE_PRECISION
    double sum = 0;
#else
    float sum = 0;
#endif



#if SIM3D

#if DEF_COEFF_COMPUTE_GPU

#if INTEGRATE_BASIS_ONE_BASIS_PER_DISPATCH
    int baseLvl, log2Aniso, templateIndex;
    float baseFreq;
    mat3 forwardRot, backwardRot;


    // data basis
    switch (b.axis) {
    case AXIS::X:
        baseLvl = b.freqLvl.y; // or z
        log2Aniso = b.freqLvl.x - baseLvl;
        forwardRot = mat3(0, 0, -1, 0, 1, 0, 1, 0, 0);
        backwardRot = mat3(0, 0, 1, 0, 1, 0, -1, 0, 0);
        break;
    case AXIS::Y:
        baseLvl = b.freqLvl.x; // or z
        log2Aniso = b.freqLvl.y - baseLvl;
        forwardRot = mat3(1, 0, 0, 0, 0, -1, 0, 1, 0);
        backwardRot = mat3(1, 0, 0, 0, 0, 1, 0, -1, 0);
        break;
    case AXIS::Z:
        baseLvl = b.freqLvl.x; // or y
        log2Aniso = b.freqLvl.z - baseLvl;
        forwardRot = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
        backwardRot = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
        break;
    }
    baseFreq = powf(2.f, baseLvl);

    if (log2Aniso == 0) { templateIndex = 0; }
    else if (log2Aniso < 0) { templateIndex = -2 * log2Aniso - 1; }
    else { templateIndex = 2 * log2Aniso; }

    auto connBasis = createPullConnection(
        &basisFlowTemplates[templateIndex]->out_vectorsMetadataTexture3D,
        &pipelineIntegrateBasisGrid_onePerDispatch->in_texBasis);
    connBasis->activate();
    pipelineIntegrateBasisGrid_onePerDispatch->in_postScalingBasis.receive(baseFreq);
    pipelineIntegrateBasisGrid_onePerDispatch->in_centerBasis.receive(b.center);
    pipelineIntegrateBasisGrid_onePerDispatch->in_halfWidthsBasis.receive(0.5f*vec3(supportBasis.right - supportBasis.left, supportBasis.top - supportBasis.bottom, supportBasis.front - supportBasis.back));
    pipelineIntegrateBasisGrid_onePerDispatch->in_forwardRotBasis.receive(forwardRot);
    pipelineIntegrateBasisGrid_onePerDispatch->in_backwardRotBasis.receive(backwardRot);


    //data vector field
    auto connField = createPullConnection(
        &velField->out_vectorsMetadataTexture3D,
        &pipelineIntegrateBasisGrid_onePerDispatch->in_texField);
    connField->activate();
    pipelineIntegrateBasisGrid_onePerDispatch->in_centerField.receive(domainCenter);
    pipelineIntegrateBasisGrid_onePerDispatch->in_halfWidthsField.receive(0.5f*vec3(supportField.right - supportField.left, supportField.top - supportField.bottom, supportField.front - supportField.back));


    // other data
    auto connTransfer = createPushConnection(
        &pipelineIntegrateBasisGrid_onePerDispatch->out_transferBuffer,
        &integrationTransferBufferGpu->in_metadataBuffer);
    connTransfer->activate();
    pipelineIntegrateBasisGrid_onePerDispatch->in_integralGridRes.receive(integralGridRes);
    pipelineIntegrateBasisGrid_onePerDispatch->in_supportInterLeftBottom.receive(vec3(supLeft, supBottom, supBack));
    pipelineIntegrateBasisGrid_onePerDispatch->in_supportInterRightTop.receive(vec3(supRight, supTop, supFront));


    unsigned int nbGroupDiv = 1;
    unsigned int it = 0;
    uvec3 nbGroups;
    do
    {
        nbGroups = uvec3(((integralGridRes + 1) - 1) / nbGroupDiv / INTEGRAL_GPU_GROUP_DIM + 1,
            ((integralGridRes + 1) - 1) / nbGroupDiv / INTEGRAL_GPU_GROUP_DIM + 1,
            ((integralGridRes + 1) - 1) / nbGroupDiv / INTEGRAL_GPU_GROUP_DIM + 1);

        pipelineIntegrateBasisGrid_onePerDispatch->in_globalIteration.receive(it++);
        pipelineIntegrateBasisGrid_onePerDispatch->goglu_nbWorkGroups.set(glm::uvec3(nbGroups.x, nbGroups.y, nbGroups.z));
        pipelineIntegrateBasisGrid_onePerDispatch->execute();
        glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

        nbGroupDiv *= INTEGRAL_GPU_GROUP_DIM;

    } while (nbGroups.x > 1 && nbGroups.y > 1 && nbGroups.z > 1);


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

    for (uint i = 0; i <= integralGridRes; i++) {
        for (uint j = 0; j <= integralGridRes; j++) {
            for (uint k = 0; k <= integralGridRes; k++) {

                vec3 p = vec3(supLeft + float(i) / integralGridRes * (supRight - supLeft),
                    supBottom + float(j) / integralGridRes * (supTop - supBottom),
                    supBack + float(k) / integralGridRes * (supFront - supBack));

                sum += ((i == 0 || i == integralGridRes) ? 0.5f : 1.f) *
                    ((j == 0 || j == integralGridRes) ? 0.5f : 1.f) *
                    ((k == 0 || k == integralGridRes) ? 0.5f : 1.f) *
                    glm::dot(translatedBasisEval(p, b.freqLvl, b.center, b.axis),
                        velField->interp(p));

            }
        }
    }

#endif

    return float(sum) * (supRight - supLeft)*(supTop - supBottom)*(supFront - supBack) / cube(integralGridRes);

#else

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

    for (uint i = 0; i <= integralGridRes; i++) {
        for (uint j = 0; j <= integralGridRes; j++) {

            vec2 p = vec2(supLeft + float(i) / integralGridRes * (supRight - supLeft),
                supBottom + float(j) / integralGridRes * (supTop - supBottom));

            sum += ((i == 0 || i == integralGridRes) ? 0.5f : 1.f) * ((j == 0 || j == integralGridRes) ? 0.5f : 1.f) *
                glm::dot(translatedBasisEval(p, b.freqLvl, b.center),
                    velField->interp(p));

        }
    }

#endif

    return float(sum) * (supRight - supLeft)*(supTop - supBottom) / sqr(integralGridRes);
#endif
}





#if !EXPLICIT_ENERGY_TRANSFER
// integrates bk . ((grad bi) . bj)  (basis I deformed by basis J, projected onto basis K (no inversion))
float AppBasisFluid::integrateBasisGradBasis(BasisFlow bi, BasisFlow bj, BasisFlow bk) {
    // compute intersection of supports
    BasisSupport supportI = bi.getSupport();
    BasisSupport supportJ = bj.getSupport();
    BasisSupport supportK = bk.getSupport();
    float supLeft = glm::max(supportI.left, glm::max(supportJ.left, supportK.left));
    float supRight = glm::min(supportI.right, glm::min(supportJ.right, supportK.right));
    float supBottom = glm::max(supportI.bottom, glm::max(supportJ.bottom, supportK.bottom));
    float supTop = glm::min(supportI.top, glm::min(supportJ.top, supportK.top));

    if (supLeft >= supRight || supBottom >= supTop
        ) {
        return 0.f;
    }

    // compute rigid motion vector
    vec2 rigidVec = matTCoeff(bi, bj);

    // compute integral as discretized sum at grid centers
#if INTEGRATION_SUM_DOUBLE_PRECISION
    double sum = 0;
#else
    float sum = 0;
#endif

#if DEF_COEFF_COMPUTE_GPU

    int minLvl;

    // data K
    auto connK = createPullConnection(
        &basisFlowTemplates[abs(bk.freqLvl.x - bk.freqLvl.y)]->out_vectorsMetadataTexture2D,
        &pipelineIntegrateBasisGradBasis->in_texBk);
    connK->activate();
    minLvl = glm::min<int>(bk.freqLvl.x, bk.freqLvl.y);
    pipelineIntegrateBasisGradBasis->in_postScalingK.receive(float(1 << minLvl));
    pipelineIntegrateBasisGradBasis->in_centerK.receive(bk.center);
    pipelineIntegrateBasisGradBasis->in_halfWidthsK.receive(0.5f*vec2(supportK.right - supportK.left, supportK.top - supportK.bottom));
    if (bk.freqLvl.x <= bk.freqLvl.y) {
        pipelineIntegrateBasisGradBasis->in_forwardRotK.receive(mat2(1, 0, 0, 1));
        pipelineIntegrateBasisGradBasis->in_backwardRotK.receive(mat2(1, 0, 0, 1));
    }
    else {
        pipelineIntegrateBasisGradBasis->in_forwardRotK.receive(mat2(0, 1, -1, 0));
        pipelineIntegrateBasisGradBasis->in_backwardRotK.receive(mat2(0, -1, 1, 0));
    }

    // data I
    auto connI = createPullConnection(
        &basisFlowTemplates[abs(bi.freqLvl.x - bi.freqLvl.y)]->out_vectorsMetadataTexture2D,
        &pipelineIntegrateBasisGradBasis->in_texBi);
    connI->activate();
    minLvl = glm::min<int>(bi.freqLvl.x, bi.freqLvl.y);
    pipelineIntegrateBasisGradBasis->in_postScalingI.receive(float(1 << minLvl));
    pipelineIntegrateBasisGradBasis->in_centerI.receive(bi.center);
    pipelineIntegrateBasisGradBasis->in_halfWidthsI.receive(0.5f*vec2(supportI.right - supportI.left, supportI.top - supportI.bottom));
    if (bi.freqLvl.x <= bi.freqLvl.y) {
        pipelineIntegrateBasisGradBasis->in_forwardRotI.receive(mat2(1, 0, 0, 1));
        pipelineIntegrateBasisGradBasis->in_backwardRotI.receive(mat2(1, 0, 0, 1));
    }
    else {
        pipelineIntegrateBasisGradBasis->in_forwardRotI.receive(mat2(0, 1, -1, 0));
        pipelineIntegrateBasisGradBasis->in_backwardRotI.receive(mat2(0, -1, 1, 0));
    }

    // data J
    auto connJ = createPullConnection(
        &basisFlowTemplates[abs(bj.freqLvl.x - bj.freqLvl.y)]->out_vectorsMetadataTexture2D,
        &pipelineIntegrateBasisGradBasis->in_texBj);
    connJ->activate();
    minLvl = glm::min<int>(bj.freqLvl.x, bj.freqLvl.y);
    pipelineIntegrateBasisGradBasis->in_postScalingJ.receive(float(1 << minLvl));
    pipelineIntegrateBasisGradBasis->in_centerJ.receive(bj.center);
    pipelineIntegrateBasisGradBasis->in_halfWidthsJ.receive(0.5f*vec2(supportJ.right - supportJ.left, supportJ.top - supportJ.bottom));
    if (bj.freqLvl.x <= bj.freqLvl.y) {
        pipelineIntegrateBasisGradBasis->in_forwardRotJ.receive(mat2(1, 0, 0, 1));
        pipelineIntegrateBasisGradBasis->in_backwardRotJ.receive(mat2(1, 0, 0, 1));
    }
    else {
        pipelineIntegrateBasisGradBasis->in_forwardRotJ.receive(mat2(0, 1, -1, 0));
        pipelineIntegrateBasisGradBasis->in_backwardRotJ.receive(mat2(0, -1, 1, 0));
    }

    // other data
    auto connTransfer = createPushConnection(
        &pipelineIntegrateBasisGradBasis->out_transferBuffer,
        &integrationTransferBufferGpu->in_metadataBuffer);
    connTransfer->activate();
    pipelineIntegrateBasisGradBasis->in_integralGridRes.receive(integralGridRes);
    pipelineIntegrateBasisGradBasis->in_supportInterLeftBottom.receive(vec2(supLeft, supBottom));
    pipelineIntegrateBasisGradBasis->in_supportInterRightTop.receive(vec2(supRight, supTop));
    pipelineIntegrateBasisGradBasis->in_epsilonGrad.receive(1.0f / (float(nbCells) * vec2(1 << bi.freqLvl) * 2.0f));
    pipelineIntegrateBasisGradBasis->in_rigidVec.receive(vec2(rigidVec));


    unsigned int nbGroupDiv = 1;
    unsigned int it = 0;
    uvec2 nbGroups;
    do
    {
        nbGroups = uvec2(((integralGridRes + 1) - 1) / nbGroupDiv / INTEGRAL_GPU_GROUP_DIM + 1,
            ((integralGridRes + 1) - 1) / nbGroupDiv / INTEGRAL_GPU_GROUP_DIM + 1);

        pipelineIntegrateBasisGradBasis->in_globalIteration.receive(it++);
        pipelineIntegrateBasisGradBasis->goglu_nbWorkGroups.set(glm::uvec3(nbGroups.x, nbGroups.y, 1));
        pipelineIntegrateBasisGradBasis->execute();
        glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

        nbGroupDiv *= INTEGRAL_GPU_GROUP_DIM;

    } while (nbGroups.x > 1 && nbGroups.y > 1);


    sum = integrationTransferBufferGpu->getCpuData(0).x;

    connK->deactivate(); delete connK;
    connI->deactivate(); delete connI;
    connJ->deactivate(); delete connJ;
    connTransfer->deactivate(); delete connTransfer;

#else

    for (uint i = 0; i <= integralGridRes; i++) {
        for (uint j = 0; j <= integralGridRes; j++) {

            vec2 p = vec2(supLeft + float(i) / integralGridRes * (supRight - supLeft),
                supBottom + float(j) / integralGridRes * (supTop - supBottom));
            vec2 valK = translatedBasisEval(p, bk.freqLvl, bk.center);
            mat2 valI = translatedBasisGradEval(p, bi.freqLvl, bi.center);
            vec2 valJ = translatedBasisEval(p, bj.freqLvl, bj.center);
            sum += ((i == 0 || i == integralGridRes) ? 0.5f : 1.f) * ((j == 0 || j == integralGridRes) ? 0.5f : 1.f) *
                glm::dot(valK, valI*(valJ - rigidVec));
        }
    }

#endif

    return float(sum) * (supRight - supLeft)*(supTop - supBottom) / sqr(integralGridRes);

#endif

}

#endif





// integrates basis(...) dot basis(...)
float IntegrateBasisBasis(BasisFlow b1, BasisFlow b2) {

    BasisSupport sup1 = b1.getSupport();
    BasisSupport sup2 = b2.getSupport();
    float supLeft = glm::max(sup1.left, sup2.left);
    float supRight = glm::min(sup1.right, sup2.right);
    float supBottom = glm::max(sup1.bottom, sup2.bottom);
    float supTop = glm::min(sup1.top, sup2.top);
#if SIM3D
    float supBack = glm::max(sup1.back, sup2.back);
    float supFront = glm::min(sup1.front, sup2.front);
#endif

    if (supLeft >= supRight || supBottom >= supTop
#if SIM3D
        || supBack >= supFront
#endif
        ) {
        return 0.f;
    }


    // compute integral as discretized sum at grid centers
#if INTEGRATION_SUM_DOUBLE_PRECISION
    double sum = 0;
#else
    float sum = 0;
#endif

#if SIM3D

#if DEF_COEFF_COMPUTE_GPU

    int baseLvl, log2Aniso, templateIndex;
    float baseFreq;
    mat3 forwardRot, backwardRot;


    // data I
    switch (b1.axis) {
    case AXIS::X:
        baseLvl = b1.freqLvl.y; // or z
        log2Aniso = b1.freqLvl.x - baseLvl;
        forwardRot = mat3(0, 0, -1, 0, 1, 0, 1, 0, 0);
        backwardRot = mat3(0, 0, 1, 0, 1, 0, -1, 0, 0);
        break;
    case AXIS::Y:
        baseLvl = b1.freqLvl.x; // or z
        log2Aniso = b1.freqLvl.y - baseLvl;
        forwardRot = mat3(1, 0, 0, 0, 0, -1, 0, 1, 0);
        backwardRot = mat3(1, 0, 0, 0, 0, 1, 0, -1, 0);
        break;
    case AXIS::Z:
        baseLvl = b1.freqLvl.x; // or y
        log2Aniso = b1.freqLvl.z - baseLvl;
        forwardRot = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
        backwardRot = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
        break;
    }
    baseFreq = powf(2.f, float(baseLvl));

    if (log2Aniso == 0) { templateIndex = 0; }
    else if (log2Aniso < 0) { templateIndex = -2 * log2Aniso - 1; }
    else { templateIndex = 2 * log2Aniso; }

    auto connI = createPullConnection(
        &basisFlowTemplates[templateIndex]->out_vectorsMetadataTexture3D,
        &pipelineIntegrateBasisBasis->in_texBi);
    connI->activate();
    pipelineIntegrateBasisBasis->in_postScalingI.receive(baseFreq);
    pipelineIntegrateBasisBasis->in_centerI.receive(b1.center);
    pipelineIntegrateBasisBasis->in_halfWidthsI.receive(0.5f*vec3(sup1.right - sup1.left, sup1.top - sup1.bottom, sup1.front - sup1.back));
    pipelineIntegrateBasisBasis->in_forwardRotI.receive(forwardRot);
    pipelineIntegrateBasisBasis->in_backwardRotI.receive(backwardRot);


    //data J
    switch (b2.axis) {
    case AXIS::X:
        baseLvl = b2.freqLvl.y; // or z
        log2Aniso = b2.freqLvl.x - baseLvl;
        forwardRot = mat3(0, 0, -1, 0, 1, 0, 1, 0, 0);
        backwardRot = mat3(0, 0, 1, 0, 1, 0, -1, 0, 0);
        break;
    case AXIS::Y:
        baseLvl = b2.freqLvl.x; // or z
        log2Aniso = b2.freqLvl.y - baseLvl;
        forwardRot = mat3(1, 0, 0, 0, 0, -1, 0, 1, 0);
        backwardRot = mat3(1, 0, 0, 0, 0, 1, 0, -1, 0);
        break;
    case AXIS::Z:
        baseLvl = b2.freqLvl.x; // or y
        log2Aniso = b2.freqLvl.z - baseLvl;
        forwardRot = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
        backwardRot = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
        break;
    }
    baseFreq = powf(2.f, float(baseLvl));

    if (log2Aniso == 0) { templateIndex = 0; }
    else if (log2Aniso < 0) { templateIndex = -2 * log2Aniso - 1; }
    else { templateIndex = 2 * log2Aniso; }

    auto connJ = createPullConnection(
        &basisFlowTemplates[templateIndex]->out_vectorsMetadataTexture3D,
        &pipelineIntegrateBasisBasis->in_texBj);
    connJ->activate();
    pipelineIntegrateBasisBasis->in_postScalingJ.receive(baseFreq);
    pipelineIntegrateBasisBasis->in_centerJ.receive(b2.center);
    pipelineIntegrateBasisBasis->in_halfWidthsJ.receive(0.5f*vec3(sup2.right - sup2.left, sup2.top - sup2.bottom, sup2.front - sup2.back));
    pipelineIntegrateBasisBasis->in_forwardRotJ.receive(forwardRot);
    pipelineIntegrateBasisBasis->in_backwardRotJ.receive(backwardRot);


    // other data
    auto connTransfer = createPushConnection(
        &pipelineIntegrateBasisBasis->out_transferBuffer,
        &integrationTransferBufferGpu->in_metadataBuffer);
    connTransfer->activate();
    pipelineIntegrateBasisBasis->in_integralGridRes.receive(integralGridRes);
    pipelineIntegrateBasisBasis->in_supportInterLeftBottom.receive(vec3(supLeft, supBottom, supBack));
    pipelineIntegrateBasisBasis->in_supportInterRightTop.receive(vec3(supRight, supTop, supFront));


    unsigned int nbGroupDiv = 1;
    unsigned int it = 0;
    uvec3 nbGroups;
    do
    {
        nbGroups = uvec3(((integralGridRes + 1) - 1) / nbGroupDiv / INTEGRAL_GPU_GROUP_DIM + 1,
            ((integralGridRes + 1) - 1) / nbGroupDiv / INTEGRAL_GPU_GROUP_DIM + 1,
            ((integralGridRes + 1) - 1) / nbGroupDiv / INTEGRAL_GPU_GROUP_DIM + 1);

        pipelineIntegrateBasisBasis->in_globalIteration.receive(it++);
        pipelineIntegrateBasisBasis->goglu_nbWorkGroups.set(glm::uvec3(nbGroups.x, nbGroups.y, nbGroups.z));
        pipelineIntegrateBasisBasis->execute();
        glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

        nbGroupDiv *= INTEGRAL_GPU_GROUP_DIM;

    } while (nbGroups.x > 1 && nbGroups.y > 1 && nbGroups.z > 1);


    sum = integrationTransferBufferGpu->getCpuData(0).x;

    connI->deactivate(); delete connI;
    connJ->deactivate(); delete connJ;
    connTransfer->deactivate(); delete connTransfer;

    return float(sum) * (supRight - supLeft)*(supTop - supBottom)*(supFront - supBack) / cube(integralGridRes);

#else
#if INTEGRATE_BASIS_BASIS_GRID_ALIGNED

    uint resX = ceil(INTEGRATE_BASIS_BASIS_GRID_ALIGNED_RESOLUTION*(domainRight - domainLeft));
    uint resY = ceil(INTEGRATE_BASIS_BASIS_GRID_ALIGNED_RESOLUTION*(domainTop - domainBottom));
    uint resZ = ceil(INTEGRATE_BASIS_BASIS_GRID_ALIGNED_RESOLUTION*(domainFront - domainBack));

    uint iMin = clamp<int>(floor((supLeft - domainLeft) / (domainRight - domainLeft)*resX), 0, resX - 1);
    uint iMax = clamp<int>(ceil((supRight - domainLeft) / (domainRight - domainLeft)*resX), 0, resX - 1);
    uint jMin = clamp<int>(floor((supBottom - domainBottom) / (domainTop - domainBottom)*resY), 0, resY - 1);
    uint jMax = clamp<int>(ceil((supTop - domainBottom) / (domainTop - domainBottom)*resY), 0, resY - 1);
    uint kMin = clamp<int>(floor((supBack - domainBack) / (domainFront - domainBack)*resZ), 0, resZ - 1);
    uint kMax = clamp<int>(ceil((supFront - domainBack) / (domainFront - domainBack)*resZ), 0, resZ - 1);

    for (uint i = iMin; i <= iMax; i++) {
        for (uint j = jMin; j <= jMax; j++) {
            for (uint k = kMin; k <= kMax; k++) {
                dvec3 p = dvec3(domainLeft + double(i + 0.5f) / resX * (domainRight - domainLeft),
                    domainBottom + double(j + 0.5f) / resY * (domainTop - domainBottom),
                    domainBack + double(k + 0.5f) / resZ * (domainFront - domainBack));
                sum += glm::dot(translatedBasisEval(p, b1.freqLvl, dvec3(b1.center), b1.axis),
                    translatedBasisEval(p, b2.freqLvl, dvec3(b2.center), b2.axis));
            }
        }
    }

    return float(sum) * (domainRight - domainLeft)*(domainTop - domainBottom)*(domainFront - domainBack)
        / (resX) / (resY) / (resZ);

#else
    for (uint i = 0; i <= integralGridRes; i++) {
        for (uint j = 0; j <= integralGridRes; j++) {
            for (uint k = 0; k <= integralGridRes; k++) {

                vec3 p = vec3(supLeft + float(i) / integralGridRes * (supRight - supLeft),
                    supBottom + float(j) / integralGridRes * (supTop - supBottom),
                    supBack + float(k) / integralGridRes * (supFront - supBack));
                sum += ((i == 0 || i == integralGridRes) ? 0.5f : 1.f) *
                    ((j == 0 || j == integralGridRes) ? 0.5f : 1.f) *
                    ((k == 0 || k == integralGridRes) ? 0.5f : 1.f) *
                    glm::dot(translatedBasisEval(vec3(p), b1.freqLvl, b1.center, b1.axis),
                        translatedBasisEval(vec3(p), b2.freqLvl, b2.center, b2.axis));
            }
        }
    }

    return float(sum) * (supRight - supLeft)*(supTop - supBottom)*(supFront - supBack) / cube(integralGridRes);

#endif

#endif



#else

#if DEF_COEFF_COMPUTE_GPU

    int minLvl;


    // data I
    auto connI = createPullConnection(
        &basisFlowTemplates[abs(b1.freqLvl.x - b1.freqLvl.y)]->out_vectorsMetadataTexture2D,
        &pipelineIntegrateBasisBasis->in_texBi);
    connI->activate();
    minLvl = glm::min<int>(b1.freqLvl.x, b1.freqLvl.y);
    pipelineIntegrateBasisBasis->in_postScalingI.receive(float(1 << minLvl));
    pipelineIntegrateBasisBasis->in_centerI.receive(b1.center);
    pipelineIntegrateBasisBasis->in_halfWidthsI.receive(0.5f*vec2(sup1.right - sup1.left, sup1.top - sup1.bottom));
    if (b1.freqLvl.x <= b1.freqLvl.y) {
        pipelineIntegrateBasisBasis->in_forwardRotI.receive(mat2(1, 0, 0, 1));
        pipelineIntegrateBasisBasis->in_backwardRotI.receive(mat2(1, 0, 0, 1));
    }
    else {
        pipelineIntegrateBasisBasis->in_forwardRotI.receive(mat2(0, 1, -1, 0));
        pipelineIntegrateBasisBasis->in_backwardRotI.receive(mat2(0, -1, 1, 0));
    }

    // data J
    auto connJ = createPullConnection(
        &basisFlowTemplates[abs(b2.freqLvl.x - b2.freqLvl.y)]->out_vectorsMetadataTexture2D,
        &pipelineIntegrateBasisBasis->in_texBj);
    connJ->activate();
    minLvl = glm::min<int>(b2.freqLvl.x, b2.freqLvl.y);
    pipelineIntegrateBasisBasis->in_postScalingJ.receive(float(1 << minLvl));
    pipelineIntegrateBasisBasis->in_centerJ.receive(b2.center);
    pipelineIntegrateBasisBasis->in_halfWidthsJ.receive(0.5f*vec2(sup2.right - sup2.left, sup2.top - sup2.bottom));
    if (b2.freqLvl.x <= b2.freqLvl.y) {
        pipelineIntegrateBasisBasis->in_forwardRotJ.receive(mat2(1, 0, 0, 1));
        pipelineIntegrateBasisBasis->in_backwardRotJ.receive(mat2(1, 0, 0, 1));
    }
    else {
        pipelineIntegrateBasisBasis->in_forwardRotJ.receive(mat2(0, 1, -1, 0));
        pipelineIntegrateBasisBasis->in_backwardRotJ.receive(mat2(0, -1, 1, 0));
    }

    // other data
    auto connTransfer = createPushConnection(
        &pipelineIntegrateBasisBasis->out_transferBuffer,
        &integrationTransferBufferGpu->in_metadataBuffer);
    connTransfer->activate();
    pipelineIntegrateBasisBasis->in_integralGridRes.receive(integralGridRes);
    pipelineIntegrateBasisBasis->in_supportInterLeftBottom.receive(vec2(supLeft, supBottom));
    pipelineIntegrateBasisBasis->in_supportInterRightTop.receive(vec2(supRight, supTop));


    unsigned int nbGroupDiv = 1;
    unsigned int it = 0;
    uvec2 nbGroups;
    do
    {
        nbGroups = uvec2(((integralGridRes + 1) - 1) / nbGroupDiv / INTEGRAL_GPU_GROUP_DIM + 1,
            ((integralGridRes + 1) - 1) / nbGroupDiv / INTEGRAL_GPU_GROUP_DIM + 1);

        pipelineIntegrateBasisBasis->in_globalIteration.receive(it++);
        pipelineIntegrateBasisBasis->goglu_nbWorkGroups.set(glm::uvec3(nbGroups.x, nbGroups.y, 1));
        pipelineIntegrateBasisBasis->execute();
        glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

        nbGroupDiv *= INTEGRAL_GPU_GROUP_DIM;

    } while (nbGroups.x > 1 && nbGroups.y > 1);


    sum = integrationTransferBufferGpu->getCpuData(0).x;

    connI->deactivate(); delete connI;
    connJ->deactivate(); delete connJ;
    connTransfer->deactivate(); delete connTransfer;

#else

    for (int i = 0; i <= integralGridRes; i++) {
        for (int j = 0; j <= integralGridRes; j++) {

            vec2 p = vec2(supLeft + float(i) / integralGridRes * (supRight - supLeft),
                supBottom + float(j) / integralGridRes * (supTop - supBottom));
            sum += ((i == 0 || i == integralGridRes) ? 0.5f : 1.f) * ((j == 0 || j == integralGridRes) ? 0.5f : 1.f) *
                glm::dot(translatedBasisEval(p, b1.freqLvl, b1.center),
                    translatedBasisEval(p, b2.freqLvl, b2.center));
        }
    }

#endif

    return float(sum) * (supRight - supLeft)*(supTop - supBottom) / sqr(integralGridRes);

#endif

}




#if EXPLICIT_TRANSPORT_ROTATION

// gives the average value of bVec over the support of bSupport
vec2 AppBasisFluid::averageBasisOnSupport(BasisFlow bVec, BasisFlow bSupport) {

    // compute intersection of supports
    BasisSupport supportVec = bVec.getSupport();
    BasisSupport supportSup = bSupport.getSupport();
    float supLeft = glm::max(supportVec.left, supportSup.left);
    float supRight = glm::min(supportVec.right, supportSup.right);
    float supBottom = glm::max(supportVec.bottom, supportSup.bottom);
    float supTop = glm::min(supportVec.top, supportSup.top);
#if SIM3D
    float supBack = glm::max(supportVec.back, supportSup.back);
    float supFront = glm::min(supportVec.front, supportSup.front);
#endif

    if (supLeft >= supRight || supBottom >= supTop
#if SIM3D
        || supBack >= supFront
#endif
        ) {
        return vec2(0);
    }

    // compute integral as discretized sum at grid centers
#if INTEGRATION_SUM_DOUBLE_PRECISION
    dvec2 sum(0);
#else
    vec2 sum(0);
#endif

#if SIM3D

#if DEF_COEFF_COMPUTE_GPU

    int baseLvl, log2Aniso, templateIndex;
    float baseFreq;
    mat3 forwardRot, backwardRot;


    // data basis i
    switch (bVec.axis) {
    case AXIS::X:
        baseLvl = bVec.freqLvl.y; // or z
        log2Aniso = bVec.freqLvl.x - baseLvl;
        forwardRot = mat3(0, 0, -1, 0, 1, 0, 1, 0, 0);
        backwardRot = mat3(0, 0, 1, 0, 1, 0, -1, 0, 0);
        break;
    case AXIS::Y:
        baseLvl = bVec.freqLvl.x; // or z
        log2Aniso = bVec.freqLvl.y - baseLvl;
        forwardRot = mat3(1, 0, 0, 0, 0, -1, 0, 1, 0);
        backwardRot = mat3(1, 0, 0, 0, 0, 1, 0, -1, 0);
        break;
    case AXIS::Z:
        baseLvl = bVec.freqLvl.x; // or y
        log2Aniso = bVec.freqLvl.z - baseLvl;
        forwardRot = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
        backwardRot = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
        break;
    }
    baseFreq = powf(2.f, float(baseLvl));

    if (log2Aniso == 0) { templateIndex = 0; }
    else if (log2Aniso < 0) { templateIndex = -2 * log2Aniso - 1; }
    else { templateIndex = 2 * log2Aniso; }

    auto connI = createPullConnection(
        &basisFlowTemplates[templateIndex]->out_vectorsMetadataTexture3D,
        &pipelineIntegrateAvgBasis->in_texBi);
    connI->activate();
    pipelineIntegrateAvgBasis->in_postScalingI.receive(baseFreq); // TODO: could be done on CPU once to save one multiplications on the GPU.
    pipelineIntegrateAvgBasis->in_centerI.receive(bVec.center);
    pipelineIntegrateAvgBasis->in_halfWidthsI.receive(0.5f*vec3(supportVec.right - supportVec.left, supportVec.top - supportVec.bottom, supportVec.front - supportVec.back));
    pipelineIntegrateAvgBasis->in_forwardRotI.receive(forwardRot);
    pipelineIntegrateAvgBasis->in_backwardRotI.receive(backwardRot);

    // other data
    auto connTransfer = createPushConnection(
        &pipelineIntegrateAvgBasis->out_transferBuffer,
        &integrationTransferBufferGpu->in_metadataBuffer);
    connTransfer->activate();
    pipelineIntegrateAvgBasis->in_integralGridRes.receive(integralGridRes);
    pipelineIntegrateAvgBasis->in_supportInterLeftBottom.receive(vec3(supLeft, supBottom, supBack));
    pipelineIntegrateAvgBasis->in_supportInterRightTop.receive(vec3(supRight, supTop, supFront));


    unsigned int nbGroupDiv = 1;
    unsigned int it = 0;
    uvec3 nbGroups;
    do
    {
        nbGroups = uvec3(((integralGridRes + 1) - 1) / nbGroupDiv / INTEGRAL_GPU_GROUP_DIM + 1,
            ((integralGridRes + 1) - 1) / nbGroupDiv / INTEGRAL_GPU_GROUP_DIM + 1,
            ((integralGridRes + 1) - 1) / nbGroupDiv / INTEGRAL_GPU_GROUP_DIM + 1);

        pipelineIntegrateAvgBasis->in_globalIteration.receive(it++);
        pipelineIntegrateAvgBasis->goglu_nbWorkGroups.set(glm::uvec3(nbGroups.x, nbGroups.y, nbGroups.z));
        pipelineIntegrateAvgBasis->execute();
        glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

        nbGroupDiv *= INTEGRAL_GPU_GROUP_DIM;

    } while (nbGroups.x > 1 && nbGroups.y > 1 && nbGroups.z > 1);


    sum = vec2(integrationTransferBufferGpu->getCpuData(0));  // vec4 -> vec3

    connI->deactivate(); delete connI;
    connTransfer->deactivate(); delete connTransfer;

#else

    for (uint i = 0; i <= integralGridRes; i++) {
        for (uint j = 0; j <= integralGridRes; j++) {
            for (uint k = 0; k <= integralGridRes; k++) {

                vec3 p = vec3(supLeft + float(i) / integralGridRes * (supRight - supLeft),
                    supBottom + float(j) / integralGridRes * (supTop - supBottom),
                    supBack + float(k) / integralGridRes * (supFront - supBack));
                sum += ((i == 0 || i == integralGridRes) ? 0.5f : 1.f) *
                    ((j == 0 || j == integralGridRes) ? 0.5f : 1.f) *
                    ((k == 0 || k == integralGridRes) ? 0.5f : 1.f) *
                    translatedBasisEval(p, bVec.freqLvl, bVec.center, bVec.axis);

            }
        }
    }

#endif

    return vec3(sum) * (supRight - supLeft)*(supTop - supBottom)*(supFront - supBack) / float(cube(integralGridRes)) / ((1.f / float(1 << bSupport.freqLvl.x))*(1.f / float(1 << bSupport.freqLvl.y))*(1.f / float(1 << bSupport.freqLvl.z)));


#else

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

    for (uint i = 0; i <= integralGridRes; i++) {
        for (uint j = 0; j <= integralGridRes; j++) {

            vec2 p = vec2(supLeft + float(i) / integralGridRes * (supRight - supLeft),
                supBottom + float(j) / integralGridRes * (supTop - supBottom));




            sum += ((i == 0 || i == integralGridRes) ? 0.5f : 1.f) * ((j == 0 || j == integralGridRes) ? 0.5f : 1.f) *
                translatedBasisEval(p, bVec.freqLvl, bVec.center);
        }
    }

#endif

    // divide by domain size to get averaged value
    return vec2(sum) * (supRight - supLeft)*(supTop - supBottom) / float(sqr(integralGridRes)) / ((1.f / float(1 << bSupport.freqLvl.x))*(1.f / float(1 << bSupport.freqLvl.y)));

#endif
}


#if SIM3D

// gives the average rotation bSupport support around its center, as caused by bVec
vec2 AppBasisFluid::averageBasisRotationOnSupport(BasisFlow bVec, BasisFlow bSupport) {

    // compute intersection of supports
    BasisSupport supportVec = bVec.getSupport();
    BasisSupport supportSup = bSupport.getSupport();
    float supLeft = glm::max(supportVec.left, supportSup.left);
    float supRight = glm::min(supportVec.right, supportSup.right);
    float supBottom = glm::max(supportVec.bottom, supportSup.bottom);
    float supTop = glm::min(supportVec.top, supportSup.top);
    float supBack = glm::max(supportVec.back, supportSup.back);
    float supFront = glm::min(supportVec.front, supportSup.front);

    if (supLeft >= supRight || supBottom >= supTop || supBack >= supFront
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

    int baseLvl, log2Aniso, templateIndex;
    float baseFreq;
    mat3 forwardRot, backwardRot;


    // data basis i
    switch (bVec.axis) {
    case AXIS::X:
        baseLvl = bVec.freqLvl.y; // or z
        log2Aniso = bVec.freqLvl.x - baseLvl;
        forwardRot = mat3(0, 0, -1, 0, 1, 0, 1, 0, 0);
        backwardRot = mat3(0, 0, 1, 0, 1, 0, -1, 0, 0);
        break;
    case AXIS::Y:
        baseLvl = bVec.freqLvl.x; // or z
        log2Aniso = bVec.freqLvl.y - baseLvl;
        forwardRot = mat3(1, 0, 0, 0, 0, -1, 0, 1, 0);
        backwardRot = mat3(1, 0, 0, 0, 0, 1, 0, -1, 0);
        break;
    case AXIS::Z:
        baseLvl = bVec.freqLvl.x; // or y
        log2Aniso = bVec.freqLvl.z - baseLvl;
        forwardRot = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
        backwardRot = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
        break;
    }
    baseFreq = powf(2.f, float(baseLvl));

    if (log2Aniso == 0) { templateIndex = 0; }
    else if (log2Aniso < 0) { templateIndex = -2 * log2Aniso - 1; }
    else { templateIndex = 2 * log2Aniso; }

    auto connI = createPullConnection(
        &basisFlowTemplates[templateIndex]->out_vectorsMetadataTexture3D,
        &pipelineIntegrateAvgBasisRotation->in_texBi);
    connI->activate();
    pipelineIntegrateAvgBasisRotation->in_postScalingI.receive(baseFreq); // TODO: could be done on CPU once to save one multiplications on the GPU.
    pipelineIntegrateAvgBasisRotation->in_centerI.receive(bVec.center);
    pipelineIntegrateAvgBasisRotation->in_centerTransported.receive(bSupport.center);
    pipelineIntegrateAvgBasisRotation->in_halfWidthsI.receive(0.5f*vec3(supportVec.right - supportVec.left, supportVec.top - supportVec.bottom, supportVec.front - supportVec.back));
    pipelineIntegrateAvgBasisRotation->in_forwardRotI.receive(forwardRot);
    pipelineIntegrateAvgBasisRotation->in_backwardRotI.receive(backwardRot);

    // other data
    auto connTransfer = createPushConnection(
        &pipelineIntegrateAvgBasisRotation->out_transferBuffer,
        &integrationTransferBufferGpu->in_metadataBuffer);
    connTransfer->activate();
    pipelineIntegrateAvgBasisRotation->in_integralGridRes.receive(integralGridRes);
    pipelineIntegrateAvgBasisRotation->in_supportInterLeftBottom.receive(vec3(supLeft, supBottom, supBack));
    pipelineIntegrateAvgBasisRotation->in_supportInterRightTop.receive(vec3(supRight, supTop, supFront));


    unsigned int nbGroupDiv = 1;
    unsigned int it = 0;
    uvec3 nbGroups;
    do
    {
        nbGroups = uvec3(((integralGridRes + 1) - 1) / nbGroupDiv / INTEGRAL_GPU_GROUP_DIM + 1,
            ((integralGridRes + 1) - 1) / nbGroupDiv / INTEGRAL_GPU_GROUP_DIM + 1,
            ((integralGridRes + 1) - 1) / nbGroupDiv / INTEGRAL_GPU_GROUP_DIM + 1);

        pipelineIntegrateAvgBasisRotation->in_globalIteration.receive(it++);
        pipelineIntegrateAvgBasisRotation->goglu_nbWorkGroups.set(glm::uvec3(nbGroups.x, nbGroups.y, nbGroups.z));
        pipelineIntegrateAvgBasisRotation->execute();
        glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

        nbGroupDiv *= INTEGRAL_GPU_GROUP_DIM;

    } while (nbGroups.x > 1 && nbGroups.y > 1 && nbGroups.z > 1);


    sum = vec2(integrationTransferBufferGpu->getCpuData(0));  // vec4 -> vec3

    connI->deactivate(); delete connI;
    connTransfer->deactivate(); delete connTransfer;

#else

    for (uint i = 0; i <= integralGridRes; i++) {
        for (uint j = 0; j <= integralGridRes; j++) {
            for (uint k = 0; k <= integralGridRes; k++) {

                vec3 p = vec3(supLeft + float(i) / integralGridRes * (supRight - supLeft),
                    supBottom + float(j) / integralGridRes * (supTop - supBottom),
                    supBack + float(k) / integralGridRes * (supFront - supBack));
                if (vecNorm(p - bSupport.center) > 1e-5) {
                    sum += ((i == 0 || i == integralGridRes) ? 0.5f : 1.f) *
                        ((j == 0 || j == integralGridRes) ? 0.5f : 1.f) *
                        ((k == 0 || k == integralGridRes) ? 0.5f : 1.f) *
                        glm::cross(
                        (p - bSupport.center) / sqr(vecNorm(p - bSupport.center)),
                            translatedBasisEval(p, bVec.freqLvl, bVec.center, bVec.axis)
                        );
                }

            }
        }
    }

#endif

    return vec3(sum) * (supRight - supLeft)*(supTop - supBottom)*(supFront - supBack) / float(cube(integralGridRes)) / ((1.f / float(1 << bSupport.freqLvl.x))*(1.f / float(1 << bSupport.freqLvl.y))*(1.f / float(1 << bSupport.freqLvl.z)));

}

#endif

#endif



void AppBasisFluid::saveCoeffsBB(string filename)
{
    if (newBBCoeffComputed) {
        ofstream file;
        file.open(filename);
        for (auto coeffPairIt = coeffsBB.begin(); coeffPairIt != coeffsBB.end(); coeffPairIt++) {
            KeyTypeBB key = coeffPairIt->first;
            float val = coeffPairIt->second;
#if SIM3D
            file << get< 0>(key) << " " <<
                get< 1>(key) << " " <<
                get< 2>(key) << " " <<
                (int)(get< 3>(key)) << " " <<
                get< 4>(key) << " " <<
                get< 5>(key) << " " <<
                get< 6>(key) << " " <<
                (int)(get< 7>(key)) << " " <<
                get< 8>(key) << " " <<
                get< 9>(key) << " " <<
                get<10>(key) << " " <<
                val << endl;
#else
            file << get<0>(key) << " " <<
                get<1>(key) << " " <<
                get<2>(key) << " " <<
                get<3>(key) << " " <<
                get<4>(key) << " " <<
                get<5>(key) << " " <<
                val << endl;
#endif
        }
        file.close();
        cout << "saved BB coefficient to " << filename << endl;
    }
}


void AppBasisFluid::loadCoeffsBB(string filename)
{
    ifstream file;
    file.open(filename);
    string line;
    stringstream ss;
    while (getline(file, line)) {
#if SIM3D
        int dataUI[8];
        float dataF[4];
        ss.clear();
        ss.str(line);
        ss >> dataUI[0] >>
            dataUI[1] >>
            dataUI[2] >>
            dataUI[3] >>
            dataUI[4] >>
            dataUI[5] >>
            dataUI[6] >>
            dataUI[7] >>
            dataF[0] >>
            dataF[1] >>
            dataF[2] >>
            dataF[3];
        KeyTypeBB key = make_tuple(
            dataUI[0],
            dataUI[1],
            dataUI[2],
            AXIS(dataUI[3]),
            dataUI[4],
            dataUI[5],
            dataUI[6],
            AXIS(dataUI[7]),
            roundToMultiple(dataF[0], coeffSnapSize),
            roundToMultiple(dataF[1], coeffSnapSize),
            roundToMultiple(dataF[2], coeffSnapSize)
        );
        float coeff = dataF[3];
#else
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
            roundToMultiple(dataF[0], coeffSnapSize),
            roundToMultiple(dataF[1], coeffSnapSize)
        );
        float coeff = dataF[2];
#endif
        coeffsBB.insert(std::pair<KeyTypeBB, float>(key, coeff));
    }
    file.close();
}



#if EXPLICIT_TRANSPORT_ROTATION

void AppBasisFluid::saveCoeffsT(string filename)
{
    if (newTCoeffComputed) {
        ofstream file;
        file.open(filename);
        for (auto coeffPairIt = coeffsT.begin(); coeffPairIt != coeffsT.end(); coeffPairIt++) {
            KeyTypeT key = coeffPairIt->first;
            vec2 val = coeffPairIt->second;
#if SIM3D
            file << get< 0>(key) << " " <<
                get< 1>(key) << " " <<
                get< 2>(key) << " " <<
                (int)(get< 3>(key)) << " " <<
                get< 4>(key) << " " <<
                get< 5>(key) << " " <<
                get< 6>(key) << " " <<
                (int)(get< 7>(key)) << " " <<
                get< 8>(key) << " " <<
                get< 9>(key) << " " <<
                get<10>(key) << " " <<
                val.x << " " <<
                val.y << " " <<
                val.z << endl;
#else
            file << get<0>(key) << " " <<
                get<1>(key) << " " <<
                get<2>(key) << " " <<
                get<3>(key) << " " <<
                get<4>(key) << " " <<
                get<5>(key) << " " <<
                val.x << " " <<
                val.y << endl;
#endif
        }
        file.close();
        cout << "saved R coefficient to " << filename << endl;
    }
}


void AppBasisFluid::loadCoeffsT(string filename)
{
    ifstream file;
    file.open(filename);
    string line;
    stringstream ss;
    while (getline(file, line)) {
#if SIM3D
        int dataUI[8];
        float dataF[6];
        ss.clear();
        ss.str(line);
        ss >> dataUI[0] >>
            dataUI[1] >>
            dataUI[2] >>
            dataUI[3] >>
            dataUI[4] >>
            dataUI[5] >>
            dataUI[6] >>
            dataUI[7] >>
            dataF[0] >>
            dataF[1] >>
            dataF[2] >>
            dataF[3] >>
            dataF[4] >>
            dataF[5];
        KeyTypeT key = make_tuple(
            dataUI[0],
            dataUI[1],
            dataUI[2],
            AXIS(dataUI[3]),
            dataUI[4],
            dataUI[5],
            dataUI[6],
            AXIS(dataUI[7]),
            roundToMultiple(dataF[0], coeffSnapSize),
            roundToMultiple(dataF[1], coeffSnapSize),
            roundToMultiple(dataF[2], coeffSnapSize)
        );
        vec3 coeff(dataF[3], dataF[4], dataF[5]);
#else
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
            roundToMultiple(dataF[0], coeffSnapSize),
            roundToMultiple(dataF[1], coeffSnapSize)
        );
        vec2 coeff(dataF[2], dataF[3]);
#endif
        coeffsT.insert(std::pair<KeyTypeT, vec2>(key, coeff));
    }
    file.close();
}


#if SIM3D

void AppBasisFluid::saveCoeffsR(string filename)
{
    if (newRCoeffComputed) {
        ofstream file;
        file.open(filename);
        for (auto coeffPairIt = coeffsR.begin(); coeffPairIt != coeffsR.end(); coeffPairIt++) {
            KeyTypeR key = coeffPairIt->first;
            vec2 val = coeffPairIt->second;
            file << get< 0>(key) << " " <<
                get< 1>(key) << " " <<
                get< 2>(key) << " " <<
                (int)(get< 3>(key)) << " " <<
                get< 4>(key) << " " <<
                get< 5>(key) << " " <<
                get< 6>(key) << " " <<
                (int)(get< 7>(key)) << " " <<
                get< 8>(key) << " " <<
                get< 9>(key) << " " <<
                get<10>(key) << " " <<
                val.x << " " <<
                val.y << " " <<
                val.z << endl;
        }
        file.close();
        cout << "saved R coefficient to " << filename << endl;
    }
}


void AppBasisFluid::loadCoeffsR(string filename)
{
    ifstream file;
    file.open(filename);
    string line;
    stringstream ss;
    while (getline(file, line)) {
        int dataUI[8];
        float dataF[6];
        ss.clear();
        ss.str(line);
        ss >> dataUI[0] >>
            dataUI[1] >>
            dataUI[2] >>
            dataUI[3] >>
            dataUI[4] >>
            dataUI[5] >>
            dataUI[6] >>
            dataUI[7] >>
            dataF[0] >>
            dataF[1] >>
            dataF[2] >>
            dataF[3] >>
            dataF[4] >>
            dataF[5];
        KeyTypeR key = make_tuple(
            dataUI[0],
            dataUI[1],
            dataUI[2],
            AXIS(dataUI[3]),
            dataUI[4],
            dataUI[5],
            dataUI[6],
            AXIS(dataUI[7]),
            roundToMultiple(dataF[0], coeffSnapSize),
            roundToMultiple(dataF[1], coeffSnapSize),
            roundToMultiple(dataF[2], coeffSnapSize)
        );
        vec3 coeff(dataF[3], dataF[4], dataF[5]);
        coeffsR.insert(std::pair<KeyTypeR, vec2>(key, coeff));
    }
    file.close();
}

#endif // #if SIM3D

#endif // #if EXPLICIT_TRANSPORT_ROTATION





#if !EXPLICIT_ENERGY_TRANSFER

void AppBasisFluid::saveCoeffsA(string filename)
{
    if (newACoeffComputed) {
        ofstream file;
        file.open(filename);
        for (auto coeffPairIt = coeffsA.begin(); coeffPairIt != coeffsA.end(); coeffPairIt++) {
            KeyTypeA key = coeffPairIt->first;
            float val = coeffPairIt->second;
#if SIM3D
            file << get< 0>(key) << " " <<
                get< 1>(key) << " " <<
                get< 2>(key) << " " <<
                (int)(get< 3>(key)) << " " <<
                get< 4>(key) << " " <<
                get< 5>(key) << " " <<
                get< 6>(key) << " " <<
                (int)(get< 7>(key)) << " " <<
                get< 8>(key) << " " <<
                get< 9>(key) << " " <<
                get<10>(key) << " " <<
                get<11>(key) << " " <<
                get<12>(key) << " " <<
                get<13>(key) << " " <<
                (int)(get<14>(key)) << " " <<
                get<15>(key) << " " <<
                get<16>(key) << " " <<
                get<17>(key) << " " <<
                val << endl;
#else
            file << get<0>(key) << " " <<
                get<1>(key) << " " <<
                get<2>(key) << " " <<
                get<3>(key) << " " <<
                get<4>(key) << " " <<
                get<5>(key) << " " <<
                get<6>(key) << " " <<
                get<7>(key) << " " <<
                get<8>(key) << " " <<
                get<9>(key) << " " <<
                val << endl;
#endif
        }
        file.close();
        cout << "saved A coefficient to " << filename << endl;
    }
}



void AppBasisFluid::loadCoeffsA(string filename)
{
    ifstream file;
    string line;
    stringstream ss;

    // count total number of lines to load
    int nbLinesTotal = 0;
    file.open(filename);
    while (getline(file, line)) {
        nbLinesTotal++;
    }
    file.close();

    file.open(filename);
    int nbLines = 0;
    while (getline(file, line)) {
#if SIM3D
        int dataUI[12];
        float dataF[7];
        ss.clear();
        ss.str(line);
        ss >> dataUI[0] >>
            dataUI[1] >>
            dataUI[2] >>
            dataUI[3] >>
            dataUI[4] >>
            dataUI[5] >>
            dataUI[6] >>
            dataUI[7] >>
            dataF[0] >>
            dataF[1] >>
            dataF[2] >>
            dataUI[8] >>
            dataUI[9] >>
            dataUI[10] >>
            dataUI[11] >>
            dataF[3] >>
            dataF[4] >>
            dataF[5] >>
            dataF[6];
        KeyTypeA key = make_tuple(
            dataUI[0],
            dataUI[1],
            dataUI[2],
            AXIS(dataUI[3]),
            dataUI[4],
            dataUI[5],
            dataUI[6],
            AXIS(dataUI[7]),
            roundToMultiple(dataF[0], coeffSnapSize),
            roundToMultiple(dataF[1], coeffSnapSize),
            roundToMultiple(dataF[2], coeffSnapSize),
            dataUI[8],
            dataUI[9],
            dataUI[10],
            AXIS(dataUI[11]),
            roundToMultiple(dataF[3], coeffSnapSize),
            roundToMultiple(dataF[4], coeffSnapSize),
            roundToMultiple(dataF[5], coeffSnapSize)
        );
        float coeff = dataF[6];
#else
        int dataUI[6];
        float dataF[5];
        ss.clear();
        ss.str(line);
        ss >> dataUI[0] >>
            dataUI[1] >>
            dataUI[2] >>
            dataUI[3] >>
            dataF[0] >>
            dataF[1] >>
            dataUI[4] >>
            dataUI[5] >>
            dataF[2] >>
            dataF[3] >>
            dataF[4];
        KeyTypeA key = make_tuple(
            dataUI[0],
            dataUI[1],
            dataUI[2],
            dataUI[3],
            roundToMultiple(dataF[0], coeffSnapSize),
            roundToMultiple(dataF[1], coeffSnapSize),
            dataUI[4],
            dataUI[5],
            roundToMultiple(dataF[2], coeffSnapSize),
            roundToMultiple(dataF[3], coeffSnapSize)
        );
        float coeff = dataF[4];
#endif
        coeffsA.insert(std::pair<KeyTypeA, float>(key, coeff));

        if (++nbLines % 10000 == 0) {
            cout << "line " << nbLines << " / " << nbLinesTotal << endl;
        }
    }
    file.close();
}



#endif



#if EXPLICIT_TRANSPORT_ROTATION

// TODO: optimize this using const&
vec2 AppBasisFluid::matTCoeff(int i, int j) {
    BasisFlow bTransported = basisFlowParams->getCpuData(i);
    BasisFlow bTransporting = basisFlowParams->getCpuData(j);
    return matTCoeff(bTransported, bTransporting);
}


// compute the integral and store it for future use
vec2 AppBasisFluid::matTCoeff(BasisFlow bTransported, BasisFlow bTransporting)
{

    if (intersectionInteriorEmpty(bTransported.getSupport(), bTransporting.getSupport())) { return vec2(0); }


    int baseLvl;
#if SIM3D
    //    int baseLvl1, baseLvl2;
    //    switch(bTransported.axis) {
    //    case AXIS::X:
    //        baseLvl1 = bTransported.freqLvl.x - bTransported.freqLvl.y; // or z
    //        break;
    //    case AXIS::Y:
    //        baseLvl1 = bTransported.freqLvl.y - bTransported.freqLvl.x; // or z
    //        break;
    //    case AXIS::Z:
    //        baseLvl1 = bTransported.freqLvl.z - bTransported.freqLvl.x; // or y
    //        break;
    //    }
    //    switch(bTransporting.axis) {
    //    case AXIS::X:
    //        baseLvl2 = bTransporting.freqLvl.x - bTransporting.freqLvl.y; // or z
    //        break;
    //    case AXIS::Y:
    //        baseLvl2 = bTransporting.freqLvl.y - bTransporting.freqLvl.x; // or z
    //        break;
    //    case AXIS::Z:
    //        baseLvl2 = bTransporting.freqLvl.z - bTransporting.freqLvl.x; // or y
    //        break;
    //    }
    //    baseLvl = glm::min(baseLvl1, baseLvl2);
    baseLvl = glm::min<int>(minVec(bTransporting.freqLvl), minVec(bTransported.freqLvl));
#else
    baseLvl = glm::min<int>(
        glm::min<int>(bTransported.freqLvl.x, bTransported.freqLvl.y),
        glm::min<int>(bTransporting.freqLvl.x, bTransporting.freqLvl.y)
        );
#endif
    float baseFreq = powf(2.f, float(baseLvl));

    // remove common frequency factors
    ivec2 normFreqLvlTransported = bTransported.freqLvl - baseLvl;
    ivec2 normFreqLvlTransporting = bTransporting.freqLvl - baseLvl;

    vec2 relativeOffset = baseFreq * vec2(
        bTransported.center.x - bTransporting.center.x,
        bTransported.center.y - bTransporting.center.y
#if SIM3D
        , bTransported.center.z - bTransporting.center.z
#endif
    );

    // snap offset to very fine grid to avoid float errors
    vec2 snappedRelativeOffset(roundToMultiple(relativeOffset.x, coeffSnapSize),
        roundToMultiple(relativeOffset.y, coeffSnapSize)
#if SIM3D
        , roundToMultiple(relativeOffset.z, coeffSnapSize)
#endif
    );

    // return coefficient if already computed, or else compute it and store it.
#if SIM3D
    KeyTypeT key = make_tuple(normFreqLvlTransported.x, normFreqLvlTransported.y, normFreqLvlTransported.z, bTransported.axis,
        normFreqLvlTransporting.x, normFreqLvlTransporting.y, normFreqLvlTransporting.z, bTransporting.axis,
        snappedRelativeOffset.x, snappedRelativeOffset.y, snappedRelativeOffset.z);
#else
    KeyTypeT key = make_tuple(normFreqLvlTransported.x, normFreqLvlTransported.y,
        normFreqLvlTransporting.x, normFreqLvlTransporting.y,
        snappedRelativeOffset.x, snappedRelativeOffset.y);
#endif
    auto coeffPair = coeffsT.find(key);
    vec2 result;

    if (coeffPair != coeffsT.end()) {
        result = coeffPair->second;
    }
    else {
        vec2 coeff;

        // compute coefficient with numerical integration
#if SIM3D
        BasisFlow bRelativeTransporting(normFreqLvlTransporting, vec2(0), bTransporting.axis);
        BasisFlow bRelativeTransported(normFreqLvlTransported, relativeOffset, bTransported.axis);
#else
        BasisFlow bRelativeTransporting(normFreqLvlTransporting, vec2(0));
        BasisFlow bRelativeTransported(normFreqLvlTransported, relativeOffset);
#endif

        coeff = averageBasisOnSupport(bRelativeTransporting, bRelativeTransported);

        coeffsT.insert(std::pair<KeyTypeT, vec2>(key, coeff));
        result = coeff;

        if (coeffsT.size() % 1000 == 0) {
            cout << "coeffs T : " << coeffsT.size() << endl;
        }

        newTCoeffComputed = true;
    }


    // scaled coefficient
#if SIM3D
    return result * baseFreq;
#else
    return result * baseFreq;
#endif

}



#endif





//TODO: optimize this using const&
float AppBasisFluid::matBBCoeff(int i, int j) {
    BasisFlow& b1 = basisFlowParams->getCpuData(i);
    BasisFlow& b2 = basisFlowParams->getCpuData(j);
    return matBBCoeff(b1, b2);
}



// compute the integral and store it for future use
float AppBasisFluid::matBBCoeff(const BasisFlow& b1, const BasisFlow& b2)
{
    if (intersectionInteriorEmpty(b1.getSupport(), b2.getSupport())) { return 0.0; }

#if ENFORCE_EXACT_ORTHOGONALITY
    if (b1.orthoGroup == b2.orthoGroup && b1.center != b2.center) {
        return 0;
    }
#endif

    int baseLvl;
#if SIM3D
    //        int baseLvl1, baseLvl2;
    //        switch(b1.axis) {
    //        case AXIS::X:
    //            baseLvl1 = b1.freqLvl.x - b1.freqLvl.y; // or z
    //            break;
    //        case AXIS::Y:
    //            baseLvl1 = b1.freqLvl.y - b1.freqLvl.x; // or z
    //            break;
    //        case AXIS::Z:
    //            baseLvl1 = b1.freqLvl.z - b1.freqLvl.x; // or y
    //            break;
    //        }
    //        switch(b2.axis) {
    //        case AXIS::X:
    //            baseLvl2 = b2.freqLvl.x - b2.freqLvl.y; // or z
    //            break;
    //        case AXIS::Y:
    //            baseLvl2 = b2.freqLvl.y - b2.freqLvl.x; // or z
    //            break;
    //        case AXIS::Z:
    //            baseLvl2 = b2.freqLvl.z - b2.freqLvl.x; // or y
    //            break;
    //        }
    //        baseLvl = glm::min(baseLvl1, baseLvl2);
    baseLvl = glm::min<int>(minVec(b1.freqLvl), minVec(b2.freqLvl));
#else
    baseLvl = glm::min<int>(
        glm::min<int>(b1.freqLvl.x, b1.freqLvl.y),
        glm::min<int>(b2.freqLvl.x, b2.freqLvl.y)
        );
#endif

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
#if SIM3D
        , b2.center.z - b1.center.z
#endif
    );



    // snap offset to very fine grid to avoid float errors
    vec2 snappedRelativeOffset = vec2(roundToMultiple(relativeOffset.x, coeffSnapSize),
        roundToMultiple(relativeOffset.y, coeffSnapSize)
#if SIM3D
        , roundToMultiple(relativeOffset.z, coeffSnapSize)
#endif
    );

    // return coefficient if already computed, or else compute it and store it.
#if SIM3D
    KeyTypeBB key = make_tuple(normFreqLvl1.x, normFreqLvl1.y, normFreqLvl1.z, b1.axis,
        normFreqLvl2.x, normFreqLvl2.y, normFreqLvl2.z, b2.axis,
        snappedRelativeOffset.x, snappedRelativeOffset.y, snappedRelativeOffset.z);
#else
    KeyTypeBB key = make_tuple(normFreqLvl1.x, normFreqLvl1.y,
        normFreqLvl2.x, normFreqLvl2.y,
        snappedRelativeOffset.x, snappedRelativeOffset.y);
#endif
    auto coeffPair = coeffsBB.find(key);
    float result;

    if (coeffPair != coeffsBB.end()) {
        result = coeffPair->second;
    }
    else {
        float coeff;

        // compute coefficient with numerical integration
#if SIM3D
        BasisFlow bRelative1(normFreqLvl1, vec2(0), b1.axis);
        BasisFlow bRelative2(normFreqLvl2, relativeOffset, b2.axis);
#else
        BasisFlow bRelative1(normFreqLvl1, vec2(0));
        BasisFlow bRelative2(normFreqLvl2, relativeOffset);
#endif
        coeff = float(integrateBasisBasis(bRelative1, bRelative2));

        coeffsBB.insert(std::pair<KeyTypeBB, float>(key, coeff));
        result = coeff;

        if (coeffsBB.size() % 1000 == 0) {
            cout << "coeffs BB : " << coeffsBB.size() << endl;
        }

        newBBCoeffComputed = true;
    }

    // scaled coefficient
#if SIM3D
    return result / baseFreq;
#else
    return result;
#endif

}




#if !EXPLICIT_ENERGY_TRANSFER

// TODO: optimize this using const&
// TODO: adjust this for it we don't use explicit transport and rotation
// TODO: add the removal of rotation when we use explicit transport and rotation but not explicit energy transfer.
float AppBasisFluid::matACoeff(int i, int j, int k) {

    BasisFlow bi = basisFlowParams->getCpuData(i);
    BasisFlow bj = basisFlowParams->getCpuData(j);
    BasisFlow bk = basisFlowParams->getCpuData(k);
    return matACoeff(bi, bj, bk);
}


// compute the integral and store it for future use
float AppBasisFluid::matACoeff(BasisFlow bi, BasisFlow bj, BasisFlow bk)
{

    //float snapSize = lengthLvl0/float(1<<maxFreqLvl)/16.0f;

    int baseLvl;
#if SIM3D
    //        int baseLvl1, baseLvl2, baseLvl3;
    //        switch(bi.axis) {
    //        case AXIS::X:
    //            baseLvl1 = bi.freqLvl.x - bi.freqLvl.y; // or z
    //            break;
    //        case AXIS::Y:
    //            baseLvl1 = bi.freqLvl.y - bi.freqLvl.x; // or z
    //            break;
    //        case AXIS::Z:
    //            baseLvl1 = bi.freqLvl.z - bi.freqLvl.x; // or y
    //            break;
    //        }
    //        switch(bj.axis) {
    //        case AXIS::X:
    //            baseLvl2 = bj.freqLvl.x - bj.freqLvl.y; // or z
    //            break;
    //        case AXIS::Y:
    //            baseLvl2 = bj.freqLvl.y - bj.freqLvl.x; // or z
    //            break;
    //        case AXIS::Z:
    //            baseLvl2 = bj.freqLvl.z - bj.freqLvl.x; // or y
    //            break;
    //        }
    //        switch(bk.axis) {
    //        case AXIS::X:
    //            baseLvl3 = bk.freqLvl.x - bk.freqLvl.y; // or z
    //            break;
    //        case AXIS::Y:
    //            baseLvl3 = bk.freqLvl.y - bk.freqLvl.x; // or z
    //            break;
    //        case AXIS::Z:
    //            baseLvl3 = bk.freqLvl.z - bk.freqLvl.x; // or y
    //            break;
    //        }
    //        baseLvl = min3(baseLvl1, baseLvl2, baseLvl3);
    baseLvl = min3<int>(minVec(bi.freqLvl), minVec(bj.freqLvl), minVec(bk.freqLvl));
#else
    baseLvl = min3<int>(
        glm::min<int>(bi.freqLvl.x, bi.freqLvl.y),
        glm::min<int>(bj.freqLvl.x, bj.freqLvl.y),
        glm::min<int>(bk.freqLvl.x, bk.freqLvl.y)
        );
#endif

    // remove common frequency factors
    ivec2 normFreqLvlI = bi.freqLvl - baseLvl;
    ivec2 normFreqLvlJ = bj.freqLvl - baseLvl;
    ivec2 normFreqLvlK = bk.freqLvl - baseLvl;

    // skip if frequency difference tolerance is not met.
    int maxFreqDiffIK, maxFreqDiffIJ;
#if SIM3D
    maxFreqDiffIK = max3<int>(
        abs(bi.freqLvl.x - bk.freqLvl.x),
        abs(bi.freqLvl.y - bk.freqLvl.y),
        abs(bi.freqLvl.z - bk.freqLvl.z)
        );
    maxFreqDiffIJ = max3<int>(
        abs(bi.freqLvl.x - bj.freqLvl.x),
        abs(bi.freqLvl.y - bj.freqLvl.y),
        abs(bi.freqLvl.z - bj.freqLvl.z)
        );
#else
    maxFreqDiffIK = glm::max<int>(
        abs(bi.freqLvl.x - bk.freqLvl.x),
        abs(bi.freqLvl.y - bk.freqLvl.y)
        );
    maxFreqDiffIJ = glm::max<int>(
        abs(bi.freqLvl.x - bj.freqLvl.x),
        abs(bi.freqLvl.y - bj.freqLvl.y)
        );
#endif
    if (maxFreqDiffIK > toleranceACoeffFreqDiff || maxFreqDiffIJ > toleranceACoeffFreqDiff) {
        return 0;
    }

    float baseFreq = powf(2.f, float(baseLvl));
#if SIM3D
    vec3 relativeOffsetJ = baseFreq * vec3(bj.center.x - bi.center.x,
        bj.center.y - bi.center.y,
        bj.center.z - bi.center.z);
    vec3 relativeOffsetK = baseFreq * vec3(bk.center.x - bi.center.x,
        bk.center.y - bi.center.y,
        bk.center.z - bi.center.z);
#else
    vec2 relativeOffsetJ = baseFreq * vec2(bj.center.x - bi.center.x,
        bj.center.y - bi.center.y);
    vec2 relativeOffsetK = baseFreq * vec2(bk.center.x - bi.center.x,
        bk.center.y - bi.center.y);
#endif

    // snap offset to very fine grid to avoid float errors
#if SIM3D
    vec3 snappedRelativeOffsetJ = vec3(roundToMultiple(relativeOffsetJ.x, coeffSnapSize),
        roundToMultiple(relativeOffsetJ.y, coeffSnapSize),
        roundToMultiple(relativeOffsetJ.z, coeffSnapSize));
    vec3 snappedRelativeOffsetK = vec3(roundToMultiple(relativeOffsetK.x, coeffSnapSize),
        roundToMultiple(relativeOffsetK.y, coeffSnapSize),
        roundToMultiple(relativeOffsetK.z, coeffSnapSize));
#else
    vec2 snappedRelativeOffsetJ = vec2(roundToMultiple(relativeOffsetJ.x, coeffSnapSize),
        roundToMultiple(relativeOffsetJ.y, coeffSnapSize));
    vec2 snappedRelativeOffsetK = vec2(roundToMultiple(relativeOffsetK.x, coeffSnapSize),
        roundToMultiple(relativeOffsetK.y, coeffSnapSize));
#endif

    // return coefficient if already computed, or else compute it and store it.
#if SIM3D
    KeyTypeA key = make_tuple(normFreqLvlI.x, normFreqLvlI.y, normFreqLvlI.z, bi.axis,
        normFreqLvlJ.x, normFreqLvlJ.y, normFreqLvlJ.z, bj.axis,
        snappedRelativeOffsetJ.x, snappedRelativeOffsetJ.y, snappedRelativeOffsetJ.z,
        normFreqLvlK.x, normFreqLvlK.y, normFreqLvlK.z, bk.axis,
        snappedRelativeOffsetK.x, snappedRelativeOffsetK.y, snappedRelativeOffsetK.z);
#else
    KeyTypeA key = make_tuple(normFreqLvlI.x, normFreqLvlI.y,
        normFreqLvlJ.x, normFreqLvlJ.y, snappedRelativeOffsetJ.x, snappedRelativeOffsetJ.y,
        normFreqLvlK.x, normFreqLvlK.y, snappedRelativeOffsetK.x, snappedRelativeOffsetK.y);
#endif
    double time1 = elapsedTime();
    auto coeffPair = coeffsA.find(key);
    double time2 = elapsedTime();
    tempTime += time2 - time1;

    float result;

    if (coeffPair != coeffsA.end()) {
        result = coeffPair->second;
    }
    else {
        float coeff;

        // compute coefficient with numerical integration
#if SIM3D
        BasisFlow tempBi(normFreqLvlI, vec2(0), bi.axis);
        BasisFlow tempBj(normFreqLvlJ, relativeOffsetJ, bj.axis);
        BasisFlow tempBk(normFreqLvlK, relativeOffsetK, bk.axis);
#else
        BasisFlow tempBi(normFreqLvlI, vec2(0));
        BasisFlow tempBj(normFreqLvlJ, relativeOffsetJ);
        BasisFlow tempBk(normFreqLvlK, relativeOffsetK);
        //            coeff = integrateBasisGradBasis(tempBi, tempBj, tempBk);  
#endif
        coeff = integrateBasisGradBasis(tempBi, tempBj, tempBk);
        //            coeff = 1.0; // TODO: remove this!

        coeffsA.insert(std::pair<KeyTypeA, float>(key, coeff));
        result = coeff;

        if (coeffsA.size() % 1000 == 0) {
            cout << "coeffs A : " << coeffsA.size() << endl;
        }

        newACoeffComputed = true;
    }

    // scaled coefficient
#if SIM3D
    return result * (baseFreq);
#else
    return result * sqr(baseFreq);
#endif

}


#endif



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
        cout << "unknown basis parameters" << endl;
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
        cout << "unknown basis parameters" << endl;
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
vec2 AppBasisFluid::translatedBasisEval(
    const vec2 p,
    const ivec2 freqLvl,
    const vec2 center) const
{
#if SAFETY_ASSERTS
    if (abs(int(freqLvl.x) - int(freqLvl.y)) > maxAnisoLvl) { cout << "translateBasisEval: frequency out of bounds." << endl; return vec2(0, 0); }
#endif

    vec2 result;

    int minLvl = glm::min<int>(freqLvl.x, freqLvl.y);
    if (freqLvl.x <= freqLvl.y) {
        result = basisFlowTemplates[freqLvl.y - freqLvl.x]->interp(
            vec2(float(1 << minLvl)*(p.x - center.x) / lengthLvl0,
                float(1 << minLvl)*(p.y - center.y) / lengthLvl0)
        );
    }
    else {
        // reverse coordinates
        result = -basisFlowTemplates[freqLvl.x - freqLvl.y]->interp(
            vec2((1 << minLvl)*(p.y - center.y) / lengthLvl0,
            (1 << minLvl)*(p.x - center.x) / lengthLvl0)
        );
        result = vec2(result.y, result.x);
    }
    return float(1 << minLvl)*result;
}


// evaluates without using stored grid, using double precision.
dvec2 TranslatedBasisEvalPrecise(const dvec2 p, const ivec2 freqLvl, const dvec2 center)
{
    dvec2 result;
    int minLvl = glm::min<int>(freqLvl.x, freqLvl.y);

    if (freqLvl.x <= freqLvl.y) {
        result = flowBasisHat(
            double(1 << minLvl)*(p - center) / double(lengthLvl0),
            freqLvl.y - minLvl
        );
    }
    else {
        // rotated
        dvec2 r = p - center;
        result = flowBasisHat(
            double(1 << minLvl)*dvec2(r.y, -r.x) / double(lengthLvl0),
            freqLvl.x - minLvl
        );
        result = dvec2(-result.y, result.x);
    }
    return double(1 << minLvl)*result;
}



mat2 TranslatedBasisGradEval(
    const vec2 p,
    const ivec2 freqLvl,
    const vec2 center)
{
#if USE_PRECISE_BASIS_EVAL
    return FMATD(translatedBasisGradEvalPrecise(dvec2(p), freqLvl, dvec2(center)));
#else
    vec2 eps = 1.f / (float(nbCells) * vec2(1 << freqLvl) * 2.0f);

    return FMATD(
        (translatedBasisEval(p + vec2(eps.x, 0), freqLvl, center) - translatedBasisEval(p - vec2(eps.x, 0), freqLvl, center)) / (2.f*eps.x),
        (translatedBasisEval(p + vec2(0, eps.y), freqLvl, center) - translatedBasisEval(p - vec2(0, eps.y), freqLvl, center)) / (2.f*eps.y)
    );
#endif

}



dmat2 TranslatedBasisGradEvalPrecise(const dvec2 p, const ivec2 freqLvl, const dvec2 center)
{
    dmat2 result;
    int minLvl = glm::min<int>(freqLvl.x, freqLvl.y);
    if (freqLvl.x <= freqLvl.y) {
        result = flowBasisHatGrad(
            double(1 << minLvl)*(p - center) / double(lengthLvl0),
            freqLvl.y - minLvl
        );
    }
    else {
        // rotated
        dvec2 r = p - center;
        result = -flowBasisHatGrad(
            double(1 << minLvl)*dvec2(r.y, r.x) / double(lengthLvl0),
            freqLvl.x - minLvl
        );
        result = dmat2(result[1][1], result[1][0], result[0][1], result[0][0]);
        //        result = -flowBasisHatGrad(
        //                    double(1<<minLvl)*(p.y-center.y)/lengthLvl0,
        //                    double(1<<minLvl)*(p.x-center.x)/lengthLvl0,
        //                    freqLvl.x-minLvl
        //                );
    }

    return double(sqr(1 << minLvl))*result;
}


#endif



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
#if SIM3D
    vec3 freq(1 << freqLvl.x, 1 << freqLvl.y, 1 << freqLvl.z);
#else
    vec2 freq(1 << freqLvl.x, 1 << freqLvl.y);
#endif
    return 0.5f / freq;
}



vec2 BasisFlow::normalizedPositionInSupport(vec2 p) {
    return (p - center) / (2.f*supportHalfSize());
}


bool BasisFlow::pointIsInSupport(vec2 p) {
    vec2 normalizedPos = normalizedPositionInSupport(p);
#if SIM3D
    return isInClosedInterval(normalizedPos.x, -0.5f, 0.5f) &&
        isInClosedInterval(normalizedPos.y, -0.5f, 0.5f) &&
        isInClosedInterval(normalizedPos.z, -0.5f, 0.5f);
#else
    return isInClosedInterval(normalizedPos.x, -0.5f, 0.5f) &&
        isInClosedInterval(normalizedPos.y, -0.5f, 0.5f);
#endif
}




