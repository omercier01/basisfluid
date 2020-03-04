
#include "Application.h"

#include <algorithm>
#include <vector>
#include <iostream>
#include <sstream>
#include <filesystem>

using namespace std;

bool Application::Init_BasisFlows() {

    std::cout << "Computing bais coefficients..." << endl;
    PrintTime();

    stringstream ss;

    ss.str("");
    ss << "Data/Coeffs" << "-" << _maxFreqLvl << "-" << _maxAnisoLvl << "-BB.txt";
    LoadCoeffsBB(ss.str());

    ss.str("");
    ss << "Data/Coeffs" << "-" << _maxFreqLvl << "-" << _maxAnisoLvl << "-T.txt";
    LoadCoeffsT(ss.str());


    // compute possible frequencies
    _freqLvls.clear();
    for (int anisoLvl = _minAnisoLvl; anisoLvl <= _maxAnisoLvl; anisoLvl++) {
        if (anisoLvl == 0) {
            for (int iFreq = _minFreqLvl; iFreq <= _maxFreqLvl; iFreq++) {
                _freqLvls.push_back(ivec2(iFreq, iFreq));
            }
        }
        else {
            for (int majorFreq = _minFreqLvl; majorFreq + anisoLvl <= _maxFreqLvl; majorFreq++) {
                _freqLvls.push_back(ivec2(majorFreq, majorFreq + anisoLvl));
                _freqLvls.push_back(ivec2(majorFreq + anisoLvl, majorFreq));
            }
        }
    }

    // order frequencies by wave number
    std::sort(
        _freqLvls.begin(),
        _freqLvls.end(),
        [](ivec2 a, ivec2 b) {return (1 << a.x)*(1 << a.x) + (1 << a.y)*(1 << a.y) < (1 << b.x)*(1 << b.x) + (1 << b.y)*(1 << b.y); });

    // set all bases of each frequency
    _basisFlowParams->resize(0);

    unsigned int nbBasesTested = 0;


    for (int iFreqLvl = 0; iFreqLvl < _freqLvls.size(); iFreqLvl++)
    {

        ivec2 freqLvl = _freqLvls[iFreqLvl];
        ivec2 freq = ivec2(powf(2.f, float(freqLvl.x)), powf(2.f, float(freqLvl.y)));

        vec2 origin = _lengthLvl0 * vec2((1.f / 4.f) / 2.f*(1.f - 1.f / freq.x), (1.f / 4.f) / 2.f*(1.f - 1.f / freq.y));
        vec2 stride = 0.5f*vec2(1.f / freq.x, 1.f / freq.y);
        vec2 extraOffsets[4] = { vec2(0,0), vec2(0.5,0), vec2(0,0.5), vec2(0.5,0.5) };

        int offsetMinX = int(ceil((_domainLeft - origin.x) / _lengthLvl0 / stride.x)) - 1;
        int offsetMaxX = int(floor((_domainRight - origin.x) / _lengthLvl0 / stride.x)) + 1;
        int offsetMinY = int(ceil((_domainBottom - origin.y) / _lengthLvl0 / stride.y)) - 1;
        int offsetMaxY = int(floor((_domainTop - origin.y) / _lengthLvl0 / stride.y)) + 1;


        unsigned int nbOffsets = sizeof(extraOffsets) / sizeof(vec2);
        vector<vector<unsigned int>> newOrthogonalBasisGroups;
        vector<unsigned int> newSameBasisTemplateGroup;
        for (unsigned int iOffset = 0; iOffset < nbOffsets; iOffset++) {
            vector<unsigned int> newGroup;
            newOrthogonalBasisGroups.push_back(newGroup);
        }

        for (int iOffsetX = offsetMinX; iOffsetX <= offsetMaxX; iOffsetX++) {
            for (int iOffsetY = offsetMinY; iOffsetY <= offsetMaxY; iOffsetY++) {

                for (unsigned int iOffset = 0; iOffset < nbOffsets; iOffset++) {

                    vec2 extraOffset = extraOffsets[iOffset];
                    vec2 center = origin + _lengthLvl0 * vec2(iOffsetX*stride.x, iOffsetY*stride.y) + _lengthLvl0 * vec2(extraOffset.x*stride.x, extraOffset.y*stride.y);

                    BasisFlow stretchedBasis_staticOnly = ComputeStretch(BasisFlow(freqLvl, center), true);

                    if (
                        AllBitsSet(stretchedBasis_staticOnly.bitFlags, INTERIOR)

                        ) {
                        _basisFlowParams->appendCpu(BasisFlow(freqLvl, center, uint(_orthogonalBasisGroupIds.size()) + iOffset));

                        // Add to basis groups
                        newOrthogonalBasisGroups[iOffset].push_back(_basisFlowParams->_nbElements - 1);
                        newSameBasisTemplateGroup.push_back(_basisFlowParams->_nbElements - 1);

                    }
                }

            }
        }

        for (unsigned int iOffset = 0; iOffset < nbOffsets; iOffset++) {
            _orthogonalBasisGroupIds.push_back(newOrthogonalBasisGroups[iOffset]);
        }
        _sameBasisTemplateGroupIds.push_back(newSameBasisTemplateGroup);
    }


    unsigned int N = _basisFlowParams->_nbElements;




    // compute basis norm squared
    BasisFlow* basisFlowParamsPointer = _basisFlowParams->getCpuDataPointer();
    for (unsigned int iBasis = 0; iBasis < N; ++iBasis) {
        BasisFlow& b = basisFlowParamsPointer[iBasis];
        b.normSquared = MatBBCoeff(b, b);
    }



    // resize matrix computation vectors
    _vecX->resize(N);
    _vecTemp->resize(N);
    _vecXForces->resize(N);
    _vecXBoundaryForces->resize(N);
    _vecB->resize(N);


    // fill basis centers acceleration structure
    for (unsigned int iBasis = 0; iBasis < N; ++iBasis) {
        BasisFlow b = _basisFlowParams->getCpuData(iBasis);
        uint idX = glm::clamp<int>(
            int(floor((b.center.x - _domainLeft) / (_domainRight - _domainLeft)*_accelBasisRes)),
            0, _accelBasisRes - 1);
        uint idY = glm::clamp<int>(
            int(floor((b.center.y - _domainBottom) / (_domainTop - _domainBottom)*_accelBasisRes)),
            0, _accelBasisRes - 1);
        _accelBasisCentersIds->getCpuData(idX, idY)->push_back(iBasis);
    }


    // precompute intersection bases (including themselves)
    _intersectingBasesIds->resize(N);
    _intersectingBasesSignificantBBIds->resize(N);
    _intersectingBasesIdsTransport->resize(N);

    for (unsigned int i = 0; i < _basisFlowParams->_nbElements; i++) {
        _intersectingBasesIds->setCpuData(i, new vector<unsigned int>);
        _intersectingBasesSignificantBBIds->setCpuData(i, new vector<unsigned int>);
        _intersectingBasesIdsTransport->setCpuData(i, new vector<unsigned int>);
    }


    basisFlowParamsPointer = _basisFlowParams->getCpuDataPointer();

    // precompute basis supports intersections
    std::vector<BasisSupport> basisSupports;
    for (unsigned int iBasis = 0; iBasis < N; ++iBasis) {
        BasisFlow& b = basisFlowParamsPointer[iBasis];
        basisSupports.push_back(b.getSupport());
    }

    // compute basis intersections and transport
    for (unsigned int iBasis1 = 0; iBasis1 < N; ++iBasis1) {
        BasisFlow& b1 = basisFlowParamsPointer[iBasis1];
        BasisSupport& b1Support = basisSupports[iBasis1];

        vec2 b1TransportLimits = b1.supportHalfSize()*_densityMultiplierBasisHalfSize*1.01f;

        for (unsigned int iBasis2 = iBasis1 + 1; iBasis2 < N; ++iBasis2) {
            BasisFlow& b2 = basisFlowParamsPointer[iBasis2];
            BasisSupport& b2Support = basisSupports[iBasis2];

            if (!IntersectionInteriorEmpty(b1Support, b2Support))
            {
                _intersectingBasesIds->getCpuData_noRefresh(iBasis1)->push_back(iBasis2);
                _intersectingBasesIds->getCpuData_noRefresh(iBasis2)->push_back(iBasis1);

                _intersectingBasesSignificantBBIds->getCpuData_noRefresh(iBasis1)->push_back(iBasis2);
                _intersectingBasesSignificantBBIds->getCpuData_noRefresh(iBasis2)->push_back(iBasis1);


                // transport
                if (
                    b2.freqLvl == b1.freqLvl &&
                    abs(b2.center.x - b1.center.x) <= b1TransportLimits.x &&
                    abs(b2.center.y - b1.center.y) <= b1TransportLimits.y
                    ) {
                    _intersectingBasesIdsTransport->getCpuData_noRefresh(iBasis1)->push_back(iBasis2);
                    _intersectingBasesIdsTransport->getCpuData_noRefresh(iBasis2)->push_back(iBasis1);
                }

            }
        }
        // include itself in intersections
        _intersectingBasesIds->getCpuData_noRefresh(iBasis1)->push_back(iBasis1);
        _intersectingBasesSignificantBBIds->getCpuData_noRefresh(iBasis1)->push_back(iBasis1);
        _intersectingBasesIdsTransport->getCpuData_noRefresh(iBasis1)->push_back(iBasis1);



        if ((iBasis1 + 1) % 1000 == 0) { std::cout << "Basis intersection: " << iBasis1 + 1 << "/" << _basisFlowParams->_nbElements << endl; }
    }
    std::cout << "Basis intersection: " << _basisFlowParams->_nbElements << "/" << _basisFlowParams->_nbElements << endl;

    // sort sets by ID number
    for (unsigned int i = 0; i < _basisFlowParams->_nbElements; ++i) {
        std::sort(_intersectingBasesSignificantBBIds->getCpuData(i)->begin(),
            _intersectingBasesSignificantBBIds->getCpuData(i)->end());
    }

    basisSupports.clear();

    // precompute decompressed T coefficients

    std::cout << "computing decompressed coefficients T..." << endl;
    _coeffsTDecompressedIntersections.clear();
    _coeffsTDecompressedIntersections.resize(N);
    for (unsigned int i = 0; i < N; i++) {
        vector<CoeffTDecompressedIntersectionInfo>& intersectionInfos = _coeffsTDecompressedIntersections[i];
        vector<unsigned int>* localIntersectingBasesIds = _intersectingBasesIds->getCpuData(i);
        for (auto itJ = localIntersectingBasesIds->begin(); itJ != localIntersectingBasesIds->end(); ++itJ) {
            vec2 coeff = MatTCoeff(i, *itJ);
            intersectionInfos.push_back(CoeffTDecompressedIntersectionInfo(*itJ, coeff));
        }

        if ((i + 1) % 1000 == 0) {
            std::cout << "decompressed T : " << i + 1 << " / " << N << endl;
        }
    }
    std::cout << "decompressed T : " << N << " / " << N << endl;


    // precompute decompressed BB coefficients
    std::cout << "computing decompressed coefficients BB..." << endl;
    _coeffsBBDecompressedIntersections.clear();
    _coeffsBBDecompressedIntersections.resize(N);
    _coeffBBExplicitTransferSum_abs.clear();
    _coeffBBExplicitTransferSum_abs.resize(N);
    _coeffBBExplicitTransferSum_sqr.clear();
    _coeffBBExplicitTransferSum_sqr.resize(N);
    for (int iRelFreq = 0; iRelFreq < _nbExplicitTransferFreqs; iRelFreq++) {
        _intersectingBasesIdsDeformation[iRelFreq]->resize(N);
        for (unsigned int i = 0; i < _basisFlowParams->_nbElements; i++) {
            _intersectingBasesIdsDeformation[iRelFreq]->setCpuData(i, new vector<CoeffBBDecompressedIntersectionInfo>);
        }
    }


    for (unsigned int i = 0; i < N; i++) {
        vector<CoeffBBDecompressedIntersectionInfo>& intersectionInfos = _coeffsBBDecompressedIntersections[i];
        vector<unsigned int>* localIntersectingBasesIds = _intersectingBasesSignificantBBIds->getCpuData(i);
        float explicitTransferTotalWeight_abs[_nbExplicitTransferFreqs] = { 0 };
        float explicitTransferTotalWeight_sqr[_nbExplicitTransferFreqs] = { 0 };


        ivec2 freqI = _basisFlowParams->getCpuData(i).freqLvl;
        for (auto it = localIntersectingBasesIds->begin(); it != localIntersectingBasesIds->end(); ++it) {
            if (*it == i) { continue; }
            float coeff = float(MatBBCoeff(i, (*it)));
            intersectionInfos.push_back(CoeffBBDecompressedIntersectionInfo((*it), coeff));

            ivec2 freqJ = _basisFlowParams->getCpuData(*it).freqLvl;
            for (int iRelFreq = 0; iRelFreq < _nbExplicitTransferFreqs; iRelFreq++) {
                if (freqJ - freqI == _explicitTransferFreqs[iRelFreq]) {

                    explicitTransferTotalWeight_abs[iRelFreq] += abs(coeff);
                    explicitTransferTotalWeight_sqr[iRelFreq] += Sqr(coeff);

                    _intersectingBasesIdsDeformation[iRelFreq]->getCpuData(i)->push_back(CoeffBBDecompressedIntersectionInfo((*it), coeff));
                }
            }
        }


        for (int iRelFreq = 0; iRelFreq < _nbExplicitTransferFreqs; iRelFreq++) {
            _coeffBBExplicitTransferSum_abs[i].coeffs[iRelFreq] = explicitTransferTotalWeight_abs[iRelFreq];
            _coeffBBExplicitTransferSum_sqr[i].coeffs[iRelFreq] = std::sqrt(explicitTransferTotalWeight_sqr[iRelFreq]);
        }

        if ((i + 1) % 1000 == 0) {
            std::cout << "Decompressed BB : " << i + 1 << " / " << N << endl;
        }
    }
    std::cout << "Decompressed BB : " << N << " / " << N << endl;


    unsigned int minNbBases = -1;
    unsigned int maxNbBases = 0;
    for (unsigned int i = 0; i < N; i++) {
        minNbBases = glm::min(minNbBases, (unsigned int)_coeffsBBDecompressedIntersections[i].size());
        maxNbBases = glm::max(maxNbBases, (unsigned int)_coeffsBBDecompressedIntersections[i].size());
    }


    // set prevBilFlags
    for (unsigned int iBasis = 0; iBasis < N; ++iBasis) {
        BasisFlow& b = basisFlowParamsPointer[iBasis];
        b.prevBitFlags = b.bitFlags;
    }


    ss.str("");
    ss << "Data/Coeffs" << "-" << _maxFreqLvl << "-" << _maxAnisoLvl << "-BB.txt";
    filesystem::create_directory("Data");
    SaveCoeffsBB(ss.str());

    ss.str("");
    ss << "Data/Coeffs" << "-" << _maxFreqLvl << "-" << _maxAnisoLvl << "-T.txt";
    SaveCoeffsT(ss.str());


    std::cout << "Basis setup done." << endl;
    PrintTime();

    return true;
}