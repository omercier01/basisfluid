
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


void Application::InverseBBMatrixMain(
    unsigned int iRow, double* vecX, double* vecB,
    BasisFlow* basisDataPointer, unsigned int basisBitMask)
{
    float tempX = float(vecB[iRow]);

    vector<CoeffBBDecompressedIntersectionInfo>& intersectionInfos = _coeffsBBDecompressedIntersections[iRow];
    for (const CoeffBBDecompressedIntersectionInfo& inter : intersectionInfos) {
        if (AllBitsSet(basisDataPointer[inter.j].bitFlags, basisBitMask)) {
            tempX -= inter.coeff * float(vecX[inter.j]);
        }
    }

    vecX[iRow] = double(tempX / float(basisDataPointer[iRow].normSquared));
}



// Solves BB * x = b using Gauss-Seidel.
// Since bases are ordered by wavenumber, doing only one step of GS is equivalent to projection by frequency layer.
// The bitMask indicates the bits that must be set for the basis to be used in the inversion.
void Application::InverseBBMatrix(
    DataBuffer1D<double>* vecX,
    DataBuffer1D<double>* vecB,
    unsigned int basisBitMask)
{
    uint n = vecX->_nbElements;

    // get references
    double* vecXPointer = vecX->getCpuDataPointer();
    double* vecTempPointer = _vecTemp->getCpuDataPointer();
    double* vecBPointer = vecB->getCpuDataPointer();
    BasisFlow* basisFlowParamsPointer = _basisFlowParams->getCpuDataPointer();


    // zero x
    for (int i = 0; i < int(n); i++) {
        vecXPointer[i] = 0.;
    }

    uint iIt = 0;
    while (iIt < _maxNbItMatBBInversion) {

        for (int iBasisGroup = 0; iBasisGroup < _orthogonalBasisGroupIds.size(); iBasisGroup++)
        {
            vector<unsigned int> ids = _orthogonalBasisGroupIds[iBasisGroup];
#pragma omp parallel for num_threads(7) shared(vecXPointer,vecBPointer,ids,basisFlowParamsPointer) firstprivate(basisBitMask,minFreq) default(none)
            for (int idid = 0; idid < ids.size(); idid++) {
                int iRow = ids[idid];
                if (
                    AllBitsSet(basisFlowParamsPointer[iRow].bitFlags, basisBitMask)
                    ) {
                    InverseBBMatrixMain(iRow, vecXPointer, vecBPointer, basisFlowParamsPointer, basisBitMask);
                }
            }
        }

        iIt++;
    }

}


dvec2 eigenLaplace(dvec2 p, dvec2 k) {
    return dvec2(
        k.y*sin(M_PI*p.x*k.x)*cos(M_PI*p.y*k.y),
        -k.x*cos(M_PI*p.x*k.x)*sin(M_PI*p.y*k.y)
    );
}


// integrates basis(...) dot a vector field defined on a grid
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
    float sum = 0;


    for (uint i = 0; i <= _integralGridRes; i++) {
        for (uint j = 0; j <= _integralGridRes; j++) {

            vec2 p = vec2(supLeft + float(i) / _integralGridRes * (supRight - supLeft),
                supBottom + float(j) / _integralGridRes * (supTop - supBottom));

            sum += ((i == 0 || i == _integralGridRes) ? 0.5f : 1.f) * ((j == 0 || j == _integralGridRes) ? 0.5f : 1.f) *
                glm::dot(TranslatedBasisEval(p, b.freqLvl, b.center),
                    velField->interp(p));

        }
    }

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
    float sum = 0;


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
    vec2 sum(0);


    for (uint i = 0; i <= _integralGridRes; i++) {
        for (uint j = 0; j <= _integralGridRes; j++) {

            vec2 p = vec2(supLeft + float(i) / _integralGridRes * (supRight - supLeft),
                supBottom + float(j) / _integralGridRes * (supTop - supBottom));

            sum += ((i == 0 || i == _integralGridRes) ? 0.5f : 1.f) * ((j == 0 || j == _integralGridRes) ? 0.5f : 1.f) *
                TranslatedBasisEval(p, bVec.freqLvl, bVec.center);
        }
    }

    // divide by domain size to get averaged value
    return vec2(sum) * (supRight - supLeft)*(supTop - supBottom) / float(Sqr(_integralGridRes)) / ((1.f / float(1 << bSupport.freqLvl.x))*(1.f / float(1 << bSupport.freqLvl.y)));

}




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
        std::cout << "saved T coefficient to " << filename << endl;
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


float Application::MatBBCoeff(int i, int j) {
    BasisFlow b1 = _basisFlowParams->getCpuData(i);
    BasisFlow b2 = _basisFlowParams->getCpuData(j);
    return MatBBCoeff(b1, b2);
}


// compute the integral and store it for future use
float Application::MatBBCoeff(const BasisFlow& b1, const BasisFlow& b2)
{
    if (IntersectionInteriorEmpty(b1.getSupport(), b2.getSupport())) { return 0.0; }

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


//evaluate basis from stored basis templates. Note that in exponent space,
//division is substraction, so lvlY-minLvl is ky/minK
vec2 Application::TranslatedBasisEval(
    const vec2 p,
    const ivec2 freqLvl,
    const vec2 center)
{

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

