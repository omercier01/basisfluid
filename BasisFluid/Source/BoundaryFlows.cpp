
#include "Application.h"
#include "Obstacles.h"

#include <glm/gtc/random.hpp>
#include <glm/ext.hpp>

using namespace glm;

//------------------------------------------------------------------------------
// stretch around boundary

// projects point _x_ onto plane with normal _n_ and member point _p_
vec2 projectOntoPlane(vec2 x, vec2 p, vec2 n) {
    return x - dot(x - p, n)*n;
}

float distanceToPlane(vec2 x, vec2 p, vec2 n) {
    return dot(x - p, n);
}

BasisFlow Application::ComputeStretch(BasisFlow b, bool staticObstaclesOnly, bool noStretch) {
    const float stretchBandRatio = STRETCH_BAND_RATIO; // band of possible corner movement has width _this_ times the basis support half size.

    BasisSupport s = b.getSupport();
    vec2 shs = b.supportHalfSize();

    b.bitFlags = SetBits(b.bitFlags, INTERIOR); // interior basis by default
    b.bitFlags = UnsetBits(b.bitFlags, DYNAMIC_BOUNDARY_PROJECTION);

    for (Obstacle* obs : _obstacles)
    {
        if (staticObstaclesOnly && obs->dynamic) { continue; }

        // if the center of the basis is inside the obstacle, we know the stretch will be too important, so we directly invalidate this basis. This prevents from running into undefined obstacle gradient cases inside the obstacle.
        if (obs->phi(b.center) < 0) {
            b.bitFlags = UnsetBits(b.bitFlags, INTERIOR);
            b.coeff = 0;
            b.coeffBoundary = 0;
        }
    }

    bool basisIsStretched = false;
    bool hasAtLeastOneCornerInside = false;
    bool hasAtLeastOneCornerOutside = false;


    // initialize stretched corners to regular corner positions
    b.stretchedCornerLB = vec2(s.left, s.bottom);
    b.stretchedCornerLT = vec2(s.left, s.top);
    b.stretchedCornerRB = vec2(s.right, s.bottom);
    b.stretchedCornerRT = vec2(s.right, s.top);


    // set basis flags and do first push and stretch

    for (Obstacle* obs : _obstacles)
    {
        if (staticObstaclesOnly && obs->dynamic) { continue; }

        // compute plane approximation of obstacles near basis center
        vec2 centerObsGrad = obs->gradPhi(b.center); // TODO: should this be not normalized, to get a better planar approximation?
        vec2 obsPlanePoint = b.center - obs->phi(b.center)*centerObsGrad;
        vec2 obsPlaneNormal = length(centerObsGrad) < 1e-6 ? centerObsGrad : normalize(centerObsGrad);

        float minOriginalDistToPlane = 99999999.f;
        float maxOriginalDistToPlane = -99999999.f;

        float minStretchedDistToPlane = 99999999.f;
        float maxStretchedDistToPlane = -99999999.f;

        vec2 originalCorner;
        vec2* stretchedCorner;
        for (int i = 0; i < 4; i++) {
            switch (i) {
            case 0:
                stretchedCorner = &b.stretchedCornerLB;
                originalCorner = vec2(s.left, s.bottom); break;
            case 1:
                stretchedCorner = &b.stretchedCornerLT;
                originalCorner = vec2(s.left, s.top); break;
            case 2:
                stretchedCorner = &b.stretchedCornerRB;
                originalCorner = vec2(s.right, s.bottom); break;
            case 3:
                stretchedCorner = &b.stretchedCornerRT;
                originalCorner = vec2(s.right, s.top); break;
            default:
                stretchedCorner = nullptr; break;
            }


            vec2 sp = projectOntoPlane(*stretchedCorner, obsPlanePoint, obsPlaneNormal);
            float dist = distanceToPlane(originalCorner, obsPlanePoint, obsPlaneNormal);

            maxOriginalDistToPlane = glm::max<float>(maxOriginalDistToPlane, dist);
            minOriginalDistToPlane = glm::min<float>(minOriginalDistToPlane, dist);

            if (dist <= 0 && // must be squished
                (
                    abs((sp - originalCorner).x) > shs.x*stretchBandRatio || abs((sp - originalCorner).y) > shs.y*stretchBandRatio
                    )
                ) {
                // stretched point too close to center, basis is invalid
                b.bitFlags = UnsetBits(b.bitFlags, INTERIOR);
                b.coeff = 0;
                b.coeffBoundary = 0;
            }
            else if (dist >= 0 && // must be stretched
                (
                    abs((sp - originalCorner).x) >= shs.x*stretchBandRatio || abs((sp - originalCorner).y) >= shs.y*stretchBandRatio
                    )
                ) {
                // corner is stretched too far away, leave it unstretched.
                maxStretchedDistToPlane = glm::max<float>(maxStretchedDistToPlane, dist);
                minStretchedDistToPlane = glm::min<float>(minStretchedDistToPlane, dist);
            }
            else
            {
                // point is stretched but not too much, keep it, and mark basis as stretched.
                basisIsStretched = true;
                *stretchedCorner = sp;

                float stretchedDist = distanceToPlane(sp, obsPlanePoint, obsPlaneNormal);
                maxStretchedDistToPlane = glm::max<float>(maxStretchedDistToPlane, stretchedDist);
                minStretchedDistToPlane = glm::min<float>(minStretchedDistToPlane, stretchedDist);
            }

            if (obs->dynamic)
            {
                if (dist <= 0) { hasAtLeastOneCornerInside = true; }
                if (dist >= 0) { hasAtLeastOneCornerOutside = true; }
            }

        }

        if (!noStretch && _nbStretchLoops > 0 && basisIsStretched)
        {
            // stretch relatively to basis center in direction orthogonal to obstacle plane
            float stretchRatio = (maxOriginalDistToPlane - minOriginalDistToPlane) / (maxStretchedDistToPlane - minStretchedDistToPlane);

            for (int i = 0; i < 4; i++) {
                switch (i) {
                case 0:
                    stretchedCorner = &b.stretchedCornerLB; break;
                case 1:
                    stretchedCorner = &b.stretchedCornerLT; break;
                case 2:
                    stretchedCorner = &b.stretchedCornerRB; break;
                case 3:
                    stretchedCorner = &b.stretchedCornerRT; break;
                }

                vec2 smc = (*stretchedCorner - b.center);
                *stretchedCorner += (stretchRatio - 1)* (smc - glm::dot(smc, obsPlaneNormal)*obsPlaneNormal);

            }

        }

    }


    if (hasAtLeastOneCornerInside && hasAtLeastOneCornerOutside
        ) {
        b.bitFlags = SetBits(b.bitFlags, DYNAMIC_BOUNDARY_PROJECTION);
    }

    b.stretched = basisIsStretched;


    if (!noStretch && AllBitsSet(b.bitFlags, INTERIOR) && b.stretched)
    {
        for (uint iStretchLoop = 0; iStretchLoop < _nbStretchLoops; iStretchLoop++)
        {

            // push
            for (Obstacle* obs : _obstacles)
            {
                if (staticObstaclesOnly && obs->dynamic) { continue; }

                // compute plane approximation of obstacles near basis center
                vec2 centerObsGrad = obs->gradPhi(b.center); // TODO: should this be not normalized, to get a better planar approximation?
                vec2 obsPlanePoint = b.center - obs->phi(b.center)*centerObsGrad;
                vec2 obsPlaneNormal = length(centerObsGrad) < 1e-6 ? centerObsGrad : normalize(centerObsGrad);

                float minOriginalDistToPlane = 99999999.f;
                float maxOriginalDistToPlane = -99999999.f;

                float minStretchedDistToPlane = 99999999.f;
                float maxStretchedDistToPlane = -99999999.f;

                vec2* stretchedCorner;
                for (int i = 0; i < 4; i++) {
                    switch (i) {
                    case 0:
                        stretchedCorner = &b.stretchedCornerLB; break;
                    case 1:
                        stretchedCorner = &b.stretchedCornerLT; break;
                    case 2:
                        stretchedCorner = &b.stretchedCornerRB; break;
                    case 3:
                        stretchedCorner = &b.stretchedCornerRT; break;
                    default:
                        stretchedCorner = nullptr; break;
                    }


                    vec2 sp = projectOntoPlane(*stretchedCorner, obsPlanePoint, obsPlaneNormal);
                    float dist = distanceToPlane(*stretchedCorner, obsPlanePoint, obsPlaneNormal);

                    maxOriginalDistToPlane = glm::max<float>(maxOriginalDistToPlane, dist);
                    minOriginalDistToPlane = glm::min<float>(minOriginalDistToPlane, dist);

                    if (dist >= 0 && // must be stretched
                        (
                            abs((sp - *stretchedCorner).x) >= shs.x*stretchBandRatio || abs((sp - *stretchedCorner).y) >= shs.y*stretchBandRatio
                            )
                        ) {
                        // corner is stretched too far away, leave it unstretched.
                        maxStretchedDistToPlane = glm::max<float>(maxStretchedDistToPlane, dist);
                        minStretchedDistToPlane = glm::min<float>(minStretchedDistToPlane, dist);
                    }
                    else
                    {
                        // point is stretched but not too much, keep it, and mark basis as stretched.
                        *stretchedCorner = sp;

                        float stretchedDist = distanceToPlane(sp, obsPlanePoint, obsPlaneNormal);
                        maxStretchedDistToPlane = glm::max<float>(maxStretchedDistToPlane, stretchedDist);
                        minStretchedDistToPlane = glm::min<float>(minStretchedDistToPlane, stretchedDist);
                    }

                }


                //stretch (not stretching on last iteration)
                if (iStretchLoop != _nbStretchLoops - 1)
                {

                    // stretch relatively to basis center in direction orthogonal to obstacle plane
                    float stretchRatio = (maxOriginalDistToPlane - minOriginalDistToPlane) / (maxStretchedDistToPlane - minStretchedDistToPlane);

                    for (int i = 0; i < 4; i++) {
                        switch (i) {
                        case 0:
                            stretchedCorner = &b.stretchedCornerLB; break;
                        case 1:
                            stretchedCorner = &b.stretchedCornerLT; break;
                        case 2:
                            stretchedCorner = &b.stretchedCornerRB; break;
                        case 3:
                            stretchedCorner = &b.stretchedCornerRT; break;
                        }

                        vec2 smc = (*stretchedCorner - b.center);
                        *stretchedCorner += (stretchRatio - 1)* (smc - glm::dot(smc, obsPlaneNormal)*obsPlaneNormal);

                    }

                }

            }
        }
    }

    return b;

}





void Application::ComputeStretches()
{
    for (unsigned int iBasis = 0; iBasis < _basisFlowParams->_nbElements; ++iBasis) {
        BasisFlow b = _basisFlowParams->getCpuData(iBasis);
        b = ComputeStretch(b);
        _basisFlowParams->setCpuData(iBasis, b);
    }

}



// Uses Newton iterations to inverse the trilinear coefficients, going from World space (p) to UV space, see http://stackoverflow.com/questions/808441/inverse-bilinear-interpolation .
vec2 Application::QuadCoord(vec2 p, BasisFlow const& b)
{
    // inverse trilinear interpolation with Newton's method
    const int maxIt = NB_NEWTON_ITERATIONS_INVERSION;
    const float tol = 1e-2f; // tolerance to check the iterations have converged (in UV space, not world space)
    vec2 c(0.5); // coordinates to find, with initial guess
    vec2 delta;

    for (int it = 0; it < maxIt; it++) {

        vec2 r =
            (1 - c.x)*(1 - c.y)*b.stretchedCornerLB
            + (1 - c.x)*(c.y)*b.stretchedCornerLT
            + (c.x)*(1 - c.y)*b.stretchedCornerRB
            + (c.x)*(c.y)*b.stretchedCornerRT
            - p;


        delta = inverse(QuadCoordInvDeriv(c, b)) * r;
        c -= delta;
    }

    if (length(delta) < tol) {
        return c;
    }
    else {
        return vec2(9999999999.f);
    }

}


// gives the Jacobian of the deformation form UV space to World space.
mat2 Application::QuadCoordInvDeriv(vec2 uv, BasisFlow const& b)
{
    mat2 J(0);
    J += glm::outerProduct(b.stretchedCornerLB, vec2((-1)*(1 - uv.y), (1 - uv.x)*(-1)));
    J += glm::outerProduct(b.stretchedCornerLT, vec2((-1)*(uv.y), (1 - uv.x)*(1)));
    J += glm::outerProduct(b.stretchedCornerRB, vec2((1)*(1 - uv.y), (uv.x)*(-1)));
    J += glm::outerProduct(b.stretchedCornerRT, vec2((1)*(uv.y), (uv.x)*(1)));

    return J;
}


// evaluates the vector field stretched around the obstacles (weighted by the basis coefficient)
vec2 Application::VecObstacle_stretch(vec2 p, BasisFlow const& b)
{

    vec2 vec;

    if (!AllBitsSet(b.bitFlags, INTERIOR) && AtLeastOneBitNotSet(b.bitFlags, DYNAMIC_BOUNDARY_PROJECTION)) {
        vec = vec2(1);
    }
    else if (!b.stretched) {
        vec = b.coeff * TranslatedBasisEval(p, b.freqLvl, b.center);
    }
    else {

        vec2 uv = QuadCoord(p, b); // World space to UV space (deformed)
        mat2 mat = QuadCoordInvDeriv(uv, b); // Jacobian of deformation from UV to World
        vec2 pos = b.center + b.supportHalfSize() * (uv * 2.f - vec2(1)); // UV space to basis space (undeformed)

        vec = b.coeff *
            mat * // to inverse deformation from World to UV
            (
                TranslatedBasisEval(pos, b.freqLvl, b.center)
                / (2.f * b.supportHalfSize()) // division to inverse deformation from UV to Basis
            );
    }

    return vec;
}



