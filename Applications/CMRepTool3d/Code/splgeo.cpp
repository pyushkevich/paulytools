#include "splgeo.h"

void SplineGeoMan::eval(double *x,double *y,double *J,double *H) {
    // Parameterization of spline
    float u = (float)x[0];
    float v = (float)x[1];

    // Patch number
    int iPatch = grid->patch(0,grid->getIndexForUV(0,u));
    int jPatch = grid->patch(1,grid->getIndexForUV(1,v));

    // Basis weights
    SMLVec4f Wu[4],Wv[4];

    // Compute the weights
    spline->basisJet(0,iPatch+3,u,Wu);
    spline->basisJet(1,jPatch+3,v,Wv);

    // Create a medial point
    MedialPoint M;

    // Interpolate the spline up to the second order
    spline->interpolateMedialPoint02(iPatch,jPatch,Wu,Wv,M);

    // Store the point
    y[0] = M.F[0];
    y[1] = M.F[1];
    y[2] = M.F[2];

    // Store the Jacobian
    for(int r=0;r<3;r++) {
        y[r] = M.F[r];
        J[SUBS2(2,r,0)] = M.Fu[r];
        J[SUBS2(2,r,1)] = M.Fv[r];
        H[SUBS3(2,2,r,0,0)] = M.Fuu[r];
        H[SUBS3(2,2,r,0,1)] = M.Fuv[r];
        H[SUBS3(2,2,r,1,0)] = M.Fuv[r];
        H[SUBS3(2,2,r,1,1)] = M.Fvv[r];
    }
}