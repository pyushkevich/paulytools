#ifndef __FourierSurface_h_
#define __FourierSurface_h_

/**
 * \class FourierSurface
 * \brief Fourier basis representation for a surface 
 */
class FourierSurface
{
public:
  /** Construct the Fourier Surface */
  FourierSurface(unsigned int nBasesU, unsigned int nBasesV, 
    unsigned int nComponents = 3);

  /** Set the coefficients of the surface */
  void SetCoefficient(unsigned int iOrderU, unsigned int iOrderV, 
    unsigned int iComponent, double xValue)
    { xCoeff[iOrderU][iOrderV][iComponent] = xValue; }

  /** Evaluate the function at a position in space */
  float Evaluate(float u, float v, 
    unsigned int iDerivativeU, unsigned int iDerivativeV, 
    unsigned int iComponent);

  /** Fit the repsesentation to a set of points */
  void FitData(unsigned int iComponent, unsigned int nPoints, 
    float *xPointU, unsigned int iStrideU, 
    float *xPointV, unsigned int iStrideV,  
    float *xValues, unsigned int iStrideX,
    unsigned int nBasesU = 0, unsigned int nBasesV = 0);

protected:  
  /** The Fourier coefficients */
  float ***xCoeff;

  /** Evaluate p-th order basis function at u or v */
  float BasisFunction(unsigned int iOrder, unsigned int iDeriv, float uValue);

  /** Order in dimensions u and v */
  unsigned int nBasesU, nBasesV;

  /** Number of components at each point */
  unsigned int nComponents;

  /** Temp storage arrays */
  float *xBasisTempU, *xBasisTempV;
};

#endif // __FourierSurface_h_
