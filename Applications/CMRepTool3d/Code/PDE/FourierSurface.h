#ifndef __FourierSurface_h_
#define __FourierSurface_h_

/**
 * \class FourierSurfaceOld
 * \brief Fourier basis representation for a surface 
 */
class FourierSurfaceOld
{
public:
  /** Construct the Fourier Surface */
  FourierSurfaceOld(unsigned int nBasesU, unsigned int nBasesV, 
    unsigned int nComponents = 3);

  /** Set the coefficients of the surface */
  void SetCoefficient(unsigned int iOrderU, unsigned int iOrderV, 
    unsigned int iComponent, double xValue)
    { xCoeff[iOrderU][iOrderV][iComponent] = xValue; }

  /** Get a coefficient */
  double GetCoefficient(unsigned int iOrderU, unsigned int iOrderV, unsigned int iComponent)
    { return xCoeff[iOrderU][iOrderV][iComponent]; }

  /** Get the number of coefficients */
  unsigned int GetNumberOfCoefficients(unsigned int dim)
    { return dim == 0 ? nBasesU : nBasesV; }

  /** Get the number of components */
  unsigned int GetNumberOfComponents()
    { return nComponents; }

  /** Evaluate the function at a position in space */
  double Evaluate(double u, double v, 
    unsigned int iDerivativeU, unsigned int iDerivativeV, 
    unsigned int iComponent);

  /** Fit the repsesentation to a set of points */
  void FitData(unsigned int iComponent, unsigned int nPoints, 
    double *xPointU, unsigned int iStrideU, 
    double *xPointV, unsigned int iStrideV,  
    double *xValues, unsigned int iStrideX,
    unsigned int nBasesU = 0, unsigned int nBasesV = 0);

protected:  
  /** The Fourier coefficients */
  double ***xCoeff;

  /** Evaluate p-th order basis function at u or v */
  double BasisFunction(unsigned int iOrder, unsigned int iDeriv, double uValue);

  /** Order in dimensions u and v */
  unsigned int nBasesU, nBasesV;

  /** Number of components at each point */
  unsigned int nComponents;

  /** Temp storage arrays */
  double *xBasisTempU, *xBasisTempV;
};
















#endif // __FourierSurface_h_
