/******************************************************************
 * MMODEL Library                                                 *
 ******************************************************************
 * Author:						Paul Yushkevich
 *
 * Date:							Feb 15, 1999
 *
 * Description					Medial model algorithms and data structures
 *
 *
 *	Sources:                Various papers on medial models at UNC,
 *                         Pythagoreah hodographs, etc.
 *
 * Dependencies:				PY Matrix, Optima, Registry libs, CLAPACK
 ******************************************************************
 * cboundariness.h
 *	---------------
 * Continous boundariness measure
 ******************************************************************/
#ifndef _CBOUNDARINESS_H_
#define _CBOUNDARINESS_H_

#include "mmodel.h"
#include <matrix/src/matrix.h>
#include <cimage/include/cimage.h>

#ifdef WIN32 
#include <map>
#else
#include <map.h>
#endif

// Begin namespace
NAMESPACE_PAULY_START





double getContinuousBoundariness(MedialInterpolant &mi,
										 cimage im,
										 double steps = 100,
										 double rho = 0.1);

double getContinuousEndness(EndInterpolant &mi,
										 cimage im,
										 double steps = 100,
										 double rho = 0.1);

// End namespace
NAMESPACE_PAULY_END

#endif