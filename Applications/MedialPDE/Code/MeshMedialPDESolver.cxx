#include "MeshMedialPDESolver.h"
#include "MedialAtomGrid.h"
#include <iomanip>
#include <vector>
#include <vnl/vnl_math.h>
#include <vnl/vnl_random.h>
#include <gsl/gsl_multimin.h>
#include <algorithm>

using namespace std;

#define SE(x) if(x != NULL) {delete x; x = NULL; }

template <class T>
void reset_ptr(T* &x)
{
  if(x != NULL)
    { delete x; x = NULL; }
}

void 
MeshMedialPDESolver
::ComputeJacobianConditionNumber()
{
  size_t i, j;

  // Compute the inverse
  double xNormInvSqr = 0.0;
  Vec ei(topology->nVertices, 0.0), ai(topology->nVertices, 0.0);
  for(i = 0; i < topology->nVertices; i++)
    {
    // Get a column of the inverse
    ei[i] = 1.0; xPardiso.Solve(ei.data_block(), ai.data_block()); ei[i] = 0.0;

    // Print "+" to indicate a solver step
    xNormInvSqr += ai.squared_magnitude();
    }

  // Compute the norm of the matrix itself
  double xNormSqr = 0.0;
  for(j = 0; j < A.GetNumberOfSparseValues(); j++)
    xNormSqr += A.GetSparseData()[j] * A.GetSparseData()[j];

  // Return the ratio of the norms
  cout << "CONDITION NUMBER: " << sqrt(xNormSqr * xNormInvSqr);
}

void
MeshMedialPDESolver
::SetMeshTopology(MeshLevel *topology, MedialAtom *inputAtoms)
{
  // We must first set the dimensions of the matrix A. The finite difference
  // equations involving the LBO are specified at internal vertices in the
  // mesh and involve all neighbors of the vertex. Equations involving the
  // gradient use Loop's approximation, which will also use the entire
  // neigborhood of the vertex (except the vertex itself, but the finite
  // difference equation does use it)
  size_t i, j, k;

  // Clean up storage
  Reset();

  // Store the mesh
  this->topology = topology;
  
  // Allocate and fill the row index array
  size_t *xRowIndex = new size_t[topology->nVertices + 1];
  xRowIndex[0] = 0;
  for(i = 0; i < topology->nVertices; i++)
    xRowIndex[i+1] = xRowIndex[i] + topology->GetVertexValence(i) + 1;

  // Get the number of sparse entries in the A matrix and allocate arrays
  size_t nSparseEntries = xRowIndex[topology->nVertices];
  size_t *xColIndex = new size_t[nSparseEntries];
  double *xSparseValues = new double[nSparseEntries];

  // Next we create a temporary array that contains indices of the neighbor
  // vertices in the vertex list (column) and in the walk around the i-th
  // vertex. This temporary array gets sorted so that we have a sorted
  // sparse matrix structure as well as a mapping from walk-index to index in
  // the sparse matrix.
  vector< pair<size_t, size_t> > xTempColIndex(nSparseEntries);

  // Allocate the arrays that map vertices and neighbor relationships into the
  // sparse matrix A. These are needed because columns in A must be sorted and
  // because A contains non-zero diagonal entries
  xMapVertexNbrToA = new size_t[topology->nbr.GetNumberOfSparseValues()];
  xMapVertexToA = new size_t[topology->nbr.GetNumberOfRows()];

  // Process the data associated with each of the vertices
  for(j = 0, i = 0; i < topology->nVertices; i++)
    {
    // Add the vertex itself to the temp index; NOID indicates for later that
    // this entry does not have a corresponding entry in the topology->nbr 
    // sparse array
    xTempColIndex[j++] = make_pair(i, NOID);

    // Create a walk about the vertex
    EdgeWalkAroundVertex it(topology, i);

    // Visit each vertex in the walk
    for(k = 0; !it.IsAtEnd(); ++it, ++k)
      xTempColIndex[j++] = 
        make_pair(it.MovingVertexId(), it.GetPositionInMeshSparseArray());

    // Sort the vertices in the range corresponding to i
    std::sort(
      xTempColIndex.begin() + xRowIndex[i], 
      xTempColIndex.begin() + xRowIndex[i + 1]);

    // Fill in the values of the column index and the vertex-A mapping
    for(k = xRowIndex[i]; k < xRowIndex[i + 1]; k++)
      {
      xColIndex[k] = xTempColIndex[k].first;
      if(xTempColIndex[k].second == NOID)
        xMapVertexToA[i] = k;
      else
        xMapVertexNbrToA[xTempColIndex[k].second] = k;
      }
    }

  // Initialize our immutable matrix
  A.SetArrays(topology->nVertices, topology->nVertices, 
    xRowIndex, xColIndex, xSparseValues);

  // This is sad, but in addition to the immutable matrix A, we have to store a 
  // PARDISO-compatible column and row index array with 1-based indexing
  xPardisoRowIndex = new int[topology->nVertices+1];
  xPardisoColIndex = new int[nSparseEntries];
  for(i = 0; i <= topology->nVertices; i++)
    xPardisoRowIndex[i] = xRowIndex[i] + 1;
  for(j = 0; j < nSparseEntries; j++)
    xPardisoColIndex[j] = xColIndex[j] + 1;

  // Initialize the right hand side and epsilon
  xRHS.set_size(topology->nVertices);
  xEpsilon.set_size(topology->nVertices);

  // Initialize the triangle area array and other such arrays
  xTriangleGeom = new TriangleGeom[topology->triangles.size()];
  xVertexGeom = new VertexGeom[topology->nVertices];

  // Initialize the tangent weights
  xLoopScheme.SetMeshLevel(topology);

  // Copy the pointer to the atoms. The atoms are passed in to this method
  // and are managed by the parent
  xAtoms = inputAtoms;

  // Set the initial values of the F field (it sucks having to do this)
  for(i = 0; i < topology->nVertices; i++)
    xAtoms[i].F = topology->IsVertexInternal(i) ? 1.0 : 0.8;

  // Initialize the weight array for gradient computation
  W = new SMLVec3d[nSparseEntries];
}

MeshMedialPDESolver
::MeshMedialPDESolver()
{
  xTriangleGeom = NULL;
  xVertexGeom = NULL;
  xMapVertexNbrToA = NULL;
  xMapVertexToA = NULL;
  topology = NULL;
  xPardisoColIndex = NULL;
  xPardisoRowIndex = NULL;
  xAtoms = NULL;
  W = NULL;
}

MeshMedialPDESolver
::~MeshMedialPDESolver()
{
  Reset();
}

void
MeshMedialPDESolver
::Reset()
{
  reset_ptr(xTriangleGeom);
  reset_ptr(xVertexGeom);
  reset_ptr(xMapVertexNbrToA);
  reset_ptr(xMapVertexToA);
  reset_ptr(xPardisoColIndex); 
  reset_ptr(xPardisoRowIndex);
  reset_ptr(W);
}

void
MeshMedialPDESolver
::ComputeMeshGeometry(bool flagGradient)
{
  size_t i, j;
  
  // Clear the vertex fan values that will be accumulated later on
  for(i = 0; i < topology->nVertices; i++)
    xVertexGeom[i].xFanArea = 0.0;
  
  // First, we precompute some triangle-wise measurements, specifically the
  // area of each of the mesh triangles. Unfortunately, this requires sqrt...
  for(i = 0; i < topology->triangles.size(); i++)
    {
    // References to the triangle data
    TriangleGeom &TG = xTriangleGeom[i];
    
    // Get the triangle and its three vertices
    SubdivisionSurface::Triangle &t = topology->triangles[i];
    const SMLVec3d &A = xAtoms[t.vertices[0]].X;
    const SMLVec3d &B = xAtoms[t.vertices[1]].X;
    const SMLVec3d &C = xAtoms[t.vertices[2]].X;

    // Compute edge vectors
    SMLVec3d BC = B - C, CA = C - A, AB = A - B;

    // Get the squared lengths of the three segments
    double a = BC.squared_magnitude();
    double b = CA.squared_magnitude();
    double c = AB.squared_magnitude();
    
    // Compute the area of the triange (use lengths)
    TG.xArea = 0.25 * sqrt(
      (a + a - b) * b + (b + b - c) * c + (c + c - a) * a); 

    // Compute the terms of the form (a + b - c)
    double qc = a + b - c, qb = a + c - b, qa = b + c - a;

    // Compute the cotangent of each angle
    double xOneOverFourArea = 0.25 / TG.xArea; 
    TG.xCotangent[0] = xOneOverFourArea * qa;
    TG.xCotangent[1] = xOneOverFourArea * qb;
    TG.xCotangent[2] = xOneOverFourArea * qc;

    // Add the weights to the fan-weight array
    xVertexGeom[t.vertices[0]].xFanArea += xTriangleGeom[i].xArea;
    xVertexGeom[t.vertices[1]].xFanArea += xTriangleGeom[i].xArea;
    xVertexGeom[t.vertices[2]].xFanArea += xTriangleGeom[i].xArea;
      
    // Now compute the derivatives of these terms with respect to each vertex
    if(flagGradient) 
      {
      // Compute the terms of the form qa * (B - C)
      SMLVec3d pa = qa * BC, pb = qb * CA, pc = qc * AB;

      // Compute the derivatives of the area
      double xOneOverEightArea = 0.5 * xOneOverFourArea;
      TG.xAreaGrad[0] = xOneOverEightArea * (pc - pb);
      TG.xAreaGrad[1] = xOneOverEightArea * (pa - pc);
      TG.xAreaGrad[2] = xOneOverEightArea * (pb - pa);

      // Compute intermediates for the derivatives of the cotangent
      double xOneOverTwoArea = xOneOverFourArea + xOneOverFourArea;
      SMLVec3d wAB = AB * xOneOverTwoArea;
      SMLVec3d wBC = BC * xOneOverTwoArea;
      SMLVec3d wCA = CA * xOneOverTwoArea;

      // Compute quantities Cot[A] / Area
      double xOneOverArea = xOneOverTwoArea + xOneOverTwoArea;
      double xCOA_A = TG.xCotangent[0] * xOneOverArea;
      double xCOA_B = TG.xCotangent[1] * xOneOverArea;
      double xCOA_C = TG.xCotangent[2] * xOneOverArea;

      // Compute the cotangent derivatives
      TG.xCotGrad[0][0] =  wAB - wCA - TG.xAreaGrad[0] * xCOA_A;
      TG.xCotGrad[0][1] =        wCA - TG.xAreaGrad[1] * xCOA_A;
      TG.xCotGrad[0][2] = -wAB       - TG.xAreaGrad[2] * xCOA_A;
      
      TG.xCotGrad[1][0] = -wBC       - TG.xAreaGrad[0] * xCOA_B;
      TG.xCotGrad[1][1] =  wBC - wAB - TG.xAreaGrad[1] * xCOA_B;
      TG.xCotGrad[1][2] =        wAB - TG.xAreaGrad[2] * xCOA_B;

      TG.xCotGrad[2][0] =        wBC - TG.xAreaGrad[0] * xCOA_C;
      TG.xCotGrad[2][1] = -wCA       - TG.xAreaGrad[1] * xCOA_C;
      TG.xCotGrad[2][2] =  wCA - wBC - TG.xAreaGrad[2] * xCOA_C;
      }

    // cout << "Triangle geometry at " << i << ": ";
    // cout << "Area = " << xTriangleGeom[i].xArea;
    // cout << " Cot = " << xTriangleGeom[i].xCotangent << endl;
    }

  // Now, compute the tangent vectors for all the points. Tangent vectors are 
  // given by Loop's formula and require iteration around vertices.
  for(i = 0; i < topology->nVertices; i++)
    {
    // Get a reference to the geometry object
    MedialAtom &a = xAtoms[i];

    // Clear the tangent for this vertex
    a.Xu = xLoopScheme.TangentX(0, i, xAtoms);
    a.Xv = xLoopScheme.TangentX(1, i, xAtoms);

    // Set the one get of the geometry descriptor
    a.G.SetOneJet(a.X.data_block(), a.Xu.data_block(), a.Xv.data_block());

    // Compute the normal vector
    a.ComputeNormalVector();
    } 
}

void
MeshMedialPDESolver
::FillNewtonMatrix(const double *phi, bool flagInputChange)
{
  size_t i, j, k;

  // At this point, the structure of the matrix A has been specified. We have
  // to specify the values. This is done one vertex at a time.
  for(i = 0; i < topology->nVertices; i++)
    {
    EdgeWalkAroundVertex it(topology, i);
    if(!it.IsOpen())
      {
      // If there is no input change since the last call (only phi has
      // changed) we don't do any computation here
      if(flagInputChange == false) continue;

      // V is an internal, LBO vertex. It, and each of its neighbors are
      // involved in the finite difference equation. Since this is a
      // differential operator, we can use the fact that all weights must add
      // up to zero 

      // Accumulator for the sum of the weights around the vertex 
      double w_accum = 0.0;

      // The scaling factor applied to cotangent of each angle
      double scale = 1.5 / xVertexGeom[i].xFanArea;

      // Add weight to each vertex
      for( ; !it.IsAtEnd(); ++it)
        {
        // What is the weight of eps_j? Under Desbrun's discretization given
        // in (2.10) of Xu's paper, coefficient of f(p_j) is 
        //   C * (Cot[\alpha_ij] + Cot[\beta_ij], 
        // where \alpha and \beta are angles facing the edge ij and C is equal
        // to 1.5 / A(p_i), the sum of triangle areas around p_i

        // Compute the weight associated with this neighbor vertex
        double cota = xTriangleGeom[it.TriangleAhead()].xCotangent[
          it.OppositeVertexIndexInTriangleAhead()];
        double cotb = xTriangleGeom[it.TriangleBehind()].xCotangent[
          it.OppositeVertexIndexInTriangleBehind()];
        double weight = scale * (cota + cotb);

        // Set the weight for the moving vertex (to be scaled)
        A.GetSparseData()[xMapVertexNbrToA[it.GetPositionInMeshSparseArray()]] = weight;

        // Update accumulators
        w_accum += weight;
        }

      // Set the diagonal entry in A
      A.GetSparseData()[xMapVertexToA[i]] = -w_accum;

      // for(j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
      //  cout << "Open Vertex " << i << ", Neighbor " 
      //    << A.GetColIndex()[j] << " = " << A.GetSparseData()[j] << endl;
      }
    else
      {
      // V is a boundary vertex for which we compute the Riemannian gradient.
      // This gradient is computed using Loop's tangent formula.
      MedialAtom &a = xAtoms[i];

      // Compute the partials of phi in the tangent directions
      double phi_1 = xLoopScheme.GetPartialDerivative(0, i, phi);
      double phi_2 = xLoopScheme.GetPartialDerivative(1, i, phi);
        
      // Multiply by the contravariant tensor to get weights
      double xi_1 = 
        a.G.xContravariantTensor[0][0] * phi_1 + 
        a.G.xContravariantTensor[1][0] * phi_2;

      double xi_2 = 
        a.G.xContravariantTensor[0][1] * phi_1 + 
        a.G.xContravariantTensor[1][1] * phi_2;

      // Add up to get the weights in the sparse matrix
      EdgeWalkAroundVertex it(topology, i);
      for(; !it.IsAtEnd(); ++it)
        {
        size_t j = xMapVertexNbrToA[it.GetPositionInMeshSparseArray()];
        A.GetSparseData()[j] = 2.0 * (
          xLoopScheme.GetNeighborWeight(0, it) * xi_1 +
          xLoopScheme.GetNeighborWeight(1, it) * xi_2);
        }

      /*
      for(j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
        {
        A.GetSparseData()[j] = 
          2.0 * (xTangentWeights[0][j] * xi_1 + xTangentWeights[1][j] * xi_2);
        }
        */

      // Finally, add the weight for the point itself (-4 \epsilon)
      A.GetSparseData()[xMapVertexToA[i]] = - 4.0 + 2.0 * (
        xLoopScheme.GetOwnWeight(0, i) * xi_1 + 
        xLoopScheme.GetOwnWeight(1, i) * xi_2);

      // for(j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
      //  cout << "Closed Vertex " << i << ", Neighbor " 
      //    << A.GetColIndex()[j] << " = " << A.GetSparseData()[j] << endl;
      }
    }
}

double
MeshMedialPDESolver
::ComputeNodeF(size_t i, const double *phi)
{
  if(topology->IsVertexInternal(i))
    {
    // To compute the laplacian of phi, simply use the weights in the corresponding
    // row of sparse matrix A.
    double lap_phi = 0.0;
    for(size_t j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
      lap_phi += A.GetSparseData()[j] * phi[A.GetColIndex()[j]];

    // Now, the right hand side has the form \Delta \phi - \rho
    return lap_phi - xAtoms[i].xLapR;
    }
  else 
    {
    // We compute the gradient of phi directly, so that it is possible to
    // call this method without the prerequisite of calling FillNewtonMatrix
    double Fu = xLoopScheme.GetPartialDerivative(0, i, phi);
    double Fv = xLoopScheme.GetPartialDerivative(1, i, phi);

    // Return value is |grad Phi|^2 - 4 \phi
    return xAtoms[i].G.SquaredGradientMagnitude(Fu, Fv) - 4.0 * phi[i];
    }
}

double
MeshMedialPDESolver
::FillNewtonRHS(const double *phi)
{
  // Loop over all vertices. This method is a little tricky because it uses all 
  // the weights already computed in matrix A to compute the right hand side.
  for(size_t i = 0; i < topology->nVertices; i++)
    xRHS[i] = - ComputeNodeF(i, phi);

  // Return the squared magnitude
  return xRHS.squared_magnitude();
}

void
MeshMedialPDESolver
::ComputeMedialAtoms(const double *soln)
{
  // Loop over all vertices in the mesh
  for(size_t i = 0; i < topology->nVertices; i++)
    {
    // Get the medial atom corresponding to a
    MedialAtom &a = xAtoms[i];

    // Set the phi in the atom
    a.F = soln[i];

    // Compute the partials of phi in the tangent directions
    a.Fu = xLoopScheme.GetPartialDerivative(0, i, soln);
    a.Fv = xLoopScheme.GetPartialDerivative(1, i, soln);

    // Compute the boundary atoms from this geometric information
    a.ComputeBoundaryAtoms(!topology->IsVertexInternal(i));
    }
}

void
MeshMedialPDESolver
::SetInputData(const SMLVec3d *X, const double *rho, const double *phi)
{
  // Fill the X and rho values of the atoms
  for(size_t i = 0; i < topology->nVertices; ++i)
    {
    xAtoms[i].X = X[i];
    xAtoms[i].xLapR = rho[i];
    if(phi)
      xAtoms[i].F = phi[i];
    }
}

int
MeshMedialPDESolver
::TestJacobianAndGradient(double *xInitSoln)
{
  size_t i;

  // Success?
  bool flagSuccess = true;

  // Compute the mesh geometry
  ComputeMeshGeometry(false);

  // Set the solution to the current values of phi (or to the initial solution
  // passed in to the method)
  vnl_vector<double> xSoln(topology->nVertices);
  for(i = 0; i < topology->nVertices; i++)
    xSoln[i] = xInitSoln == NULL ? xAtoms[i].F : xInitSoln[i];

  // Compute the Jacobian matrix
  FillNewtonMatrix(xSoln.data_block(), true);

  // Test whether the Jacobian matrix is really Jacobian
  double eps = 0.0001;
  for(i = 0; i < topology->nVertices; i++)
    {
    for(SparseMat::RowIterator r = A.Row(i); !r.IsAtEnd(); ++r)
      {
      size_t j = r.Column(); double xj = xSoln[j];

      // Compute central difference gradient
      xSoln[j] = xj + eps;
      double f1 = ComputeNodeF(i, xSoln.data_block());
      xSoln[j] = xj - eps;
      double f2 = ComputeNodeF(i, xSoln.data_block());
      xSoln[j] = xj;

      double cd = (f1 - f2) * 0.5 / eps;
      double ad = r.Value();
      if(fabs(ad - cd) > eps)
        {
        flagSuccess = false;
        cout << "Jacobian test failed: " << ad << " <> " << cd << endl;
        }
      }
    }

  // Compute finite difference gradient of 0.5 F.F
  Vec cdGrad(topology->nVertices, 0.0);
  for(i = 0; i < topology->nVertices; i++)
    {
    // Compute the central difference
    double xi = xSoln[i];
    xSoln[i] = xi + eps; double f1 = 0.5 * this->FillNewtonRHS(xSoln.data_block());
    xSoln[i] = xi - eps; double f2 = 0.5 * this->FillNewtonRHS(xSoln.data_block());
    xSoln[i] = xi;
    cdGrad[i] = (f1 - f2) * 0.5 / eps;
    }

  // Compute the gradient using Jacobian
  this->FillNewtonRHS(xSoln.data_block());
  Vec adGrad = - A.MultiplyTransposeByVector(xRHS);

  // Report the difference
  double xGradDiff = (adGrad - cdGrad).inf_norm();
  if(xGradDiff > eps)
    {
    flagSuccess = false;
    cout << "Gradient test failed: " << xGradDiff << endl;
    }

  return flagSuccess ? 0 : -1;
}

/*
void
MeshMedialPDESolver
::BacktrackNRC(
  const Vec &xold,          // Starting point for Newton's method
  double fold,              // Value of the function 0.5 F.F at this point
  const Vec &g,             // The gradient direction at xold
  const Vec &p,             // The Newton step direction
  Vec &x,                   // Output: the new point along Newton direction
  double &f)                // Output: the value of the function at the new point
{


}
*/

void
MeshMedialPDESolver
::SolveEquation(double *xInitSoln, bool flagGradient)
{
  size_t i;

  // Compute the mesh geometry
  ComputeMeshGeometry(flagGradient);

  // Set the solution to the current values of phi (or to the initial solution
  // passed in to the method)
  vnl_vector<double> xSoln(topology->nVertices);
  for(i = 0; i < topology->nVertices; i++)
    xSoln[i] = xInitSoln == NULL ? xAtoms[i].F : xInitSoln[i];

  // Flag indicating that Newton's method converged (without tricks)
  bool flagConverge = false;

  // Run for up to 50 iterations
  for(i = 0; i < 20; i++)
    {
    // First, fill in the A matrix (only a part of it during most iterations)
    FillNewtonMatrix(xSoln.data_block(), i == 0);

    // Now, fill in the right hand side 
    double xStartingRHS = FillNewtonRHS(xSoln.data_block());

    // If the right hand side is close enough to zero, we can exit
    if(xRHS.inf_norm() < 1.0e-10) 
      {
      // Set the converge flag
      flagConverge = true;

      // Exit the loop
      break;
      }

    // During the first iteration perform factorization
    if(i == 0)
      xPardiso.SymbolicFactorization(
        topology->nVertices, xPardisoRowIndex, xPardisoColIndex, A.GetSparseData());

    // In the subsequent iterations, only do the numeric factorization and solve
    xPardiso.NumericFactorization(A.GetSparseData());
    xPardiso.Solve(xRHS.data_block(), xEpsilon.data_block());

    // Add epsilon to the current guess
    xSoln += xEpsilon;
    }

  // If the iteration failed to converge, we may need to backtrack
  if(!flagConverge)
    {
    cout << "{Conv. Fail. " << xRHS.inf_norm() << "}" << endl;
    throw MedialModelException("Convergence Failure");
    }

  // If Newton fails to converge revert to gradient descent
  /*
  if(xPostSolveRHS > 1e-10) 
    {
    cout << "REVERTING TO GRADIENT DESCENT" << endl;

    // Set up the optimizer with target functions
    gsl_multimin_function_fdf my_func;
    my_func.f = &MeshMedialPDESolver::gsl_f;
    my_func.df = &MeshMedialPDESolver::gsl_df;
    my_func.fdf = &MeshMedialPDESolver::gsl_fdf;
    my_func.n = topology->nVertices;
    my_func.params = this;

    // Create the initial solution
    gsl_vector *my_start = gsl_vector_alloc(my_func.n);
    memcpy(my_start->data, xSoln.data_block(), sizeof(double) * my_func.n);

    // Create conjugate gradient minimizer
    gsl_multimin_fdfminimizer *my_min = 
      gsl_multimin_fdfminimizer_alloc( gsl_multimin_fdfminimizer_conjugate_pr, my_func.n);

    // Set up the parameters of the optimizer
    gsl_multimin_fdfminimizer_set(my_min, &my_func, my_start, 0.1, 1e-7);

    // Iterate until we converge
    while(0 == gsl_multimin_fdfminimizer_iterate(my_min))
      {
      cout << gsl_multimin_fdfminimizer_minimum(my_min) << endl;
      }
    cout << endl;

    // Report
    cout << "Conj Grad Final Value: " << gsl_multimin_fdfminimizer_minimum(my_min);
    
    // Get the best solution
    xSoln.copy_in(gsl_multimin_fdfminimizer_x(my_min)->data);
    }
  */

  // Compute the medial atoms
  ComputeMedialAtoms(xSoln.data_block());
}

double
MeshMedialPDESolver::gsl_f(const gsl_vector *x, void *params)
{
  // Get the pointer to the solver
  MeshMedialPDESolver *solver = static_cast<MeshMedialPDESolver *>(params);

  // Compute the function for this data
  return solver->FillNewtonRHS(x->data);
}

void
MeshMedialPDESolver::gsl_fdf(const gsl_vector *v, void *params, double *f, gsl_vector *df)
{
  // Get the pointer to the solver
  MeshMedialPDESolver *solver = static_cast<MeshMedialPDESolver *>(params);

  // Compute the Jacobian matrix and the function value
  solver->FillNewtonMatrix(v->data, false);
  *f = solver->FillNewtonRHS(v->data);

  // Compute the gradient vector
  Vec grad = -2.0 * solver->A.MultiplyTransposeByVector(solver->xRHS);

  // Copy the gradient in
  std::copy(grad.data_block(), grad.data_block() + grad.size(), df->data);
}

void
MeshMedialPDESolver::gsl_df(const gsl_vector *v, void *params, gsl_vector *df)
{
  double dummy;
  MeshMedialPDESolver::gsl_fdf(v, params, &dummy, df);
}

void
MeshMedialPDESolver
::ComputeRHSGradientMatrix()
{
  // Fill the elements of matrix W
  size_t i, j, k;

  // At this point, the structure of the matrix W has been specified. We have
  // to specify the values. This is done one vertex at a time.
  for(i = 0; i < topology->nVertices; i++)
    {
    // Get a reference to the atom, which contains the current solution
    MedialAtom &a = xAtoms[i];

    // Walk around the vertex
    EdgeWalkAroundVertex it(topology, i);
    if(!it.IsOpen())
      {
      // Compute the laplacian of phi (TODO: is it needed?)
      double xLapPhi = 0.0;
      for(j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
        xLapPhi += A.GetSparseData()[j] * xAtoms[A.GetColIndex()[j]].F;
      
      // The first component of the partial derivative is the derivative 
      // of the term (xFanArea). It is relatively easy to compute.
      double xRhoOverFanArea = xLapPhi / xVertexGeom[i].xFanArea;

      // Scaling factor for cotangent derivatives
      double xScale = 1.5 / xVertexGeom[i].xFanArea;
      
      // Reference to the weight for the fixed vertex
      SMLVec3d &wii = W[xMapVertexToA[i]];
      wii.fill(0.0);

      // Get the value of phi at the center
      double phiFixed = a.F;

      for( ; !it.IsAtEnd(); ++it)
        {
        // The weight vector for the moving vertex
        SMLVec3d &wij = W[xMapVertexNbrToA[it.GetPositionInMeshSparseArray()]];

        // Get the triangle index in front and behind
        size_t ta = it.TriangleAhead(), tb = it.TriangleBehind();
        size_t va = it.MovingVertexIndexInTriangleAhead();
        size_t vb = it.MovingVertexIndexInTriangleBehind();
        size_t qa = it.OppositeVertexIndexInTriangleAhead();
        size_t qb = it.OppositeVertexIndexInTriangleBehind();
        size_t pa = it.FixedVertexIndexInTriangleAhead();
        size_t pb = it.FixedVertexIndexInTriangleBehind();

        // Get the phi's at the surrounding vertices
        double phiAhead = xAtoms[it.VertexIdAhead()].F;
        double phiBehind = xAtoms[it.VertexIdBehind()].F;
        double phiMoving = xAtoms[it.MovingVertexId()].F;

        // Assign the first component of the weight
        wij = (phiMoving - phiFixed) *
          (xTriangleGeom[ta].xCotGrad[qa][va] + xTriangleGeom[tb].xCotGrad[qb][vb]);
        wij += (phiAhead - phiFixed) * xTriangleGeom[ta].xCotGrad[va][va];
        wij += (phiBehind - phiFixed) * xTriangleGeom[tb].xCotGrad[vb][vb];

        // Scale the vector accumulated so far
        wij *= xScale;

        // Assign the second half of the value to the vector
        wij -= xRhoOverFanArea * 
          (xTriangleGeom[ta].xAreaGrad[va] + xTriangleGeom[tb].xAreaGrad[vb]);

        // The cental value can be computed based on the fact that all weights must
        // add up to zero, so the sum of the derivatives should also be zero
        // wii += xScale * (phiMoving - phiFixed) *
        //  (xTriangleGeom[ta].xCotGrad[qa][pa] + xTriangleGeom[tb].xCotGrad[qb][pb]);
        // wii -= xRhoOverFanArea * xTriangleGeom[ta].xAreaGrad[pa];
        wii -= wij;
        }
      }
    else
      {
      // In this step, we must compute the weights of nu_j (j in Nhd(i)) as
      // they contribute to the right hand side of the variational PDE. To see
      // the derivation, look in the March 27, 2006 notebook (red), page with three
      // stars in the top right corner.
      
      // Get references to the differential geometry
      double (*g)[2] = a.G.xContravariantTensor;

      // Compute the tensor h
      double p1 = g[0][0] * a.Fu + g[1][0] * a.Fv;
      double p2 = g[0][1] * a.Fu + g[1][1] * a.Fv;
      double h11 = - 2 * p1 * p1, h12 = - 2 * p1 * p2, h22 = - 2 * p2 * p2;
      
      // Compute the vectors Y1 and Y2
      SMLVec3d Y1 = h11 * a.Xu + h12 * a.Xv;
      SMLVec3d Y2 = h12 * a.Xu + h22 * a.Xv;

      // Add up to get the weights in the sparse matrix
      for(EdgeWalkAroundVertex it(topology, i); !it.IsAtEnd(); ++it)
        W[xMapVertexNbrToA[it.GetPositionInMeshSparseArray()]] = 
          xLoopScheme.GetNeighborWeight(0, it) * Y1 +
          xLoopScheme.GetNeighborWeight(1, it) * Y2;

      // Finally, add the weight for the point itself (-4 \epsilon)
      W[xMapVertexToA[i]] = 
        xLoopScheme.GetOwnWeight(0, i) * Y1 + 
        xLoopScheme.GetOwnWeight(1, i) * Y2;
      }
    }
}

void
MeshMedialPDESolver
::ComputeGradient(vector<MedialAtom *> xVariations)
{
  size_t i, q, j;

  // Compute the weight matrix for the right hand side computations
  this->ComputeRHSGradientMatrix();

  // Each variational derivative is a linear system Ax = b. The matrix A
  // should be set correctly from the last call to Solve() and all we have to
  // change are the vectors B.
  vnl_vector<double> rhs(topology->nVertices, 0.0);
  vnl_vector<double> soln(topology->nVertices, 0.0);

  // Precompute common terms for atom derivatives
  MedialAtom::DerivativeTerms *dt = new MedialAtom::DerivativeTerms[topology->nVertices];
  for(i = 0; i < topology->nVertices; i++)
    xAtoms[i].ComputeCommonDerivativeTerms(dt[i]);

  // Compute the derivative for each variation
  for(q = 0; q < xVariations.size(); q++)
    {
    // Get the derivative atoms for this variation
    MedialAtom *dAtoms = xVariations[q];

    // Repeat for each atom
    for(i = 0; i < topology->nVertices; i++) 
      {
      // Set up the initial right hand side value
      rhs[i] = (topology->IsVertexInternal(i)) ? dAtoms[i].xLapR : 0.0;

      // Compute the rest of the right hand side
      for(size_t j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
        rhs[i] -= dot_product(W[j], dAtoms[A.GetColIndex()[j]].X);
      }

    // Solve the partial differential equation (dPhi/dVar)
    xPardiso.Solve(rhs.data_block(), soln.data_block());

    // Compute each atom
    for(i = 0; i < topology->nVertices; i++) 
      {
      // For each atom, compute F, Fu, Fv, Xu and Xv
      MedialAtom &da = dAtoms[i];
      da.F = soln[i];
      da.Fu = xLoopScheme.GetPartialDerivative(0, i, soln.data_block());
      da.Fv = xLoopScheme.GetPartialDerivative(1, i, soln.data_block());
      da.Xu = xLoopScheme.TangentX(0, i, dAtoms);
      da.Xv = xLoopScheme.TangentX(1, i, dAtoms);
      da.Xuu = da.Xuv = da.Xvv = SMLVec3d(0.0); 

      // Compute the metric tensor derivatives of the atom
      xAtoms[i].ComputeMetricTensorDerivatives(dAtoms[i]);

      // Compute the derivatives of the boundary nodes
      xAtoms[i].ComputeBoundaryAtomDerivatives(dAtoms[i], dt[i]);
      }
    }

  // Delete junk
  delete dt;
}


int
MeshMedialPDESolver
::TestPartialDerivatives()
{
  bool flagSuccess = true;
  size_t i, j, k;
  
  // Epsilon for central differences
  double xEpsilon = 0.00001;
  double z = 0.5 / xEpsilon;

  // Create a pair of atom arrays
  MedialAtom *atom1 = new MedialAtom[topology->nVertices];
  MedialAtom *atom2 = new MedialAtom[topology->nVertices];
  
  // Create a pair of alternative solutions
  MeshMedialPDESolver m1, m2;
  m1.SetMeshTopology(topology, atom1);
  m2.SetMeshTopology(topology, atom2);

  // Compute the solution in the current model
  this->SolveEquation(NULL, true);

  // Compute the gradient terms
  this->ComputeRHSGradientMatrix();

  // Random generator
  vnl_random rnd;

  // Create a pair of offsets (random variations)
  vector<SMLVec3d> varX(topology->nVertices);
  vector<double> varRho(topology->nVertices);

  // Create the solution vector (current solution)
  vnl_vector<double> phi(topology->nVertices);

  // Populate the variations
  for(i = 0; i < topology->nVertices; i++)
    {
    // Set the offsets
    varX[i][0] = rnd.drand32(-1.0, 1.0);
    varX[i][1] = rnd.drand32(-1.0, 1.0);
    varX[i][2] = rnd.drand32(-1.0, 1.0);
    varRho[i] = rnd.drand32(-1.0, 1.0);
    
    // Set phi
    phi[i] = xAtoms[i].F;

    // Apply to the atoms
    m1.xAtoms[i].X = xAtoms[i].X + xEpsilon * varX[i];
    m1.xAtoms[i].xLapR = xAtoms[i].xLapR + xEpsilon * varRho[i];
    m1.xAtoms[i].F = phi[i];

    m2.xAtoms[i].X = xAtoms[i].X - xEpsilon * varX[i];
    m2.xAtoms[i].xLapR = xAtoms[i].xLapR - xEpsilon * varRho[i];
    m2.xAtoms[i].F = phi[i];
    }

  // Compute the geometry in all three cases
  m1.ComputeMeshGeometry(false);
  m2.ComputeMeshGeometry(false);

  // Compare the gradient with the results
  for(j = 0; j < topology->triangles.size(); j++)
    {
    // Get the triangles
    TriangleGeom &T = xTriangleGeom[j];
    TriangleGeom &T1 = m1.xTriangleGeom[j];
    TriangleGeom &T2 = m2.xTriangleGeom[j];

    // Look at the triangle areas
    double cdArea = z * (T1.xArea - T2.xArea);
    double adArea = 
      dot_product(T.xAreaGrad[0], varX[topology->triangles[j].vertices[0]]) +
      dot_product(T.xAreaGrad[1], varX[topology->triangles[j].vertices[1]]) +
      dot_product(T.xAreaGrad[2], varX[topology->triangles[j].vertices[2]]);
    if(fabs(cdArea - adArea) > 2.0 * xEpsilon)
      {
      flagSuccess = false;
      cout << "j = " << std::setw(5) << j 
        << "; cdArea = " << std::setw(12) << cdArea 
        << "; adArea = " << std::setw(12) << adArea 
        << "; delta = " << std::setw(12) << fabs(cdArea - adArea) << endl;
      }

    // Look at the cotangents of the angles
    for(k = 0; k < 3; k++)
      {
      double cdCot = z * (T1.xCotangent[k] - T2.xCotangent[k]);
      double adCot =
        dot_product(T.xCotGrad[k][0], varX[topology->triangles[j].vertices[0]]) +
        dot_product(T.xCotGrad[k][1], varX[topology->triangles[j].vertices[1]]) +
        dot_product(T.xCotGrad[k][2], varX[topology->triangles[j].vertices[2]]);

      // We have to allow for some slack with the Cot derivatives because for
      // small angles the third derivative can be quite large and the central
      // difference approximation can be inaccurate
      if(fabs(cdCot - adCot) > 2.0 * xEpsilon * std::max(1.0, fabs(T.xCotangent[k])))
        {
        flagSuccess = false;
        cout << "j = " << std::setw(5) << j 
          << "; cdCot = " << std::setw(12) << cdCot 
          << "; adCot = " << std::setw(12) << adCot 
          << "; delta = " << std::setw(12) << fabs(cdCot - adCot) << endl;
        }
      }
    }

  // Access matrix A in the offset models
  size_t *xRowIndex = A.GetRowIndex();
  size_t *xColIndex = A.GetColIndex();
  double *A1 = m1.A.GetSparseData();
  double *A2 = m2.A.GetSparseData();

  // Perform the gradient computation
  m1.FillNewtonMatrix(phi.data_block(), true);
  m2.FillNewtonMatrix(phi.data_block(), true);

  // Now, compute the derivative of the laplacian at every vertex
  for(i = 0; i < topology->nVertices; i++)
    {
    if(topology->IsVertexInternal(i))
      {
      double L1 = 0, L2 = 0, cdLap = 0, adLap = 0;
      for(j = xRowIndex[i]; j < xRowIndex[i+1]; j++)
        {
        L1 += A1[j] * phi[xColIndex[j]];
        L2 += A2[j] * phi[xColIndex[j]];
        adLap += dot_product(W[j], varX[xColIndex[j]]);
        }
      cdLap = z * (L1 - L2);
      
      if(fabs(cdLap - adLap) > 2.0 * xEpsilon)
        {
        flagSuccess = false;
        cout << "i = " << std::setw(5) << i 
          << "; cdLap = " << std::setw(12) << cdLap 
          << "; adLap = " << std::setw(12) << adLap 
          << "; delta = " << std::setw(12) << fabs(cdLap - adLap) << endl;
        }
      }
    else
      {
      double GM1 = 0, GM2 = 0, cdGM = 0, adGM = 0;
      double Fu = 0, Fv = 0;
      for(j = xRowIndex[i]; j < xRowIndex[i+1]; j++) 
        {
        adGM += dot_product(W[j], varX[xColIndex[j]]);
        }

      // Compute gradient magnitude of phi on both subdivision surfaces
      GM1 = m1.xAtoms[i].G.SquaredGradientMagnitude(xAtoms[i].Fu, xAtoms[i].Fv);
      GM2 = m2.xAtoms[i].G.SquaredGradientMagnitude(xAtoms[i].Fu, xAtoms[i].Fv);
      cdGM = z * (GM1 - GM2);

      if(fabs(cdGM - adGM) > 2.0 * xEpsilon)
        {
        flagSuccess = false;
        cout << "i = " << std::setw(5) << i 
          << "; cdGM = " << std::setw(12) << cdGM 
          << "; adGM = " << std::setw(12) << adGM 
          << "; delta = " << std::setw(12) << fabs(cdGM - adGM) << endl;
        }
      }
    }

  // Clean up
  delete atom1; delete atom2;

  return flagSuccess ? 0 : -1;
}
