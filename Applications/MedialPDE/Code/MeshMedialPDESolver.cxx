#include "MeshMedialPDESolver.h"
#include <vector>
#include <vnl/vnl_math.h>
#include <vnl/vnl_random.h>
#include "PardisoInterface.h"

using namespace std;

#define SE(x) if(x != NULL) {delete x; x = NULL; }

template <class T>
void reset_ptr(T* &x)
{
  if(x != NULL)
    { delete x; x = NULL; }
}

void
MeshMedialModel
::SetMeshTopology(MeshLevel *topology)
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
  xAtoms = new MedialAtom[topology->nVertices];

  // Initialize the tangent weights
  xTangentWeights[0] = new double[nSparseEntries]; 
  xTangentWeights[1] = new double[nSparseEntries];
  fill_n(xTangentWeights[0], nSparseEntries, 0.0);
  fill_n(xTangentWeights[1], nSparseEntries, 0.0);

  // Compute tangent weight at each vertex
  for(i = 0; i < topology->nVertices; i++)
    {
    // Create an iterator around the vertex
    EdgeWalkAroundVertex it(topology, i);

    // Get the valence of the vertex
    size_t n = it.Valence();

    // If the vertex is internal, weights are easy to compute
    if(!it.IsOpen()) 
      {
      for(j = 0; !it.IsAtEnd(); ++it, ++j)
        {
        size_t idxA = xMapVertexNbrToA[it.GetPositionInMeshSparseArray()];
        xTangentWeights[0][idxA] = cos(j * vnl_math::pi * 2.0 / n);
        xTangentWeights[1][idxA] = sin(j * vnl_math::pi * 2.0 / n);
        }
      }
    else 
      {
      // Get the index of the first neighbor
      size_t idxFirst = xMapVertexNbrToA[it.GetPositionInMeshSparseArray()];

      // Move the index to the last neighbor
      it.GoToLastEdge();
      size_t idxLast = xMapVertexNbrToA[it.GetPositionInMeshSparseArray()];

      // We can now set the along-the-edge tangent vector
      xTangentWeights[1][idxFirst] = 1.0;
      xTangentWeights[1][idxLast] = -1.0;

      // Now we branch according to the valence
      if(n == 2)
        {
        xTangentWeights[0][idxFirst] = -1.0;
        xTangentWeights[0][idxLast] = -1.0;
        xTangentWeights[0][xMapVertexToA[i]] = 2.0;
        }
      else if (n == 3)
        {
        // Move iterator to the previous (middle) vertex
        it.GoToFirstEdge(); ++it;

        // Get the corresponding A index
        size_t idxMiddle = xMapVertexNbrToA[it.GetPositionInMeshSparseArray()];

        // Set the weights
        xTangentWeights[0][idxMiddle] = -1.0;
        xTangentWeights[0][xMapVertexToA[i]] = 1.0;
        }
      else
        {
        // Compute the angle theta
        double theta = vnl_math::pi / (n - 1);
        xTangentWeights[0][idxFirst] = xTangentWeights[0][idxLast] = sin(theta);

        // Assign the weights to intermediate vertices
        it.GoToFirstEdge(); ++it;
        for(j = 1; j < n - 1; ++j, ++it)
          {
          size_t idxInt = xMapVertexNbrToA[it.GetPositionInMeshSparseArray()];
          xTangentWeights[0][idxInt] = (2.0 * cos(theta) - 2.0) * sin(theta * j);
          }
        }
      }
    }

  // Initialize the weight array for gradient computation
  W = new SMLVec3d[nSparseEntries];

  // Last step: set the atom iteration context
  xIterationContext = new SubdivisionSurfaceMedialIterationContext(topology);
}

MeshMedialModel
::MeshMedialModel()
{
  xTriangleGeom = NULL;
  xVertexGeom = NULL;
  xMapVertexNbrToA = NULL;
  xMapVertexToA = NULL;
  topology = NULL;
  xTangentWeights[0] = xTangentWeights[1] = NULL;
  xPardisoColIndex = NULL;
  xPardisoRowIndex = NULL;
  W = NULL;
}

MeshMedialModel
::~MeshMedialModel()
{
  Reset();
}

void
MeshMedialModel
::Reset()
{
  reset_ptr(xTriangleGeom);
  reset_ptr(xVertexGeom);
  reset_ptr(xAtoms);
  reset_ptr(xMapVertexNbrToA);
  reset_ptr(xMapVertexToA);
  reset_ptr(xTangentWeights[0]); 
  reset_ptr(xTangentWeights[1]);
  reset_ptr(xPardisoColIndex); 
  reset_ptr(xPardisoRowIndex);
  reset_ptr(xIterationContext);
  reset_ptr(W);
}

void
MeshMedialModel
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
    a.Xu.fill(0.0);
    a.Xv.fill(0.0);

    // Computing the tangent is easy. We use the sparse matrix structure from A to do it
    size_t *xRowIndex = A.GetRowIndex(), *xColIndex = A.GetColIndex();
    for(j = xRowIndex[i]; j < xRowIndex[i+1]; j++)
      {
      const SMLVec3d &xn = xAtoms[xColIndex[j]].X;
      a.Xu += xn * xTangentWeights[0][j];
      a.Xv += xn * xTangentWeights[1][j];
      }
    
    /*
    cout << "Vertex " << (i+1) << " Tangent 1 = ";
    for(j = xRowIndex[i]; j < xRowIndex[i+1]; j++)
      cout << xTangentWeights[0][j] << " X[" << xColIndex[j] << "] + ";
    cout << endl;
    
    cout << "Vertex " << (i+1) << " Tangent 2 = ";
    for(j = xRowIndex[i]; j < xRowIndex[i+1]; j++)
      cout << xTangentWeights[1][j] << " X[" << xColIndex[j] << "] + ";
    cout << endl;

    // cout << "Vertex " << i << ", Tangent 1: " << G.t1 << endl;
    // cout << "Vertex " << i << ", Tangent 2: " << G.t2 << endl;
    SMLVec3d nrm = vnl_cross_3d(G.t1, G.t2);
    nrm /= nrm.magnitude();
    cout << "Vertex " << i << ", Normal Vector: " << nrm << endl;
    */

    // Now we can compute the covariant metric tensor
    a.G.xCovariantTensor[0][0] = dot_product(a.Xu, a.Xu);
    a.G.xCovariantTensor[1][1] = dot_product(a.Xv, a.Xv);
    a.G.xCovariantTensor[0][1] = 
      a.G.xCovariantTensor[1][0] = dot_product(a.Xu, a.Xv);

    // The determinant of the tensor
    a.G.g = 
      a.G.xCovariantTensor[0][0] * a.G.xCovariantTensor[1][1] -
      a.G.xCovariantTensor[0][1] * a.G.xCovariantTensor[1][0];
    a.G.gInv = 1.0 / a.G.g;

    // Compute the contravariant tensor
    a.G.xContravariantTensor[0][0] = a.G.gInv * a.G.xCovariantTensor[1][1];
    a.G.xContravariantTensor[1][1] = a.G.gInv * a.G.xCovariantTensor[0][0];
    a.G.xContravariantTensor[0][1] = 
      a.G.xContravariantTensor[1][0] = - a.G.gInv * a.G.xCovariantTensor[0][1];

    // Compute the normal vector
    a.ComputeNormalVector();
    } 
}

void
MeshMedialModel
::FillNewtonMatrix(const double *phi)
{
  size_t i, j, k;

  // At this point, the structure of the matrix A has been specified. We have
  // to specify the values. This is done one vertex at a time.
  for(i = 0; i < topology->nVertices; i++)
    {
    EdgeWalkAroundVertex it(topology, i);
    if(!it.IsOpen())
      {
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
      double phi_1 = 0.0, phi_2 = 0.0;
      for(j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
        {
        double phi_i = phi[A.GetColIndex()[j]];
        phi_1 += xTangentWeights[0][j] * phi_i;
        phi_2 += xTangentWeights[1][j] * phi_i;
        }

      // Multiply by the contravariant tensor to get weights
      double xi_1 = 
        a.G.xContravariantTensor[0][0] * phi_1 + 
        a.G.xContravariantTensor[1][0] * phi_2;

      double xi_2 = 
        a.G.xContravariantTensor[0][1] * phi_1 + 
        a.G.xContravariantTensor[1][1] * phi_2;

      // Add up to get the weights in the sparse matrix
      for(j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
        {
        A.GetSparseData()[j] = 
          2.0 * (xTangentWeights[0][j] * xi_1 + xTangentWeights[1][j] * xi_2);
        }

      // Finally, add the weight for the point itself (-4 \epsilon)
      A.GetSparseData()[xMapVertexToA[i]] += -4.0;

      // for(j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
      //  cout << "Closed Vertex " << i << ", Neighbor " 
      //    << A.GetColIndex()[j] << " = " << A.GetSparseData()[j] << endl;
      }
    }
}

void
MeshMedialModel
::FillNewtonRHS(const double *phi)
{
  size_t i, j, k;

  // Loop over all vertices. This method is a little tricky because it uses all 
  // the weights already computed in matrix A to compute the right hand side.
  for(i = 0; i < topology->nVertices; i++)
    {
    if(topology->IsVertexInternal(i))
      {
      // To compute the laplacian of phi, simply use the weights in the corresponding
      // row of sparse matrix A.
      double lap_phi = 0.0;
      for(j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
        lap_phi += A.GetSparseData()[j] * phi[A.GetColIndex()[j]];

      // Now, the right hand side has the form \rho - \Delta \phi
      xRHS[i] = xAtoms[i].xLapR - lap_phi;
      }
    else 
      {
      // For an internal vertex, dot product of phi with A gives 
      // 2 \| \Nabla \phi \|^2 - 4 \phi. So if we let that be Z and take
      // 2 phi - 0.5 Z, we get the rhs
      double Z = 0.0;
      for(j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
        Z += A.GetSparseData()[j] * phi[A.GetColIndex()[j]];
      xRHS[i] = 2.0 * phi[i] - 0.5 * Z;
      }
    }
}

void
MeshMedialModel
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
    double phi_1 = 0.0, phi_2 = 0.0;
    for(size_t j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
      {
      double phi_i = soln[A.GetColIndex()[j]];
      phi_1 += xTangentWeights[0][j] * phi_i;
      phi_2 += xTangentWeights[1][j] * phi_i;
      }

    // Multiply by the contravariant tensor to get weights
    double xi_1 = 
      a.G.xContravariantTensor[0][0] * phi_1 + 
      a.G.xContravariantTensor[1][0] * phi_2;

    double xi_2 = 
      a.G.xContravariantTensor[0][1] * phi_1 + 
      a.G.xContravariantTensor[1][1] * phi_2;

    // Compute the gradient vector
    SMLVec3d gradPhi = xi_1 * a.Xu + xi_2 * a.Xv;

    // Compute the medial atom from the gradient vector
    a.ComputeBoundaryAtoms(gradPhi, !topology->IsVertexInternal(i));
    }
}

void
MeshMedialModel
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

void
MeshMedialModel
::SolveEquation(bool flagGradient)
{
  size_t i;

  // Compute the mesh geometry
  ComputeMeshGeometry(flagGradient);

  // Initialize the solver
  UnsymmetricRealPARDISO xPardiso;

  // Set the solution to the current values of phi
  vnl_vector<double> xSoln(topology->nVertices);
  for(i = 0; i < topology->nVertices; i++)
    xSoln[i] = xAtoms[i].F;

  // Run for up to 50 iterations
  for(i = 0; i < 40; i++)
    {
    // First, fill in the A matrix
    FillNewtonMatrix(xSoln.data_block());

    // Now, fill in the right hand side 
    FillNewtonRHS(xSoln.data_block());

    // Check if the right hand side is close enough to zero that we can terminate
    cout << "Iteration: " << i << ", error: " << xRHS.inf_norm() << endl;
    // cout << "Sparse Mat: " << A << endl;
    // cout << "RHS: " << xRHS << endl;

    if(xRHS.inf_norm() < 1.0e-10) break;

    // During the first iteration perform factorization
    if(i == 0)
      xPardiso.SymbolicFactorization(
        topology->nVertices, xPardisoRowIndex, xPardisoColIndex, A.GetSparseData());

    // In the subsequent iterations, only do the numeric factorization and solve
    xPardiso.NumericFactorization(A.GetSparseData());
    xPardiso.Solve(xRHS.data_block(), xEpsilon.data_block());

    // Add epsilon to the current guess
    xSoln += xEpsilon;

    // cout << "Epsilon: " << xEpsilon << endl;
    // cout << "Solution: " << xSoln << endl;

    // Perform backtracking (later)
    }

  // Compute the medial atoms
  ComputeMedialAtoms(xSoln.data_block());
}

void
MeshMedialModel
::ComputeRHSGradientMatrix(const double *phi)
{
  // Fill the elements of matrix W
  size_t i, j, k;

  // At this point, the structure of the matrix W has been specified. We have
  // to specify the values. This is done one vertex at a time.
  for(i = 0; i < topology->nVertices; i++)
    {
    EdgeWalkAroundVertex it(topology, i);
    if(!it.IsOpen())
      {
      // Compute the laplacian of phi (TODO: is it needed?)
      double xLapPhi = 0.0;
      for(j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
        xLapPhi += A.GetSparseData()[j] * phi[A.GetColIndex()[j]];
      
      // The first component of the partial derivative is the derivative 
      // of the term (xFanArea). It is relatively easy to compute.
      double xRhoOverFanArea = xLapPhi / xVertexGeom[i].xFanArea;

      // Scaling factor for cotangent derivatives
      double xScale = 1.5 / xVertexGeom[i].xFanArea;
      
      // Reference to the weight for the fixed vertex
      SMLVec3d &wii = W[xMapVertexToA[i]];
      wii.fill(0.0);

      // Get the value of phi at the center
      double phiFixed = phi[i];

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
        double phiAhead = phi[it.VertexIdAhead()];
        double phiBehind = phi[it.VertexIdBehind()];
        double phiMoving = phi[it.MovingVertexId()];

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
      }
    }
}

void
MeshMedialModel
::ComputeGradient(vector<MedialAtom *> dAtoms)
{
}


void
MeshMedialModel
::TestPartialDerivatives()
{
  bool flagSuccess = true;
  size_t i, j, k;
  
  // Epsilon for central differences
  double xEpsilon = 0.0001;
  double z = 0.5 / xEpsilon;
  
  // Create a pair of alternative solutions
  MeshMedialModel m1, m2;
  m1.SetMeshTopology(topology);
  m2.SetMeshTopology(topology);

  // Random generator
  vnl_random rnd;

  // Create a pair of offsets (random variations)
  vector<SMLVec3d> varX(topology->nVertices);
  vector<double> varRho(topology->nVertices);
  vnl_vector<double> phi(topology->nVertices);
  for(i = 0; i < topology->nVertices; i++)
    {
    // Set the offsets
    varX[i][0] = rnd.drand32(-1.0, 1.0);
    varX[i][1] = rnd.drand32(-1.0, 1.0);
    varX[i][2] = rnd.drand32(-1.0, 1.0);
    varRho[i] = rnd.drand32(-1.0, 1.0);
    phi[i] = rnd.drand32(0.25, 1.0);

    // varX[i] = (i == 239) ? SMLVec3d(1., 0., 0.) : SMLVec3d(0.);

    // Apply to the atoms
    m1.xAtoms[i].X = xAtoms[i].X + xEpsilon * varX[i];
    m2.xAtoms[i].X = xAtoms[i].X - xEpsilon * varX[i];
    m1.xAtoms[i].xLapR = xAtoms[i].xLapR + xEpsilon * varRho[i];
    m2.xAtoms[i].xLapR = xAtoms[i].xLapR - xEpsilon * varRho[i];
    }

  // Compute the geometry in all three cases
  m1.ComputeMeshGeometry(false);
  m2.ComputeMeshGeometry(false);
  this->ComputeMeshGeometry(true);

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
      cout << "j = " << j << "; cdArea = " << cdArea 
        << "; adArea = " << adArea << endl;
      }

    // Look at the cotangents of the angles
    for(k = 0; k < 3; k++)
      {
      double cdCot = z * (T1.xCotangent[k] - T2.xCotangent[k]);
      double adCot =
        dot_product(T.xCotGrad[k][0], varX[topology->triangles[j].vertices[0]]) +
        dot_product(T.xCotGrad[k][1], varX[topology->triangles[j].vertices[1]]) +
        dot_product(T.xCotGrad[k][2], varX[topology->triangles[j].vertices[2]]);
      if(fabs(cdCot - adCot) > 2.0 * xEpsilon)
        {
        flagSuccess = false;
        cout << "j = " << j << "; k = " << k << "; cdCot = " << cdCot 
          << "; adCot = " << adCot << endl;
        }
      }
    }

  // Access matrix A in the offset models
  size_t *xRowIndex = A.GetRowIndex();
  size_t *xColIndex = A.GetColIndex();
  double *A1 = m1.A.GetSparseData();
  double *A2 = m2.A.GetSparseData();

  // Perform the gradient computation
  this->FillNewtonMatrix(phi.data_block());
  this->ComputeRHSGradientMatrix(phi.data_block());
  m1.FillNewtonMatrix(phi.data_block());
  m2.FillNewtonMatrix(phi.data_block());

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
        cout << "i = " << i << "; cdLap = " << cdLap 
          << "; adLap = " << adLap << endl;
        }
      }
    }

}
