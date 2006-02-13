#include "MeshMedialPDESolver.h"
#include <vector>
#include <vnl/vnl_math.h>
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
MeshMedialPDESolver
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
      xTangentWeights[0][idxFirst] = 1.0;
      xTangentWeights[0][idxLast] = -1.0;

      // Now we branch according to the valence
      if(n == 2)
        {
        xTangentWeights[1][idxFirst] = 1.0;
        xTangentWeights[1][idxLast] = 1.0;
        xTangentWeights[1][xMapVertexToA[i]] = -2.0;
        }
      else if (n == 3)
        {
        // Move iterator to the previous (middle) vertex
        it.GoToFirstEdge(); ++it;

        // Get the corresponding A index
        size_t idxMiddle = xMapVertexNbrToA[it.GetPositionInMeshSparseArray()];

        // Set the weights
        xTangentWeights[1][idxMiddle] = 1.0;
        xTangentWeights[1][xMapVertexToA[i]] = -1.0;
        }
      else
        {
        // Compute the angle theta
        double theta = vnl_math::pi / (n - 1);
        xTangentWeights[1][idxFirst] = xTangentWeights[1][idxLast] = sin(theta);

        // Assign the weights to intermediate vertices
        it.GoToFirstEdge(); ++it;
        for(j = 1; j < n - 1; ++j, ++it)
          {
          size_t idxInt = xMapVertexNbrToA[it.GetPositionInMeshSparseArray()];
          xTangentWeights[1][idxInt] = (2.0 * cos(theta) - 2.0) * sin(theta * j);
          }
        }
      }
    }
}

MeshMedialPDESolver
::MeshMedialPDESolver()
{
  xTriangleGeom = NULL;
  xVertexGeom = NULL;
  xMapVertexNbrToA = NULL;
  xMapVertexToA = NULL;
  topology = NULL;
  xTangentWeights[0] = xTangentWeights[1] = NULL;
  xPardisoColIndex = NULL;
  xPardisoRowIndex = NULL;
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
  reset_ptr(xTangentWeights[0]); 
  reset_ptr(xTangentWeights[1]);
  reset_ptr(xPardisoColIndex); 
  reset_ptr(xPardisoRowIndex);
}

void
MeshMedialPDESolver
::ComputeMeshGeometry(const SMLVec3d *X)
{
  size_t i, j;
  
  // Clear the vertex fan values that will be accumulated later on
  for(i = 0; i < topology->nVertices; i++)
    xVertexGeom[i].xFanArea = 0.0;
  
  // First, we precompute some triangle-wise measurements, specifically the
  // area of each of the mesh triangles. Unfortunately, this requires sqrt...
  for(i = 0; i < topology->triangles.size(); i++)
    {
    // Get the triangle and its three vertices
    SubdivisionSurface::Triangle &t = topology->triangles[i];
    const SMLVec3d &X0 = X[t.vertices[0]];
    const SMLVec3d &X1 = X[t.vertices[1]];
    const SMLVec3d &X2 = X[t.vertices[2]];

    // Get the squared lengths of the three segments
    double l0 = squared_distance(X1, X2);
    double l1 = squared_distance(X2, X0);
    double l2 = squared_distance(X0, X1);
    
    // Compute the area of the triange (use lengths)
    xTriangleGeom[i].xArea = 0.25 * sqrt(
      (l0 + l0 - l1) * l1 + (l1 + l1 - l2) * l2 + (l2 + l2 - l0) * l0); 

    // Add the weights to the fan-weight array
    xVertexGeom[t.vertices[0]].xFanArea += xTriangleGeom[i].xArea;
    xVertexGeom[t.vertices[1]].xFanArea += xTriangleGeom[i].xArea;
    xVertexGeom[t.vertices[2]].xFanArea += xTriangleGeom[i].xArea;
    
    // Compute the cotangent of each angle
    xTriangleGeom[i].xCotangent = SMLVec3d(l2 + l1 - l0, l0 + l2 - l1, l1 + l0 - l2);
    xTriangleGeom[i].xCotangent *= 0.25 / xTriangleGeom[i].xArea;

    cout << "Triangle geometry at " << i << ": ";
    cout << "Area = " << xTriangleGeom[i].xArea;
    cout << " Cot = " << xTriangleGeom[i].xCotangent << endl;
    }

  // Now, compute the tangent vectors for all the points. Tangent vectors are 
  // given by Loop's formula and require iteration around vertices.
  for(i = 0; i < topology->nVertices; i++)
    {
    // Get a reference to the geometry object
    VertexGeom &G = xVertexGeom[i];

    // Clear the tangent for this vertex
    G.t1.fill(0.0);
    G.t2.fill(0.0);

    // Computing the tangent is easy. We use the sparse matrix structure from A to do it
    size_t *xRowIndex = A.GetRowIndex(), *xColIndex = A.GetColIndex();
    for(j = xRowIndex[i]; j < xRowIndex[i+1]; j++)
      {
      const SMLVec3d &xn = X[xColIndex[j]];
      G.t1 += xn * xTangentWeights[0][j];
      G.t2 += xn * xTangentWeights[1][j];
      }

    cout << "Vertex " << (i+1) << " Tangent 1 = ";
    for(j = xRowIndex[i]; j < xRowIndex[i+1]; j++)
      cout << xTangentWeights[0][j] << " X[" << xColIndex[j] + 1 << "] + ";
    cout << endl;
    
    cout << "Vertex " << (i+1) << " Tangent 2 = ";
    for(j = xRowIndex[i]; j < xRowIndex[i+1]; j++)
      cout << xTangentWeights[1][j] << " X[" << xColIndex[j] + 1 << "] + ";
    cout << endl;

    // cout << "Vertex " << i << ", Tangent 1: " << G.t1 << endl;
    // cout << "Vertex " << i << ", Tangent 2: " << G.t2 << endl;
    SMLVec3d nrm = vnl_cross_3d(G.t1, G.t2);
    nrm /= nrm.magnitude();
    cout << "Vertex " << i << ", Normal Vector: " << nrm << endl;

    // Now we can compute the covariant metric tensor
    G.gCovariant[0][0] = dot_product(G.t1, G.t1);
    G.gCovariant[1][1] = dot_product(G.t2, G.t2);
    G.gCovariant[0][1] = G.gCovariant[1][0] = dot_product(G.t1, G.t2);

    // The contravariant tensor as well
    double gInv = 1.0 / (G.gCovariant[0][0] * G.gCovariant[1][1] - 
                         G.gCovariant[0][1] * G.gCovariant[0][1]);
    G.gContravariant[0][0] = gInv * G.gCovariant[1][1];
    G.gContravariant[1][1] = gInv * G.gCovariant[0][0];
    G.gContravariant[0][1] = G.gContravariant[1][0] = - gInv * G.gCovariant[0][1];
    }  
}

void
MeshMedialPDESolver
::FillNewtonMatrix(const SMLVec3d *X, const double *phi)
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
        double cota = 
          xTriangleGeom[it.TriangleAhead()].xCotangent[it.OppositeVertexIndexInTriangleAhead()];
        double cotb = 
          xTriangleGeom[it.TriangleBehind()].xCotangent[it.OppositeVertexIndexInTriangleBehind()];
        double weight = scale * (cota + cotb);

        // Set the weight for the moving vertex (to be scaled)
        A.GetSparseData()[xMapVertexNbrToA[it.GetPositionInMeshSparseArray()]] = weight;

        // Update accumulators
        w_accum += weight;
        }

      // Set the diagonal entry in A
      A.GetSparseData()[xMapVertexToA[i]] = -w_accum;

      for(j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
        cout << "Open Vertex " << i << ", Neighbor " << A.GetColIndex()[j] << " = " << A.GetSparseData()[j] << endl;
      }
    else
      {
      // V is a boundary vertex for which we compute the Riemannian gradient.
      // This gradient is computed using Loop's tangent formula.
      VertexGeom &G = xVertexGeom[i];

      // Compute the partials of phi in the tangent directions
      double phi_1 = 0.0, phi_2 = 0.0;
      for(j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
        {
        double phi_i = phi[A.GetColIndex()[j]];
        phi_1 += xTangentWeights[0][j] * phi_i;
        phi_2 += xTangentWeights[1][j] * phi_i;
        }

      // Multiply by the contravariant tensor to get weights
      double xi_1 = G.gContravariant[0][0] * phi_1 + G.gContravariant[1][0] * phi_2;
      double xi_2 = G.gContravariant[1][0] * phi_1 + G.gContravariant[1][1] * phi_2;

      // Add up to get the weights in the sparse matrix
      for(j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
        {
        A.GetSparseData()[j] = 
          2.0 * (xTangentWeights[0][j] * xi_1 + xTangentWeights[1][j] * xi_2);
        }

      // Finally, add the weight for the point itself (-4 \epsilon)
      A.GetSparseData()[xMapVertexToA[i]] += -4.0;

      for(j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
        cout << "Closed Vertex " << i << ", Neighbor " << A.GetColIndex()[j] << " = " << A.GetSparseData()[j] << endl;
      }
    }
}

void
MeshMedialPDESolver
::FillNewtonRHS(const SMLVec3d *X, const double *rho, const double *phi)
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
      xRHS[i] = rho[i] - lap_phi;
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
MeshMedialPDESolver
::SolveEquation(const SMLVec3d *X, const double *rho, const double *phi, double *soln)
{
  // Compute the mesh geometry
  ComputeMeshGeometry(X);

  // Initialize the solver
  UnsymmetricRealPARDISO xPardiso;

  // Set the solution to phi
  vnl_vector<double> xSoln(phi, topology->nVertices);

  // Run for up to 50 iterations
  for(size_t i = 0; i < 10; i++)
    {
    // First, fill in the A matrix
    FillNewtonMatrix(X, xSoln.data_block());

    // Now, fill in the right hand side 
    FillNewtonRHS(X, rho, xSoln.data_block());

    // Check if the right hand side is close enough to zero that we can terminate
    cout << "Iteration: " << i << ", error: " << xRHS.inf_norm() << endl;
    cout << "Sparse Mat: " << A << endl;
    cout << "RHS: " << xRHS << endl;

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

    cout << "Epsilon: " << xEpsilon << endl;
    cout << "Solution: " << xSoln << endl;

    // Perform backtracking (later)
    }
}
