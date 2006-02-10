#include "MeshMedialPDESolver.h"
#include <vector>
#include <vnl/vnl_math.h>

using namespace std;

#define LCOS(i,n) cos((i) * 2.0 * vnl_math::pi / (n))
#define LSIN(i,n) cos((i) * 2.0 * vnl_math::pi / (n))

void
MeshMedialPDESolver
::MeshMedialPDESolver()
{
  // Initialize the cosine table if needed
  if(xCosTable.rows() == 0)
    {
    xCosTable.set_size(10, 10);
    xSinTable.set_size(10, 10);
    for(size_t i = 0; i < 10; i++)
      for(size_t j = 0; j < 10; j++)
        {
        xCosTable[i][j] = cos(j * 2.0 * vnl_math::pi / (i + 1)); 
        xSinTable[i][j] = sin(j * 2.0 * vnl_math::pi / (i + 1)); 
        }
    }
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
    xRowIndex[i+1] = xRowIndex[i] + topology->GetVertexValence() + 1;

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
  xMapVertexToA = new size_t[topology->nbr.GetNumberOfRow()];

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

  // Initialize the triangle area array and other such arrays
  xTriangleGeom = new TriangleGeom[topology->triangles.size()];
  xVertexGeom = new VertexGeom[topology->nVertices];
}

MeshMedialPDESolver
::MeshMedialPDESolver()
{
  xTriangleGeom = NULL;
  xVertexGeom = NULL;
  xMapVertexNbrToA = NULL;
  xMapVertexToA = NULL;
  topology = NULL;
}

void
MeshMedialPDESolver
::Reset()
{
  // Delete all pointers
  if(xTriangleGeom) delete xTriangleGeom; xTriangleGeom = NULL;
  if(xVertexGeom) delete xVertexGeom; xVertexGeom = NULL;
  if(xMapVertexNbrToA) delete xMapVertexNbrToA; xMapVertexNbrToA = NULL;
  if(xMapVertexToA) delete xMapVertexToA; xMapVertexToA = NULL;
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
  for(i = 0, i < topology->triangles.size(); i++)
    {
    // Get the triangle and its three vertices
    Triangle &t = topology->triangles[i];
    SMLVec3d &X0 = X[t->vertices[0]];
    SMLVec3d &X1 = X[t->vertices[1]];
    SMLVec3d &X2 = X[t->vertices[2]];

    // Get the squared lengths of the three segments
    double l0 = (X2-X1).magnitude_sqr();
    double l1 = (X2-X0).magnitude_sqr();
    double l2 = (X1-X0).magnitude_sqr();
    
    // Compute the area of the triange (use lengths)
    xTriangleGeom[i].xArea = 0.25 * sqrt(
      (l0 + l0 - l1) * ll + (l1 + l1 - l2) * l2 + (l2 + l2 - l0) * l0); 

    // Add the weights to the fan-weight array
    xVertexGeom[t.vertices[0]].xFanArea += xTriangleGeom[i].xArea[i];
    xVertexGeom[t.vertices[1]].xFanArea += xTriangleGeom[i].xArea[i];
    xVertexGeom[t.vertices[2]].xFanArea += xTriangleGeom[i].xArea[i];
    
    // Compute the cotangent of each angle
    xTriangleGeom[i].xCotangent[i] = SMLVec3d(l2 + l1 - l0, l0 + l2 - l1, l1 + l0 - l2);
    xTriangleGeom[i].xCotangent[i] *= 0.25 / xTriangleGeom[i].xArea[i];
    }

  // Now, compute the tangent vectors for all the points. Tangent vectors are 
  // given by Loop's formula and require iteration around vertices.
  for(i = 0; i < topology->nVertices; i++)
    {
    // Clear the tangent for this vertex
    xVertexGeom[i].t1.fill(0.0);
    xVertexGeom[i].t2.fill(0.0);
    
    // Get the walk around the vertex
    EdgeWalkAroundVertex it(topology, i);
    size_t n = it->Valence();
    if(!it.IsOpen())
      {
      for(j = 0; !it.IsAtEnd(); ++it, ++j)
        {
        xVertexGeom[i].t1 += X[it.MovingVertexId()] * xCosTable[n-1][j];
        xVertexGeom[i].t2 += X[it.MovingVertexId()] * xSinTable[n-1][j];
        }
      }
    
    }
  
  
}

void
MeshMedialPDESolver
::SolveEquation(const SMLVec3d *X, const double *rho, double *soln)
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
      double scale = 1.5 / xVertexFanArea[i];
      
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
          xCotangent[it.TriangleAhead()][it.OppositeVertexIndexInTriangleAhead()];
        double cotb = 
          xCotangent[it.TriangleBehind()][it.OppositeVertexIndexInTriangleBehind()];
        double weight = scale * (cota + cotb);

        // Set the weight for the moving vertex (to be scaled)
        A.GetSparseData()[xMapVertexNbrToA[it.MovingVertexId()]] = weight;

        // Update accumulators
        w_accum += weight;
        }

      // Set the diagonal entry in A
      A.GetSparseData()[xMapVertexToA[i]] = -w_accum;
      }
    else
      {
      // V is a boundary vertex for which we compute the Riemannian gradient.
      // This gradient is computed using Loop's tangent formula.
      
      }
    }


}
