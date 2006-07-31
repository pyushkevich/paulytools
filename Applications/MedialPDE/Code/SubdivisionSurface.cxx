#include "SubdivisionSurface.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vnl/vnl_vector_fixed.h"
#include <string>
#include <list>
#include <map>

#ifndef vtkFloatingPointType
#define vtkFloatingPointType float
#endif

using namespace std;

SubdivisionSurface::Triangle::Triangle()
{
  for(size_t i = 0; i < 3; i++)
    {
    vertices[i] = NOID;
    neighbors[i] = NOID;
    nedges[i] = -1;
    }
}

void 
SubdivisionSurface
::SetOddVertexWeights(MutableSparseMatrix &W, 
  const MeshLevel *parent, MeshLevel *child, size_t t, size_t v)
{
  // Weight constants
  const static double W_INT_EDGE_CONN = 3.0 / 8.0;
  const static double W_INT_FACE_CONN = 1.0 / 8.0;
  const static double W_BND_EDGE_CONN = 1.0 / 2.0;
  
  // Get the parent and child triangles
  const Triangle &tp = parent->triangles[t / 4];
  Triangle &tc = child->triangles[t];

  // Get the child vertex index
  size_t ivc = tc.vertices[v];

  // Find out if this is a boundary vertex. Since the child triangle is a center
  // triangle, we must check if the parent triangle has a neighbor accross edge v
  if(tp.neighbors[v] != NOID)
    {
    // Get the neighbor of the parent triangle
    const Triangle &topp = parent->triangles[tp.neighbors[v]];
    
    // Internal triangle. It's weights are 3/8 for the edge-connected parent 
    // vertices and 1/8 for the face-connected vertices
    W(ivc,tp.vertices[(v+1)%3]) = W_INT_EDGE_CONN;
    W(ivc,tp.vertices[(v+2)%3]) = W_INT_EDGE_CONN;
    W(ivc,tp.vertices[v]) = W_INT_FACE_CONN;
    W(ivc,topp.vertices[tp.nedges[v]]) = W_INT_FACE_CONN;
    }
  else
    {
    // Only the edge-connected vertices are involved
    W(ivc,tp.vertices[(v+1)%3]) = W_BND_EDGE_CONN;
    W(ivc,tp.vertices[(v+2)%3]) = W_BND_EDGE_CONN;
    }
}

inline bool NeighborPairPredicate (
  const pair<size_t, SubdivisionSurface::NeighborInfo> &p1,
  const pair<size_t, SubdivisionSurface::NeighborInfo> &p2)
{
  return p1.first < p2.first;  
}

void SubdivisionSurface
::ComputeWalks(MeshLevel *mesh)
{
  // Create a temporary mutable structure that represents walk info
  typedef std::pair<size_t, NeighborInfo> Entry;
  typedef std::list<Entry> Row;
  std::vector<Row> walks(mesh->nVertices);

  // Loop over all triangles, all vertices
  for(size_t t = 0; t < mesh->triangles.size(); t++) for(size_t v = 0; v < 3; v++)
    {
    // This is the index of the current vertex
    size_t ivtx = mesh->triangles[t].vertices[v];

    // Check if the walk around the vertex has already been generated
    if(walks[ivtx].size() == 0)
      {
      // The current position in the walk
      size_t tWalk = t; short vWalk = v;

      // The position at which the walk will loop around
      size_t tLoop = t;
      
      // Walk until reaching a NOID triangle or looping around
      // cout << "Walk around vtx. " << ivtx << endl;
      do
        {
        // Get a reference to the current triangle
        Triangle &T = mesh->triangles[tWalk];

        // Represent the next triangle in the walk
        size_t tNext = T.neighbors[ror(vWalk)];
        short vNext = ror(T.nedges[ror(vWalk)]);
        
        // Put the current edge in the triangle into the list
        // cout << "Fwd walk: visiting vtx. " << T.vertices[rol(vWalk)] <<
        //    "; behind: (" << tWalk << "," << vWalk <<
        //    "); ahead: (" << tNext << "," << vNext << ")" <<  endl;
        walks[ivtx].push_back( make_pair(
            T.vertices[rol(vWalk)], 
            NeighborInfo(tNext, vNext, tWalk, vWalk)));
              
        // Update the position in the walk
        tWalk = tNext; vWalk = vNext;
        } 
      while(tWalk != NOID && tWalk != tLoop);

      // Now, if we hit a NOID, we can need to walk in the opposite direction,
      // starting from the same triangle as before
      if(tWalk == NOID)
        {
        // Reset the position in the walk
        tWalk = t; vWalk = v;

        // Walk again, in a symmetrical loop
        do
          {
          // Get a reference to the current triangle
          Triangle &T = mesh->triangles[tWalk];

          // Represent the next triangle in the walk
          size_t tNext = T.neighbors[rol(vWalk)];
          short vNext = rol(T.nedges[rol(vWalk)]);
          
          // Put the current edge in the triangle into the list
          // cout << "Rev walk: visiting vtx. " << T.vertices[ror(vWalk)] <<
          //   "; behind: (" << tNext << "," << vNext <<
          //   "); ahead: (" << tWalk << "," << vWalk << ")" << endl;
          walks[ivtx].push_front( make_pair(
              T.vertices[ror(vWalk)], 
              NeighborInfo(tWalk, vWalk, tNext, vNext)));

          // Update the position in the walk
          tWalk = tNext; vWalk = vNext;
          } 
        while(tWalk != NOID);
        }
      else
        {
        // Rotate the walk so that the smallest member is the first element
        // (for consistency with MMA code)
        rotate(walks[ivtx].begin(), 
          min_element(walks[ivtx].begin(), walks[ivtx].end(), &NeighborPairPredicate), 
          walks[ivtx].end());
        }
      }
    }

  // Now we have visited all the vertices and computed a walk around each one.
  // All that is left is to transfer this result into a sparse matrix
  mesh->nbr.SetFromSTL(walks, mesh->triangles.size());
}

/** 
 * This method should be called for child-level triangles whose index in the
 * parent triangle matches the vertex id (v).
 */
void SubdivisionSurface
::SetEvenVertexWeights(MutableSparseMatrix &W,
  const MeshLevel *parent, MeshLevel *child, size_t t, size_t v)
{
  // Weight constants
  const static double W_BND_EDGE_CONN = 1.0 / 8.0;
  const static double W_BND_SELF = 3.0 / 4.0;
  
  // Get the child vertex index
  size_t ivc = child->triangles[t].vertices[v];

  // Get the corresponding parent triangle and vertex index
  size_t tp = t / 4;
  size_t ivp = parent->triangles[tp].vertices[v];

  // Get the vertex object for this vertex
  EdgeWalkAroundVertex it(parent, ivp);
  if(it.IsOpen()) 
    {
    // Add the point itself
    W(ivc, ivp) = W_BND_SELF;

    // Add the starting point
    W(ivc, it.MovingVertexId()) = W_BND_EDGE_CONN;

    // Get the ending point
    it.GoToLastEdge();
    W(ivc, it.MovingVertexId()) = W_BND_EDGE_CONN;
    }
  else
    {
    // Compute the beta constant
    int n = it.Valence();
    double beta = (n > 3) ? 3.0 / (8.0 * n) : 3.0 / 16.0;

    // Assign beta to each of the edge-adjacent vertices
    for( ; !it.IsAtEnd(); ++it)
      W(ivc, it.MovingVertexId()) = beta;

    // Assign the balance to the coincident vertex
    W(ivc, ivp) = 1.0 - n * beta;
    }
}

void 
SubdivisionSurface
::RecursiveAssignVertexLabel(
  MeshLevel *mesh, size_t t, size_t v, size_t id)
{
  if(mesh->triangles[t].vertices[v] != NOID)
    return;
  else
    mesh->triangles[t].vertices[v] = id;
  if(mesh->triangles[t].neighbors[(v+1) % 3] != NOID)
    RecursiveAssignVertexLabel(mesh, mesh->triangles[t].neighbors[(v+1) % 3],
                               (mesh->triangles[t].nedges[(v+1) % 3] + 1) % 3, id);
  if(mesh->triangles[t].neighbors[(v+2) % 3] != NOID)
    RecursiveAssignVertexLabel(mesh, mesh->triangles[t].neighbors[(v+2) % 3],
                               (mesh->triangles[t].nedges[(v+2) % 3] + 2) % 3, id);  
}

void SubdivisionSurface::Subdivide(const MeshLevel *parent, MeshLevel *child)
{
  // Get the numbers of triangles before and after
  size_t ntParent = parent->triangles.size();
  size_t ntChild = 4 * ntParent;
  size_t i, j, k;

  // Connect the child to the parent
  child->parent = parent;

  // Initialize the number of vertices in the new mesh
  child->nVertices = 0;

  // Subdivide each triangle into four
  child->triangles.resize(ntChild);
  for (i = 0; i < ntParent; i++)
    {
    // Get pointers to the four ctren
    const Triangle &pt = parent->triangles[i];
    Triangle *ct[4] = {
      &child->triangles[4*i], &child->triangles[4*i+1], 
      &child->triangles[4*i+2], &child->triangles[4*i+3]};

    // Set the neighbors within this triangle
    for (j = 0; j < 3; j++)
      {
      // Assign the neighborhoods within the pt triangle
      ct[j]->neighbors[j] = 4*i + 3;
      ct[3]->neighbors[j] = 4*i + j;
      ct[j]->nedges[j] = j;
      ct[3]->nedges[j] = j;

      // Assign neighborhoods outside the pt triangle
      if (pt.neighbors[(j+1) % 3] != NOID)
        {
        ct[j]->neighbors[(j+1) % 3] = 
        pt.neighbors[(j+1) % 3] * 4 + ((pt.nedges[(j+1) % 3] + 1) % 3);      
        ct[j]->nedges[(j+1) % 3] = pt.nedges[(j+1) % 3];
        }
      if (pt.neighbors[(j+2) % 3] != NOID)
        {
        ct[j]->neighbors[(j+2) % 3] = 
        pt.neighbors[(j+2) % 3] * 4 + ((pt.nedges[(j+2) % 3] + 2) % 3);
        ct[j]->nedges[(j+2) % 3] = pt.nedges[(j+2) % 3];
        }
      }
    }

  // Compute the number of edges in the pt graph
  size_t neParent = 0;
  for(i = 0; i < ntParent; i++) for(j = 0; j < 3; j++)
    neParent += (parent->triangles[i].neighbors[j] == NOID) ? 2 : 1;
  neParent >>= 1;

  // Assign vertex ids and weights. Under this scheme, the vertex ids
  // of the pt and ct should match (for even vertices) and odd
  // vertices appear at the end of the list
  
  // Create a mutable sparse matrix representation for the weights
  MutableSparseMatrix W(parent->nVertices + neParent, parent->nVertices);
  
  // Visit each of the even vertices, assigning it an id and weights
  for(i = 0; i < ntParent; i++) for(j = 0; j < 3; j++)
    if (child->triangles[4 * i + j].vertices[j] == NOID)
      {
      RecursiveAssignVertexLabel(child, 4*i+j, j, parent->triangles[i].vertices[j]);
      SetEvenVertexWeights(W, parent, child, 4*i+j, j);
      }
  
  // Visit each of the odd vertices, assigning it an id, and weights
  child->nVertices = parent->nVertices;
  for(i = 0; i < ntParent; i++) for(j = 0; j < 3; j++)
    if (child->triangles[4 * i + 3].vertices[j] == NOID)
      {
      RecursiveAssignVertexLabel(child, 4*i+3, j, child->nVertices++);
      SetOddVertexWeights(W, parent, child, 4*i+3, j);
      }

  // Copy the sparse matrix into immutable form
  child->weights.SetFromVNL(W);

  // If the parent's parent is not NULL, we need to multiply the sparse
  // matrices of the parent and child
  if(parent->parent)
    ImmutableSparseMatrix<double>::Multiply(
      child->weights, child->weights, parent->weights);

  // Compute the walks in the child mesh
  ComputeWalks(child);
}

void
SubdivisionSurface::
RecursiveSubdivide(const MeshLevel *src, MeshLevel *dst, size_t n)
{
  // N has to be greater than 1
  if(n == 0) 
    { *dst = *src; return; }
  else if(n == 1)
    { Subdivide(src, dst); return; }

  // Create an array of intermedial mesh levels
  MeshLevel *temp = new MeshLevel[n - 1];

  // Subdivide the intermediate levels
  const MeshLevel *parent = src; 
  for(size_t i = 0; i < n-1; i++)
    {
    Subdivide(parent, temp + i);
    parent = temp + i;
    }

  // Subdivide the last level
  Subdivide(parent, dst);

  // Delete the intermediates
  delete[] temp;

  // Set the parent pointer in dst to src (bypass intermediates)
  dst->parent = src;
}

void SubdivisionSurface::ExportLevelToVTK(MeshLevel &src, vtkPolyData *mesh)
{
  // The important thing with importing and exporting mesh levels is that the
  // order of triangles and vertices must be maintained.

  // Allocate an array of cells
  mesh->Allocate(src.triangles.size());

  // Add all of the triangles in the mesh
  for(size_t i = 0; i < src.triangles.size(); i++)
    {
    vtkIdType pts[3];
    pts[0] = src.triangles[i].vertices[0];
    pts[1] = src.triangles[i].vertices[1];
    pts[2] = src.triangles[i].vertices[2];
    mesh->InsertNextCell(VTK_TRIANGLE, 3, pts);
    }
}

void SubdivisionSurface::ImportLevelFromVTK(vtkPolyData *mesh, MeshLevel &dest)
{
  size_t i, j, k;

  // Prepare the mesh
  mesh->BuildCells();

  // Typedefs
  typedef pair<size_t, size_t> HalfEdge;
  typedef pair<size_t, short> TriangleRep;
  typedef map<HalfEdge, TriangleRep> TriangleMap;

  // Get the number of triangles in the mesh
  size_t nTriangles = mesh->GetNumberOfCells();
  dest.triangles.resize(nTriangles);

  // Set the number of vertices in the mesh
  dest.nVertices = mesh->GetNumberOfPoints();

  // Initialize the triangle map
  TriangleMap tmap;

  // For each triangle, compute the neighbor. This can be done by first enumerating
  // all the edges in the mesh. For each edge there will be one or two triangles
  for(i = 0; i < nTriangles; i++)
    {
    // Get the points from the current triangle
    vtkIdType npts, *pts;
    mesh->GetCellPoints(i, npts, pts);

    // If the number of points is not 3, return with an exception
    if(npts != 3) throw string("Mesh contains cells other than triangles");

    // Associate each half-edge with a triangle
    for(j = 0; j < 3; j++)
      {
      // Set the vertices in each triangle
      dest.triangles[i].vertices[j] = pts[j];

      // Create the key and value
      HalfEdge he(pts[(j+1) % 3], pts[(j+2) % 3]);
      TriangleRep trep(i, (short) j);

      // Insert the half-edge and check for uniqueness
      pair<TriangleMap::iterator, bool> rc = tmap.insert(make_pair(he, trep));
      if(!rc.second)
        throw string("Half-edge appears twice in the mesh, that is illegal!");
      }
    }

  // Take a pass through all the half-edges. For each, set the corresponding triangle's
  // neighbor and neighbor index
  for(TriangleMap::iterator it = tmap.begin(); it != tmap.end(); ++it)
    {
    TriangleRep &trep = it->second;
    HalfEdge opposite(it->first.second, it->first.first);
    TriangleMap::iterator itopp = tmap.find(opposite);
    if(itopp != tmap.end())
      {
      dest.triangles[trep.first].neighbors[trep.second] = itopp->second.first;
      dest.triangles[trep.first].nedges[trep.second] = itopp->second.second;
      }
    }

  // Set the mesh's parent to NULL
  dest.parent = NULL;

  // Finally, compute the walks in this mesh
  ComputeWalks(&dest);
}

void SubdivisionSurface
::ApplySubdivision(vtkPolyData *src, vtkPolyData *target, MeshLevel &m)
{
  size_t i, j, k;

  // Initialize the target mesh
  vtkPoints *tpoints = vtkPoints::New();
  vtkCellArray *tcells = vtkCellArray::New();
  
  // Add the points to the output mesh
  for(i = 0; i < m.nVertices; i++)
    {
    // Output point
    vnl_vector_fixed<vtkFloatingPointType, 3> p(0.0);
    
    // Iterate over columns
    typedef ImmutableSparseMatrix<double>::RowIterator IteratorType;
    for(IteratorType it = m.weights.Row(i); !it.IsAtEnd(); ++it)
      {
      vnl_vector_fixed<vtkFloatingPointType, 3> x(src->GetPoints()->GetPoint(it.Column()));
      p += x * (vtkFloatingPointType) it.Value();
      }
    
    // Add the point to the output mesh
    tpoints->InsertNextPoint(p.data_block());
    }

  // Add the cells to the mesh
  for(i = 0; i < m.triangles.size(); i++)
    {
    Triangle &t = m.triangles[i];
    vtkIdType ids[3];
    ids[0] = t.vertices[0];
    ids[1] = t.vertices[1];
    ids[2] = t.vertices[2];
    tcells->InsertNextCell(3, ids);
    }

  // Set the points in the mesh
  target->SetPoints(tpoints);
  target->SetPolys(tcells);
  tpoints->Delete();
  tcells->Delete();
}

void SubdivisionSurface
::ApplySubdivision(SMLVec3d *xSrc, double *rhoSrc, SMLVec3d *xTrg, double *rhoTrg, MeshLevel &m)
{
  size_t i, j, k;

  // Add the points to the output mesh
  for(i = 0; i < m.nVertices; i++)
    {
    // Output point
    xTrg[i].fill(0.0); rhoTrg[i] = 0;
    
    // Iterate over columns
    typedef ImmutableSparseMatrix<double>::RowIterator IteratorType;
    for(IteratorType it = m.weights.Row(i); !it.IsAtEnd(); ++it)
      {
      xTrg[i] += it.Value() * xSrc[it.Column()];
      rhoTrg[i] += it.Value() * rhoSrc[it.Column()];
      }
    }
}

void SubdivisionSurface
::ApplySubdivision(const double *xsrc, double *xdst, size_t nComp, MeshLevel &m)
{  
  typedef ImmutableSparseMatrix<double>::RowIterator IteratorType;

  // Get the total number of values
  size_t n = nComp * m.nVertices;

  // Set the target array to zero
  std::fill_n(xdst, n, 0.0);

  // Iterate over the output vertices
  for(size_t iv = 0; iv < m.nVertices; iv++)
    {
    size_t i = iv * nComp;

    // Iterate over contributing input vertices
    for(IteratorType it = m.weights.Row(iv); !it.IsAtEnd(); ++it)
      {
      double w = it.Value();
      size_t j = it.Column() * nComp;

      // Iterate over the component
      for(size_t k = 0; k < nComp; k++)
        xdst[i+k] = w * xsrc[j+k];
      }
    }
}
  
bool SubdivisionSurface::CheckMeshLevel (MeshLevel &mesh)
{
  // Check the following rules for all triangles in the mesh
  // 1. T[T[i].nbr[j]].nbr[T[i].ne[j]] == i for all i, j    (i am my neighbors neighbor)
  // 2. T[T[i].nbr[j]].ne[T[i].ne[j]] == j for all i, j    (i am my neighbors neighbor)
  // 3. T[i].v[j] == T[T[i].nbr[j+k]].v[T[i].ne[j+k]+k] for all i, j, k=(1, 2), 
  // with modulo 3 addition
  size_t nerr = 0;
  for(size_t i = 0; i < mesh.triangles.size(); i++)
    {
    Triangle &t = mesh.triangles[i];
    for(size_t j = 0; j < 3; j++)
      {
      if(t.neighbors[j] != NOID) 
        {
        Triangle &tn = mesh.triangles[t.neighbors[j]];
        if(tn.neighbors[t.nedges[j]] != i)
          {
          cout << "Error " << nerr++ << 
            " Rule 1 violated for i = " << i << " and j = " << j << endl;
          }
        if(tn.nedges[t.nedges[j]] != j)
          {
          cout << "Error " << nerr++ << 
            " Rule 2 violated for i = " << i << " and j = " << j << endl;
          }
        }
      for(size_t k = 1; k < 3; k++)
        {
        if(t.neighbors[(j+k) % 3] != NOID) 
          {
          Triangle &tk = mesh.triangles[t.neighbors[(j+k) % 3]];
          if(t.vertices[j] != tk.vertices[(t.nedges[(j+k) % 3] + k ) % 3])
            {
            cout << "Error " << nerr++ << 
              " Rule 3 violated for i = " << i << ", j = " << j << " and k = " << k << endl;
            }
          }
        }
      }
    }

  return (nerr == 0);
}

// We need to instantiate sparse matrix of NeighborInfo objects
#include "SparseMatrix.txx"
template class ImmutableSparseArray<SubdivisionSurface::NeighborInfo>;

