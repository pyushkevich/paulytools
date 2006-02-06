#include "SubdivisionSurface.h"
#include "vtkPolyData.h"
#include <string>
#include <map>

using namespace std;

// Not an ID
const size_t SubdivisionSurface::NOID = 0xFFFFFFFF;

SubdivisionSurface::Triangle::Triangle()
{
  for(size_t i = 0; i < 3; i++)
    {
    vertices[i] = NOID;
    neighbors[i] = NOID;
    nedges[i] = NOID;
    }
}

void SubdivisionSurface::RecursiveAssignVertexLabel(MeshLevel &mesh, size_t t, size_t v, size_t id)
{
  if(mesh.triangles[t].vertices[v] != NOID)
    return;
  else
    mesh.triangles[t].vertices[v] = id;
  if(mesh.triangles[t].neighbors[(v+1) % 3] != NOID)
    RecursiveAssignVertexLabel(mesh, mesh.triangles[t].neighbors[(v+1) % 3],
                               (mesh.triangles[t].nedges[(v+1) % 3] + 1) % 3, id);
  if(mesh.triangles[t].neighbors[(v+2) % 3] != NOID)
    RecursiveAssignVertexLabel(mesh, mesh.triangles[t].neighbors[(v+2) % 3],
                               (mesh.triangles[t].nedges[(v+2) % 3] + 2) % 3, id);  
}

void SubdivisionSurface::Subdivide(MeshLevel &src, MeshLevel &dst)
{
  // Get the numbers of triangles before and after
  size_t ntParent = src.triangles.size();
  size_t ntChild = 4 * ntParent;
  size_t i, j, k;

  // Initialize the number of vertices in the new mesh
  dst.nVertices = 0;

  // Subdivide each triangle into four
  dst.triangles.resize(ntChild);
  for (i = 0; i < ntParent; i++)
    {
    // Get pointers to the four children
    Triangle &parent = src.triangles[i];
    Triangle *child[4] = {
      &dst.triangles[4*i], &dst.triangles[4*i+1], 
      &dst.triangles[4*i+2], &dst.triangles[4*i+3]};

    // Set the neighbors within this triangle
    for (j = 0; j < 3; j++)
      {
      // Assign the neighborhoods within the parent triangle
      child[j]->neighbors[j] = 4*i + 3;
      child[3]->neighbors[j] = 4*i + j;
      child[j]->nedges[j] = j;
      child[3]->nedges[j] = j;

      // Assign neighborhoods outside the parent triangle
      if (parent.neighbors[(j+1) % 3] != NOID)
        {
        child[j]->neighbors[(j+1) % 3] = 
        parent.neighbors[(j+1) % 3] * 4 + ((parent.nedges[(j+1) % 3] + 1) % 3);      
        child[j]->nedges[(j+1) % 3] = parent.nedges[(j+1) % 3];
        }
      if (parent.neighbors[(j+2) % 3] != NOID)
        {
        child[j]->neighbors[(j+2) % 3] = 
        parent.neighbors[(j+2) % 3] * 4 + ((parent.nedges[(j+2) % 3] + 2) % 3);
        child[j]->nedges[(j+2) % 3] = parent.nedges[(j+2) % 3];
        }
      }
    }

  // Now that the triangle neighborhoods are in place, we need to assign vertices to each of 
  // the new triangles. Since each vertex can be shared by multiple triangles, we must be
  // quite careful in doing this.
  for (i = 0; i < ntChild; i++) for (j = 0; j < 3; j++)
    if (dst.triangles[i].vertices[j] == NOID)
      RecursiveAssignVertexLabel(dst, i, j, dst.nVertices++);

  // Now, we have a second-level mesh that should be valid. We should test the validity separately
}

struct ImportEdge
{
  size_t v1, v2;
  size_t t1, t2;
};

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
