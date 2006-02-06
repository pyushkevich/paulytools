#include <vector>
#include "Registry.h"

class vtkPolyData;

using namespace std;

class SubdivisionSurface 
{
public:

  struct Triangle
  {
    size_t vertices[3];
    size_t neighbors[3];
    short nedges[3];

    // Initializes to dummy values
    Triangle();
  };

  struct MeshLevel
  {
    typedef vector<Triangle>::iterator TriangleIt;
    vector<Triangle> triangles;
    size_t nVertices;
  };

  /** Subdivide a mesh level once */
  void Subdivide(MeshLevel &src, MeshLevel &dst);

  /** Import a mesh from a VTK mesh that's been wrapped in a half-edge structure */
  void ImportLevelFromVTK(vtkPolyData *, MeshLevel &dest);

  /** Load a mesh level from a registry */
  void LoadMeshLevel(Registry &registry);

  /** Test the correctness of a mesh level */
  bool CheckMeshLevel(MeshLevel &mesh);


private:

  // Vertex assignment function (visit vertices to assign labels)
  void RecursiveAssignVertexLabel(MeshLevel &mesh, size_t t, size_t v, size_t id);

  const static size_t NOID;
};

