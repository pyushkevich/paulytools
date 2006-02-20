#ifndef __CartesianGridMedialIterationContext_h_
#define __CartesianGridMedialIterationContext_h_

class CartesianGridMedialIterationContext
{
public:
  typedef SubdivisionSurface::MeshLevel MeshLevel;
  
  SubdivisionSurfaceMedialIterationContext(MeshLevel &inMesh)
    : mesh(inMesh), n(inMesh.nVertices)
    {
    // Generate the mapping from atom index / boundary side into 
    // boundary index
    xBndMap = new size_t[n + 1];
    xBndMap[0] = 0;
    for(size_t i = 0; i < n; i++)
      xBndMap[i+1] = xBndMap[i] + mesh.IsVertexInternal(i) ? 2 : 1;    
    };

  ~SubdivisionSurfaceMedialIterationContext()
    { delete xBndMap; }
    
  bool IsEdgeAtom(size_t i)
    { return xBndMap[i+1] == 1 + xBndMap[i]; }
  
  size_t GetNumberOfAtoms()
    { return n; }
  
  size_t GetAtomIndexInTriangle(size_t iTri, size_t jVert)
    { return mesh.triangles[iTri].vertices[jVert]; }

  size_t GetBoundaryPointIndex(size_t iAtom, size_t iSide)
    { return xBndMap[iAtom + iSide] - iSide; }

  size_t GetInternalPointIndex(size_t iAtom, size_t iSide, size_t iDepth)
    { 
    return iDepth == 0 ? iAtom : 
      n + (iDepth - 1) * xBndMap[n] + GetBoundaryPointIndex(iAtom, iSide);
    }
      
  size_t GetProfileIntervalIndex(size_t iBoundary, size_t iDepth)
    { return iDepth * xBndMap[n] + iBoundary; }

private:
  MeshLevel &mesh;

  // This is an array of form 0 2 4 5 6 8 ... 
  // difference of 2 between x[i+1] and x[i] indicates i is internal.
  // Then the boundary index for atom i, side j (0,1) is x[i+j]-j
  size_t *xBndMap, n;
};


#endif
