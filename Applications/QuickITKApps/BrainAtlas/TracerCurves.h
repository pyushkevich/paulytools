/** 
 * TracerCurves
 * A class representing the curve network used by the BrainTracer
 */
#ifndef _TracerCurves_h_
#define _TracerCurves_h_

#include <list>
#include <map>
#include <iostream>
#include <fstream>
#include <set>

#include <vtkSystemIncludes.h>
#include <vnl/vnl_vector_fixed.h>

using namespace std;

/**
 * A control point is a point selected by the user
 * A span is an interval between two control points that is traced
 * on the surface of the mesh
 * A curve is a list of connected soanbs:
 */

class TracerCurves {
public:

  /** Pointer to a mesh vertex */
  typedef vtkIdType MeshVertex;
  typedef list<MeshVertex> MeshCurve;
  typedef unsigned int IdType;
  typedef list<IdType> IdList;

  /** Representation of a 3D vector */
  typedef vnl_vector_fixed<double, 3> Vec;

  /** Used to report curve names in sorted order */
  typedef pair<string, IdType> StringIdPair;
  
  /** Structure representing a control point */
  struct ControlPoint 
    {
    /** The index of the control point as a vertex in the mesh */
    MeshVertex iVertex;

    /** The actual position of the vertex in space */
    Vec xVertex;

    /** Set of links to which the control point belongs */
    set<IdType> links;
    
    ControlPoint() {}
    ControlPoint(MeshVertex inId, const Vec &inX) 
      : iVertex(inId), xVertex(inX) {}
  };

  /** Structure representing a span between two control points */
  struct Link {
    IdType ctlStart, ctlEnd;
    MeshCurve points;
    IdType curve;
    
    Link() : ctlStart(0), ctlEnd(0), curve(0) {}
    Link(IdType cpStart, IdType cpEnd, IdType inCurve)
      : ctlStart(cpStart), ctlEnd(cpEnd), curve(inCurve) {}
  };

  /** Structure representing a curve */
  struct Curve {
    string name;
    list<IdType> links, controls;
    MeshCurve points;

    Curve() {}
    Curve(const char *inName) 
      : name(inName) {}
  };

  /** Add a control point to the list */
  IdType AddControlPoint(MeshVertex id, const Vec &x)
    { return AddControlPoint( GenerateId(m_Controls), id, x); }

  /** Add a link between two control points, appending it to a given curve */
  IdType AddLink( IdType iStartCtl, IdType iEndCtl, IdType iCurve, MeshCurve &path )
    { return AddLink( GenerateId(m_Links), iStartCtl, iEndCtl, iCurve, path); }
  
  /** Add an empty curve to the collection */
  IdType AddCurve( const char * name )
    { return AddCurve( GenerateId(m_Curves), name); }

  /** Get a chain of mesh vertices associated with a curve */
  const MeshCurve & GetCurveVertices(IdType iCurve) const
    { return m_Curves.find(iCurve)->second.points; }

  /** Get a chain of control points associated with a curve */
  const list<IdType> & GetCurveControls(IdType iCurve) const
    { return m_Curves.find(iCurve)->second.controls; }

  /** Check if a curve exists */
  bool IsCurvePresent(IdType id) 
    { return m_Curves.find(id) != m_Curves.end(); }

  /** Get the number of curves */
  unsigned int GetNumberOfCurves() const
    { return m_Curves.size(); }

  /** Get the curve name */
  const char * GetCurveName(IdType iCurve) const
    { return m_Curves.find(iCurve)->second.name.c_str(); }

  /** Change the name of a curve */
  void SetCurveName(IdType iCurve, const char *name)
    { m_Curves[iCurve].name = name; }

  /** A sorting predicate for pairs */
  struct StringIdPred : public binary_function<StringIdPair,StringIdPair,bool>
    {
    bool operator()(const StringIdPair &p1, const StringIdPair &p2)
      { return p1.second > p2.second; }
    };
    
  /** Get the list of curves sorted by name */
  void GetAlphabeticCurveList(list<StringIdPair> &target) const
    {
    CurveMap::const_iterator it;
    
    target.clear();
    for(it = m_Curves.begin(); it!=m_Curves.end(); it++)
      target.push_back(StringIdPair(it->second.name, it->first));
    target.sort();
    }

  /** Get a list of curves in any order */
  void GetCurveIdList(list<IdType> &target) const
    {
    CurveMap::const_iterator it;
    target.clear();
    for(it = m_Curves.begin(); it!=m_Curves.end(); it++)
      target.push_back(it->first);
    }

  /** Delete a curve */
  void DeleteCurve(IdType iCurve);

  /** Delete the lasp point on a curve */
  void DeleteLastControlPointInCurve(IdType iCurve);
  
  /** Save the data structure to a file */
  void SaveAsText(std::ostream &sout);

  /** Read the structure from a file (optionally append existing curves) */
  bool LoadAsText(std::istream &sin, bool clear = true);

  /** Clean the data by removing all curves, points, and links */
  void RemoveAllData()
    {
    m_Curves.clear();
    m_Controls.clear();
    m_Links.clear();
    }

  /** Get the vertex id at a given control point */
  MeshVertex GetControlPointVertexId(IdType iControl) const
    { return m_Controls.find(iControl)->second.iVertex; }

  /** Get the vertex position at a given control point */
  const Vec &GetControlPointPosition(IdType iControl) const
    { return m_Controls.find(iControl)->second.xVertex; }

  /** Get the vertex position at a given control point */
  void SetControlPointPosition(IdType iControl, const Vec &x)
    { m_Controls.find(iControl)->second.xVertex = x; }

private:

  /** A list of curves */
  typedef map<IdType, Curve> CurveMap;
  CurveMap m_Curves;
  
  /** A list of control points */
  typedef map<IdType, ControlPoint> ControlMap;
  ControlMap m_Controls;

  /** A list of links */
  typedef map<IdType, Link> LinkMap;
  LinkMap m_Links;

  /** An ID generator for map */
  template <class T> static IdType GenerateId(const map<IdType, T> &xMap)
    { 
    IdType id = xMap.size() + 1;
    while(xMap.find(id) != xMap.end())
      id = (rand() << 16) + rand();
    return id; 
    }
  
  /** Add a control point to the list */
  IdType AddControlPoint(IdType id, MeshVertex v, const Vec &x);

  /** Add a link between two control points, appending it to a given curve */
  IdType AddLink(IdType id, IdType iStartCtl, IdType iEndCtl, 
    IdType iCurve, list<MeshVertex> path );
  
  /** Add an empty curve to the collection */
  IdType AddCurve(IdType id, const char * name );

  /** Clean up links and controls that do not belong to any curves */
  void CleanUpDeadLinksAndControls();

  /** Rebuild the point list in a curve */
  void RebuildCurvePoints(IdType iCurve);
};


#endif // _TracerCurves_h_
