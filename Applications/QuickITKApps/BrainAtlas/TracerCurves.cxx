#include "TracerCurves.h"
#include <sstream>
#include <cassert>

TracerCurves::IdType
TracerCurves
::AddControlPoint(IdType id, MeshVertex v)
{
  // Insert the control point with the new id
  m_Controls[id] = ControlPoint(v);

  // Return the id
  return id;
}

TracerCurves::IdType
TracerCurves
::AddLink(IdType id, IdType iStartCtl, IdType iEndCtl, IdType iCurve, list<MeshVertex> path)
{
  // Make sure correct parameters are passed in
  assert(m_Curves.find(iCurve) != m_Curves.end());
  assert(m_Controls.find(iStartCtl) != m_Controls.end());
  assert(m_Controls.find(iEndCtl) != m_Controls.end());

  // Make sure that the curve actually ends with the control being added
  assert(m_Curves[iCurve].links.size() == 0 ||
    m_Links[m_Curves[iCurve].links.back()].ctlEnd == iStartCtl);

  // Insert the control point with the new id
  m_Links[id] = Link(iStartCtl,iEndCtl,iCurve);
  m_Links[id].points = path;

  // Associate the control points with the links
  m_Controls[iStartCtl].links.insert(id);
  m_Controls[iEndCtl].links.insert(id);

  // Add the link to the curve
  m_Curves[iCurve].links.push_back(id);

  // Add the control points and mesh points to the curve
  
  // If the curve is empty, insert the vertex corresponding to the 
  // first control point
  if(m_Curves[iCurve].points.size() == 0)
    {
    m_Curves[iCurve].points.push_back(m_Controls[iStartCtl].vertex);
    m_Curves[iCurve].controls.push_back(iStartCtl);
    }

  // Insert the intermediate mesh points
  m_Curves[iCurve].points.insert(
    m_Curves[iCurve].points.end(),path.begin(),path.end());
    
  // Insert the end control point
  m_Curves[iCurve].points.push_back(m_Controls[iEndCtl].vertex);
  m_Curves[iCurve].controls.push_back(iEndCtl);

  // Return the id
  return id;
}

TracerCurves::IdType
TracerCurves
::AddCurve(IdType id, const char *name)
{
  // Generate the id for the control point
  m_Curves[id] = Curve(name);

  // Return the id
  return id;
}

void 
TracerCurves
::SaveAsText(ostream &fout)
{
  // Write the header
  fout << "# GeeLab BrainTracer Curve File " << endl;
  fout << "# Version 1.0 " << endl;
  fout << "# Please don't edit the above two lines!" << endl;
  fout << "# " << endl;
  fout << "# Format Specification: " << endl;
  fout << "#    CONTROL id vertex-id " << endl;
  fout << "#    CURVE   id name" << endl;
  fout << "#    LINK    id start-ctl end-ctl curve num-vertices v1 .. vn" << endl;
  fout << "#" << endl;

  // Write the contols
  ControlMap::const_iterator itControl = m_Controls.begin();
  while(itControl != m_Controls.end())
    {
    const ControlPoint &cp = itControl->second;
    fout << "CONTROL\t" << itControl->first << "\t" << cp.vertex << endl;
    itControl++;
    }

  // Write the curves
  CurveMap::const_iterator itCurve = m_Curves.begin();
  while(itCurve != m_Curves.end())
    {
    const Curve &c = itCurve->second;
    fout << "CURVE\t" << itCurve->first << "\t\"" << c.name << "\"" << endl;
    itCurve++;
    }

  // Write the links
  LinkMap::const_iterator itLink = m_Links.begin();
  while(itLink != m_Links.end())
    {
    const Link &l = itLink->second;
    fout << "LINK\t" << itLink->first << "\t" << l.ctlStart << "\t"
      << l.ctlEnd << "\t" << l.curve << "\t" << l.points.size();
   
    // Write out the connecting points
    MeshCurve::const_iterator it = l.points.begin();
    while(it != l.points.end())
      {
      fout << " " << *it; 
      it++;
      }

    fout << endl;
    itLink++;
    }
}

bool
TracerCurves
::LoadAsText(istream &fin, bool clear)
{
  // Create a copy of the data or blank data
  map<IdType, Curve> oldCurves = m_Curves;
  map<IdType, ControlPoint> oldControls = m_Controls;
  map<IdType, Link> oldLinks = m_Links;

  // Clear the data if necessary
  if(!clear)
    {
    m_Controls.clear();
    m_Links.clear();
    m_Curves.clear();
    }
  
  try 
    {
    // String holding the current line in the file, and the starting key word
    string line, key;

    // Read each line of the file separately
    for(unsigned int iLine=0;!fin.eof();iLine++)
      {
      // Read the line into a string
      getline(fin,line);

      // Check if the line is a comment or a blank line
      if(line[0] == '#' || line.length() == 0)
        continue;

      // Create a stream to parse that string
      istringstream iss(line), tss(line);
      iss.exceptions(ios_base::badbit | ios_base::failbit);

      // We need a keyword before proceeding
      iss >> key;
      if(key == "CONTROL")
        {
        // This is a control point
        IdType id; iss >> id;
        MeshVertex vtx; iss >> vtx;

        // Add the control point
        AddControlPoint(id, vtx);
        }
      else if(key == "CURVE")
        {
        // This is a curve defition
        IdType id; iss >> id;
        string name; iss >> name;

        // Add the curve
        AddCurve(id, name.c_str());
        }
      else if(key == "LINK")
        {
        // Read the link
        IdType id, iStart, iEnd, iCurve, iLength;
        iss >> id; iss >> iStart; iss >> iEnd; iss >> iCurve; iss >> iLength;
        
        // Create a list to store the path
        MeshCurve path;
        for(unsigned int i = 0;i < iLength;i++)
          {
          MeshVertex v; iss >> v;
          path.push_back(v);
          }

        // Add the link
        AddLink(id,iStart,iEnd,iCurve,path);
        }
      else
        continue;
      }

    return true;
    }
  catch(...)
    {
    // Restore the data
    m_Curves = oldCurves;
    m_Controls = oldControls;
    m_Links = oldLinks;
    return false;
    }
}

