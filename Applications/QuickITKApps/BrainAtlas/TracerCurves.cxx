#include "TracerCurves.h"
#include <sstream>
#include <cassert>

string StringEncode(const string &source)
{
  ostringstream oss;
  for(unsigned int i=0;i<source.size();i++)
    {
    char c = source[i];
    if((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
      (c >= '0' && c <= '9') || (c == '_'))
      oss << c;
    else
      oss << '\\' << std::setfill('0') << std::oct << 
        std::setw(3) << ((int)c);
    }
  return oss.str();
}

string StringDecode(const string &source)
{
  ostringstream oss;
  for(unsigned int i=0;i<source.size();i++)
    {
    char c = source[i];
    if(c != '\\')
      oss << c;
    else if(i+3 < source.size())
      oss << (char)(
        64 * (int)(source[++i] - '0') + 
         8 * (int)(source[++i] - '0') + 
         1 * (int)(source[++i] - '0'));
    }
  return oss.str();
}

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
::RebuildCurvePoints(IdType iCurve)
{
  // Go through and rebuild the points in the given curve
  Curve &xCurve = m_Curves[iCurve];

  // Clear the list of points
  xCurve.points.clear();

  // Create iterators for links and control points
  IdList::const_iterator itCtl = xCurve.controls.begin();
  IdList::const_iterator itLink = xCurve.links.begin();

  // Alternate between control points and links
  while(itCtl != xCurve.controls.end())
    {
    // Add the current control to the list
    xCurve.points.push_back( m_Controls[*itCtl].vertex );

    // If there is no link, stop
    if(itLink == xCurve.links.end()) break;

    // Add the contents of the current link to the control
    MeshCurve &xLinkPoints = m_Links[*itLink].points;
    xCurve.points.insert(xCurve.points.end(), 
      xLinkPoints.begin(), xLinkPoints.end());

    // Increment the iterators
    ++itLink; ++itCtl;
    }
}

void 
TracerCurves
::DeleteCurve(IdType iCurve)
{
  // Delete the curve
  m_Curves.erase(iCurve);
  
  // Clean up links and control points
  CleanUpDeadLinksAndControls();
}

void 
TracerCurves
::DeleteLastControlPointInCurve(IdType iCurve)
{
  // Check input validity
  assert(m_Curves.find(iCurve) != m_Curves.end());
  assert(m_Curves[iCurve].controls.size() > 0);
 
  // Delete the last control point
  m_Curves[iCurve].controls.pop_back();

  // Delete the last link if there is one
  if(m_Curves[iCurve].links.size())
    m_Curves[iCurve].links.pop_back();

  // Clean up dead controls and links
  CleanUpDeadLinksAndControls();

  // Rebuild the points in the curve
  RebuildCurvePoints(iCurve);
}

void 
TracerCurves
::CleanUpDeadLinksAndControls()
{
  // Create tables in order to remove links and control points that 
  // belong to zero curves
  typedef set< IdType > CountSet;
  CountSet xLiveLinks, xLiveControls;
  
  // Populate the tables
  CurveMap::const_iterator itCurve = m_Curves.begin();
  while(itCurve != m_Curves.end())
    {
    // Get a curve reference
    const Curve &xCurve = itCurve->second;
    
    // Add the links in that curve to list of live links
    xLiveLinks.insert( xCurve.links.begin(), xCurve.links.end() );
    
    // Add the controls in that curve to list of live controls
    xLiveControls.insert( xCurve.controls.begin(), xCurve.controls.end() );
    
    // On to the next curve
    ++itCurve;
    }
  
  // Delete all links that don't belong to any curve
  LinkMap::iterator itLink = m_Links.begin();
  while(itLink != m_Links.end())
    {
    // Check if link is live
    if(xLiveLinks.find(itLink->first) == xLiveLinks.end())
      {
      // Remove references to this link in control points
      m_Controls[itLink->second.ctlStart].links.erase(itLink->first);
      m_Controls[itLink->second.ctlEnd].links.erase(itLink->first);

      // Erase the link from the link map
      m_Links.erase(itLink);
      }
    
    // Go to the next element  
    ++itLink;
    }

  // Delete all the controls that don't belong to any curve
  ControlMap::iterator itControl = m_Controls.begin();
  while(itControl != m_Controls.end())
    {
    // Check if control is live
    if(xLiveControls.find(itControl->first) == xLiveControls.end())
      {
      m_Controls.erase(itControl);
      }
    
    // Go to the next element
    ++itControl;
    }
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
    fout << "CURVE\t" << itCurve->first << " " << 
      StringEncode(c.name) << endl;
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
        AddCurve(id, StringDecode(name).c_str());
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

