#ifndef __EventListenerModel_h_
#define __EventListenerModel_h_

#include <set>

/** Event object, parent to all events */
template <class TSource>
class EventObject
{
public:
  /** Create a new abstract event */
  EventObject(TSource *source)
    { m_Source = source; }

  /** Virtual destructor */
  virtual ~EventObject() {}

  /** Get the source of the event */
  TSource *GetSource()
    { return m_Source; }

protected:
  /** Event source */
  TSource *m_Source;
};

/** Macros for setting up an event source */
#define pyAddListenerMacro(name,interface) \
virtual void Add##name(interface *listener) \
{ m_##name.insert(listener); } 

#define pyRemoveListenerMacro(name,interface) \
virtual void Remove##name(interface *listener) \
{ m_##name.erase(listener); }

#define pyBroadcastEventMacro(name, interface, function, event) \
virtual void Broadcast##function(event *evt) \
{ \
  std::set<interface *>::iterator it = m_##name.begin(); \
  while(it != m_##name.end()) \
    (*it++)->function(evt); \
}

#define pySetWithNofityMacro(name, type, function, event) \
virtual void Set##name(const type &in##name) \
{ \
  if(in##name != m_##name) \
    { \
    event evt(this); \
    Broadcast##function(&evt); \
    m_##name = in##name; \
    } \
}
    
#define pyListenerArrayMacro(name,interface) \
std::set<interface *> m_##name;

#endif // __EventListenerModel_h_
