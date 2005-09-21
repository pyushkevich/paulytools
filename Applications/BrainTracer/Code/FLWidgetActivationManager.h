#ifndef _FLWidgetActivationManager_h_
#define _FLWidgetActivationManager_h_

#include "FL/Fl_Widget.H"
#include "FL/Fl_Menu_Item.H"
#include <map>
#include <set>

using namespace std;

template<typename TFlag>
class FLWidgetActivationManager {
public:
  // FLWidgetActivationManager();
  // ~FLWidgetActivationManager();

  void AddWidget(Fl_Widget *widget, TFlag flag)
    {
    // FlagData &fdata = GetFlagData(flag);
    m_Flags[flag].Widgets.insert(widget);
    widget->deactivate();
    }

  void AddMenuItem(Fl_Menu_Item *menu, TFlag flag)
    {
    // FlagData &fdata = GetFlagData(flag);
    m_Flags[flag].MenuItems.insert(menu);
    menu->deactivate();
    }

  /** Flag A implies flag B */
  void SetFlagImplies(TFlag flagA, TFlag flagB)
    {
    //FlagData &fda = GetFlagData(flagA);
    //FlagData &fdb = GetFlagData(flagB);
    
    m_Flags[flagA].Implies.insert(flagB);
    m_Flags[flagB].ImpliedBy.insert(flagA);
    }

  void UpdateFlag(TFlag flag, bool on)
    {
    // Update the flag
    m_Flags[flag].State = on;

    if(on)
      {
      // Activate all the widgets
      typename set<Fl_Widget *>::iterator itWidget = m_Flags[flag].Widgets.begin();
      while(itWidget != m_Flags[flag].Widgets.end())
        (*itWidget++)->activate();

      // Activate all the widgets
      typename set<Fl_Menu_Item *>::iterator itMenu = m_Flags[flag].MenuItems.begin();
      while(itMenu != m_Flags[flag].MenuItems.end())
        (*itMenu++)->activate();

      // Impty all the other flags
      typename set<TFlag>::iterator it = m_Flags[flag].Implies.begin();
      while(it != m_Flags[flag].Implies.end())
        UpdateFlag(*(it++), on);
      }
    else
      {
      // Deactivate all the widgets
      typename set<Fl_Widget *>::iterator itWidget = m_Flags[flag].Widgets.begin();
      while(itWidget != m_Flags[flag].Widgets.end())
        (*itWidget++)->deactivate();

      // Activate all the widgets
      typename set<Fl_Menu_Item *>::iterator itMenu = m_Flags[flag].MenuItems.begin();
      while(itMenu != m_Flags[flag].MenuItems.end())
        (*itMenu++)->deactivate();

      // Unset all the flags implied by this one
      typename set<TFlag>::iterator it = m_Flags[flag].ImpliedBy.begin();
      while(it != m_Flags[flag].ImpliedBy.end())
        UpdateFlag(*(it++), on);
      }
    }

private:
  
  struct FlagData
    {
    set<TFlag> Implies;
    set<TFlag> ImpliedBy;
    set<Fl_Widget *> Widgets;
    set<Fl_Menu_Item *> MenuItems;
    bool State;

    FlagData() { State = false; }
    };

  typedef map<TFlag, FlagData> FlagMap;
  typedef typename FlagMap::iterator FlagIterator;

  FlagMap m_Flags;
  
  FlagData &GetFlagData(TFlag flag)
    {
    FlagIterator it = m_Flags.find(flag);
    if(it == m_Flags.end())
      {
      m_Flags[flag] = FlagData();
      return m_Flags[flag];
      }
    else return it->second;
    }

};

#endif
