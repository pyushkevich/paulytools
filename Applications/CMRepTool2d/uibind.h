#pragma warning( disable : 4786 )  

#ifndef __UIBIND_H__
#define __UIBIND_H__

#include <map>
#include <string>
#include <Fl/Fl_Widget.H>
#include <Fl/Fl_Button.H>
#include <Fl/Fl_Menu_.H>
#include <Fl/Fl_Valuator.H>
#include <Fl/Fl_Group.H>
#include <Fl/Fl_Light_Button.H>

using namespace std;

class UIBinding {
public:
	string key;

	virtual void setWidget(string s) = 0;
	virtual string readWidget() = 0;
	virtual ~UIBinding() {};
};

// A way to bind widgets to settings...
class UIBinder {
private:
	map <string,UIBinding *> keyMap;
	map <Fl_Widget *,UIBinding *> widMap;
public:
	void bindRadios(string key,Fl_Group *container,Fl_Button *deflt = NULL);
	void bindChoice(string key,Fl_Menu_ *menu,int deflt = -1);
	void bindColor(string key,Fl_Light_Button *widget);
	void bind(string key,Fl_Valuator *widget);
	void bind(string key,Fl_Button *widget);

	// This method is called whenever a widget calls back
	bool readWidget(Fl_Widget *w,string &key,string &value);
	void setWidget(string key,string value);
};

extern UIBinder uiBind;

#endif //__UIBIND_H__
