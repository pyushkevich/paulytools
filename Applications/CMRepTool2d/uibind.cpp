#pragma warning( disable : 4786 )  

#include <list>
#include <string>

#include <Fl/Fl_Widget.H>
#include <Fl/Fl_Button.H>
#include <Fl/Fl_Light_Button.H>
#include <Fl/Fl_Menu.H>
#include <Fl/Fl_Menu_.H>
#include <Fl/Fl_Valuator.H>
#include <Fl/Fl_Group.H>
#include <Fl/Fl.H>

#include "uibind.h"

using namespace std;

class ValuatorBinding : public UIBinding {
	Fl_Valuator *widget;
public:
	ValuatorBinding(string inKey,Fl_Valuator *inWidget) {
		key = inKey;
		widget = inWidget;
	}
	
	virtual void setWidget(string value) {
		widget->value(atof(value.c_str()));
	}

	virtual string readWidget() {
		char out[20];
		sprintf(out,"%lg",widget->value());
		return string(out); 
	}
};

class ColorButtonBinding : public UIBinding {
	Fl_Light_Button *widget;
public:
	ColorButtonBinding(string inKey,Fl_Light_Button *inWidget) {
		key = inKey;
		widget = inWidget;
	}
	
	virtual void setWidget(string value) {
		// double r,g,b;
		// sscanf(value.c_str(),"%lg,%lg,%lg",&r,&g,&b);
		widget->color(atoi(value.c_str()));
	}

	virtual string readWidget() {
		char out[20];

		// unsigned char r,g,b;

		sprintf(out,"%d",(int)Fl::get_color(widget->color()));
		return string(out); 
	}
};

class CheckBinding : public UIBinding {
	Fl_Button *widget;
public:
	CheckBinding(string inKey,Fl_Button *inWidget) {
		key = inKey;
		widget = inWidget;
	}
	
	virtual void setWidget(string value) {
		widget->value(atoi(value.c_str()));
	}

	virtual string readWidget() {
		char out[20];
		sprintf(out,"%d",widget->value());
		return string(out); 
	}
};

class RadioBinding : public UIBinding {
	list<Fl_Button *> widgets;
	Fl_Button *deflt;
public:
	RadioBinding(string inKey) {
		key = inKey;
		deflt = NULL;
	}
	
	void add(Fl_Button *inWidget,bool isDefault) {
		widgets.push_back(inWidget);
		if(isDefault) {
			deflt = inWidget;
		}
	}

	virtual void setWidget(string value) {
		bool found = false;

		list<Fl_Button *>::iterator j;
		for(j = widgets.begin();j!=widgets.end();j++) {
			Fl_Button *widget = *j;
			string bval((char *)widget->user_data());

			if(bval == value) {
				widget->value(1);
				found = true;
			}
			else {
				widget->value(0);
			}
		}

		if(!found && deflt!=NULL)
			deflt->value(1);
	}

	virtual string readWidget() {
		list<Fl_Button *>::iterator j;
		for(j = widgets.begin();j!=widgets.end();j++) {
			Fl_Button *widget = *j;
			if(widget->value()) {
				return string((char *)widget->user_data());
			}
		}
		return "";
	}
};

class ChoiceBinding : public UIBinding {
	Fl_Menu_ *menu;
	int deflt;
public:
	ChoiceBinding(string inKey,Fl_Menu_ *inMenu,int inDefault = -1) {
		key = inKey;
		menu = inMenu;
		deflt = inDefault;
	}
	
	virtual void setWidget(string value) {
		bool found = false;

		const Fl_Menu_Item *mi = menu->menu();
		for(int i=0;i<menu->size();i++) {			
			if(value == string((char *)mi->user_data())) {
				menu->value(i);
				found = true;
				break;
			}
			mi++;
		}

		if(!found && deflt >= 0)
			menu->value(deflt);
	}

	virtual string readWidget() {
		if(menu->mvalue()) 
			return string((char *)menu->mvalue()->user_data());
		return "";
	}
};

void UIBinder::bind(string key,Fl_Valuator *widget) {
	ValuatorBinding *vb = new ValuatorBinding(key,widget);
	keyMap[key] = vb;
	widMap[widget] = vb;
}

void UIBinder::bindColor(string key,Fl_Light_Button *widget) {
	ColorButtonBinding *vb = new ColorButtonBinding(key,widget);
	keyMap[key] = vb;
	widMap[widget] = vb;
}

void UIBinder::bind(string key,Fl_Button *widget) {
	CheckBinding *vb = new CheckBinding(key,widget);
	keyMap[key] = vb;
	widMap[widget] = vb;
}

void UIBinder::bindRadios(string key,Fl_Group *container,Fl_Button *deflt) {
	RadioBinding *rb = new RadioBinding(key);
	keyMap[key] = rb;
	
	// Go through the group and pull out all the buttons
	for(int i=0;i<container->children();i++) {
		Fl_Widget *widget = container->child(i);
		if(widget->type() == 'f') {
			// Fl_Button *button = dynamic_cast<Fl_Button*>(widget);
			Fl_Button *button = (Fl_Button*) widget;
			if(button) {	
				rb->add(button,(button==deflt));
				widMap[button] = rb;
			}
		}
	}
}

void UIBinder::bindChoice(string key,Fl_Menu_ *menu,int deflt) {
	ChoiceBinding *cb = new ChoiceBinding(key,menu,deflt);
	keyMap[key] = cb;
	widMap[menu] = cb;
}

bool UIBinder::readWidget(Fl_Widget *w,string &key,string &value) {
	UIBinding *binding = widMap[w];
	if(binding) {
		key = binding->key;
		value = binding->readWidget();
		return true;
	}
	return false;
}

void UIBinder::setWidget(string key,string value) {
	UIBinding *binding = keyMap[key];
	if(binding) {
		binding->setWidget(value);
	}
}


