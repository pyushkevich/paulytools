// generated by Fast Light User Interface Designer (fluid) version 1.00

#include "DevToolApp.H"

inline void BasicWindow::cb_Quit_i(Fl_Menu_*, void*) {
  exit(0);
}
void BasicWindow::cb_Quit(Fl_Menu_* o, void* v) {
  ((BasicWindow*)(o->parent()->user_data()))->cb_Quit_i(o,v);
}

inline void BasicWindow::cb_Edit_i(Fl_Menu_*, void*) {
  app->EditPreferences();
}
void BasicWindow::cb_Edit(Fl_Menu_* o, void* v) {
  ((BasicWindow*)(o->parent()->user_data()))->cb_Edit_i(o,v);
}

inline void BasicWindow::cb_Tags_i(Fl_Menu_*, void*) {
  app->DoTagsBrowser();
}
void BasicWindow::cb_Tags(Fl_Menu_* o, void* v) {
  ((BasicWindow*)(o->parent()->user_data()))->cb_Tags_i(o,v);
}

inline void BasicWindow::cb_Make_i(Fl_Menu_*, void*) {
  app->DoMakeParser();
}
void BasicWindow::cb_Make(Fl_Menu_* o, void* v) {
  ((BasicWindow*)(o->parent()->user_data()))->cb_Make_i(o,v);
}

inline void BasicWindow::cb_Split_i(Fl_Menu_*, void*) {
  app->Split();
}
void BasicWindow::cb_Split(Fl_Menu_* o, void* v) {
  ((BasicWindow*)(o->parent()->user_data()))->cb_Split_i(o,v);
}

inline void BasicWindow::cb_Merge_i(Fl_Menu_*, void*) {
  app->Merge();
}
void BasicWindow::cb_Merge(Fl_Menu_* o, void* v) {
  ((BasicWindow*)(o->parent()->user_data()))->cb_Merge_i(o,v);
}

inline void BasicWindow::cb_Tags1_i(Fl_Menu_*, void*) {
  app->HelpTagsBrowser();
}
void BasicWindow::cb_Tags1(Fl_Menu_* o, void* v) {
  ((BasicWindow*)(o->parent()->user_data()))->cb_Tags1_i(o,v);
}

inline void BasicWindow::cb_Make1_i(Fl_Menu_*, void*) {
  app->HelpMakeParser();
}
void BasicWindow::cb_Make1(Fl_Menu_* o, void* v) {
  ((BasicWindow*)(o->parent()->user_data()))->cb_Make1_i(o,v);
}

Fl_Menu_Item BasicWindow::menu_[] = {
 {"File", 0,  0, 0, 64, 0, 0, 12, 0},
 {"Quit", 0,  (Fl_Callback*)BasicWindow::cb_Quit, 0, 0, 0, 0, 12, 0},
 {0},
 {"Settings", 0,  0, 0, 64, 0, 0, 12, 0},
 {"Edit ..", 0,  (Fl_Callback*)BasicWindow::cb_Edit, 0, 128, 0, 0, 12, 0},
 {"Save", 0,  0, 0, 0, 0, 0, 12, 0},
 {"Save As ..", 0,  0, 0, 128, 0, 0, 12, 0},
 {"Save As Default", 0,  0, 0, 0, 0, 0, 12, 0},
 {0},
 {"Utilities", 0,  0, 0, 64, 0, 0, 12, 0},
 {"Tags Browser", 0,  (Fl_Callback*)BasicWindow::cb_Tags, 0, 0, 0, 0, 12, 0},
 {"Make Parser", 0,  (Fl_Callback*)BasicWindow::cb_Make, 0, 0, 0, 0, 12, 0},
 {"Multifile Search", 0,  0, 0, 16, 0, 0, 12, 0},
 {0},
 {"Windows", 0,  0, 0, 64, 0, 0, 12, 0},
 {"Split", 0,  (Fl_Callback*)BasicWindow::cb_Split, 0, 0, 0, 0, 12, 0},
 {"Merge", 0,  (Fl_Callback*)BasicWindow::cb_Merge, 0, 0, 0, 0, 12, 0},
 {0},
 {"Help", 0,  0, 0, 64, 0, 0, 12, 0},
 {"Tags Browser ..", 0,  (Fl_Callback*)BasicWindow::cb_Tags1, 0, 0, 0, 0, 12, 0},
 {"Make Parser ..", 0,  (Fl_Callback*)BasicWindow::cb_Make1, 0, 0, 0, 0, 12, 0},
 {"Multifile Search ..", 0,  0, 0, 16, 0, 0, 12, 0},
 {0},
 {0}
};

BasicWindow::BasicWindow():Fl_Double_Window(0,0) {
  Fl_Group* w;
  { Fl_Group* o = mGroup = new Fl_Group(0, 0, 421, 310);
    w = o;
    o->user_data((void*)(this));
    { Fl_Menu_Bar* o = new Fl_Menu_Bar(0, 0, 420, 25);
      o->box(FL_THIN_UP_BOX);
      o->menu(menu_);
    }
    { Fl_Group* o = mMainGroup = new Fl_Group(0, 25, 420, 285);
      { Fl_Tile* o = mTile = new Fl_Tile(0, 25, 420, 285);
        { Fl_Box* o = new Fl_Box(40, 90, 330, 140, "label");
          o->hide();
          Fl_Group::current()->resizable(o);
        }
        o->end();
        Fl_Group::current()->resizable(o);
      }
      o->end();
      Fl_Group::current()->resizable(o);
    }
    o->end();
  }
  resize(x(),y(),mGroup->w(),mGroup->h());
  resizable(mGroup);
}
DevToolsApp* app;

inline void PreferencesWindow::cb_OK_i(Fl_Button*, void*) {
  Done(1);
}
void PreferencesWindow::cb_OK(Fl_Button* o, void* v) {
  ((PreferencesWindow*)(o->parent()->parent()->user_data()))->cb_OK_i(o,v);
}

inline void PreferencesWindow::cb_Cancel_i(Fl_Button*, void*) {
  Done(0);
}
void PreferencesWindow::cb_Cancel(Fl_Button* o, void* v) {
  ((PreferencesWindow*)(o->parent()->parent()->user_data()))->cb_Cancel_i(o,v);
}

PreferencesWindow::PreferencesWindow():Fl_Window(0,0) {
  Fl_Group* w;
  { Fl_Group* o = mGroup = new Fl_Group(0, 0, 221, 174);
    w = o;
    o->user_data((void*)(this));
    { Fl_Group* o = mButtonGroup = new Fl_Group(-1, 148, 222, 27);
      { Fl_Button* o = new Fl_Button(3, 151, 49, 20, "OK");
        o->box(FL_THIN_UP_BOX);
        o->labelsize(12);
        o->callback((Fl_Callback*)cb_OK);
      }
      { Fl_Button* o = new Fl_Button(55, 151, 49, 20, "Cancel");
        o->box(FL_THIN_UP_BOX);
        o->labelsize(12);
        o->callback((Fl_Callback*)cb_Cancel);
      }
      o->end();
    }
    { Fl_Box* o = new Fl_Box(105, 0, 115, 149, "label");
      o->hide();
      Fl_Group::current()->resizable(o);
    }
    o->end();
  }
  resize(x(),y(),mGroup->w(),mGroup->h());
  resizable(mGroup);
}
