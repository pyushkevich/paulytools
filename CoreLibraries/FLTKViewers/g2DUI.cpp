// generated by Fast Light User Interface Designer (fluid) version 1.00

#include "g2DUI.h"

inline void Graph2DUI::cb_inpRes_i(Fl_Slider*, void*) {
  outRes->value(1 << ((int)inpRes->value()));
}
void Graph2DUI::cb_inpRes(Fl_Slider* o, void* v) {
  ((Graph2DUI*)(o->parent()->parent()->parent()->parent()->user_data()))->cb_inpRes_i(o,v);
}

inline void Graph2DUI::cb_bSeeGrid_i(Fl_Check_Button*, void*) {
  graph->showGrid(bSeeGrid->value());
}
void Graph2DUI::cb_bSeeGrid(Fl_Check_Button* o, void* v) {
  ((Graph2DUI*)(o->parent()->parent()->parent()->parent()->user_data()))->cb_bSeeGrid_i(o,v);
}

inline void Graph2DUI::cb_bSeeAxis_i(Fl_Check_Button*, void*) {
  graph->showAxis(bSeeAxis->value());
}
void Graph2DUI::cb_bSeeAxis(Fl_Check_Button* o, void* v) {
  ((Graph2DUI*)(o->parent()->parent()->parent()->parent()->user_data()))->cb_bSeeAxis_i(o,v);
}

inline void Graph2DUI::cb_bSeeLabels_i(Fl_Check_Button*, void*) {
  graph->showLabels(bSeeLabels->value());
}
void Graph2DUI::cb_bSeeLabels(Fl_Check_Button* o, void* v) {
  ((Graph2DUI*)(o->parent()->parent()->parent()->parent()->user_data()))->cb_bSeeLabels_i(o,v);
}

inline void Graph2DUI::cb_bNext_i(Fl_Button*, void*) {
  graph->next();
}
void Graph2DUI::cb_bNext(Fl_Button* o, void* v) {
  ((Graph2DUI*)(o->parent()->parent()->parent()->parent()->user_data()))->cb_bNext_i(o,v);
}

inline void Graph2DUI::cb_bBack_i(Fl_Button*, void*) {
  graph->back();
}
void Graph2DUI::cb_bBack(Fl_Button* o, void* v) {
  ((Graph2DUI*)(o->parent()->parent()->parent()->parent()->user_data()))->cb_bBack_i(o,v);
}

inline void Graph2DUI::cb_bZoomOut_i(Fl_Button*, void*) {
  graph->zoomOut();
}
void Graph2DUI::cb_bZoomOut(Fl_Button* o, void* v) {
  ((Graph2DUI*)(o->parent()->parent()->parent()->parent()->user_data()))->cb_bZoomOut_i(o,v);
}

inline void Graph2DUI::cb_bZoomIn_i(Fl_Button*, void*) {
  graph->zoomIn();
}
void Graph2DUI::cb_bZoomIn(Fl_Button* o, void* v) {
  ((Graph2DUI*)(o->parent()->parent()->parent()->parent()->user_data()))->cb_bZoomIn_i(o,v);
}

inline void Graph2DUI::cb_bCompute_i(Fl_Button*, void*) {
  graph->compute();
}
void Graph2DUI::cb_bCompute(Fl_Button* o, void* v) {
  ((Graph2DUI*)(o->parent()->parent()->parent()->parent()->user_data()))->cb_bCompute_i(o,v);
}

inline void Graph2DUI::cb_bResetPalette_i(Fl_Button*, void*) {
  winPalette->reset();
}
void Graph2DUI::cb_bResetPalette(Fl_Button* o, void* v) {
  ((Graph2DUI*)(o->parent()->parent()->parent()->user_data()))->cb_bResetPalette_i(o,v);
}

Fl_Window* Graph2DUI::makeWindow() {
  Fl_Window* w;
  { Fl_Window* o = winMain = w = new Fl_Window(602, 532, "Graph 2D - Plotting functions in 2D Domain");
    o->color(32);
    o->user_data((void*)(this));
    { Fl_Group* o = new Fl_Group(465, 5, 135, 450);
      o->box(FL_UP_BOX);
      { Fl_Group* o = new Fl_Group(465, 5, 135, 440);
        { Fl_Group* o = groupX = new Fl_Group(470, 25, 125, 50, "X Coordinate:");
          o->align(5);
          { Fl_Value_Input* o = inpXMin = new Fl_Value_Input(505, 30, 90, 20, "from:");
            o->labelsize(12);
            o->textsize(10);
          }
          { Fl_Value_Input* o = inpXMax = new Fl_Value_Input(505, 50, 90, 20, "to:");
            o->labelsize(12);
            o->textsize(10);
          }
          o->end();
        }
        { Fl_Group* o = groupY = new Fl_Group(470, 95, 125, 50, "Y Coordinate:");
          o->align(5);
          { Fl_Value_Input* o = inpYMin = new Fl_Value_Input(505, 100, 90, 20, "from:");
            o->labelsize(12);
            o->textsize(10);
          }
          { Fl_Value_Input* o = inpYMax = new Fl_Value_Input(505, 120, 90, 20, "to:");
            o->labelsize(12);
            o->textsize(10);
          }
          o->end();
        }
        { Fl_Group* o = new Fl_Group(470, 165, 125, 40, "Resolution:");
          o->align(5);
          { Fl_Value_Output* o = outRes = new Fl_Value_Output(485, 170, 40, 20);
            o->value(256);
            o->textsize(12);
          }
          { Fl_Slider* o = inpRes = new Fl_Slider(530, 170, 65, 20);
            o->type(1);
            o->selection_color(62);
            o->minimum(6);
            o->maximum(10);
            o->step(1);
            o->value(8);
            o->callback((Fl_Callback*)cb_inpRes);
          }
          o->end();
        }
        { Fl_Group* o = new Fl_Group(470, 220, 110, 65, "Display Options:");
          o->align(5);
          { Fl_Check_Button* o = bSeeGrid = new Fl_Check_Button(490, 240, 85, 20, "show grid");
            o->down_box(FL_DIAMOND_DOWN_BOX);
            o->selection_color(62);
            o->labelsize(12);
            o->callback((Fl_Callback*)cb_bSeeGrid);
          }
          { Fl_Check_Button* o = bSeeAxis = new Fl_Check_Button(490, 220, 85, 20, "show axis");
            o->down_box(FL_DIAMOND_DOWN_BOX);
            o->value(1);
            o->selection_color(62);
            o->labelsize(12);
            o->callback((Fl_Callback*)cb_bSeeAxis);
          }
          { Fl_Check_Button* o = bSeeLabels = new Fl_Check_Button(490, 260, 85, 20, "show x,y");
            o->down_box(FL_DIAMOND_DOWN_BOX);
            o->selection_color(62);
            o->labelsize(12);
            o->callback((Fl_Callback*)cb_bSeeLabels);
          }
          o->end();
        }
        { Fl_Group* o = new Fl_Group(470, 305, 120, 135, "Zoom Control:");
          o->align(5);
          { Fl_Button* o = bNext = new Fl_Button(535, 415, 30, 20, "@->");
            o->labeltype(FL_SYMBOL_LABEL);
            o->labelcolor(130);
            o->callback((Fl_Callback*)cb_bNext);
            o->deactivate();
          }
          { Fl_Button* o = bBack = new Fl_Button(495, 415, 30, 20, "@<-");
            o->labeltype(FL_SYMBOL_LABEL);
            o->labelcolor(130);
            o->callback((Fl_Callback*)cb_bBack);
            o->deactivate();
          }
          { Fl_Button* o = bZoomOut = new Fl_Button(480, 310, 100, 20, "Zoom Out (2:1)");
            o->labelsize(12);
            o->callback((Fl_Callback*)cb_bZoomOut);
          }
          { Fl_Button* o = bZoomIn = new Fl_Button(480, 335, 100, 20, "Zoom In (1:2)");
            o->labelsize(12);
            o->callback((Fl_Callback*)cb_bZoomIn);
          }
          { Fl_Button* o = bCompute = new Fl_Button(490, 380, 80, 25, "Compute");
            o->labelfont(3);
            o->callback((Fl_Callback*)cb_bCompute);
          }
          o->end();
        }
        o->end();
      }
      { Fl_Group* o = new Fl_Group(465, 445, 135, 10);
        o->end();
        Fl_Group::current()->resizable(o);
      }
      o->end();
    }
    { Fl_Group* o = new Fl_Group(5, 460, 455, 70);
      o->box(FL_UP_BOX);
      { Graph2DPalette* o = winPalette = new Graph2DPalette(90, 465, 360, 60);
        o->box(FL_DOWN_BOX);
        o->color(4);
        o->selection_color(4);
        o->align(133);
        Fl_Group::current()->resizable(o);
      }
      { Fl_Group* o = new Fl_Group(10, 495, 75, 30, "Intensity Windowing");
        o->align(133);
        { Fl_Button* o = bResetPalette = new Fl_Button(20, 500, 55, 20, "Reset");
          o->labelsize(12);
          o->callback((Fl_Callback*)cb_bResetPalette);
        }
        o->end();
      }
      o->end();
    }
    { Graph2DWindow* o = winGraph = new Graph2DWindow(5, 5, 455, 450);
      o->box(FL_UP_BOX);
      o->color(220);
      Fl_Group::current()->resizable(o);
    }
    { Fl_Group* o = new Fl_Group(465, 460, 135, 70);
      o->box(FL_UP_BOX);
      { Fl_Value_Output* o = outX = new Fl_Value_Output(500, 465, 95, 20, "x:");
        o->labelsize(12);
        o->textsize(10);
      }
      { Fl_Value_Output* o = outY = new Fl_Value_Output(500, 485, 95, 20, "y:");
        o->labelsize(12);
        o->textsize(10);
      }
      { Fl_Value_Output* o = outFxy = new Fl_Value_Output(500, 505, 95, 20, "f(x,y):");
        o->labelsize(12);
        o->textsize(10);
      }
      o->end();
    }
    o->end();
  }
  { Fl_Window* o = winComputing = w = new Fl_Window(322, 75, "Graph 2D");
    o->box(FL_ENGRAVED_BOX);
    o->user_data((void*)(this));
    { Fl_Slider* o = sldProgress = new Fl_Slider(5, 45, 310, 20);
      o->type(3);
      o->box(FL_THIN_DOWN_BOX);
      o->selection_color(137);
      o->labeltype(FL_ENGRAVED_LABEL);
      o->align(5);
    }
    { Fl_Group* o = new Fl_Group(10, 25, 300, 15, "Computing Graph.  Please Wait...");
      o->labelfont(3);
      o->align(5);
      o->end();
    }
    o->set_modal();
    o->clear_border();
    o->end();
  }
  return w;
}
