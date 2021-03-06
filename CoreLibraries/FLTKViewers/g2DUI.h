// generated by Fast Light User Interface Designer (fluid) version 1.00

#include "Graph2DWindow.h"
#include "graph2d.h"
#include <FL/Fl.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Slider.H>
#include <FL/Fl_Value_Input.H>
#include <FL/Fl_Value_Output.H>
#include <FL/Fl_Window.H>

class Graph2DUI {
public:
  Fl_Window* makeWindow();
  Fl_Window *winMain;
  Fl_Group *groupX;
  Fl_Value_Input *inpXMin;
  Fl_Value_Input *inpXMax;
  Fl_Group *groupY;
  Fl_Value_Input *inpYMin;
  Fl_Value_Input *inpYMax;
  Fl_Value_Output *outRes;
  Fl_Slider *inpRes;
private:
  inline void cb_inpRes_i(Fl_Slider*, void*);
  static void cb_inpRes(Fl_Slider*, void*);
public:
  Fl_Check_Button *bSeeGrid;
private:
  inline void cb_bSeeGrid_i(Fl_Check_Button*, void*);
  static void cb_bSeeGrid(Fl_Check_Button*, void*);
public:
  Fl_Check_Button *bSeeAxis;
private:
  inline void cb_bSeeAxis_i(Fl_Check_Button*, void*);
  static void cb_bSeeAxis(Fl_Check_Button*, void*);
public:
  Fl_Check_Button *bSeeLabels;
private:
  inline void cb_bSeeLabels_i(Fl_Check_Button*, void*);
  static void cb_bSeeLabels(Fl_Check_Button*, void*);
public:
  Fl_Button *bNext;
private:
  inline void cb_bNext_i(Fl_Button*, void*);
  static void cb_bNext(Fl_Button*, void*);
public:
  Fl_Button *bBack;
private:
  inline void cb_bBack_i(Fl_Button*, void*);
  static void cb_bBack(Fl_Button*, void*);
public:
  Fl_Button *bZoomOut;
private:
  inline void cb_bZoomOut_i(Fl_Button*, void*);
  static void cb_bZoomOut(Fl_Button*, void*);
public:
  Fl_Button *bZoomIn;
private:
  inline void cb_bZoomIn_i(Fl_Button*, void*);
  static void cb_bZoomIn(Fl_Button*, void*);
public:
  Fl_Button *bCompute;
private:
  inline void cb_bCompute_i(Fl_Button*, void*);
  static void cb_bCompute(Fl_Button*, void*);
public:
  Graph2DPalette *winPalette;
  Fl_Button *bResetPalette;
private:
  inline void cb_bResetPalette_i(Fl_Button*, void*);
  static void cb_bResetPalette(Fl_Button*, void*);
public:
  Graph2DWindow *winGraph;
  Fl_Value_Output *outX;
  Fl_Value_Output *outY;
  Fl_Value_Output *outFxy;
  Fl_Window *winComputing;
  Fl_Slider *sldProgress;
  Graph2D *graph;
};
