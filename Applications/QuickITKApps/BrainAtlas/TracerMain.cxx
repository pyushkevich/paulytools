/**
 * Tracer application: trace curves on the grey matter surface
 */

#include "TracerData.h"
#include "TracerUserInterfaceLogic.h"

int main(int argc, char *argv[])
{
  // Initialize FLTK
  Fl::visual(FL_DOUBLE|FL_INDEX);
  Fl::gl_visual(FL_RGB);  
  Fl::background(236,233,216);
  
  // Create the main window
  TracerUserInterfaceLogic *ui = new TracerUserInterfaceLogic();
  ui->MakeWindow();

  // Show the windows
  ui->ShowWindows();

  // Create the tracer data
  TracerData *data = new TracerData();

  // Assign the data to the UI
  ui->SetTracerData(data);

  // Run the FL loop
  Fl::run();

  // Clean up
  delete ui;
  delete data;
}
