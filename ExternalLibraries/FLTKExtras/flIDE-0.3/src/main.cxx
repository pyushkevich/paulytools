#include <FL/fl_ask.H>
#include "DevToolApp.H"
#include "DevTool.H"
#include "MakeParser.H"
#include "TagsBrowser.H"


DevToolsApp::DevToolsApp(void)
{
	mTagsBrowser = 0;
	mMakeParser = 0;	
	mMakeParserWindow = 0;
	mTagsBrowserWindow = 0;
	mSharedWindow = 0;
	mEmptyWindow = 0;	
}

BasicWindow* DevToolsApp::DoDevTool(DevTool* tool)
{
	bool created=false;
	BasicWindow* thisToolsWindow;
	
	if (mEmptyWindow) {
		thisToolsWindow = mEmptyWindow;
		mEmptyWindow = 0;
	}else{
		thisToolsWindow = new BasicWindow;
		created=true;
	}
	thisToolsWindow->size(		
		tool->mMainGroup->w(),
		thisToolsWindow->mMainGroup->y()+tool->mMainGroup->h()
	);
	tool->mMainGroup->position(
		thisToolsWindow->mMainGroup->x(),
		thisToolsWindow->mMainGroup->y());
	thisToolsWindow->mMainGroup->add(tool->mMainGroup);

	if (created) {
		thisToolsWindow->show();
	}else{
		thisToolsWindow->redraw();
	}
	return thisToolsWindow;
}

void DevToolsApp::DoTagsBrowser(void)
{
	if (mSharedWindow) return;
	if (mTagsBrowser==0) {
		mTagsBrowser=new TagsBrowser;
		mTagsBrowserWindow = DoDevTool(mTagsBrowser);
		mTagsBrowserWindow->label("Tags Browser");
	}
}

void DevToolsApp::DoMakeParser()
{
	if (mSharedWindow) return;
	if (mMakeParser==0) {
		mMakeParser=new MakeParser;
		mMakeParserWindow = DoDevTool(mMakeParser);
		mMakeParserWindow->label("Make Parser");
	}
}

void DevToolsApp::HelpMakeParser(void)
{
	fl_alert(
		"Not implemented yet. TODO: add Help window to MakeParser.fl (from flmake),"
		"and call it from main.cxx: DevToolsApp::HelpMakeParser"
	);
}

void DevToolsApp::HelpTagsBrowser(void)
{
	fl_alert(
		"Not implemented yet. TODO: add Help window to TagsBrowser.fl (from fltags),"
		"and call it from main.cxx: DevToolsApp::HelpTagsBrowser"
	);
}

void DevToolsApp::Merge(void)
{
	if (
		mTagsBrowser && mMakeParser &&
		mSharedWindow == 0) 
	{
		mSharedWindow = mTagsBrowserWindow;
		mSharedWindow->label("flIDE");

		Fl_Group* g= new		Fl_Group(	
				mSharedWindow->mTile->x(),mSharedWindow->mTile->y(),
				mSharedWindow->mTile->w(),mSharedWindow->mTile->h()/2);

		mTagsBrowser->mMainGroup->resize(g->x(),g->y(),g->w(),g->h()-2);

		g->resizable(mTagsBrowser->mMainGroup);
		Fl_Box* b2=mSharedWindow->mSplitter = 
			new	Fl_Box(g->x()+4,g->y()+g->h()-2,g->w()-8,2);
		b2->box(FL_THIN_DOWN_BOX);

		g->add(mTagsBrowser->mMainGroup);
		g->add(b2);
		g->end();
	
		mMakeParser->mMainGroup->resize(
			g->x(),g->y()+g->h(),
			g->w(),mSharedWindow->mTile->h()-g->h());

		mSharedWindow->mTile->add(g);				
		mSharedWindow->mTile->add(mMakeParser->mMainGroup);				
		delete mMakeParserWindow;

		mTagsBrowserWindow = 0;
		mMakeParserWindow = 0;

		mTagsBrowser->mMainGroup->box(FL_FLAT_BOX);
		mMakeParser->mMainGroup->box(FL_FLAT_BOX);

		mSharedWindow->redraw();		
	}
}

void DevToolsApp::Split(void)
{
	if (
		mTagsBrowser && mMakeParser &&
		mSharedWindow) {
			mTagsBrowserWindow = mSharedWindow;
			mMakeParserWindow = new BasicWindow;

		mTagsBrowserWindow->label("Tags Browser");
		mMakeParserWindow->label("Make Parser");

		mTagsBrowserWindow->mMainGroup->add(mTagsBrowser->mMainGroup);
		mMakeParserWindow->mMainGroup->add(mMakeParser->mMainGroup);
		mTagsBrowser->mMainGroup->resize(
			mTagsBrowserWindow->mMainGroup->x(),mTagsBrowserWindow->mMainGroup->y(),
			mTagsBrowserWindow->mMainGroup->w(),mTagsBrowserWindow->mMainGroup->h());
		mMakeParser->mMainGroup->resize(
			mMakeParserWindow->mMainGroup->x(),mMakeParserWindow->mMainGroup->y(),
			mMakeParserWindow->mMainGroup->w(),mMakeParserWindow->mMainGroup->h());

		mTagsBrowserWindow->mTile->remove(mTagsBrowserWindow->mSplitter->parent());
		mTagsBrowserWindow->mTile->remove(mTagsBrowserWindow->mSplitter);

		mMakeParserWindow->show();

		mTagsBrowserWindow->redraw();
		
		mSharedWindow = 0;
	}
}

void PreferencesWindow::Init(void)
{
/*	int hh = mPreferences->h();
	g1->end();
	
	mPreferences->size(g1->w(),g1->h());
	mPreferences->add(g1);
*/      done=false;
	Fl_Group* g1=app->mTagsBrowser->Preferences();
	Fl_Box* b;
	g1->add(
		b=new Fl_Box(
			g1->x()+g1->w(),
			g1->y()+g1->h(),
			0,0));

	int hh = g1->h();
	
	g1->resizable(b);
	g1->end();

	Fl_Group* g2=app->mMakeParser->Preferences();
	g2->add(
		b=new Fl_Box(
			g2->x()+g2->w(),
			g2->y()+g2->h(),
			0,0));
	g2->resizable(b);
	g2->end();


	size(g1->w(),g1->h()+g2->h()+mButtonGroup->h());
	add(g1);
	add(g2);
	g2->position(0,hh);
}

void PreferencesWindow::Done(int ok)
{
 done=true;
 if(!ok) return;
}

void DevToolsApp::EditPreferences(void)
{
	PreferencesWindow* w = new PreferencesWindow;
	w->Init();
	w->show();
	while (!w->done && w->visible()) {
		Fl::wait();
	}
	delete w;

}

main(int argc, char **argv) 
{
	app = new DevToolsApp;
	
	app->DoMakeParser();
	app->DoTagsBrowser();

	app->mMakeParser->Init(argc,argv);
	app->mTagsBrowser->Init(argc,argv);

	Fl::run();
}
