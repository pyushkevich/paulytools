// generated by Fast Light User Interface Designer (fluid) version 1.00

#ifndef MakeParser_H
#define MakeParser_H
#include <FL/Fl.H>
#include <stdio.h>
#include <sys/types.h>
#include <malloc.h>
#include <regex.h>
#include <FL/Fl.H>
#include <FL/fl_draw.H>
#include <FL/Fl_Hold_Browser.H>
#include "DevTool.H"

class NiceBrowser : public Fl_Hold_Browser {
public:
  NiceBrowser(int x,int y,int w,int h,const char *l=0);
  void item_draw(void* v, int x, int y, int w, int h) const;
};
#include <FL/Fl_Box.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Tabs.H>
#include <FL/Fl_Input.H>

class MakeParser : public DevTool {
public:
  MakeParser();
private:
  Fl_Group *mLeftMakeMenu;
public:
  Fl_Button *mPrevError;
private:
  inline void cb_mPrevError_i(Fl_Button*, void*);
  static void cb_mPrevError(Fl_Button*, void*);
public:
  Fl_Button *mNextError;
private:
  inline void cb_mNextError_i(Fl_Button*, void*);
  static void cb_mNextError(Fl_Button*, void*);
public:
  Fl_Button *mPrevWarning;
private:
  inline void cb_mPrevWarning_i(Fl_Button*, void*);
  static void cb_mPrevWarning(Fl_Button*, void*);
public:
  Fl_Button *mNextWarning;
private:
  inline void cb_mNextWarning_i(Fl_Button*, void*);
  static void cb_mNextWarning(Fl_Button*, void*);
  Fl_Group *mTopMakeMenu;
public:
  Fl_Button *mMakeButton;
private:
  inline void cb_mMakeButton_i(Fl_Button*, void*);
  static void cb_mMakeButton(Fl_Button*, void*);
public:
  Fl_Button *mMakeCleanButton;
private:
  inline void cb_mMakeCleanButton_i(Fl_Button*, void*);
  static void cb_mMakeCleanButton(Fl_Button*, void*);
public:
  Fl_Button *mStopButton;
private:
  inline void cb_mStopButton_i(Fl_Button*, void*);
  static void cb_mStopButton(Fl_Button*, void*);
public:
  Fl_Button *mRunButton;
private:
  inline void cb_mRunButton_i(Fl_Button*, void*);
  static void cb_mRunButton(Fl_Button*, void*);
public:
  Fl_Tabs *mTabbedOutput;
  Fl_Group *mErrorTab;
  NiceBrowser *mErrorBrowser;
private:
  inline void cb_mErrorBrowser_i(NiceBrowser*, void*);
  static void cb_mErrorBrowser(NiceBrowser*, void*);
public:
  Fl_Group *mRunTab;
  NiceBrowser *mRunTargetsBrowser;
private:
  inline void cb_mRunTargetsBrowser_i(NiceBrowser*, void*);
  static void cb_mRunTargetsBrowser(NiceBrowser*, void*);
public:
  Fl_Group* Preferences();
private:
  void SelectErrorLine();
public:
  void Init(int argc,char** argv);
private:
  void Run();
  void Make();
  void MakeClean();
  void StopMake();
  void NextErrorMsg();
  void PrevErrorMsg();
  void NextWarningMsg();
  void PrevWarningMsg();
public:
  void MakeThread(void);
private:
  void RunMake(char* cmdin);
  void NextMessage(const char*);
  void PrevMessage(const char*);
public:
  void LookForExecutable(const char*);
  char* LookForArg(const char* arg,const char* makeoutline);
private:
  const char* SplitLineRegExp(char* intext,char* filename,char *text,int* linenumber,int *msgtype,int concatdir);
public:
  void InitRegex(void);
private:
  int CompileRegEx(regex_t *preg,const char* regex);
  int TestRegSet(const char *buff,int nmatch, regmatch_t pmatch[]);
public:
  void InitMessageHeaders(void);
private:
  int SplitRunTargetsLine(const char *text,char *exename,char *params);
  void ChangeRunTargetsParams(void);
  void RunWithConsole(void);
  void RunInDebugger(void);
  int FindTokInBrowser(Fl_Browser *mBrowser,const char *identTok);
public:
  char mMakeArgs[1024];
  char mCurrDir[1024];
  char mInvEditor[1024];
private:
  int mNPreg;
  int mWNPreg;
  int mINPreg;
  regex_t mPreg[64];
  char *mMsgHeaders[64];
  int mNumMsgHeaders;
  int mEditOnMsgBrowse;
  int mCleanRunTargets;
};
#endif
