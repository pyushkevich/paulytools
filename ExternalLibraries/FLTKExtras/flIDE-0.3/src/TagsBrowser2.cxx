#if defined(WIN32)
#define popen _popen
#endif

#include <FL/Fl_ToggleTree.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Radio_Button.H>
#include <FL/Fl_Scroll.H>
#include <FL/Fl.H>

#include <FL/folder_small.xpm>
#include <FL/file_small.xpm>
#include "function_small.xpm"
#include "point_small.xpm"
#include "box_small.xpm"
#include "boxes_small.xpm"
#include "global_small.xpm"

#include "TagsBrowser.H"
#include "FL/Fl_ToggleTreeResizer.H"

#include <stdio.h>
#include <stdlib.h>

const int MAX_SEARCH_ENTRY = 25;

Fl_Pixmap* fileSmall = 0;
Fl_Pixmap* folderSmall = 0;
Fl_Pixmap* pointSmall = 0;
Fl_Pixmap* functionSmall = 0;
Fl_Pixmap* boxSmall = 0;
Fl_Pixmap* boxesSmall = 0;
Fl_Pixmap* globalSmall = 0;

class TagNodeInfo
{
public:
	TagNodeInfo(int t)
	{
		type=t;
	}
	int type;
};

#define kTagNodeGlobal		1
#define kTagNodeGroup			2
#define kTagNodeMembers		10
#define kTagNodeMethods 	11
#define kTagNodeSource  	12
#define kTagNodeRef				20
#define kTagNodeCategory	50
#define kTagNodeDeleting	128
#define kTagNodeSources		256

class GroupNodeInfo:public TagNodeInfo
{
public:
	GroupNodeInfo():TagNodeInfo(kTagNodeGroup){
		methods = 0;
		members = 0;
	}
	Fl_ToggleNode* methods;
	Fl_ToggleNode* members;
};

class SourceNodeInfo:public TagNodeInfo
{
public:
	SourceNodeInfo():TagNodeInfo(kTagNodeSource){}
};

class MembersNodeInfo:public TagNodeInfo
{
public:
	MembersNodeInfo():TagNodeInfo(kTagNodeMembers){}
};

class MethodsNodeInfo:public TagNodeInfo
{
public:
	MethodsNodeInfo():TagNodeInfo(kTagNodeMethods){}
};

class RefNodeInfo:public TagNodeInfo
{
public:
	RefNodeInfo(char* Filename,int Line):TagNodeInfo(kTagNodeRef){
		filename=Filename;
		line=Line;
	}
	char* filename;
	int line;
};

int LineSortFunction(Fl_Node* a, Fl_Node* b);
int NameSortFunction(Fl_Node* a, Fl_Node* b);

void TagsBrowser::SortByLine(void)
{
	mTree->sort_tree(LineSortFunction);
	mTree->redraw();
}

void TagsBrowser::SortByName(void)
{
	mTree->sort_tree(NameSortFunction);
	mTree->redraw();
}

int LineSortFunction(Fl_Node* a, Fl_Node* b) {
	TagNodeInfo* ia=(TagNodeInfo*) ((Fl_ToggleNode*)a)->data();
	TagNodeInfo* ib=(TagNodeInfo*) ((Fl_ToggleNode*)b)->data();
	
	if (ia->type==kTagNodeGlobal) return -1;
	if (ib->type==kTagNodeGlobal) return 1;
	if (ia->type==kTagNodeMembers) return -1;
	if (ib->type==kTagNodeMembers) return 1;
		
	if (ia->type==kTagNodeRef) {
		RefNodeInfo* ra=(RefNodeInfo*) ia;
		RefNodeInfo* rb=(RefNodeInfo*) ib;
		int c=strcmp(rb->filename,ra->filename);
		if (c==0) {
			if (ra->line<rb->line) return -1;
			if (ra->line>rb->line) return 1;
			return 0;
		}
		return c;
	}
	if (ia->type>kTagNodeCategory) {
		if (ia->type>ib->type) return 1;
		if (ib->type>ia->type) return -1;
		return 0;
	}		
		
	return strcmp(((Fl_ToggleNode*)a)->label(), ((Fl_ToggleNode*)b)->label());
}

int NameSortFunction(Fl_Node* a, Fl_Node* b) {
	TagNodeInfo* ia=(TagNodeInfo*) ((Fl_ToggleNode*)a)->data();
	TagNodeInfo* ib=(TagNodeInfo*) ((Fl_ToggleNode*)b)->data();
	
	if (ia->type==kTagNodeGlobal) return -1;
	if (ib->type==kTagNodeGlobal) return 1;
	if (ia->type==kTagNodeMembers) return -1;
	if (ib->type==kTagNodeMembers) return 1;
	if (ia->type>kTagNodeCategory) {
		if (ia->type>ib->type) return 1;
		if (ib->type>ia->type) return -1;
		return 0;
	}		
	return strcmp(((Fl_ToggleNode*)a)->label(), ((Fl_ToggleNode*)b)->label());
}

Fl_ToggleNode* catnodes[8];
int ncats=0;

class TagTreeLoader
{
public:
	Fl_ToggleTree* mTree;
	char** mFileNames;
	int mFileNamesCount;		
	TagTreeLoader(Fl_ToggleTree* t)
	{
		mTree=t;
		mFileNames=0;
		mFileNamesCount=0;
	}

	Fl_ToggleNode* new_group(char* str,Fl_ToggleNode* cat,Fl_Pixmap* pixmap=folderSmall)
	{
		mTree->traverse_start(cat);

		Fl_ToggleNode* n=new Fl_ToggleNode(str,1,pixmap);
		GroupNodeInfo* i = new GroupNodeInfo;
		mTree->add_sub(n);
		mTree->close(n);
		n->data(i);
		if (strcmp(cat->label(),"class")==0) {
			i->methods=mTree->add_sub("Methods",1,boxesSmall);
			MethodsNodeInfo* i2 = new MethodsNodeInfo;
			i->methods->data(i2);
			mTree->close(i->methods);
			i->members=mTree->add_next("Members",1,boxesSmall);
			MembersNodeInfo* i3 = new MembersNodeInfo;
			i->members->data(i3);
			mTree->close(i->members);
		}else{
			i->members=n;
			i->methods=n;
		}
		
		return n;
	};

	void new_global(void)
	{
		if (get_cat("global")) return;
		Fl_ToggleNode* global=new_cat("global",globalSmall);
		GroupNodeInfo* i = new GroupNodeInfo;
		i->type=1;
		global->data(i);

		i->methods=mTree->add_sub("Functions",1,boxesSmall);
		MethodsNodeInfo* i2 = new MethodsNodeInfo;
		i->methods->data(i2);
		mTree->close(i->methods);
		i->members=mTree->add_next("Variables",1,boxesSmall);
		MembersNodeInfo* i3 = new MembersNodeInfo;
		i->members->data(i3);
		mTree->close(i->members);
	};
    // F.C add-on
        void new_sources(void)
	{
	    static GroupNodeInfo* i = new GroupNodeInfo;
		if (!get_cat("sources"))
		{
			Fl_ToggleNode* sources=new_cat("sources",folderSmall);
			i->type= kTagNodeSources;
			sources->data(i);
			i->methods=mTree->add_sub("Include",1,folderSmall);
			SourceNodeInfo* i2 = new SourceNodeInfo;
			i->methods->data(i2);
			mTree->close(i->methods);
			i->members=mTree->add_next("Source",1,folderSmall);
			SourceNodeInfo* i3 = new SourceNodeInfo;
			i->members->data(i3);
			mTree->close(i->members);
		}
		Fl_ToggleNode* node = NULL;
		for (int n=0; n<mFileNamesCount;n++)
		{
		    char * s = name_only(mFileNames[n]);
		    size_t l = strlen(s);
		    if (s && l>2)
		    {
			// is it an std include file ?
			if (strcmp(s+l-2,".H")==0 || strcmp(s+l-2,".h")==0 
			    || strcmp(s+l-3,".hh")==0 || strcmp(s+l-4,".hpp")==0 )
			    node = i->methods;
			else
			    node = i->members;
			mTree->traverse_start(node);
			Fl_ToggleNode* nn=new Fl_ToggleNode(s,0,fileSmall);
			mTree->add_sub(nn);
			RefNodeInfo* ni=new RefNodeInfo(mFileNames[n],0);
			nn->data(ni);
		    }
		}
	};

	Fl_ToggleNode* get_cat(char* catstr) {
		int cati;
		
		for (cati=0;cati<ncats;cati++) {
			if (strcmp(catstr,catnodes[cati]->label())==0) {
				return catnodes[cati];
			}
		}
		return 0;	
	}

	Fl_ToggleNode* new_cat(char* catstr,Fl_Pixmap* pixmap=folderSmall)
	{
		Fl_ToggleNode* cat;
		if (ncats) {
			mTree->traverse_start(catnodes[ncats-1]);
		}else{
			mTree->traverse_start();
		}
		cat=new Fl_ToggleNode(catstr,1,pixmap);
		mTree->add_next(cat);
		GroupNodeInfo* i = new GroupNodeInfo;
		i->type=kTagNodeCategory+ncats;
		cat->data(i);
		mTree->close(cat);		
		catnodes[ncats++]=cat;
	
		return cat;
	}
	
	Fl_ToggleNode* get_or_new_cat(char* catstr)
	{
		Fl_ToggleNode* cat=get_cat(catstr);
		if (cat) return cat;
		return new_cat(catstr);	
	}
	
	void load(char** files,int nfiles)
	{
		new_global();
	
		get_or_new_cat("class");
		get_or_new_cat("struct");
		get_or_new_cat("enum");
		
		char* cmd=0;
		int i;
		if (nfiles) {
			int l=0;
			for (i=0;i<nfiles;i++) {
				l+=1+strlen(files[i]);
			}
			char* cmdctags="ctags -n -i cdefgmpstuv -f -";
			l+=strlen(cmdctags);
			cmd=new char[l+1];
			strcpy(cmd,cmdctags);
			for (i=0;i<nfiles;i++) {
				strcat(cmd," ");
				strcat(cmd,files[i]);
			}
		}
		
		char line[256];
		FILE* ctagsout=0;
		
		if (nfiles) {
#ifdef VERBOSE
			printf("running %s\n",cmd);
#endif
			ctagsout=popen(cmd,"r");			
			delete cmd;
		}else{
#ifdef VERBOSE
			printf("opening file \"tags\"");
#endif
			ctagsout=fopen("tags","r");			
		}
		if (ctagsout) while (fgets(line,255,ctagsout)) {
			char* name=strtok(line,"\t");
			char* file=strtok(0,"\t");
			char* line=strtok(0,"\t");
			char* type=strtok(0,"\t");
			char* tmp1=strtok(0,"\n");
			char* owner=0;
			Fl_ToggleNode* cat=catnodes[0];
						
			{
			    char* catstr=strtok(tmp1,":");
			    owner=strtok(0,"\0");
			    if (catstr) {
				cat=get_or_new_cat(catstr);
				if (owner==0) {
				    owner="unknown";
				}
			    }									
			}

			if (type && (type[0]=='f' || type[0]=='m' || type[0]=='e')) {
			    Fl_ToggleNode * curr=0;
				
			    if (owner) {
				mTree->traverse_start(cat);
				curr = mTree->traverse_forward();
	  			while (curr) {
				    TagNodeInfo* ni=(TagNodeInfo*) curr->data();
				    if (ni->type==2 || ni->type==(kTagNodeGroup|kTagNodeDeleting)) {
					if (strcmp(curr->label(),owner) == 0) {
					    ni->type=kTagNodeGroup;
					    break;
					}
				    }
				    if (ni->type==kTagNodeCategory) {
					curr=0;
				    }else{
		    			curr = mTree->traverse_forward();
				    }
				}
			    }
			    if (curr == 0) {
				if (owner) {
				    curr = new_group(owner,cat,boxSmall);
				}else {
				    curr = catnodes[0];
				}
			    }
			    if (type[0]=='f') {
				new_method((GroupNodeInfo*) curr->data(),
					   name,file,atoi(line));
			    }else{
				new_member((GroupNodeInfo*) curr->data(),
					   name,file,atoi(line));
			    }
			}	

		}
		// add sources at the end...
		new_sources();
	}
	
    char * name_only(char *filename)
    {
	if (filename )
	{
	    char *s;
	    size_t l = strlen(filename);
	    if (l)
		for (s=filename+l-1;s>=filename;s--) 
		    if (*s == '/' || *s =='\\' || *s==':') return (s+1);
	}
	return filename;
    }

    char* mFileNamestr(char* filename)
	{
	    for (int i=0;i<mFileNamesCount;i++) {
		if (strcmp(filename,mFileNames[i])==0) return mFileNames[i];
		}
		
		mFileNames=(char**) realloc(mFileNames,(mFileNamesCount+1)*sizeof(char*));
		filename=strdup(filename);
		mFileNames[mFileNamesCount++]=filename;
		return filename;
	}
	void new_method(GroupNodeInfo* ci,char* str,char* filename,int line)
	{
		mTree->traverse_start(ci->methods);
		Fl_ToggleNode* n=new Fl_ToggleNode(str,0,functionSmall);
		mTree->add_sub(n);
		RefNodeInfo* ni=new RefNodeInfo(mFileNamestr(filename),line);
		n->data(ni);
	}
	void new_member(GroupNodeInfo* ci,char* str,char* filename,int line)
	{
		mTree->traverse_start(ci->members);
		Fl_ToggleNode* n=new Fl_ToggleNode(str,0,pointSmall);
		mTree->add_sub(n);
		RefNodeInfo* ni=new RefNodeInfo(mFileNamestr(filename),line);
		n->data(ni);
	}
};

void TagsBrowser::cb_tree(Fl_Widget* t) {
	((TagsBrowser*) (t->parent()->parent()->user_data()))->TreeSelect();
}

void TagsBrowser::TreeSelect(void) {
  if (mTree->selected()) {
	TagNodeInfo* ti=(TagNodeInfo*) mTree->selected()->data();
	if (ti->type==kTagNodeRef) {
		RefNodeInfo* ri=(RefNodeInfo*) ti;
        char cmd[1024];
		int i=0,j=0;
		while (mInvokeEditorStr[i]) {
			if (mInvokeEditorStr[i]=='%') {
				i++;
				switch (mInvokeEditorStr[i]) {
				case '%':
					cmd[j]='%';
					break;
				case 'f':
					cmd[j]='\0';
					strcat(&cmd[j],ri->filename);
					j+=strlen(ri->filename);
					break;
				case 'n':
					cmd[j]='\0';
					{
						char tmp[8];
						sprintf(tmp,"%d",ri->line);
						strcat(&cmd[j],tmp);
						j+=strlen(tmp);
					}
					break;
				default:
					cmd[j++]='%';
					cmd[j]=mInvokeEditorStr[i];
				}
			}else{
				cmd[j]=mInvokeEditorStr[i];
				j++;
			}		
			i++;
		}
		cmd[j]='\0';
#ifdef VERBOSE
		printf("calling editor: %s %s\n",mInvokeEditorStr,cmd);
#endif
		system(cmd);
	}
  }
}

#define kByName 0

char prefsfilename[1024];
/*
void TagsBrowser::Init(void)
{
	sprintf(prefsfilename,"%s/.fltags",getenv("HOME"));

	mSettingsParser.Add(
		"make_ctags",&mMakeCTagsStr,
			strdup("make ctags"));
	mSettingsParser.Add(
		"invoke_ctags",&mInvokeCTagsStr,
			strdup("ctags -n -i cdefgmpstuv -f -"));
	mSettingsParser.Add(
		"invoke_editor",&mInvokeEditorStr,
			strdup("nc -noask +%n %f"));
	mSettingsParser.Add(
		"sort_type",&mSortType,0);

	mSettingsParser.Parse(prefsfilename);

	InitTree();
}
*/

void TagsBrowser::Load(char** files,int nfiles)
{
	mFiles = files;
	mFilesCount = nfiles;
	
	_Load(mFiles,mFilesCount);
		
}

void TagsBrowser::_Load(char** files,int nfiles)
{
    TagTreeLoader tl(mTree);
  
    tl.load(files,nfiles);

    if (mSortType==kByName)
	SortByName();
    else
	SortByLine();	
}

void TagsBrowser::MakeTagsAndReload(void)
{
	system(mMakeCTagsStr);
	
	Reload(1);
}

void TagsBrowser::Reload(int force_load_tags)
{
	Fl_ToggleNode *curr;
	curr=(Fl_ToggleNode*) mTree->traverse_start();
	
	while (curr) {
		TagNodeInfo* ti=(TagNodeInfo*) curr->data();
		if (ti->type==kTagNodeGroup) ti->type=kTagNodeGroup|kTagNodeDeleting;
		else if (ti->type==kTagNodeRef) mTree->remove(curr);
		curr=(Fl_ToggleNode*) mTree->traverse_forward();
	}

	if (force_load_tags) {
		_Load(0,0);
	}else{
		_Load(mFiles,mFilesCount);
	}

	curr=(Fl_ToggleNode*) mTree->traverse_start();
	
	while (curr) {
		TagNodeInfo* ti=(TagNodeInfo*) curr->data();

		if (ti->type==(kTagNodeGroup|kTagNodeDeleting)){
			mTree->remove(curr);
		}
		curr=(Fl_ToggleNode*) mTree->traverse_forward();
	}

	mTree->unselect();

	scroll->redraw();
}

void TagsBrowser::InitTree(void)
{
	mTree = new Fl_ToggleTree(
		scroll->x()+2,
		scroll->y()+2,
		scroll->w()-scroll->scrollbar.w()-3, 0); 
	scroll->add(mTree);

	Fl_ToggleTreeResizer* ttr = new Fl_ToggleTreeResizer(scroll,mTree);

	((Fl_Group*)(scroll->parent()))->add(ttr);

	mTree->callback(cb_tree);
}

void TagsBrowser::EditPreferences(void)
{
}

void TagsBrowser::SavePreferences(void)
{
}

void TagsBrowser::About(void)
{
}

void TagsBrowser::Search(void)
{
}

int CalcScore(const char* str1,const char* str2)
{
	int i=0;
	int score=0;
	int a=20;
	while (str1[i] && str2[i]) {
		int dif=str1[i]-str2[i];
		if (dif<0) dif=-dif;
		score+=((127-dif)<<a);
		a-=4;
		i++;
	} 
	return score;
}

void TagsBrowser::Init(int argc,char** argv) 
{
	mMakeCTagsStr="make tags";
	mInvokeCTagsStr="ctags -n -i cdefgmpstuv -f -";
	mInvokeEditorStr="nc -noask +%n %f";

	if (fileSmall==0) {
	  fileSmall = new Fl_Pixmap(file_small);
	  folderSmall = new Fl_Pixmap(folder_small);
	  pointSmall = new Fl_Pixmap(point_small);
	  functionSmall = new Fl_Pixmap(function_small);
	 	boxSmall = new Fl_Pixmap(box_small);
	  boxesSmall = new Fl_Pixmap(boxes_small);
	  globalSmall = new Fl_Pixmap(global_small);
	}

	InitTree();
	
//  gui.Load(argv+1,argc-1);

//  gui.Run();
}

