#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkUnaryFunctorImageFilter.h>
#include <string>
#include <vector>
#include <set>
#include <list>
#include <map>
#include <iostream>
#include <algorithm>

// Regular expression support
#include <itksys/RegularExpression.hxx>

using namespace std;


static const unsigned int slice_dim = 2;

// Mapper functions for thresholding body, head, tail
template<class TInput, class TOutput> class PickFunctor 
{
public:
  TInput x;
  TOutput operator()(const TInput &in)
    { if(in == x) return 1; else return 0; }
  bool operator == (const PickFunctor<TInput, TOutput> &dummy)
    { return dummy.x == x; }
  bool operator != (const PickFunctor<TInput, TOutput> &dummy)
    { return dummy.x != x; }
};

class Exc : public exception {
protected:
  string m_SimpleMessage;

public:
  Exc();
  Exc(const char *message, ...) : exception()
    {
    char buffer[1024];
    va_list args;
    va_start(args, message);
    vsprintf(buffer,message,args);
    va_end (args);
    m_SimpleMessage = buffer;
    }

  virtual ~Exc() throw() {};

  virtual const char * what() const throw() { return m_SimpleMessage.c_str(); }

  operator const char *() { return this->what(); }
};

int usage()
{
  printf(
    "subfield_slice_rules: generate exclusion priors from segmentation and rules\n"
    "usage:\n"
    "  subfield_slice_rules input_seg.nii rules.txt output_pattern.nii\n"
    "or\n"
    "  subfield_slice_rules --check-rules rules.txt\n"
    "parameters:\n"
    "  input_seg.nii        Input multi-label segmentation\n"
    "  rules.txt            Text file containing rules\n"
    "  output_pattern.nii   Output pattern (printf format, e.g. out_%%02d.nii.gz\n"
    "rules file format:\n"
    "The rules file consists of lines of text, each being one of three possible\n"
    "declarations: CLASS, GROUP, RULE, and RANGE_RULE.\n"
    "\n"
    "CLASS declaration\n"
    "=================\n"
    "\n"
    "The CLASS declaration states that a set of labels forms a named class:\n"
    "\n"
    "    CLASS class_name = labels\n"
    "\n"
    "where class_name is an identifier, and labels is a number or a comma separated\n"
    "list of numbers. For example, \n"
    "\n"
    "    CLASS body = 1,2,3,4\n"
    "\n"
    "Each slice will be assigned to one of the classes specified, based on the\n"
    "prevalence of voxels with that class in the slice. If there are no voxels of any\n"
    "class in a slice, the slice is assigned the NULL class.\n"
    "\n"
    "GROUP declaration\n"
    "=================\n"
    "\n"
    "The GROUP declaration states that a list of classes are exclusive, i.e., a slice\n"
    "can only belong to one of the classes in the group. There can be multiple groups\n"
    "\n"
    "    GROUP class_names [BACKGROUND fraction]\n"
    "\n"
    "where class_names is a comma-separated list. Within each group, slices will be\n"
    "assigned to one of the classes. However, putting classes in a group does\n"
    "not automatically set up any exclusion rules -- these rules must be supplied\n"
    "separately (see below).\n"
    "\n"
    "The optional BACKGROUND keyword followed by the fraction (between 0.0 and 1.0)\n"
    "specified that some slices in the group may be assigned the special BACKGROUND\n"
    "label if the number of foreground voxels in that slice is smaller than a certain\n"
    "threshold. Foreground voxels are voxels assigned any of the classes in the group, \n"
    "and the threshold is calculated as the fraction times the median number of foreground\n"
    "voxels for the group (slices that have zero foreground voxels are not used in the\n"
    "background calculation). \n"
    "\n"
    "Using the background keyword makes it possible to mask out slices where there are\n"
    "very few straddler foreground voxels. This is useful when the manual segmentation is\n"
    "restricted to somewhat arbitrary slice boundaries.\n"
    "\n"
    "RULE declaration\n"
    "================\n"
    "\n"
    "The RULE declaration states that if a slice belongs to a certain class, then\n"
    "voxels in that slice are not allowed to have certain labels:\n"
    "\n"
    "    RULE class_name EXCLUDES labels\n"
    "\n"
    "where class_name must be declared in one of the CLASS declarations and labels is\n"
    "a single label, single class, or a comma-separated list of labels and classes.\n"
    "For example:\n"
    "\n"
    "    RULE body EXCLUDES head, tail, 8\n"
    "\n"
    "means that voxels in slices belonging to the body class are not allowed to have\n"
    "labels in classes head, tail, and the label 8.\n"
    "\n"
    "RANGE_RULE declaration\n"
    "======================\n"
    "\n"
    "The RANGE_RULE declaration applies exclusions based on the position of the slice\n"
    "relative to the extents of a class. It has the following format:\n"
    "\n"
    "    RANGE_RULE range_start TO range_end EXCLUDES labels\n"
    "\n"
    "Where range_start and range_end are specifications in the form. NO SPACES!\n"
    "\n"
    "    <class_name:<FIRST|LAST>[<+|->offset]|END>\n"
    "\n"
    "For example:\n"
    "\n"
    "    RANGE_RULE body:FIRST+1 TO body:LAST-1 EXCLUDES 9\n"
    "    RANGE_RULE END TO body:LAST EXCLUDES 8\n"
    "\n"
    "EXAMPLE OF EXCLUDING STRADDLERS\n"
    "===============================\n"
    "\n"
    "This example rule file demonstrates how to 'clean up' straddler slices in the \n"
    "extra-hippocampal structures.\n"
    "\n"
    "CLASS ehmtl 9,10,11,12\n"
    "GROUP ehmtl BACKGROUND 0.25\n"
    "RULE BACKGROUND excludes ehmtl\n"
    "\n"
    "A class for extrahippocampal structures is first established. Next, we specify that\n"
    "a slice is classified as ehmtl only if there are more than 0.25 * N voxels in that\n"
    "slice labeled as ehmtl, where N is the median number of ehmtl voxels per slice that\n"
    "includes the ehmtl label. For instance, if there are eight slices and the number of\n"
    "ehmtl voxels per slice are 0,0,10,100,120,95,20,0 then N=95, and the cutoff will be\n"
    "0.25 * 95. The slices with 10 and 20 voxels will be marked as background. Last the\n"
    "RULE specifies that any slice marked as background should exclude the ehmtl class.\n"
    "\n"
    "EXAMPLE OF A COMPLEX RULE FOR HEAD/BODY/TAIL SEGMENTATION\n"
    "=========================================================\n"
    "\n"
    "# Three classes into which we partition the slices\n"
    "CLASS body = 1,2,3,4\n"
    "CLASS head = 5\n"
    "CLASS tail = 6\n"
    "CLASS sub = 8\n"
    "CLASS ercp = 9\n"
    "CLASS erca = 10\n"
    "CLASS prc = 11,12\n"
    "\n"
    "# These classes form a single mutually exclusive group\n"
    "GROUP body, head, tail\n"
    "GROUP erca, ercp\n"
    "\n"
    "# Additional exclusion rules applied\n"
    "RULE head EXCLUDES body tail\n"
    "RULE body EXCLUDES head tail\n"
    "RULE tail EXCLUDES head body\n"
    "\n"
    "# Additional range rules for ERC and PHG\n"
    "RANGE_RULE END TO head:LAST-2 EXCLUDES ercp\n"
    "RANGE_RULE head:LAST+2 TO END EXCLUDES erca,ercp,prc\n"
    "RANGE_RULE head:LAST-1 TO head:LAST+1 EXCLUDES erca\n"
    " \n");

  return -1;
}

using namespace std;
using namespace itksys;

typedef set<int> LabelSet;
typedef LabelSet::iterator LabelSetIter;

typedef map<string, LabelSet> ClassMap;
typedef set<string> ClassGroup;

enum RelPos { NA, FIRST, LAST };

struct RangeLimit
{
  bool end;
  RelPos pos;
  string cls;
  int offset;
  RangeLimit() : end(true), pos(NA), offset(0) {}
};

struct Rule
{
  // What is excluded
  LabelSet rhs;

  // The class that is matched, or "" if a range rule
  string cls;

  // The lower and upper range specs
  RangeLimit rl, ru;
};

struct Group
{
  // List of classes 
  ClassGroup clsgrp;

  // Background level (for cutting off straddlers)
  double bkg;
};

struct RuleCollection
{
  // Set of classes in the rules
  ClassMap cls;

  // Set of groups formed by the classes
  list<Group> grp;

  // Set of rules
  list<Rule> rules;
};

list<string> ReadCommaSeparatedList(const char *text)
{
  list<string> matches;
  RegularExpression re(" *([a-zA-Z0-9]+) *,*");
  while(re.find(text))
    {
    matches.push_back(re.match(1));
    text += re.end();
    }
  return matches;
}

set<int> ReadCommaSeparatedIntegerSet(const char *text)
{
  set<int> matches;
  RegularExpression re(" *([0-9]+) *,*");
  while(re.find(text))
    {
    matches.insert(atoi(re.match(1).c_str()));
    text += re.end();
    }
  return matches;
}

set<int> ReadCommaSeparatedIntegerAndClassSet(const char *text, ClassMap &cls)
{
  set<int> matches;
  RegularExpression re(" *([a-zA-Z0-9]+) *,*");
  while(re.find(text))
    {
    string match = re.match(1);
    if(cls.find(match) != cls.end())
      {
      for(LabelSet::iterator it = cls[match].begin(); it != cls[match].end(); it++)
        matches.insert(*it);
      }
    else
      {
      int k = atoi(match.c_str());
      if(k == 0)
        throw Exc("Bad label expression %s in rule right hand side %s", match.c_str(), text);
      matches.insert(k);
      }
    text += re.end();
    }
  return matches;
}

RangeLimit ReadLimitDesc(const char *text)
{
  RegularExpression reLim("([a-zA-Z]+):(FIRST|LAST)([+-][0-9]+|)");
  RegularExpression reEnd("^END$");
  RangeLimit rl;

  if(reEnd.find(text))
    {
    rl.end = true;
    }
  else if(reLim.find(text))
    {
    rl.end = false;
    rl.cls = reLim.match(1);
    rl.pos = (reLim.match(2) == "FIRST" ? FIRST : LAST);
    rl.offset = atoi(reLim.match(3).c_str());  
    }
  else throw Exc("Unable to parse limit descriptor %s", text);
  return rl;
}


RuleCollection parseRules(const char *rulesFile)
{
  // The output collection
  RuleCollection rc;

  // Open the file
  FILE *f = fopen(rulesFile, "rt");
  if(!f)
    throw Exc("Unable to read rules file %s", rulesFile);

  // Read each line
  char buffer[4096];
  while(fgets(buffer, 4096, f))
    {
    // Get rid of the newline
    char *pendl = strchr(buffer,'\n');
    if(pendl) *pendl = 0;

    // Replace tabs with spaces to shorten regex
    replace(buffer, buffer+strlen(buffer), '\t', ' ');

    // Skip empty lines and comment lines
    if(RegularExpression("^ *$").find(buffer) || RegularExpression("^ *#").find(buffer))
      continue;

    // Recognizable commands
    RegularExpression reClass("^ *CLASS +([a-zA-Z]+) *= *(.*) *$");
    RegularExpression reRule("^ *RULE +([a-zA-Z]+) +EXCLUDES +(.*) *$");
    RegularExpression reRange("^ *RANGE_RULE *([a-zA-Z0-9:\\+\\-]+) +TO +([a-zA-Z0-9:\\+\\-]+) +EXCLUDES +(.*)$");
    RegularExpression reGroup("^ *GROUP +(.*)$");
    RegularExpression reGroupBkg("^ *GROUP +(.*) +BACKGROUND +(.*)$");
    
    // Check if this is a class line
    if(reClass.find(buffer))
      {
      // Found a class declaration
      string clsName = reClass.match(1);

      // Get the list of labels
      LabelSet labels = ReadCommaSeparatedIntegerSet(reClass.match(2).c_str());
      if(labels.size() == 0)
        throw Exc("Bad label specification %s", reClass.match(2).c_str());

      // Save to the list of classes
      rc.cls[clsName] = labels;
      }

    else if(reGroupBkg.find(buffer))
      {
      // Get the names
      list<string> classes = ReadCommaSeparatedList(reGroupBkg.match(1).c_str());
      if(classes.size() == 0)
        throw Exc("Bad class list specification %s", reGroupBkg.match(1).c_str());

      // Get the background key
      double background = atof(reGroupBkg.match(2).c_str());

      // Store the group spec
      Group grp;
      grp.bkg = background;
      for(list<string>::iterator it = classes.begin(); it!=classes.end();it++)
        {
        if(rc.cls.find(*it) == rc.cls.end())
          throw Exc("Bad class specification %s in GROUP command %s", it->c_str(), buffer);
        grp.clsgrp.insert(*it);
        }


      rc.grp.push_back(grp);
      }

    // Check if this is a group line
    else if(reGroup.find(buffer))
      {
      // Get the names
      list<string> classes = ReadCommaSeparatedList(reGroup.match(1).c_str());
      if(classes.size() == 0)
        throw Exc("Bad class list specification %s", reGroup.match(1).c_str());

      // Store the group spec
      Group grp;
      grp.bkg = 0.0;
      for(list<string>::iterator it = classes.begin(); it!=classes.end();it++)
        {
        if(rc.cls.find(*it) == rc.cls.end())
          throw Exc("Bad class specification %s in GROUP command %s", it->c_str(), buffer);
        grp.clsgrp.insert(*it);
        }
      rc.grp.push_back(grp);
      }

    // Check if this is a rule line
    else if(reRule.find(buffer))
      {
      Rule rule;
      rule.cls = reRule.match(1);
      rule.rhs = ReadCommaSeparatedIntegerAndClassSet(reRule.match(2).c_str(), rc.cls);
      rc.rules.push_back(rule);
      }

    // Check if this is a range rule line
    else if(reRange.find(buffer))
      {
      Rule rule;
      rule.rhs = ReadCommaSeparatedIntegerAndClassSet(reRange.match(3).c_str(), rc.cls);
      string lim1 = reRange.match(1), lim2 = reRange.match(2);
      
      // Read the range descriptor
      rule.rl = ReadLimitDesc(lim1.c_str());
      rule.ru = ReadLimitDesc(lim2.c_str());
      rc.rules.push_back(rule);
      }

    else 
      {
      throw Exc("Unable to parse command %s", buffer);
      }
    }
  
  fclose(f);
  return rc;
}

int main(int argc, char *argv[])
{
  typedef itk::Image<short, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  string fnRules, fnInput, fnOutPattern;
  bool checkOnly = false;

  if(argc == 3 && !strcmp(argv[1], "--check-rules"))
    {
    // Special case - checking rules
    checkOnly = true;
    fnRules = argv[2];
    }
  else if(argc == 4)
    {
    fnInput = argv[1];
    fnRules = argv[2];
    fnOutPattern = argv[3];
    }
  else return usage();

  // Parse the rules
  RuleCollection rc;
  try
    {
    rc = parseRules(fnRules.c_str());
    printf("Success parsing rules file. The rules follow:\n");
    }
  catch(Exc &exc)
    {
    cout << "EXCEPTION parsing rules: " << exc.what() << endl;
    return -1;
    }

  // Describe the rules
  printf("  [CLASSES]\n");
  for(ClassMap::iterator itmap = rc.cls.begin(); itmap!=rc.cls.end(); itmap++)
    {
    printf("  %10s   :", itmap->first.c_str());
    for(LabelSet::iterator itset = itmap->second.begin(); itset != itmap->second.end(); itset++)
      printf(" %d", *itset);
    printf("\n");
    }

  printf("  [GROUPS]\n");
  int igrp = 0;
  for(list<Group>::iterator it = rc.grp.begin(); it!=rc.grp.end(); it++)
    {
    printf("  %10d   :", ++igrp);
    for(ClassGroup::iterator itg = it->clsgrp.begin(); itg!=it->clsgrp.end();itg++)
      printf(" %s", itg->c_str());
    printf(" BACKGROUND %f", it->bkg);
    printf("\n");
    }

  printf("  [RULES]\n");
  int irule = 0;
  for(list<Rule>::iterator it = rc.rules.begin(); it!=rc.rules.end(); it++)
    {
    printf("  %10d   :", ++irule);
    Rule r = *it;
    if(r.cls.size())
      {
      printf(" %s EXCLUDES", r.cls.c_str());
      }
    else
      {
      printf(" RANGE ");
      if(r.rl.end)
        printf(" END ");
      else
        printf(" %s:%s%c%d ", r.rl.cls.c_str(), (r.rl.pos == FIRST ? "FIRST" : "LAST"), r.rl.offset < 0 ? '-' : '+', abs(r.rl.offset));
      printf("TO");
      if(r.ru.end)
        printf(" END ");
      else
        printf(" %s:%s%c%d ", r.ru.cls.c_str(), (r.ru.pos == FIRST ? "FIRST" : "LAST"), r.ru.offset < 0 ? '-' : '+', abs(r.ru.offset));
      printf("EXCLUDES");
      }
    for(LabelSet::iterator lit = r.rhs.begin(); lit != r.rhs.end(); lit++)
      printf(" %d", *lit);
    printf("\n");
    }

  // We now actually have the rules!
  if(checkOnly)
    {
    return 0;
    }
    
  // Read reference segmentation
  ReaderType::Pointer fReader = ReaderType::New();
  fReader->SetFileName(fnInput.c_str());
  try 
    {
    fReader->Update();
    }
  catch(exception &exc)
    {
    cerr << "Exception reading input image " << fnInput << ": " << exc.what() << endl;
    return -1;
    }
  ImageType::Pointer ref = fReader->GetOutput();

  // Get the union of all the classes covered by the exclusions
  LabelSet allex;
  for(list<Rule>::iterator it = rc.rules.begin(); it!=rc.rules.end(); it++)
    {
    Rule r = *it;
    for(LabelSetIter qt = r.rhs.begin(); qt != r.rhs.end(); qt++)
      {
      allex.insert(*qt);
      }
    }


  // The slicing dimension, by default 2. Could be specified from command line
  int slice_dim = 2;
  int ns = ref->GetBufferedRegion().GetSize()[slice_dim]; 

  // For each slice, we need to pick the dominant class in each group. For that we
  // count the number of voxels of each label in each slice. This is done using a 
  // simple list of maps
  typedef map<int, long> Histogram;
  vector<Histogram> sliceHist(ns);

  // Build the slice histograms
  typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
  for(IteratorType it(ref, ref->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    int slice = it.GetIndex()[slice_dim];
    int label = it.Get();
    Histogram::iterator qt = sliceHist[slice].find(label);
    if(qt == sliceHist[slice].end())
      sliceHist[slice][label] = 1;
    else
      qt->second++;
    }

  // A list of classes assigned to each slice. A slice can have multiple
  // classes if there are multiple groups
  vector< set<string> > sliceCls(ns);

  // We will compute the extent of each of the classes
  map<string, int> clsFirst, clsLast;

  // For each group, we will compute the median number of voxels per slice
  // assigned to that group. This is used for setting the background 
  std::vector<int> grpMedVox;

  for(list<Group>::iterator it = rc.grp.begin(); it != rc.grp.end(); it++)
    {
    std::vector<int> vox_count;
    for(int i = 0; i < ns; i++)
      {
      int slice_vox_count = 0;

      // Find the number of voxels in this slice that belong to one of the classes in the group
      for(ClassGroup::iterator qt = it->clsgrp.begin(); qt != it->clsgrp.end(); qt++)
        {
        LabelSet clab = rc.cls[*qt];
        for(LabelSetIter gt = clab.begin(); gt != clab.end(); gt++)
          {
          slice_vox_count += sliceHist[i][*gt];
          }
        }

      // If this number is greater than zero, include it for median computation
      if(slice_vox_count > 0)
        vox_count.push_back(slice_vox_count);
      }

    // Find the median value
    std::sort(vox_count.begin(), vox_count.end());
    int median = vox_count[vox_count.size() / 2];
    grpMedVox.push_back(median);
    }

  // Within each group, within each slice, find the dominant class
  for(int i = 0; i < ns; i++)
    {
    int igrp = 0;
    for(list<Group>::iterator it = rc.grp.begin(); it != rc.grp.end(); it++, igrp++)
      {

      // Find the dominant class for this group in this slice
      string bestClass = "BACKGROUND";
      int bestClassSize = (int) (grpMedVox[igrp] * it->bkg);
      for(ClassGroup::iterator qt = it->clsgrp.begin(); qt != it->clsgrp.end(); qt++)
        {
        // Set of labels in that class
        LabelSet clab = rc.cls[*qt];
        int myClassSize = 0;
        for(LabelSetIter gt = clab.begin(); gt != clab.end(); gt++)
          myClassSize += sliceHist[i][*gt];
        if(myClassSize > bestClassSize)
          {
          bestClassSize = myClassSize;
          bestClass = *qt;
          }
        }

      // Store this information
      if(bestClassSize > 0)
        {
        // Set the class for the slice
        sliceCls[i].insert(bestClass);

        // Update the extents of the class
        if(clsFirst.find(bestClass) == clsFirst.end())
          clsFirst[bestClass] = i;

        clsLast[bestClass] = i;
        }
      }
    }

  // Print the classification for all the slices
  printf("Slice Classification:\n");
  igrp = 0;
  for(list<Group>::iterator it = rc.grp.begin(); it != rc.grp.end(); it++)
    {
    printf("  Group %d\n", ++igrp);
    for(ClassGroup::iterator qt = it->clsgrp.begin(); qt != it->clsgrp.end(); qt++)
      {
      printf("    %10s: |", qt->c_str());
      for(int i = 0; i < ns; i++)
        printf(sliceCls[i].count(*qt) ? "+" : " ");
      printf("|\n");
      }
    }

  // Now, we can apply the rules for each slice and each excluded label. We just need
  // a boolean for each label as to whether it is excluded or not
  vector< set<int> > sliceExcl(ns);

  // Go over all the slices
  for(int i = 0; i < ns; i++)
    {
    // Apply all the rules
    for(list<Rule>::iterator it = rc.rules.begin(); it != rc.rules.end(); it++)
      {
      Rule rule = *it;
      bool met = false;

      // Normal exclusion rule?
      if(rule.cls.size() && sliceCls[i].count(rule.cls))
        {
        met = true;
        }

      // Range rule
      else if(rule.cls.size() == 0) 
        {
        // Check if we are above or equal to the lower bound
        int lpos = rule.rl.end ? 0 : 
          (rule.rl.pos == FIRST ? clsFirst[rule.rl.cls] : clsLast[rule.rl.cls]) + rule.rl.offset;

        int upos = rule.ru.end ? ns : 
          (rule.ru.pos == FIRST ? clsFirst[rule.ru.cls] : clsLast[rule.ru.cls]) + rule.ru.offset;

        met = (i >= lpos) && (i <= upos);
        }

      if(met)
        {
        // Apply all of the rule exclusions to this slice
        sliceExcl[i].insert(rule.rhs.begin(), rule.rhs.end());
        }
      }
    }

  // Print a map of slices and exclusions
  printf("Slice Exclusions:\n");
  for(LabelSetIter it = allex.begin(); it != allex.end(); it++)
    {
    printf("    %10d: |", *it);
    for(int i = 0; i < ns; i++)
      printf(sliceExcl[i].count(*it) ? "*" : " ");
    printf("|\n");
    }

  // Generate an image for each of the exclusion classes
  for(LabelSetIter it = allex.begin(); it!=allex.end(); ++it)
    {
    // Create the exclusion image
    ImageType::Pointer exim = ImageType::New();
    exim->CopyInformation(ref);
    exim->SetRegions(ref->GetLargestPossibleRegion());
    exim->Allocate();
    exim->FillBuffer(0);

    // Iterate over the exclusion image
    for(itk::ImageRegionIteratorWithIndex<ImageType> qt(exim, exim->GetBufferedRegion());
      !qt.IsAtEnd(); ++qt)
      {
      int s = (int) qt.GetIndex()[slice_dim];
      qt.Set(sliceExcl[s].count(*it));
      }

    // Write the exclusion image
    char fname[4096];
    sprintf(fname, fnOutPattern.c_str(), *it);
    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer fWriter = WriterType::New();
    fWriter->SetInput(exim);
    fWriter->SetFileName(fname);
    fWriter->Update();
    }
}
