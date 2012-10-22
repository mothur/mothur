//
//  classifysharedcommand.h
//  Mothur
//
//  Created by Abu Zaher Md. Faridee on 8/13/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__classifysharedcommand__
#define __Mothur__classifysharedcommand__

#include "command.hpp"
#include "inputdata.h"

class ClassifySharedCommand : public Command {
public:
  ClassifySharedCommand();
  ClassifySharedCommand(string);
  ~ClassifySharedCommand() {};
  
  vector<string> setParameters();
  string getCommandName()			{ return "classify.shared";     }
   string getCommandCategory()		{ return "OTU-Based Approaches";		}
  string getOutputFileNameTag(string, string);
  string getHelpString();
  string getCitation() { return "http://www.mothur.org/wiki/Classify.shared\n"; }
  string getDescription()		{ return "description"; }
  int execute();
  
  void help() { m->mothurOut(getHelpString()); }

private:
    bool abort;
    string outputDir;
    vector<string> outputNames, Groups;
  
    string sharedfile, designfile, otupersplit, splitcriteria;
    set<string> labels;
    bool allLines;
  
    int processors;
    bool useTiming;

    GroupMap designMap;
  
    int numDecisionTrees;
    string treeSplitCriterion, optimumFeatureSubsetSelectionCriteria;
    bool doPruning, discardHighErrorTrees;
    double pruneAggressiveness, highErrorTreeDiscardThreshold, featureStandardDeviationThreshold;
    
    void processSharedAndDesignData(vector<SharedRAbundVector*> lookup);
};

#endif /* defined(__Mothur__classifysharedcommand__) */
