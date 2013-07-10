//
//  classifysvmsharedcommand.h
//  Mothur
//
//  Created by Joshua Lynch on 6/28/2013.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//
//  This class is based on ClassifySharedCommand
//

#ifndef __Mothur__classifysvmsharedcommand__
#define __Mothur__classifysvmsharedcommand__

#include "command.hpp"
#include "inputdata.h"
#include "svm.hpp"

class ClassifySvmSharedCommand : public Command {
public:
  ClassifySvmSharedCommand();
  ClassifySvmSharedCommand(string);
  ~ClassifySvmSharedCommand() {};
  
  vector<string> setParameters();
  string getCommandName()			{ return "classifysvm.shared";     }
  string getCommandCategory()		{ return "OTU-Based Approaches";		}  
  string getHelpString();	
  string getOutputPattern(string);
  string getCitation() { return "http://www.mothur.org/wiki/ClassifySvm.shared\n"; }
  string getDescription()		{ return "implements the support vector machine machine learning algorithm to identify OTUs that can be used to differentiate between various groups of samples"; }
  int execute();
  
  void help() { m->mothurOut(getHelpString()); }

  static void readSharedAndDesignFiles(const std::string&, const std::string&, LabeledObservationVector&);

private:
    bool abort;
    string outputDir;
    vector<string> outputNames, Groups;
  
    string sharedfile, designfile;
    set<string> labels;
    bool allLines;
  
    int processors;
    bool useTiming;

    GroupMap designMap;
  
    //int numDecisionTrees;
    //string treeSplitCriterion, optimumFeatureSubsetSelectionCriteria;
    //bool doPruning, discardHighErrorTrees;
    //double pruneAggressiveness, highErrorTreeDiscardThreshold, featureStandardDeviationThreshold;
    
    void processSharedAndDesignData(vector<SharedRAbundVector*> lookup);
};

#endif /* defined(__Mothur__classifysvmsharedcommand__) */
