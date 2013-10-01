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

class ClassifySvmSharedCommand : public Command, ExternalSvmTrainingInterruption {
public:
  ClassifySvmSharedCommand();
  ClassifySvmSharedCommand(string);
  ~ClassifySvmSharedCommand() throw() {};
  
  vector<string> setParameters();
  string getCommandName()			{ return "classifysvm.shared";     }
  string getCommandCategory()		{ return "OTU-Based Approaches";		}  
  string getHelpString();	
  string getOutputPattern(string);
  string getCitation() { return "http://www.mothur.org/wiki/ClassifySvm.shared\n"; }
  string getDescription()		{ return "implements the support vector machine machine learning algorithm to identify OTUs that can be used to differentiate between various groups of samples"; }
  int execute();
  
  void help() { m->mothurOut(getHelpString()); }

  void readSharedAndDesignFiles(const std::string&, const std::string&, LabeledObservationVector&, FeatureVector&);
  void readSharedRAbundVectors(vector<SharedRAbundVector*>&, GroupMap&, LabeledObservationVector&, FeatureVector&);

  bool interruptTraining() { return m->control_pressed; }

  std::vector<double>& getSmocList() { return smocList; }
  const KernelParameterRangeMap& getKernelParameterRangeMap() { return kernelParameterRangeMap; }

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

    int evaluationFoldCount;
    int trainingFoldCount;
    std::vector<double> smocList;
    KernelParameterRangeMap kernelParameterRangeMap;

    //int numDecisionTrees;
    //string treeSplitCriterion, optimumFeatureSubsetSelectionCriteria;
    //bool doPruning, discardHighErrorTrees;
    //double pruneAggressiveness, highErrorTreeDiscardThreshold, featureStandardDeviationThreshold;

    void processSharedAndDesignData(vector<SharedRAbundVector*> lookup);
    void trainSharedAndDesignData(vector<SharedRAbundVector*> lookup);

    void getParameterValue(int& target, std::string pstring, int defaultvalue) {
        if (pstring == "not found" or pstring == "") {
            target = defaultvalue;
        }
        else {
            m->mothurConvert(pstring, target);
        }
    }


};

#endif /* defined(__Mothur__classifysvmsharedcommand__) */
