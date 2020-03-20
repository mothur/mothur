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
#include "designmap.h"

class ClassifySvmSharedCommand : public Command {
public:
  
  ClassifySvmSharedCommand(string);
  ~ClassifySvmSharedCommand() {};
  
  vector<string> setParameters();
  string getCommandName()			{ return "classify.svm";     }
  string getCommandCategory()		{ return "OTU-Based Approaches";		}  
  string getHelpString();	
  string getOutputPattern(string);
  string getCitation()              { return "http://www.mothur.org/wiki/Classify.svm\n"; }
  string getDescription()		    { return "implements the support vector machine machine learning algorithm to identify OTUs that can be used to differentiate between various groups of samples"; }
  int execute();
  
  void help() { m->mothurOut(getHelpString()); }

  void readSharedAndDesignFiles(const string&, const string&, LabeledObservationVector&, FeatureVector&);
  void readSharedRAbundVectors(vector<SharedRAbundVector*>&, DesignMap&, LabeledObservationVector&, FeatureVector&, vector<string>);

  vector<double>& getSmocList() { return smocList; }
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

    DesignMap designMap;
    
    // mode is either "rfe" or "classify"
    string mode;

    int evaluationFoldCount;
    int trainingFoldCount;
    vector<double> smocList;
    KernelParameterRangeMap kernelParameterRangeMap;

    string transformName;

    int verbosity;

    double stdthreshold;


    void processSharedAndDesignData(vector<SharedRAbundVector*> lookup, vector<string>);
    void trainSharedAndDesignData(vector<SharedRAbundVector*> lookup, vector<string>);

    void getParameterValue(int& target, string pstring, int defaultvalue) {
        if (pstring == "not found" or pstring == "") {
            target = defaultvalue;
        }
        else {
            util.mothurConvert(pstring, target);
        }
    }


};

#endif /* defined(__Mothur__classifysvmsharedcommand__) */
