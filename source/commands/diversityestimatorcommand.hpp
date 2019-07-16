//
//  diversityestimatorcommand.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/4/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef diversityestimatorcommand_hpp
#define diversityestimatorcommand_hpp

#include "command.hpp"
#include "inputdata.h"
#include "validcalculator.h"


//******************************************************

class EstimatorSingleCommand : public Command {
public:
    EstimatorSingleCommand(string);
    EstimatorSingleCommand();
    ~EstimatorSingleCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "estimator.single";            }
    string getCommandCategory()		{ return "OTU-Based Approaches";		}
    
    string getHelpString();
    string getCommonQuestions();
    string getOutputPattern(string);
    string getCitation() { return "http://www.mothur.org/wiki/Estimator.single"; }
    string getDescription()		{ return "This command implements the diversity estimators from https://github.com/chrisquince/DiversityEstimates"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    
    bool abort, allLines, burnSampleSet, burnSet, createSampling, itersSet;
    string label, calc, outputDir, sharedfile, listfile, rabundfile, sabundfile, format, inputfile, samplefile;
    double freq, sigmaAlpha, sigmaBeta, sigmaS, sigmaN, coverage;
    int iters, burn, burnSample, fitIters;
    vector<string> outputNames;
    set<string> labels; //holds labels to be used
    vector<string> groups, rarefactCalcs, abundCalcs, smallBurn;
    map<string, vector<mcmcSample> > sampling;
    map<string, vector<mcmcSample> > ::iterator it;
    set<string> samplingCalcs;
    map<string, string> calcToSamplingCalc;
    map<string, string> ::iterator itCalcSample;
    
    int fillSampling(int, int, bool filldNu=false);
    int processSingleSample();
    int processSharedFile();
    int processShared(SharedRAbundVectors*& shared, vector<ofstream*>& out, string fileRoot);
    int processSingle(SAbundVector*&, string, vector<ofstream*>&, string);
    
    int runRarefactCalcs(int numSeqs, string groupName, ofstream& out);
    vector<string> runSamplingCalcs(SAbundVector*&, string);
    vector<double> runAbundCalcs(SAbundVector*&, string groupName);
    
};
//*******************************************************

#endif /* diversityestimatorcommand_hpp */
