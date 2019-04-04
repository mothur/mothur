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

class EstimatorSingleCommand : public Command {
public:
    EstimatorSingleCommand(string);
    EstimatorSingleCommand();
    ~EstimatorSingleCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "estimator.single";            }
    string getCommandCategory()		{ return "OTU-Based Approaches";		}
    
    string getHelpString();
    string getOutputPattern(string);
    string getCitation() { return "http://www.mothur.org/wiki/Estimator.single"; }
    string getDescription()		{ return "This command implements the diversity estimators from https://github.com/chrisquince/DiversityEstimates"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    
    bool abort, allLines;
    string label, calc, outputDir, sharedfile, listfile, rabundfile, sabundfile, format, inputfile;
    float freq;
    vector<string> outputNames;
    set<string> labels; //holds labels to be used
    vector<string>  Estimators, groups;
    map<string, string> variables;
    
    vector<string> parseSharedFile(string);
    
    int process(SAbundVector*&);
    int runErarefaction(SAbundVector*&);
    
};


#endif /* diversityestimatorcommand_hpp */
