//
//  sparcccommand.h
//  Mothur
//
//  Created by SarahsWork on 5/10/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef Mothur_sparcccommand_h
#define Mothur_sparcccommand_h

#include "command.hpp"
#include "inputdata.h"
#include "calcsparcc.h"

/**************************************************************************************************/

class SparccCommand : public Command {
public:
    SparccCommand(string);
    SparccCommand();
    ~SparccCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "sparcc";			}
    string getCommandCategory()		{ return "OTU-Based Approaches";		}
    
    string getOutputPattern(string);
    //commmand category choices: Sequence Processing, OTU-Based Approaches, Hypothesis Testing, Phylotype Analysis, General, Clustering and Hidden
	string getHelpString();
    string getCitation() { return "Friedman J, Alm EJ (2012) Inferring Correlation Networks from Genomic Survey Data. PLoS Comput Biol 8(9): e1002687. doi:10.1371/journal.pcbi.1002687 http://www.mothur.org/wiki/Sparcc"; }
    string getDescription()		{ return "Calculates correlations between OTUs using a method that is insensitive to the use of relative abundance data"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    bool abort, allLines;
    string outputDir, sharedfile, normalizeMethod;
    int numSamplings, maxIterations, numPermutations, processors;
    set<string> labels;
    vector<string> Groups;
    vector<string> outputNames;
    
    int process(SharedRAbundVectors*&);
    vector<vector<float> > createProcesses(vector<vector<float> >&, vector<vector<float> >&);
    //vector<vector<float> > driver(vector<vector<float> >&, vector<vector<float> >&, int);
    //vector<vector<float> > shuffleSharedVector(vector<vector<float> >&);
};

/**************************************************************************************************/

#endif
