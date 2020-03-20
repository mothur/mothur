//
//  primerdesigncommand.h
//  Mothur
//
//  Created by Sarah Westcott on 1/18/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef Mothur_primerdesigncommand_h
#define Mothur_primerdesigncommand_h

#include "command.hpp"
#include "listvector.hpp"
#include "inputdata.h"
#include "sequence.hpp"
#include "alignment.hpp"
#include "needlemanoverlap.hpp"

/**************************************************************************************************/

class PrimerDesignCommand : public Command {
public:
    PrimerDesignCommand(string);
    ~PrimerDesignCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "primer.design";		}
    string getCommandCategory()		{ return "OTU-Based Approaches";		} 
    
    string getOutputPattern(string);
	string getHelpString();	
    string getCitation() { return "http://www.mothur.org/wiki/Primer.design"; }
    string getDescription()		{ return "identify sequence fragments that are specific to particular OTUs"; }
    
    int execute(); 
    void help() { m->mothurOut(getHelpString()); }	
    
private:
    
    bool abort, allLines, large;
    int cutoff, pdiffs, length, processors, alignedLength;
    string outputDir, listfile, otulabel, namefile, countfile, fastafile, label;
    double minTM, maxTM;
    vector<string> outputNames;

    char getBase(vector<unsigned int> counts, int size);
    ListVector* getListVector();
    set<string> getPrimer(Sequence);
    int findMeltingPoint(string primer, double&, double&);
    set<int> createProcesses(string, vector<double>&, vector<double>&, set<string>&, vector<Sequence>&, int, vector<string>&);
    map<string, int> getSequenceBinAssignments(ListVector* list, map<string, int>& nameMap);
    vector<Sequence> createProcessesConSeqs(map<string, int>&, map<string, int>&, vector<string>&);
    int findIndex(string binLabel, vector<string> binLabels);
    
};

/**************************************************************************************************/

#endif
