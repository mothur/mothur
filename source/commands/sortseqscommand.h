#ifndef Mothur_sortseqscommand_h
#define Mothur_sortseqscommand_h


//
//  sortseqscommand.h
//  Mothur
//
//  Created by Sarah Westcott on 2/3/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//



#include "command.hpp"
#include "counttable.h"

class SortSeqsCommand : public Command {
	
public:
	
    SortSeqsCommand(string);	
    ~SortSeqsCommand(){}
	
    vector<string> setParameters();
    string getCommandName()			{ return "sort.seqs";				}
    string getCommandCategory()		{ return "Sequence Processing";		}
    
	string getHelpString();	
    string getOutputPattern(string);	
    string getCitation() { return "http://www.mothur.org/wiki/Sort.seqs"; }
    string getDescription()		{ return "puts sequences from a fasta, name, group, quality, flow or taxonomy file in the same order"; }
    
    int execute(); 
    void help() { m->mothurOut(getHelpString()); }	
	
	
private:
    map<string, int> names;
    string accnosfile, fastafile, namefile, taxfile, qualfile, flowfile;
    bool abort, large;
    vector<string> outputNames;
    
    int readFasta();
    int readFlow();
    int readName();
    int readTax();
    int readQual();
    
};

#endif


