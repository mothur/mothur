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

class SortSeqsCommand : public Command {
	
public:
	
    SortSeqsCommand(string);	
    SortSeqsCommand();
    ~SortSeqsCommand(){}
	
    vector<string> setParameters();
    string getCommandName()			{ return "sort.seqs";				}
    string getCommandCategory()		{ return "Sequence Processing";		}
    string getHelpString();	
    string getCitation() { return "http://www.mothur.org/wiki/Sort.seqs"; }
    string getDescription()		{ return "puts sequences from a fasta, name, group, quality or taxonomy file in the same order"; }
    
    int execute(); 
    void help() { m->mothurOut(getHelpString()); }	
	
	
private:
    map<string, int> names;
    string accnosfile, fastafile, namefile, groupfile, taxfile, qualfile, outputDir;
    bool abort, large;
    vector<string> outputNames;
    
    int readFasta();
    int readName();
    int readGroup();
    int readAccnos();
    int readTax();
    int readQual();
    
};

#endif


