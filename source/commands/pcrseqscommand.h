#ifndef Mothur_pcrseqscommand_h
#define Mothur_pcrseqscommand_h

//
//  pcrseqscommand.h
//  Mothur
//
//  Created by Sarah Westcott on 3/14/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//


#include "command.hpp"
#include "sequence.hpp"
#include "trimoligos.h"
#include "alignment.hpp"
#include "needlemanoverlap.hpp"
#include "counttable.h"
#include "oligos.h"
#include "removeseqscommand.h"

class PcrSeqsCommand : public Command {
public:
	PcrSeqsCommand(string);
	PcrSeqsCommand();
	~PcrSeqsCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "pcr.seqs";	}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Pcr.seqs"; }
	string getDescription()		{ return "pcr.seqs"; }
    
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
    bool abort, keepprimer, keepdots, fileAligned, pairedOligos;
	string fastafile, oligosfile, taxfile, groupfile, namefile, countfile, ecolifile, outputDir, nomatch;
	int start, end, processors, length, pdiffs, rdiffs;
    vector<string> outputNames;
    
    int writeAccnos(set<string>, string);
    Sequence readEcoli();
	long long createProcesses(string, string, string, set<string>&);
    int adjustDots(string goodFasta, map<string, vector<int> > locations, int pstart, int pend, bool, bool);
};

/**************************************************************************************************/

#endif
