#ifndef SEQSUMMARYCOMMAND_H
#define SEQSUMMARYCOMMAND_H

/*
 *  seqcoordcommand.h
 *  Mothur
 *
 *  Created by Pat Schloss on 5/30/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "sequence.hpp"

/**************************************************************************************************/

class SeqSummaryCommand : public Command {
public:
	SeqSummaryCommand(string);
	SeqSummaryCommand();
	~SeqSummaryCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "summary.seqs";			}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Summary.seqs"; }
	string getDescription()		{ return "summarize the quality of sequences in an unaligned or aligned fasta file"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }		
private:
	bool abort;
	string fastafile, outputDir, namefile, countfile;
	int processors;
	vector<string> outputNames;
};

#endif

/**************************************************************************************************/


