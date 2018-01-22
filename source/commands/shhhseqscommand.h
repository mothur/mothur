#ifndef SHHHSEQSCOMMAND_H
#define SHHHSEQSCOMMAND_H

/*
 *  shhhseqscommand.h
 *  Mothur
 *
 *  Created by westcott on 11/8/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "myseqdist.h"
#include "seqnoise.h"
#include "sequenceparser.h"
#include "deconvolutecommand.h"
#include "clustercommand.h"

//**********************************************************************************************************************

class ShhhSeqsCommand : public Command {
	
public:
	ShhhSeqsCommand(string);
	ShhhSeqsCommand();
	~ShhhSeqsCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "shhh.seqs";	}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Schloss PD, Gevers D, Westcott SL (2011).  Reducing the effects of PCR amplification and sequencing artifacts on 16S rRNA-based studies.  PLoS ONE.  6:e27310.\nQuince C, Lanzen A, Davenport RJ, Turnbaugh PJ (2011).  Removing noise from pyrosequenced amplicons.  BMC Bioinformatics  12:38.\nhttp://www.mothur.org/wiki/Shhh.seqs"; }
	string getDescription()		{ return "shhh.seqs"; }
	
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }		
	
private:
	bool abort;
	string outputDir, fastafile, namefile, groupfile;
	int processors;
	double sigma;
	vector<string> outputNames;
	
    int readData(correctDist*, seqNoise&, vector<string>&, vector<string>&, vector<string>&, vector<int>&);
	vector<string> createProcessesGroups(string, string, string);
	int deconvoluteResults(string, string);
};

/**************************************************************************************************/

#endif
