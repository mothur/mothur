#ifndef PRECLUSTERCOMMAND_H
#define PRECLUSTERCOMMAND_H


/*
 *  preclustercommand.h
 *  Mothur
 *
 *  Created by westcott on 12/21/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"
#include "sequence.hpp"
#include "sequenceparser.h"
#include "sequencecountparser.h"
#include "alignment.hpp"
#include "gotohoverlap.hpp"
#include "needlemanoverlap.hpp"
#include "blastalign.hpp"
#include "noalign.hpp"
#include "filters.h"
#include "getseqscommand.h"


//************************************************************/
class PreClusterCommand : public Command {

public:
	PreClusterCommand(string);
	~PreClusterCommand(){}

	vector<string> setParameters();
	string getCommandName()			{ return "pre.cluster";				}
	string getCommandCategory()		{ return "Sequence Processing";		}

	string getHelpString();
    string getOutputPattern(string);
	string getCitation() { return "Schloss PD, Gevers D, Westcott SL (2011).  Reducing the effects of PCR amplification and sequencing artifacts on 16S rRNA-based studies.  PLoS ONE.  6:e27310.\nhttp://www.mothur.org/wiki/Pre.cluster"; }
	string getDescription()		{ return "implements a pseudo-single linkage algorithm with the goal of removing sequences that are likely due to pyrosequencing errors"; }

	int execute();
	void help() { m->mothurOut(getHelpString()); }

private:
    int diffs, length, processors;
    float match, misMatch, gapOpen, gapExtend, alpha, delta, error_rate, indel_prob, max_indels;
    vector<float> error_dist;
    bool abort, bygroup;
    string fastafile, outputDir, countfile, pc_method, align_method, align;
    vector<string> outputNames;
    
    void createProcessesGroups(map<string, vector<string> >&, vector<string>, string, string);
    string mergeGroupCounts(string, string);
    void printFasta(string newFastaFileName, string accnosFile);
};


/**************************************************************************************************/

#endif
