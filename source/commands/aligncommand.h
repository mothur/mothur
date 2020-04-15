#ifndef ALIGNCOMMAND_H
#define ALIGNCOMMAND_H
/*
 *  aligncommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/15/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "database.hpp"
#include "alignment.hpp"
#include "alignmentdb.h"
#include "sequence.hpp"
#include "gotohoverlap.hpp"
#include "needlemanoverlap.hpp"
#include "blastalign.hpp"
#include "noalign.hpp"
#include "nast.hpp"
#include "nastreport.hpp"

//test
class AlignCommand : public Command {
	
#ifdef UNIT_TEST
    friend class TestAlignSeqsIntegration;
#endif
    
public:
	AlignCommand(string);	
	~AlignCommand();
	
	vector<string> setParameters();
	string getCommandName()			{ return "align.seqs";			}
	string getCommandCategory()		{ return "Sequence Processing"; }
	
	string getHelpString();
    string getCommonQuestions();
    string getOutputPattern(string);	
	string getCitation() { return "DeSantis TZ, Jr., Hugenholtz P, Keller K, Brodie EL, Larsen N, Piceno YM, Phan R, Andersen GL (2006). NAST: a multiple sequence alignment server for comparative analysis of 16S rRNA genes. Nucleic Acids Res 34: W394-9.\nSchloss PD (2009). A high-throughput DNA sequence aligner for microbial ecology studies. PLoS ONE 4: e8230.\nSchloss PD (2010). The effects of alignment quality, distance calculation method, sequence filtering, and region on the analysis of 16S rRNA gene-based studies. PLoS Comput Biol 6: e1000844.\nhttp://www.mothur.org/wiki/Align.seqs http://www.mothur.org/wiki/Align.seqs"; }
	string getDescription()		{ return "align sequences"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
protected:
	
	long long createProcesses(string, string, string, string, vector<long long>&);
	void appendReportFiles(string, string);
		
	string fastafile, templateFileName, distanceFileName, search, align;
	float match, misMatch, gapOpen, gapExtend, threshold;
	int processors, kmerSize;
	vector<string> outputNames;
	
	bool abort, flip, calledHelp, save;
    AlignmentDB* templateDB;

};

#endif
