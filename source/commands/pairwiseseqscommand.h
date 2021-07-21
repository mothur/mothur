#ifndef PAIRWISESEQSCOMMAND_H
#define PAIRWISESEQSCOMMAND_H

/*
 *  pairwiseseqscommand.h
 *  Mothur
 *
 *  Created by westcott on 10/20/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "database.hpp"
#include "alignment.hpp"
#include "validcalculator.h"
#include "calculator.h"
#include "sequencedb.h"
#include "sequence.hpp"

#include "gotohoverlap.hpp"
#include "needlemanoverlap.hpp"
#include "blastalign.hpp"
#include "noalign.hpp"

#include "ignoregaps.h"
#include "eachgapdist.h"
#include "eachgapignore.h"
#include "onegapdist.h"
#include "onegapignore.h"
#include "writer.h"

class PairwiseSeqsCommand : public Command {
	
public:
	PairwiseSeqsCommand(string);	
	~PairwiseSeqsCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "pairwise.seqs";		}
	string getCommandCategory()		{ return "Sequence Processing"; }
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Needleman SB, Wunsch CD (1970). A general method applicable to the search for similarities in the amino acid sequence of two proteins. J Mol Biol 48: 443-53. [ for needleman ]\nGotoh O (1982). An improved algorithm for matching biological sequences. J Mol Biol 162: 705-8. [ for gotoh ] \nhttp://www.mothur.org/wiki/Pairwise.seqs"; }
	string getDescription()		{ return "calculates pairwise distances from an unaligned fasta file"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	SequenceDB alignDB;
	
	void createProcesses(string);
    bool sanityCheck();
    
    bool abort, countends, compress, fitCalc;
	string fastaFileName, align, calc,  output, oldfastafile, column;
	float match, misMatch, gapOpen, gapExtend, cutoff;
	int processors, longestBase, numDistsBelowCutoff;
	vector<string> Estimators, outputNames;
    
    vector< vector< int > > kmerDB; //kmerDB[0] = vector<int> maxKmers long, contains kmer counts
    vector< int > lengths;
	
	
};

#endif

