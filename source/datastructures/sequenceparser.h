#ifndef SEQUENCEPARSER_H
#define SEQUENCEPARSER_H

/*
 *  sequenceParser.h
 *  Mothur
 *
 *  Created by westcott on 9/9/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */


#include "mothur.h"
#include "mothurout.h"
#include "sequence.hpp"
#include "groupmap.h"

/* This class reads a fasta and group file with a namesfile as optional and parses the data by group. 
 
	Note: The sum of all the groups unique sequences will be larger than the original number of unique sequences. 
	This is because when we parse the name file we make a unique for each group instead of 1 unique for all
	groups. 
 
 */

class SequenceParser {
	
	public:
	
		SequenceParser(string, string);			//group, fasta - file mismatches will set m->control_pressed = true
		SequenceParser(string, string, string);	//group, fasta, name  - file mismatches will set m->control_pressed = true
		~SequenceParser();
		
		//general operations
		int getNumGroups();
		vector<string> getNamesOfGroups();	
		bool isValidGroup(string);  //return true if string is a valid group
		
		int getNumSeqs(string);		//returns the number of unique sequences in a specific group
		vector<Sequence> getSeqs(string); //returns unique sequences in a specific group
		map<string, string> getNameMap(string); //returns seqName -> namesOfRedundantSeqs separated by commas for a specific group - the name file format, but each line is parsed by group.
		
		int getSeqs(string, string, string, string, bool); //prints unique sequences in a specific group to a file - group, filename, uchimeFormat=false, tag(/ab= or ;size=), tag2(/ or ;)
		int getNameMap(string, string); //print seqName -> namesOfRedundantSeqs separated by commas for a specific group - group, filename
		
		map<string, string> getAllSeqsMap(){ return allSeqsMap; }  //returns map where the key=sequenceName and the value=representativeSequence - helps us remove duplicates after group by group processing
	private:
	
		GroupMap* groupMap;
		MothurOut* m;
	
		int numSeqs;
		map<string, string> allSeqsMap;
		map<string, vector<Sequence> > seqs; //a vector for each group
		map<string, map<string, string> > nameMapPerGroup; //nameMap for each group
};

#endif

