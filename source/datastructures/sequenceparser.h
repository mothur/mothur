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


#include "utils.hpp"
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
	
		SequenceParser(string, string, vector<string>);			//group, fasta, groups (if blanks then all) - file mismatches will set m->setControl_pressed(true)
		SequenceParser(string, string, string, vector<string>);	//group, fasta, name, groups (if blanks then all) - file mismatches will set m->setControl_pressed(true)
		~SequenceParser();
		
		//general operations
		int getNumGroups();
		vector<string> getNamesOfGroups();	
		
		int getNumSeqs(string);		//returns the number of unique sequences in a specific group
		vector<Sequence> getSeqs(string); //returns unique sequences in a specific group
        bool fillWeighted(vector< seqPNode* >&, string group, int& length); //return true for aligned and false for unaligned
		map<string, string> getNameMap(string); //returns seqName -> namesOfRedundantSeqs separated by commas for a specific group - the name file format, but each line is parsed by group.
		
		int getSeqs(string, string, string, string, long long&, bool); //prints unique sequences in a specific group to a file - group, filename, uchimeFormat=false, tag(/ab= or ;size=), tag2(/ or ;)
		int getNameMap(string, string); //print seqName -> namesOfRedundantSeqs separated by commas for a specific group - group, filename
		
        map<string, string> getAllSeqsMap();  //returns map where the key=sequenceName and the value=representativeSequence - helps us remove duplicates after group by group processing
    
	private:
		MothurOut* m;
        GroupMap groupMap;
        Utils util;
        bool hasName;
        map<string, string> allSeqsMap;
        map<string, string> nameMap;
	
        vector<Sequence> seqs; //unique
        map<string, int> groupIndexMap; //maps sample name to index in namesOfGroups and groupsToSeqs
        vector< vector<int> > groupToSeqs; //maps group -> vector of sequences indexes that belong to the group. (groupToSeqs[0] contains sequence info for sample namesOfGroups[0]).  groupToSeqs[0] = <1,4,7,8> means seqs[1], seqs[4], seq[7], seqs[8] belong to group namesOfGroups[0]. namesOfGroups in same order as groupToSeqs
        vector<string> namesOfGroups; //namesOfGroups in same order as groupToSeqs
    
        int readFasta(string fastafile);
};

#endif

