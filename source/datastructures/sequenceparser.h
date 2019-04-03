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
#include "splitgroupscommand.h"

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
		
        vector<string> getFiles(string);  //returns fasta and count file a specific group.
        map<string, vector<string> > getFiles() { return groupToFiles; } //returns all files groupName - > vector of groups files (fasta, optionalName, group);
    
	private:
		MothurOut* m;
        Utils util;
        bool hasName;
    
        map<string, vector<string> > groupToFiles; //groupName -> fasta, name, group or groupName -> fasta, group
        vector<string> namesOfGroups; //namesOfGroups in same order as groupToSeqs;
 };   

#endif

