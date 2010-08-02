#ifndef PHYLODIVERSITYCOMMAND_H
#define PHYLODIVERSITYCOMMAND_H

/*
 *  phylodiversitycommand.h
 *  Mothur
 *
 *  Created by westcott on 4/30/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "treemap.h"
#include "globaldata.hpp"

class PhyloDiversityCommand : public Command {
	
	public:
		PhyloDiversityCommand(string);
		~PhyloDiversityCommand();
		int execute();	
		void help();
	
	private:
		GlobalData* globaldata;
		
		float freq;
		int iters;  
		bool abort, rarefy;
		string groups, outputDir;
		vector<string> Groups, outputNames; //holds groups to be used, and outputFile names
		
		void printData(set<int>&, map< string, vector<float> >&, string);
		float calcBranchLength(Tree*, int);
};

#endif

