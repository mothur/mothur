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
		PhyloDiversityCommand();
		~PhyloDiversityCommand();
		vector<string> getRequiredParameters();
		vector<string> getValidParameters();
		vector<string> getRequiredFiles();
		map<string, vector<string> > getOutputFiles() { return outputTypes; }
		int execute();	
		void help();
	
	private:
		GlobalData* globaldata;
		
		float freq;
		int iters, processors;  
		bool abort, rarefy, summary, collect, scale;
		string groups, outputDir;
		vector<string> Groups, outputNames; //holds groups to be used, and outputFile names
		map<string, vector<string> > outputTypes;
		
		void printData(set<int>&, map< string, vector<float> >&, ofstream&, int);
		void printSumData(map< string, vector<float> >&, ofstream&, int);
		float calcBranchLength(Tree*, int, map< string, set<int> >&);
		int driver(Tree*, map< string, vector<float> >&, map<string, vector<float> >&, int, int, vector<int>&, set<int>&, ofstream&, ofstream&, bool);
		int createProcesses(vector<int>&, Tree*, map< string, vector<float> >&, map<string, vector<float> >&, int, int, vector<int>&, set<int>&, ofstream&, ofstream&);

};

#endif

