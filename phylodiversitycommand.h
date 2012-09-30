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
#include "counttable.h"
#include "sharedutilities.h"
#include "tree.h"

class PhyloDiversityCommand : public Command {
	
	public:
		PhyloDiversityCommand(string);
		PhyloDiversityCommand();
		~PhyloDiversityCommand(){}
	
		vector<string> setParameters();
		string getCommandName()			{ return "phylo.diversity";			}
		string getCommandCategory()		{ return "Hypothesis Testing";		}
		string getOutputFileNameTag(string, string);
	string getHelpString();	
		string getCitation() { return "Faith DP (1994). Phylogenetic pattern and the quantification of organismal biodiversity. Philos Trans R Soc Lond B Biol Sci 345: 45-58. \nhttp://www.mothur.org/wiki/Phylo.diversity"; }
		string getDescription()		{ return "phylo.diversity"; }

		int execute();
		void help() { m->mothurOut(getHelpString()); }
private:
		CountTable* ct;
		float freq;
		int iters, processors, numUniquesInName;  
		bool abort, rarefy, summary, collect, scale;
		string groups, outputDir, treefile, groupfile, namefile, countfile;
		vector<string> Groups, outputNames; //holds groups to be used, and outputFile names
		
        map<string, int> getRootForGroups(Tree* t);
		int readNamesFile();
		void printData(set<int>&, map< string, vector<float> >&, ofstream&, int);
		void printSumData(map< string, vector<float> >&, ofstream&, int);
        vector<float> calcBranchLength(Tree*, int, vector< map<string, bool> >&, map<string, int>);
		int driver(Tree*, map< string, vector<float> >&, map<string, vector<float> >&, int, int, vector<int>&, set<int>&, ofstream&, ofstream&, bool);
		int createProcesses(vector<int>&, Tree*, map< string, vector<float> >&, map<string, vector<float> >&, int, int, vector<int>&, set<int>&, ofstream&, ofstream&);

};

#endif

