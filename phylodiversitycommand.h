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
#include "readtree.h"
#include "sharedutilities.h"


class PhyloDiversityCommand : public Command {
	
	public:
		PhyloDiversityCommand(string);
		PhyloDiversityCommand();
		~PhyloDiversityCommand(){}
	
		vector<string> setParameters();
		string getCommandName()			{ return "phylo.diversity";			}
		string getCommandCategory()		{ return "Hypothesis Testing";		}
		string getHelpString();	
		string getCitation() { return "http://www.mothur.org/wiki/Phylo.diversity"; }
	
		int execute();
		void help() { m->mothurOut(getHelpString()); }
private:
		ReadTree* read;
		TreeMap* tmap;
		float freq;
		int iters, processors, numUniquesInName;  
		bool abort, rarefy, summary, collect, scale;
		string groups, outputDir, treefile, groupfile, namefile;
		vector<string> Groups, outputNames; //holds groups to be used, and outputFile names
		map<string, string> nameMap;
		
		int readNamesFile();
		void printData(set<int>&, map< string, vector<float> >&, ofstream&, int);
		void printSumData(map< string, vector<float> >&, ofstream&, int);
		vector<float> calcBranchLength(Tree*, int, map< string, set<int> >&);
		int driver(Tree*, map< string, vector<float> >&, map<string, vector<float> >&, int, int, vector<int>&, set<int>&, ofstream&, ofstream&, bool);
		int createProcesses(vector<int>&, Tree*, map< string, vector<float> >&, map<string, vector<float> >&, int, int, vector<int>&, set<int>&, ofstream&, ofstream&);

};

#endif

