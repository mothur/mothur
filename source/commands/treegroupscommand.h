#ifndef TREEGROUPCOMMAND_H
#define TREEGROUPCOMMAND_H

/*
 *  treegroupscommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
 
#include "command.hpp"
#include "inputdata.h"
#include "groupmap.h"
#include "validcalculator.h"
#include "tree.h"
#include "counttable.h"
#include "readmatrix.hpp"
#include "readcolumn.h"
#include "readphylip.h"
#include "sharedsobscollectsummary.h"
#include "sharedchao1.h"
#include "sharedace.h"
#include "sharednseqs.h"
#include "sharedjabund.h"
#include "sharedsorabund.h"
#include "sharedjclass.h"
#include "sharedsorclass.h"
#include "sharedjest.h"
#include "sharedsorest.h"
#include "sharedthetayc.h"
#include "sharedthetan.h"
#include "sharedkstest.h"
#include "whittaker.h"
#include "sharedochiai.h"
#include "sharedanderbergs.h"
#include "sharedkulczynski.h"
#include "sharedkulczynskicody.h"
#include "sharedlennon.h"
#include "sharedmorisitahorn.h"
#include "sharedbraycurtis.h"
#include "whittaker.h"
#include "odum.h"
#include "canberra.h"
#include "structeuclidean.h"
#include "structchord.h"
#include "hellinger.h"
#include "manhattan.h"
#include "structpearson.h"
#include "soergel.h"
#include "spearman.h"
#include "structkulczynski.h"
#include "structchi2.h"
#include "speciesprofile.h"
#include "hamming.h"
#include "gower.h"
#include "memchi2.h"
#include "memchord.h"
#include "memeuclidean.h"
#include "mempearson.h"
#include "sharedrjsd.h"
#include "sharedjsd.h"



/* This command create a tree file for each similarity calculator at distance level, using various calculators to find the similiarity between groups. 
	The user can select the lines or labels they wish to use as well as the groups they would like included.
	They can also use as many or as few calculators as they wish. */
	
/**************************************************************************************************************/

class TreeGroupCommand : public Command {
	
public:
	TreeGroupCommand(string);	
	TreeGroupCommand();
	~TreeGroupCommand();
	
	vector<string> setParameters();
	string getCommandName()			{ return "tree.shared";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Tree.shared"; }
	string getDescription()		{ return "generate a tree file that describes the dissimilarity among groups"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	string lastLabel;
	string format, groupNames, filename, sharedfile, countfile, inputfile;
	int numGroups, subsampleSize, iters, processors;
	ofstream out;
	float precision, cutoff;
    vector<string> Treenames;

	bool abort, allLines, subsample, withReplacement;
	set<string> labels; //holds labels to be used
	string phylipfile, columnfile, namefile, calc, groups, label, outputDir;
	vector<string>  Estimators, Groups, outputNames; //holds estimators to be used
	
    int createProcesses(SharedRAbundVectors*& thisLookup, CountTable&);
    void printSims(ostream&, vector< vector<double> >&, vector<string>);
    int makeSimsShared(InputData&, SharedRAbundVectors*&, CountTable& );
    vector< vector<double> > makeSimsDist(SparseDistanceMatrix*, int);
};

/**************************************************************************************************************/
	
#endif


