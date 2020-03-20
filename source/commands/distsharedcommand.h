#ifndef MATRIXOUTPUTCOMMAND_H
#define MATRIXOUTPUTCOMMAND_H

/*
 *  distsharedcommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/20/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */ 
#include "command.hpp"
#include "inputdata.h"
#include "groupmap.h"
#include "validcalculator.h"
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
#include "sharedjsd.h"
#include "sharedrjsd.h"


// aka. dist.shared()

/* This command create a tree file for each similarity calculator at distance level, using various calculators to find the similiarity between groups. 
	The user can select the labels they wish to use as well as the groups they would like included.
	They can also use as many or as few calculators as they wish. */
	

class DistSharedCommand : public Command {
	
public:
	DistSharedCommand(string);
	~DistSharedCommand();
	
	vector<string> setParameters();
	string getCommandName()			{ return "dist.shared";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Dist.shared"; }
	string getDescription()		{ return "generate a distance matrix that describes the dissimilarity among multiple groups"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
    string exportFileName, output, sharedfile;
	int numGroups, processors, iters, subsampleSize;
	ofstream out;

	bool abort, allLines, subsample, withReplacement;
	set<string> labels; //holds labels to be used
	string outputFile, calc, groups, label, outputDir, mode;
	vector<string>  Estimators, Groups, outputNames; //holds estimators to be used
	
    int createProcesses(SharedRAbundVectors*&);
	int driver(vector<SharedRAbundVector*>&, vector< vector<seqDist> >&, vector<Calculator*>);
    void printDists(ostream&, vector< vector<double> >&, vector<string>);

};
	
/**************************************************************************************************/
struct distSharedData {
    SharedRAbundVectors* thisLookup;
    vector< vector< vector<seqDist> > > calcDistsTotals;  //each iter, one for each calc, then each groupCombos dists. this will be used to make .dist files
    vector< vector< vector<double> > > matrices; //for each calculator a square matrix to represent the distances, only filled by main thread
    vector<string>  Estimators;
    long long numIters;
	MothurOut* m;
    int count, subsampleSize;
    bool mainThread, subsample, withReplacement;
	
	distSharedData(){}
	distSharedData(long long st, bool mt, bool su, int subsize, bool wr, vector<string> est, SharedRAbundVectors* lu) {
        m = MothurOut::getInstance();
		numIters = st;
        Estimators = est;
        thisLookup = lu;
        count = 0;
        mainThread = mt;
        subsample = su;
        subsampleSize = subsize;
        withReplacement = wr;
	}
};
/**************************************************************************************************/

#endif

