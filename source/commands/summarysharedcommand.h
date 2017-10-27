#ifndef SUMMARYSHAREDCOMMAND_H
#define SUMMARYSHAREDCOMMAND_H
/*
 *  summarysharedcommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "command.hpp"
#include "inputdata.h"
#include "calculator.h"
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
//#include "sharedjackknife.h"
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

class SummarySharedCommand : public Command {

public:
	SummarySharedCommand(string);
	SummarySharedCommand();
	~SummarySharedCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "summary.shared";			}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Summary.shared"; }
	string getDescription()		{ return "generate a summary file containing calculator values for each line in the OTU data and for all possible comparisons between groups"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	vector<linePair> lines;
	
	bool abort, allLines, mult, all, createPhylip, subsample;
	set<string> labels; //holds labels to be used
	string label, calc, groups, sharedfile, output;
	vector<string>  Estimators, Groups, outputNames, sumCalculatorsNames;
	
	string format, outputDir;
	int numGroups, processors, subsampleSize, iters, numCalcs;
	int process(vector<SharedRAbundVector*>, string, string, vector<string>);
    int printSims(ostream&, vector< vector<double> >&, vector<string>);
    int runCalcs(vector<SharedRAbundVector*>&, string, string, vector< vector<seqDist>  >&);

};

/**************************************************************************************************/
//custom data structure for threads to use.
//main process handling the calcs that can do more than 2 groups
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct summarySharedData {
    vector<SharedRAbundVector*> thisLookup;
    vector< vector<seqDist> > calcDists;
    vector<string>  Estimators;
	unsigned long long start;
	unsigned long long end;
	MothurOut* m;
	string sumFile, sumAllFile;
    int count;
    bool main, mult;
    
	summarySharedData(){}
	summarySharedData(string sf, string sfa, MothurOut* mout, unsigned long long st, unsigned long long en, vector<string> est, vector<SharedRAbundVector*> lu, bool mai, bool mu) {
		sumFile = sf;
        sumAllFile = sfa;
		m = mout;
		start = st;
		end = en;
        Estimators = est;
        thisLookup = lu;
        count=0;
        main = mai;
        mult = mu;
	}
    ~summarySharedData() {  }
};
/**************************************************************************************************/

#endif
