//
//  getmetacommunitycommand.h
//  Mothur
//
//  Created by SarahsWork on 4/9/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef Mothur_getmetacommunitycommand_h
#define Mothur_getmetacommunitycommand_h

#include "command.hpp"
#include "inputdata.h"
#include "qFinderDMM.h"
#include "pam.h"
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
#include "sharedjackknife.h"
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

/**************************************************************************************************/

class GetMetaCommunityCommand : public Command {
public:
    GetMetaCommunityCommand(string);
    GetMetaCommunityCommand();
    ~GetMetaCommunityCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "get.communitytype";		}
    string getCommandCategory()		{ return "OTU-Based Approaches";         }
    
    string getOutputPattern(string);
    
	string getHelpString();
    string getCitation() { return "Holmes I, Harris K, Quince C (2012) Dirichlet Multinomial Mixtures: Generative Models for Microbial Metagenomics. PLoS ONE 7(2): e30126. doi:10.1371/journal.pone.0030126 http://www.mothur.org/wiki/get.communitytype"; }
    string getDescription()		{ return "Assigns samples to bins using a Dirichlet multinomial mixture model"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    struct linePair {
		unsigned long long start;
		unsigned long long end;
		linePair(unsigned long long i, unsigned long long j) : start(i), end(j) {}
	};
    bool abort, allLines, subsample;
    string outputDir;
    vector<string> outputNames;
    string sharedfile, method, calc;
    int minpartitions, maxpartitions, optimizegap, processors, iters, subsampleSize;
    vector<string> Groups, Estimators;
    set<string> labels;
    
    vector<vector<double> > generateDistanceMatrix(vector<SharedRAbundVector*>& lookup);
    int driver(vector<SharedRAbundVector*> thisLookup, vector< vector<seqDist> >& calcDists, Calculator*);
    int processDriver(vector<SharedRAbundVector*>&, vector<int>&, string, vector<string>, vector<string>, vector<string>, int);
    int createProcesses(vector<SharedRAbundVector*>&);
    vector<double> generateDesignFile(int, map<string,string>);
    int generateSummaryFile(int, map<string,string>, vector<double>);

};

/**************************************************************************************************/
struct summaryData {
    
    string name;
    double refMean, difference;
    vector<double> partMean, partLCI, partUCI;
    
};

#endif
