//
//  clusterfitcommand.hpp
//  Mothur
//
//  Created by Sarah Westcott on 1/22/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef clusterfitcommand_hpp
#define clusterfitcommand_hpp

#include "command.hpp"
#include "listvector.hpp"
#include "cluster.hpp"
#include "counttable.h"
#include "optifitcluster.hpp"
#include "optirefmatrix.hpp"
#include "calculator.h"
#include "distancecommand.h"
#include "aligncommand.h"
#include "filterseqscommand.h"
#include "deconvolutecommand.h"
#include "listseqscommand.h"
#include "getdistscommand.h"
#include "getseqscommand.h"

class ClusterFitCommand : public Command {
    
public:
    ClusterFitCommand(string);
    ~ClusterFitCommand();
    
    vector<string> setParameters();
    string getCommandName()			{ return "cluster.fit";		}
    string getCommandCategory()		{ return "Clustering";      }
    
    string getHelpString();
    string getOutputPattern(string);
    string getCitation() { return "\nhttp://www.mothur.org/wiki/Cluster.fit"; }
    string getDescription()		{ return "fit your sequences into existing OTUs"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    bool abort, sim, print_start, selfReference, printref;
    string refdistfile, reffastafile, refnamefile, refcountfile, reflistfile, refNameOrCount;
    string namefile, refformat, distfile, countfile, fastafile, columnfile, nameOrCount, accnosfile;
    string comboDistFile;
    
    string method, fileroot, tag,  inputDir, metric, initialize, metricName, criteria, refWeight;
    double cutoff, stableMetric;
    float adjust, fitPercent;
    int precision, length, maxIters, processors, denovoIters;
    vector<string> outputNames, listFiles;
    unsigned long loops;
    
    ListVector* runUserRefOptiCluster(OptiData*&, ClusterMetric*&, map<string, int>&, string);
    string runRefOptiCluster(OptiData*&, ClusterMetric*&, ListVector*&, map<string, int>&, string);
    string runDenovoOptiCluster(OptiData*&, ClusterMetric*&, map<string, int>&, string);
    ListVector* clusterRefs(OptiData*& refsMatrix, ClusterMetric*&);
    void createReferenceNameCount();
    string calcDists();
    string runSensSpec(OptiData*& matrix, ClusterMetric*& userMetric, ListVector*& list, map<string, int>& counts);
    string runSensSpec(string distFileName, string dupsFile, string dupsFormat, ClusterMetric*&, string);
    void outputSteps(string outputName, bool& printHeaders, double tp, double tn, double fp, double fn, vector<double> results, long long numBins, double fittp, double fittn, double fitfp, double fitfn, vector<double> fitresults, long long numFitBins, int, bool, int);
};

#endif /* clusterfitcommand_hpp */
