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
    ClusterFitCommand();
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
    bool abort, sim, print_start, selfReference;
    string refcolumnfile, refphylipfile, refdistfile, reffastafile, refnamefile, refcountfile, reflistfile, refNameOrCount;
    string namefile, refformat, distfile, countfile, fastafile, columnfile, nameOrCount;
    string comboDistFile;
    
    string method, fileroot, tag, outputDir, inputDir, metric, initialize, metricName;
    double cutoff, stableMetric;
    float adjust;
    int precision, length, maxIters, processors;
    vector<string> outputNames;
    unsigned long loops;
    
    string runRefOptiCluster(OptiData*&, ClusterMetric*&, ListVector*&, map<string, int>&, string);
    void createReferenceNameCount();
    string calcDists();
    void runSensSpec(string listFileName, string distFileName, string dupsFile, string dupsFormat);
    void outputSteps(string outputName, bool printHeaders, long long tp, long long tn, long long fp, long long fn, vector<double> results, long long numBins, long long fittp, long long fittn, long long fitfp, long long fitfn, vector<double> fitresults, long long numFitBins, int);
};

#endif /* clusterfitcommand_hpp */
