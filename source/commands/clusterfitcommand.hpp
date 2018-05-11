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
    string refcolumnfile, reffastafile, refnamefile, refcountfile, reflistfile, refNameOrCount;
    string namefile, format, distfile, countfile, fastafile, columnfile, nameOrCount;
    string comboDistFile;
    
    string method, fileroot, tag, outputDir, inputDir, metric, initialize, metricName;
    double cutoff, stableMetric;
    float adjust;
    int precision, length, maxIters, processors;
    vector<string> outputNames;
    unsigned long loops;
    //long long truePositives, falsePositives, trueNegatives, falseNegatives;
    map<string, int> counts;
    
    int runOptiCluster(ListVector*&);
    void createReferenceNameCount();
    string calcDists();
};

#endif /* clusterfitcommand_hpp */
