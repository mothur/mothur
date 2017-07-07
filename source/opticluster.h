//
//  opticluster.h
//  Mothur
//
//  Created by Sarah Westcott on 4/20/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__opticluster__
#define __Mothur__opticluster__

#include "cluster.hpp"
#include "optimatrix.h"
#include "calculator.h"
#include "mcc.hpp"
#include "sensitivity.hpp"
#include "specificity.hpp"
#include "fdr.hpp"
#include "npv.hpp"
#include "ppv.hpp"
#include "f1score.hpp"
#include "tp.hpp"
#include "fp.hpp"
#include "fpfn.hpp"
#include "tptn.hpp"
#include "tn.hpp"
#include "fn.hpp"
#include "accuracy.hpp"

/***********************************************************************/

class OptiCluster : public Cluster {

public:
    
#ifdef UNIT_TEST
    friend class TestOptiCluster;
    OptiCluster() : Cluster() { m = MothurOut::getInstance(); truePositives = 0; trueNegatives = 0; falseNegatives = 0; falsePositives = 0; } //for testing class
    void setVariables(OptiMatrix* mt, ClusterMetric* met) { matrix = mt; metric = met; }
#endif
    
    OptiCluster(OptiMatrix* mt, ClusterMetric* met, long long ns) : Cluster() {
        m = MothurOut::getInstance(); matrix = mt; metric = met; truePositives = 0; trueNegatives = 0; falseNegatives = 0; falsePositives = 0; numSingletons = ns;
    }
    ~OptiCluster() {}
    bool updateDistance(PDistCell& colCell, PDistCell& rowCell) { return false; } //inheritance compliant
    string getTag() { string tag = "opti_" + metric->getName(); return tag; }
    long long getNumBins();
    int initialize(double&, bool, string);  //randomize and place in "best" OTUs
    bool update(double&); //returns whether list changed and MCC
    vector<double> getStats( long long&,  long long&,  long long&,  long long&);
    ListVector* getList();
    
protected:
    MothurOut* m;
    map<int, int> seqBin; //sequence# -> bin#
    OptiMatrix* matrix;
    vector<int> randomizeSeqs;
    vector< vector<int> > bins; //bin[0] -> seqs in bin[0]
    ClusterMetric* metric;
    long long truePositives, trueNegatives, falsePositives, falseNegatives, numSeqs, insertLocation, totalPairs, numSingletons;
    
    int findInsert();
    vector<long long> getCloseFarCounts(int seq, int newBin);
};

#endif /* defined(__Mothur__opticluster__) */
