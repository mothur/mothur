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

/***********************************************************************/

class OptiCluster : public Cluster {
    
public:
    OptiCluster(OptiMatrix* mt, string met, int i, double st) {
        matrix = mt; metric = met; truePositives = 0; trueNegatives = 0; falseNegatives = 0; falsePositives = 0; maxIters = i; stableMetric = st;
    }
    ~OptiCluster() { delete matrix; }
    bool updateDistance(PDistCell& colCell, PDistCell& rowCell) { return false; } //inheritance compliant
    string getTag() { return("opti"); }
    
    int initialize();  //randomize and place in "best" OTUs
    bool update(double&); //returns whether list changed and MCC
    ListVector* getList();
    
private:
    map<int, int> seqBin;
    ListVector* list;
    OptiMatrix* matrix;
    vector< vector<int> > bins; //bin[0] -> seqs in bin[0]
    string metric;
    double listVectorMetric, stableMetric;
    int truePositives, trueNegatives, falsePositives, falseNegatives, maxIters, numSeqs, insertLocation;
    
    double calcMCC(double, double, double, double);
    double moveAdjustTFValues(int bin, int seq, int newBin);
    int eraseIndex(vector<int>& otus, int value);
    
};

#endif /* defined(__Mothur__opticluster__) */
