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
    
#ifdef UNIT_TEST
    friend class TestOptiCluster;
    OptiCluster() : Cluster() { m = MothurOut::getInstance(); truePositives = 0; trueNegatives = 0; falseNegatives = 0; falsePositives = 0; } //for testing class
    void setVariables(OptiMatrix* mt, string met) { matrix = mt; metric = met; }
#endif
    
public:
    OptiCluster(OptiMatrix* mt, string met, double ns) : Cluster() {
        m = MothurOut::getInstance(); matrix = mt; metric = met; truePositives = 0; trueNegatives = 0; falseNegatives = 0; falsePositives = 0; numSingletons = ns;
    }
    ~OptiCluster() {}
    bool updateDistance(PDistCell& colCell, PDistCell& rowCell) { return false; } //inheritance compliant
    string getTag() { return("opti"); }
    int initialize(double&, bool);  //randomize and place in "best" OTUs
    bool update(double&); //returns whether list changed and MCC
    vector<double> getStats();
    ListVector* getList();
    
private:
    MothurOut* m;
    map<int, int> seqBin;
    ListVector* list;
    OptiMatrix* matrix;
    vector<int> randomizeSeqs;
    vector< vector<int> > bins; //bin[0] -> seqs in bin[0]
    string metric;
    double truePositives, trueNegatives, falsePositives, falseNegatives, numSeqs, insertLocation, totalPairs, numSingletons;
    
    double calcMCC(double, double, double, double);
    double calcSens(double, double, double, double);
    double calcSpec(double, double, double, double);
    double calcTPTN(double, double, double, double);
    double calcTP2TN(double, double, double, double);
    double calcFPFN(double tp, double tn, double fp, double fn);
    double moveAdjustTFValues(int bin, int seq, int newBin, double&, double&, double&, double&);
};

#endif /* defined(__Mothur__opticluster__) */
