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
    string getTag() { string tag = "opti_" + metric; return tag; }
    int initialize(double&, bool);  //randomize and place in "best" OTUs
    bool update(double&); //returns whether list changed and MCC
    vector<double> getStats(unsigned long long&, unsigned long long&, unsigned long long&, unsigned long long&);
    ListVector* getList();
    
private:
    MothurOut* m;
    map<int, int> seqBin; //sequence# -> bin#
    OptiMatrix* matrix;
    vector<int> randomizeSeqs;
    vector< vector<int> > bins; //bin[0] -> seqs in bin[0]
    string metric;
    unsigned long long truePositives, trueNegatives, falsePositives, falseNegatives, numSeqs, insertLocation, totalPairs, numSingletons;
    
    double calcMCC(unsigned long long, unsigned long long, unsigned long long, unsigned long long);
    double calcSens(unsigned long long, unsigned long long, unsigned long long, unsigned long long);
    double calcSpec(unsigned long long, unsigned long long, unsigned long long, unsigned long long);
    double calcTPTN(unsigned long long, unsigned long long, unsigned long long, unsigned long long);
    double calcTP(unsigned long long, unsigned long long, unsigned long long, unsigned long long);
    double calcTN(unsigned long long, unsigned long long, unsigned long long, unsigned long long);
    double calcFP(unsigned long long, unsigned long long, unsigned long long, unsigned long long);
    double calcFN(unsigned long long, unsigned long long, unsigned long long, unsigned long long);
    double calcFPFN(unsigned long long, unsigned long long, unsigned long long, unsigned long long);
    double calcF1Score(unsigned long long, unsigned long long, unsigned long long, unsigned long long);
    double calcAccuracy(unsigned long long, unsigned long long, unsigned long long, unsigned long long);
    double calcPPV(unsigned long long, unsigned long long, unsigned long long, unsigned long long);
    double calcNPV(unsigned long long, unsigned long long, unsigned long long, unsigned long long);
    double calcFDR(unsigned long long, unsigned long long, unsigned long long, unsigned long long);
    double moveAdjustTFValues(int bin, int seq, int newBin, unsigned long long&, unsigned long long&, unsigned long long&, unsigned long long&);
};

#endif /* defined(__Mothur__opticluster__) */
