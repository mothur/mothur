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
    int initialize(double&, bool, string);  //randomize and place in "best" OTUs
    bool update(double&); //returns whether list changed and MCC
    vector<double> getStats( long long&,  long long&,  long long&,  long long&);
    ListVector* getList();
    
private:
    MothurOut* m;
    map<int, int> seqBin; //sequence# -> bin#
    OptiMatrix* matrix;
    vector<int> randomizeSeqs;
    vector< vector<int> > bins; //bin[0] -> seqs in bin[0]
    string metric;
    long long truePositives, trueNegatives, falsePositives, falseNegatives, numSeqs, insertLocation, totalPairs, numSingletons;
    
    int findInsert();
    double calcMCC(long long, long long, long long, long long);
    double calcSens( long long,  long long,  long long,  long long);
    double calcSpec( long long,  long long,  long long,  long long);
    double calcTPTN( long long,  long long,  long long,  long long);
    double calcTP( long long,  long long,  long long,  long long);
    double calcTN( long long,  long long,  long long,  long long);
    double calcFP( long long,  long long,  long long,  long long);
    double calcFN( long long,  long long,  long long,  long long);
    double calcFPFN( long long,  long long,  long long,  long long);
    double calcF1Score( long long,  long long,  long long,  long long);
    double calcAccuracy( long long,  long long,  long long,  long long);
    double calcPPV( long long,  long long,  long long,  long long);
    double calcNPV( long long,  long long,  long long,  long long);
    double calcFDR( long long,  long long,  long long,  long long);
    double moveAdjustTFValues(int bin, int seq, int newBin,  long long&,  long long&,  long long&,  long long&);
};

#endif /* defined(__Mothur__opticluster__) */
