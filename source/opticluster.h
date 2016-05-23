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
    OptiCluster(OptiMatrix* mt, string met, double mc) : Cluster() {
        m = MothurOut::getInstance(); matrix = mt; metric = met; truePositives = 0; trueNegatives = 0; falseNegatives = 0; falsePositives = 0; metricCutoff = mc;
    }
    ~OptiCluster() {}
    bool updateDistance(PDistCell& colCell, PDistCell& rowCell) { return false; } //inheritance compliant
    string getTag() { return("opti"); }
    
    int initialize(double&);  //randomize and place in "best" OTUs
    bool update(double&); //returns whether list changed and MCC
    ListVector* getList();
    
private:
    MothurOut* m;
    map<int, int> seqBin;
    ListVector* list;
    OptiMatrix* matrix;
    vector< vector<int> > bins; //bin[0] -> seqs in bin[0]
    string metric;
    int truePositives, trueNegatives, falsePositives, falseNegatives, numSeqs, insertLocation;
    double metricCutoff;
    
    double calcMCC(double, double, double, double);
    int removeDups(set<int>& left, set<int>& right);
    double moveAdjustTFValues(int bin, int seq, int newBin, double&, double&, double&, double&);
    int eraseIndex(vector<int>& otus, int value);
    
};

#endif /* defined(__Mothur__opticluster__) */
