//
//  optifitcluster.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/10/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef optifitcluster_hpp
#define optifitcluster_hpp

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
#include "opticluster.h"

/***********************************************************************/

class OptiFitCluster : public Cluster {
    
public:
    
    OptiFitCluster(OptiData* mt, ClusterMetric* met, long long ns);
    ~OptiFitCluster() {}
    
    int initialize(double& value, bool randomize, vector<vector< string > > existingBins, vector<string>, string);
    bool update(double&); //returns whether list changed and MCC
    bool updateDistance(PDistCell& colCell, PDistCell& rowCell) { return false; } //inheritance compliant
    string getTag() { string tag = "optifit_" + metric->getName(); return tag; }
    long long getNumBins();
    long long getNumFitBins();
    
    vector<double>  getStats( long long&,  long long&,  long long&,  long long&);  //combo stats
    vector<double>  getFitStats( long long&,  long long&,  long long&,  long long&); //fitted seqs stats
    
    ListVector* getList();
    ListVector* getFittedList(long long&);
    
protected:
    MothurOut* m;
    Utils util;
    ClusterMetric* metric;
    OptiData* matrix;
    
    map<int, int> seqBin; //sequence# -> bin#
    vector<int> randomizeSeqs;
    vector< vector<int> > bins; //bin[0] -> seqs in bin[0]
    map<int, string> binLabels; //for fitting - maps binNumber to existing reference label
    set<int> fitSeqs; //matrix indexes for movable sequences
    long long maxRefBinNumber;
    bool closed;
    
    long long fittruePositives, fittrueNegatives, fitfalsePositives, fitfalseNegatives, numFitSeqs, insertLocation, numFitSingletons;
    long long combotruePositives, combotrueNegatives, combofalsePositives, combofalseNegatives, numComboSeqs, numComboSingletons;
    
    int findInsert();
    vector<long long> getCloseFarCounts(int seq, int newBin);
    vector<long long> getCloseFarFitCounts(int seq, int newBin);
    void clusterUnfitted(OptiData*);
    
};


#endif /* optifitcluster_hpp */
