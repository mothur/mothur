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
    
    OptiFitCluster(OptiData* mt, ClusterMetric* met, long long ns, string crit);
    ~OptiFitCluster() {}
    
    int initialize(double& value, bool randomize, vector<vector< string > > existingBins, vector<string>, string, bool);
    bool update(double&); //returns whether list changed and MCC
    bool updateDistance(PDistCell& colCell, PDistCell& rowCell) { return false; } //inheritance compliant
    string getTag() { string tag = "optifit_" + metric->getName(); return tag; }
    long long getNumBins();
    long long getNumFitBins();
    
    vector<double>  getStats( long long&,  long long&,  long long&,  long long&);  //combo stats
    vector<double>  getFitStats( long long&,  long long&,  long long&,  long long&); //fitted seqs stats
    
    ListVector* getList();
    ListVector* getFittedList(string, bool);
    
protected:
    MothurOut* m;
    Utils util;
    ClusterMetric* metric;
    OptiData* matrix;
    string criteria;
    
    map<long long, long long> seqBin; //sequence# -> bin#
    vector<long long> randomizeSeqs;
    vector< vector<long long> > bins; //bin[0] -> seqs in bin[0]
    map<long long, string> binLabels; //for fitting - maps binNumber to existing reference label
    set<long long> fitSeqs; //matrix indexes for movable sequences
    long long maxRefBinNumber;
    bool closed, denovo;
    
    long long fittruePositives, fittrueNegatives, fitfalsePositives, fitfalseNegatives, numFitSeqs, insertLocation, numFitSingletons; 
    long long combotruePositives, combotrueNegatives, combofalsePositives, combofalseNegatives, numComboSeqs, numComboSingletons;
    
    int findInsert();
    vector<long long> getCloseFarCounts(long long seq, long long newBin);
    vector<long long> getCloseFarFitCounts(long long seq, long long newBin);
    ListVector* clusterUnfitted(OptiData*, string);
    
};


#endif /* optifitcluster_hpp */
