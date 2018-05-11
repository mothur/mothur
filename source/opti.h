//
//  opti.h
//  Mothur
//
//  Created by Sarah Westcott on 5/10/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef opti_h
#define opti_h

class Opti : public Cluster {
    
public:
    
#ifdef UNIT_TEST
    friend class TestOpti;
    Opti() : Cluster() { m = MothurOut::getInstance(); } //for testing class
    void setVariables(OptiMatrix* mt, ClusterMetric* met) { matrix = mt; metric = met; }
#endif
    
    OptiCluster(OptiMatrix* mt, ClusterMetric* met, long long ns) : Cluster() {
        m = MothurOut::getInstance(); matrix = mt; metric = met; truePositives = 0; trueNegatives = 0; falseNegatives = 0; falsePositives = 0; numSeqs = 0; numSingletons = ns; numRefSeqs = 0; numRefSingletons = 0; fitCalc=false; reftruePositives = 0; reffalsePositives = 0; reffalseNegatives = 0; reftrueNegatives = 0;
    }
    ~OptiCluster() {}
    bool updateDistance(PDistCell& colCell, PDistCell& rowCell) { return false; } //inheritance compliant
    string getTag() { string tag = "opti_" + metric->getName(); return tag; }
    long long getNumBins();
    int initialize(double&, bool, string);  //randomize and place in "best" OTUs
    int initialize(double& value, bool randomize, vector<vector< string > > existingBins, vector<string>);
    bool update(double&); //returns whether list changed and MCC
    vector<double> getStats( long long&,  long long&,  long long&,  long long&);
    vector< vector<double> >  getFitStats( long long&,  long long&,  long long&,  long long&);  //results[0][...] ref results, results[1][...] Users results, results[2][...] ref+users results
    ListVector* getList();
    ListVector* getList(set<string>&);
    
protected:
    MothurOut* m;
    Utils util;
    map<int, int> seqBin; //sequence# -> bin#
    OptiMatrix* matrix;
    vector<int> randomizeSeqs;
    vector< vector<int> > bins; //bin[0] -> seqs in bin[0]
    map<int, string> binLabels; //for fitting - maps binNumber to existing reference label
    //set<string> immovableNames;
    //set<int> namesSeqs; //matrix indexes for movable sequences
    ClusterMetric* metric;
    long long truePositives, trueNegatives, falsePositives, falseNegatives, numSeqs, insertLocation, numSingletons;
    long long reftruePositives, reftrueNegatives, reffalsePositives, reffalseNegatives, numRefSeqs, numRefSingletons;
    bool removeTrainers, fitCalc;
    
    int findInsert();
    vector<long long> getCloseFarCounts(int seq, int newBin);
    
};


#endif /* opti_h */
