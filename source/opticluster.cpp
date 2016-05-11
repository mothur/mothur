//
//  opticluster.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/20/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#include "opticluster.h"

/***********************************************************************/
//randomly assign sequences to OTUs
int OptiCluster::initialize() {
    try {
        numSeqs = matrix->getNumSeqs();
        
        bins.resize(numSeqs); //place seqs in own bin
        
        for (int i = 1; i < numSeqs; i++) { bins[i].push_back(i); }
        
        random_shuffle(bins.begin(), bins.end());
        
        //maps randomized sequences to bins
        for (int i = 0; i < bins[0].size(); i++) {
            seqBin[i] = bins[i][0];
        }
        
        bins[numSeqs].push_back(-1);
        seqBin[numSeqs] = -1;
        insertLocation = numSeqs;
    
        if (metric == "mcc") {
            truePositives = 0;
            falsePositives = 0;
            falseNegatives = numSeqs;
            trueNegatives = numSeqs * (numSeqs-1)/2 - (falsePositives + falseNegatives + truePositives);
            
            listVectorMetric = calcMCC(truePositives, trueNegatives, falsePositives, falseNegatives);
        }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "initialize");
        exit(1);
    }
}
/***********************************************************************/
/* for each sequence with mutual information (close)
* remove from current OTU and calculate MCC when sequence forms its own OTU or joins one of the other OTUs where there is a sequence within the `threshold` (no need to calculate MCC if the paired sequence is already in same OTU and no need to try every OTU - just those where there's a close sequence) 
 * keep or move the sequence to the OTU where the `metric` is the largest - flip a coin on ties */
bool OptiCluster::update(double& bestBinMetric) {
    try {
        //int iters = 0;
        //while ((listVectorMetric < stableMetric) && (iters < maxIters)) {
        
            //for each sequence (singletons removed on read)
            for (int i = 0; i < bins.size(); i++) {
                if ( m->control_pressed) { break; }
                
                map<double, int> otuMetric; //maps otu to metric to find "best" otu
                
                //remove from old otu and create singleton
                map<int, int>::iterator it = seqBin.find(i);
                
                if ((bins[it->second].size()) == 1) {} //already singleton
                else {
                    double singleMetric = moveAdjustTFValues(it->second, i, -1);
                    otuMetric[singleMetric] = -1;
                }
                
                //merge into each "close" otu
                for (int j = 0; j < (matrix->getCloseSeqs(i)).size(); j++) {
                    double newMetric = moveAdjustTFValues(it->second, i, j);
                    otuMetric[newMetric] = it->second;
                }
                
                //choose "best" otu - which is highest value or last item in map.
                bestBinMetric = (otuMetric.rend())->first;
                int bestBin = (otuMetric.rend())->second;
                
                if (bestBin == -1) {  bestBin = insertLocation;  }
                
                //move seq from i to j
                bins[bestBin].push_back(i); //add seq to bestbin
                bins[i].erase(remove(bins[i].begin(), bins[i].end(), bestBin), bins[i].end()); //remove from old bin i
                
                it = seqBin.find(i); //find Otu to be plucked from
                
                if (bins[it->second].size() == 0) { seqBin[it->second] = -1; insertLocation = it->second; } //set flag if old bin is empty.

                //update seqBins
                seqBin[i] = bestBin; //set new OTU location
            }
       // }
        
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "update");
        exit(1);
    }
}
/***********************************************************************/
double OptiCluster::moveAdjustTFValues(int bin, int seq, int newBin) {
    try {
        
        vector<int> seqs; seqs = bins[bin];
        int closeCount = 0;
        for (int i = 0; i < seqs.size(); i++) { //how many close sequences are in the old bin?
            if (matrix->isClose(seq, seqs[i])) {  closeCount++;  }
        }
        
        //move out of old bin into singleton
        int farApart = (seqs.size()-closeCount);
        falsePositives-=farApart; truePositives-=closeCount; falseNegatives+=closeCount;
        
        double result = 0.0;
        
        //making a singleton bin. Close but we are forcing apart.
        if (newBin == -1) { }
        else { //merging a bin
            vector<int> newSeqs; newSeqs = bins[newBin];
            int newCloseCount = 0;
            for (int i = 0; i < newSeqs.size(); i++) { //how many close sequences are in the old bin?
                if (matrix->isClose(seq, newSeqs[i])) {  newCloseCount++;  }
            }
            farApart = (newSeqs.size()-newCloseCount);
            truePositives+=closeCount; falsePositives+=farApart; falseNegatives-=closeCount;
        }
        
        trueNegatives = numSeqs * (numSeqs-1)/2 - (falsePositives + falseNegatives + truePositives);
        
        if (metric == "mcc") { result = calcMCC(truePositives, trueNegatives, falsePositives, falseNegatives);  }
        
        return result;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "moveAdjustTFValues");
        exit(1);
    }
}
/***********************************************************************/
double OptiCluster::calcMCC(double tp, double tn, double fp, double fn) {
    try {
        
        double p = tp + fn;
        double n = fp + tn;
        double pPrime = tp + fp;
        double nPrime = tn + fn;
        
        double matthewsCorrCoef = (tp * tn - fp * fn) / sqrt(p * n * pPrime * nPrime);
        if(p == 0 || n == 0){	matthewsCorrCoef = 0;	}
        
        return matthewsCorrCoef;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcMCC");
        exit(1);
    }
}
/***********************************************************************/
ListVector* OptiCluster::getList() {
    try {
        ListVector* list = new ListVector();
        
        
        
        
        return list;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "getList");
        exit(1);
    }
}
/***********************************************************************/
