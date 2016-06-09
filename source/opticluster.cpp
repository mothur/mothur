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
int OptiCluster::initialize(double& value) {
    try {
        numSeqs = matrix->getNumSeqs();
        
        bins.resize(numSeqs); //place seqs in own bin
        
        for (int i = 0; i < numSeqs; i++) { bins[i].push_back(i); }
        
        //maps randomized sequences to bins
        for (int i = 0; i < numSeqs; i++) {
            seqBin[i] = bins[i][0];
            randomizeSeqs.push_back(i);
        }
        
        random_shuffle(randomizeSeqs.begin(), randomizeSeqs.end());
        
        vector<int> temp;
        bins.push_back(temp);
        seqBin[numSeqs] = -1;
        insertLocation = numSeqs;
    
        truePositives = 0;
        falsePositives = 0;
        falseNegatives = 0;
        trueNegatives = 0;
        //for each sequence (singletons removed on read)
        for (map<int, int>::iterator it = seqBin.begin(); it != seqBin.end(); it++) {
            if (it->second == -1) { }
            else {
                int numCloseSeqs = (matrix->getCloseSeqs(it->first)).size(); //does not include self
                falseNegatives += numCloseSeqs;
            }
        }
        falseNegatives /= 2; //square matrix
        trueNegatives = numSeqs * (numSeqs-1)/2 - (falsePositives + falseNegatives + truePositives); //since everyone is a singleton no one clusters together. True negative = num far apart
        totalPairs = trueNegatives + truePositives + falseNegatives + falsePositives;
        
        value = 0;
        if (metric == "mcc")        { value = calcMCC(truePositives, trueNegatives, falsePositives, falseNegatives);    }
        else if (metric == "sens")  { value = calcSens(truePositives, trueNegatives, falsePositives, falseNegatives);   }
        else if (metric == "spec")  { value = calcSpec(truePositives, trueNegatives, falsePositives, falseNegatives);   }
        else if (metric == "tptn")  { value = calcTPTN(truePositives, trueNegatives, falsePositives, falseNegatives);   }
        else if (metric == "tp2tn") { value = calcTP2TN(truePositives, trueNegatives, falsePositives, falseNegatives);  }
        else if (metric == "fpfn")  { value = calcFPFN(truePositives, trueNegatives, falsePositives, falseNegatives);   }
       
        return value;
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
bool OptiCluster::update(double& listMetric) {
    try {

        //for each sequence (singletons removed on read)
        for (int i = 0; i < randomizeSeqs.size(); i++) {
            
            if (m->control_pressed) { break; }
            
            map<int, int>::iterator it = seqBin.find(randomizeSeqs[i]);
            
            int seqNumber = it->first;
            int binNumber = it->second;
            
            if (binNumber == -1) { }
            else {
                
                double tn, tp, fp, fn;
                double bestMetric = -1;
                vector< vector<int> > bestMetricsTPValues;
                //already singleton
                if ((bins[binNumber].size()) == 1) {
                    tn = trueNegatives; tp = truePositives; fp = falsePositives; fn = falseNegatives;
                    vector<int> temp;
                    double singleMetric;
                    if (metric == "mcc") {
                        singleMetric = calcMCC(tp, tn, fp, fn);
                        temp.push_back(binNumber); temp.push_back(tp); temp.push_back(tn); temp.push_back(fp); temp.push_back(fn);
                    }else if (metric == "sens") {
                        singleMetric = calcSens(tp, tn, fp, fn);
                        temp.push_back(binNumber); temp.push_back(tp); temp.push_back(tn); temp.push_back(fp); temp.push_back(fn);
                    }else if (metric == "spec") {
                        singleMetric = calcSpec(tp, tn, fp, fn);
                        temp.push_back(binNumber); temp.push_back(tp); temp.push_back(tn); temp.push_back(fp); temp.push_back(fn);
                    }else if (metric == "tptn") {
                        singleMetric = calcTPTN(tp, tn, fp, fn);
                        temp.push_back(binNumber); temp.push_back(tp); temp.push_back(tn); temp.push_back(fp); temp.push_back(fn);
                    }else if (metric == "tp2tn") {
                        singleMetric = calcTP2TN(tp, tn, fp, fn);
                        temp.push_back(binNumber); temp.push_back(tp); temp.push_back(tn); temp.push_back(fp); temp.push_back(fn);
                    }else if (metric == "fpfn") {
                        singleMetric = calcFPFN(tp, tn, fp, fn);
                        temp.push_back(binNumber); temp.push_back(tp); temp.push_back(tn); temp.push_back(fp); temp.push_back(fn);
                    }
                    bestMetric = singleMetric;
                    bestMetricsTPValues.push_back(temp);
                }else {
                    //make a singleton
                    tn = trueNegatives; tp = truePositives; fp = falsePositives; fn = falseNegatives;
                    double singleMetric = moveAdjustTFValues(binNumber, seqNumber, -1, tp, tn, fp, fn);
                    vector<int> temp;
                    temp.push_back(-1); temp.push_back(tp); temp.push_back(tn); temp.push_back(fp); temp.push_back(fn);
                    bestMetric = singleMetric;
                    bestMetricsTPValues.push_back(temp);
                }
                
                
                set<int> binsToTry; //find bins to search eliminating dups
                for (int j = 0; j < (matrix->getCloseSeqs(seqNumber)).size(); j++) { binsToTry.insert(seqBin[matrix->getCloseSeqs(seqNumber)[j]]); }
                
                //merge into each "close" otu
                for (set<int>::iterator it = binsToTry.begin(); it != binsToTry.end(); it++) {
                    tn = trueNegatives; tp = truePositives; fp = falsePositives; fn = falseNegatives;
                    double newMetric = moveAdjustTFValues(binNumber, seqNumber, *it, tp, tn, fp, fn);
                    vector<int> temp;
                    temp.push_back(*it); temp.push_back(tp); temp.push_back(tn); temp.push_back(fp); temp.push_back(fn);
                    if (newMetric > bestMetric) { //new best
                        bestMetric = newMetric;
                        bestMetricsTPValues.clear();
                        bestMetricsTPValues.push_back(temp);
                    }else if (newMetric == bestMetric) { //tie
                        bestMetricsTPValues.push_back(temp);
                    }
                }
                
                //choose "best" otu - which is highest value or last item in map. - shuffle for ties
                random_shuffle(bestMetricsTPValues.begin(), bestMetricsTPValues.end());
                int bestBin = bestMetricsTPValues[0][0];
                truePositives = bestMetricsTPValues[0][1];
                trueNegatives = bestMetricsTPValues[0][2];
                falsePositives = bestMetricsTPValues[0][3];
                falseNegatives = bestMetricsTPValues[0][4];
                
                if (bestBin == -1) {  bestBin = insertLocation;  }
                
                if (bestBin != binNumber) {
                    //move seq from i to j
                    bins[bestBin].push_back(seqNumber); //add seq to bestbin
                    bins[binNumber].erase(remove(bins[binNumber].begin(), bins[binNumber].end(), seqNumber), bins[binNumber].end()); //remove from old bin i
                }
                if (bins[binNumber].size() == 0) { insertLocation = binNumber;  } //set flag if old bin is empty.
                
                //update seqBins
                seqBin[seqNumber] = bestBin; //set new OTU location
            }
        }
        
        if (metric == "mcc")        { listMetric = calcMCC(truePositives, trueNegatives, falsePositives, falseNegatives);    }
        else if (metric == "sens")  { listMetric = calcSens(truePositives, trueNegatives, falsePositives, falseNegatives);   }
        else if (metric == "spec")  { listMetric = calcSpec(truePositives, trueNegatives, falsePositives, falseNegatives);   }
        else if (metric == "tptn")  { listMetric = calcTPTN(truePositives, trueNegatives, falsePositives, falseNegatives);   }
        else if (metric == "tp2tn") { listMetric = calcTP2TN(truePositives, trueNegatives, falsePositives, falseNegatives);  }
        else if (metric == "fpfn")  { listMetric = calcFPFN(truePositives, trueNegatives, falsePositives, falseNegatives);   }
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "update");
        exit(1);
    }
}
/***********************************************************************/
double OptiCluster::moveAdjustTFValues(int bin, int seq, int newBin, double& tp, double& tn, double& fp, double& fn) {
    try {
        if (bin == newBin) {
            if (metric == "mcc") {
                return calcMCC(tp, tn, fp, fn);
            }else if (metric == "sens") {
                return calcSens(tp, tn, fp, fn);
            }else if (metric == "spec") {
                return calcSpec(tp, tn, fp, fn);
            }else if (metric == "tptn") {
                return calcTPTN(tp, tn, fp, fn);
            }else if (metric == "tp2tn") {
                return calcTP2TN(tp, tn, fp, fn);
            }else if (metric == "fpfn") {
                return calcFPFN(tp, tn, fp, fn);
            }
        }
 
        int cCount = 0; int fCount = 0;
        for (int i = 0; i < bins[bin].size(); i++) { //how many close sequences are in the old bin?
            if (seq == bins[bin][i]) {}
            else if (matrix->isClose(seq, bins[bin][i])) {  cCount++;   }
            else { fCount++;  }
        }

        //move out of old bin
        fn+=cCount; tn+=fCount; fp-=fCount; tp-=cCount;
        
        //making a singleton bin. Close but we are forcing apart.
        if (newBin == -1) {
        }else { //merging a bin
            int ncCount = 0; int nfCount = 0;
            for (int i = 0; i < bins[newBin].size(); i++) { //how many close sequences are in the old bin?
                if (seq == bins[newBin][i]) {}
                else if (matrix->isClose(seq, bins[newBin][i])) { ncCount++; }
                else { nfCount++;  }
            }
   
            //move into new bin
            fn-=ncCount; tn-=nfCount;  tp+=ncCount; fp+=nfCount;
        }

        double result = 0.0;
        if (metric == "mcc") {
            result =  calcMCC(tp, tn, fp, fn);
        }else if (metric == "sens") {
            result =  calcSens(tp, tn, fp, fn);
        }else if (metric == "spec") {
            result =  calcSpec(tp, tn, fp, fn);
        }else if (metric == "tptn") {
            result =  calcTPTN(tp, tn, fp, fn);
        }else if (metric == "tp2tn") {
            result =  calcTP2TN(tp, tn, fp, fn);
        }else if (metric == "fpfn") {
            result =  calcFPFN(tp, tn, fp, fn);
        }

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
        
        double matthewsCorrCoef = ((tp * tn) - (fp * fn)) / sqrt(p * n * pPrime * nPrime);
        if(p == 0 || n == 0 || pPrime == 0 || nPrime == 0){	matthewsCorrCoef = 0;	}
        
        return matthewsCorrCoef;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcMCC");
        exit(1);
    }
}
/***********************************************************************/
double OptiCluster::calcSens(double tp, double tn, double fp, double fn) {
    try {
        
        double p = tp + fn;
        double sensitivity = tp / p;
        
        if(p == 0)	{	sensitivity = 0;	}
        
        return sensitivity;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcSens");
        exit(1);
    }
}
/***********************************************************************/
double OptiCluster::calcSpec(double tp, double tn, double fp, double fn) {
    try {
        double n = fp + tn;
        double specificity = tn / n;
        
        if(n == 0)			{	specificity = 0;	}
        
        return specificity;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcSpec");
        exit(1);
    }
}
/***********************************************************************/
double OptiCluster::calcTPTN(double tp, double tn, double fp, double fn) {
    try {
        
        double p = tp + tn;
        
        double tptn = p / (double)(tp + tn + fp + fn);
        
        return tptn;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcTPTN");
        exit(1);
    }
}
/***********************************************************************/
double OptiCluster::calcTP2TN(double tp, double tn, double fp, double fn) {
    try {
        
        double p = tp + (2*tn);
        
        double tptn = p / (double)(tp + tn + fp + fn);
        
        return tptn;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcTP2TN");
        exit(1);
    }
}

/***********************************************************************/
double OptiCluster::calcFPFN(double tp, double tn, double fp, double fn) {
    try {
        
        double p = fp + fn;
        
        double tptn = 1.0 - (p / (double)(tp + tn + fp + fn)); //minimize
        
        return tptn;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcFPFN");
        exit(1);
    }
}
/***********************************************************************/
vector<double> OptiCluster::getStats() {
    try {
        
        double tp = (double) truePositives;
        double fp = (double) falsePositives;
        double tn = (double) trueNegatives;
        double fn = (double) falseNegatives;
        
        double p = tp + fn;
        double n = fp + tn;
        double pPrime = tp + fp;
        double nPrime = tn + fn;
        
        double sensitivity = tp / p;
        double specificity = tn / n;
        double positivePredictiveValue = tp / pPrime;
        double negativePredictiveValue = tn / nPrime;
        double falseDiscoveryRate = fp / pPrime;
        
        double accuracy = (tp + tn) / (p + n);
        double matthewsCorrCoef = (tp * tn - fp * fn) / sqrt(p * n * pPrime * nPrime);	if(p == 0 || n == 0){	matthewsCorrCoef = 0;	}
        double f1Score = 2.0 * tp / (p + pPrime);
        
        
        if(p == 0)			{	sensitivity = 0;	matthewsCorrCoef = 0;	}
        if(n == 0)			{	specificity = 0;	matthewsCorrCoef = 0;	}
        if(p + n == 0)		{	accuracy = 0;								}
        if(p + pPrime == 0)	{	f1Score = 0;								}
        if(pPrime == 0)		{	positivePredictiveValue = 0;	falseDiscoveryRate = 0;	matthewsCorrCoef = 0;	}
        if(nPrime == 0)		{	negativePredictiveValue = 0;	matthewsCorrCoef = 0;							}
        
        vector<double> results;
        
        results.push_back(truePositives); results.push_back(trueNegatives); results.push_back(falsePositives); results.push_back(falseNegatives);
        results.push_back(sensitivity); results.push_back(specificity); results.push_back(positivePredictiveValue); results.push_back(negativePredictiveValue);
        results.push_back(falseDiscoveryRate); results.push_back(accuracy); results.push_back(matthewsCorrCoef); results.push_back(f1Score);
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "getStats");
        exit(1);
    }
}
/***********************************************************************/
ListVector* OptiCluster::getList() {
    try {
        ListVector* list = new ListVector();
        ListVector* singleton = matrix->getListSingle();
        
        if (singleton != NULL) { //add in any sequences above cutoff in read. Removing these saves clustering time.
            for (int i = 0; i < singleton->getNumBins(); i++) {
                if (singleton->get(i) != "") {
                    list->push_back(singleton->get(i));
                }
            }
            delete singleton;
        }

        for (int i = 0; i < bins.size(); i++) {
            if (bins[i].size() != 0) {
                string otu = "";
                otu += matrix->getName(bins[i][0]);
                
                for (int j = 1; j < bins[i].size(); j++) {
                    otu += "," + matrix->getName(bins[i][j]);
                }
                list->push_back(otu);
            }
        }
        return list;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "getList");
        exit(1);
    }
}
/***********************************************************************/
