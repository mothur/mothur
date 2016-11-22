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
int OptiCluster::initialize(double& value, bool randomize, string initialize) {
    try {
        numSeqs = matrix->getNumSeqs();
        truePositives = 0;
        falsePositives = 0;
        falseNegatives = 0;
        trueNegatives = 0;
        
        bins.resize(numSeqs); //place seqs in own bin
        
        vector<int> temp;
        bins.push_back(temp);
        seqBin[numSeqs] = -1;
        insertLocation = numSeqs;
        
        if (initialize == "singleton") {
            
            //put everyone in own bin
            for (int i = 0; i < numSeqs; i++) { bins[i].push_back(i); }
            
            //maps randomized sequences to bins
            for (int i = 0; i < numSeqs; i++) {
                seqBin[i] = bins[i][0];
                randomizeSeqs.push_back(i);
            }
            
            if (randomize) { random_shuffle(randomizeSeqs.begin(), randomizeSeqs.end()); }
            
            //for each sequence (singletons removed on read)
            for (map<int, int>::iterator it = seqBin.begin(); it != seqBin.end(); it++) {
                if (it->second == -1) { }
                else {
                    long long numCloseSeqs = (matrix->getCloseSeqs(it->first)).size(); //does not include self
                    falseNegatives += numCloseSeqs;
                }
            }
            falseNegatives /= 2; //square matrix
            trueNegatives = numSeqs * (numSeqs-1)/2 - (falsePositives + falseNegatives + truePositives); //since everyone is a singleton no one clusters together. True negative = num far apart
            totalPairs = trueNegatives + truePositives + falseNegatives + falsePositives;
        }else {
            
            //put everyone in first bin
            for (int i = 0; i < numSeqs; i++) {
                bins[0].push_back(i);
                seqBin[i] = 0;
                randomizeSeqs.push_back(i);
            }
            
            if (randomize) { random_shuffle(randomizeSeqs.begin(), randomizeSeqs.end()); }
            
            //for each sequence (singletons removed on read)
            for (map<int, int>::iterator it = seqBin.begin(); it != seqBin.end(); it++) {
                if (it->second == -1) { }
                else {
                    long long numCloseSeqs = (matrix->getCloseSeqs(it->first)).size(); //does not include self
                    truePositives += numCloseSeqs;
                }
            }
            truePositives /= 2; //square matrix
            falsePositives = numSeqs * (numSeqs-1)/2 - (trueNegatives + falseNegatives + truePositives);
        }
        
        value = 0;
        if (metric == "mcc")             { value = calcMCC(truePositives, trueNegatives, falsePositives, falseNegatives);        }
        else if (metric == "sens")       { value = calcSens(truePositives, trueNegatives, falsePositives, falseNegatives);       }
        else if (metric == "spec")       { value = calcSpec(truePositives, trueNegatives, falsePositives, falseNegatives);       }
        else if (metric == "tptn")       { value = calcTPTN(truePositives, trueNegatives, falsePositives, falseNegatives);       }
        else if (metric == "tp")         { value = calcTP(truePositives, trueNegatives, falsePositives, falseNegatives);         }
        else if (metric == "tn")         { value = calcTN(truePositives, trueNegatives, falsePositives, falseNegatives);         }
        else if (metric == "fp")         { value = calcFP(truePositives, trueNegatives, falsePositives, falseNegatives);         }
        else if (metric == "fn")         { value = calcFN(truePositives, trueNegatives, falsePositives, falseNegatives);         }
        else if (metric == "fpfn")       { value = calcFPFN(truePositives, trueNegatives, falsePositives, falseNegatives);       }
        else if (metric == "f1score")    { value = calcF1Score(truePositives, trueNegatives, falsePositives, falseNegatives);    }
        else if (metric == "accuracy")   { value = calcAccuracy(truePositives, trueNegatives, falsePositives, falseNegatives);   }
        else if (metric == "ppv")        { value = calcPPV(truePositives, trueNegatives, falsePositives, falseNegatives);        }
        else if (metric == "npv")        { value = calcNPV(truePositives, trueNegatives, falsePositives, falseNegatives);        }
        else if (metric == "fdr")        { value = calcFDR(truePositives, trueNegatives, falsePositives, falseNegatives);        }
       
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
                
                long long tn, tp, fp, fn;
                double bestMetric = -1;
                int bestBin, bestTp, bestTn, bestFn, bestFp;
                tn = trueNegatives; tp = truePositives; fp = falsePositives; fn = falseNegatives;
                
                //this calculation is used to save time in the move and adjust function.
                //we don't need to calculate the cost of moving out of our current bin each time we
                //test moving into a new bin. Only calc this once per iteration.
                long long currentBinTP, currentBinTN, currentBinFP, currentBinFN;
                currentBinFN = falseNegatives; currentBinFP = falsePositives; currentBinTN = trueNegatives; currentBinTP = truePositives;
                
                long long cCount = 0;  long long fCount = 0;
                for (int i = 0; i < bins[binNumber].size(); i++) { //how many close sequences are in the old bin?
                    if (seqNumber == bins[binNumber][i]) {}
                    else if (!matrix->isClose(seqNumber, bins[binNumber][i])) {  fCount++;   }
                    else { cCount++;  }
                }
                //move out of old bin
                currentBinFN+=cCount; currentBinTN+=fCount; currentBinFP-=fCount; currentBinTP-=cCount;
                
                //metric in current bin
                bestMetric = calcScoreCurrentBin(tp, tn, fp, fn); bestBin = binNumber; bestTp = tp; bestTn = tn; bestFp = fp; bestFn = fn;
                
                //if not already singleton, then calc value if singleton was created
                if (!((bins[binNumber].size()) == 1)) {
                    //make a singleton
                    double singleMetric = moveAdjustTFValues(binNumber, seqNumber, -1, currentBinTP, currentBinTN, currentBinFP, currentBinFN);
                    bestBin = -1; bestTp = tp; bestTn = tn; bestFp = fp; bestFn = fn;
                    bestMetric = singleMetric;
                }
                
                set<int> binsToTry; //find bins to search eliminating dups
                for (int j = 0; j < (matrix->getCloseSeqs(seqNumber)).size(); j++) {  binsToTry.insert(seqBin[matrix->getCloseSeqs(seqNumber)[j]]); }
                
                //merge into each "close" otu
                for (set<int>::iterator it = binsToTry.begin(); it != binsToTry.end(); it++) {
                    tn = trueNegatives; tp = truePositives; fp = falsePositives; fn = falseNegatives;
                    double newMetric = moveAdjustTFValues(binNumber, seqNumber, *it, currentBinTP, currentBinTN, currentBinFP, currentBinFN);
                    
                    //new best
                    if (newMetric > bestMetric) { bestMetric = newMetric; bestBin = (*it); bestTp = tp; bestTn = tn; bestFp = fp; bestFn = fn; }
                }
                
                bool usedInsert = false;
                if (bestBin == -1) {  bestBin = insertLocation;  usedInsert = true;  }
                
                if (bestBin != binNumber) {
                    truePositives = bestTp; trueNegatives = bestTn; falsePositives = bestFp; falseNegatives = bestFn;
            
                    //move seq from i to j
                    bins[bestBin].push_back(seqNumber); //add seq to bestbin
                    bins[binNumber].erase(remove(bins[binNumber].begin(), bins[binNumber].end(), seqNumber), bins[binNumber].end()); //remove from old bin i
                }
                
                if (usedInsert) {
                    if (bins[binNumber].size() == 0) { insertLocation = binNumber;  } //set flag if old bin is empty.
                    insertLocation = findInsert();
                }
                
                //update seqBins
                seqBin[seqNumber] = bestBin; //set new OTU location
            }
        }
        
        if (metric == "mcc")             { listMetric = calcMCC(truePositives, trueNegatives, falsePositives, falseNegatives);       }
        else if (metric == "sens")       { listMetric = calcSens(truePositives, trueNegatives, falsePositives, falseNegatives);      }
        else if (metric == "spec")       { listMetric = calcSpec(truePositives, trueNegatives, falsePositives, falseNegatives);      }
        else if (metric == "tptn")       { listMetric = calcTPTN(truePositives, trueNegatives, falsePositives, falseNegatives);      }
        else if (metric == "tp")         { listMetric = calcTP(truePositives, trueNegatives, falsePositives, falseNegatives);        }
        else if (metric == "tn")         { listMetric = calcTN(truePositives, trueNegatives, falsePositives, falseNegatives);        }
        else if (metric == "fp")         { listMetric = calcFP(truePositives, trueNegatives, falsePositives, falseNegatives);        }
        else if (metric == "fn")         { listMetric = calcFN(truePositives, trueNegatives, falsePositives, falseNegatives);        }
        else if (metric == "f1score")    { listMetric = calcF1Score(truePositives, trueNegatives, falsePositives, falseNegatives);   }
        else if (metric == "accuracy")   { listMetric = calcAccuracy(truePositives, trueNegatives, falsePositives, falseNegatives);  }
        else if (metric == "ppv")        { listMetric = calcPPV(truePositives, trueNegatives, falsePositives, falseNegatives);       }
        else if (metric == "npv")        { listMetric = calcNPV(truePositives, trueNegatives, falsePositives, falseNegatives);       }
        else if (metric == "fdr")        { listMetric = calcFDR(truePositives, trueNegatives, falsePositives, falseNegatives);       }
        else if (metric == "fpfn")       { listMetric = calcFPFN(truePositives, trueNegatives, falsePositives, falseNegatives);      }
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "update");
        exit(1);
    }
}
/***********************************************************************/
double OptiCluster::moveAdjustTFValues(int bin, int seq, int newBin,  long long& tp,  long long& tn,  long long& fp,  long long& fn) {
    try {
        
        //making a singleton bin. Close but we are forcing apart.
        if (newBin == -1) {
        }else { //merging a bin
            long long ncCount = 0;  long long nfCount = 0;
            for (int i = 0; i < bins[newBin].size(); i++) { //how many close sequences are in the old bin?
                if (seq == bins[newBin][i]) {}
                else if (!matrix->isClose(seq, bins[newBin][i])) { nfCount++; }
                else { ncCount++;  }
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
        }else if (metric == "tp") {
            result =  calcTP(tp, tn, fp, fn);
        }else if (metric == "tn") {
            result =  calcTN(tp, tn, fp, fn);
        }else if (metric == "fp") {
            result =  calcFP(tp, tn, fp, fn);
        }else if (metric == "fn") {
            result =  calcFN(tp, tn, fp, fn);
        }else if (metric == "f1score") {
            result =  calcF1Score(tp, tn, fp, fn);
        }else if (metric == "accuracy") {
            result =  calcAccuracy(tp, tn, fp, fn);
        }else if (metric == "ppv") {
            result =  calcPPV(tp, tn, fp, fn);
        }else if (metric == "npv") {
            result =  calcNPV(tp, tn, fp, fn);
        }else if (metric == "fdr") {
            result =  calcFDR(tp, tn, fp, fn);
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
double OptiCluster::calcScoreCurrentBin( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
        
        if (metric == "mcc") {
            return calcMCC(tp, tn, fp, fn);
        }else if (metric == "sens") {
            return calcSens(tp, tn, fp, fn);
        }else if (metric == "spec") {
            return calcSpec(tp, tn, fp, fn);
        }else if (metric == "tptn") {
            return calcTPTN(tp, tn, fp, fn);
        }else if (metric == "tp") {
            return calcTP(tp, tn, fp, fn);
        }else if (metric == "tn") {
            return calcTN(tp, tn, fp, fn);
        }else if (metric == "fp") {
            return calcFP(tp, tn, fp, fn);
        }else if (metric == "fn") {
            return calcFN(tp, tn, fp, fn);
        }else if (metric == "f1score") {
            return calcF1Score(tp, tn, fp, fn);
        }else if (metric == "accuracy") {
            return calcAccuracy(tp, tn, fp, fn);
        }else if (metric == "ppv") {
            return calcPPV(tp, tn, fp, fn);
        }else if (metric == "npv") {
            return calcNPV(tp, tn, fp, fn);
        }else if (metric == "fdr") {
            return calcFDR(tp, tn, fp, fn);
        }else if (metric == "fpfn") {
            return calcFPFN(tp, tn, fp, fn);
        }
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcScoreCurrentBin");
        exit(1);
    }
}

/***********************************************************************/
double OptiCluster::calcMCC( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
        
        long long p = tp + fn;
        long long n = fp + tn;
        long long pPrime = tp + fp;
        long long nPrime = tn + fn;
        
        double matthewsCorrCoef = ((tp * tn) - (fp * fn)) / (double) sqrt(p * n * pPrime * nPrime);
        if(p == 0 || n == 0 || pPrime == 0 || nPrime == 0){	matthewsCorrCoef = 0;	}
        
        return matthewsCorrCoef;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcMCC");
        exit(1);
    }
}
/***********************************************************************/
double OptiCluster::calcSens( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
        
         long long p = tp + fn;
        double sensitivity = tp / (double) p;
        
        if(p == 0)	{	sensitivity = 0;	}
        
        return sensitivity;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcSens");
        exit(1);
    }
}
/***********************************************************************/
double OptiCluster::calcSpec( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
         long long n = fp + tn;
        double specificity = tn / (double) n;
        
        if(n == 0)			{	specificity = 0;	}
        
        return specificity;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcSpec");
        exit(1);
    }
}
/***********************************************************************/
double OptiCluster::calcTPTN( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
        
         long long p = tp + tn;
        
        double tptn = p / (double)(tp + tn + fp + fn);
        
        return tptn;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcTPTN");
        exit(1);
    }
}
/***********************************************************************/
double OptiCluster::calcTP( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
        
        double tpmax = tp / (double)(tp + tn + fp + fn);
        
        return tpmax;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcTP");
        exit(1);
    }
}
/***********************************************************************/
double OptiCluster::calcTN( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
        
        double tnmax = tn / (double)(tp + tn + fp + fn);
        
        return tnmax;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcTN");
        exit(1);
    }
}
/***********************************************************************/
double OptiCluster::calcFP( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
        
        double fpmin = fp / (double)(tp + tn + fp + fn);
        
        return (1.0 - fpmin);
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcFP");
        exit(1);
    }
}
/***********************************************************************/
double OptiCluster::calcFN( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
        
        double fnmin = fn / (double)(tp + tn + fp + fn);
        
        return (1.0 - fnmin);
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcFN");
        exit(1);
    }
}

/***********************************************************************/
double OptiCluster::calcFPFN( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
        
         long long p = fp + fn;
        
        double fpfn = 1.0 - (p / (double)(tp + tn + fp + fn)); //minimize
        
        return fpfn;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcFPFN");
        exit(1);
    }
}
/***********************************************************************/
double OptiCluster::calcF1Score( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
        
         long long p = tp + fn;
         long long pPrime = tp + fp;
        double f1Score = 2.0 * tp / (double) (p + pPrime);
        
        if(p + pPrime == 0)	{	f1Score = 0;	}
        
        return f1Score;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcF1Score");
        exit(1);
    }
}
/***********************************************************************/
double OptiCluster::calcAccuracy( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
        
         long long p = tp + fn;
         long long n = fp + tn;
        double accuracy = (tp + tn) / (double) (p + n);
        if(p + n == 0)		{	accuracy = 0;								}
        return accuracy;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcAccuracy");
        exit(1);
    }
}
/***********************************************************************/
double OptiCluster::calcPPV( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
        
         long long pPrime = tp + fp;
        double positivePredictiveValue = tp / (double) pPrime;
        
        if(pPrime == 0)		{	positivePredictiveValue = 0;		}
        
        return positivePredictiveValue;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcPPV");
        exit(1);
    }
}
/***********************************************************************/
double OptiCluster::calcNPV( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
        
         long long nPrime = tn + fn;
        double negativePredictiveValue = tn / (double) nPrime;
        
        if(nPrime == 0)		{	negativePredictiveValue = 0;		}
        
        return negativePredictiveValue;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcNPV");
        exit(1);
    }
}
/***********************************************************************/
double OptiCluster::calcFDR( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
        
         long long pPrime = tp + fp;
        double falseDiscoveryRate = fp / (double) pPrime;
        
        if(pPrime == 0)		{	falseDiscoveryRate = 0;		}
        
        return (1.0-falseDiscoveryRate);
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "calcFDR");
        exit(1);
    }
}
/***********************************************************************/
vector<double> OptiCluster::getStats( long long& tp,  long long& tn,  long long& fp,  long long& fn) {
    try {
        long long singletn = matrix->getNumSingletons() + numSingletons;
        long long tempnumSeqs = numSeqs + singletn;
        
        tp = truePositives;
        fp = falsePositives;
        fn = falseNegatives;
        tn = tempnumSeqs * (tempnumSeqs-1)/2 - (falsePositives + falseNegatives + truePositives); //adds singletons to tn
        
         long long p = tp + fn;
         long long n = fp + tn;
         long long pPrime = tp + fp;
         long long nPrime = tn + fn;
        
        double sensitivity = tp /(double) p;
        double specificity = tn / (double)n;
        double positivePredictiveValue = tp / (double)pPrime;
        double negativePredictiveValue = tn / (double)nPrime;
        double falseDiscoveryRate = fp / (double)pPrime;
        
        double accuracy = (tp + tn) / (double)(p + n);
        double matthewsCorrCoef = (tp * tn - fp * fn) / (double) sqrt(p * n * pPrime * nPrime);	if(p == 0 || n == 0){	matthewsCorrCoef = 0;	}
        double f1Score = 2.0 * tp / (double)(p + pPrime);
        
        
        if(p == 0)			{	sensitivity = 0;	matthewsCorrCoef = 0;	}
        if(n == 0)			{	specificity = 0;	matthewsCorrCoef = 0;	}
        if(p + n == 0)		{	accuracy = 0;								}
        if(p + pPrime == 0)	{	f1Score = 0;								}
        if(pPrime == 0)		{	positivePredictiveValue = 0;	falseDiscoveryRate = 0;	matthewsCorrCoef = 0;	}
        if(nPrime == 0)		{	negativePredictiveValue = 0;	matthewsCorrCoef = 0;							}
        
        vector<double> results;
        
        //results.push_back(truePositives); results.push_back(tn); results.push_back(falsePositives); results.push_back(falseNegatives);
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
                string otu = matrix->getName(bins[i][0]);
                
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
long long OptiCluster::getNumBins() {
    try {
        long long singletn = matrix->getNumSingletons();
        
        for (int i = 0; i < bins.size(); i++) {
            if (bins[i].size() != 0) {
                singletn++;
            }
        }
        
        return singletn;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "getNumBins");
        exit(1);
    }
}

/***********************************************************************/
int OptiCluster::findInsert() {
    try {
        
        //for each sequence (singletons removed on read)
        for (int i = 0; i < bins.size(); i++) {
            
            if (m->control_pressed) { break; }
            
            if (bins[i].size() == 0) { return i;  }
        }
        
        return -1;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "getList");
        exit(1);
    }
}

/***********************************************************************/
