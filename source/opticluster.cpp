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
        Utils util;
        
        if (initialize == "singleton") {
            
            //put everyone in own bin
            for (int i = 0; i < numSeqs; i++) { bins[i].push_back(i); }
            
            //maps randomized sequences to bins
            for (int i = 0; i < numSeqs; i++) {
                seqBin[i] = bins[i][0];
                randomizeSeqs.push_back(i);
            }
            
            if (randomize) { util.mothurRandomShuffle(randomizeSeqs); }
            
            //for each sequence (singletons removed on read)
            for (map<int, int>::iterator it = seqBin.begin(); it != seqBin.end(); it++) {
                if (it->second == -1) { }
                else {
                    long long numCloseSeqs = (matrix->getNumClose(it->first)); //does not include self
                    falseNegatives += numCloseSeqs;
                }
            }
            falseNegatives /= 2; //square matrix
            trueNegatives = numSeqs * (numSeqs-1)/2 - (falsePositives + falseNegatives + truePositives); //since everyone is a singleton no one clusters together. True negative = num far apart
        }else {
            
            //put everyone in first bin
            for (int i = 0; i < numSeqs; i++) {
                bins[0].push_back(i);
                seqBin[i] = 0;
                randomizeSeqs.push_back(i);
            }
            
            if (randomize) { util.mothurRandomShuffle(randomizeSeqs); }
            
            //for each sequence (singletons removed on read)
            for (map<int, int>::iterator it = seqBin.begin(); it != seqBin.end(); it++) {
                if (it->second == -1) { }
                else {
                    long long numCloseSeqs = (matrix->getNumClose(it->first)); //does not include self
                    truePositives += numCloseSeqs;
                }
            }
            truePositives /= 2; //square matrix
            falsePositives = numSeqs * (numSeqs-1)/2 - (trueNegatives + falseNegatives + truePositives);
        }
        
        value = metric->getValue(truePositives, trueNegatives, falsePositives, falseNegatives);
        
        return value;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "initialize");
        exit(1);
    }
}
/***********************************************************************/
//populate preset OTUs, randomize sequences to "fit" into these OTUs
int OptiCluster::initialize(double& value, bool randomize, vector<vector< string > > existingBins, vector<string> bls) {
    try {
        for (int i = 0; i < existingBins.size(); i++) { for (int j = 0; j < existingBins[i].size(); j++) { immovableNames.insert(existingBins[i][j]); }  }
        vector< vector< int> > translatedBins;
        vector<int> indexSeqs = matrix->getNumSeqs(existingBins, translatedBins); //number of seqs not in existingBins, the seqs to "fit" into OTUs
        for (int i = 0; i < numSeqs; i++) { namesSeqs.insert(indexSeqs[i]); }
        truePositives = 0;
        falsePositives = 0;
        falseNegatives = 0;
        trueNegatives = 0;
        numSeqs = namesSeqs.size();
        removeTrainers = true;
        fitCalc=true;
        
        Utils util;
        
        int binNumber = 0;
        int placeHolderIndex = -1;
        for (int i = 0; i < translatedBins.size(); i++) {
            if (translatedBins[i].size() != 0) { //remove any bins that contain seqeunces above cutoff, better to be a singleton than in a bin with no close seqs
                binLabels[binNumber] = bls[i];
                bins.push_back(translatedBins[i]);
                for (int j = 0; j < translatedBins[i].size(); j++) {
                    if (translatedBins[i][j] < 0) {  translatedBins[i][j] = placeHolderIndex; placeHolderIndex--; }
                    seqBin[translatedBins[i][j]] = binNumber;
                }
                binNumber++;
            }
        }
        
        
        
        //randomly assigns fit seqs to an OTU
        for (int i = 0; i < numSeqs; i++) {
            randomizeSeqs.push_back(indexSeqs[i]);
            
            set<int> binsToTry;
            set<int> closeSeqs = matrix->getCloseSeqs(indexSeqs[i]);
            for (set<int>::iterator itClose = closeSeqs.begin(); itClose != closeSeqs.end(); itClose++) { if (namesSeqs.count(*itClose) == 0) {  binsToTry.insert(seqBin[*itClose]);  } }
            
            bins[i].push_back(indexSeqs[i]);
            seqBin[indexSeqs[i]] = i;
            
            
        }
        
        if (randomize) { util.mothurRandomShuffle(randomizeSeqs); }
        
        //for each sequence (singletons removed on read)
        for (map<int, int>::iterator it = seqBin.begin(); it != seqBin.end(); it++) {
            if (it->second == -1) { }
            else {
                long long numCloseSeqs = (matrix->getNumClose(it->first)); //does not include self
                falseNegatives += numCloseSeqs;
            }
        }
        falseNegatives /= 2; //square matrix
        trueNegatives = numSeqs * (numSeqs-1)/2 - (falsePositives + falseNegatives + truePositives); //since everyone is a singleton no one clusters together. True negative = num far apart
        
        
        
        //add insert location
        seqBin[bins.size()] = -1;
        insertLocation = bins.size();
        vector<int> temp;
        bins.push_back(temp);
        
        value = metric->getValue(truePositives, trueNegatives, falsePositives, falseNegatives);
        
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
            
            if (m->getControl_pressed()) { break; }
            
            map<int, int>::iterator it = seqBin.find(randomizeSeqs[i]);
            
            int seqNumber = it->first;
            int binNumber = it->second;
            
            if (binNumber == -1) { }
            else {
                
                long long tn, tp, fp, fn;
                double bestMetric = -1;
                long long bestBin, bestTp, bestTn, bestFn, bestFp;
                tn = trueNegatives; tp = truePositives; fp = falsePositives; fn = falseNegatives;
                
                //close / far count in current bin
                vector<long long> results = getCloseFarCounts(seqNumber, binNumber);
                long long cCount = results[0];  long long fCount = results[1];
                
                //metric in current bin
                bestMetric = metric->getValue(tp, tn, fp, fn); bestBin = binNumber; bestTp = tp; bestTn = tn; bestFp = fp; bestFn = fn;
                
                //if not already singleton, then calc value if singleton was created
                if (!((bins[binNumber].size()) == 1)) {
                    //make a singleton
                    //move out of old bin
                    fn+=cCount; tn+=fCount; fp-=fCount; tp-=cCount;
                    double singleMetric = metric->getValue(tp, tn, fp, fn);
                    if (singleMetric > bestMetric) {
                        bestBin = -1; bestTp = tp; bestTn = tn; bestFp = fp; bestFn = fn;
                        bestMetric = singleMetric;
                    }
                }
                
                set<int> binsToTry;
                set<int> closeSeqs = matrix->getCloseSeqs(seqNumber);
                for (set<int>::iterator itClose = closeSeqs.begin(); itClose != closeSeqs.end(); itClose++) {
                    if (fitCalc) { //we only want to try reference bins
                        if (namesSeqs.count(*itClose) == 0) {  binsToTry.insert(seqBin[*itClose]);  }
                    }else {  binsToTry.insert(seqBin[*itClose]);  }
                }
                
                if (binsToTry.size() != 0) {
                    //merge into each "close" otu
                    for (set<int>::iterator it = binsToTry.begin(); it != binsToTry.end(); it++) {
                        tn = trueNegatives; tp = truePositives; fp = falsePositives; fn = falseNegatives;
                        fn+=cCount; tn+=fCount; fp-=fCount; tp-=cCount; //move out of old bin
                        results = getCloseFarCounts(seqNumber, *it);
                        fn-=results[0]; tn-=results[1];  tp+=results[0]; fp+=results[1]; //move into new bin
                        double newMetric = metric->getValue(tp, tn, fp, fn); //score when sequence is moved
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
                    
                    if (usedInsert) { insertLocation = findInsert(); }
                    
                    //update seqBins
                    seqBin[seqNumber] = bestBin; //set new OTU location
                }else { cout << seqNumber << " has no bins to try\n"; }
            }
        }
        
        listMetric = metric->getValue(truePositives, trueNegatives, falsePositives, falseNegatives);
        
        if (m->getDebug()) { ListVector* list = getList(); list->print(cout); delete list; }
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "update");
        exit(1);
    }
}
/***********************************************************************/
vector<long long> OptiCluster::getCloseFarCounts(int seq, int newBin) {
    try {
        vector<long long> results; results.push_back(0); results.push_back(0);
        
        if (newBin == -1) { }  //making a singleton bin. Close but we are forcing apart.
        else { //merging a bin
            for (int i = 0; i < bins[newBin].size(); i++) {
                if (seq == bins[newBin][i]) {} //ignore self
                else if (!matrix->isClose(seq, bins[newBin][i])) { results[1]++; }  //this sequence is "far away" from sequence i - above the cutoff
                else { results[0]++;  }  //this sequence is "close" to sequence i - distance between them is less than cutoff
            }
        }
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "getCloseFarCounts");
        exit(1);
    }
}

/***********************************************************************/
vector<double> OptiCluster::getStats(long long& tp,  long long& tn,  long long& fp,  long long& fn) {
    try {
        if (fitCalc) { return getFitStats(tp,tn,fp,fn); }
        
        long long singletn = matrix->getNumSingletons() + numSingletons;
        long long tempnumSeqs = numSeqs + singletn;
        
        tp = truePositives;
        fp = falsePositives;
        fn = falseNegatives;
        tn = tempnumSeqs * (tempnumSeqs-1)/2 - (falsePositives + falseNegatives + truePositives); //adds singletons to tn
        
        vector<double> results;
        
        Sensitivity sens;   double sensitivity = sens.getValue(tp, tn, fp, fn); results.push_back(sensitivity);
        Specificity spec;   double specificity = spec.getValue(tp, tn, fp, fn); results.push_back(specificity);
        PPV ppv;            double positivePredictiveValue = ppv.getValue(tp, tn, fp, fn); results.push_back(positivePredictiveValue);
        NPV npv;            double negativePredictiveValue = npv.getValue(tp, tn, fp, fn); results.push_back(negativePredictiveValue);
        FDR fdr;            double falseDiscoveryRate = fdr.getValue(tp, tn, fp, fn); results.push_back(falseDiscoveryRate);
        Accuracy acc;       double accuracy = acc.getValue(tp, tn, fp, fn); results.push_back(accuracy);
        MCC mcc;            double matthewsCorrCoef = mcc.getValue(tp, tn, fp, fn); results.push_back(matthewsCorrCoef);
        F1Score f1;         double f1Score = f1.getValue(tp, tn, fp, fn); results.push_back(f1Score);
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "getStats");
        exit(1);
    }
}
/***********************************************************************/
vector<double> OptiCluster::getFitStats(long long& tp,  long long& tn,  long long& fp,  long long& fn) {
    try {
        long long singletn = matrix->getNumFitSingletons() + numSingletons;
        long long tempnumSeqs = numSeqs + singletn;
        
        tp = truePositives;
        fp = falsePositives;
        fn = falseNegatives;
        tn = tempnumSeqs * (tempnumSeqs-1)/2 - (falsePositives + falseNegatives + truePositives); //adds singletons to tn
        
        vector<double> results;
        
        Sensitivity sens;   double sensitivity = sens.getValue(tp, tn, fp, fn); results.push_back(sensitivity);
        Specificity spec;   double specificity = spec.getValue(tp, tn, fp, fn); results.push_back(specificity);
        PPV ppv;            double positivePredictiveValue = ppv.getValue(tp, tn, fp, fn); results.push_back(positivePredictiveValue);
        NPV npv;            double negativePredictiveValue = npv.getValue(tp, tn, fp, fn); results.push_back(negativePredictiveValue);
        FDR fdr;            double falseDiscoveryRate = fdr.getValue(tp, tn, fp, fn); results.push_back(falseDiscoveryRate);
        Accuracy acc;       double accuracy = acc.getValue(tp, tn, fp, fn); results.push_back(accuracy);
        MCC mcc;            double matthewsCorrCoef = mcc.getValue(tp, tn, fp, fn); results.push_back(matthewsCorrCoef);
        F1Score f1;         double f1Score = f1.getValue(tp, tn, fp, fn); results.push_back(f1Score);
        
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
            vector<int> thisBin;
            for (int j = 0; j < bins[i].size(); j++) {  if (bins[i][j] >= 0) { thisBin.push_back(bins[i][j]); } }
                
            if (thisBin.size() != 0) {

                string otu = matrix->getName(thisBin[0]);
                
                for (int j = 1; j < thisBin.size(); j++) { otu += "," + matrix->getName(thisBin[j]); }
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
ListVector* OptiCluster::getList(set<string>& unfittedSeqs) {
    try {
        ListVector* list = new ListVector();
        vector<string> newLabels;
        
        if (removeTrainers) {
            map<int, string> newBins;
            for (int i = 0; i < randomizeSeqs.size(); i++) { //build otus
                
                if (m->getControl_pressed()) { break; }
                
                map<int, int>::iterator it = seqBin.find(randomizeSeqs[i]);
                
                int seqNumber = it->first;
                int binNumber = it->second;
                
                map<int, string>::iterator itBinLabels = binLabels.find(binNumber); //do we have a label for this bin.  If the seq maps to existing bin then we should, otherwise we couldn't "fit" this sequence
                
                if (itBinLabels != binLabels.end()) {
                    map<int, string>::iterator itBin = newBins.find(binNumber); // have we seen this otu yet?
                    
                    if (itBin == newBins.end()) { //create bin
                        newBins[binNumber] = matrix->getName(seqNumber);
                    }else { //append bin
                        newBins[binNumber] += "," + matrix->getName(seqNumber);
                    }
                }else { unfittedSeqs.insert(matrix->getName(seqNumber)); }
            }
            
            for (map<int, string>::iterator itBin = newBins.begin(); itBin != newBins.end(); itBin++) { list->push_back(itBin->second); newLabels.push_back(binLabels[itBin->first]);   }
            
            list->setLabels(newLabels);
            
            ListVector* singleton = matrix->getListSingle();
            
            if (singleton != NULL) { //add in any sequences above cutoff in read. Removing these saves clustering time.
                for (int i = 0; i < singleton->getNumBins(); i++) {
                    string bin = singleton->get(i);
                    if (bin != "") {
                        vector<string> names; util.splitAtComma(bin, names);
                        if (immovableNames.count(names[0]) == 0) { //you are not a reference sequence
                            list->push_back(singleton->get(i));
                        }
                    }
                }
                delete singleton;
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
        
        long long singletn = 0;
        
        if (fitCalc) {
            singletn = matrix->getNumFitSingletons();
            
            for (int i = 0; i < bins.size(); i++) { // we only want to count the bins with movable names
                bool containsMovable = false;
                for (int j = 0; j < bins[i].size(); j++) {
                    if (bins[i][j] >= 0) {
                        if (immovableNames.count(matrix->getName(bins[i][j])) == 0) { containsMovable = true; break; }
                    }
                }
                if (containsMovable) { singletn++; }
            }
            
        }else {
            singletn = matrix->getNumSingletons();
        
            for (int i = 0; i < bins.size(); i++) { if (bins[i].size() != 0) { singletn++; } }
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
        
        //initially there are bins for each sequence (excluding singletons removed on read)
        for (int i = 0; i < bins.size(); i++) {
            
            if (m->getControl_pressed()) { break; }
            
            if (bins[i].size() == 0) { return i;  } //this bin is empty
        }
        
        return -1;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiCluster", "findInsert");
        exit(1);
    }
}

/***********************************************************************/
