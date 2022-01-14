//
//  optifitcluster.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/10/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "optifitcluster.hpp"


/***********************************************************************/
OptiFitCluster::OptiFitCluster(OptiData* mt, ClusterMetric* met, long long ns) : Cluster(), matrix(mt), metric(met), numComboSingletons(ns) {
    m = MothurOut::getInstance();
    maxRefBinNumber = 0;
    closed = false;

    numFitSeqs = 0;  fittruePositives = 0; fitfalsePositives = 0; fitfalseNegatives = 0; fittrueNegatives = 0; numFitSingletons = 0;
    numComboSeqs = 0; numComboSingletons = 0; combotruePositives = 0; combofalsePositives = 0; combofalseNegatives = 0; combotrueNegatives = 0;
}
/***********************************************************************/
int OptiFitCluster::initialize(double& value, bool randomize, vector<vector< string > > existingBins, vector<string> bls, string meth, bool denov) {
    try {
        double reftruePositives, reftrueNegatives, reffalsePositives, reffalseNegatives, numRefSeqs;
        numRefSeqs = 0; reftruePositives = 0; reffalsePositives = 0; reffalseNegatives = 0; reftrueNegatives = 0;
        
        if (meth == "closed") { closed = true; }
        denovo = denov;
        
        vector< vector< long long> > translatedBins;
        randomizeSeqs = matrix->getTranslatedBins(existingBins, translatedBins); //otus in existingBins, otus with matrix names
        
        int binNumber = 0;
        int placeHolderIndex = -1;
        for (long long i = 0; i < translatedBins.size(); i++) {
            binLabels[binNumber] = bls[i];
            bins.push_back(translatedBins[i]);
            numRefSeqs += translatedBins[i].size();
            
            for (int j = 0; j < translatedBins[i].size(); j++) {
                for (int k = 0; k < j; k++) {
                    if (translatedBins[i][j] < 0) { //no dists in matrix
                        translatedBins[i][j] = placeHolderIndex; placeHolderIndex--;
                        reffalsePositives++;
                    }else { //j has distances in the matrix, but is it close to k?
                        if (matrix->isClose(translatedBins[i][j], translatedBins[i][k])) {
                            reftruePositives++;
                        }else { reffalsePositives++; }
                    }
                }
                seqBin[translatedBins[i][j]] = binNumber;
            }
            binNumber++;
        }
        
        maxRefBinNumber = binNumber;
        reffalseNegatives = matrix->getNumRefDists() - reftruePositives; //number of distance in matrix for reference seqs - reftruePositives
        reftrueNegatives = numRefSeqs * (numRefSeqs-1)/2 - (reffalsePositives + reffalseNegatives + reftruePositives);
        
        //add fit seqs as singletons
        int numRefBins = translatedBins.size();
        numFitSingletons = 0;
        //put every fit seq in own bin
        for (long long i = 0; i < randomizeSeqs.size(); i++) {
            vector<long long> thisBin;
            thisBin.push_back(randomizeSeqs[i]);
            bins.push_back(thisBin);
            seqBin[randomizeSeqs[i]] = numRefBins+i;
            
            long long numCloseSeqs = (matrix->getNumFitClose(randomizeSeqs[i])); //does not include self
            fitfalseNegatives += numCloseSeqs;
            if (numCloseSeqs == 0) { numFitSingletons++; } //you are a singletons counted by the matrix as a fitSingleton, but you are not removed because you have ref dists we want to use in the fitting. Don't want to count you twice in stats output.
        }
        numFitSeqs = randomizeSeqs.size();
       
        fitfalseNegatives /= 2; //square matrix
        fittrueNegatives = numFitSeqs * (numFitSeqs-1)/2 - (fitfalsePositives + fitfalseNegatives + fittruePositives); //since everyone is a singleton no one clusters together. True negative = num far apart
        
        numComboSeqs = numRefSeqs + randomizeSeqs.size();
        
        combofalseNegatives = matrix->getNumDists() - reftruePositives; //number of distance in matrix for reference seqs - reftruePositives
        combotrueNegatives = numComboSeqs * (numComboSeqs-1)/2 - (reffalsePositives + reffalseNegatives + reftruePositives);
        combotruePositives = reftruePositives;
        combofalsePositives = reffalsePositives;
        
        double comboValue = metric->getValue(combotruePositives, combotrueNegatives, combofalsePositives, combofalseNegatives);
        
        //add insert location
        seqBin[bins.size()] = -1;
        insertLocation = bins.size();
        vector<long long> temp;
        bins.push_back(temp);
        
        if (randomize) { util.mothurRandomShuffle(randomizeSeqs); }
        
        value = comboValue;
        
        return value;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiFitCluster", "initialize");
        exit(1);
    }
}
/***********************************************************************/
/* for each sequence with mutual information (close)
 * remove from current OTU and calculate MCC when sequence forms its own OTU or joins one of the other OTUs where there is a sequence within the `threshold` (no need to calculate MCC if the paired sequence is already in same OTU and no need to try every OTU - just those where there's a close sequence)
 * keep or move the sequence to the OTU where the `metric` is the largest - flip a coin on ties */
bool OptiFitCluster::update(double& listMetric) {
    try {
        //for each sequence (singletons removed on read)
        for (int i = 0; i < randomizeSeqs.size(); i++) {
            
            if (m->getControl_pressed()) { break; }
            
            map<long long, long long>::iterator it = seqBin.find(randomizeSeqs[i]);
            
            int seqNumber = it->first;
            int binNumber = it->second;
            
            if (binNumber == -1) { }
            else {
                vector<long long> bestBin; bestBin.resize(2, binNumber);
                vector<double> tn; tn.push_back(fittrueNegatives); tn.push_back(combotrueNegatives);
                vector<double> tp; tp.push_back(fittruePositives); tp.push_back(combotruePositives);
                vector<double> fp; fp.push_back(fitfalsePositives); fp.push_back(combofalsePositives);
                vector<double> fn; fn.push_back(fitfalseNegatives); fn.push_back(combofalseNegatives);
                vector<double> bestMetric; bestMetric.resize(2, -1);  //bestMetric[0] = fitSeqs alone, bestMetric[1] = combo or ref and fit
                
                vector<double> bestTp; bestTp.resize(2, 0);
                vector<double> bestTn; bestTn.resize(2, 0);
                vector<double> bestFp; bestFp.resize(2, 0);
                vector<double> bestFn; bestFn.resize(2, 0);
            
                //close / far count in current bin
                vector<double> results = getCloseFarCounts(seqNumber, binNumber);
                double combocCount = results[0];  double combofCount = results[1];
                
                //close / far count in current bin for fit seqs
                vector<double> fitresults = getCloseFarFitCounts(seqNumber, binNumber);
                double fitcCount = fitresults[0];  double fitfCount = fitresults[1];
                
                //fit metrics in current bin
                bestMetric[0] = metric->getValue(tp[0], tn[0], fp[0], fn[0]);
                bestTp[0] = tp[0]; bestTn[0] = tn[0]; bestFp[0] = fp[0]; bestFn[0] = fn[0];
                
                //combo metric in current bin
                bestMetric[1] = metric->getValue(tp[1], tn[1], fp[1], fn[1]);
                bestTp[1] = tp[1]; bestTn[1] = tn[1]; bestFp[1] = fp[1]; bestFn[1] = fn[1];
                
                //if not already singleton, then calc value if singleton was created
                if (!((bins[binNumber].size()) == 1)) {
                    //make a singleton
                    fn[0]+=fitcCount; tn[0]+=fitfCount; fp[0]-=fitfCount; tp[0]-=fitcCount;
                    fn[1]+=combocCount; tn[1]+=combofCount; fp[1]-=combofCount; tp[1]-=combocCount;
                    
                    double singleFitMetric = metric->getValue(tp[0], tn[0], fp[0], fn[0]);
                    double singleComboMetric = metric->getValue(tp[1], tn[1], fp[1], fn[1]);
                    if ((singleFitMetric > bestMetric[0]) || (singleComboMetric > bestMetric[1])) {
                        bestBin[1] = -1; bestTp[1] = tp[1]; bestTn[1] = tn[1]; bestFp[1] = fp[1]; bestFn[1] = fn[1];
                        bestMetric[1] = singleComboMetric;
                        
                        bestBin[0] = -1; bestTp[0] = tp[0]; bestTn[0] = tn[0]; bestFp[0] = fp[0]; bestFn[0] = fn[0];
                        bestMetric[0] = singleFitMetric;
                    }
                }
                
                set<long long> binsToTry;
                set<long long> closeSeqs = matrix->getCloseRefSeqs(seqNumber);
                for (set<long long>::iterator itClose = closeSeqs.begin(); itClose != closeSeqs.end(); itClose++) { binsToTry.insert(seqBin[*itClose]); }
                
                //merge into each "close" otu
                vector<vector<double> > ties; vector<vector<double> > ties0;
                for (set<long long>::iterator it = binsToTry.begin(); it != binsToTry.end(); it++) {
                    //reset tn, tp,fp,fn values to original bin
                    tn[0] = fittrueNegatives; tp[0] = fittruePositives; fp[0] = fitfalsePositives; fn[0] = fitfalseNegatives;
                    tn[1] = combotrueNegatives; tp[1] = combotruePositives; fp[1] = combofalsePositives; fn[1] = combofalseNegatives;
                    
                    //move out of old bin
                    fn[0]+=fitcCount; tn[0]+=fitfCount; fp[0]-=fitfCount; tp[0]-=fitcCount;
                    fn[1]+=combocCount; tn[1]+=combofCount; fp[1]-=combofCount; tp[1]-=combocCount;
                    
                    results = getCloseFarCounts(seqNumber, *it); //results[0] = close count, results[1] = far count
                    fn[1]-=results[0]; tn[1]-=results[1];  tp[1]+=results[0]; fp[1]+=results[1]; //move into new bin
                    
                    results = getCloseFarFitCounts(seqNumber, *it);
                    fn[0]-=results[0]; tn[0]-=results[1];  tp[0]+=results[0]; fp[0]+=results[1]; //move into new bin - only consider fit seqs
                    
                    double newComboMetric = metric->getValue(tp[1], tn[1], fp[1], fn[1]); //score when sequence is moved
                    double newFitMetric = metric->getValue(tp[0], tn[0], fp[0], fn[0]); //score when sequence is moved
                    //new best
                    if (newComboMetric > bestMetric[1]) {
                        ties.clear(); ties0.clear();
                        bestMetric[1] = newComboMetric; bestBin[1] = (*it); bestTp[1] = tp[1]; bestTn[1] = tn[1]; bestFp[1] = fp[1]; bestFn[1] = fn[1];
                        bestMetric[0] = newFitMetric; bestBin[0] = (*it); bestTp[0] = tp[0]; bestTn[0] = tn[0]; bestFp[0] = fp[0]; bestFn[0] = fn[0];
                        vector<double> tie; tie.push_back(bestMetric[1]); tie.push_back(bestBin[1]); tie.push_back(bestTp[1]);
                        tie.push_back(bestTn[1]); tie.push_back(bestFp[1]); tie.push_back(bestFn[1]); ties.push_back(tie);
                        vector<double> tie0; tie0.push_back(bestMetric[0]); tie0.push_back(bestBin[0]); tie0.push_back(bestTp[0]);
                        tie0.push_back(bestTn[0]); tie0.push_back(bestFp[0]); tie0.push_back(bestFn[0]); ties0.push_back(tie0);

                        
                    }else if (newComboMetric == bestMetric[1]) {
                        bestMetric[1] = newComboMetric; bestBin[1] = (*it); bestTp[1] = tp[1]; bestTn[1] = tn[1]; bestFp[1] = fp[1]; bestFn[1] = fn[1];
                        bestMetric[0] = newFitMetric; bestBin[0] = (*it); bestTp[0] = tp[0]; bestTn[0] = tn[0]; bestFp[0] = fp[0]; bestFn[0] = fn[0];
                        vector<double> tie; tie.push_back(bestMetric[1]); tie.push_back(bestBin[1]); tie.push_back(bestTp[1]);
                        tie.push_back(bestTn[1]); tie.push_back(bestFp[1]); tie.push_back(bestFn[1]); ties.push_back(tie);
                        vector<double> tie0; tie0.push_back(bestMetric[0]); tie0.push_back(bestBin[0]); tie0.push_back(bestTp[0]);
                        tie0.push_back(bestTn[0]); tie0.push_back(bestFp[0]); tie0.push_back(bestFn[0]); ties0.push_back(tie0);

                    }
                }
                
                if (ties.size() > 1) {
                    int randomTie = util.getRandomIndex((int)ties.size()-1);
                    bestMetric[1] = ties[randomTie][0]; bestBin[1] = ties[randomTie][1]; bestTp[1] = ties[randomTie][2]; bestTn[1] = ties[randomTie][3]; bestFp[1] = ties[randomTie][4]; bestFn[1] = ties[randomTie][5];
                    bestMetric[0] = ties0[randomTie][0]; bestBin[0] = ties0[randomTie][1]; bestTp[0] = ties0[randomTie][2]; bestTn[0] = ties0[randomTie][3]; bestFp[0] = ties0[randomTie][4]; bestFn[0] = ties0[randomTie][5];
                }
                
                //how to choose the best bin if they differ????
                long long newBin = bestBin[1];
                
                bool usedInsert = false;
                if (newBin == -1) {  newBin = insertLocation;  usedInsert = true;  }
                
                if (newBin != binNumber) {
                    combotruePositives = bestTp[1]; combotrueNegatives = bestTn[1]; combofalsePositives = bestFp[1]; combofalseNegatives = bestFn[1];
                    fittruePositives = bestTp[0]; fittrueNegatives = bestTn[0]; fitfalsePositives = bestFp[0]; fitfalseNegatives = bestFn[0];
                    
                    //move seq from i to j
                    bins[newBin].push_back(seqNumber); //add seq to bestbin
                    bins[binNumber].erase(remove(bins[binNumber].begin(), bins[binNumber].end(), seqNumber), bins[binNumber].end()); //remove from old bin i
                }
                
                if (usedInsert) { insertLocation = findInsert(); }
                
                //update seqBins
                seqBin[seqNumber] = newBin; //set new OTU location
            }
        }
        
        listMetric = metric->getValue(combotruePositives, combotrueNegatives, combofalsePositives, combofalseNegatives); 
        
       
        if (m->getDebug()) { ListVector* list = getList(); list->print(cout); delete list; }
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "OptiFitCluster", "update");
        exit(1);
    }
}
/***********************************************************************/
vector<double> OptiFitCluster::getCloseFarCounts(long long seq, long long newBin) {
    try {
        vector<double> results; results.push_back(0); results.push_back(0); //results[0] = close count, results[1] = far count
        
        if (newBin == -1) { }  //making a singleton bin. Close but we are forcing apart.
        else { //merging a bin
            for (long long i = 0; i < bins[newBin].size(); i++) {
                if (seq == bins[newBin][i]) {} //ignore self
                else if (!matrix->isClose(seq, bins[newBin][i])) { results[1]++; }  //this sequence is "far away" from sequence i - above the cutoff
                else { results[0]++;  }  //this sequence is "close" to sequence i - distance between them is less than cutoff
            }
        }
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiFitCluster", "getCloseFarCounts");
        exit(1);
    }
}
/***********************************************************************/
vector<double> OptiFitCluster::getCloseFarFitCounts(long long seq, long long newBin) {
    try {
        vector<double> results; results.push_back(0); results.push_back(0); //results[0] = close count, results[1] = far count
        
        if (newBin == -1) { }  //making a singleton bin. Close but we are forcing apart.
        else { //merging a bin
            for (long long i = 0; i < bins[newBin].size(); i++) {
                
                if (seq == bins[newBin][i]) {} //ignore self
                else {
                    bool isFit = true;
                    bool closeFit = matrix->isCloseFit(seq, bins[newBin][i], isFit);
                    if (closeFit) { //you are close if you are fit and close
                        results[0]++;
                    }else if (isFit) { results[1]++; } //this sequence is "far away" and fit - above the cutoff
                }
            }
        }
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiFitCluster", "getCloseFarCounts");
        exit(1);
    }
}

/***********************************************************************/
vector<double> OptiFitCluster::getStats(double& tp,  double& tn,  double& fp,  double& fn) {
    try {
        double singletn = 0;
        if (!closed) { singletn = matrix->getNumSingletons(); }
        double tempnumSeqs = numComboSeqs + singletn;
        
        tp = combotruePositives;
        fp = combofalsePositives;
        fn = combofalseNegatives;
        tn = tempnumSeqs * (tempnumSeqs-1)/2 - (combofalsePositives + combofalseNegatives + combotruePositives); //adds singletons to tn
        
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
vector<double> OptiFitCluster::getFitStats(double& tp,  double& tn,  double& fp,  double& fn) {
    try {
        double singletn = 0;
        if (!closed) { singletn = matrix->getNumFitTrueSingletons(); }
        double tempnumSeqs = numFitSeqs + singletn; //numFitSingletons are reads that are selected as the fit seqs, that have dists to reference but not dists to other fit seqs. They are included
        
        tp = fittruePositives;
        fp = fitfalsePositives;
        fn = fitfalseNegatives;
        tn = tempnumSeqs * (tempnumSeqs-1)/2 - (fitfalsePositives + fitfalseNegatives + fittruePositives); //adds singletons to tn
        
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
        m->errorOut(e, "OptiFitCluster", "getFitStats");
        exit(1);
    }
}
/***********************************************************************/
ListVector* OptiFitCluster::getList() {
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
        m->errorOut(e, "OptiFitCluster", "getList");
        exit(1);
    }
}
/***********************************************************************/
ListVector* OptiFitCluster::getFittedList(string label, bool includerefs) {
    try {
        ListVector* list = new ListVector();
        
        map<long long, string> newBins;
        set<long long> unFitted;
        long long numListSeqs = 0;
        for (long long i = 0; i < randomizeSeqs.size(); i++) { //build otus
            
            if (m->getControl_pressed()) { break; }
            
            map<long long, long long>::iterator it = seqBin.find(randomizeSeqs[i]);
            
            long long seqNumber = it->first;
            long long binNumber = it->second;
            
            map<long long, string>::iterator itBinLabels = binLabels.find(binNumber); //do we have a label for this bin.  If the seq maps to existing bin then we should, otherwise we couldn't "fit" this sequence
            
            if (itBinLabels != binLabels.end()) {
                numListSeqs++;
                map<long long, string>::iterator itBin = newBins.find(binNumber); // have we seen this otu yet?
                
                if (itBin == newBins.end()) { //create bin
                    newBins[binNumber] = matrix->getName(seqNumber);
                }else { //append bin
                    newBins[binNumber] += "," + matrix->getName(seqNumber);
                }
            }else { unFitted.insert(seqNumber); }
        }
        
        if (denovo || includerefs) { //add in refs
            vector<long long> refs = matrix->getRefSeqs();
            
            for (long long i = 0; i < refs.size(); i++) {
                if (m->getControl_pressed()) { break; }
                map<long long, long long>::iterator it = seqBin.find(refs[i]);
                
                long long seqNumber = it->first;
                long long binNumber = it->second;
                
                map<long long, string>::iterator itBin = newBins.find(binNumber); // have we seen this otu yet?
                
                if (itBin == newBins.end()) { //create bin
                    newBins[binNumber] = matrix->getName(seqNumber);
                }else { //append bin
                    newBins[binNumber] += "," + matrix->getName(seqNumber);
                }
            }
        }
        
        //numFitSeqs does not include any kind of singleton
        long long numUnFitted = (numFitSeqs + numFitSingletons - numListSeqs); //getNumFitTrueSingletons are fit reads that have no dists in the matrix. This can be confusing, think of it like this: there are true singletons, meaning we don't care if you are a ref or fit and you have no dists below the cutoff. This means you will be in your own OTU no matter what we do. There are fitSingletons, meaning you are a fit sequence and have no dists below the cutoff that coorespond to other fit seqs ( NOTE: you may or may not have dists to ref seqs or you could be a true singleton or a just a singleton because of the references chosen).
        
        long long numSingletonBins = 0;
        if ((label != "") && (numUnFitted != 0)) {
            
            m->mothurOut("\nFitted " + toString(numListSeqs) + " sequences to " + toString(newBins.size()) + " existing OTUs.\n");
            
            if (!closed) { //cluster the unfitted seqs separately
                m->mothurOut(toString(numUnFitted) + " sequences were unable to be fitted existing OTUs, excluding singletons.\n");
                
                m->mothurOut("\n**************** Clustering the unfitted sequences ****************\n");
                
                OptiData* unFittedMatrix = matrix->extractMatrixSubset(unFitted);
                
                ListVector* unfittedList = clusterUnfitted(unFittedMatrix, label); //unfittedList includes unfitted singletons
                
                if (unfittedList != NULL) {
                    
                    m->mothurOut("The unfitted sequences clustered into " + toString(unfittedList->getNumBins()) + " new OTUs.\n"); //+unFittedMatrix->getNumSingletons()+ matrix->getNumFitSingletons()
                    
                    for (int i = 0; i < unfittedList->getNumBins(); i++) {
                        string bin = unfittedList->get(i);
                        if (bin != "") { list->push_back(unfittedList->get(i)); }
                    }
                    delete unfittedList;
                }
                delete unFittedMatrix;
                
                m->mothurOut("\n*******************************************************************\n\n");
                
                //add in fit singletons 
                ListVector* singleton = matrix->getFitListSingle();
                
                if (singleton != NULL) { //add in any sequences above cutoff in read. Removing these saves clustering time.
                    for (int i = 0; i < singleton->getNumBins(); i++) {
                        if (m->getControl_pressed()) { break; }
                        if (singleton->get(i) != "") { list->push_back(singleton->get(i)); }
                    }
                    numSingletonBins += singleton->getNumBins();
                    delete singleton;
                }
                
            }else {
                m->mothurOut("\nSequences that were unable to be fitted existing OTUs will be listed in the *.optifit_scrap.accnos file.\n");
                unfittedNames = matrix->getNames(unFitted);
                
                //add in fit singletons
                ListVector* singleton = matrix->getFitListSingle();
                
                if (singleton != NULL) { //add in any sequences above cutoff in read. Removing these saves clustering time.
                    for (int i = 0; i < singleton->getNumBins(); i++) {
                        if (m->getControl_pressed()) { break; }
                        if (singleton->get(i) != "") { unfittedNames.insert(singleton->get(i)); }
                    }
                    delete singleton;
                }
            }
        }else {
            if (label != "") {
                m->mothurOut("\nFitted all " + toString(list->getNumSeqs()) + " sequences to existing OTUs. \n");
            }
        }
        
        vector<string> newLabels = list->getLabels();
        
        for (map<long long, string>::iterator itBin = newBins.begin(); itBin != newBins.end(); itBin++) {
            list->push_back(itBin->second);
            newLabels.push_back(binLabels[itBin->first]);
        }
        
        list->setLabels(newLabels);
    
        return list;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiFitCluster", "getFittedList");
        exit(1);
    }
}
/***********************************************************************/
ListVector* OptiFitCluster::clusterUnfitted(OptiData* unfittedMatrix, string label) {
    try {
        ListVector* list = NULL;
        
        OptiCluster cluster(unfittedMatrix, metric, 0);
        
        int iters = 0;
        double listVectorMetric = 0; //worst state
        double delta = 1;
        
        cluster.initialize(listVectorMetric, true, "singleton");
        
        long long numBins = cluster.getNumBins();
        m->mothurOut("\n\niter\ttime\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");
        
        double tp, tn, fp, fn;
        vector<double> results = cluster.getStats(tp, tn, fp, fn);
        m->mothurOut("0\t0\t" + label + "\t" + toString(numBins) + "\t"+ label + "\t" + toString(tp) + "\t" + toString(tn) + "\t" + toString(fp) + "\t" + toString(fn) + "\t");
       
        for (int i = 0; i < results.size(); i++) { m->mothurOut(toString(results[i]) + "\t");  }
        m->mothurOutEndLine();
        
        while ((delta > 0.0001) && (iters < 100)) {
            
            long start = time(NULL);
            
            if (m->getControl_pressed()) { break; }
            double oldMetric = listVectorMetric;
            
            cluster.update(listVectorMetric);
            
            delta = abs(oldMetric - listVectorMetric);
            iters++;
            
            results = cluster.getStats(tp, tn, fp, fn);
            numBins = cluster.getNumBins();
            
            m->mothurOut(toString(iters) + "\t" + toString(time(NULL) - start) + "\t" + label + "\t" + toString(numBins) + "\t" + label + "\t"+ toString(tp) + "\t" + toString(tn) + "\t" + toString(fp) + "\t" + toString(fn) + "\t");
            
            for (int i = 0; i < results.size(); i++) { m->mothurOut(toString(results[i]) + "\t");  }
            m->mothurOutEndLine();
            
        }
        m->mothurOutEndLine(); m->mothurOutEndLine();
        
        if (m->getControl_pressed()) { return list; }
        
        list = cluster.getList();
        list->setLabel(label);
        
        return list;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiFitCluster", "clusterUnfitted");
        exit(1);
    }
}

/***********************************************************************/
long long OptiFitCluster::getNumBins() {
    try {
        
        long long singletn = 0;
        
        singletn = matrix->getNumSingletons();
            
        for (int i = 0; i < bins.size(); i++) { if (bins[i].size() != 0) { singletn++; } }
        
        return singletn;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiFitCluster", "getNumBins");
        exit(1);
    }
}
/***********************************************************************/
long long OptiFitCluster::getNumFitBins() {
    try {
        ListVector* list = getFittedList("", false);
        
        int numBins = 0;
        if (list != NULL) {
            numBins = list->getNumBins();
            delete list;
        }
        
        return numBins;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiFitCluster", "getNumFitBins");
        exit(1);
    }
}
/***********************************************************************/
int OptiFitCluster::findInsert() {
    try {
        
        //initially there are bins for each sequence (excluding singletons removed on read)
        for (int i = 0; i < bins.size(); i++) {
            
            if (m->getControl_pressed()) { break; }
            
            if (bins[i].size() == 0) { return i;  } //this bin is empty
        }
        
        return -1;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiFitCluster", "findInsert");
        exit(1);
    }
}

/***********************************************************************/
