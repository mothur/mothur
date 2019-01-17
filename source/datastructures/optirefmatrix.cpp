//
//  optirefmatrix.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/3/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "optirefmatrix.hpp"

#include "counttable.h"

/***********************************************************************/
OptiRefMatrix::OptiRefMatrix(string distFile, string distFormat, string dupsFile, string dupsFormat, double c, float fP, string refWeight) : OptiData(c) {
    
    numFitSingletons = 0;
    numRefSingletons = 0;
    numSingletons = 0;
    numBetweenDists = 0;
    numFitDists = 0;
    numRefDists = 0;
    numFitSeqs = 0;
    refWeightMethod = refWeight;
    
    fitPercent = fP / 100.0;
    if (fitPercent < 0.001) { fitPercent = 0.10; m->mothurOut("[WARNING]: fit percentage must be between 0.001 (0.1%) and 1.0 (100%). Setting to 0.10 or 10%. \n"); } //minumum of 0.1%
    else if (fitPercent > 100.0) {  m->mothurOut("[ERROR]: fit percentage must be between 0.0001 and 100.0\n"); m->setControl_pressed(true); }
    
    square = false;
    
    set<string> noRefNamesSet;
    readFiles(distFile, distFormat, dupsFile, dupsFormat, noRefNamesSet);
}
/***********************************************************************/
OptiRefMatrix::OptiRefMatrix(string distFile, string distFormat, string dupsFile, string dupsFormat, double c, string accnosfile) : OptiData(c) {
    
    numFitSingletons = 0;
    numRefSingletons = 0;
    numSingletons = 0;
    numBetweenDists = 0;
    numFitDists = 0;
    numRefDists = 0;
    numFitSeqs = 0;
    refWeightMethod = "accnos";
    
    set<string> accnosRefFileNames = util.readAccnos(accnosfile);
    
    square = false;
    
    readFiles(distFile, distFormat, dupsFile, dupsFormat, accnosRefFileNames);
}

/***********************************************************************/
OptiRefMatrix::OptiRefMatrix(string d, string nc, string f, string df, double c, string fit, string fitnc, string fitf, string fitdf, string betweend, string betweendf) : OptiData(c) {
    
    string refdistfile, refnamefile, refcountfile, refformat, refdistformat, fitdistfile, fitnamefile, fitcountfile, fitformat, fitdistformat, betweendistfile, betweendistformat;
    
    refdistfile = d; refdistformat = df; refformat = f; fitdistfile = fit; fitdistformat = fitdf; fitformat = fitf; betweendistfile = betweend; betweendistformat = betweendf;
    
    numFitSingletons = 0;
    numRefSingletons = 0;
    numSingletons = 0;
    numBetweenDists = 0;
    numFitDists = 0;
    numRefDists = 0;
    numFitSeqs = 0;
    
    fitPercent = 0;
    refWeightMethod = "none";
    
    square = false;
    
    if (refformat == "name") { refnamefile = nc; refcountfile = ""; }
    else if (refformat == "count") { refcountfile = nc; refnamefile = ""; }
    else { refcountfile = ""; refnamefile = ""; }
    
    if (fitformat == "name") { fitnamefile = fitnc; fitcountfile = ""; }
    else if (fitformat == "count") { fitcountfile = fitnc; fitnamefile = ""; }
    else { fitcountfile = ""; fitnamefile = ""; }
    
    readFiles(refdistfile, refnamefile, refcountfile, refformat, refdistformat, fitdistfile, fitnamefile, fitcountfile, fitformat, fitdistformat, betweendistfile, betweendistformat);
}
/***********************************************************************/
//Since we are extracting a subset of the seqs some reads that may not have been singletons
OptiData* OptiRefMatrix::extractRefMatrix() {
    try {
        set<long long> seqs; for (long long i = 0; i < isRef.size(); i++) { if (isRef[i]) { seqs.insert(i); } }
            
        vector<string> subsetNameMap;
        vector<string> subsetSingletons;
        vector< set<long long> > subsetCloseness;
        map<long long, long long> thisNameMap;
        map<long long, long long> nonSingletonNameMap;
        vector<bool> singleton; singleton.resize(seqs.size(), true);
        int count = 0;
        
        for (set<long long>::iterator it = seqs.begin(); it != seqs.end(); it++) {
            long long seqNum = *it;
            thisNameMap[seqNum] = count;
            nonSingletonNameMap[count] = seqNum;
            
            set<long long> thisSeqsCloseSeqs = getCloseSeqs(seqNum);
            for (set<long long>::iterator itClose = thisSeqsCloseSeqs.begin(); itClose != thisSeqsCloseSeqs.end(); itClose++) {
                
                if (m->getControl_pressed()) { break; }
                
                long long thisSeq = *itClose;
                
                //is this seq in the set of unfitted?
                if (seqs.count(thisSeq) != 0) { singleton[thisNameMap[seqNum]] = false; }
            }
            count++;
        }
        
        int nonSingletonCount = 0;
        for (long long i = 0; i < singleton.size(); i++) {
            if (!singleton[i]) { //if you are not a singleton
                nonSingletonNameMap[i] = nonSingletonCount;
                nonSingletonCount++;
            }else { seqs.erase(nonSingletonNameMap[i]);  subsetSingletons.push_back(getName(nonSingletonNameMap[i])); } //remove from unfitted
        }
        singleton.clear();
        
        subsetCloseness.resize(nonSingletonCount);
        for (set<long long>::iterator it = seqs.begin(); it != seqs.end(); it++) {
            
            if (m->getControl_pressed()) { break; }
            
            long long seqNum = *it;
            
            set<long long> thisSeqsCloseSeqs = getCloseSeqs(seqNum);
            set<long long> thisSeqsCloseUnFittedSeqs;
            for (set<long long>::iterator itClose = thisSeqsCloseSeqs.begin(); itClose != thisSeqsCloseSeqs.end(); itClose++) {
                
                if (m->getControl_pressed()) { break; }
                
                long long thisSeq = *itClose;
                
                //is this seq in the set of unfitted?
                if (seqs.count(thisSeq) != 0) { thisSeqsCloseUnFittedSeqs.insert(nonSingletonNameMap[thisNameMap[thisSeq]]); }
            }
            
            if (!thisSeqsCloseUnFittedSeqs.empty()) {
                subsetCloseness[nonSingletonNameMap[thisNameMap[seqNum]]] = thisSeqsCloseUnFittedSeqs;
                subsetNameMap.push_back(getName(seqNum));
            }
            
        }
        
        for (int i = 0; i < isSingleRef.size(); i++) { if (isSingleRef[i]) { subsetSingletons.push_back(singletons[i]); } }
            
        OptiData* unfittedMatrix = new OptiMatrix(subsetCloseness, subsetNameMap, subsetSingletons, cutoff);
        
        return unfittedMatrix;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "extractRefMatrix");
        exit(1);
    }
}
/***********************************************************************/
//given matrix indexes of seqs, pull out their dists and create optimatrix
OptiData* OptiRefMatrix::extractMatrixSubset(set<long long> & seqs) {
    try {
        vector<string> subsetNameMap;
        vector<string> subsetSingletons;
        vector< set<long long> > subsetCloseness;
        map<long long, long long> thisNameMap;
        map<long long, long long> nonSingletonNameMap;
        vector<bool> singleton; singleton.resize(seqs.size(), true);
        int count = 0;
        
        for (set<long long>::iterator it = seqs.begin(); it != seqs.end(); it++) {
            long long seqNum = *it;
            thisNameMap[seqNum] = count;
            nonSingletonNameMap[count] = seqNum;
            
            set<long long> thisSeqsCloseSeqs = getCloseSeqs(seqNum);
            for (set<long long>::iterator itClose = thisSeqsCloseSeqs.begin(); itClose != thisSeqsCloseSeqs.end(); itClose++) {
                
                if (m->getControl_pressed()) { break; }
                
                long long thisSeq = *itClose;
                
                //is this seq in the set of unfitted?
                if (seqs.count(thisSeq) != 0) { singleton[thisNameMap[seqNum]] = false; }
            }
            count++;
        }
        
        int nonSingletonCount = 0;
        for (long long i = 0; i < singleton.size(); i++) {
            if (!singleton[i]) { //if you are a singleton
                nonSingletonNameMap[i] = nonSingletonCount;
                nonSingletonCount++;
            }else { seqs.erase(nonSingletonNameMap[i]);  subsetSingletons.push_back(getName(nonSingletonNameMap[i])); } //remove from unfitted
        }
        singleton.clear();
        
        subsetCloseness.resize(nonSingletonCount);
        for (set<long long>::iterator it = seqs.begin(); it != seqs.end(); it++) {
            
            if (m->getControl_pressed()) { break; }
            
            long long seqNum = *it;
            
            set<long long> thisSeqsCloseSeqs = getCloseSeqs(seqNum);
            set<long long> thisSeqsCloseUnFittedSeqs;
            for (set<long long>::iterator itClose = thisSeqsCloseSeqs.begin(); itClose != thisSeqsCloseSeqs.end(); itClose++) {
                
                if (m->getControl_pressed()) { break; }
                
                long long thisSeq = *itClose;
                
                //is this seq in the set of unfitted?
                if (seqs.count(thisSeq) != 0) { thisSeqsCloseUnFittedSeqs.insert(nonSingletonNameMap[thisNameMap[thisSeq]]); }
            }
            
            if (!thisSeqsCloseUnFittedSeqs.empty()) {
                subsetCloseness[nonSingletonNameMap[thisNameMap[seqNum]]] = thisSeqsCloseUnFittedSeqs;
                subsetNameMap.push_back(getName(seqNum));
            }
            
        }
        
        OptiData* unfittedMatrix = new OptiMatrix(subsetCloseness, subsetNameMap, subsetSingletons, cutoff);
        
        return unfittedMatrix;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "extractMatrixSubset");
        exit(1);
    }
}
/***********************************************************************/
vector<long long> OptiRefMatrix::getTranslatedBins(vector<vector<string> > & binNames, vector<vector<long long> > & fixedBins) {
    try {
        fixedBins.clear();
        
        map<string, long long> nameIndexes;
        set<string> unique;
        for (long long i = 0; i < nameMap.size(); i++) { //vector of string representing the sequences in the matrix from the name file.
            vector<string> thisSeqsReps; util.splitAtComma(nameMap[i], thisSeqsReps); //split redundant names
            if (i < closeness.size()) {  nameIndexes[thisSeqsReps[0]] = i;  } //this is a sequence with distances in the matrix
            if (thisSeqsReps.size() == 1) { //you are unique
                unique.insert(thisSeqsReps[0]);
            }
        }
        
        for (long long i = 0; i < singletons.size(); i++) {
            if (isSingleRef[i]) {
                vector<string> thisSeqsReps; util.splitAtComma(singletons[i], thisSeqsReps); //split redundant names
                nameIndexes[thisSeqsReps[0]] = -1;
                if (thisSeqsReps.size() == 1) { unique.insert(thisSeqsReps[0]); }
            }
        }
        
        for (long long i = 0; i < binNames.size(); i++) { //for each OTU
            vector<long long> thisBinsSeqs;
            for (long long j = 0; j < binNames[i].size(); j++) { //for each sequence
                map<string, long long>::iterator it = nameIndexes.find(binNames[i][j]);
                
                if (it == nameIndexes.end()) { }//not in distance matrix, but needs a value in fixedBins. 2 reasons for making it here: you are a redundant name in the listfile, you do not have any distances
                else { thisBinsSeqs.push_back(it->second);  } //"name" of sequence in matrix
            }
            fixedBins.push_back(thisBinsSeqs);
        }
        
        return (getFitSeqs());
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "getTranslatedBins");
        exit(1);
    }
}
/***********************************************************************/
//assumes that i is a fitSeq
bool OptiRefMatrix::isCloseFit(long long i, long long toFind, bool& isFit){
    try {
        if (i < 0) { return false; }
        else if (i > closeness.size()) { m->mothurOut("[ERROR]: index is not valid.\n"); m->setControl_pressed(true); return false; }
        
        bool found = false;
        if (!isRef[toFind]) { //are you a fit seq
            if (closeness[i].count(toFind) != 0) {  //are you close
                found = true;
            }
            isFit = true;
        }else { isFit = false;  }
        return found;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "isCloseFit");
        exit(1);
    }
}
/***********************************************************************/
//does not include singletons, only reads in closeness
vector<long long> OptiRefMatrix::getRefSeqs() {
    try {
        vector<long long> refSeqsIndexes;
        for (long long i = 0; i < isRef.size(); i++) {
            if (isRef[i]) { refSeqsIndexes.push_back(i); }
        }
        return refSeqsIndexes;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "getRefSeqs");
        exit(1);
    }
}
/***********************************************************************/
vector<string> OptiRefMatrix::getRefSingletonNames() {
    try {
        vector<string> refSeqsNames;
    
        for (long long i = 0; i < isSingleRef.size(); i++) {
            if (isSingleRef[i]) { refSeqsNames.push_back(singletons[i]); }
        }
        
        return refSeqsNames;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "getRefSingletonNames");
        exit(1);
    }
}
/***********************************************************************/
vector<long long> OptiRefMatrix::getFitSeqs() {
    try {
        vector<long long> fitSeqsIndexes;
        for (long long i = 0; i < isRef.size(); i++) {
            if (!isRef[i]) { fitSeqsIndexes.push_back(i);  }
        }
        return fitSeqsIndexes;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "getFitSeqs");
        exit(1);
    }
}
/***********************************************************************/
long long OptiRefMatrix::getNumFitTrueSingletons() {
    try {
        long long numFitTrueSingletons = 0;
        
        for (long long i = 0; i < isSingleRef.size(); i++) {
            if (!isSingleRef[i]) { numFitTrueSingletons++; }
        }
        
        return numFitTrueSingletons;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiData", "getNumFitTrueSingletons");
        exit(1);
    }
}
/***********************************************************************/
long long OptiRefMatrix::getNumFitClose(long long index) {
    try {
        long long numClose = 0;
        
        if (index < 0) { }
        else if (index > closeness.size()) { m->mothurOut("[ERROR]: index is not valid.\n"); m->setControl_pressed(true);  }
        else {
            //reference seqs all have indexes less than refEnd
            for (set<long long>::iterator it = closeness[index].begin(); it != closeness[index].end(); it++) {
                if (!isRef[*it]) {  numClose++; } //you are a fit seq
            }
        }
        
        return numClose;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiData", "getNumClose");
        exit(1);
    }
}
/***********************************************************************/
long long OptiRefMatrix::getNumRefClose(long long index) {
    try {
        long long numClose = 0;
        
        if (index < 0) { }
        else if (index > closeness.size()) { m->mothurOut("[ERROR]: index is not valid.\n"); m->setControl_pressed(true);  }
        else {
            //reference seqs all have indexes less than refEnd
            for (set<long long>::iterator it = closeness[index].begin(); it != closeness[index].end(); it++) {
                if (isRef[*it]) {  numClose++; } //you are a ref seq
            }
        }
        
        return numClose;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiData", "getNumClose");
        exit(1);
    }
}
/***********************************************************************/
set<long long> OptiRefMatrix::getCloseFitSeqs(long long index){
    try {
        set<long long> closeSeqs;
        
        if (index < 0) { }
        else if (index > closeness.size()) { m->mothurOut("[ERROR]: index is not valid.\n"); m->setControl_pressed(true);  } //
        else {
            //reference seqs all have indexes less than refEnd
            for (set<long long>::iterator it = closeness[index].begin(); it != closeness[index].end(); it++) {
                if (!isRef[*it]) {  closeSeqs.insert(*it); } //you are a fit seq
            }
        }
        
        return closeSeqs;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "getCloseFitSeqs");
        exit(1);
    }
}
/***********************************************************************/
set<long long> OptiRefMatrix::getCloseRefSeqs(long long index){
    try {
        set<long long> closeSeqs;
        
        if (index < 0) { }
        else if (index > closeness.size()) { m->mothurOut("[ERROR]: index is not valid.\n"); m->setControl_pressed(true);  }
        else {
            //reference seqs all have indexes less than refEnd
            for (set<long long>::iterator it = closeness[index].begin(); it != closeness[index].end(); it++) {
                if (isRef[*it]) { closeSeqs.insert(*it); } //you are a ref seq
            }
        }
        
        return closeSeqs;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "getCloseFitSeqs");
        exit(1);
    }
}
/***********************************************************************/
ListVector* OptiRefMatrix::getFitListSingle() {
    try {
        ListVector* singlelist = NULL;
        
        if (singletons.size() == 0) { }
        else {
            singlelist = new ListVector();
            
            for (int i = 0; i < isSingleRef.size(); i++) { if (!isRef[i]) { singlelist->push_back(singletons[i]); } }
        }
        
        return singlelist;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "getFitListSingle");
        exit(1);
    }
}
/***********************************************************************/
void OptiRefMatrix::randomizeRefs() {
    try {
        long long totalSeqs = (isRef.size()+isSingleRef.size());
        long long numToSelect = totalSeqs * fitPercent;
        long long refSingletonCutoff = isRef.size();
        long long singleSize = isSingleRef.size();
        
        //select sequences to be reference
        set<long long> fitSeqsIndexes;
        if (weights.size() != 0) {  fitSeqsIndexes = subsample.getWeightedSample(weights, numToSelect);  } //you have weighted selection
        else {
            long long numSelected = 0;
            while (numSelected < numToSelect) {
                if (m->getControl_pressed()) { break; }
                fitSeqsIndexes.insert(util.getRandomIndex(totalSeqs-1)); //no repeats
                numSelected = fitSeqsIndexes.size();
            }
        }
        
        //initilize isRef to true
        isRef.clear(); isRef.resize(refSingletonCutoff, true);
        isSingleRef.clear(); isSingleRef.resize(singleSize, true);
        
        //set isRef values
        for (set<long long>::iterator it = fitSeqsIndexes.begin(); it != fitSeqsIndexes.end(); it++) {
            if (m->getControl_pressed()) { break; }
            
            long long thisSeq = *it;
            if (thisSeq < refSingletonCutoff) { //you are a non singleton seq in the closeness
                isRef[thisSeq] = false;
            }else { //thisSeq is a singleton
                isSingleRef[thisSeq-refSingletonCutoff] = false;
            }
        }
        
        //find number of fitDists, refDists and between dists
        calcCounts();
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "randomizeRefs");
        exit(1);
    }

}
/***********************************************************************/
//for denovo method
int OptiRefMatrix::readFiles(string distFile, string distFormat, string dupsFile, string dupsFormat, set<string>& optionalRefNames) {
    try {
        string namefile, countfile;
        if (dupsFormat == "name") { namefile = dupsFile; countfile = ""; }
        else if (dupsFormat == "count") { countfile = dupsFile; namefile = ""; }
        else { countfile = ""; namefile = ""; }
        
        map<string, long long> nameAssignment;
        if (namefile != "") { util.readNames(namefile, nameAssignment); }
        else  {
            CountTable ct; ct.readTable(countfile, false, true);
            map<string, int> temp = ct.getNameMap();
            for (map<string, int>::iterator it = temp.begin(); it!= temp.end(); it++) {  nameAssignment[it->first] = it->second; }
        }
        
        //select sequences to be reference
        set<long long> fitSeqsIndexes;
        long long count = 0;
        for (map<string, long long>::iterator it = nameAssignment.begin(); it!= nameAssignment.end(); it++) {
            if (refWeightMethod == "abundance")          { weights[count] = it->second; }
            else if (refWeightMethod == "connectivity")  { weights[count] = 1;          } //initialize
            else if (refWeightMethod == "accnos") { //fill fit indexes
                if (optionalRefNames.count(it->first) == 0) { //you are not a reference sequence
                    fitSeqsIndexes.insert(count); //add as fit seq
                }
            }
            it->second = count; count++;
            nameMap.push_back(it->first);
            nameAssignment[it->first] = it->second;
        }
        
        //read file to find singletons
        vector<bool> singleton; singleton.resize(count, true);
        map<long long, long long> singletonIndexSwap;
        
        if (distFormat == "column")        {  singletonIndexSwap = readColumnSingletons(singleton, distFile, nameAssignment);           }
        else if (distFormat == "phylip")   {  singletonIndexSwap = readPhylipSingletons(singleton, distFile, count, nameAssignment);    }
        
        int nonSingletonCount = 0;
        for (int i = 0; i < singleton.size(); i++) {
            if (!singleton[i]) { //if you are not a singleton
                singletonIndexSwap[i] = nonSingletonCount;
                nonSingletonCount++;
            }else { singletons.push_back(nameMap[i]); }
        }
        numSingletons = singletons.size();
        closeness.resize(nonSingletonCount);
        
        map<string, string> names;
        if (namefile != "") {
            //update names for reference
            util.readNames(namefile, names);
            for (int i = 0; i < numSingletons; i++) {
                singletons[i] = names[singletons[i]];
            }
        }
        
        //read reference file distances
        bool hasName = false;
        if (namefile != "") { hasName = true; }
        if (distFormat == "column")        {  readColumn(distFile, hasName, names, nameAssignment, singletonIndexSwap);     }
        else if (distFormat == "phylip")   {  readPhylip(distFile, hasName, names, nameAssignment, singletonIndexSwap);     }
        
        
        //randomly select the "fit" seqs
        long long numToSelect = nameAssignment.size() * fitPercent;
        if (weights.size() != 0) {  fitSeqsIndexes = subsample.getWeightedSample(weights, numToSelect);  } //you have weighted selection
        else {
            if (refWeightMethod == "accnos") { } //fitIndexes are filled above
            else { //randomly select references
                long long numSelected = 0;
                long long totalSeqs = nameAssignment.size();
                while (numSelected < numToSelect) {
                    if (m->getControl_pressed()) { break; }
                    fitSeqsIndexes.insert(util.getRandomIndex(totalSeqs-1)); //no repeats
                    numSelected = fitSeqsIndexes.size();
                }
            }
        }
        
        //flag reference seqs singleton or not
        for (int i = 0; i < singleton.size(); i++) {
            if (!singleton[i]) { //if you are not a singleton
                
                if (fitSeqsIndexes.count(i) != 0) { //you are a fit seq
                    isRef.push_back(false);
                }else { isRef.push_back(true);  } //its a reference
            }else {
                if (fitSeqsIndexes.count(i) != 0) { //you are a fit seq singleton
                    isSingleRef.push_back(false);
                }else { isSingleRef.push_back(true); } //its a singleton reference
            }
        }
        singleton.clear();
        
        //find number of fitDists, refDists and between dists
        calcCounts();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "readFiles");
        exit(1);
    }
}
/***********************************************************************/
//for reading reference and fit files separately, reference method
int OptiRefMatrix::readFiles(string refdistfile, string refnamefile, string refcountfile, string refformat, string refdistformat, string fitdistfile, string fitnamefile, string fitcountfile, string fitformat, string fitdistformat, string betweendistfile, string betweendistformat){
    try {
        map<string, long long> nameAssignment;
        if (refnamefile != "") { util.readNames(refnamefile, nameAssignment); }
        else  {
            CountTable ct; ct.readTable(refcountfile, false, true);
            map<string, int> temp = ct.getNameMap();
            for (map<string, int>::iterator it = temp.begin(); it!= temp.end(); it++) {  nameAssignment[it->first] = it->second; }
        }
        
        long long count = 0;
        for (map<string, long long>::iterator it = nameAssignment.begin(); it!= nameAssignment.end(); it++) {
            it->second = count; count++;
            nameMap.push_back(it->first);
            nameAssignment[it->first] = it->second;
        }
        
        long long refCount = count;
        vector<bool> singleton; singleton.resize(count, true); //resize will only set new elements to true
        map<long long, long long> refSingletonIndexSwap; //index into
        if (refdistformat == "column")        {  refSingletonIndexSwap = readColumnSingletons(singleton, refdistfile, nameAssignment);          }
        else if (refdistformat == "phylip")   {  refSingletonIndexSwap = readPhylipSingletons(singleton, refdistfile, count, nameAssignment);   }
        
        //read fit file to find singletons
        map<long long, long long> fitSingletonIndexSwap;
        map<string, long long> fitnameAssignment;
        if (fitnamefile != "") { util.readNames(fitnamefile, fitnameAssignment); }
        else  {
            CountTable ct; ct.readTable(fitcountfile, false, true);
            map<string, int> temp = ct.getNameMap();
            for (map<string, int>::iterator it = temp.begin(); it!= temp.end(); it++) {  fitnameAssignment[it->first] = it->second; }
        }
        
        for (map<string, long long>::iterator it = fitnameAssignment.begin(); it!= fitnameAssignment.end(); it++) {
            it->second = count; count++;
            nameMap.push_back(it->first);
            nameAssignment[it->first] = it->second;
        }

        singleton.resize(count, true);
        if (fitdistformat == "column")        {  fitSingletonIndexSwap = readColumnSingletons(singleton, fitdistfile, nameAssignment);          }
        else if (fitdistformat == "phylip")   {  fitSingletonIndexSwap = readPhylipSingletons(singleton, fitdistfile, count, nameAssignment);   }

        fitPercent = ((count-refCount) / (float) count);
        
        //read bewtween file to update singletons
        readColumnSingletons(singleton, betweendistfile, nameAssignment);
        
        long long nonSingletonCount = 0;
        map<long long, long long> singletonIndexSwap;
        for (long long i = 0; i < refCount; i++) {
            if (!singleton[i]) { //if you are not a singleton
                singletonIndexSwap[i] = nonSingletonCount;
                isRef.push_back(true);
                nonSingletonCount++;
            }else {
                singletons.push_back(nameMap[i]);
                isSingleRef.push_back(true);
            }
        }
        refSingletonIndexSwap.clear();
        
        for (long long i = refCount; i < singleton.size(); i++) {
            if (!singleton[i]) { //if you are not a singleton
                singletonIndexSwap[i] = nonSingletonCount;
                isRef.push_back(false);
                nonSingletonCount++;
            }else {
                singletons.push_back(nameMap[i]);
                isSingleRef.push_back(false);
            }
        }
        singleton.clear();
        fitSingletonIndexSwap.clear();
        
        numSingletons = singletons.size();
        closeness.resize(nonSingletonCount);
        
        map<string, string> names;
        if (refnamefile != "") { util.readNames(refnamefile, names); }
        
        if (fitnamefile != "") {
            map<string, string> fitnames;
            util.readNames(fitnamefile, fitnames);
            
            names.insert(fitnames.begin(), fitnames.end()); //copy fit names into names
        }
        
        if ((fitnamefile != "") || (refnamefile != "")) {
            for (int i = 0; i < singletons.size(); i++) {
                map<string, string>::iterator it = names.find(singletons[i]);
                if (it != names.end()) { //update singletons
                    singletons[i] = it->second;
                }
            }
        }
        
        //read reference file distances
        bool refHasName = false;
        if (refnamefile != "") { refHasName = true; }
        if (refdistformat == "column")        {  readColumn(refdistfile, refHasName, names, nameAssignment, singletonIndexSwap);     }
        else if (refdistformat == "phylip")   {  readPhylip(refdistfile, refHasName, names, nameAssignment, singletonIndexSwap);     }
        
        
        //read fit distances
        bool fitHasName = false;
        if (fitnamefile != "") { fitHasName = true; }
        if (fitdistformat == "column")        {  readColumn(fitdistfile, fitHasName, names, nameAssignment, singletonIndexSwap);     }
        else if (fitdistformat == "phylip")   {  readPhylip(fitdistfile, fitHasName, names, nameAssignment, singletonIndexSwap);     }
        
        
        //read in between distances
        bool hasName = fitHasName;
        if (!hasName && refHasName) { hasName = true; } //if either the ref or fit has a name file then set hasName
        if (betweendistformat == "column")        {  readColumn(betweendistfile, hasName, names, nameAssignment, singletonIndexSwap);     }
        else if (betweendistformat == "phylip")   {  readPhylip(betweendistfile, hasName, names, nameAssignment, singletonIndexSwap);     }
        
        //find number of fitDists, refDists and between dists
        calcCounts();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "readFiles");
        exit(1);
    }
}
/***********************************************************************/
map<long long, long long> OptiRefMatrix::readColumnSingletons(vector<bool>& singleton, string distFile, map<string, long long>& nameAssignment){
    try {
        
        ifstream fileHandle;
        util.openInputFile(distFile, fileHandle);
        
        string firstName, secondName;
        double distance;
        map<long long, long long> singletonIndexSwap;
        
        while(fileHandle){  //let's assume it's a triangular matrix...
            
            fileHandle >> firstName; util.gobble(fileHandle);
            fileHandle >> secondName; util.gobble(fileHandle);
            fileHandle >> distance;	util.gobble(fileHandle); // get the row and column names and distance
            
            if (m->getDebug()) { cout << firstName << '\t' << secondName << '\t' << distance << endl; }
            
            if (m->getControl_pressed()) {  break; }
            
            if (util.isEqual(distance,-1)) { distance = 1000000; }
            
            if(distance <= cutoff){
                map<string,long long>::iterator itA = nameAssignment.find(firstName);
                map<string,long long>::iterator itB = nameAssignment.find(secondName);
                
                if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the name or count file, please correct\n"); exit(1);  }
                if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the name or count file, please correct\n"); exit(1);  }
                
                long long indexA = (itA->second);
                long long indexB = (itB->second);
                singleton[indexA] = false;
                singleton[indexB] = false;
                singletonIndexSwap[indexA] = indexA;
                singletonIndexSwap[indexB] = indexB;
            }
        }
        fileHandle.close();
        
        return singletonIndexSwap;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "readColumnSingletons");
        exit(1);
    }
}
/***********************************************************************/

map<long long, long long> OptiRefMatrix::readPhylipSingletons(vector<bool>& singleton, string distFile, long long& count, map<string, long long>& nameAssignment){
    try {
        float distance;
        long long nseqs;
        string name;
        map<long long, long long> singletonIndexSwap;
        
        ifstream fileHandle;
        string numTest;
        
        util.openInputFile(distFile, fileHandle);
        fileHandle >> numTest >> name;
        nameMap.push_back(name);
        singletonIndexSwap[0] = 0;
        nameAssignment[name] = 0;
        
        if (!util.isContainingOnlyDigits(numTest)) { m->mothurOut("[ERROR]: expected a number and got " + numTest + ", quitting.\n"); m->setControl_pressed(true); return singletonIndexSwap; }
        else { convert(numTest, nseqs); }
        
        //square test
        char d;
        while((d=fileHandle.get()) != EOF){
            if(isalnum(d)){ square = true; fileHandle.putback(d); for(int i=0;i<nseqs;i++){ fileHandle >> distance;  } break; }
            if(d == '\n'){ square = false; break; }
        }
        
        singleton.resize((count+nseqs), true);
        if(square == 0){
            
            for(long long i=1;i<nseqs;i++){
                if (m->getControl_pressed()) {  break; }
                
                fileHandle >> name; nameMap.push_back(name); singletonIndexSwap[i] = i;  nameAssignment[name] = i;
                
                for(long long j=0;j<i;j++){
                    
                    fileHandle >> distance;
                    
                    if (util.isEqual(distance,-1)) { distance = 1000000; }
                    
                    if(distance <= cutoff){
                        singleton[i] = false;
                        singleton[j] = false;
                    }
                }
            }
        }else{
            for(long long i=1;i<nseqs;i++){
                if (m->getControl_pressed()) {  break; }
                
                fileHandle >> name; nameMap.push_back(name); singletonIndexSwap[i] = i; nameAssignment[name] = i;
                
                for(long long j=0;j<nseqs;j++){
                    fileHandle >> distance;
                    
                    if (util.isEqual(distance,-1)) { distance = 1000000; }
                    
                    if(distance <= cutoff && j < i){
                        singleton[i] = false;
                        singleton[j] = false;
                    }
                }
            }
        }
        fileHandle.close();
    
        count += nseqs;
        
        return singletonIndexSwap;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "readPhylipSingletons");
        exit(1);
    }
}
/***********************************************************************/
int OptiRefMatrix::readPhylip(string distFile, bool hasName, map<string, string>& names, map<string, long long>& nameAssignment, map<long long, long long>& singletonIndexSwap){
    try {
        long long nseqs;
        string name;
        double distance;
        
        ifstream in; string numTest;
        util.openInputFile(distFile, in);
        
        in >> numTest >> name;
        
        if (hasName) { name = names[name]; } //redundant names
        nameMap[singletonIndexSwap[0]] = name;
        
        
        if (!util.isContainingOnlyDigits(numTest)) { m->mothurOut("[ERROR]: expected a number and got " + numTest + ", quitting.\n"); m->setControl_pressed(true); return 0; }
        else { convert(numTest, nseqs); }
        
        //square test
        char d;
        while((d=in.get()) != EOF){
            if(isalnum(d)){ square = true; in.putback(d); for(int i=0;i<nseqs;i++){ in >> distance;  } break; }
            if(d == '\n'){ square = false; break; }
        }
        
        string line = "";
        if(!square){
            
            for(long long i=1;i<nseqs;i++){
                
                if (m->getControl_pressed()) {  break; }
                
                in >> name; util.gobble(in);
                
                if (hasName) { name = names[name]; } //redundant names
                nameMap[singletonIndexSwap[i]] = name;
                
                for(long long j=0;j<i;j++){
                    
                    in >> distance; util.gobble(in);
                    
                    if (util.isEqual(distance,-1)) { distance = 1000000; }
                    
                    if(distance <= cutoff){
                        if (refWeightMethod == "connectivity") { //count dists
                            weights[i]++; weights[j]++;
                        }
                        long long newB = singletonIndexSwap[j];
                        long long newA = singletonIndexSwap[i];
                        closeness[newA].insert(newB);
                        closeness[newB].insert(newA);
                    }
                }
            }
        }else{
            for(long long i=0;i<nseqs;i++){ in >> distance;  } util.gobble(in);
            
            for(long long i=1;i<nseqs;i++){
                if (m->getControl_pressed()) {  break; }
                
                in >> name; util.gobble(in);
                
                if (hasName) { name = names[name]; } //redundant names
                nameMap[singletonIndexSwap[i]] = name;
                
                for(long long j=0;j<nseqs;j++){
                    in >> distance; util.gobble(in);
                    
                    if (util.isEqual(distance,-1)) { distance = 1000000; }
                    
                    if(distance <= cutoff && j < i){
                        if (refWeightMethod == "connectivity") { //count dists
                            weights[i]++; weights[j]++;
                        }
                        long long newB = singletonIndexSwap[j];
                        long long newA = singletonIndexSwap[i];
                        closeness[newA].insert(newB);
                        closeness[newB].insert(newA);
                    }
                }
            }
        }
        in.close();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "readPhylip");
        exit(1);
    }
}
/***********************************************************************/

int OptiRefMatrix::readColumn(string distFile, bool hasName, map<string, string>& names, map<string, long long>& nameAssignment, map<long long, long long>& singletonIndexSwap){
    try {
        string firstName, secondName;
        double distance;
        
        ifstream in; util.openInputFile(distFile, in);
        
        while(in){  //let's assume it's a triangular matrix...
            
            in >> firstName; util.gobble(in);
            in >> secondName; util.gobble(in);
            in >> distance;	util.gobble(in); // get the row and column names and distance
            
            if (m->getDebug()) { cout << firstName << '\t' << secondName << '\t' << distance << endl; }
            
            if (m->getControl_pressed()) {  in.close();   return 0; }
            
            if (util.isEqual(distance,-1)) { distance = 1000000; }
            
            if(distance <= cutoff){
                map<string,long long>::iterator itA = nameAssignment.find(firstName);
                map<string,long long>::iterator itB = nameAssignment.find(secondName);
                
                if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the name or count file, please correct\n"); exit(1);  }
                if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the name or count file, please correct\n"); exit(1);  }
                
                long long indexA = (itA->second);
                long long indexB = (itB->second);
                
                if (refWeightMethod == "connectivity") { //count dists
                    weights[indexA]++; weights[indexB]++;
                }
                
                long long newB = singletonIndexSwap[indexB];
                long long newA = singletonIndexSwap[indexA];
                closeness[newA].insert(newB);
                closeness[newB].insert(newA);
                
                if (hasName) {
                    map<string, string>::iterator itName1 = names.find(firstName);
                    map<string, string>::iterator itName2 = names.find(secondName);
                    
                    if (itName1 != names.end()) { firstName = itName1->second;  } //redundant names
                    if (itName2 != names.end()) { secondName = itName2->second;  } //redundant names
                }
                
                nameMap[newA] = firstName;
                nameMap[newB] = secondName;
            }
        }
        in.close();
        
        return 1;
        
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "readColumn");
        exit(1);
    }
}
/***********************************************************************/

void OptiRefMatrix::calcCounts(){
    try {
        //find number of fitDists, refDists and between dists
        numRefDists = 0;
        numFitDists = 0;
        numBetweenDists = 0;
        numFitSingletons = 0;
        numFitSeqs = 0;
        numRefSingletons = 0;
        
        for (long long i = 0; i < closeness.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            bool thisSeqIsRef = isRef[i];
            long long thisSeqsNumRefDists = 0;
            long long thisSeqsNumFitDists = 0;
            
            for (set<long long>::iterator it = closeness[i].begin(); it != closeness[i].end(); it++) {
                long long newB = *it;
                
                if ((thisSeqIsRef) && (isRef[newB])) {  thisSeqsNumRefDists++; } //both refs
                else if ((thisSeqIsRef) && (!isRef[newB])) { numBetweenDists++; } // ref to fit dist
                else if ((!thisSeqIsRef) && (isRef[newB])) { numBetweenDists++; } // fit to ref dist
                else if ((!thisSeqIsRef) && (!isRef[newB])) { thisSeqsNumFitDists++; } // both fit
            }
            
            //a refSingleton or Fitsingleton may not be a true singleton (no valid dists in matrix), but may be a refSeq with no distances to other refs but distances to fitseqs. a fitsingleton may have dists to refs but no dists to other fitseqs.
            
            //you are a ref with no refdists, so you are a refsingleton
            if ((thisSeqIsRef) && (thisSeqsNumRefDists == 0)) {  numRefSingletons++; }
            else if ((!thisSeqIsRef) && (thisSeqsNumFitDists == 0)) {  numFitSingletons++; }
            else if ((!thisSeqIsRef) && (thisSeqsNumFitDists != 0)) {  numFitSeqs++; }
            
            numRefDists += thisSeqsNumRefDists;
            numFitDists += thisSeqsNumFitDists;
        }
        
        //counted twice
        numRefDists /= 2;
        numFitDists /= 2;
        numBetweenDists /= 2;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "calcCounts");
        exit(1);
    }
}
/***********************************************************************/


