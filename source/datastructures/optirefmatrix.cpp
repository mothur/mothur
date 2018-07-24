//
//  optirefmatrix.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/3/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "optirefmatrix.hpp"
#include "progress.hpp"
#include "counttable.h"

/***********************************************************************/
OptiRefMatrix::OptiRefMatrix(string distFile, string distFormat, string dupsFile, string dupsFormat, double c, float fP) : OptiData(c) {
    
    numFitSingletons = 0;
    numRefSingletons = 0;
    numSingletons = 0;
    numBetweenDists = 0;
    numFitDists = 0;
    numRefDists = 0;
    numFitSeqs = 0;
    
    fitPercent = fP / 100.0;
    if (fitPercent < 0.001) { fitPercent = 0.10; m->mothurOut("[WARNING]: fit percentage must be between 0.001 (0.1%) and 1.0 (100%). Setting to 0.10 or 10%. \n"); } //minumum of 0.1%
    else if (fitPercent > 100.0) {  m->mothurOut("[ERROR]: fit percentage must be between 0.0001 and 100.0\n"); m->setControl_pressed(true); }
    
    square = false;
    
    readFiles(distFile, distFormat, dupsFile, dupsFormat);
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
        
        for (long long i = 0; i < numRefSingletons; i++) {
            vector<string> thisSeqsReps; util.splitAtComma(singletons[i], thisSeqsReps); //split redundant names
            nameIndexes[thisSeqsReps[0]] = -1;
            if (thisSeqsReps.size() == 1) { unique.insert(thisSeqsReps[0]); }
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
        long long numSelected = 0;
        long long refSingletonCutoff = isRef.size();
        long long singleSize = isSingleRef.size();
        
        set<long long> fitSeqsIndexes;
        while (numSelected < numToSelect) {
            if (m->getControl_pressed()) { break; }
            
            fitSeqsIndexes.insert(util.getRandomIndex(totalSeqs-1)); //no repeats
            
            numSelected = fitSeqsIndexes.size();
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
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "randomizeRefs");
        exit(1);
    }

}
/***********************************************************************/
int OptiRefMatrix::readFiles(string distFile, string distFormat, string dupsFile, string dupsFormat) {
    try {
        string namefile, countfile;
        if (dupsFormat == "name") { namefile = dupsFile; countfile = ""; }
        else if (dupsFormat == "count") { countfile = dupsFile; namefile = ""; }
        else { countfile = ""; namefile = ""; }
        
        //read file to find singletons
        vector<bool> singleton;
        map<long long, long long> singletonIndexSwap;
        map<string, long long> nameAssignment;
        long long count = 0;
        
        if (distFormat == "column")        {  singletonIndexSwap = readColumnSingletons(singleton, namefile, countfile, distFile, count, nameAssignment);     }
        else if (distFormat == "phylip")   {  singletonIndexSwap = readPhylipSingletons(singleton, distFile, count, nameAssignment);                          }
        
        //randomly select the "fit" seqs
        long long numToSelect = nameAssignment.size() * fitPercent;
        long long numSelected = 0;
        long long totalSeqs = nameAssignment.size();
        set<long long> fitSeqsIndexes;
        while (numSelected < numToSelect) {
            if (m->getControl_pressed()) { break; }
            
            fitSeqsIndexes.insert(util.getRandomIndex(totalSeqs-1)); //no repeats
            
            numSelected = fitSeqsIndexes.size();
        }
                
        int nonSingletonCount = 0;
        for (int i = 0; i < singleton.size(); i++) {
            if (!singleton[i]) { //if you are not a singleton
                singletonIndexSwap[i] = nonSingletonCount;
                
                if (fitSeqsIndexes.count(i) != 0) { //you are a fit seq
                    isRef.push_back(false);
                }else { isRef.push_back(true);  } //its a reference

                nonSingletonCount++;
            }else {
                singletons.push_back(nameMap[i]);
                
                if (fitSeqsIndexes.count(i) != 0) { //you are a fit seq singleton
                    isSingleRef.push_back(false);
                }else { isSingleRef.push_back(true); } //its a singleton reference
            }
        }
        
        singleton.clear();
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
        
        //find number of fitDists, refDists and between dists
        numRefDists = 0;
        numFitDists = 0;
        numBetweenDists = 0;
        numRefSingletons = 0;
        numFitSingletons = 0;
        numFitSeqs = 0;
        
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


        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "readFiles");
        exit(1);
    }
}
/***********************************************************************/
//lets create separate reads for single and combined files.  Will produce some dup code
int OptiRefMatrix::readFiles(string refdistfile, string refnamefile, string refcountfile, string refformat, string refdistformat, string fitdistfile, string fitnamefile, string fitcountfile, string fitformat, string fitdistformat, string betweendistfile, string betweendistformat){
    try {
        vector<bool> singleton;
        //read reference file to find singletons
        map<long long, long long> refSingletonIndexSwap;
        map<string, long long> nameAssignment;
        long long count = 0;
        
        if (refdistformat == "column")        {  refSingletonIndexSwap = readColumnSingletons(singleton, refnamefile, refcountfile, refdistfile, count, nameAssignment);     }
        else if (refdistformat == "phylip")   {  refSingletonIndexSwap = readPhylipSingletons(singleton, refdistfile, count, nameAssignment);                                }
        
        long long refSeqsEnd = count;
        
        numRefSingletons = 0;
        for (long long i = 0; i < refSeqsEnd; i++) { if (singleton[i]) { numRefSingletons++; } }
        
        //read fit file to find singletons
        map<long long, long long> fitSingletonIndexSwap;
        
        if (fitdistformat == "column")        {  fitSingletonIndexSwap = readColumnSingletons(singleton, fitnamefile, fitcountfile, fitdistfile, count, nameAssignment);     }
        else if (fitdistformat == "phylip")   {  fitSingletonIndexSwap = readPhylipSingletons(singleton, fitdistfile, count, nameAssignment);                                }
        
        numFitSingletons = 0;
        for (long long i = refSeqsEnd; i < singleton.size(); i++) {  if (singleton[i]) { numFitSingletons++; }  }

        fitPercent = ((count-refSeqsEnd) / (float) count);
        if (fitPercent < 0.001) { fitPercent = 0.001; } //minumum of 0.1%
        
        //read bewtween file to update singletons
        readColumnSingletons(singleton, betweendistfile, nameAssignment);
        
        numFitSeqs = 0;
        for (long long i = refSeqsEnd; i < singleton.size(); i++) {  if (!singleton[i]) { numFitSeqs++; }  }
        
        long long nonSingletonCount = 0;
        map<long long, long long> singletonIndexSwap;
        for (long long i = 0; i < refSeqsEnd; i++) {
            if (!singleton[i]) { //if you are not a singleton
                singletonIndexSwap[i] = nonSingletonCount;
                nonSingletonCount++;
            }else { singletons.push_back(nameMap[i]); }
        }
        refSingletonIndexSwap.clear();
        long long refEnd = nonSingletonCount; // reference sequences are stored in beginning of closeness, fit seqs stored after
        long long refSingletonsEnd = singletons.size();
        
        for (long long i = refSeqsEnd; i < singleton.size(); i++) {
            if (!singleton[i]) { //if you are not a singleton
                singletonIndexSwap[i] = nonSingletonCount;
                nonSingletonCount++;
            }else { singletons.push_back(nameMap[i]); }
        }
        singleton.clear();
        fitSingletonIndexSwap.clear();
        
        numSingletons = singletons.size();
        closeness.resize(nonSingletonCount);
        isRef.resize(nonSingletonCount, true);
        for (long long i = refEnd; i < nonSingletonCount; i++) { isRef[i] = false; } //you are a fitseq
        isSingleRef.resize(numSingletons, true);
        for (long long i = refSingletonsEnd; i < numSingletons; i++) { isSingleRef[i] = false; } //you are a fitseq
        
        map<string, string> names;
        if (refnamefile != "") {
            //update names for reference
            util.readNames(refnamefile, names);
            for (int i = 0; i < numRefSingletons; i++) {
                singletons[i] = names[singletons[i]];
            }
        }
        
        if (fitnamefile != "") {
            
            map<string, string> fitnames;
            //update names for fit seqs
            util.readNames(fitnamefile, fitnames);
            for (long long i = numRefSingletons; i < numSingletons; i++) {
                singletons[i] = fitnames[singletons[i]];
            }
            
            //copy fit names into names
            names.insert(fitnames.begin(), fitnames.end());
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
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "readFiles");
        exit(1);
    }
}
/***********************************************************************/
//when running in combined mode the singleton parameter must have room for all the reads
//reading the in between file, so we already know "who" is in the file and that it's in column format since mothur created it.
void OptiRefMatrix::readColumnSingletons(vector<bool>& singleton, string distFile, map<string, long long>& nameAssignment){
    try {
        string firstName, secondName;
        double distance;
        
        ifstream fileHandle;
        util.openInputFile(distFile, fileHandle);
        
        while(fileHandle){  //let's assume it's a triangular matrix...
            
            fileHandle >> firstName; util.gobble(fileHandle);
            fileHandle >> secondName; util.gobble(fileHandle);
            fileHandle >> distance;	util.gobble(fileHandle); // get the row and column names and distance
            
            if (m->getDebug()) { cout << firstName << '\t' << secondName << '\t' << distance << endl; }
            
            if (m->getControl_pressed()) {  break; }
            
            if (distance == -1) { distance = 1000000; }
            
            if(distance < cutoff){
                map<string,long long>::iterator itA = nameAssignment.find(firstName);
                map<string,long long>::iterator itB = nameAssignment.find(secondName);
                
                if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the name or count file, please correct\n"); exit(1);  }
                if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the name or count file, please correct\n"); exit(1);  }
                
                long long indexA = (itA->second);
                long long indexB = (itB->second);
                
                singleton[indexA] = false;
                singleton[indexB] = false;
            }
        }
        fileHandle.close();
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "readColumnSingletons");
        exit(1);
    }
}
/***********************************************************************/
map<long long, long long> OptiRefMatrix::readColumnSingletons(vector<bool>& singleton, string namefile, string countfile, string distFile, long long& count, map<string, long long>& nameAssignment){
    try {
        bool readingRef = false;
        if (count == 0) { readingRef = true; }
        
        map<string, long long> thisFilesNameAssignment;
        if (namefile != "") { util.readNames(namefile, thisFilesNameAssignment); }
        else  {
            CountTable ct; ct.readTable(countfile, false, true);
            map<string, int> temp = ct.getNameMap();
            for (map<string, int>::iterator it = temp.begin(); it!= temp.end(); it++) {  thisFilesNameAssignment[it->first] = it->second; }
        }
            
        for (map<string, long long>::iterator it = thisFilesNameAssignment.begin(); it!= thisFilesNameAssignment.end(); it++) {
            it->second = count; count++;
            nameMap.push_back(it->first);
            nameAssignment[it->first] = it->second; //add thisFilesNameAssignment to nameAssignment
        }
        
        ifstream fileHandle;
        util.openInputFile(distFile, fileHandle);
        
        string firstName, secondName;
        double distance;
        map<long long, long long> singletonIndexSwap;
        singleton.resize(count, true); //resize will only set new elements to true
        
        while(fileHandle){  //let's assume it's a triangular matrix...
            
            fileHandle >> firstName; util.gobble(fileHandle);
            fileHandle >> secondName; util.gobble(fileHandle);
            fileHandle >> distance;	util.gobble(fileHandle); // get the row and column names and distance
            
            if (m->getDebug()) { cout << firstName << '\t' << secondName << '\t' << distance << endl; }
            
            if (m->getControl_pressed()) {  break; }
            
            if (distance == -1) { distance = 1000000; }
            
            if(distance < cutoff){
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
        bool readingRef = false;
        if (count == 0) { readingRef = true; }
        
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
                    
                    if (distance == -1) { distance = 1000000; }
                    
                    if(distance < cutoff){
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
                    
                    if (distance == -1) { distance = 1000000; }
                    
                    if(distance < cutoff && j < i){
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
        
        ifstream in;
        util.openInputFile(distFile, in);
        
        in >> nseqs >> name;
        
        if (hasName) { name = names[name]; } //redundant names
        nameMap[singletonIndexSwap[0]] = name;
        
        string line = "";
        if(!square){
            long long index = 0;
            
            for(long long i=1;i<nseqs;i++){
                
                if (m->getControl_pressed()) {  break; }
                
                in >> name; util.gobble(in);
                
                if (hasName) { name = names[name]; } //redundant names
                nameMap[singletonIndexSwap[i]] = name;
                
                for(long long j=0;j<i;j++){
                    
                    in >> distance; util.gobble(in);
                    
                    if (distance == -1) { distance = 1000000; } 
                    
                    if(distance < cutoff){
                        long long newB = singletonIndexSwap[j];
                        long long newA = singletonIndexSwap[i];
                        closeness[newA].insert(newB);
                        closeness[newB].insert(newA);
                        
                        if ((isRef[newA]) && (isRef[newB])) {  numRefDists++; } //both refs
                        else if ((isRef[newA]) && (!isRef[newB])) { numBetweenDists++; } // ref to fit dist
                        else if ((!isRef[newA]) && (isRef[newB])) { numBetweenDists++; } // fit to ref dist
                        else if ((!isRef[newA]) && (!isRef[newB])) { numFitDists++; } // both fit
                        
                    }
                    index++;
                }
            }
        }else{
            long long index = nseqs;
            
            for(long long i=0;i<nseqs;i++){ in >> distance;  } util.gobble(in);
            
            for(long long i=1;i<nseqs;i++){
                if (m->getControl_pressed()) {  break; }
                
                in >> name; util.gobble(in);
                
                if (hasName) { name = names[name]; } //redundant names
                nameMap[singletonIndexSwap[i]] = name;
                
                for(long long j=0;j<nseqs;j++){
                    in >> distance; util.gobble(in);
                    
                    if (distance == -1) { distance = 1000000; }
                    
                    if(distance < cutoff && j < i){
                        long long newB = singletonIndexSwap[j];
                        long long newA = singletonIndexSwap[i];
                        closeness[newA].insert(newB);
                        closeness[newB].insert(newA);
                        
                        if ((isRef[newA]) && (isRef[newB])) {  numRefDists++; } //both refs
                        else if ((isRef[newA]) && (!isRef[newB])) { numBetweenDists++; } // ref to fit dist
                        else if ((!isRef[newA]) && (isRef[newB])) { numBetweenDists++; } // fit to ref dist
                        else if ((!isRef[newA]) && (!isRef[newB])) { numFitDists++; } // both fit

                    }
                    index++;
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
            
            if (distance == -1) { distance = 1000000; }
            
            if(distance < cutoff){
                map<string,long long>::iterator itA = nameAssignment.find(firstName);
                map<string,long long>::iterator itB = nameAssignment.find(secondName);
                
                if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the name or count file, please correct\n"); exit(1);  }
                if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the name or count file, please correct\n"); exit(1);  }
                
                long long indexA = (itA->second);
                long long indexB = (itB->second);
                
                long long newB = singletonIndexSwap[indexB];
                long long newA = singletonIndexSwap[indexA];
                closeness[newA].insert(newB);
                closeness[newB].insert(newA);
                
                if ((isRef[newA]) && (isRef[newB])) {  numRefDists++; } //both refs
                else if ((isRef[newA]) && (!isRef[newB])) { numBetweenDists++; } // ref to fit dist
                else if ((!isRef[newA]) && (isRef[newB])) { numBetweenDists++; } // fit to ref dist
                else if ((!isRef[newA]) && (!isRef[newB])) { numFitDists++; } // both fit
                
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


