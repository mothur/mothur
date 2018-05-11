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
OptiRefMatrix::OptiRefMatrix(string d, string nc, string f, string df, double c, string fit, string fitnc, string fitf, string fitdf, string betweend, string betweendf) : refdistfile(d), refdistformat(df), refformat(f), fitdistfile(fit), fitdistformat(fitdf), fitformat(fitf), betweendistfile(betweend), betweendistformat(betweendf), OptiData(c) {
    
    numFitSingletons = 0;
    numRefSingletons = 0;
    numSingletons = 0;
    numBetweenDists = 0;
    numFitDists = 0;
    numRefDists = 0;
    
    square = false;
    
    if (refformat == "name") { refnamefile = nc; refcountfile = ""; }
    else if (refformat == "count") { refcountfile = nc; refnamefile = ""; }
    else { refcountfile = ""; refnamefile = ""; }
    
    if (fitformat == "name") { fitnamefile = fitnc; fitcountfile = ""; }
    else if (fitformat == "count") { fitcountfile = fitnc; fitnamefile = ""; }
    else { fitcountfile = ""; fitnamefile = ""; }
    
    readFiles();
}
/***********************************************************************/
vector<int> OptiRefMatrix::getNumSeqs(vector<vector<string> > & binNames, vector<vector<int> > & fixedBins) {
    try {
        fixedBins.clear();
        
        map<string, int> nameIndexes;
        set<string> unique;
        for (int i = 0; i < nameMap.size(); i++) { //vector of string representing the name file.
            vector<string> thisSeqsReps; util.splitAtComma(nameMap[i], thisSeqsReps); //split redundant names
            if (i < closeness.size()) {  nameIndexes[thisSeqsReps[0]] = i;  } //this is a sequence with distances in the matrix
            if (thisSeqsReps.size() == 1) { //you are unique
                unique.insert(thisSeqsReps[0]);
            }
        }
        
        for (int i = 0; i < binNames.size(); i++) { //for each OTU
            vector<int> thisBinsSeqs;
            for (int j = 0; j < binNames[i].size(); j++) { //for each sequence
                map<string, int>::iterator it = nameIndexes.find(binNames[i][j]);
                
                if (it == nameIndexes.end()) {//not in distance matrix, but needs a value in fixedBins. 2 reasons for making it here: you are a redundant name in the listfile, you do not have any distances below the cutoff
                    if (unique.count(binNames[i][j]) == 0) { } //you are redundant seq in list file because namefile was used. You should be edited out of the listfile.
                    else { thisBinsSeqs.push_back(-1); } //you are unique, but have no distances so add placeholder
                }else { thisBinsSeqs.push_back(it->second);  nameIndexes.erase(it); } //"name" of sequence in matrix
            }
            fixedBins.push_back(thisBinsSeqs);
        }
        
        return (getFitSeqs());
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "getNumSeqs");
        exit(1);
    }
}
/***********************************************************************/
bool OptiRefMatrix::isCloseFit(int i, int toFind){
    try {
        if (i > closeness.size()) { m->mothurOut("[ERROR]: index is not valid.\n"); m->setControl_pressed(true); return false; }
        else if (index < 0) { return false; }
        
        bool found = false;
        if (toFind >= refEnd) { //are you a fit seq
            if (closeness[i].count(toFind) != 0) {  //are you close
                found = true;
            }
        }
        return found;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "isCloseFit");
        exit(1);
    }
}
/***********************************************************************/
vector<int> OptiRefMatrix::getRefSeqs() {
    try {
        vector<int> refSeqsIndexes;
        for (int i = 0; i < refEnd; i++) { refSeqsIndexes.push_back(i); }
        return refSeqsIndexes;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "getRefSeqs");
        exit(1);
    }
}
/***********************************************************************/
vector<int> OptiRefMatrix::getFitSeqs() {
    try {
        vector<int> fitSeqsIndexes;
        for (int i = refEnd; i < closeness.size(); i++) { fitSeqsIndexes.push_back(i); }
        return fitSeqsIndexes;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "getFitSeqs");
        exit(1);
    }
}
/***********************************************************************/
int OptiRefMatrix::getNumFitClose(int index) {
    try {
        int numClose = 0;
        
        if (index > closeness.size()) { m->mothurOut("[ERROR]: index is not valid.\n"); m->setControl_pressed(true);  }
        else if (index < 0) { }
        else {
            //reference seqs all have indexes less than refEnd
            for (set<int>::iterator it = closeness[index].begin(); it != closeness[index].end(); it++) {
                if ((*it) >= refEnd) {  numClose++; } //you are a fit seq
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
int OptiRefMatrix::getNumRefClose(int index) {
    try {
        int numClose = 0;
        
        if (index > closeness.size()) { m->mothurOut("[ERROR]: index is not valid.\n"); m->setControl_pressed(true);  }
        else if (index < 0) { }
        else {
            //reference seqs all have indexes less than refEnd
            for (set<int>::iterator it = closeness[index].begin(); it != closeness[index].end(); it++) {
                if ((*it) < refEnd) {  numClose++; } //you are a ref seq
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
set<int> OptiRefMatrix::getCloseFitSeqs(int index){
    try {
        set<int> closeSeqs;
        
        if (index > closeness.size()) { m->mothurOut("[ERROR]: index is not valid.\n"); m->setControl_pressed(true);  }
        else if (index < 0) { } //
        else {
            //reference seqs all have indexes less than refEnd
            for (set<int>::iterator it = closeness[index].begin(); it != closeness[index].end(); it++) {
                if ((*it) >= refEnd) {  closeSeqs.insert(*it); } //you are a fit seq
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
set<int> OptiRefMatrix::getCloseRefSeqs(int index){
    try {
        set<int> closeSeqs;
        
        if (index > closeness.size()) { m->mothurOut("[ERROR]: index is not valid.\n"); m->setControl_pressed(true);  }
        else if (index < 0) { } //
        else {
            //reference seqs all have indexes less than refEnd
            for (set<int>::iterator it = closeness[index].begin(); it != closeness[index].end(); it++) {
                if ((*it) < refEnd) {  closeSeqs.insert(*it); } //you are a ref seq
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
            
            for (int i = numRefSingletons; i < singletons.size(); i++) { singlelist->push_back(singletons[i]); }
        }
        
        return singlelist;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "getFitListSingle");
        exit(1);
    }
}
/***********************************************************************/
//lets create separate reads for single and combined files.  Will produce some dup code
int OptiRefMatrix::readFiles(){
    try {
        vector<bool> singleton;
        //read reference file to find singletons
        map<int, int> refSingletonIndexSwap;
        map<string, int> nameAssignment;
        int count = 0;
        
        if (refdistformat == "column")        {  refSingletonIndexSwap = readColumnSingletons(singleton, refnamefile, refcountfile, refdistfile, count, nameAssignment);     }
        else if (refdistformat == "phylip")   {  refSingletonIndexSwap = readPhylipSingletons(singleton, refdistfile, count, nameAssignment);                                }
        
        int numRefSeqs = count+1;
        
        //read fit file to find singletons
        map<int, int> fitSingletonIndexSwap;
        
        if (fitdistformat == "column")        {  fitSingletonIndexSwap = readColumnSingletons(singleton, fitnamefile, fitcountfile, fitdistfile, count, nameAssignment);     }
        else if (fitdistformat == "phylip")   {  fitSingletonIndexSwap = readPhylipSingletons(singleton, fitdistfile, count, nameAssignment);                                }
        
        int numFitSeqs = (count+1) - numRefSeqs;
        
        //read bewtween file to update singletons
        readColumnSingletons(singleton, betweendistfile, nameAssignment);
        
        int nonSingletonCount = 0;
        map<int, int> singletonIndexSwap;
        for (int i = 0; i < numRefSeqs; i++) {
            if (!singleton[i]) { //if you are not a singleton
                singletonIndexSwap[i] = nonSingletonCount;
                nonSingletonCount++;
            }else { singletons.push_back(nameMap[i]); }
        }
        refSingletonIndexSwap.clear();
        numRefSingletons = singletons.size();
        refEnd = nonSingletonCount; // reference sequences are stored in beginning of closeness, fit seqs stored after
        
        int numSeqs = numFitSeqs+numRefSeqs;
        for (int i = numRefSeqs; i < numSeqs; i++) {
            if (!singleton[i]) { //if you are not a singleton
                singletonIndexSwap[i] = nonSingletonCount;
                nonSingletonCount++;
            }else { singletons.push_back(nameMap[i]); }
        }
        singleton.clear();
        fitSingletonIndexSwap.clear();
        
        numSingletons = singletons.size();
        numFitSingletons = numSingletons - numRefSingletons;
        closeness.resize(nonSingletonCount);
        
        map<string, string> names;
        if (refnamefile != "") {
            
            //update names for reference
            util.readNames(refnamefile, names);
            for (int i = 0; i < numRefSeqs; i++) {
                singletons[i] = names[singletons[i]];
            }
        }
        
        if (fitnamefile != "") {
            
            map<string, string> fitnames;
            //update names for fit seqs
            util.readNames(fitnamefile, fitnames);
            for (int i = numRefSeqs; i < numSeqs; i++) {
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
        if (fitdistformat == "column")        {  readColumn(fitdistfile, hasName, names, nameAssignment, singletonIndexSwap);     }
        else if (fitdistformat == "phylip")   {  readPhylip(fitdistfile, hasName, names, nameAssignment, singletonIndexSwap);     }
        
    }
    catch(exception& e) {
        m->errorOut(e, "OptiRefMatrix", "readFiles");
        exit(1);
    }
}
/***********************************************************************/
//when running in combined mode the singleton parameter must have room for all the reads
//reading the in between file, so we already know "who" is in the file and that it's in column format since mothur created it.
void OptiRefMatrix::readColumnSingletons(vector<bool>& singleton, string distFile, map<string, int>& nameAssignment){
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
                map<string,int>::iterator itA = nameAssignment.find(firstName);
                map<string,int>::iterator itB = nameAssignment.find(secondName);
                
                if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the name or count file, please correct\n"); exit(1);  }
                if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the name or count file, please correct\n"); exit(1);  }
                
                int indexA = (itA->second);
                int indexB = (itB->second);
                
                singleton[indexA] = false;
                singleton[indexB] = false;
                numBetweenDists++;
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
map<int, int> OptiRefMatrix::readColumnSingletons(vector<bool>& singleton, string namefile, string countfile, string distFile, int& count, map<string, int>& nameAssignment){
    try {
        bool readingRef = false;
        if (count == 0) { readingRef = true; }
        
        map<string, int> thisFilesNameAssignment;
        if (namefile != "") { thisFilesNameAssignment = util.readNames(namefile); }
        else  {  CountTable ct; ct.readTable(countfile, false, true); thisFilesNameAssignment = ct.getNameMap(); }
            
        for (map<string, int>::iterator it = thisFilesNameAssignment.begin(); it!= thisFilesNameAssignment.end(); it++) {
            it->second = count; count++;
            nameMap.push_back(it->first);
            nameAssignment[it->first] = it->second; //add thisFilesNameAssignment to nameAssignment
        }
        
        ifstream fileHandle;
        util.openInputFile(distFile, fileHandle);
        
        string firstName, secondName;
        double distance;
        map<int, int> singletonIndexSwap;
        singleton.resize(count, true); //resize will only set new elements to true
        
        while(fileHandle){  //let's assume it's a triangular matrix...
            
            fileHandle >> firstName; util.gobble(fileHandle);
            fileHandle >> secondName; util.gobble(fileHandle);
            fileHandle >> distance;	util.gobble(fileHandle); // get the row and column names and distance
            
            if (m->getDebug()) { cout << firstName << '\t' << secondName << '\t' << distance << endl; }
            
            if (m->getControl_pressed()) {  break; }
            
            if (distance == -1) { distance = 1000000; }
            
            if(distance < cutoff){
                map<string,int>::iterator itA = nameAssignment.find(firstName);
                map<string,int>::iterator itB = nameAssignment.find(secondName);
                
                if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the name or count file, please correct\n"); exit(1);  }
                if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the name or count file, please correct\n"); exit(1);  }
                
                int indexA = (itA->second);
                int indexB = (itB->second);
                singleton[indexA] = false;
                singleton[indexB] = false;
                singletonIndexSwap[indexA] = indexA;
                singletonIndexSwap[indexB] = indexB;
                
                if (readingRef) { numRefDists++; }
                else { numFitDists++; }
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

map<int, int> OptiRefMatrix::readPhylipSingletons(vector<bool>& singleton, string distFile, int& count, map<string, int>& nameAssignment){
    try {
        bool readingRef = false;
        if (count == 0) { readingRef = true; }
        
        float distance;
        int nseqs;
        string name;
        map<int, int> singletonIndexSwap;
        
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
            
            for(int i=1;i<nseqs;i++){
                if (m->getControl_pressed()) {  break; }
                
                fileHandle >> name; nameMap.push_back(name); singletonIndexSwap[i] = i;  nameAssignment[name] = i;
                
                for(int j=0;j<i;j++){
                    
                    fileHandle >> distance;
                    
                    if (distance == -1) { distance = 1000000; }
                    
                    if(distance < cutoff){
                        singleton[i] = false;
                        singleton[j] = false;
                        
                        if (readingRef) { numRefDists++; }
                        else { numFitDists++; }
                    }
                }
            }
        }else{
            for(int i=1;i<nseqs;i++){
                if (m->getControl_pressed()) {  break; }
                
                fileHandle >> name; nameMap.push_back(name); singletonIndexSwap[i] = i; nameAssignment[name] = i;
                
                for(int j=0;j<nseqs;j++){
                    fileHandle >> distance;
                    
                    if (distance == -1) { distance = 1000000; }
                    
                    if(distance < cutoff && j < i){
                        singleton[i] = false;
                        singleton[j] = false;
                        
                        if (readingRef) { numRefDists++; }
                        else { numFitDists++; }
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
int OptiRefMatrix::readPhylip(string distFile, bool hasName, map<string, string>& names, map<string, int>& nameAssignment, map<int, int>& singletonIndexSwap){
    try {
        int nseqs;
        string name;
        double distance;
        
        ifstream in;
        util.openInputFile(distFile, in);
        
        in >> nseqs >> name;
        
        if (hasName) { name = names[name]; } //redundant names
        nameMap[singletonIndexSwap[0]] = name;
        
        string line = "";
        if(!square){
            int index = 0;
            
            for(int i=1;i<nseqs;i++){
                
                if (m->getControl_pressed()) {  break; }
                
                in >> name; util.gobble(in);
                
                if (hasName) { name = names[name]; } //redundant names
                nameMap[singletonIndexSwap[i]] = name;
                
                for(int j=0;j<i;j++){
                    
                    in >> distance; util.gobble(in);
                    
                    if (distance == -1) { distance = 1000000; } 
                    
                    if(distance < cutoff){
                        int newB = singletonIndexSwap[j];
                        int newA = singletonIndexSwap[i];
                        closeness[newA].insert(newB);
                        closeness[newB].insert(newA);
                    }
                    index++;
                }
            }
        }else{
            int index = nseqs;
            
            for(int i=0;i<nseqs;i++){ in >> distance;  } util.gobble(in);
            
            for(int i=1;i<nseqs;i++){
                if (m->getControl_pressed()) {  break; }
                
                in >> name; util.gobble(in);
                
                if (hasName) { name = names[name]; } //redundant names
                nameMap[singletonIndexSwap[i]] = name;
                
                for(int j=0;j<nseqs;j++){
                    in >> distance; util.gobble(in);
                    
                    if (distance == -1) { distance = 1000000; }
                    
                    if(distance < cutoff && j < i){
                        int newB = singletonIndexSwap[j];
                        int newA = singletonIndexSwap[i];
                        closeness[newA].insert(newB);
                        closeness[newB].insert(newA);
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

int OptiRefMatrix::readColumn(string distFile, bool hasName, map<string, string>& names, map<string, int>& nameAssignment, map<int, int>& singletonIndexSwap){
    try {
        string firstName, secondName;
        double distance;
        
        ifstream in;
        util.openInputFile(distFile, in);
        
        while(in){  //let's assume it's a triangular matrix...
            
            in >> firstName; util.gobble(in);
            in >> secondName; util.gobble(in);
            in >> distance;	util.gobble(in); // get the row and column names and distance
            
            if (m->getDebug()) { cout << firstName << '\t' << secondName << '\t' << distance << endl; }
            
            if (m->getControl_pressed()) {  in.close();   return 0; }
            
            if (distance == -1) { distance = 1000000; }
            
            if(distance < cutoff){
                map<string,int>::iterator itA = nameAssignment.find(firstName);
                map<string,int>::iterator itB = nameAssignment.find(secondName);
                
                if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the name or count file, please correct\n"); exit(1);  }
                if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the name or count file, please correct\n"); exit(1);  }
                
                int indexA = (itA->second);
                int indexB = (itB->second);
                
                int newB = singletonIndexSwap[indexB];
                int newA = singletonIndexSwap[indexA];
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


