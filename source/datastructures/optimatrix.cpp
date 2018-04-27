//
//  optimatrix.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/20/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#include "optimatrix.h"
#include "progress.hpp"
#include "counttable.h"

/***********************************************************************/

OptiMatrix::OptiMatrix(string d, string df, double c, bool s) : distFile(d), distFormat(df), cutoff(c), sim(s) {
    m = MothurOut::getInstance();
    countfile = ""; namefile = "";
    numMovableSingletons = 0;
    
    setBlastVariables(5, 0.10, true);
    if (distFormat == "phylip")         { readPhylip();     }
    else if (distFormat == "column")    { readColumn();     }
    else if (distFormat == "blast")     {   readBlast();    }
}
/***********************************************************************/
OptiMatrix::OptiMatrix(string d, string nc, string f, string df, double c, bool s) : distFile(d), distFormat(df), format(f), cutoff(c), sim(s) {
    m = MothurOut::getInstance();
    numMovableSingletons = 0;
    
    if (format == "name") { namefile = nc; countfile = ""; }
    else if (format == "count") { countfile = nc; namefile = ""; }
    else { countfile = ""; namefile = ""; }
    
    setBlastVariables(5, 0.10, true);
    if (distFormat == "phylip")         { readPhylip();     }
    else if (distFormat == "column")    { readColumn();     }
    else if (distFormat == "blast")     {   readBlast();    }
    
}
/***********************************************************************/
int OptiMatrix::readFile(string d, string nc, string f, string df, double c, bool s)  {
    distFile = d; format = f; cutoff = c; sim = s; distFormat = df;
    numMovableSingletons = 0;
    
    if (format == "name") { namefile = nc; countfile = ""; }
    else if (format == "count") { countfile = nc; namefile = ""; }
    else { countfile = ""; namefile = ""; }
    
    setBlastVariables(5, 0.10, true);
    if (distFormat == "phylip")         { readPhylip();     }
    else if (distFormat == "column")    { readColumn();     }
    else if (distFormat == "blast")     {   readBlast();    }
    
    return 0;
}
/***********************************************************************/
ListVector* OptiMatrix::getListSingle() {
    try {
        ListVector* singlelist = NULL;
        
        if (singletons.size() == 0) { }
        else {
            singlelist = new ListVector();
            
            for (int i = 0; i < singletons.size(); i++) { singlelist->push_back(singletons[i]); }
        }
        
        return singlelist;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiMatrix", "getListSingle");
        exit(1);
    }
}
/***********************************************************************/
long int OptiMatrix::print(ostream& out) {
    try {
        long int count = 0;
        for (int i = 0; i < closeness.size(); i++) {
            for(set<int>::iterator it = closeness[i].begin(); it != closeness[i].end(); it++){
                out << *it << '\t';
                count++;
            }
            out << endl;
        }
        out << endl;
        return count;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiMatrix", "getName");
        exit(1);
    }
}
/***********************************************************************/
//maps unique name to index in distance matrix
//used by sensspec to get translate the list file name to the index name for closeness shirt
map<string, int> OptiMatrix::getNameIndexMap() {
    try {
        map<string, int> nameIndexes;
        for (int i = 0; i < nameMap.size(); i++) {
            vector<string> thisBinsSeqs; util.splitAtComma(nameMap[i], thisBinsSeqs);
            if (i < closeness.size()) {  nameIndexes[thisBinsSeqs[0]] = i;  }
        }
        
        return nameIndexes;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiMatrix", "getNameIndexMap");
        exit(1);
    }
}
/***********************************************************************/
vector<int> OptiMatrix::getNumSeqs(vector<vector<string> > & binNames, vector<vector<int> > & fixedBins) {
    try {
        fixedBins.clear();
        vector<int> movableSeqs;
        
        map<string, int> nameIndexes;
        set<string> unique;
        for (int i = 0; i < nameMap.size(); i++) { //vector of string representing the name file.
            vector<string> thisBinsSeqs; util.splitAtComma(nameMap[i], thisBinsSeqs); //split redundant names
            if (i < closeness.size()) {  nameIndexes[thisBinsSeqs[0]] = i;  } //this is a sequence with distances in the matrix
            if (thisBinsSeqs.size() == 1) { //you are unique
                unique.insert(thisBinsSeqs[0]);
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
        
        for (map<string, int>::iterator it = nameIndexes.begin(); it != nameIndexes.end(); it++) { movableSeqs.push_back(it->second); }
        
        set<string> immovableNames;
        for (int i = 0; i < binNames.size(); i++) { for (int j = 0; j < binNames[i].size(); j++) { immovableNames.insert(binNames[i][j]); }  }
        getNumSingletons(immovableNames);
        
        return movableSeqs;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiMatrix", "getNumSeqs");
        exit(1);
    }
}
/***********************************************************************/
int OptiMatrix::getNumSingletons(set<string> immovables) {
    try {
        numMovableSingletons = 0;
        for (int i = 0; i < singletons.size(); i++) {
            vector<string> theseSeqs; util.splitAtComma(singletons[i], theseSeqs);
            for (int j = 0; j < theseSeqs.size(); j++) {
                if (immovables.count(theseSeqs[j]) == 0) { //you are not immovable
                    numMovableSingletons++;
                }
            }
        }
        
        return numMovableSingletons;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiMatrix", "getNumSingletons");
        exit(1);
    }
}
/***********************************************************************/
string OptiMatrix::getName(int index) {
    try {
        //return toString(index);
        if (index > closeness.size()) { m->mothurOut("[ERROR]: index is not valid.\n"); m->setControl_pressed(true); return ""; }
        string name = nameMap[index];
        return name;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiMatrix", "getName");
        exit(1);
    }
}
/***********************************************************************/
string OptiMatrix::getOverlapName(int index) {
    try {
        //return toString(index);
        if (index > blastOverlap.size()) { m->mothurOut("[ERROR]: index is not valid.\n"); m->setControl_pressed(true); return ""; }
        string name = overlapNameMap[index];
        return name;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiMatrix", "getOverlapName");
        exit(1);
    }
}
/***********************************************************************/
bool OptiMatrix::isClose(int i, int toFind){
    try {
        if (i < 0) { return false; } //if we are running cluster fit and you are a reference sequence with no dists in the matrix
        
        bool found = false;
        if (closeness[i].count(toFind) != 0) { found = true; }
        return found;
        
    }
    catch(exception& e) {
        m->errorOut(e, "OptiMatrix", "isClose");
        exit(1);
    }
}
/***********************************************************************/
long long OptiMatrix::getNumDists(){
    try {
        long long foundDists = 0;
        for (int i = 0; i < closeness.size(); i++) { foundDists += closeness[i].size(); }
        return foundDists;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiMatrix", "isClose");
        exit(1);
    }
}
/***********************************************************************/

int OptiMatrix::readPhylip(){
    try {
        nameMap.clear();
        float distance;
        int square, nseqs;
        string name;
        map<int, int> singletonIndexSwap;
        
        ifstream fileHandle;
        string numTest;
        
        Utils util; util.openInputFile(distFile, fileHandle);
        fileHandle >> numTest >> name;
        nameMap.push_back(name);
        singletonIndexSwap[0] = 0;
        
        if (!util.isContainingOnlyDigits(numTest)) { m->mothurOut("[ERROR]: expected a number and got " + numTest + ", quitting."); m->mothurOutEndLine(); exit(1); }
        else { convert(numTest, nseqs); }
        
        //square test
        char d;
        while((d=fileHandle.get()) != EOF){
            if(isalnum(d)){ square = 1; fileHandle.putback(d); for(int i=0;i<nseqs;i++){ fileHandle >> distance;  } break; }
            if(d == '\n'){ square = 0; break; }
        }
        
        vector<bool> singleton; singleton.resize(nseqs, true);
        ///////////////////// Read to eliminate singletons ///////////////////////
        if(square == 0){
            
            for(int i=1;i<nseqs;i++){
                if (m->getControl_pressed()) {  fileHandle.close();  return 0; }
                
                fileHandle >> name; nameMap.push_back(name); singletonIndexSwap[i] = i;
                
                for(int j=0;j<i;j++){
                    
                    fileHandle >> distance;
                    
                    if (distance == -1) { distance = 1000000; } else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                    
                    if(distance < cutoff){
                        singleton[i] = false;
                        singleton[j] = false;
                    }
                }
            }
        }else{
            for(int i=1;i<nseqs;i++){
                if (m->getControl_pressed()) {  fileHandle.close();   return 0; }
                
                fileHandle >> name; nameMap.push_back(name); singletonIndexSwap[i] = i;
                
                for(int j=0;j<nseqs;j++){
                    fileHandle >> distance;
                    
                    if (distance == -1) { distance = 1000000; } else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                    
                    if(distance < cutoff && j < i){
                        singleton[i] = false;
                        singleton[j] = false;
                    }
                }
            }
        }
        fileHandle.close();
        //////////////////////////////////////////////////////////////////////////
       
        int nonSingletonCount = 0;
        for (int i = 0; i < singleton.size(); i++) {
            if (!singleton[i]) { //if you are not a singleton
                singletonIndexSwap[i] = nonSingletonCount;
                nonSingletonCount++;
            }else { singletons.push_back(nameMap[i]); }
        }
        singleton.clear();

        closeness.resize(nonSingletonCount);
        
        map<string, string> names;
        if (namefile != "") {
            util.readNames(namefile, names);
            for (int i = 0; i < singletons.size(); i++) {
                singletons[i] = names[singletons[i]];
            }
        }
        
        Progress* reading;
        ifstream in;
        
        util.openInputFile(distFile, in);
        in >> nseqs >> name;
        
        if (namefile != "") { name = names[name]; } //redundant names
        nameMap[singletonIndexSwap[0]] = name;
        
        string line = "";
        if(square == 0){
            
            reading = new Progress("Reading matrix:     ", nseqs * (nseqs - 1) / 2);
            int index = 0;
            
            for(int i=1;i<nseqs;i++){
                
                if (m->getControl_pressed()) {  in.close();  delete reading; return 0; }
                
                in >> name; util.gobble(in);
                
                if (namefile != "") { name = names[name]; } //redundant names
                nameMap[singletonIndexSwap[i]] = name;
                
                for(int j=0;j<i;j++){
                    
                    in >> distance; util.gobble(in);
                    
                    if (distance == -1) { distance = 1000000; } else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                    
                    if(distance < cutoff){
                        int newB = singletonIndexSwap[j];
                        int newA = singletonIndexSwap[i];
                        closeness[newA].insert(newB);
                        closeness[newB].insert(newA);
                    }
                    index++; reading->update(index);
                }
            }
        }else{
            reading = new Progress("Reading matrix:     ", nseqs * nseqs);
            int index = nseqs;
            
            for(int i=0;i<nseqs;i++){ in >> distance;  } util.gobble(in);
            
            for(int i=1;i<nseqs;i++){
                if (m->getControl_pressed()) {  in.close();  delete reading; return 0; }
                
                in >> name; util.gobble(in);
                
                if (namefile != "") { name = names[name]; } //redundant names
                nameMap[singletonIndexSwap[i]] = name;
                
                for(int j=0;j<nseqs;j++){
                    in >> distance; util.gobble(in);

                    if (distance == -1) { distance = 1000000; } else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                    
                    if(distance < cutoff && j < i){
                        int newB = singletonIndexSwap[j];
                        int newA = singletonIndexSwap[i];
                        closeness[newA].insert(newB);
                        closeness[newB].insert(newA);
                    }
                    index++; reading->update(index);
                }
            }
        }
        in.close();
        reading->finish();
        delete reading;

        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "OptiMatrix", "readPhylip");
        exit(1);
    }
}
/***********************************************************************/

int OptiMatrix::readColumn(){
    try {
        Utils util;
        map<string, int> nameAssignment;
        if (namefile != "") { nameAssignment = util.readNames(namefile); }
        else  {  CountTable ct; ct.readTable(countfile, false, true); nameAssignment = ct.getNameMap(); }
        int count = 0;
        for (map<string, int>::iterator it = nameAssignment.begin(); it!= nameAssignment.end(); it++) {
            it->second = count; count++;
            nameMap.push_back(it->first);
        }
        
        string firstName, secondName;
        float distance;
        
        ///////////////////// Read to eliminate singletons ///////////////////////
        ifstream fileHandle;
        util.openInputFile(distFile, fileHandle);
        vector<bool> singleton; singleton.resize(nameAssignment.size(), true);
        map<int, int> singletonIndexSwap;
        while(fileHandle){  //let's assume it's a triangular matrix...
            
            fileHandle >> firstName; util.gobble(fileHandle);
            fileHandle >> secondName; util.gobble(fileHandle);
            fileHandle >> distance;	util.gobble(fileHandle); // get the row and column names and distance
            
            if (m->getDebug()) { cout << firstName << '\t' << secondName << '\t' << distance << endl; }
            
            if (m->getControl_pressed()) {  fileHandle.close();   return 0; }
            
            map<string,int>::iterator itA = nameAssignment.find(firstName);
            map<string,int>::iterator itB = nameAssignment.find(secondName);
            
            if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the name or count file, please correct\n"); exit(1);  }
            if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the name or count file, please correct\n"); exit(1);  }
            
            if (distance == -1) { distance = 1000000; }
            else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
            
            if(distance < cutoff){
                int indexA = (itA->second);
                int indexB = (itB->second);
                singleton[indexA] = false;
                singleton[indexB] = false;
                singletonIndexSwap[indexA] = indexA;
                singletonIndexSwap[indexB] = indexB;
            }
        }
        fileHandle.close();
        //////////////////////////////////////////////////////////////////////////
        
        int nonSingletonCount = 0;
        for (int i = 0; i < singleton.size(); i++) {
            if (!singleton[i]) { //if you are a singleton
                singletonIndexSwap[i] = nonSingletonCount;
                nonSingletonCount++;
            }else { singletons.push_back(nameMap[i]); }
        }
        singleton.clear();
        
        closeness.resize(nonSingletonCount);
        
        map<string, string> names;
        if (namefile != "") {
            util.readNames(namefile, names);
            for (int i = 0; i < singletons.size(); i++) {
                singletons[i] = names[singletons[i]];
            }
        }
        
        ifstream in;
        util.openInputFile(distFile, in);
        
        while(in){  //let's assume it's a triangular matrix...
            
            in >> firstName; util.gobble(in);
            in >> secondName; util.gobble(in);
            in >> distance;	util.gobble(in); // get the row and column names and distance
            
            if (m->getDebug()) { cout << firstName << '\t' << secondName << '\t' << distance << endl; }
            
            if (m->getControl_pressed()) {  in.close();   return 0; }

            map<string,int>::iterator itA = nameAssignment.find(firstName);
            map<string,int>::iterator itB = nameAssignment.find(secondName);
            
            if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the name or count file, please correct\n"); exit(1);  }
            if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the name or count file, please correct\n"); exit(1);  }
            
            if (distance == -1) { distance = 1000000; }
            else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
            
            if(distance < cutoff){
                int indexA = (itA->second);
                int indexB = (itB->second);
                
                int newB = singletonIndexSwap[indexB];
                int newA = singletonIndexSwap[indexA];
                closeness[newA].insert(newB);
                closeness[newB].insert(newA);
                
                if (namefile != "") {
                    firstName = names[firstName];  //redundant names
                    secondName = names[secondName]; //redundant names
                }
                
                nameMap[newA] = firstName;
                nameMap[newB] = secondName;
            }
        }
        in.close();
        nameAssignment.clear();
        
        return 1;
        
    }
    catch(exception& e) {
        m->errorOut(e, "OptiMatrix", "readColumn");
        exit(1);
    }
}
/***********************************************************************/
int OptiMatrix::readBlast(){
    try {
        Utils util;
        map<string, int> nameAssignment;
        if (namefile != "") { nameAssignment = util.readNames(namefile); }
        else if (countfile != "") {  CountTable ct; ct.readTable(countfile, false, true); nameAssignment = ct.getNameMap(); }
        else { readBlastNames(nameAssignment);  }
        int count = 0;
        for (map<string, int>::iterator it = nameAssignment.begin(); it!= nameAssignment.end(); it++) {
            it->second = count; count++;
            nameMap.push_back(it->first);
            overlapNameMap.push_back(it->first);
        }
        
        m->mothurOut("Reading Blast File... "); cout.flush();
        
        string firstName, secondName, eScore, currentRow; currentRow = "";
        string repeatName = "";
        float distance, thisoverlap, refScore;
        float percentId;
        float numBases, mismatch, gap, startQuery, endQuery, startRef, endRef, score, lengthThisSeq;
        map<string, float> thisRowsBlastScores;
        
        ///////////////////// Read to eliminate singletons ///////////////////////
        ifstream fileHandle;
        util.openInputFile(distFile, fileHandle);
        
        map<int, int> singletonIndexSwap;
        map<int, int> blastSingletonIndexSwap;
        vector<bool> singleton; singleton.resize(nameAssignment.size(), true);
        vector<bool> overlapSingleton; overlapSingleton.resize(nameAssignment.size(), true);
        vector< map<string,float> > dists;  dists.resize(nameAssignment.size());
        
        if (!fileHandle.eof()) {
            //read in line from file
            fileHandle >> firstName >> secondName >> percentId >> numBases >> mismatch >> gap >> startQuery >> endQuery >> startRef >> endRef >> eScore >> score;
            util.gobble(fileHandle);
            
            currentRow = firstName;
            lengthThisSeq = numBases;
            repeatName = firstName + secondName;
            
            if (firstName == secondName) {   refScore = score;  }
            else{
                //convert name to number
                map<string,int>::iterator itA = nameAssignment.find(firstName);
                map<string,int>::iterator itB = nameAssignment.find(secondName);
                if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the names file, please correct\n"); exit(1);  }
                if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the names file, please correct\n"); exit(1);  }
                
                thisRowsBlastScores[secondName] = score;
                
                //calc overlap score
                thisoverlap = 1.0 - (percentId * (lengthThisSeq - startQuery) / endRef / 100.0 - penalty);
                
                //if there is a valid overlap, add it
                if ((startRef <= length) && ((endQuery+length) >= lengthThisSeq) && (thisoverlap < cutoff)) {
                    int indexA = (itA->second);
                    int indexB = (itB->second);
                    overlapSingleton[indexA] = false;
                    overlapSingleton[indexB] = false;
                    blastSingletonIndexSwap[indexA] = indexA;
                    blastSingletonIndexSwap[indexB] = indexB;
                }
            }
        }else { m->mothurOut("Error in your blast file, cannot read."); m->mothurOutEndLine(); exit(1); }
        
        
        while(fileHandle){  //let's assume it's a triangular matrix...
            
            if (m->getControl_pressed()) { fileHandle.close(); return 0; }
            
            //read in line from file
            fileHandle >> firstName >> secondName >> percentId >> numBases >> mismatch >> gap >> startQuery >> endQuery >> startRef >> endRef >> eScore >> score;
            util.gobble(fileHandle);
            
            string temp = firstName + secondName; //to check if this file has repeat lines, ie. is this a blast instead of a blscreen file
            
            //if this is a new pairing
            if (temp != repeatName) {
                repeatName = temp;
                
                if (currentRow == firstName) {
                    if (firstName == secondName) {  refScore = score; }
                    else{
                        //convert name to number
                        map<string,int>::iterator itA = nameAssignment.find(firstName);
                        map<string,int>::iterator itB = nameAssignment.find(secondName);
                        if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the names file, please correct\n"); exit(1);  }
                        if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the names file, please correct\n"); exit(1);  }
                        
                        //save score
                        thisRowsBlastScores[secondName] = score;
                        
                        //calc overlap score
                        thisoverlap = 1.0 - (percentId * (lengthThisSeq - startQuery) / endRef / 100.0 - penalty);
                        
                        //if there is a valid overlap, add it
                        if ((startRef <= length) && ((endQuery+length) >= lengthThisSeq) && (thisoverlap < cutoff)) {
                            int indexA = (itA->second);
                            int indexB = (itB->second);
                            overlapSingleton[indexA] = false;
                            overlapSingleton[indexB] = false;
                            blastSingletonIndexSwap[indexA] = indexA;
                            blastSingletonIndexSwap[indexB] = indexB;
                        }
                    } //end else
                }else { //end row
                    //convert blast scores to distance and add cell to sparse matrix if we can
                    map<string, float>::iterator it;
                    map<string, float>::iterator itDist;
                    for(it=thisRowsBlastScores.begin(); it!=thisRowsBlastScores.end(); it++) {
                        distance = 1.0 - (it->second / refScore);
                        
                        //do we already have the distance calculated for b->a
                        map<string,int>::iterator itA = nameAssignment.find(currentRow);
                        map<string,int>::iterator itB = nameAssignment.find(it->first);
                        itDist = dists[itB->second].find(itA->first);
                        
                        //if we have it then compare
                        if (itDist != dists[itB->second].end()) {
                            
                            //if you want the minimum blast score ratio, then pick max distance
                            if(minWanted) {	 distance = max(itDist->second, distance);  }
                            else{	distance = min(itDist->second, distance);  }
                            
                            //is this distance below cutoff
                            if (distance < cutoff) {
                                int indexA = (itA->second);
                                int indexB = (itB->second);
                                singleton[indexA] = false;
                                singleton[indexB] = false;
                                singletonIndexSwap[indexA] = indexA;
                                singletonIndexSwap[indexB] = indexB;
                            }
                            //not going to need this again
                            dists[itB->second].erase(itDist);
                        }else { //save this value until we get the other ratio
                            dists[itA->second][it->first] = distance;
                        }
                    }
                    //clear out last rows info
                    thisRowsBlastScores.clear();
                    
                    currentRow = firstName;
                    lengthThisSeq = numBases;
                    
                    //add this row to thisRowsBlastScores
                    if (firstName == secondName) {   refScore = score;  }
                    else{ //add this row to thisRowsBlastScores
                        
                        //convert name to number
                        map<string,int>::iterator itA = nameAssignment.find(firstName);
                        map<string,int>::iterator itB = nameAssignment.find(secondName);
                        if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the names file, please correct\n"); exit(1);  }
                        if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the names file, please correct\n"); exit(1);  }
                        
                        thisRowsBlastScores[secondName] = score;
                        
                        //calc overlap score
                        thisoverlap = 1.0 - (percentId * (lengthThisSeq - startQuery) / endRef / 100.0 - penalty);
                        
                        //if there is a valid overlap, add it
                        if ((startRef <= length) && ((endQuery+length) >= lengthThisSeq) && (thisoverlap < cutoff)) {
                            int indexA = (itA->second);
                            int indexB = (itB->second);
                            overlapSingleton[indexA] = false;
                            overlapSingleton[indexB] = false;
                            blastSingletonIndexSwap[indexA] = indexA;
                            blastSingletonIndexSwap[indexB] = indexB;
                        }
                    }
                }//end if current row
            }//end if repeat
        }
        fileHandle.close();
        
        //convert blast scores to distance and add cell to sparse matrix if we can
        map<string, float>::iterator it;
        map<string, float>::iterator itDist;
        for(it=thisRowsBlastScores.begin(); it!=thisRowsBlastScores.end(); it++) {
            distance = 1.0 - (it->second / refScore);
            
            //do we already have the distance calculated for b->a
            map<string,int>::iterator itA = nameAssignment.find(currentRow);
            map<string,int>::iterator itB = nameAssignment.find(it->first);
            itDist = dists[itB->second].find(itA->first);
            
            //if we have it then compare
            if (itDist != dists[itB->second].end()) {
                
                //if you want the minimum blast score ratio, then pick max distance
                if(minWanted) {	 distance = max(itDist->second, distance);  }
                else{	distance = min(itDist->second, distance);  }
                
                //is this distance below cutoff
                if (distance < cutoff) {
                    int indexA = (itA->second);
                    int indexB = (itB->second);
                    singleton[indexA] = false;
                    singleton[indexB] = false;
                    singletonIndexSwap[indexA] = indexA;
                    singletonIndexSwap[indexB] = indexB;
                }
                //not going to need this again
                dists[itB->second].erase(itDist);
            }else { //save this value until we get the other ratio
                dists[itA->second][it->first] = distance;
            }
        }
        //clear out info
        thisRowsBlastScores.clear();
        dists.clear();

        //////////////////////////////////////////////////////////////////////////
        int nonSingletonCount = 0;
        for (int i = 0; i < singleton.size(); i++) {
            if (!singleton[i]) { //if you are a singleton
                singletonIndexSwap[i] = nonSingletonCount;
                nonSingletonCount++;
            }else { singletons.push_back(nameMap[i]); }
        }
        singleton.clear();
        
        int overlapNonSingletonCount = 0;
        for (int i = 0; i < overlapSingleton.size(); i++) {
            if (!overlapSingleton[i]) { //if you are a singleton
                blastSingletonIndexSwap[i] = overlapNonSingletonCount;
                overlapNonSingletonCount++;
            }
        }
        overlapSingleton.clear();

        ifstream in;
        util.openInputFile(distFile, in);
        
        dists.resize(nameAssignment.size());
        closeness.resize(nonSingletonCount);
        blastOverlap.resize(overlapNonSingletonCount);
        
        map<string, string> names;
        if (namefile != "") {
            util.readNames(namefile, names);
            for (int i = 0; i < singletons.size(); i++) {
                singletons[i] = names[singletons[i]];
            }
        }
        
        m->mothurOut(" halfway ... "); cout.flush();
        
        if (!in.eof()) {
            //read in line from file
            in >> firstName >> secondName >> percentId >> numBases >> mismatch >> gap >> startQuery >> endQuery >> startRef >> endRef >> eScore >> score;
            util.gobble(fileHandle);
            
            currentRow = firstName;
            lengthThisSeq = numBases;
            repeatName = firstName + secondName;
            
            if (firstName == secondName) {   refScore = score;  }
            else{
                //convert name to number
                map<string,int>::iterator itA = nameAssignment.find(firstName);
                map<string,int>::iterator itB = nameAssignment.find(secondName);
                if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the names file, please correct\n"); exit(1);  }
                if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the names file, please correct\n"); exit(1);  }
                
                thisRowsBlastScores[secondName] = score;
                
                if (namefile != "") {
                    firstName = names[firstName];  //redundant names
                    secondName = names[secondName]; //redundant names
                }
                
                nameMap[singletonIndexSwap[itA->second]] = firstName;
                nameMap[singletonIndexSwap[itB->second]] = secondName;
                
                //calc overlap score
                thisoverlap = 1.0 - (percentId * (lengthThisSeq - startQuery) / endRef / 100.0 - penalty);
                
                //if there is a valid overlap, add it
                if ((startRef <= length) && ((endQuery+length) >= lengthThisSeq) && (thisoverlap < cutoff)) {
                    int indexA = (itA->second);
                    int indexB = (itB->second);
                    
                    int newB = blastSingletonIndexSwap[indexB];
                    int newA = blastSingletonIndexSwap[indexA];
                    blastOverlap[newA].insert(newB);
                    blastOverlap[newB].insert(newA);
                    
                    overlapNameMap[newA] = firstName;
                    overlapNameMap[newB] = secondName;
                }
            }
        }else { m->mothurOut("Error in your blast file, cannot read."); m->mothurOutEndLine(); exit(1); }
        
        
        while(in){  //let's assume it's a triangular matrix...
            
            if (m->getControl_pressed()) { fileHandle.close(); return 0; }
            
            //read in line from file
            in >> firstName >> secondName >> percentId >> numBases >> mismatch >> gap >> startQuery >> endQuery >> startRef >> endRef >> eScore >> score;
            util.gobble(fileHandle);
            
            string temp = firstName + secondName; //to check if this file has repeat lines, ie. is this a blast instead of a blscreen file
            
            //if this is a new pairing
            if (temp != repeatName) {
                repeatName = temp;
                
                if (currentRow == firstName) {
                    if (firstName == secondName) {  refScore = score; }
                    else{
                        //convert name to number
                        map<string,int>::iterator itA = nameAssignment.find(firstName);
                        map<string,int>::iterator itB = nameAssignment.find(secondName);
                        if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the names file, please correct\n"); exit(1);  }
                        if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the names file, please correct\n"); exit(1);  }
                        
                        //save score
                        thisRowsBlastScores[secondName] = score;
                        
                        if (namefile != "") {
                            firstName = names[firstName];  //redundant names
                            secondName = names[secondName]; //redundant names
                        }
                        
                        nameMap[singletonIndexSwap[itA->second]] = firstName;
                        nameMap[singletonIndexSwap[itB->second]] = secondName;
                        
                        //calc overlap score
                        thisoverlap = 1.0 - (percentId * (lengthThisSeq - startQuery) / endRef / 100.0 - penalty);
                        
                        //if there is a valid overlap, add it
                        if ((startRef <= length) && ((endQuery+length) >= lengthThisSeq) && (thisoverlap < cutoff)) {
                            int indexA = (itA->second);
                            int indexB = (itB->second);
                            
                            int newB = blastSingletonIndexSwap[indexB];
                            int newA = blastSingletonIndexSwap[indexA];
                            blastOverlap[newA].insert(newB);
                            blastOverlap[newB].insert(newA);
                            
                            overlapNameMap[newA] = firstName;
                            overlapNameMap[newB] = secondName;                        }
                    } //end else
                }else { //end row
                    //convert blast scores to distance and add cell to sparse matrix if we can
                    map<string, float>::iterator it;
                    map<string, float>::iterator itDist;
                    for(it=thisRowsBlastScores.begin(); it!=thisRowsBlastScores.end(); it++) {
                        distance = 1.0 - (it->second / refScore);
                        
                        //do we already have the distance calculated for b->a
                        map<string,int>::iterator itA = nameAssignment.find(currentRow);
                        map<string,int>::iterator itB = nameAssignment.find(it->first);
                        itDist = dists[itB->second].find(itA->first);
                        
                        //if we have it then compare
                        if (itDist != dists[itB->second].end()) {
                            
                            //if you want the minimum blast score ratio, then pick max distance
                            if(minWanted) {	 distance = max(itDist->second, distance);  }
                            else{	distance = min(itDist->second, distance);  }
                            
                            //is this distance below cutoff
                            if (distance < cutoff) {
                                int indexA = (itA->second);
                                int indexB = (itB->second);
                                
                                int newB = singletonIndexSwap[indexB];
                                int newA = singletonIndexSwap[indexA];
                                closeness[newA].insert(newB);
                                closeness[newB].insert(newA);
                            }
                            //not going to need this again
                            dists[itB->second].erase(itDist);
                        }else { //save this value until we get the other ratio
                            dists[itA->second][it->first] = distance;
                        }
                    }
                    //clear out last rows info
                    thisRowsBlastScores.clear();
                    
                    currentRow = firstName;
                    lengthThisSeq = numBases;
                    
                    //add this row to thisRowsBlastScores
                    if (firstName == secondName) {   refScore = score;  }
                    else{ //add this row to thisRowsBlastScores
                        
                        //convert name to number
                        map<string,int>::iterator itA = nameAssignment.find(firstName);
                        map<string,int>::iterator itB = nameAssignment.find(secondName);
                        if(itA == nameAssignment.end()){  m->mothurOut("AError: Sequence '" + firstName + "' was not found in the names file, please correct\n"); exit(1);  }
                        if(itB == nameAssignment.end()){  m->mothurOut("BError: Sequence '" + secondName + "' was not found in the names file, please correct\n"); exit(1);  }
                        
                        thisRowsBlastScores[secondName] = score;
                        
                        //calc overlap score
                        thisoverlap = 1.0 - (percentId * (lengthThisSeq - startQuery) / endRef / 100.0 - penalty);
                        
                        //if there is a valid overlap, add it
                        if ((startRef <= length) && ((endQuery+length) >= lengthThisSeq) && (thisoverlap < cutoff)) {
                            int indexA = (itA->second);
                            int indexB = (itB->second);
                            
                            int newB = blastSingletonIndexSwap[indexB];
                            int newA = blastSingletonIndexSwap[indexA];
                            blastOverlap[newA].insert(newB);
                            blastOverlap[newB].insert(newA);
                            
                            overlapNameMap[newA] = firstName;
                            overlapNameMap[newB] = secondName;
                        }
                    }
                }//end if current row
            }//end if repeat
        }
        in.close();
        
        //convert blast scores to distance and add cell to sparse matrix if we can
        for(it=thisRowsBlastScores.begin(); it!=thisRowsBlastScores.end(); it++) {
            distance = 1.0 - (it->second / refScore);
            
            //do we already have the distance calculated for b->a
            map<string,int>::iterator itA = nameAssignment.find(currentRow);
            map<string,int>::iterator itB = nameAssignment.find(it->first);
            itDist = dists[itB->second].find(itA->first);
            
            //if we have it then compare
            if (itDist != dists[itB->second].end()) {
                
                //if you want the minimum blast score ratio, then pick max distance
                if(minWanted) {	 distance = max(itDist->second, distance);  }
                else{	distance = min(itDist->second, distance);  }
                
                //is this distance below cutoff
                if (distance < cutoff) {
                    int indexA = (itA->second);
                    int indexB = (itB->second);
                    
                    int newB = singletonIndexSwap[indexB];
                    int newA = singletonIndexSwap[indexA];
                    closeness[newA].insert(newB);
                    closeness[newB].insert(newA);
                }
                //not going to need this again
                dists[itB->second].erase(itDist);
            }else { //save this value until we get the other ratio
                dists[itA->second][it->first] = distance;
            }
        }
        //clear out info
        thisRowsBlastScores.clear();
        dists.clear();
        nameAssignment.clear();
        
        m->mothurOut(" done.\n"); 
        
        return 1;
        
    }
    catch(exception& e) {
        m->errorOut(e, "OptiMatrix", "readBlast");
        exit(1);
    }
}
/*********************************************************************************************/
int OptiMatrix::readBlastNames(map<string, int>& nameAssignment) {
    try {
        m->mothurOut("Reading names... "); cout.flush();
        
        string name, hold, prevName;
        int num = 0;
        
        ifstream in;
        Utils util; util.openInputFile(distFile, in);
        
        //read first line
        in >> prevName;
        
        for (int i = 0; i < 11; i++) {  in >> hold;  }
        util.gobble(in);
        
        //save name in nameMap
        nameAssignment[prevName] = num; num++;
        
        map<string, int>::iterator it;
        while (!in.eof()) {
            if (m->getControl_pressed()) { in.close(); return 0; }
            
            //read line
            in >> name;
            
            for (int i = 0; i < 11; i++) {  in >> hold;  }
            util.gobble(in);
            
            //is this a new name?
            if (name != prevName) {
                prevName = name;
                
                it = nameAssignment.find(name);
                if (it != nameAssignment.end()) { m->mothurOut("[ERROR]: trying to exact names from blast file, and I found dups.  Are you sequence names unique? quitting.\n"); m->setControl_pressed(true); }
                else {
                    nameAssignment[name] = num; num++;
                }
            }
        }
        
        in.close();
        
        if (m->getControl_pressed()) { return 0; }
        
        m->mothurOut(toString(num) + " names read.\n");
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "OptiMatrix", "readBlastNames");
        exit(1);
    }
} 

/***********************************************************************/

