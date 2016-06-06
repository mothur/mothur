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

OptiMatrix::OptiMatrix(string d, double c, bool s) : distFile(d), cutoff(c), sim(s) {
    m = MothurOut::getInstance();
    countfile = ""; namefile = "";
    distFormat = findDistFormat(distFile);
    
    if (distFormat == "phylip") { readPhylip(); }
    else { readColumn();  }
}
/***********************************************************************/
OptiMatrix::OptiMatrix(string d, string nc, string f, double c, bool s) : distFile(d), format(f), cutoff(c), sim(s) {
    m = MothurOut::getInstance();
    
    if (format == "name") { namefile = nc; countfile = ""; }
    else { countfile = nc; namefile = ""; }
    
    distFormat = findDistFormat(distFile);
    
    if (distFormat == "phylip") { readPhylip(); }
    else { readColumn();  }
}
/***********************************************************************/
ListVector* OptiMatrix::getListSingle() {
    try {
        ListVector* singlelist = NULL;
        
        if (singletons.size() == 0) { }
        else {
            singlelist = new ListVector();
            
            for (int i = 0; i < singletons.size(); i++) {
                string otu = singletons[i];
                singlelist->push_back(otu);
            }
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
            for(int j=0; j < closeness[i].size();j++){
                out << closeness[i][j] << '\t';
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
string OptiMatrix::getName(int index) {
    try {
        map<int, string>::iterator it = nameMap.find(index);
        
        string name = "not found";
        if (it == nameMap.end()) { m->mothurOut("[ERROR]: cannot find name for index " + toString(index) + "\n"); m-> control_pressed = true; }
        else { name = it->second; }
        
        return name;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiMatrix", "getName");
        exit(1);
    }
}
/***********************************************************************/
//assumes sorted optimatrix
bool OptiMatrix::isClose(int i, int toFind){
    try {
        // Returns index of toFind in sortedArray, or -1 if not found
        int low = 0;
        int high = closeness[i].size() - 1;
        int mid;
        
        int l = closeness[i][low];
        int h = closeness[i][high];
        
        while (l <= toFind && h >= toFind) {
            mid = (low + high)/2;
            
            int m = closeness[i][mid];
            
            if (m < toFind) {
                l = closeness[i][low = mid + 1];
            } else if (m > toFind) {
                h = closeness[i][high = mid - 1];
            } else {
                return true;
            }
        }
        
        if (closeness[i][low] == toFind) {
            return true;
        }else{
            return false; // Not found
        }
    }
    catch(exception& e) {
        m->errorOut(e, "OptiMatrix", "isClose");
        exit(1);
    }
}
/***********************************************************************/

string OptiMatrix::findDistFormat(string distFile){
    try {
        string fileFormat = "column";
        
        ifstream fileHandle;
        string numTest;
        
        m->openInputFile(distFile, fileHandle);
        fileHandle >> numTest;
        fileHandle.close();
        
        if (m->isContainingOnlyDigits(numTest)) { fileFormat = "phylip"; }
        
        return fileFormat;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiMatrix", "findDistFormat");
        exit(1);
    }
}
/***********************************************************************/

int OptiMatrix::readPhylip(){
    try {
        float distance;
        int square, nseqs;
        string name;
        int count = 0;
        vector< map<int, string> > temp; temp.resize(nseqs);
        map<int, string> tempNameMap;

        
        ifstream fileHandle;
        string numTest;
        
        m->openInputFile(distFile, fileHandle);
        fileHandle >> numTest >> name;
        
        if (!m->isContainingOnlyDigits(numTest)) { m->mothurOut("[ERROR]: expected a number and got " + numTest + ", quitting."); m->mothurOutEndLine(); exit(1); }
        else { convert(numTest, nseqs); }
        
        //map shorten name to real name - space saver
        tempNameMap[count] = name;
        
        //square test
        char d;
        while((d=fileHandle.get()) != EOF){
            
            if(isalnum(d)){
                square = 1;
                fileHandle.putback(d);
                for(int i=0;i<nseqs;i++){
                    fileHandle >> distance;
                }
                break;
            }
            if(d == '\n'){
                square = 0;
                break;
            }
        }
        
        Progress* reading;
       
        if(square == 0){
            
            reading = new Progress("Reading matrix:     ", nseqs * (nseqs - 1) / 2);
            
            int index = 0;
            
            for(int i=1;i<nseqs;i++){
                if (m->control_pressed) {  fileHandle.close();  delete reading; return 0; }
                
                fileHandle >> name;
                tempNameMap[i] = name;
                
                for(int j=0;j<i;j++){
                    
                    if (m->control_pressed) { delete reading; fileHandle.close(); return 0;  }
                    
                    fileHandle >> distance;
                    
                    if (distance == -1) { distance = 1000000; }
                    else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                    
                    if(distance < cutoff){
                        map<int, string> tempA;
                        tempA[j] = name;
                        temp[i] = tempA;
                        temp[j][i] = tempNameMap[j];
                    }
                    index++;
                    reading->update(index);
                }
            }
        }
        else{
            
            reading = new Progress("Reading matrix:     ", nseqs * nseqs);
            
            int index = nseqs;
            
            for(int i=1;i<nseqs;i++){
                fileHandle >> name;
                
                tempNameMap[i] = name;
                
                //list->push_back(toString(i));
                
                for(int j=0;j<nseqs;j++){
                    fileHandle >> distance;
                    
                    if (m->control_pressed) {  fileHandle.close();  delete reading; return 0; }
                    
                    if (distance == -1) { distance = 1000000; }
                    else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                    
                    if(distance < cutoff && j < i){
                        map<int, string> tempA;
                        tempA[j] = name;
                        temp[i] = tempA;
                        temp[j][i] = tempNameMap[j];
                    }
                    index++;
                    reading->update(index);
                }
            }
        }
        
        map<string, string> names;
        if (namefile != "") { m->readNames(namefile, names); }
        
        count = 0;
        map<int, int> closenessIndexMap;
        for (int i = 0; i < temp.size(); i++) {
            //add to singletons +1 singletons count
            if (temp[i].size() == 0) {
                //add new Index singleton to nameMap
                map<int, string>::iterator it = tempNameMap.find(i);
                if (it != tempNameMap.end()) {
                    if (namefile != "") {
                        map<string, string>::iterator it2 = names.find(it->second);
                        if (it2 != names.end()) { //set name = names in namefile
                            singletons.push_back(it2->second);
                        }else{  m->mothurOut("[ERROR]: cannot find " + it->second + " int your name file, please correct.\n"); m->control_pressed = true; }
                    }else{
                        singletons.push_back(it->second);
                    }
                }
            }else {
                int newIndex = closeness.size();
                
                vector<int> thisClose;
                for(map<int, string>:: iterator it = temp[i].begin(); it != temp[i].end(); it++) {
                    thisClose.push_back(it->first);
                }
                closeness.push_back(thisClose);
                closenessIndexMap[i] = closeness.size()-1;
                
                //add new Index singleton to nameMap
                map<int, string>::iterator it = tempNameMap.find(i);
                if (namefile != "") {
                    map<string, string>::iterator it2 = names.find(it->second);
                    if (it2 != names.end()) { //set name = names in namefile
                        nameMap[newIndex] = it2->second;
                    }else{  m->mothurOut("[ERROR]: cannot find " + it->second + " int your name file, please correct.\n"); m->control_pressed = true; }
                }else{
                    nameMap[newIndex] = it->second;
                }
            
                temp[i].clear();
            }
        }
        tempNameMap.clear();
        
        
        for (int i = 0; i < closeness.size(); i++) {
            for (int j = 0; j < closeness[i].size(); j++) {
                closeness[i][j] = closenessIndexMap[closeness[i][j]];
            }
            sort(closeness[i].begin(), closeness[i].end());
        }
        
        if (m->control_pressed) {  fileHandle.close();  delete reading; return 0; }
        
        reading->finish();
        delete reading;
        
        //list->setLabel("0");
        fileHandle.close();
        
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
        map<string, int> nameAssignment;
        if (namefile != "") { nameAssignment = m->readNames(namefile); }
        else  { CountTable ct; ct.readTable(countfile, false, true); nameAssignment = ct.getNameMap(); }
        int count = 0; for (map<string, int>::iterator it = nameAssignment.begin(); it!= nameAssignment.end(); it++) { it->second = count; count++; }
        
        string firstName, secondName;
        float distance;
        map<string, int> indexMap;
        //list = new ListVector();
        
        ifstream fileHandle;
        m->openInputFile(distFile, fileHandle);
        
        vector< map<int, string> > temp; temp.resize(nameAssignment.size());
        vector<string> tempNameMap; tempNameMap.resize(nameAssignment.size(), "");
        
        while(fileHandle){  //let's assume it's a triangular matrix...
            
            fileHandle >> firstName; m->gobble(fileHandle);
            fileHandle >> secondName; m->gobble(fileHandle);
            fileHandle >> distance;	// get the row and column names and distance
            
            if (m->debug) { cout << firstName << '\t' << secondName << '\t' << distance << endl; }
            
            if (m->control_pressed) {  fileHandle.close();   return 0; }
            
            map<string,int>::iterator itA = nameAssignment.find(firstName);
            map<string,int>::iterator itB = nameAssignment.find(secondName);
            
            if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the name or count file, please correct\n"); exit(1);  }
            if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the name or count file, please correct\n"); exit(1);  }
            
            if (distance == -1) { distance = 1000000; }
            else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
            
            int indexA = (itA->second);
            int indexB = (itB->second);
            
            if(distance < cutoff){
                temp[indexA][indexB] = secondName;
                temp[indexB][indexA] = firstName;
                tempNameMap[indexA] = firstName;
                tempNameMap[indexB] = secondName;
            }
            m->gobble(fileHandle);
        }
        
        map<string, string> names;
        if (namefile != "") { m->readNames(namefile, names); }
        
        map<int, int> closenessIndexMap;
        for (int i = 0; i < temp.size(); i++) {
            if (temp[i].size() == 0) {
            }else {
                string newName = tempNameMap[i];
                int newIndex = closeness.size();
                
                vector<int> thisClose;
                for(map<int, string>:: iterator it = temp[i].begin(); it != temp[i].end(); it++) {
                    thisClose.push_back(it->first);
                }
                closeness.push_back(thisClose);
                closenessIndexMap[i] = closeness.size()-1;
                
                //add new Index singleton to nameMap
                if (namefile != "") {
                    map<string, string>::iterator it2 = names.find(newName);
                    if (it2 != names.end()) { //set name = names in namefile
                        nameMap[newIndex] = it2->second;
                    }else{  m->mothurOut("[ERROR]: cannot find " + it2->second + " int your name file, please correct.\n"); m->control_pressed = true; }
                }else{
                    nameMap[newIndex] = newName;
                }
                temp[i].clear();
            }
        }
        tempNameMap.clear();
        
        if (nameAssignment.size() != nameMap.size()) { //you have singletons to remove
            for(map<int, string>::iterator it = nameMap.begin(); it != nameMap.end(); it++) {  //remove non singletons
                nameAssignment.erase(it->second);
            }
            for(map<string, int>::iterator it = nameAssignment.begin(); it != nameAssignment.end(); it++) {
                if (namefile == "") { singletons.push_back(it->first); }
                else {  map<string, string>::iterator it2 = names.find(it->first);
                    if (it2 != names.end()) { //set name = names in namefile
                        singletons.push_back(it2->second);
                    }else{  m->mothurOut("[ERROR]: cannot find " + it->first + " in your name file, please correct.\n"); m->control_pressed = true; }
                }
            }
        }

        
        for (int i = 0; i < closeness.size(); i++) {
            for (int j = 0; j < closeness[i].size(); j++) {
                closeness[i][j] = closenessIndexMap[closeness[i][j]];
            }
            sort(closeness[i].begin(), closeness[i].end());
        }
        
        if (m->control_pressed) {  fileHandle.close();   return 0; }
        
        fileHandle.close();
        
        //list->setLabel("0");
        
        return 1;
        
    }
    catch(exception& e) {
        m->errorOut(e, "OptiMatrix", "readColumn");
        exit(1);
    }
}
/***********************************************************************/
