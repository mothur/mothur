//
//  optimatrix.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/20/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#include "optimatrix.h"
#include "progress.hpp"
#include "nameassignment.hpp"

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
            singlelist = new ListVector(singletons.size());
            
            for (int i = 0; i < singletons.size(); i++) {
                string otu = getName(singletons[i]);
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
string OptiMatrix::getListSingleFile() {
    try {
        string filename = "none";
        
        if (singletons.size() == 0) { return "none"; }
        else {
            string filename = distFile + ".extra.temp";
            ofstream out;
            m->openOutputFile(filename, out);
            if (format == "count") { out << out << "Representative_Sequence\ttotal" << endl;  }
            
            for (int i = 0; i < singletons.size(); i++) {
                string otu = getName(singletons[i]);
                if (format == "name") {  out << singletons[i] << '\t' << otu << endl;  }
                else  {  out << singletons[i] << '\t' << singletons[i] << endl;   }
            }
        }
        
        return filename;
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
        
        ifstream fileHandle;
        string numTest;
        
        m->openInputFile(distFile, fileHandle);
        fileHandle >> numTest >> name;
        
        if (!m->isContainingOnlyDigits(numTest)) { m->mothurOut("[ERROR]: expected a number and got " + numTest + ", quitting."); m->mothurOutEndLine(); exit(1); }
        else { convert(numTest, nseqs); }
        
        //map shorten name to real name - space saver
        nameMap[count] = name;
        //list = new ListVector();
        //list->push_back(toString(count)); //list without singletons above cutoff
        
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
        
        vector< map<int, string> > temp; temp.resize(nseqs);
        map<int, string> tempNameMap;

        if(square == 0){
            
            reading = new Progress("Reading matrix:     ", nseqs * (nseqs - 1) / 2);
            
            int index = 0;
            
            for(int i=1;i<nseqs;i++){
                if (m->control_pressed) {  fileHandle.close();  delete reading; return 0; }
                
                fileHandle >> name;
                tempNameMap[i] = name;
                
                //list->push_back(toString(i));
                
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
                int newIndex = temp.size()+count;
                singletons.push_back(newIndex);
                count++;
                
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
        NameAssignment* nameAssignment;
        if (namefile != "") { nameAssignment = new NameAssignment(namefile);    }
        else                { nameAssignment = new NameAssignment(countfile);   }
        nameAssignment->readMap();
        
        string firstName, secondName;
        float distance;
        map<string, int> indexMap;
        //list = new ListVector();
        
        ifstream fileHandle;
        m->openInputFile(distFile, fileHandle);
        
        vector< map<int, string> > temp; temp.resize(nameAssignment->size());
        map<int, string> tempNameMap;
        
        while(fileHandle){  //let's assume it's a triangular matrix...
            
            fileHandle >> firstName; m->gobble(fileHandle);
            fileHandle >> secondName; m->gobble(fileHandle);
            fileHandle >> distance;	// get the row and column names and distance
            
            if (m->debug) { cout << firstName << '\t' << secondName << '\t' << distance << endl; }
            
            if (m->control_pressed) {  fileHandle.close();   return 0; }
            
            map<string,int>::iterator itA = nameAssignment->find(firstName);
            map<string,int>::iterator itB = nameAssignment->find(secondName);
            
            if(itA == nameAssignment->end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the names file, please correct\n"); exit(1);  }
            if(itB == nameAssignment->end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the names file, please correct\n"); exit(1);  }
            
            if (distance == -1) { distance = 1000000; }
            else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
            
            int indexA = (itA->second)-1;
            int indexB = (itB->second)-1;
            
            if(distance < cutoff){
                temp[indexA][indexB] = secondName;
                temp[indexB][indexA] = firstName;
                tempNameMap[indexA] = firstName;
                tempNameMap[indexB] = secondName;
            }
            m->gobble(fileHandle);
        }
        delete nameAssignment;
        
        map<string, string> names;
        if (namefile != "") { m->readNames(namefile, names); }
        
        map<int, int> closenessIndexMap;
        int count = temp.size();
        for (int i = 0; i < temp.size(); i++) {
            //add to singletons +1 singletons count
            if (temp[i].size() == 0) {
                int newIndex = count;
                singletons.push_back(newIndex);
                count++;
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
