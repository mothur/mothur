//
//  optimatrix.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/20/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#include "optimatrix.h"
#include "progress.hpp"

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
//assumes sorted optimatrix
bool OptiMatrix::isClose(int i, int j){
    try {
        int low = 0;
        int high = closeness[i].size() - 1;
        int mid = 0;
        
        int l = closeness[i][low];
        int h = closeness[i][high];
        
        while (l <= j && h >= j) {
            mid = low + ((int)(int)(high - low)*(int)(j - l))/((int)(h-l));
            
            int m = closeness[i][mid];
            
            if (m < j) {
                l = closeness[i][low = mid + 1];
            } else if (m > j) {
                h = closeness[i][high = mid - 1];
            } else {
                return mid;
            }
        }
        
        if (closeness[i][low] == j) {
            return true;
        }else{
            return false; // Not found
        }
        
        return false;
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
        list = new ListVector();
        list->push_back(toString(count)); //list without singletons above cutoff
        
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
        
        vector< map<int, string> > temp;
        map<int, string> tempNameMap;

        if(square == 0){
            
            reading = new Progress("Reading matrix:     ", nseqs * (nseqs - 1) / 2);
            
            int index = 0;
            
            for(int i=1;i<nseqs;i++){
                if (m->control_pressed) {  fileHandle.close();  delete reading; return 0; }
                
                fileHandle >> name;
                tempNameMap[i] = name;
                
                list->push_back(toString(i));
                
                for(int j=0;j<i;j++){
                    
                    if (m->control_pressed) { delete reading; fileHandle.close(); return 0;  }
                    
                    fileHandle >> distance;
                    
                    if (distance == -1) { distance = 1000000; }
                    else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                    
                    if(distance < cutoff){
                        map<int, string> tempA; tempA[j] = name;
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
                
                list->push_back(toString(i));
                
                for(int j=0;j<nseqs;j++){
                    fileHandle >> distance;
                    
                    if (m->control_pressed) {  fileHandle.close();  delete reading; return 0; }
                    
                    if (distance == -1) { distance = 1000000; }
                    else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                    
                    if(distance < cutoff && j < i){
                        map<int, string> tempA; tempA[j] = name;
                        temp[i] = tempA;
                        temp[j][i] = tempNameMap[j];
                    }
                    index++;
                    reading->update(index);
                }
            }
        }
        
        count = 0;
        for (int i = 0; i < temp.size(); i++) {
            //add to singletons +1 singletons count
            if (temp[i].size() == 0) {
                int newIndex = temp.size()+count;
                singletons.push_back(newIndex);
                count++;
                
                //add new Index singleton to nameMap
                map<int, string>::iterator it = tempNameMap.find(i);
                nameMap[newIndex] = it->second;
                
                tempNameMap.erase(it);
            }else {
                int newIndex = closeness.size();
                
                vector<int> thisClose;
                for(map<int, string>:: iterator it = temp[i].begin(); it != temp[i].end(); it++) {
                    thisClose.push_back(it->first);
                }
                closeness.push_back(thisClose);
                
                //add new Index singleton to nameMap
                map<int, string>::iterator it = tempNameMap.find(i);
                nameMap[newIndex] = it->second;
                
                tempNameMap.erase(it);
                
                temp[i].clear();
            }
        }

        
        for (int i = 0; i < closeness.size(); i++) {  sort(closeness[i].begin(), closeness[i].end());  }
        
        if (m->control_pressed) {  fileHandle.close();  delete reading; return 0; }
        
        reading->finish();
        delete reading;
        
        list->setLabel("0");
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
        
        string firstName, secondName;
        float distance;
        map<string, int> indexMap;
        list = new ListVector();
        
        ifstream fileHandle;
        m->openInputFile(distFile, fileHandle);
        
        vector< map<int, string> > temp;
        map<int, string> tempNameMap;
        
        while(fileHandle){  //let's assume it's a triangular matrix...
            
            fileHandle >> firstName; m->gobble(fileHandle);
            fileHandle >> secondName; m->gobble(fileHandle);
            fileHandle >> distance;	// get the row and column names and distance
            
            if (m->debug) { cout << firstName << '\t' << secondName << '\t' << distance << endl; }
            
            if (m->control_pressed) {  fileHandle.close();   return 0; }
            
            map<string,int>::iterator itA = indexMap.find(firstName);
            map<string,int>::iterator itB = indexMap.find(secondName);
            
            int indexI, indexJ; bool newA = false;  bool newB = false;
            if(itA == indexMap.end()){
                newA = true;
                indexMap[firstName] = indexMap.size();
                indexI = indexMap.size();
                tempNameMap[indexI] = firstName;
                list->push_back(toString(indexI));
            }else { indexI = itA->second; }
            
            if(itB == indexMap.end()){
                newB = true;
                indexMap[secondName] = indexMap.size();
                indexJ = indexMap.size();
                tempNameMap[indexJ] = secondName;
                list->push_back(toString(indexJ));
            }else { indexJ = itB->second; }
            
            if (distance == -1) { distance = 1000000; }
            else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
            
            if(distance < cutoff){
                if (newA) {
                    map<int, string> tempA; tempA[indexJ] = firstName;
                    temp[indexI] = tempA;
                }else {  temp[indexI][indexJ] = firstName; }
                if (newB) {
                    map<int, string> tempB; tempB[indexI] = secondName;
                    temp[indexJ] = tempB;
                }else {  temp[indexJ][indexI] = secondName; }
            }
            m->gobble(fileHandle);
        }
        
        int count = 0;
        for (int i = 0; i < temp.size(); i++) {
            //add to singletons +1 singletons count
            if (temp[i].size() == 0) {
                int newIndex = temp.size()+count;
                singletons.push_back(newIndex);
                count++;
                
                //add new Index singleton to nameMap
                map<int, string>::iterator it = tempNameMap.find(i);
                nameMap[newIndex] = it->second;
                
                tempNameMap.erase(it);
            }else {
                int newIndex = closeness.size();
                
                vector<int> thisClose;
                for(map<int, string>:: iterator it = temp[i].begin(); it != temp[i].end(); it++) {
                    thisClose.push_back(it->first);
                }
                closeness.push_back(thisClose);
                
                //add new Index singleton to nameMap
                map<int, string>::iterator it = tempNameMap.find(i);
                nameMap[newIndex] = it->second;
                
                tempNameMap.erase(it);
                
                temp[i].clear();
            }
        }
        
        for (int i = 0; i < closeness.size(); i++) {  sort(closeness[i].begin(), closeness[i].end());  }
        
        if (m->control_pressed) {  fileHandle.close();   return 0; }
        
        fileHandle.close();
        
        list->setLabel("0");
        
        return 1;
        
    }
    catch(exception& e) {
        m->errorOut(e, "OptiMatrix", "readColumn");
        exit(1);
    }
}
/***********************************************************************/
