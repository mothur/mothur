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

/***********************************************************************

OptiMatrix::OptiMatrix(string d, string df, double c, bool s) : distFile(d), distFormat(df), cutoff(c), sim(s), OptiData(c){
    countfile = ""; namefile = "";
    
    if (distFormat == "phylip")         { readPhylip();     }
    else if (distFormat == "column")    { readColumn();     }
}
/***********************************************************************/
OptiMatrix::OptiMatrix(string d, string nc, string f, string df, double c, bool s) : distFile(d), distFormat(df), format(f), sim(s), OptiData(c) {
    
    if (format == "name") { namefile = nc; countfile = ""; }
    else if (format == "count") { countfile = nc; namefile = ""; }
    else { countfile = ""; namefile = ""; }
    
    if (distFormat == "phylip")         { readPhylip();     }
    else if (distFormat == "column")    { readColumn();     }
    
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
            
            if (distance == -1) { distance = 1000000; }
            else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
            
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
            
            if (distance == -1) { distance = 1000000; }
            else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
            
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


