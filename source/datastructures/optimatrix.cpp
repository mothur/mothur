//
//  optimatrix.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/20/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#include "optimatrix.h"
#include "counttable.h"

/***********************************************************************/
OptiMatrix::OptiMatrix(vector< set<long long> > close, vector<string> name, vector<string> single, double c) : OptiData(c) {
    closeness = close;
    nameMap = name;
    singletons = single;
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
        long long nseqs;
        bool square = false;
        string name;
        map<long long, long long> singletonIndexSwap;
        
        ifstream fileHandle;
        string numTest;
        
        Utils util; util.openInputFile(distFile, fileHandle);
        fileHandle >> numTest >> name;
        nameMap.push_back(name);
        singletonIndexSwap[0] = 0;
        
        if (!util.isContainingOnlyDigits(numTest)) { m->mothurOut("[ERROR]: expected a number and got " + numTest + ", quitting.\n");  exit(1); }
        else { convert(numTest, nseqs); }
        
        //square test
        char d;
        while((d=fileHandle.get()) != EOF){
            if(isalnum(d)){ square = true; fileHandle.putback(d); for(int i=0;i<nseqs;i++){ fileHandle >> distance;  } break; }
            if(d == '\n'){ square = false; break; }
        }
        
        vector<bool> singleton; singleton.resize(nseqs, true);
        ///////////////////// Read to eliminate singletons ///////////////////////
        if(!square){
            
            for(long long i=1;i<nseqs;i++){
                if (m->getControl_pressed()) {  fileHandle.close();  return 0; }
                
                fileHandle >> name; nameMap.push_back(name); singletonIndexSwap[i] = i;
                
                for(long long j=0;j<i;j++){
                    
                    fileHandle >> distance;
                    
                    if (util.isEqual(distance,-1)) { distance = 1000000; } else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                    
                    if(distance <= cutoff){
                        singleton[i] = false;
                        singleton[j] = false;
                    }
                }
            }
        }else{
            for(long long i=1;i<nseqs;i++){
                if (m->getControl_pressed()) {  fileHandle.close();   return 0; }
                
                fileHandle >> name; nameMap.push_back(name); singletonIndexSwap[i] = i;
                
                for(long long j=0;j<nseqs;j++){
                    fileHandle >> distance;
                    
                    if (util.isEqual(distance,-1)) { distance = 1000000; } else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                    
                    if(distance <= cutoff && j < i){
                        singleton[i] = false;
                        singleton[j] = false;
                    }
                }
            }
        }
        fileHandle.close();
        //////////////////////////////////////////////////////////////////////////
       
        long long nonSingletonCount = 0;
        for (long long i = 0; i < singleton.size(); i++) {
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
            for (long long i = 0; i < singletons.size(); i++) {
                singletons[i] = names[singletons[i]];
            }
        }
        
        ifstream in; util.openInputFile(distFile, in);
        
        in >> nseqs >> name;
        
        if (namefile != "") { name = names[name]; } //redundant names
        nameMap[singletonIndexSwap[0]] = name;
        
        string line = "";
        if(!square){
            
            int index = 0;
            
            for(long long i=1;i<nseqs;i++){
                
                if (m->getControl_pressed()) {  in.close();   return 0; }
                
                in >> name; gobble(in);
                
                if (namefile != "") { name = names[name]; } //redundant names
                nameMap[singletonIndexSwap[i]] = name;
                
                for(long long j=0;j<i;j++){
                    
                    in >> distance; gobble(in);
                    
                    if (util.isEqual(distance,-1)) { distance = 1000000; } else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                    
                    if(distance <= cutoff){
                        long long newB = singletonIndexSwap[j];
                        long long newA = singletonIndexSwap[i];
                        closeness[newA].insert(newB);
                        closeness[newB].insert(newA);
                    }
                    index++; 
                }
            }
        }else{
            long long index = nseqs;
            
            for(long long i=0;i<nseqs;i++){ in >> distance;  } gobble(in);
            
            for(long long i=1;i<nseqs;i++){
                if (m->getControl_pressed()) {  in.close();   return 0; }
                
                in >> name; gobble(in);
                
                if (namefile != "") { name = names[name]; } //redundant names
                nameMap[singletonIndexSwap[i]] = name;
                
                for(long long j=0;j<nseqs;j++){
                    in >> distance; gobble(in);

                    if (util.isEqual(distance,-1)) { distance = 1000000; } else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                    
                    if(distance <= cutoff && j < i){
                        long long newB = singletonIndexSwap[j];
                        long long newA = singletonIndexSwap[i];
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
        m->errorOut(e, "OptiMatrix", "readPhylip");
        exit(1);
    }
}
/***********************************************************************/

int OptiMatrix::readColumn(){
    try {
        Utils util;
        map<string, long long> nameAssignment;
        if (namefile != "") { util.readNames(namefile, nameAssignment); }
        else  {
            CountTable ct; ct.readTable(countfile, false, true);
            map<string, int> temp = ct.getNameMap();
            for (map<string, int>::iterator it = temp.begin(); it!= temp.end(); it++) {  nameAssignment[it->first] = it->second; }
        }
        long long count = 0;
        for (map<string, long long>::iterator it = nameAssignment.begin(); it!= nameAssignment.end(); it++) {
            it->second = count; count++;
            nameMap.push_back(it->first);
        }
        
        string firstName, secondName;
        float distance;
        
        ///////////////////// Read to eliminate singletons ///////////////////////
        ifstream fileHandle;
        util.openInputFile(distFile, fileHandle);
        vector<bool> singleton; singleton.resize(nameAssignment.size(), true);
        map<long long, long long> singletonIndexSwap;
        while(fileHandle){  //let's assume it's a triangular matrix...
            
            fileHandle >> firstName; gobble(fileHandle);
            fileHandle >> secondName; gobble(fileHandle);
            fileHandle >> distance;	gobble(fileHandle); // get the row and column names and distance
            
            if (m->getDebug()) { cout << firstName << '\t' << secondName << '\t' << distance << endl; }
            
            if (m->getControl_pressed()) {  fileHandle.close();   return 0; }
            
            if (util.isEqual(distance,-1)) { distance = 1000000; }
            else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
            
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
        //////////////////////////////////////////////////////////////////////////
        
        long long nonSingletonCount = 0;
        for (long long i = 0; i < singleton.size(); i++) {
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
            for (long long i = 0; i < singletons.size(); i++) {
                singletons[i] = names[singletons[i]];
            }
        }
        
        ifstream in; util.openInputFile(distFile, in);
        
        while(in){  //let's assume it's a triangular matrix...
            
            in >> firstName; gobble(in);
            in >> secondName; gobble(in);
            in >> distance;	gobble(in); // get the row and column names and distance
            
            if (m->getDebug()) { cout << firstName << '\t' << secondName << '\t' << distance << endl; }
            
            if (m->getControl_pressed()) {  in.close();   return 0; }
            
            if (util.isEqual(distance,-1)) { distance = 1000000; }
            else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
            
            if(distance <= cutoff){
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


