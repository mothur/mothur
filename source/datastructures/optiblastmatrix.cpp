//
//  optiblastmatrix.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/10/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "optiblastmatrix.hpp"

/***********************************************************************/
OptiBlastMatrix::OptiBlastMatrix(string d, string nc, string f, bool s, double c, int l, float p, bool min) : OptiData(c), length(l), penalty(p), minWanted(min) {
    
    m = MothurOut::getInstance();
    
    distFile = d; format = f; sim = s; 
    
    if (format == "name") { namefile = nc; countfile = ""; }
    else if (format == "count") { countfile = nc; namefile = ""; }
    else { countfile = ""; namefile = ""; }
    
    readBlast();
}
/***********************************************************************/
string OptiBlastMatrix::getOverlapName(long long index) {
    try {
        if (index > blastOverlap.size()) { m->mothurOut("[ERROR]: index is not valid.\n"); m->setControl_pressed(true); return ""; }
        string name = overlapNameMap[index];
        return name;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiBlastMatrix", "getOverlapName");
        exit(1);
    }
}
/***********************************************************************/
int OptiBlastMatrix::readBlast(){
    try {
        Utils util;
        map<string, long long> nameAssignment;
        if (namefile != "") { util.readNames(namefile, nameAssignment); }
        else if (countfile != "") {
            CountTable ct; ct.readTable(countfile, false, true);
            map<string, int> temp = ct.getNameMap();
            for (map<string, int>::iterator it = temp.begin(); it!= temp.end(); it++) {  nameAssignment[it->first] = it->second; }
        }
        else { readBlastNames(nameAssignment);  }
        int count = 0;
        for (map<string, long long>::iterator it = nameAssignment.begin(); it!= nameAssignment.end(); it++) {
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
                thisRowsBlastScores[secondName] = score;
                
                //calc overlap score
                thisoverlap = 1.0 - (percentId * (lengthThisSeq - startQuery) / endRef / 100.0 - penalty);
                
                //if there is a valid overlap, add it
                if ((startRef <= length) && ((endQuery+length) >= lengthThisSeq) && (thisoverlap <= cutoff)) {
                    //convert name to number
                    map<string,long long>::iterator itA = nameAssignment.find(firstName);
                    map<string,long long>::iterator itB = nameAssignment.find(secondName);
                    if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the names file, please correct\n"); exit(1);  }
                    if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the names file, please correct\n"); exit(1);  }

                    int indexA = (itA->second);
                    int indexB = (itB->second);
                    overlapSingleton[indexA] = false;
                    overlapSingleton[indexB] = false;
                    blastSingletonIndexSwap[indexA] = indexA;
                    blastSingletonIndexSwap[indexB] = indexB;
                }
            }
        }else { m->mothurOut("Error in your blast file, cannot read.\n");  exit(1); }
        
        
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
                        //save score
                        thisRowsBlastScores[secondName] = score;
                        
                        //calc overlap score
                        thisoverlap = 1.0 - (percentId * (lengthThisSeq - startQuery) / endRef / 100.0 - penalty);
                        
                        //if there is a valid overlap, add it
                        if ((startRef <= length) && ((endQuery+length) >= lengthThisSeq) && (thisoverlap <= cutoff)) {
                            //convert name to number
                            map<string,long long>::iterator itA = nameAssignment.find(firstName);
                            map<string,long long>::iterator itB = nameAssignment.find(secondName);
                            if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the names file, please correct\n"); exit(1);  }
                            if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the names file, please correct\n"); exit(1);  }
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
                        map<string,long long>::iterator itA = nameAssignment.find(currentRow);
                        map<string,long long>::iterator itB = nameAssignment.find(it->first);
                        itDist = dists[itB->second].find(itA->first);
                        
                        //if we have it then compare
                        if (itDist != dists[itB->second].end()) {
                            
                            //if you want the minimum blast score ratio, then pick max distance
                            if(minWanted) {	 distance = max(itDist->second, distance);  }
                            else{	distance = min(itDist->second, distance);  }
                            
                            //is this distance below cutoff
                            if (distance <= cutoff) {
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
                        thisRowsBlastScores[secondName] = score;
                        
                        //calc overlap score
                        thisoverlap = 1.0 - (percentId * (lengthThisSeq - startQuery) / endRef / 100.0 - penalty);
                        
                        //if there is a valid overlap, add it
                        if ((startRef <= length) && ((endQuery+length) >= lengthThisSeq) && (thisoverlap <= cutoff)) {
                            //convert name to number
                            map<string,long long>::iterator itA = nameAssignment.find(firstName);
                            map<string,long long>::iterator itB = nameAssignment.find(secondName);
                            if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the names file, please correct\n"); exit(1);  }
                            if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the names file, please correct\n"); exit(1);  }
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
            map<string,long long>::iterator itA = nameAssignment.find(currentRow);
            map<string,long long>::iterator itB = nameAssignment.find(it->first);
            itDist = dists[itB->second].find(itA->first);
            
            //if we have it then compare
            if (itDist != dists[itB->second].end()) {
                
                //if you want the minimum blast score ratio, then pick max distance
                if(minWanted) {	 distance = max(itDist->second, distance);  }
                else{	distance = min(itDist->second, distance);  }
                
                //is this distance below cutoff
                if (distance <= cutoff) {
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
        
        ifstream in; util.openInputFile(distFile, in);
        
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
                map<string,long long>::iterator itA = nameAssignment.find(firstName);
                map<string,long long>::iterator itB = nameAssignment.find(secondName);
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
                if ((startRef <= length) && ((endQuery+length) >= lengthThisSeq) && (thisoverlap <= cutoff)) {
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
        }else { m->mothurOut("Error in your blast file, cannot read.\n");  exit(1); }
        
        
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
                        map<string,long long>::iterator itA = nameAssignment.find(firstName);
                        map<string,long long>::iterator itB = nameAssignment.find(secondName);
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
                        if ((startRef <= length) && ((endQuery+length) >= lengthThisSeq) && (thisoverlap <= cutoff)) {
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
                        map<string,long long>::iterator itA = nameAssignment.find(currentRow);
                        map<string,long long>::iterator itB = nameAssignment.find(it->first);
                        itDist = dists[itB->second].find(itA->first);
                        
                        //if we have it then compare
                        if (itDist != dists[itB->second].end()) {
                            
                            //if you want the minimum blast score ratio, then pick max distance
                            if(minWanted) {	 distance = max(itDist->second, distance);  }
                            else{	distance = min(itDist->second, distance);  }
                            
                            //is this distance below cutoff
                            if (distance <= cutoff) {
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
                        map<string,long long>::iterator itA = nameAssignment.find(firstName);
                        map<string,long long>::iterator itB = nameAssignment.find(secondName);
                        if(itA == nameAssignment.end()){  m->mothurOut("AError: Sequence '" + firstName + "' was not found in the names file, please correct\n"); exit(1);  }
                        if(itB == nameAssignment.end()){  m->mothurOut("BError: Sequence '" + secondName + "' was not found in the names file, please correct\n"); exit(1);  }
                        
                        thisRowsBlastScores[secondName] = score;
                        
                        //calc overlap score
                        thisoverlap = 1.0 - (percentId * (lengthThisSeq - startQuery) / endRef / 100.0 - penalty);
                        
                        //if there is a valid overlap, add it
                        if ((startRef <= length) && ((endQuery+length) >= lengthThisSeq) && (thisoverlap <= cutoff)) {
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
            map<string,long long>::iterator itA = nameAssignment.find(currentRow);
            map<string,long long>::iterator itB = nameAssignment.find(it->first);
            itDist = dists[itB->second].find(itA->first);
            
            //if we have it then compare
            if (itDist != dists[itB->second].end()) {
                
                //if you want the minimum blast score ratio, then pick max distance
                if(minWanted) {	 distance = max(itDist->second, distance);  }
                else{	distance = min(itDist->second, distance);  }
                
                //is this distance below cutoff
                if (distance <= cutoff) {
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
        m->errorOut(e, "OptiBlastMatrix", "readBlast");
        exit(1);
    }
}
/*********************************************************************************************/
int OptiBlastMatrix::readBlastNames(map<string, long long>& nameAssignment) {
    try {
        m->mothurOut("Reading names... "); cout.flush();
        
        string name, hold, prevName;
        int num = 0;
        
        ifstream in; Utils util; util.openInputFile(distFile, in);
        
        //read first line
        in >> prevName;
        
        for (int i = 0; i < 11; i++) {  in >> hold;  }
        util.gobble(in);
        
        //save name in nameMap
        nameAssignment[prevName] = num; num++;
        
        map<string, long long>::iterator it;
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
                
                else { nameAssignment[name] = num; num++; }
            }
        }
        
        in.close();
        
        if (m->getControl_pressed()) { return 0; }
        
        m->mothurOut(toString(num) + " names read.\n");
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "OptiBlastMatrix", "readBlastNames");
        exit(1);
    }
}
/***********************************************************************/



