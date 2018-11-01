/*
 *  readblast.cpp
 *  Mothur
 *
 *  Created by westcott on 12/10/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "readblast.h"


//********************************************************************************************************************
//sorts lowest to highest
inline bool compareOverlap(seqDist left, seqDist right){
	return (left.dist < right.dist);	
} 
/*********************************************************************************************/
ReadBlast::ReadBlast(string file, float c, float p, int l, bool ms) : blastfile(file), cutoff(c), penalty(p), length(l), minWanted(ms) {
	try {
		m = MothurOut::getInstance();
		matrix = NULL;
	}
	catch(exception& e) {
		m->errorOut(e, "ReadBlast", "ReadBlast");
		exit(1);
	}
} 
/*********************************************************************************************/
//assumptions about the blast file: 
//1. if duplicate lines occur the first line is always best and is chosen
//2. blast scores are grouped together, ie. a a .... score, a b .... score, a c ....score...
int ReadBlast::read(NameAssignment* nameMap) {
	try {
	
		//if the user has not given a names file read names from blastfile
		if (nameMap->size() == 0) { readNames(nameMap);  }
		int nseqs = nameMap->size();
		
		if (m->getControl_pressed()) { return 0; }

		ifstream fileHandle;
		util.openInputFile(blastfile, fileHandle);
		
		string firstName, secondName, eScore, currentRow;
		string repeatName = "";
		int count = 1;
		float distance, thisoverlap, refScore;
		float percentId; 
		float numBases, mismatch, gap, startQuery, endQuery, startRef, endRef, score, lengthThisSeq;
		
		ofstream outDist;
		ofstream outOverlap;
		
		//create objects needed for read
        matrix = new SparseDistanceMatrix();
        matrix->resize(nseqs);
		
		
		if (m->getControl_pressed()) { fileHandle.close(); delete matrix; return 0; }
		
		//this is used to quickly find if we already have a distance for this combo
		vector< map<int,float> > dists;  dists.resize(nseqs);  //dists[0][1] = distance from seq0 to seq1
		map<int, float> thisRowsBlastScores;
		
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
				map<string,int>::iterator itA = nameMap->find(firstName);
				map<string,int>::iterator itB = nameMap->find(secondName);
				if(itA == nameMap->end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the names file, please correct\n"); exit(1);  }
				if(itB == nameMap->end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the names file, please correct\n"); exit(1);  }
				
				thisRowsBlastScores[itB->second] = score;
				
				//calc overlap score
				thisoverlap = 1.0 - (percentId * (lengthThisSeq - startQuery) / endRef / 100.0 - penalty);
				
				//if there is a valid overlap, add it
				if ((startRef <= length) && ((endQuery+length) >= lengthThisSeq) && (thisoverlap <= cutoff)) {
                    seqDist overlapValue(itA->second, itB->second, thisoverlap);
                    overlap.push_back(overlapValue);
				}
			}
		}else { m->mothurOut("Error in your blast file, cannot read."); m->mothurOutEndLine(); exit(1); }

       
		//read file
		while(!fileHandle.eof()){  
		
			if (m->getControl_pressed()) { fileHandle.close(); delete matrix;  return 0; }
			
			//read in line from file
			fileHandle >> firstName >> secondName >> percentId >> numBases >> mismatch >> gap >> startQuery >> endQuery >> startRef >> endRef >> eScore >> score;
			
			util.gobble(fileHandle);
			
			string temp = firstName + secondName; //to check if this file has repeat lines, ie. is this a blast instead of a blscreen file
			
			//if this is a new pairing
			if (temp != repeatName) {
				repeatName = temp; 
				
				if (currentRow == firstName) {
						
					if (firstName == secondName) { 
						refScore = score;  
						count++; 
					}else{
						//convert name to number
						map<string,int>::iterator itA = nameMap->find(firstName);
						map<string,int>::iterator itB = nameMap->find(secondName);
						if(itA == nameMap->end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the names file, please correct\n"); exit(1);  }
						if(itB == nameMap->end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the names file, please correct\n"); exit(1);  }
						
						//save score
						thisRowsBlastScores[itB->second] = score;
							
						//calc overlap score
						thisoverlap = 1.0 - (percentId * (lengthThisSeq - startQuery) / endRef / 100.0 - penalty);
						
						//if there is a valid overlap, add it
						if ((startRef <= length) && ((endQuery+length) >= lengthThisSeq) && (thisoverlap <= cutoff)) {
                            seqDist overlapValue(itA->second, itB->second, thisoverlap);
                            overlap.push_back(overlapValue);
						}
					} //end else
				}else { //end row
					//convert blast scores to distance and add cell to sparse matrix if we can
					map<int, float>::iterator it;
					map<int, float>::iterator itDist;
					for(it=thisRowsBlastScores.begin(); it!=thisRowsBlastScores.end(); it++) {  
						distance = 1.0 - (it->second / refScore);
		
						
						//do we already have the distance calculated for b->a
						map<string,int>::iterator itA = nameMap->find(currentRow);
						itDist = dists[it->first].find(itA->second);
						
						//if we have it then compare
						if (itDist != dists[it->first].end()) {
	
							//if you want the minimum blast score ratio, then pick max distance
							if(minWanted) {	 distance = max(itDist->second, distance);  }
							else{	distance = min(itDist->second, distance);  }

							//is this distance below cutoff
							if (distance <= cutoff) {
                                if (itA->second < it->first) {
                                    PDistCell value(it->first, distance);
                                    matrix->addCell(itA->second, value);
                                }else {
                                    PDistCell value(itA->second, distance);
                                    matrix->addCell(it->first, value);
                                }
							}
							//not going to need this again
							dists[it->first].erase(itDist);
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
						map<string,int>::iterator itA = nameMap->find(firstName);
						map<string,int>::iterator itB = nameMap->find(secondName);
						if(itA == nameMap->end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the names file, please correct\n"); exit(1);  }
						if(itB == nameMap->end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the names file, please correct\n"); exit(1);  }
						
						thisRowsBlastScores[itB->second] = score;
						
						//calc overlap score
						thisoverlap = 1.0 - (percentId * (lengthThisSeq - startQuery) / endRef / 100.0 - penalty);
						
						//if there is a valid overlap, add it
						if ((startRef <= length) && ((endQuery+length) >= lengthThisSeq) && (thisoverlap <= cutoff)) {
                            seqDist overlapValue(itA->second, itB->second, thisoverlap);
                            overlap.push_back(overlapValue);
						}
					}
				}//end if current row
			}//end if repeat
		}//end while
		
		//get last rows info stored
		//convert blast scores to distance and add cell to sparse matrix if we can
		map<int, float>::iterator it;
		map<int, float>::iterator itDist;
		for(it=thisRowsBlastScores.begin(); it!=thisRowsBlastScores.end(); it++) {  
			distance = 1.0 - (it->second / refScore);
			
			//do we already have the distance calculated for b->a
			map<string,int>::iterator itA = nameMap->find(currentRow);
			itDist = dists[it->first].find(itA->second);
			
			//if we have it then compare
			if (itDist != dists[it->first].end()) {
				//if you want the minimum blast score ratio, then pick max distance
				if(minWanted) {	 distance = max(itDist->second, distance);  }
				else{	distance = min(itDist->second, distance);  }
				
				//is this distance below cutoff
				if (distance <= cutoff) {
                    if (itA->second < it->first) {
                        PDistCell value(it->first, distance);
                        matrix->addCell(itA->second, value);
                    }else {
                        PDistCell value(itA->second, distance);
                        matrix->addCell(it->first, value);
                    }
				}
				//not going to need this again
				dists[it->first].erase(itDist);
			}else { //save this value until we get the other ratio
				dists[itA->second][it->first] = distance;
			}
		}
		//clear out info
		thisRowsBlastScores.clear();
		dists.clear();
		
		if (m->getControl_pressed()) {  fileHandle.close(); delete matrix;  return 0; }
		
        sort(overlap.begin(), overlap.end(), compareOverlap);
 
		if (m->getControl_pressed()) {  fileHandle.close(); delete matrix;  return 0; }
		
		
		
		fileHandle.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ReadBlast", "read");
		exit(1);
	}
} 
/*********************************************************************************************/
int ReadBlast::readNames(NameAssignment* nameMap) {
	try {
		m->mothurOut("Reading names... "); cout.flush();
		
		string name, hold, prevName;
		int num = 1;
		
		ifstream in;
		util.openInputFile(blastfile, in);
		
		//ofstream outName;
		//util.openOutputFile((blastfile + ".tempOutNames"), outName);
		
		//read first line
		in >> prevName;
	
		for (int i = 0; i < 11; i++) {  in >> hold;  }
		util.gobble(in);
				
		//save name in nameMap
		nameMap->push_back(prevName);
		
		while (!in.eof()) {
			if (m->getControl_pressed()) { in.close(); return 0; }
			
			//read line
			in >> name;
	
			for (int i = 0; i < 11; i++) {  in >> hold;  }
			util.gobble(in);
			
			//is this a new name?
			if (name != prevName) {
				prevName = name;
                
                if (nameMap->get(name) != -1) { m->mothurOut("[ERROR]: trying to exact names from blast file, and I found dups.  Are you sequence names unique? quitting.\n"); m->setControl_pressed(true); }
                else {
                    nameMap->push_back(name);
                }
                //outName << name << '\t' << name << endl;
				num++;
			}
		}
	
		in.close();
		
		//write out names file
		//string outNames = util.getRootName(blastfile) + "names";
		//ofstream out;
		//util.openOutputFile(outNames, out);
		//nameMap->print(out);
		//out.close();
		
		if (m->getControl_pressed()) { return 0; }
		
		m->mothurOut(toString(num) + " names read."); m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ReadBlast", "readNames");
		exit(1);
	}
} 
/*********************************************************************************************/



