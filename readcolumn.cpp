/*
 *  readcolumn.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/21/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "readcolumn.h"
#include "progress.hpp"

/***********************************************************************/

ReadColumnMatrix::ReadColumnMatrix(string df) : distFile(df){
	
	successOpen = m->openInputFile(distFile, fileHandle);
	sim = false;
	
}
/***********************************************************************/

ReadColumnMatrix::ReadColumnMatrix(string df, bool s) : distFile(df){
	
	successOpen = m->openInputFile(distFile, fileHandle);
	sim = s;
}

/***********************************************************************/

int ReadColumnMatrix::read(NameAssignment* nameMap){
	try {		

		string firstName, secondName;
		float distance;
		int nseqs = nameMap->size();

		list = new ListVector(nameMap->getListVector());
	
		Progress* reading = new Progress("Reading matrix:     ", nseqs * nseqs);

		int lt = 1;
		int refRow = 0;	//we'll keep track of one cell - Cell(refRow,refCol) - and see if it's transpose
		int refCol = 0; //shows up later - Cell(refCol,refRow).  If it does, then its a square matrix

		//need to see if this is a square or a triangular matrix...
	
		while(fileHandle && lt == 1){  //let's assume it's a triangular matrix...

		
			fileHandle >> firstName >> secondName >> distance;	// get the row and column names and distance
			
			if (m->control_pressed) {  fileHandle.close();  delete reading; return 0; }
	
			map<string,int>::iterator itA = nameMap->find(firstName);
			map<string,int>::iterator itB = nameMap->find(secondName);
				
			if(itA == nameMap->end()){
				cerr << "AAError: Sequence '" << firstName << "' was not found in the names file, please correct\n"; exit(1);
			}
			if(itB == nameMap->end()){
				cerr << "ABError: Sequence '" << secondName << "' was not found in the names file, please correct\n"; exit(1);
			}
//if (((itA->second == 8) && (itB->second == 1588)) || ((itA->second == 1588) && (itB->second == 8))) { cout << "found it" << endl; }

			if (distance == -1) { distance = 1000000; }
			else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
			
			if(distance < cutoff && itA != itB){
				if(itA->second > itB->second){
					PCell value(itA->second, itB->second, distance);
			
					if(refRow == refCol){		// in other words, if we haven't loaded refRow and refCol...
						refRow = itA->second;
						refCol = itB->second;
						D->addCell(value);
					}
					else if(refRow == itA->second && refCol == itB->second){
						lt = 0;
					}
					else{
						D->addCell(value);
					}
				}
				else if(itA->second < itB->second){
					PCell value(itB->second, itA->second, distance);
			
					if(refRow == refCol){		// in other words, if we haven't loaded refRow and refCol...
						refRow = itA->second;
						refCol = itB->second;
						D->addCell(value);
					}
					else if(refRow == itB->second && refCol == itA->second){
						lt = 0;
					}
					else{
						D->addCell(value);
					}
				}
				reading->update(itA->second * nseqs);
			}
			m->gobble(fileHandle);
		}

		if(lt == 0){  // oops, it was square
	
			fileHandle.close();  //let's start over
			D->clear();  //let's start over
		   
			m->openInputFile(distFile, fileHandle);  //let's start over

			while(fileHandle){
				fileHandle >> firstName >> secondName >> distance;
				
				if (m->control_pressed) {  fileHandle.close();  delete reading; return 0; }
		
				map<string,int>::iterator itA = nameMap->find(firstName);
				map<string,int>::iterator itB = nameMap->find(secondName);
				
				if(itA == nameMap->end()){
					cerr << "BError: Sequence '" << firstName << "' was not found in the names file, please correct\n";
				}
				if(itB == nameMap->end()){
					cerr << "BError: Sequence '" << secondName << "' was not found in the names file, please correct\n";
				}
				
				if (distance == -1) { distance = 1000000; }
				else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
				
				if(distance < cutoff && itA->second > itB->second){
					PCell value(itA->second, itB->second, distance);
					D->addCell(value);
					reading->update(itA->second * nseqs);
				}
		
				m->gobble(fileHandle);
			}
		}
		
		if (m->control_pressed) {  fileHandle.close();  delete reading; return 0; }
		
		reading->finish();
		fileHandle.close();

		list->setLabel("0");
		
		return 1;

	}
	catch(exception& e) {
		m->errorOut(e, "ReadColumnMatrix", "read");
		exit(1);
	}
}

/***********************************************************************/

ReadColumnMatrix::~ReadColumnMatrix(){
	//delete D;
	//delete list;
}


