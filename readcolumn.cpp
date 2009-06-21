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
	
	successOpen = openInputFile(distFile, fileHandle);
	
}

/***********************************************************************/

void ReadColumnMatrix::read(NameAssignment* nameMap){
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
	
			map<string,int>::iterator itA = nameMap->find(firstName);
			map<string,int>::iterator itB = nameMap->find(secondName);
			
			if(itA == nameMap->end()){
				cerr << "AAError: Sequence '" << firstName << "' was not found in the names file, please correct\n";
			}
			if(itB == nameMap->end()){
				cerr << "ABError: Sequence '" << secondName << "' was not found in the names file, please correct\n";
			}

			if (distance == -1) { distance = 1000000; }
			
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
			gobble(fileHandle);
		}

		if(lt == 0){  // oops, it was square
			fileHandle.close();  //let's start over
			D->clear();  //let's start over
		   
			openInputFile(distFile, fileHandle);  //let's start over

			while(fileHandle){
				fileHandle >> firstName >> secondName >> distance;
		
				map<string,int>::iterator itA = nameMap->find(firstName);
				map<string,int>::iterator itB = nameMap->find(secondName);
				
				if(itA == nameMap->end()){
					cerr << "BError: Sequence '" << firstName << "' was not found in the names file, please correct\n";
				}
				if(itB == nameMap->end()){
					cerr << "BError: Sequence '" << secondName << "' was not found in the names file, please correct\n";
				}
				
				if (distance == -1) { distance = 1000000; }
				
				if(distance < cutoff && itA->second > itB->second){
					PCell value(itA->second, itB->second, distance);
					D->addCell(value);
					reading->update(itA->second * nseqs);
				}
		
				gobble(fileHandle);
			}
		}

		reading->finish();
		fileHandle.close();

		list->setLabel("0");

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadColumnMatrix class Function read. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadColumnMatrix class function read. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

/***********************************************************************/

ReadColumnMatrix::~ReadColumnMatrix(){
	//delete D;
	//delete list;
}


