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
		
				if(nameMap->count(firstName)==0){
					cerr << "AError: Sequence '" << firstName << "' was not found in the names file, please correct\n";
				}
				if(nameMap->count(secondName)==0){
					cerr << "AError: Sequence '" << secondName << "' was not found in the names file, please correct\n";
				}
				
				if (distance == -1) { distance = 1000000; }
				
				if(distance < cutoff && nameMap->get(firstName) != nameMap->get(secondName)){
					if(nameMap->get(firstName) > nameMap->get(secondName)){
						PCell value(nameMap->get(firstName), nameMap->get(secondName), distance);
				
						if(refRow == refCol){		// in other words, if we haven't loaded refRow and refCol...
							refRow = nameMap->get(firstName);
							refCol = nameMap->get(secondName);
							D->addCell(value);
						}
						else if(refRow == nameMap->get(firstName) && refCol == nameMap->get(secondName)){
							lt = 0;
						}
						else{
							D->addCell(value);
						}
					}
					else if(nameMap->get(firstName) < nameMap->get(secondName)){
						PCell value(nameMap->get(secondName), nameMap->get(firstName), distance);
				
						if(refRow == refCol){		// in other words, if we haven't loaded refRow and refCol...
							refRow = nameMap->get(firstName);
							refCol = nameMap->get(secondName);
							D->addCell(value);
						}
						else if(refRow == nameMap->get(secondName) && refCol == nameMap->get(firstName)){
							lt = 0;
						}
						else{
							D->addCell(value);
						}
					}
					reading->update(nameMap->get(firstName) * nseqs);
				}
				gobble(fileHandle);
			}

			if(lt == 0){  // oops, it was square
				fileHandle.close();  //let's start over
				D->clear();  //let's start over
			   
				openInputFile(distFile, fileHandle);  //let's start over

				while(fileHandle){
					fileHandle >> firstName >> secondName >> distance;
			
					if(nameMap->count(firstName)==0){
						cerr << "BError: Sequence '" << firstName << "' was not found in the names file, please correct\n";
					}
					if(nameMap->count(secondName)==0){
						cerr << "BError: Sequence '" << secondName << "' was not found in the names file, please correct\n";
					}
					
					if (distance == -1) { distance = 1000000; }
					
					if(distance < cutoff && nameMap->get(firstName) > nameMap->get(secondName)){
						PCell value(nameMap->get(firstName), nameMap->get(secondName), distance);
						D->addCell(value);
						reading->update(nameMap->get(firstName) * nseqs);
					}
			
					gobble(fileHandle);
				}
			}
		//	else if(lt == 0){
		//		while(fileHandle){
		//			fileHandle >> firstName >> secondName >> distance;
		//			
		//			if(nameMap->count(firstName)==0){
		//				cerr << "CError: Sequence '" << firstName << "' was not found in the names file, please correct\n";
		//			}
		//			if(nameMap->count(secondName)==0){
		//				cerr << "CError: Sequence '" << secondName << "' was not found in the names file, please correct\n";
		//			}
		//			if (distance == -1) { distance = 1000000; }
		
		//			if(distance < cutoff && (*nameMap)[firstName].second < (*nameMap)[secondName].second){
		////				cout << (*nameMap)[secondName] << ' ' << (*nameMap)[firstName] << ' ' << distance << endl;
		//				D->addCell(Cell((*nameMap)[secondName].second, (*nameMap)[firstName].second, distance));
		//				reading->update((*nameMap)[secondName].second * nseqs);
		//			}
		//
		//			gobble(fileHandle);
		//		}
		//	}	
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
	delete D;
	delete list;
}


