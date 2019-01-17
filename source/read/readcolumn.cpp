/*
 *  readcolumn.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/21/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "readcolumn.h"


/***********************************************************************/

ReadColumnMatrix::ReadColumnMatrix(string df) : distFile(df){
	
	successOpen = util.openInputFile(distFile, fileHandle);
	sim = false;
	
}
/***********************************************************************/

ReadColumnMatrix::ReadColumnMatrix(string df, bool s) : distFile(df){
	
	successOpen = util.openInputFile(distFile, fileHandle);
	sim = s;
}

/***********************************************************************/

int ReadColumnMatrix::read(NameAssignment* nameMap){
	try {		

		string firstName, secondName;
		float distance;
		int nseqs = nameMap->size();
        DMatrix->resize(nseqs);
		list = new ListVector(nameMap->getListVector());
	
        int lt = 1;
		int refRow = 0;	//we'll keep track of one cell - Cell(refRow,refCol) - and see if it's transpose
		int refCol = 0; //shows up later - Cell(refCol,refRow).  If it does, then its a square matrix

		//need to see if this is a square or a triangular matrix...
	
		while(fileHandle && lt == 1){  //let's assume it's a triangular matrix...

		
			fileHandle >> firstName; util.gobble(fileHandle);
            fileHandle >> secondName; util.gobble(fileHandle);
            fileHandle >> distance;	// get the row and column names and distance
            
            if (m->getDebug()) { cout << firstName << '\t' << secondName << '\t' << distance << endl; }
			
			if (m->getControl_pressed()) {  fileHandle.close();   return 0; }
	
            map<string,int>::iterator itA = nameMap->find(firstName);
            map<string,int>::iterator itB = nameMap->find(secondName);
            
            if(itA == nameMap->end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the names file, please correct\n"); exit(1);  }
            if(itB == nameMap->end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the names file, please correct\n"); exit(1);  }

			if (util.isEqual(distance, -1)) { distance = 1000000; }
			else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
			
			if(distance <= cutoff && itA != itB){
				if(itA->second > itB->second){
                    PDistCell value(itA->second, distance);
                    
                    
					if(refRow == refCol){		// in other words, if we haven't loaded refRow and refCol...
						refRow = itA->second;
						refCol = itB->second;
						DMatrix->addCell(itB->second, value);
					}
					else if(refRow == itA->second && refCol == itB->second){
						lt = 0;
					}
					else{
						DMatrix->addCell(itB->second, value);
					}
				}
				else if(itA->second < itB->second){
					PDistCell value(itB->second, distance);
			
					if(refRow == refCol){		// in other words, if we haven't loaded refRow and refCol...
						refRow = itA->second;
						refCol = itB->second;
						DMatrix->addCell(itA->second, value);
					}
					else if(refRow == itB->second && refCol == itA->second){
						lt = 0;
					}
					else{
						DMatrix->addCell(itA->second, value);
					}
				}
			}
			util.gobble(fileHandle);
		}

		if(lt == 0){  // oops, it was square
	
			fileHandle.close();  //let's start over
			DMatrix->clear();  //let's start over
		   
			util.openInputFile(distFile, fileHandle);  //let's start over

			while(fileHandle){
				fileHandle >> firstName; util.gobble(fileHandle);
                fileHandle >> secondName; util.gobble(fileHandle);
                fileHandle >> distance;	// get the row and column names and distance
				
				if (m->getControl_pressed()) {  fileHandle.close();   return 0; }
		
				map<string,int>::iterator itA = nameMap->find(firstName);
				map<string,int>::iterator itB = nameMap->find(secondName);
				
				if(itA == nameMap->end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the names file, please correct\n"); exit(1);  }
				if(itB == nameMap->end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the names file, please correct\n"); exit(1);  }
				
				if (util.isEqual(distance, -1)) { distance = 1000000; }
				else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
				
				if(distance <= cutoff && itA->second > itB->second){
                    PDistCell value(itA->second, distance);
					DMatrix->addCell(itB->second, value);
				}
		
				util.gobble(fileHandle);
			}
		}
		
		if (m->getControl_pressed()) {  fileHandle.close();   return 0; }
		
		
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

int ReadColumnMatrix::read(CountTable* countTable){
	try {		
        
		string firstName, secondName;
		float distance;
		int nseqs = countTable->size();
        
        DMatrix->resize(nseqs);
		list = new ListVector(countTable->getListVector());
        
		int lt = 1;
		int refRow = 0;	//we'll keep track of one cell - Cell(refRow,refCol) - and see if it's transpose
		int refCol = 0; //shows up later - Cell(refCol,refRow).  If it does, then its a square matrix
        
		//need to see if this is a square or a triangular matrix...
               
		while(fileHandle && lt == 1){  //let's assume it's a triangular matrix...
            
            
			fileHandle >> firstName; util.gobble(fileHandle);
            fileHandle >> secondName; util.gobble(fileHandle);
            fileHandle >> distance;	// get the row and column names and distance
            
			if (m->getControl_pressed()) {  fileHandle.close();   return 0; }
            
			int itA = countTable->get(firstName);
			int itB = countTable->get(secondName);
            
            if (m->getControl_pressed()) { exit(1); }
            
			if (util.isEqual(distance, -1)) { distance = 1000000; }
			else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
			
			if(distance <= cutoff && itA != itB){
				if(itA > itB){
                    PDistCell value(itA, distance);
                    
                    
					if(refRow == refCol){		// in other words, if we haven't loaded refRow and refCol...
						refRow = itA;
						refCol = itB;
						DMatrix->addCell(itB, value);
					}
					else if(refRow == itA && refCol == itB){
						lt = 0;
					}
					else{
						DMatrix->addCell(itB, value);
					}
				}
				else if(itA < itB){
					PDistCell value(itB, distance);
                    
					if(refRow == refCol){		// in other words, if we haven't loaded refRow and refCol...
						refRow = itA;
						refCol = itB;
						DMatrix->addCell(itA, value);
					}
					else if(refRow == itB && refCol == itA){
						lt = 0;
					}
					else{
						DMatrix->addCell(itA, value);
					}
				}
			}
			util.gobble(fileHandle);
		}
        
		if(lt == 0){  // oops, it was square
            
			fileHandle.close();  //let's start over
			DMatrix->clear();  //let's start over
            
			util.openInputFile(distFile, fileHandle);  //let's start over
            
			while(fileHandle){
				fileHandle >> firstName; util.gobble(fileHandle);
                fileHandle >> secondName; util.gobble(fileHandle);
                fileHandle >> distance;	// get the row and column names and distance
				
				if (m->getControl_pressed()) {  fileHandle.close();   return 0; }
                
				int itA = countTable->get(firstName);
                int itB = countTable->get(secondName);
                
                
                if (m->getControl_pressed()) { exit(1); }
				
				if (util.isEqual(distance, -1)) { distance = 1000000; }
				else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
				
				if(distance <= cutoff && itA > itB){
                    PDistCell value(itA, distance);
					DMatrix->addCell(itB, value);
				}
                
				util.gobble(fileHandle);
			}
		}
		
		if (m->getControl_pressed()) {  fileHandle.close();   return 0; }
		
		
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
ReadColumnMatrix::~ReadColumnMatrix(){}
/***********************************************************************/

