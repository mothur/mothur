/*
 *  fullmatrix.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/6/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "fullmatrix.h"

/**************************************************************************/
//This constructor reads a distance matrix file and stores the data in the matrix.
FullMatrix::FullMatrix(ifstream& filehandle) {
	try{
		globaldata = GlobalData::getInstance();
		groupmap = globaldata->gGroupmap;
		
		string name, group;
		filehandle >> numSeqs >> name;
		
		//make the matrix filled with zeros
		matrix.resize(numSeqs); 
		for(int i = 0; i < numSeqs; i++) {
			matrix[i].resize(numSeqs, 0);
		}
		
		group = groupmap->getGroup(name);
		if(group == "not found") {	cout << "Error: Sequence '" << name << "' was not found in the group file, please correct." << endl; exit(1); }
		index[0].groupname = group; 
		index[0].seqName = name;
		
		//determine if matrix is square or lower triangle
		//if it is square read the distances for the first sequence
		char d;
		while((d=filehandle.get()) != EOF){
			
			//is d a number meaning its square
			if(isalnum(d)){ 
				square = true;
				filehandle.putback(d);
				for(int i=0;i<numSeqs;i++){
					filehandle >> matrix[0][i];
				}
				break;
			}
			
			//is d a line return meaning its lower triangle
			if(d == '\n'){
				square = false;
				break;
			}
		}
		
		//read rest of matrix
		if (square == true) { readSquareMatrix(filehandle); }
		else { readLTMatrix(filehandle); }
		
		
		
	printMatrix(cout);
		//sort sequences so they are gathered in groups for processing
		sortGroups(0, numSeqs-1);
		cout << "after sort" << endl;
	printMatrix(cout);
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the FullMatrix class Function FullMatrix. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the FullMatrix class function FullMatrix. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/**************************************************************************/
void FullMatrix::readSquareMatrix(ifstream& filehandle) {
	try {
	
		Progress* reading;
		reading = new Progress("Reading matrix:    ", numSeqs * numSeqs);
		
		int count = 0;
		float distance;
		string group, name;
		
		for(int i=1;i<numSeqs;i++){
			filehandle >> name;		
			
			group = groupmap->getGroup(name);
			index[i].groupname = group;
			index[i].seqName = name;
			
			if(group == "not found") {	cout << "Error: Sequence '" << name << "' was not found in the group file, please correct." << endl; exit(1); }
				
			for(int j=0;j<numSeqs;j++){
				filehandle >> distance;
					
				matrix[i][j] = distance;
				count++;
				reading->update(count);
			}
		}
		reading->finish();
		delete reading;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the FullMatrix class Function readSquareMatrix. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the FullMatrix class function readSquareMatrix. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

} 
/**************************************************************************/
void FullMatrix::readLTMatrix(ifstream& filehandle) {
	try {
		Progress* reading;
		reading = new Progress("Reading matrix:    ", numSeqs * (numSeqs - 1) / 2);
		
		int count = 0;
		float distance;
		string group, name;
		
		for(int i=1;i<numSeqs;i++){
			filehandle >> name;		
						
			group = groupmap->getGroup(name);
			index[i].groupname = group;
			index[i].seqName = name;
	
			if(group == "not found") {	cout << "Error: Sequence '" << name << "' was not found in the group file, please correct." << endl;  exit(1); }
				
			for(int j=0;j<i;j++){
				filehandle >> distance;
					
				matrix[i][j] = distance;  matrix[j][i] = distance;
				count++;
				reading->update(count);
			}
		}
		reading->finish();
		delete reading;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the FullMatrix class Function readLTMatrix. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the FullMatrix class function readLTMatrix. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

/**************************************************************************/
void FullMatrix::sortGroups(int low, int high){
	try{
	
		int i = low;
		int j = high;
		int y = 0;
		string name;
		
		/* compare value */
		//what group does this row belong to
		string z = index[(low + high) / 2].groupname;

		/* partition */
		do {
			/* find member above ... */
			while(index[i].groupname < z) i++;

			/* find element below ... */
			while(index[j].groupname > z) j--;
			
			if(i <= j) {
				/* swap rows*/
				for (int h = 0; h < numSeqs; h++) {
					y = matrix[i][h];
					matrix[i][h] = matrix[j][h]; 
					matrix[j][h] = y;
				}
				
				/* swap columns*/
				for (int b = 0; b < numSeqs; b++) {
					y = matrix[b][i];
					matrix[b][i] = matrix[b][j]; 
					matrix[b][j] = y;
				}
				
				//swap map elements
				z = index[i].groupname;
				index[i].groupname = index[j].groupname;
				index[j].groupname = z;
				
				name = index[i].seqName;
				index[i].seqName = index[j].seqName;
				index[j].seqName = name;

				
				i++; 
				j--;
cout << "swapping rows " << i << " " << j << endl;
printMatrix(cout); cout << endl;
			}
		} while(i <= j);

		/* recurse */
		if(low < j) 
		sortGroups(low, j);

		if(i < high) 
		sortGroups(i, high); 

	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the FullMatrix class Function sortGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the FullMatrix class function sortGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

/**************************************************************************/	
int FullMatrix::getNumSeqs(){ return numSeqs; }
/**************************************************************************/
//print out matrix
void FullMatrix::printMatrix(ostream& out) {
	try{
		for (int i = 0; i < numSeqs; i++) {
			out << "row " << i << " group = " << index[i].groupname << " name = " << index[i].seqName << endl;
			for (int j = 0; j < numSeqs; j++) {
				out << matrix[i][j] << " ";
			}
			out << endl;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the FullMatrix class Function printMatrix. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the FullMatrix class function printMatrix. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
/**************************************************************************/

