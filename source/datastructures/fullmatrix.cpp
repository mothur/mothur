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
FullMatrix::FullMatrix(ifstream& filehandle, GroupMap* g, bool s) : groupmap(g), sim(s) {
	try{
		m = MothurOut::getInstance();
		
		string name, group;
		
		filehandle >> numSeqs >> name;
	
		//make the matrix filled with zeros
		matrix.resize(numSeqs); 
		for(int i = 0; i < numSeqs; i++) {
			matrix[i].resize(numSeqs, 0.0);
		}
		group = groupmap->getGroup(name);
		if(group == "not found") {	m->mothurOut("Error: Sequence '" + name + "' was not found in the group file, please correct."); m->mothurOutEndLine(); exit(1); }
		index.resize(numSeqs);
		index[0].seqName = name;
		index[0].groupName = group;
		
		//determine if matrix is square or lower triangle
		//if it is square read the distances for the first sequence
		char d;
		bool square;
		while((d=filehandle.get()) != EOF){
			
			//is d a number meaning its square
			if(isalnum(d)){ 
				square = true;
				filehandle.putback(d);
				
				for(int i=0;i<numSeqs;i++){
					filehandle >> matrix[0][i];
					if (sim) {  matrix[0][i] = 1.0 - matrix[0][i];  }
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
		if (square ) {  readSquareMatrix(filehandle); }
		else {  readLTMatrix(filehandle); }
		
		filehandle.close();
		
		if (!m->getControl_pressed()) { sortGroups(0, numSeqs-1); }
				
	}
	catch(exception& e) {
		m->errorOut(e, "FullMatrix", "FullMatrix");
		exit(1);
	}
}
/**************************************************************************/
int FullMatrix::readSquareMatrix(ifstream& filehandle) {
	try {
	
		Progress* reading;
		reading = new Progress("Reading matrix:     ", numSeqs * numSeqs);
		
		int count = 0;
		
		string group, name;
	
		for(int i=1;i<numSeqs;i++){
			filehandle >> name;		
			
			group = groupmap->getGroup(name);
			index[i].seqName = name;
			index[i].groupName = group;
			
			if(group == "not found") {	m->mothurOut("Error: Sequence '" + name + "' was not found in the group file, please correct."); m->mothurOutEndLine(); exit(1); }
				
			for(int j=0;j<numSeqs;j++){
				if (m->getControl_pressed()) { delete reading;  return 0; }
				
				filehandle >> matrix[i][j];
				if (sim) {  matrix[i][j] = 1.0 - matrix[i][j];  }
				
				count++;
				reading->update(count);
			}
		}
		
		if (m->getControl_pressed()) { delete reading;  return 0; }
		
		reading->finish();
		delete reading;
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "FullMatrix", "readSquareMatrix");
		exit(1);
	}
} 
/**************************************************************************/
int FullMatrix::readLTMatrix(ifstream& filehandle) {
	try {
		
		Progress* reading;
		reading = new Progress("Reading matrix:     ", numSeqs * (numSeqs - 1) / 2);
		
		int count = 0;
		float distance;

		string group, name;
	
		for(int i=1;i<numSeqs;i++){
			filehandle >> name;		
					
			group = groupmap->getGroup(name);
			index[i].seqName = name;
			index[i].groupName = group;
	
			if(group == "not found") {	m->mothurOut("Error: Sequence '" + name + "' was not found in the group file, please correct."); m->mothurOutEndLine();  exit(1); }
				
			for(int j=0;j<i;j++){
				if (m->getControl_pressed()) { delete reading;  return 0; }
				
				filehandle >> distance;
				if (sim) {  distance = 1.0 - distance;  }
				
				matrix[i][j] = distance;  matrix[j][i] = distance;
				
				count++;
				reading->update(count);
			}
		}
		
		if (m->getControl_pressed()) { delete reading;  return 0; }
		
		reading->finish();
		delete reading;
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "FullMatrix", "readLTMatrix");
		exit(1);
	}
}

/**************************************************************************/

void FullMatrix::sortGroups(int low, int high){
	try{
		
		if (low < high) {
			int i = low+1;
			int j = high;
			int pivot = (low+high) / 2;
			
			swapRows(low, pivot);  //puts pivot in final spot
			
			/* compare value */
			//what group does this row belong to
			string key = index[low].groupName;
			
			/* partition */
			while(i <= j) {
				/* find member above ... */
				while((i <= high) && (index[i].groupName <= key))	{  i++;  }  
				
				/* find element below ... */
				while((j >= low) && (index[j].groupName > key))		{  j--;  } 
								
				if(i < j) {
					swapRows(i, j);
				}
			} 
			
			swapRows(low, j);
			
			/* recurse */
			sortGroups(low, j-1);
			sortGroups(j+1, high); 
		}
	
	}
	catch(exception& e) {
		m->errorOut(e, "FullMatrix", "sortGroups");
		exit(1);
	}
}

/**************************************************************************/	
void FullMatrix::swapRows(int i, int j) {
	try {
	
		float y;
		string z, name;
		
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
		z = index[i].groupName;
		index[i].groupName = index[j].groupName;
		index[j].groupName = z;
		
		name = index[i].seqName;
		index[i].seqName = index[j].seqName;
		index[j].seqName = name;
		
		
	}
	catch(exception& e) {
		m->errorOut(e, "FullMatrix", "swapRows");
		exit(1);
	}
}
/**************************************************************************/	

float FullMatrix::get(int i, int j){	return matrix[i][j];		}

/**************************************************************************/	

vector<string> FullMatrix::getGroups(){	return groups;		}

/**************************************************************************/	

vector<int> FullMatrix::getSizes(){	return sizes;		}

/**************************************************************************/	

int FullMatrix::getNumGroups(){	return groups.size();		}

/**************************************************************************/	

int FullMatrix::getNumSeqs(){	return numSeqs;		}

/**************************************************************************/

void FullMatrix::printMatrix(ostream& out) {
	try{
		for (int i = 0; i < numSeqs; i++) {
			out << "row " << i << " group = " << index[i].groupName << " name = " << index[i].seqName << endl;
			for (int j = 0; j < numSeqs; j++) {
				out << i << '\t' << j << '\t' << matrix[i][j] << endl;
			}
			out << endl;
		}
		
		for (int i = 0; i < numSeqs; i++) {  out << i << '\t' <<  index[i].seqName << endl;  }
	}
	catch(exception& e) {
		m->errorOut(e, "FullMatrix", "printMatrix");
		exit(1);
	}
}

/**************************************************************************/

