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
		
		//sort sequences so they are gathered in groups for processing
		sortGroups(0, numSeqs-1);
			
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
		
		string group, name;
		
		for(int i=1;i<numSeqs;i++){
			filehandle >> name;		
			
			group = groupmap->getGroup(name);
			index[i].groupname = group;
			index[i].seqName = name;
			
			if(group == "not found") {	cout << "Error: Sequence '" << name << "' was not found in the group file, please correct." << endl; exit(1); }
				
			for(int j=0;j<numSeqs;j++){
				filehandle >> matrix[i][j];
				
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
		float y = 0;
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
void FullMatrix::setBounds(){
	try{
		numGroups = globaldata->gGroupmap->namesOfGroups.size();
		
		//sort globaldata->gGroupmap.namesOfGroups so that it will match the matrix
		sort(globaldata->gGroupmap->namesOfGroups.begin(), globaldata->gGroupmap->namesOfGroups.end());
		
		//one for each comparision
		//minsForRows.resize(numGroups*numGroups);
		
		/*************************************************/
		//find where in matrix each group starts and stops
		/*************************************************/
		bounds.resize(numGroups);
		
		bounds[0] = 0;
		bounds[numGroups] = numSeqs;

		//for each group find bounds of subgroup/comparison
		for (int i = 1; i < numGroups; i++) {
			getBounds(bounds[i], globaldata->gGroupmap->namesOfGroups[i-1]);
		}
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the FullMatrix class Function getMinsForRowsVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the FullMatrix class function getMinsForRowsVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
/**************************************************************************/
vector<float> FullMatrix::getMins(int x) {
	try{	
		//clear out old data
		minsForRows.clear();
		
		/************************************************************/
		//fill the minsForRows vector for the box the user wants
		/************************************************************/
		int count = 0;
		int lowBoundx = bounds[0]; //where first group starts
		int lowBoundy = bounds[0]; 
		int highBoundx = bounds[1]; //where second group starts
		int highBoundy = bounds[1]; 
		
		int countx = 1;  //index in bound
		int county = 1;	//index in bound
		
		//find the bounds for the box the user wants
		for (int i = 0; i < (numGroups * numGroups); i++) {
		
			//are you at the box?
			if (count == x) { break; }
			else { count++; }
			
			//move to next box
			if (county < numGroups) {
				county++;
				highBoundy = bounds[county];
				lowBoundy = bounds[county-1];
			}else{ //you are moving to a new row of "boxes"
				county = 1;
				countx++;
				highBoundx = bounds[countx];
				lowBoundx = bounds[countx-1];
				highBoundy = bounds[county];
				lowBoundy = bounds[county-1];
			}
		}
				
		//each row in the box
		for (int x = lowBoundx; x < highBoundx; x++) {
			float min4Row = 100000.0;
			//each entry in that row
			for (int y = lowBoundy; y < highBoundy; y++) {
				//if you are not on the diagonal and you are less than previous minimum
				if ((x != y) && (matrix[x][y] < min4Row)) {
					min4Row = matrix[x][y];
				}
			}
			//save minimum value for that row in minsForRows vector of vectors
			minsForRows.push_back(min4Row);
		}
			
		return minsForRows;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the FullMatrix class Function getMins. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the FullMatrix class function getMins. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/**************************************************************************/
void FullMatrix::getBounds(int& higher, string group) {
	try{
		bool gotLower = false;
		
		//for each group find bounds of subgroup/comparison
		for (it = index.begin(); it != index.end(); it++) {
			if (it->second.groupname == group) {
				gotLower = true; 
			}else if ((gotLower == true) && (it->second.groupname != group)) {  higher = it->first; break; }
		}
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the FullMatrix class Function getBounds. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the FullMatrix class function getBounds. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

/**************************************************************************/
//print out matrix
void FullMatrix::printMinsForRows(ostream& out) {
	try{
		for (int j = 0; j < minsForRows.size(); j++) {
			out << minsForRows[j] << " ";
		}
		out << endl;

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
//shuffles the sequences in the 2 groups passed in.
void FullMatrix::shuffle(string groupA, string groupB){
	try{
		vector<int> rows2Swap;
		vector<int> shuffled;
		float y = 0;
		string name = "";
		
			
		/********************************/
		//save rows you want to randomize
		/********************************/
		//go through the matrix map to find the rows from groups you want to randomize
		for (it = index.begin(); it != index.end(); it++) {
			//is this row from group A or B?
			if ((it->second.groupname == groupA) || (it->second.groupname == groupB)) {
				rows2Swap.push_back(it->first);
				shuffled.push_back(it->first);
			}
		}
		
		//randomize rows to shuffle in shuffled
		random_shuffle(shuffled.begin(), shuffled.end());
		
		/***************************************/
		//swap rows and columns to randomize box
		/***************************************/
		for (int i = 0; i < shuffled.size(); i++) {

			//record the swaps you are making so you can undo them in restore function
			restoreIndex[i].a = shuffled[i];
			restoreIndex[i].b = rows2Swap[i];
			
			/* swap rows*/
			for (int h = 0; h < numSeqs; h++) {
				y = matrix[shuffled[i]][h];
				matrix[shuffled[i]][h] = matrix[rows2Swap[i]][h]; 
				matrix[rows2Swap[i]][h] = y;
			}
				
			/* swap columns */
			for (int b = 0; b < numSeqs; b++) {
				y = matrix[b][shuffled[i]];
				matrix[b][shuffled[i]] = matrix[b][rows2Swap[i]]; 
				matrix[b][rows2Swap[i]] = y;
			}
				
			//swap map elements
			name = index[shuffled[i]].seqName;
			index[shuffled[i]].seqName = index[rows2Swap[i]].seqName;
			index[rows2Swap[i]].seqName = name;

		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the FullMatrix class Function shuffle. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the FullMatrix class function shuffle. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
} 
/**************************************************************************/
//unshuffles the matrix.
void FullMatrix::restore(){
	try{
		float y = 0;
		string name = "";

		//reverse iterate through swaps and undo them to restore original matrix and index map.
		for(it2 = restoreIndex.rbegin(); it2 != restoreIndex.rend(); it2++) {
			/* swap rows */

			for (int h = 0; h < numSeqs; h++) {
				y = matrix[it2->second.a][h];
				matrix[it2->second.a][h] = matrix[it2->second.b][h]; 
				matrix[it2->second.b][h] = y;
			}
			
			/* swap columns */
			for (int b = 0; b < numSeqs; b++) {
				y = matrix[b][it2->second.a];
				matrix[b][it2->second.a] = matrix[b][it2->second.b]; 
				matrix[b][it2->second.b] = y;
			}
			
				
			//swap map elements
			name = index[it2->second.a].seqName;
			index[it2->second.a].seqName = index[it2->second.b].seqName;
			index[it2->second.b].seqName = name;

		}

		//clear restore for next shuffle
		restoreIndex.clear();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the FullMatrix class Function restore. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the FullMatrix class function restore. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}  
/**************************************************************************/
void FullMatrix::getDist(vector<float>& distances) {
	try{
		map<float, float> dist;	 //holds the distances for the integral form
		map<float, float>::iterator it;

		/************************************************************/
		//fill the minsForRows vectors for each group the user wants
		/************************************************************/
		int lowBoundx = bounds[0]; //where first group starts
		int lowBoundy = bounds[0]; 
		int highBoundx = bounds[1]; //where second group starts
		int highBoundy = bounds[1]; 
		
		int countx = 1;  //index in bound
		int county = 1;	//index in bound
		
		//go through each "box" in the matrix
		for (int i = 0; i < (numGroups * numGroups); i++) {
			//each row in the box
			for (int x = lowBoundx; x < highBoundx; x++) {
				float min4Row = 100000.0;
				//each entry in that row
				for (int y = lowBoundy; y < highBoundy; y++) {
					//if you are not on the diagonal and you are less than previous minimum
					if ((x != y) && (matrix[x][y] < min4Row)){
						min4Row = matrix[x][y];
					}
				}
				//save minimum value 
				dist[min4Row] = min4Row;
			}
			
			//****** reset bounds to process next "box" ********
			//if you still have more "boxes" in that row
			if (county < numGroups) {
				county++;
				highBoundy = bounds[county];
				lowBoundy = bounds[county-1];
			}else{ //you are moving to a new row of "boxes"
				county = 1;
				countx++;
				highBoundx = bounds[countx];
				lowBoundx = bounds[countx-1];
				highBoundy = bounds[county];
				lowBoundy = bounds[county-1];
			}
		}

		//store distances in users vector
		for (it = dist.begin(); it != dist.end(); it++) {
			distances.push_back(it->first);
		}
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the FullMatrix class Function restore. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the FullMatrix class function restore. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

