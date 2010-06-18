/*
 *  cluster.cpp
 *  
 *
 *  Created by Pat Schloss on 8/14/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

#include "cluster.hpp"
#include "rabundvector.hpp"
#include "listvector.hpp"
#include "sparsematrix.hpp"

/***********************************************************************/

Cluster::Cluster(RAbundVector* rav, ListVector* lv, SparseMatrix* dm, float c, string f) :
rabund(rav), list(lv), dMatrix(dm), method(f)
{
	try {
/*
	cout << "sizeof(MatData): " << sizeof(MatData) << endl;
	cout << "sizeof(PCell*): " << sizeof(PCell*) << endl;

	int nCells = dMatrix->getNNodes();
	time_t start = time(NULL);

	MatVec matvec = MatVec(nCells); 
	int i = 0;
	for (MatData currentCell = dMatrix->begin(); currentCell != dMatrix->end(); currentCell++) {
		matvec[i++] = currentCell;
	}
	for (i= matvec.size();i>0;i--) {
		dMatrix->rmCell(matvec[i-1]);
	}
	MatData it = dMatrix->begin(); 
	while (it != dMatrix->end()) { 
		it = dMatrix->rmCell(it);
	}
	cout << "Time to remove " << nCells << " cells: " << time(NULL) - start << " seconds" << endl;
    exit(0);
	MatData it = dMatrix->begin();
	cout << it->row << "/" << it->column << "/" << it->dist << endl;
	dMatrix->rmCell(dMatrix->begin());
	cout << it->row << "/" << it->column << "/" << it->dist << endl;
	exit(0);
*/

	// Create a data structure to quickly access the PCell information
	// for a certain sequence. It consists of a vector of lists, where 
	// a list contains pointers (iterators) to the all distances related
	// to a certain sequence. The Vector is accessed via the index of a 
	// sequence in the distance matrix.
//ofstream outtemp;
//string temp = "temp";
//openOutputFile(temp, outtemp);	
//cout << lv->size() << endl;
	seqVec = vector<MatVec>(lv->size());
	for (MatData currentCell = dMatrix->begin(); currentCell != dMatrix->end(); currentCell++) {
//outtemp << currentCell->row << '\t' << currentCell->column  << '\t' << currentCell->dist << endl;
		seqVec[currentCell->row].push_back(currentCell);
		seqVec[currentCell->column].push_back(currentCell);
	}
//outtemp.close();
	mapWanted = false;  //set to true by mgcluster to speed up overlap merge
	
	//save so you can modify as it changes in average neighbor
	cutoff = c;
	m = MothurOut::getInstance();
	
	}
	catch(exception& e) {
		m->errorOut(e, "Cluster", "Cluster");
		exit(1);
	}
}

/***********************************************************************/

void Cluster::getRowColCells() {
	try {
		PCell* smallCell = dMatrix->getSmallestCell();	//find the smallest cell - this routine should probably not be in the SpMat class
	
		smallRow = smallCell->row;		// get its row
		smallCol = smallCell->column;	// get its column
		smallDist = smallCell->dist;	// get the smallest distance
	//cout << "small row = " << smallRow << "small col = " << smallCol << "small dist = " << smallDist << endl;
	 
		rowCells = seqVec[smallRow];    // all distances related to the row index
		colCells = seqVec[smallCol];    // all distances related to the column index
		nRowCells = rowCells.size();
		nColCells = colCells.size();
//cout << "num rows = " << nRowCells << "num col = " << nColCells << endl;
		
		//for (int i = 0; i < nColCells; i++) { cout << colCells[i]->row << '\t' << colCells[i]->column << endl;  }
		//for (int i = 0; i < nRowCells; i++) { cout << rowCells[i]->row << '\t' << rowCells[i]->column << endl;  }
	}
	catch(exception& e) {
		m->errorOut(e, "Cluster", "getRowColCells");
		exit(1);
	}

}
/***********************************************************************/
// Remove the specified cell from the seqVec and from the sparse
// matrix
void Cluster::removeCell(const MatData& cell, int vrow, int vcol, bool rmMatrix){
	try {
	
		ull drow = cell->row;
			ull dcol = cell->column;
			if (((vrow >=0) && (drow != smallRow)) ||
				((vcol >=0) && (dcol != smallCol))) {
				ull dtemp = drow;
				drow = dcol;
				dcol = dtemp;
			}

			ull crow;
			ull ccol;
			int nCells;
			if (vrow < 0) {
				nCells = seqVec[drow].size();
				for (vrow=0; vrow<nCells;vrow++) {
					crow = seqVec[drow][vrow]->row;
					ccol = seqVec[drow][vrow]->column;
					if (((crow == drow) && (ccol == dcol)) ||
						((ccol == drow) && (crow == dcol))) {
						break;
					}
				}
			}

			seqVec[drow].erase(seqVec[drow].begin()+vrow);
			if (vcol < 0) {
				nCells = seqVec[dcol].size();
				for (vcol=0; vcol<nCells;vcol++) {
					crow = seqVec[dcol][vcol]->row;
					ccol = seqVec[dcol][vcol]->column;
					if (((crow == drow) && (ccol == dcol)) ||
						((ccol == drow) && (crow == dcol))) {
						break;
					}
				}
			}
		
			seqVec[dcol].erase(seqVec[dcol].begin()+vcol);
		
			if (rmMatrix) {
			//cout << " removing = " << cell->row << '\t' << cell->column  << '\t' << cell->dist << endl;
				dMatrix->rmCell(cell);
		//	cout << "done" << endl;
			}
		
	}
	catch(exception& e) {
		m->errorOut(e, "Cluster", "removeCell");
		exit(1);
	}
}
/***********************************************************************/

void Cluster::clusterBins(){
	try {
	//	cout << smallCol << '\t' << smallRow << '\t' << smallDist << '\t' << rabund->get(smallRow) << '\t' << rabund->get(smallCol);

		rabund->set(smallCol, rabund->get(smallRow)+rabund->get(smallCol));	
		rabund->set(smallRow, 0);	
		rabund->setLabel(toString(smallDist));

	//	cout << '\t' << rabund->get(smallRow) << '\t' << rabund->get(smallCol) << endl;
	}
	catch(exception& e) {
		m->errorOut(e, "Cluster", "clusterBins");
		exit(1);
	}


}

/***********************************************************************/

void Cluster::clusterNames(){
	try {
	//	cout << smallCol << '\t' << smallRow << '\t' << smallDist << '\t' << list->get(smallRow) << '\t' << list->get(smallCol);
		if (mapWanted) {  updateMap();  }
		
		list->set(smallCol, list->get(smallRow)+','+list->get(smallCol));
		list->set(smallRow, "");	
		list->setLabel(toString(smallDist));
	
	//	cout << '\t' << list->get(smallRow) << '\t' << list->get(smallCol) << endl;
    }
	catch(exception& e) {
		m->errorOut(e, "Cluster", "clusterNames");
		exit(1);
	}

}

/***********************************************************************/
//This function clusters based on the method of the derived class
//At the moment only average and complete linkage are covered, because
//single linkage uses a different approach.
void Cluster::update(double& cutOFF){
	try {
		getRowColCells();	
//cout << "got rowcells" << endl;

		vector<int> foundCol(nColCells, 0);

		int search;
		bool changed;

		// The vector has to be traversed in reverse order to preserve the index
		// for faster removal in removeCell()
		for (int i=nRowCells-1;i>=0;i--) {
			if (!((rowCells[i]->row == smallRow) && (rowCells[i]->column == smallCol))) {
				if (rowCells[i]->row == smallRow) {
					search = rowCells[i]->column;
				} else {
					search = rowCells[i]->row;
				}
				
				bool merged = false;
				for (int j=0;j<nColCells;j++) {
					if (!((colCells[j]->row == smallRow) && (colCells[j]->column == smallCol))) { //if you are not hte smallest distance
						if (colCells[j]->row == search || colCells[j]->column == search) {
							foundCol[j] = 1;
							merged = true;
							changed = updateDistance(colCells[j], rowCells[i]);
							// If the cell's distance changed and it had the same distance as 
							// the smallest distance, invalidate the mins vector in SparseMatrix
							if (changed) {
								if (colCells[j]->vectorMap != NULL) {
									*(colCells[j]->vectorMap) = NULL;
									colCells[j]->vectorMap = NULL;
								}
							}
							break;
						}
					}		
				}
				//if not merged it you need it for warning 
				if ((!merged) && (method == "average" || method == "weighted")) {  
					//m->mothurOut("Warning: trying to merge cell " + toString(rowCells[i]->row+1) + " " + toString(rowCells[i]->column+1) + " distance " + toString(rowCells[i]->dist) + " with value above cutoff. Results may vary from using cutoff at cluster command instead of read.dist."); m->mothurOutEndLine(); 
					if (cutOFF > rowCells[i]->dist) {  
						cutOFF = rowCells[i]->dist;  
						//m->mothurOut("changing cutoff to " + toString(cutOFF));  m->mothurOutEndLine(); 
					}

				}
				removeCell(rowCells[i], i , -1);  
				
			}
		}
		clusterBins();
		clusterNames();

		// Special handling for singlelinkage case, not sure whether this
		// could be avoided
		for (int i=nColCells-1;i>=0;i--) {
			if (foundCol[i] == 0) {
				if (method == "average" || method == "weighted") {
					if (!((colCells[i]->row == smallRow) && (colCells[i]->column == smallCol))) {
						//m->mothurOut("Warning: merging cell " + toString(colCells[i]->row+1) + " " + toString(colCells[i]->column+1) + " distance " + toString(colCells[i]->dist) + " value above cutoff. Results may vary from using cutoff at cluster command instead of read.dist."); m->mothurOutEndLine();
						if (cutOFF > colCells[i]->dist) {  
							cutOFF = colCells[i]->dist;  
							//m->mothurOut("changing cutoff to " + toString(cutOFF));  m->mothurOutEndLine(); 
						}
					}
				}
				removeCell(colCells[i], -1, i);
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "Cluster", "update");
		exit(1);
	}
}
/***********************************************************************/
void Cluster::setMapWanted(bool f)  {  
	try {
		mapWanted = f;
		
		//initialize map
		for (int i = 0; i < list->getNumBins(); i++) {
			
			//parse bin 
			string names = list->get(i);
			while (names.find_first_of(',') != -1) { 
				//get name from bin
				string name = names.substr(0,names.find_first_of(','));
				//save name and bin number
				seq2Bin[name] = i;
				names = names.substr(names.find_first_of(',')+1, names.length());
			}
			
			//get last name
			seq2Bin[names] = i;
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "Cluster", "setMapWanted");
		exit(1);
	}
}
/***********************************************************************/
void Cluster::updateMap() {
try {
		//update location of seqs in smallRow since they move to smallCol now
		string names = list->get(smallRow);
		while (names.find_first_of(',') != -1) { 
			//get name from bin
			string name = names.substr(0,names.find_first_of(','));
			//save name and bin number
			seq2Bin[name] = smallCol;
			names = names.substr(names.find_first_of(',')+1, names.length());
		}
			
		//get last name
		seq2Bin[names] = smallCol;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Cluster", "updateMap");
		exit(1);
	}
}
/***********************************************************************/



