
#include "mothur.h"
#include "cluster.hpp"


/***********************************************************************/

SingleLinkage::SingleLinkage(RAbundVector* rav, ListVector* lv, SparseMatrix* dm, float c, string s) :
Cluster(rav, lv, dm, c, s)
{}


/***********************************************************************/
//This function returns the tag of the method.
string SingleLinkage::getTag() {
	return("nn");
}

/***********************************************************************/
//This function clusters based on the single linkage method.
void  SingleLinkage::update(){
	try {
		getRowColCells();	
	
		vector<bool> deleted(nRowCells, false);
		int rowInd;
		int search;
		bool changed;

		// The vector has to be traversed in reverse order to preserve the index
		// for faster removal in removeCell()
		for (int i=nRowCells-1;i>=0;i--) {
			if ((rowCells[i]->row == smallRow) && (rowCells[i]->column == smallCol)) {
				rowInd = i;   // The index of the smallest distance cell in rowCells
			} else {
				if (rowCells[i]->row == smallRow) {
					search = rowCells[i]->column;
				} else {
					search = rowCells[i]->row;
				}
		
				for (int j=0;j<nColCells;j++) {
					if (!((colCells[j]->row == smallRow) && (colCells[j]->column == smallCol))) {
						if (colCells[j]->row == search || colCells[j]->column == search) {
							changed = updateDistance(colCells[j], rowCells[i]);
							// If the cell's distance changed and it had the same distance as 
							// the smallest distance, invalidate the mins vector in SparseMatrix
							if (changed) {
								if (colCells[j]->vectorMap != NULL) {
									*(colCells[j]->vectorMap) = NULL;
									colCells[j]->vectorMap = NULL;
								}
							}
							removeCell(rowCells[i], i , -1);
							deleted[i] = true;
							break;
						}
					}
				}
				if (!deleted[i]) {
					// Assign the cell to the new cluster 
					// remove the old cell from seqVec and add the cell
					// with the new row and column assignment again
					removeCell(rowCells[i], i , -1, false);
					if (search < smallCol){
						rowCells[i]->row = smallCol;
						rowCells[i]->column = search;
					} else {
						rowCells[i]->row = search;
						rowCells[i]->column = smallCol;
					}
					seqVec[rowCells[i]->row].push_back(rowCells[i]);
					seqVec[rowCells[i]->column].push_back(rowCells[i]);
				}
			}	
		}
		clusterBins();
		clusterNames();
		// remove also the cell with the smallest distance
		removeCell(rowCells[rowInd], -1 , -1);
	}
	catch(exception& e) {
		errorOut(e, "SingleLinkage", "update");
		exit(1);
	}
}


/***********************************************************************/
//This function updates the distance based on the nearest neighbor method.
bool SingleLinkage::updateDistance(MatData& colCell, MatData& rowCell) {
	try {
		bool changed = false;
		if (colCell->dist > rowCell->dist) {
			colCell->dist = rowCell->dist;
		}
		return(changed);
	}
	catch(exception& e) {
		errorOut(e, "SingleLinkage", "updateDistance");
		exit(1);
	}
}
/***********************************************************************/
