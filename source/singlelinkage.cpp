
#include "cluster.hpp"


/***********************************************************************/

SingleLinkage::SingleLinkage(RAbundVector* rav, ListVector* lv, SparseDistanceMatrix* dm, float c, string s, float a) :
Cluster(rav, lv, dm, c, s, a)
{}


/***********************************************************************/
//This function returns the tag of the method.
string SingleLinkage::getTag() {
	return("nn");
}

/***********************************************************************
//This function clusters based on the single linkage method.
void SingleLinkage::update(double& cutOFF){
	try {
		smallCol = dMatrix->getSmallestCell(smallRow);
        nColCells = dMatrix->seqVec[smallCol].size();
        nRowCells = dMatrix->seqVec[smallRow].size();	
	
		vector<bool> deleted(nRowCells, false);
		int rowInd;
		int search;
		bool changed;

		// The vector has to be traversed in reverse order to preserve the index
		// for faster removal in removeCell()
		for (int i=nRowCells-1;i>=0;i--) {
                if (dMatrix->seqVec[smallRow][i].index == smallCol) {
                    rowInd = i;   // The index of the smallest distance cell in rowCells
                } else {
                    search = dMatrix->seqVec[smallRow][i].index;
                    
                    for (int j=0;j<nColCells;j++) {
                        if (dMatrix->seqVec[smallCol][j].index != smallRow) { //if you are not the small cell
                            if (dMatrix->seqVec[smallCol][j].index == search) {
                                changed = updateDistance(dMatrix->seqVec[smallCol][j], dMatrix->seqVec[smallRow][i]);
                                dMatrix->updateCellCompliment(smallCol, j);
                                dMatrix->rmCell(smallRow, i);
                                deleted[i] = true;
                                break;
                            }
                        }
                    }
                    if (!deleted[i]) {
                        // Assign the cell to the new cluster 
                        // remove the old cell from seqVec and add the cell
                        // with the new row and column assignment again
                        float distance =  dMatrix->seqVec[smallRow][i].dist;
                        dMatrix->rmCell(smallRow, i);
                        if (search < smallCol){
                            PDistCell value(smallCol, distance);
                            dMatrix->addCell(search, value);
                        } else {
                            PDistCell value(search, distance);
                            dMatrix->addCell(smallCol, value);
                        }
                        sort(dMatrix->seqVec[smallCol].begin(), dMatrix->seqVec[smallCol].end(), compareIndexes);
                        sort(dMatrix->seqVec[search].begin(), dMatrix->seqVec[search].end(), compareIndexes); 
                    }
                }
		}
		clusterBins();
		clusterNames();
		// remove also the cell with the smallest distance

		dMatrix->rmCell(smallRow, rowInd);
	}
	catch(exception& e) {
		m->errorOut(e, "SingleLinkage", "update");
		exit(1);
	}
}


/***********************************************************************/
//This function updates the distance based on the nearest neighbor method.
bool SingleLinkage::updateDistance(PDistCell& colCell, PDistCell& rowCell) {
	try {
		bool changed = false;
		if (colCell.dist > rowCell.dist) {
			colCell.dist = rowCell.dist;
		}
		return(changed);
	}
	catch(exception& e) {
		m->errorOut(e, "SingleLinkage", "updateDistance");
		exit(1);
	}
}
/***********************************************************************/
