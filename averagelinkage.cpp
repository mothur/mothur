#ifndef AVERAGE_H
#define AVERAGE_H


#include "mothur.h"
#include "cluster.hpp"
#include "rabundvector.hpp"
#include "sparsematrix.hpp"

/* This class implements the average UPGMA, average neighbor clustering algorithm */

/***********************************************************************/

AverageLinkage::AverageLinkage(RAbundVector* rav, ListVector* lv, SparseMatrix* dm) :
Cluster(rav, lv, dm)
{}

/***********************************************************************/
//THis function clusters based on the average method
void AverageLinkage::update(){
	try{
		getRowColCells();		
	
		vector<int> found(nColCells, 0);
	
		int rowBin = rabund->get(smallRow);
		int colBin = rabund->get(smallCol);
		int totalBin = rowBin + colBin;
	
		for(int i=1;i<nRowCells;i++){
		
			int search;
		
			if(rowCells[i]->row == smallRow){
				search = rowCells[i]->column;
			}
			else{
				search = rowCells[i]->row;
			}
		
			for(int j=1;j<nColCells;j++){
			
				if(colCells[j]->row == search || colCells[j]->column == search){
					colCells[j]->dist = (colBin * colCells[j]->dist + rowBin * rowCells[i]->dist) / totalBin;		
				
					found[j] = 1;
				
					if(colCells[j]->vectorMap != NULL){
						*(colCells[j]->vectorMap) = NULL;
						colCells[j]->vectorMap = NULL;
					}
				
					break;
				}
			
			}			
			dMatrix->rmCell(rowCells[i]);
		}	
	
		clusterBins();
		clusterNames();
	
		for(int i=0;i<nColCells;i++){
			if(found[i] == 0){
				dMatrix->rmCell(colCells[i]);
			}
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the AverageLinkage class Function update. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the AverageLinkage class function update. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

	
}

/***********************************************************************/

#endif


