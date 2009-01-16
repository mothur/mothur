
#include "cluster.hpp"
#include <exception>

/***********************************************************************/

CompleteLinkage::CompleteLinkage(RAbundVector* rav, ListVector* lv, SparseMatrix* dm) :
Cluster(rav, lv, dm)
{}

/***********************************************************************/
//This function clusters based on the furthest neighbor method.
void CompleteLinkage::update(){
	try {
		getRowColCells();	
	
		vector<int> found(nColCells, 0);
	
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
				
					if(colCells[j]->dist < rowCells[i]->dist){
						colCells[j]->dist = rowCells[i]->dist;
					
						if(colCells[j]->vectorMap != NULL){
							*(colCells[j]->vectorMap) = NULL;
							colCells[j]->vectorMap = NULL;
						}
					
					}
				
					found[j] = 1;
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
		cout << "Standard Error: " << e.what() << " has occurred in the CompleteLinkage class Function update. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the CompleteLinkage class function update. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

/***********************************************************************/
