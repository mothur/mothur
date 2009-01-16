
#include "cluster.hpp"
#include <exception>

/***********************************************************************/

SingleLinkage::SingleLinkage(RAbundVector* rav, ListVector* lv, SparseMatrix* dm) :
Cluster(rav, lv, dm)
{}

/***********************************************************************/
//This function clusters based on nearest method.
void SingleLinkage::update(){
	try {
		getRowColCells();		
	
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
				
					if(colCells[j]->dist > rowCells[i]->dist){
						colCells[j]->dist = rowCells[i]->dist;
					
						if(colCells[j]->vectorMap != NULL){
							*(colCells[j]->vectorMap) = NULL;
							colCells[j]->vectorMap = NULL;
						}
					
					}
					dMatrix->rmCell(rowCells[i]);
					break;
				}
			}			
		
			if(search < smallCol){
				rowCells[i]->row = smallCol;
				rowCells[i]->column = search;
			}
			else{
				rowCells[i]->row = search;
				rowCells[i]->column = smallCol;
			}
		
		}	
		clusterBins();
		clusterNames();
		dMatrix->rmCell(rowCells[0]);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SingleLinkage class Function update. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SingleLinkage class function update. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/
