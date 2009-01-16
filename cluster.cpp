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
#include <exception>

/***********************************************************************/

Cluster::Cluster(RAbundVector* rav, ListVector* lv, SparseMatrix* dm) :
rabund(rav), list(lv), dMatrix(dm)
{
}

/***********************************************************************/

void Cluster::getRowColCells(){
	try {
		PCell* smallCell = dMatrix->getSmallestCell();	//find the smallest cell - this routine should probably not be in the SpMat class
	
		smallRow = smallCell->row;		//get its row
		smallCol = smallCell->column;	//get its column
		smallDist = smallCell->dist;	//get the smallest distance
	
		rowCells.clear();
		colCells.clear();
		
		for(MatData currentCell=dMatrix->begin();currentCell!=dMatrix->end();currentCell++){
		
			if(&*currentCell == smallCell){				//put the smallest cell first
				rowCells.insert(rowCells.begin(), currentCell);
				colCells.insert(colCells.begin(), currentCell);
			}
			else if(currentCell->row == smallRow){
				rowCells.push_back(currentCell);
			}
			else if(currentCell->column == smallRow){
				rowCells.push_back(currentCell);
			}
			else if(currentCell->row == smallCol){
				colCells.push_back(currentCell);
			}
			else if(currentCell->column == smallCol){
				colCells.push_back(currentCell);
			}
		}
	
		nRowCells = rowCells.size();
		nColCells = colCells.size();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Cluster class Function getRowColCells. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Cluster class function getRowColCells. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
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
		cout << "Standard Error: " << e.what() << " has occurred in the Cluster class Function clusterBins. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Cluster class function clusterBins. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	


}

/***********************************************************************/

void Cluster::clusterNames(){
	try {
	//	cout << smallCol << '\t' << smallRow << '\t' << smallDist << '\t' << list->get(smallRow) << '\t' << list->get(smallCol);

		list->set(smallCol, list->get(smallRow)+','+list->get(smallCol));
		list->set(smallRow, "");	
		list->setLabel(toString(smallDist));
	
	//	cout << '\t' << list->get(smallRow) << '\t' << list->get(smallCol) << endl;
    }
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Cluster class Function clusterNames. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Cluster class function clusterNames. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/

