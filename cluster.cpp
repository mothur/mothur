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
		errorOut(e, "Cluster", "getRowColCells");
		exit(1);
	}
}

/***********************************************************************/

void Cluster::clusterBins(){
	try {
	
		rabund->set(smallCol, rabund->get(smallRow)+rabund->get(smallCol));	
		rabund->set(smallRow, 0);	
		rabund->setLabel(toString(smallDist));

	}
	catch(exception& e) {
		errorOut(e, "Cluster", "clusterBins");
		exit(1);
	}
}

/***********************************************************************/

void Cluster::clusterNames(){
	try {
	
		list->set(smallCol, list->get(smallRow)+','+list->get(smallCol));
		list->set(smallRow, "");	
		list->setLabel(toString(smallDist));
	
    }
	catch(exception& e) {
		errorOut(e, "Cluster", "clusterNames");
		exit(1);
	}
}

/***********************************************************************/

