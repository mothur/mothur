#ifndef AVERAGE_H
#define AVERAGE_H


#include "mothur.h"
#include "cluster.hpp"
#include "rabundvector.hpp"
#include "sparsematrix.hpp"

/* This class implements the average UPGMA, average neighbor clustering algorithm */

/***********************************************************************/

AverageLinkage::AverageLinkage(RAbundVector* rav, ListVector* lv, SparseMatrix* dm, float c, string s) :
	Cluster(rav, lv, dm, c, s)
{
	saveRow = -1;
	saveCol = -1;
}


/***********************************************************************/
//This function returns the tag of the method.
string AverageLinkage::getTag() {
	return("an");
}


/***********************************************************************/
//This function updates the distance based on the average linkage method.
bool AverageLinkage::updateDistance(MatData& colCell, MatData& rowCell) {
	try {
		if ((saveRow != smallRow) || (saveCol != smallCol)) {
			rowBin = rabund->get(smallRow);
			colBin = rabund->get(smallCol);
			totalBin = rowBin + colBin;
			saveRow = smallRow;
			saveCol = smallCol;
		}
		
		colCell->dist = (colBin * colCell->dist + rowBin * rowCell->dist) / totalBin;
		
		return(true);
	}
	catch(exception& e) {
		errorOut(e, "AverageLinkage", "updateDistance");
		exit(1);
	}
}

/***********************************************************************/


/***********************************************************************/

#endif
