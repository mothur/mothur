#ifndef WEIGHTEDLINKAGE_H
#define WEIGHTEDLINKAGE_H


#include "mothur.h"
#include "cluster.hpp"
#include "rabundvector.hpp"
#include "sparsematrix.hpp"

/* This class implements the WPGMA, weighted average neighbor clustering algorithm */

/***********************************************************************/

WeightedLinkage::WeightedLinkage(RAbundVector* rav, ListVector* lv, SparseMatrix* dm, float c, string s) :
	Cluster(rav, lv, dm, c, s)
{
	saveRow = -1;
	saveCol = -1;
}


/***********************************************************************/
//This function returns the tag of the method.
string WeightedLinkage::getTag() {
	return("wn");
}


/***********************************************************************/
//This function updates the distance based on the average linkage method.
bool WeightedLinkage::updateDistance(MatData& colCell, MatData& rowCell) {
	try {
		if ((saveRow != smallRow) || (saveCol != smallCol)) {
//			rowBin = rabund->get(smallRow);
//			colBin = rabund->get(smallCol);
//			totalBin = rowBin + colBin;
			saveRow = smallRow;
			saveCol = smallCol;
		}
		
		colCell->dist = (colCell->dist + rowCell->dist) / 2.0;
		
		return(true);
	}
	catch(exception& e) {
		m->errorOut(e, "WeightedLinkage", "updateDistance");
		exit(1);
	}
}

/***********************************************************************/


/***********************************************************************/

#endif
