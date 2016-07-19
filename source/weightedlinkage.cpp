#ifndef WEIGHTEDLINKAGE_H
#define WEIGHTEDLINKAGE_H

#include "cluster.hpp"

/* This class implements the WPGMA, weighted average neighbor clustering algorithm */

/***********************************************************************/

WeightedLinkage::WeightedLinkage(RAbundVector* rav, ListVector* lv, SparseDistanceMatrix* dm, float c, string s, float a) :
	Cluster(rav, lv, dm, c, s, a)
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
bool WeightedLinkage::updateDistance(PDistCell& colCell, PDistCell& rowCell) {
	try {
		if ((saveRow != smallRow) || (saveCol != smallCol)) {
//			rowBin = rabund->get(smallRow);
//			colBin = rabund->get(smallCol);
//			totalBin = rowBin + colBin;
			saveRow = smallRow;
			saveCol = smallCol;
		}
		
		colCell.dist = (colCell.dist + rowCell.dist) / 2.0;
		
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
