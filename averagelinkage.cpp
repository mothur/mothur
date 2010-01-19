#ifndef AVERAGE_H
#define AVERAGE_H


#include "mothur.h"
#include "cluster.hpp"
#include "rabundvector.hpp"
#include "sparsematrix.hpp"

/* This class implements the average UPGMA, average neighbor clustering algorithm */

/***********************************************************************/

AverageLinkage::AverageLinkage(RAbundVector* rav, ListVector* lv, SparseMatrix* dm, float c) :
	Cluster(rav, lv, dm, c)
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
		
		float oldColCell = colCell->dist;
		
		colCell->dist = (colBin * colCell->dist + rowBin * rowCell->dist) / totalBin;
		
		//warn user if merge with value above cutoff produces a value below cutoff
		if ((colCell->dist < cutoff) && ((oldColCell > cutoff) || (rowCell->dist > cutoff)) ) {
			mothurOut("Warning: merging " + toString(oldColCell) + " with " + toString(rowCell->dist) + ", new value = " + toString(colCell->dist) + ". Results will differ from those if cutoff was used in the read.dist command."); mothurOutEndLine();
		}

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
