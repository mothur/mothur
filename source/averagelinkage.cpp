#ifndef AVERAGE_H
#define AVERAGE_H

//test
#include "mothur.h"
#include "cluster.hpp"
#include "rabundvector.hpp"
#include "sparsedistancematrix.h"

/* This class implements the average UPGMA, average neighbor clustering algorithm */

/***********************************************************************/

AverageLinkage::AverageLinkage(RAbundVector* rav, ListVector* lv, SparseDistanceMatrix* dm, float c, string s, float a) :
Cluster(rav, lv, dm, c, s, a)
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
bool AverageLinkage::updateDistance(PDistCell& colCell, PDistCell& rowCell) {
	try {
		if ((saveRow != smallRow) || (saveCol != smallCol)) {
			rowBin = rabund->get(smallRow);
			colBin = rabund->get(smallCol);
			totalBin = rowBin + colBin;
			saveRow = smallRow;
			saveCol = smallCol;
		}
		//cout << "colcell.dist = " << colCell.dist << '\t' << smallRow << '\t' << smallCol << '\t' << rowCell.dist << endl;
		colCell.dist = (colBin * colCell.dist + rowBin * rowCell.dist) / totalBin;
        
		return(true);
	}
	catch(exception& e) {
		m->errorOut(e, "AverageLinkage", "updateDistance");
		exit(1);
	}
}

/***********************************************************************/


/***********************************************************************/

#endif
