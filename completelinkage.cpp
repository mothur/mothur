
#include "cluster.hpp"

/***********************************************************************/

CompleteLinkage::CompleteLinkage(RAbundVector* rav, ListVector* lv, SparseMatrix* dm) :
	Cluster(rav, lv, dm)
{}

/***********************************************************************/
//This function returns the tag of the method.
string CompleteLinkage::getTag() {
	return("fn");
}


/***********************************************************************/
//This function updates the distance based on the furthest neighbor method.
bool CompleteLinkage::updateDistance(MatData& colCell, MatData& rowCell) {
	try {
		bool changed = false;
		if (colCell->dist < rowCell->dist) {
			colCell->dist = rowCell->dist;
			changed = true;
		}	
		return(changed);
	}
	catch(exception& e) {
		errorOut(e, "CompleteLinkage", "updateDistance");
		exit(1);
	}
}

/***********************************************************************/
