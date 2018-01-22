
#include "cluster.hpp"


/***********************************************************************/

SingleLinkage::SingleLinkage(RAbundVector* rav, ListVector* lv, SparseDistanceMatrix* dm, float c, string s, float a) :
Cluster(rav, lv, dm, c, s, a)
{}


/***********************************************************************/
//This function returns the tag of the method.
string SingleLinkage::getTag() {
	return("nn");
}
/***********************************************************************/
//This function updates the distance based on the nearest neighbor method.
bool SingleLinkage::updateDistance(PDistCell& colCell, PDistCell& rowCell) {
	try {
		bool changed = false;
		if (colCell.dist > rowCell.dist) {
			colCell.dist = rowCell.dist;
		}
		return(changed);
	}
	catch(exception& e) {
		m->errorOut(e, "SingleLinkage", "updateDistance");
		exit(1);
	}
}
/***********************************************************************/
