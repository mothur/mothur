/*
 *  parsimony.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/26/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "parsimony.h"

/**************************************************************************************************/
EstOutput Parsimony::getValues(Tree* t) {
	try {
		data.resize(1,0);
			
		int score = 0;
			
		for(int i=t->getNumLeaves();i<t->getNumNodes();i++){
			int lc = t->tree[i].getLChild();
			int rc = t->tree[i].getRChild();
			
			 if(t->tree[i].pGroups.size() > t->tree[rc].pGroups.size() || t->tree[i].pGroups.size() > t->tree[lc].pGroups.size()){
				score++;
			}
		} 
		
		data[0] = score;
		
		return data;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Parsimony class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Parsimony class function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
/**************************************************************************************************/
