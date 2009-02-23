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
		globaldata = GlobalData::getInstance();
		
		data.resize(1,0);
			
		int score = 0;
					
		for(int i=t->getNumLeaves();i<t->getNumNodes();i++){
			int lc = t->tree[i].getLChild();
			int rc = t->tree[i].getRChild();
			
			int iSize = t->tree[i].pGroups.size();
			int rcSize = t->tree[rc].pGroups.size();
			int lcSize = t->tree[lc].pGroups.size();
			
			//add in all the groups the users wanted
			for (it = t->tree[i].pGroups.begin(); it != t->tree[i].pGroups.end(); it++) {
				if (inUsersGroups(it->first, globaldata->Groups) != true) {  iSize--;  }
			}
			//add in all the groups the users wanted
			for (it = t->tree[rc].pGroups.begin(); it != t->tree[rc].pGroups.end(); it++) {
				if (inUsersGroups(it->first, globaldata->Groups) != true) {  rcSize--;  }
			}
			
			//add in all the groups the users wanted
			for (it = t->tree[lc].pGroups.begin(); it != t->tree[lc].pGroups.end(); it++) {
				if (inUsersGroups(it->first, globaldata->Groups) != true) {  lcSize--;  }
			}
			
			//if isize are 0 then that branch is to be ignored
			if (iSize == 0) { }
			else if ((rcSize == 0) || (lcSize == 0)) { }
			//if you have more groups than either of your kids then theres been a change.
			else if(iSize > rcSize || iSize > lcSize){
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

