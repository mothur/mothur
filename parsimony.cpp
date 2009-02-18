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
			
			int iSize = 0;
			int rcSize = 0;
			int lcSize = 0;

			//add in all the groups the users wanted
			for (it = t->tree[i].pGroups.begin(); it != t->tree[i].pGroups.end(); it++) {
				if (inUsersGroups(it->first) == true) {  iSize++;  }
			}

			//if that leaves no groups give it 1 so it will cause no change to parent
			if (iSize == 0) { iSize++; }
			
			//add in all the groups the users wanted
			for (it = t->tree[rc].pGroups.begin(); it != t->tree[rc].pGroups.end(); it++) {

				if (inUsersGroups(it->first) == true) {  rcSize++;  }
			}
			
			//if that leaves no groups give it 1 so it will cause no change to parent
			if (rcSize == 0) { rcSize++; }

				
			//add in all the groups the users wanted
			for (it = t->tree[lc].pGroups.begin(); it != t->tree[lc].pGroups.end(); it++) {

				if (inUsersGroups(it->first) == true) {  lcSize++;  }
			}
			
			//if that leaves no groups give it 1 so it will cause no change to parent
			if (lcSize == 0) { lcSize++; }


			//if you have more groups than either of your kids then theres been a change.
			 if(iSize > rcSize || iSize > lcSize){
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

bool Parsimony::inUsersGroups(string groupname) {
	try {
		for (int i = 0; i < globaldata->Groups.size(); i++) {
			if (groupname == globaldata->Groups[i]) { return true; }
		}
		return false;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Parsimony class Function inUsersGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Parsimony class function inUsersGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/**************************************************************************************************/
