/*
 *  unweighted.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "unweighted.h"

/**************************************************************************************************/

EstOutput Unweighted::getValues(Tree* t) {
	try {
		globaldata = GlobalData::getInstance();
		
		//clear out old values
		data.resize(1,0); 
		
		double UniqueBL=0.0000;  //a branch length is unique if it's chidren are from the same group
		double totalBL = 0.00;	//all branch lengths
		double UW = 0.00;		//Unweighted Value = UniqueBL / totalBL;
		
		map<string, int>::iterator it;  //iterator to traverse pgroups
		map<string, int> copyIpcount;
		
		for(int i=t->getNumLeaves();i<t->getNumNodes();i++){
		
			int lc = t->tree[i].getLChild();  //lc = vector index of left child
			int rc = t->tree[i].getRChild();  //rc = vector index of right child
			
			/**********************************************************************/
			//This section adds in all lengths that are non leaf
			
			copyIpcount = t->tree[i].pcount;
			for (it = copyIpcount.begin(); it != copyIpcount.end(); it++) {
				if (inUsersGroups(it->first, globaldata->Groups) != true) {	copyIpcount.erase(it->first);	}
			}
			
			//if i's children are from the same group then i's pcount size will be 1 
			//if copyIpcount.size() = 0 they are from a branch that is entirely from a group the user doesn't want
			if (copyIpcount.size() == 0) { }
			else if ((t->tree[i].getBranchLength() != -1) && (copyIpcount.size() == 1)) {  UniqueBL += t->tree[i].getBranchLength();	}
			
			//add i's BL to total if it is from the groups the user wants
			if ((t->tree[i].getBranchLength() != -1) && (copyIpcount.size() != 0)) {  
				totalBL += t->tree[i].getBranchLength(); 
			}
			
			/**********************************************************************/
			//This section adds in all lengths that are leaf
			
			//if i's chidren are leaves
			if (t->tree[rc].getRChild() == -1) {
				//if rc is a valid group and rc has a BL
				if ((inUsersGroups(t->tree[rc].getGroup(), globaldata->Groups) == true) && (t->tree[rc].getBranchLength() != -1)) {
					UniqueBL += t->tree[rc].getBranchLength();
					totalBL += t->tree[rc].getBranchLength(); 
				}
			}
			
			if (t->tree[lc].getLChild() == -1) {
				//if lc is a valid group and lc has a BL
				if ((inUsersGroups(t->tree[lc].getGroup(), globaldata->Groups) == true) && (t->tree[lc].getBranchLength() != -1)) {
					UniqueBL += t->tree[lc].getBranchLength();
					totalBL += t->tree[lc].getBranchLength(); 
				}
			}
			
			/**********************************************************************/
		}
		
		UW = (UniqueBL / totalBL);  
	
		if (isnan(UW) || isinf(UW)) { UW = 0; }
	
		data[0] = UW;
		
		return data;
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Unweighted class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Unweighted class function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

