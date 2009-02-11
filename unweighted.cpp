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
	
		//clear out old values
		data.resize(1,0); 
		penalty.resize(t->getNumLeaves(), 0);
		
		map<string,double> unique;  //group, total of all branch lengths of nodes with that group.
		double shared = 0.0000;
		double UW=0.0000;
		
		//add up the branch lengths for each group. 
		for(int i=0;i<t->getNumLeaves();i++){
			if(t->tree[i].pGroups.size() > 0){
				unique[t->tree[i].pGroups.begin()->first] += t->tree[i].getBranchLength();
			}
		}
		
		//for each non-leaf node
		for(int i=t->getNumLeaves();i<t->getNumNodes();i++){
		
			int lc = t->tree[i].getLChild();  //lc = vector index of left child
			int rc = t->tree[i].getRChild();  //rc = vector index of right child
			
			//get penalty values
			if(t->tree[rc].pGroups.size() == 0 || t->tree[lc].pGroups.size() == 0){
				penalty.push_back(penalty[t->tree[rc].getIndex()]+penalty[t->tree[lc].getIndex()]);
			}
			else if(t->tree[i].pGroups.size() > t->tree[rc].pGroups.size() || t->tree[i].pGroups.size() > t->tree[lc].pGroups.size()){
				penalty.push_back(penalty[t->tree[rc].getIndex()]+penalty[t->tree[lc].getIndex()]+1);
			}
			else{
				penalty.push_back(penalty[t->tree[rc].getIndex()]+penalty[t->tree[lc].getIndex()]);
			}

			//not sure when this would ever be true??? if your parent is root could be, but pGroups.size() should never be 0.
			if(t->tree[i].getParent() == -1 && (t->tree[lc].pGroups.size() == 0 || t->tree[rc].pGroups.size() == 0)){
				shared -= 1; 
			}
			else if(penalty[i] != 0 && t->tree[i].pGroups.size() != 0){
				shared += t->tree[i].getBranchLength();
			}
			else if( t->tree[i].pGroups.size() != 0){
				unique[t->tree[i].pGroups.begin()->first] += t->tree[i].getBranchLength();	    
			}
		}
    
		map<string,double>::iterator pos;
		for(pos=unique.begin();pos!=unique.end();pos++){
			if(pos->first!="xxx"){	    
				UW += unique[pos->first];
			}
		}
	
		UW /= (UW + shared);
	
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

/**************************************************************************************************/