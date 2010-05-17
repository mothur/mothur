/*
 *  phylodiversity.cpp
 *  Mothur
 *
 *  Created by westcott on 4/30/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "phylodiversity.h"

/**************************************************************************************************/
EstOutput PhyloDiversity::getValues(Tree* t, vector<int> treeNodes) {
    try {
		
		map<string, float> DScore;
		float totalLength = 0.0;
		data.clear();
		
		//initialize Dscore
		for (int i=0; i<globaldata->Groups.size(); i++) {		DScore[globaldata->Groups[i]] = 0.0;	}
	
		/********************************************************/
		//calculate a D value for each group 
		for(int v=0;v<treeNodes.size();v++){
				
			if (m->control_pressed) { return data; }
			
			//is this node from a sequence which is in one of the users groups
			if (inUsersGroups(t->tree[treeNodes[v]].getGroup(), globaldata->Groups) == true) {
					
				//calc the branch length
				//while you aren't at root
				float sum = 0.0;
				int index = treeNodes[v];

				while(t->tree[index].getParent() != -1){
					
					//if you have a BL
					if(t->tree[index].getBranchLength() != -1){
						sum += abs(t->tree[index].getBranchLength());
					}
					index = t->tree[index].getParent();
				}
					
				//get last breanch length added
				if(t->tree[index].getBranchLength() != -1){
					sum += abs(t->tree[index].getBranchLength());
				}
					
				//for each group in the groups update the total branch length accounting for the names file
				vector<string> groups = t->tree[treeNodes[v]].getGroup();
				for (int j = 0; j < groups.size(); j++) {
					int numSeqsInGroupJ = 0;
					map<string, int>::iterator it;
					it = t->tree[treeNodes[v]].pcount.find(groups[j]);
					if (it != t->tree[treeNodes[v]].pcount.end()) { //this leaf node contains seqs from group j
						numSeqsInGroupJ = it->second;
					}

					//add branch length to total for group
					DScore[groups[j]] += (numSeqsInGroupJ * sum);
				}
			}
		}
		
	
		for (int i=0; i<globaldata->Groups.size(); i++) {   
			//if (groupTotals[globaldata->Groups[i]] != 0.0) { //avoid divide by zero error
				//float percent = DScore[globaldata->Groups[i]] / groupTotals[globaldata->Groups[i]];
			float percent = DScore[globaldata->Groups[i]]; 
			data.push_back(percent);  
			//}else {  data.push_back(0.0);  }
		}
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloDiversity", "getValues");
		exit(1);
	}
}
/**************************************************************************************************/
void PhyloDiversity::setTotalGroupBranchLengths(Tree* t) {
    try {
		
		groupTotals.clear();
		
		//initialize group totals
		for (int i=0; i<globaldata->Groups.size(); i++) {		groupTotals[globaldata->Groups[i]] = 0.0;	}
		
		
		/********************************************************/
		//calculate a D value for each group 
		for(int v=0;v<t->getNumLeaves();v++){
			
			//is this node from a sequence which is in one of the users groups
			if (inUsersGroups(t->tree[v].getGroup(), globaldata->Groups) == true) {
			
				//calc the branch length
				int index = v;
				float sum = 0.0;
				
				while(t->tree[index].getParent() != -1){  //while you aren't at root
						
					//if you have a BL
					if(t->tree[index].getBranchLength() != -1){
						sum += abs(t->tree[index].getBranchLength());
					}
					index = t->tree[index].getParent();
				}
					
				//get last breanch length added
				if(t->tree[index].getBranchLength() != -1){
					sum += abs(t->tree[index].getBranchLength());
				}
							
				//account for the names file
				vector<string> groups = t->tree[v].getGroup();
				for (int j = 0; j < groups.size(); j++) {
					int numSeqsInGroupJ = 0;
					map<string, int>::iterator it;
					it = t->tree[v].pcount.find(groups[j]);
					if (it != t->tree[v].pcount.end()) { //this leaf node contains seqs from group j
						numSeqsInGroupJ = it->second;
					}

					//add branch length to total for group
					groupTotals[groups[j]] += (numSeqsInGroupJ * sum);
				}//end for
			}//end if
		}//end for
		
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloDiversity", "setTotalGroupBranchLengths");
		exit(1);
	}
}
/**************************************************************************************************/


