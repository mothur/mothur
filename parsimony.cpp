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
		vector<string> groups;
		
		copyTree = new Tree();
		
		//if the users enters no groups then give them the score of all groups
		int numGroups = globaldata->Groups.size();
		
		//calculate number of comparsions
		int numComp = 0;
		for (int r=0; r<numGroups; r++) { 
			for (int l = r+1; l < numGroups; l++) {
				numComp++;
			}
		}

		//numComp+1 for AB, AC, BC, ABC
		data.resize(numComp+1,0);
		
		int count = 0;
		for (int a=0; a<numGroups; a++) { 
			for (int l = 0; l < a; l++) {
				int score = 0;
				
				//groups in this combo
				groups.push_back(globaldata->Groups[a]); groups.push_back(globaldata->Groups[l]);
				
				//copy users tree so that you can redo pgroups 
				copyTree->getCopy(t);

				//create pgroups that reflect the groups the user want to use
				for(int i=copyTree->getNumLeaves();i<copyTree->getNumNodes();i++){
					copyTree->tree[i].pGroups = (copyTree->mergeUserGroups(i, groups));
				}
		
				for(int i=copyTree->getNumLeaves();i<copyTree->getNumNodes();i++){
				
					if (m->control_pressed) { return data; }
					
					int lc = copyTree->tree[i].getLChild();
					int rc = copyTree->tree[i].getRChild();
			
					int iSize = copyTree->tree[i].pGroups.size();
					int rcSize = copyTree->tree[rc].pGroups.size();
					int lcSize = copyTree->tree[lc].pGroups.size();
		
					//if isize are 0 then that branch is to be ignored
					if (iSize == 0) { }
					else if ((rcSize == 0) || (lcSize == 0)) { }
					//if you have more groups than either of your kids then theres been a change.
					else if(iSize > rcSize || iSize > lcSize){
						score++;
					}
				} 
				
				data[count] = score;
				count++;
				groups.clear();
			}
		}
		
		if (numComp != 1) {
			if (numGroups == 0) {
				//get score for all users groups
				for (int i = 0; i < tmap->namesOfGroups.size(); i++) {
					if (tmap->namesOfGroups[i] != "xxx") {
						groups.push_back(tmap->namesOfGroups[i]);
					}
				}
			}else {
				for (int i = 0; i < globaldata->Groups.size(); i++) {
					groups.push_back(globaldata->Groups[i]);
				}
			}
			
			//copy users tree so that you can redo pgroups 
			copyTree->getCopy(t);
			int score = 0;
		
			//create pgroups that reflect the groups the user want to use
			for(int i=copyTree->getNumLeaves();i<copyTree->getNumNodes();i++){
				copyTree->tree[i].pGroups = (copyTree->mergeUserGroups(i, groups));
			}
		
//			map<string,int>::iterator it;
			
			for(int i=copyTree->getNumLeaves();i<copyTree->getNumNodes();i++){
			
				if (m->control_pressed) { return data; }
				
				int lc = copyTree->tree[i].getLChild();
				int rc = copyTree->tree[i].getRChild();
			
				int iSize = copyTree->tree[i].pGroups.size();
				int rcSize = copyTree->tree[rc].pGroups.size();
				int lcSize = copyTree->tree[lc].pGroups.size();
				
					
				//if isize are 0 then that branch is to be ignored
				if (iSize == 0) { }
				else if ((rcSize == 0) || (lcSize == 0)) { }
				//if you have more groups than either of your kids then theres been a change.
				else if(iSize > rcSize || iSize > lcSize){
					score++;
				}
			} 
		
			data[count] = score;

		}
		
		delete copyTree;
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Parsimony", "getValues");
		exit(1);
	}
}

/**************************************************************************************************/

