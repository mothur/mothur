/*
 *  weighted.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "weighted.h"

/**************************************************************************************************/

EstOutput Weighted::getValues(Tree* t) {
    try {
		globaldata = GlobalData::getInstance();
		int numGroups;
		vector<double> D;
		
		//if the user has not entered specific groups to analyze then do them all
		if (globaldata->Groups.size() == 0) {
			numGroups = tmap->getNumGroups();
		}else {
			numGroups = globaldata->Groups.size();
		}
		
		//calculate number of comparisons i.e. with groups A,B,C = AB, AC, BC = 3;
		int n = 1;
		int count = 0;
		for (int i=1; i<numGroups; i++) { 
			for (int l = n; l < numGroups; l++) {	
				//initialize weighted scores
				if (globaldata->Groups.size() == 0) {
					WScore[tmap->namesOfGroups[i-1]+tmap->namesOfGroups[l]] = 0.0;
				}else {
					WScore[globaldata->Groups[i-1]+globaldata->Groups[l]] = 0.0;
				}
				
				D.push_back(0.0000); //initialize a spot in D for each combination
				
				/********************************************************/
				//calculate a D value for each group combo
				for(int v=0;v<t->getNumLeaves();v++){
					int index = v;
					double sum = 0.0000;
		
					//while you aren't at root
					while(t->tree[index].getParent() != -1){
							
						//if you have a BL
						if(t->tree[index].getBranchLength() != -1){
							sum += t->tree[index].getBranchLength();
						}
						index = t->tree[index].getParent();
					}
						
					//get last breanch length added
					if(t->tree[index].getBranchLength() != -1){
						sum += t->tree[index].getBranchLength();
					}
						
					if (globaldata->Groups.size() == 0) {
						//is this sum from a sequence which is in one of the users groups
						if (inUsersGroups(t->tree[v].getGroup(), tmap->namesOfGroups) == true) {
							//is this sum from a sequence which is in this groupCombo
							if ((t->tree[v].getGroup() == tmap->namesOfGroups[i-1]) || (t->tree[v].getGroup() == tmap->namesOfGroups[l])) {
								sum /= (double)tmap->seqsPerGroup[t->tree[v].getGroup()];
								D[count] += sum; 
							}
						}
					}else {
						//is this sum from a sequence which is in one of the users groups
						if (inUsersGroups(t->tree[v].getGroup(), globaldata->Groups) == true) {
							//is this sum from a sequence which is in this groupCombo
							if ((t->tree[v].getGroup() == globaldata->Groups[i-1]) || (t->tree[v].getGroup() == globaldata->Groups[l])) {
								sum /= (double)tmap->seqsPerGroup[t->tree[v].getGroup()];
								D[count] += sum; 
							}
						}
					}
				}
				/*********************************************************/
				count++;
			}
			n++;
		}
		
		data.clear(); //clear out old values
	
		for(int i=0;i<t->getNumNodes();i++){
			//calculate weighted score for each of the group comb i.e. with groups A,B,C = AB, AC, BC.
			n = 1;
			for (int b=1; b<numGroups; b++) { 
				for (int l = n; l < numGroups; l++) {
					//calculate a u value for each combo
					double u;
					//the user has not entered specific groups
					if (globaldata->Groups.size() == 0) {
						//does this node have descendants from group b-1
						it = t->tree[i].pcount.find(tmap->namesOfGroups[b-1]);
						//if it does u = # of its descendants with a certain group / total number in tree with a certain group
						if (it != t->tree[i].pcount.end()) {
							u = (double) t->tree[i].pcount[tmap->namesOfGroups[b-1]] / (double) tmap->seqsPerGroup[tmap->namesOfGroups[b-1]];
						}else { u = 0.00; }
		
						//does this node have descendants from group l
						it = t->tree[i].pcount.find(tmap->namesOfGroups[l]);
						//if it does subtract their percentage from u
						if (it != t->tree[i].pcount.end()) {
							u -= (double) t->tree[i].pcount[tmap->namesOfGroups[l]] / (double) tmap->seqsPerGroup[tmap->namesOfGroups[l]];
						}
						
						u = abs(u) * t->tree[i].getBranchLength();
					
						//save groupcombs u value
						WScore[tmap->namesOfGroups[b-1]+tmap->namesOfGroups[l]] += u;
						
					//the user has entered specific groups	
					}else {
						//does this node have descendants from group b-1
						it = t->tree[i].pcount.find(globaldata->Groups[b-1]);
						//if it does u = # of its descendants with a certain group / total number in tree with a certain group
						if (it != t->tree[i].pcount.end()) {
							u = (double) t->tree[i].pcount[globaldata->Groups[b-1]] / (double) tmap->seqsPerGroup[globaldata->Groups[b-1]];
						}else { u = 0.00; }
		
						//does this node have descendants from group l
						it = t->tree[i].pcount.find(globaldata->Groups[l]);
						//if it does subtract their percentage from u
						if (it != t->tree[i].pcount.end()) {
							u -= (double) t->tree[i].pcount[globaldata->Groups[l]] / (double) tmap->seqsPerGroup[globaldata->Groups[l]];
						}
						
						u = abs(u) * t->tree[i].getBranchLength();
					
						//save groupcombs u value
						WScore[globaldata->Groups[b-1]+globaldata->Groups[l]] += u;
					}
				/*********************************************************/
				}
				n++;
			}
		}
  
		//calculate weighted score for each group combination
		double UN;	
		n = 1;
		count = 0;
		for (int i=1; i<numGroups; i++) { 
			for (int l = n; l < numGroups; l++) {
				//the user has not entered specific groups
				if (globaldata->Groups.size() == 0) {
					UN = (WScore[tmap->namesOfGroups[i-1]+tmap->namesOfGroups[l]] / D[count]);
				}else {//they have entered specific groups
					UN = (WScore[globaldata->Groups[i-1]+globaldata->Groups[l]] / D[count]);
				}
				if (isnan(UN) || isinf(UN)) { UN = 0; } 
				data.push_back(UN);
				count++;
			}
			n++;
		}
		return data;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Weighted class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Weighted class function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

/**************************************************************************************************/
EstOutput Weighted::getValues(Tree* t, string groupA, string groupB) { 
 try {
		globaldata = GlobalData::getInstance();
		
		data.clear(); //clear out old values
		
		//initialize weighted score
		WScore[(groupA+groupB)] = 0.0;
		float D = 0.0;
		
						
		/********************************************************/
		//calculate a D value for the group combo
		for(int v=0;v<t->getNumLeaves();v++){
			int index = v;
			double sum = 0.0000;
		
			//while you aren't at root
			while(t->tree[index].getParent() != -1){
							
				//if you have a BL
				if(t->tree[index].getBranchLength() != -1){
					sum += t->tree[index].getBranchLength();
				}
				index = t->tree[index].getParent();
			}
						
			//get last breanch length added
			if(t->tree[index].getBranchLength() != -1){
				sum += t->tree[index].getBranchLength();
			}
						
			if ((t->tree[v].getGroup() == groupA) || (t->tree[v].getGroup() == groupB)) {
				sum /= (double)tmap->seqsPerGroup[t->tree[v].getGroup()];
				D += sum; 
			}
		}
		/********************************************************/				
		
		//calculate u for the group comb 
		for(int i=0;i<t->getNumNodes();i++){
			double u;
			//does this node have descendants from groupA
			it = t->tree[i].pcount.find(groupA);
			//if it does u = # of its descendants with a certain group / total number in tree with a certain group
			if (it != t->tree[i].pcount.end()) {
				u = (double) t->tree[i].pcount[groupA] / (double) tmap->seqsPerGroup[groupA];
			}else { u = 0.00; }
		
						
			//does this node have descendants from group l
			it = t->tree[i].pcount.find(groupB);
			//if it does subtract their percentage from u
			if (it != t->tree[i].pcount.end()) {
				u -= (double) t->tree[i].pcount[groupB] / (double) tmap->seqsPerGroup[groupB];
			}
						
			u = abs(u) * t->tree[i].getBranchLength();
					
			//save groupcombs u value
			WScore[(groupA+groupB)] += u;
		}
		
		/********************************************************/
		
		//calculate weighted score for the group combination
		double UN;	
		UN = (WScore[(groupA+groupB)] / D);
		
		if (isnan(UN) || isinf(UN)) { UN = 0; } 
		data.push_back(UN);
				
		return data; 
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Weighted class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Weighted class function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}









