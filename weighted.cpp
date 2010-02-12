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
		
		numGroups = globaldata->Groups.size();
		
		//calculate number of comparisons i.e. with groups A,B,C = AB, AC, BC = 3;
		int count = 0;
		for (int i=0; i<numGroups; i++) { 
			for (int l = i+1; l < numGroups; l++) {	
				//initialize weighted scores
				WScore[globaldata->Groups[i]+globaldata->Groups[l]] = 0.0;
				
				vector<string> groups; groups.push_back(globaldata->Groups[i]); groups.push_back(globaldata->Groups[l]);
				
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
							sum += abs(t->tree[index].getBranchLength());
						}
						index = t->tree[index].getParent();
					}
						
					//get last breanch length added
					if(t->tree[index].getBranchLength() != -1){
						sum += abs(t->tree[index].getBranchLength());
					}
						
					//is this sum from a sequence which is in one of the users groups
					if (inUsersGroups(t->tree[v].getGroup(), globaldata->Groups) == true) {
						//is this sum from a sequence which is in this groupCombo
						if (inUsersGroups(t->tree[v].getGroup(), groups)) {
							int numSeqsInGroupI, numSeqsInGroupL;
							
							map<string, int>::iterator it;
							it = t->tree[v].pcount.find(groups[0]);
							if (it != t->tree[v].pcount.end()) { //this leaf node contains seqs from group i
								numSeqsInGroupI = it->second;
							}else{ numSeqsInGroupI = 0; }
							
							it = t->tree[v].pcount.find(groups[1]);
							if (it != t->tree[v].pcount.end()) { //this leaf node contains seqs from group l
								numSeqsInGroupL = it->second;
							}else{ numSeqsInGroupL = 0; }
							
							double weightedSum = ((numSeqsInGroupI * sum) / (double)tmap->seqsPerGroup[groups[0]]) + ((numSeqsInGroupL * sum) / (double)tmap->seqsPerGroup[groups[1]]);

							//sum /= (double)tmap->seqsPerGroup[t->tree[v].getGroup()];
							
							D[count] += weightedSum; 
						}
					}
				}
				/*********************************************************/
				count++;
			}
		}
		
		data.clear(); //clear out old values
	
		for(int i=0;i<t->getNumNodes();i++){
			//calculate weighted score for each of the group comb i.e. with groups A,B,C = AB, AC, BC.
			for (int b=0; b<numGroups; b++) { 
				for (int l = b+1; l < numGroups; l++) {
					//calculate a u value for each combo
					double u;
					//does this node have descendants from group b-1
					it = t->tree[i].pcount.find(globaldata->Groups[b]);
					//if it does u = # of its descendants with a certain group / total number in tree with a certain group
					if (it != t->tree[i].pcount.end()) {
						u = (double) t->tree[i].pcount[globaldata->Groups[b]] / (double) tmap->seqsPerGroup[globaldata->Groups[b]];
					}else { u = 0.00; }
		
					//does this node have descendants from group l
					it = t->tree[i].pcount.find(globaldata->Groups[l]);
					//if it does subtract their percentage from u
					if (it != t->tree[i].pcount.end()) {
						u -= (double) t->tree[i].pcount[globaldata->Groups[l]] / (double) tmap->seqsPerGroup[globaldata->Groups[l]];
					}
						
					u = abs(u * t->tree[i].getBranchLength());
					
					//save groupcombs u value
					WScore[globaldata->Groups[b]+globaldata->Groups[l]] += u;
				/*********************************************************/
				}
			}
		}
  
		//calculate weighted score for each group combination
		double UN;	
		count = 0;
		for (int i=0; i<numGroups; i++) { 
			for (int l = i+1; l < numGroups; l++) {
				UN = (WScore[globaldata->Groups[i]+globaldata->Groups[l]] / D[count]);
		
				if (isnan(UN) || isinf(UN)) { UN = 0; } 
				data.push_back(UN);
				count++;
			}
		}
		return data;
	}
	catch(exception& e) {
		errorOut(e, "Weighted", "getValues");
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
		
		vector<string> groups; groups.push_back(groupA); groups.push_back(groupB);
						
		/********************************************************/
		//calculate a D value for the group combo
		for(int v=0;v<t->getNumLeaves();v++){
			int index = v;
			double sum = 0.0000;
		
			//while you aren't at root
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
						
			if (inUsersGroups(t->tree[v].getGroup(), groups)) {
				int numSeqsInGroupI, numSeqsInGroupL;
							
				map<string, int>::iterator it;
				it = t->tree[v].pcount.find(groups[0]);
				if (it != t->tree[v].pcount.end()) { //this leaf node contains seqs from group i
					numSeqsInGroupI = it->second;
				}else{ numSeqsInGroupI = 0; }
				
				it = t->tree[v].pcount.find(groups[1]);
				if (it != t->tree[v].pcount.end()) { //this leaf node contains seqs from group l
					numSeqsInGroupL = it->second;
				}else{ numSeqsInGroupL = 0; }
				
				double weightedSum = ((numSeqsInGroupI * sum) / (double)tmap->seqsPerGroup[groups[0]]) + ((numSeqsInGroupL * sum) / (double)tmap->seqsPerGroup[groups[1]]);
				
				//sum /= (double)tmap->seqsPerGroup[t->tree[v].getGroup()];
				
				D += weightedSum; 
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
						
			u = abs(u * t->tree[i].getBranchLength());
					
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
		errorOut(e, "Weighted", "getValues");
		exit(1);
	}
}









