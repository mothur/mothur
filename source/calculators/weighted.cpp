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
Weighted::Weighted(bool r, vector<string> G) : includeRoot(r), Groups(G) {
    try {
        int numGroups = Groups.size();
        
        //calculate number of comparisons i.e. with groups A,B,C = AB, AC, BC = 3;
        for (int i=0; i<numGroups; i++) {
            for (int l = 0; l < i; l++) {
                vector<string> groups; groups.push_back(Groups[i]); groups.push_back(Groups[l]);
                namesOfGroupCombos.push_back(groups);
            }
        }
    }
    catch(exception& e) {
        m->errorOut(e, "Weighted", "getValues");
        exit(1);
    }
}
/**************************************************************************************************/

EstOutput Weighted::getValues(Tree* t, int p, string o) {
    try {
		processors = p; outputDir = o;
        
        CountTable* ct = t->getCountTable();
        Treenames = t->getTreeNames();
		
		if (m->getControl_pressed()) { return data; }
		
        return (createProcesses(t, ct));
	}
	catch(exception& e) {
		m->errorOut(e, "Weighted", "getValues");
		exit(1);
	}
}
/***********************************************************************/
struct weightedData {
    int start;
    int num, count;
    MothurOut* m;
    EstOutput results;
    vector< vector<string> > namesOfGroupCombos;
    Tree* t;
    CountTable* ct;
    bool includeRoot;
    
    
    weightedData(){}
    weightedData(int st, int en, vector< vector<string> > ngc, Tree* tree, CountTable* count, bool ir) {
        m = MothurOut::getInstance();
        start = st;
        num = en;
        namesOfGroupCombos = ngc;
        t = tree;
        ct = count;
        includeRoot = ir;
        count = 0;
    }
};
/**************************************************************************************************/
double getLengthToRoot(Tree* t, bool includeRoot, int v, string groupA, string groupB, map< vector<string>, set<int> >& rootForGrouping) {
    MothurOut* m; m = MothurOut::getInstance();
    try {
        double sum = 0.0;
        int index = v;
        
        //you are a leaf
        if(t->tree[index].getBranchLength() != -1){	sum += abs(t->tree[index].getBranchLength());	}
        double tempTotal = 0.0;
        index = t->tree[index].getParent();
        
        vector<string> grouping; grouping.push_back(groupA); grouping.push_back(groupB);
        
        rootForGrouping[grouping].insert(index);
        
        //while you aren't at root
        while(t->tree[index].getParent() != -1){
            
            if (m->getControl_pressed()) {  return sum; }
            
            int parent = t->tree[index].getParent();
            
            if (includeRoot) { //add everyone
                if(t->tree[index].getBranchLength() != -1){	sum += abs(t->tree[index].getBranchLength());	}
            }else {
                
                //am I the root for this grouping? if so I want to stop "early"
                //does my sibling have descendants from the users groups?
                int lc = t->tree[parent].getLChild();
                int rc = t->tree[parent].getRChild();
                
                int sib = lc;
                if (lc == index) { sib = rc; }
                
                map<string, int>::iterator itGroup;
                int pcountSize = 0;
                itGroup = t->tree[sib].pcount.find(groupA);
                if (itGroup != t->tree[sib].pcount.end()) { pcountSize++;  }
                itGroup = t->tree[sib].pcount.find(groupB);
                if (itGroup != t->tree[sib].pcount.end()) { pcountSize++;  }
                
                //if yes, I am not the root so add me
                if (pcountSize != 0) {
                    if (t->tree[index].getBranchLength() != -1) {
                        sum += abs(t->tree[index].getBranchLength()) + tempTotal;
                        tempTotal = 0.0;
                    }else {
                        sum += tempTotal;
                        tempTotal = 0.0;
                    }
                    rootForGrouping[grouping].clear();
                    rootForGrouping[grouping].insert(parent);
                }else { //if no, I may be the root so add my br to tempTotal until I am proven innocent
                    if (t->tree[index].getBranchLength() != -1) {
                        tempTotal += abs(t->tree[index].getBranchLength());
                    }
                }
            }
            
            index = parent;	
        }
        
        //get all nodes above the root to add so we don't add their u values above
        index = *(rootForGrouping[grouping].begin());
        
        while(t->tree[index].getParent() != -1){
            int parent = t->tree[index].getParent();
            rootForGrouping[grouping].insert(parent);
            index = parent;
        }
        return sum;
    }
    catch(exception& e) {
        m->errorOut(e, "Weighted", "getBranchLengthSums");
        exit(1);
    }
}
/**************************************************************************************************/
int driverWeighted(weightedData* params) {
 try {
		vector<double> D;
		params->count = 0;
        map<string, double> WScore;
        map< vector<string>, set<int> > rootForGrouping;
     
		for (int h = params->start; h < (params->start+params->num); h++) {
		
            if (params->m->getControl_pressed()) { break; }
		
			//initialize weighted score
			string groupA = params->namesOfGroupCombos[h][0];
			string groupB = params->namesOfGroupCombos[h][1];
			
			set<int> validBranches;
            WScore[groupA+groupB] = 0.0;
			D.push_back(0.0000); //initialize a spot in D for each combination
			
			//adding the wieghted sums from group i
			for (int j = 0; j < params->t->groupNodeInfo[groupA].size(); j++) { //the leaf nodes that have seqs from group i
				map<string, int>::iterator it = params->t->tree[params->t->groupNodeInfo[groupA][j]].pcount.find(groupA);
				int numSeqsInGroupI = it->second;
				
				double sum = getLengthToRoot(params->t, params->includeRoot, params->t->groupNodeInfo[groupA][j], groupA, groupB, rootForGrouping);
				double weightedSum = ((numSeqsInGroupI * sum) / (double)params->ct->getGroupCount(groupA));
			
				D[params->count] += weightedSum;
			}
			
			//adding the wieghted sums from group l
			for (int j = 0; j < params->t->groupNodeInfo[groupB].size(); j++) { //the leaf nodes that have seqs from group l
				map<string, int>::iterator it = params->t->tree[params->t->groupNodeInfo[groupB][j]].pcount.find(groupB);
				int numSeqsInGroupL = it->second;
				
				double sum = getLengthToRoot(params->t, params->includeRoot, params->t->groupNodeInfo[groupB][j], groupA, groupB, rootForGrouping);
				double weightedSum = ((numSeqsInGroupL * sum) / (double)params->ct->getGroupCount(groupB));
			
				D[params->count] += weightedSum;
			}
			params->count++;
		}
	 
		//calculate u for the group comb 
		for (int h = params->start; h < (params->start+params->num); h++) {
            
			string groupA = params->namesOfGroupCombos[h][0];
			string groupB = params->namesOfGroupCombos[h][1];
			
			//calculate u for the group comb 
			for(int i=0;i<params->t->getNumNodes();i++){
				
                if (params->m->getControl_pressed()) { break; }
				
				double u;
				//int pcountSize = 0;
				//does this node have descendants from groupA
				map<string, int>::iterator it = params->t->tree[i].pcount.find(groupA);
				//if it does u = # of its descendants with a certain group / total number in tree with a certain group
				if (it != params->t->tree[i].pcount.end()) {
					u = (double) params->t->tree[i].pcount[groupA] / (double) params->ct->getGroupCount(groupA);
				}else { u = 0.00; }
				
				
				//does this node have descendants from group l
				it = params->t->tree[i].pcount.find(groupB);
				
				//if it does subtract their percentage from u
				if (it != params->t->tree[i].pcount.end()) {
					u -= (double) params->t->tree[i].pcount[groupB] / (double) params->ct->getGroupCount(groupB);
				}
				
				if (params->includeRoot) {
					if (params->t->tree[i].getBranchLength() != -1) {
						u = abs(u * params->t->tree[i].getBranchLength());
						WScore[(groupA+groupB)] += u; 
					}
				}else {
					//if this is not the root then add it
					if (rootForGrouping[params->namesOfGroupCombos[h]].count(i) == 0) {
						if (params->t->tree[i].getBranchLength() != -1) {
							u = abs(u * params->t->tree[i].getBranchLength());
							WScore[(groupA+groupB)] += u; 
						}
					}
				}
			}
			
		}
		
		/********************************************************/
		//calculate weighted score for the group combination
		double UN;	
		params->count = 0;
		for (int h = params->start; h < (params->start+params->num); h++) {
			UN = (WScore[params->namesOfGroupCombos[h][0]+params->namesOfGroupCombos[h][1]] / D[params->count]);
			if (isnan(UN) || isinf(UN)) { UN = 0; } 
			params->results.push_back(UN);
			params->count++;
		}
	}
	catch(exception& e) {
		params->m->errorOut(e, "Weighted", "driver");
		exit(1);
	}
}
/**************************************************************************************************/
EstOutput Weighted::getValues(Tree* t, string groupA, string groupB) { 
 try {
		
		EstOutput data;
        CountTable* ct = t->getCountTable();
        map< vector<string>, set<int> > rootForGrouping;
		
		if (m->getControl_pressed()) { return data; }
		
		//initialize weighted score
		map<string, double> WScore; WScore[(groupA+groupB)] = 0.0;
		double D = 0.0;
		set<int> validBranches;
		
		vector<string> groups; groups.push_back(groupA); groups.push_back(groupB);
		
		//adding the wieghted sums from group i
		for (int j = 0; j < t->groupNodeInfo[groups[0]].size(); j++) { //the leaf nodes that have seqs from group i
			map<string, int>::iterator it = t->tree[t->groupNodeInfo[groups[0]][j]].pcount.find(groups[0]);
			int numSeqsInGroupI = it->second;
			
			double sum = getLengthToRoot(t, includeRoot, t->groupNodeInfo[groups[0]][j], groups[0], groups[1], rootForGrouping);
			double weightedSum = ((numSeqsInGroupI * sum) / (double)ct->getGroupCount(groups[0]));
		
			D += weightedSum;
		}
		
		//adding the wieghted sums from group l
		for (int j = 0; j < t->groupNodeInfo[groups[1]].size(); j++) { //the leaf nodes that have seqs from group l
			map<string, int>::iterator it = t->tree[t->groupNodeInfo[groups[1]][j]].pcount.find(groups[1]);
			int numSeqsInGroupL = it->second;
			
			double sum = getLengthToRoot(t, includeRoot, t->groupNodeInfo[groups[1]][j], groups[0], groups[1], rootForGrouping);
			double weightedSum = ((numSeqsInGroupL * sum) / (double)ct->getGroupCount(groups[1]));
		
			D += weightedSum;
		}
				
		//calculate u for the group comb 
		for(int i=0;i<t->getNumNodes();i++){
		 
			if (m->getControl_pressed()) { return data; }
			
			double u;
			//int pcountSize = 0;
			//does this node have descendants from groupA
			map<string, int>::iterator it =  t->tree[i].pcount.find(groupA);
			//if it does u = # of its descendants with a certain group / total number in tree with a certain group
			if (it != t->tree[i].pcount.end()) {
				u = (double) t->tree[i].pcount[groupA] / (double) ct->getGroupCount(groupA);
			}else { u = 0.00; }
			
			
			//does this node have descendants from group l
			it = t->tree[i].pcount.find(groupB);
			//if it does subtract their percentage from u
			if (it != t->tree[i].pcount.end()) {
				u -= (double) t->tree[i].pcount[groupB] / (double) ct->getGroupCount(groupB);
			}
			
			if (includeRoot) {
				if (t->tree[i].getBranchLength() != -1) {
					u = abs(u * t->tree[i].getBranchLength());
					WScore[(groupA+groupB)] += u;
				}
			}else{
				//if this is not the root then add it
				if (rootForGrouping[groups].count(i) == 0) {
					if (t->tree[i].getBranchLength() != -1) {
						u = abs(u * t->tree[i].getBranchLength());
						WScore[(groupA+groupB)] += u;
					}
				}
			}
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
		m->errorOut(e, "Weighted", "getValues");
		exit(1);
	}
}
/**************************************************************************************************/

EstOutput Weighted::createProcesses(Tree* t, CountTable* ct) {
    try {
        vector<linePair> lines;
        int remainingPairs = namesOfGroupCombos.size();
        if (remainingPairs < processors) { processors = remainingPairs; }
        int startIndex = 0;
        for (int remainingProcessors = processors; remainingProcessors > 0; remainingProcessors--) {
            int numPairs = remainingPairs; //case for last processor
            if (remainingProcessors != 1) { numPairs = ceil(remainingPairs / remainingProcessors); }
            lines.push_back(linePair(startIndex, numPairs)); //startIndex, numPairs
            startIndex = startIndex + numPairs;
            remainingPairs = remainingPairs - numPairs;
        }
        
        //create array of worker threads
        vector<thread*> workerThreads;
        vector<weightedData*> data;
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            CountTable* copyCount = new CountTable();
            copyCount->copy(ct);
            Tree* copyTree = new Tree(copyCount, Treenames);
            copyTree->getCopy(t);
            
            weightedData* dataBundle = new weightedData(lines[i+1].start, lines[i+1].end, namesOfGroupCombos, copyTree, copyCount, includeRoot);
            data.push_back(dataBundle);
            
            workerThreads.push_back(new thread(driverWeighted, dataBundle));
        }
        
        weightedData* dataBundle = new weightedData(lines[0].start, lines[0].end, namesOfGroupCombos, t, ct, includeRoot);
        driverWeighted(dataBundle);
        EstOutput results = dataBundle->results;
        delete dataBundle;
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            
            for (int j = 0; j < data[i]->results.size(); j++) {  results.push_back(data[i]->results[j]);  }
            if (data[i]->count != data[i]->num) { //you didn't complete your tasks
                m->mothurOut("[ERROR]: thread " + toString(i+1) + " failed to complete it's tasks, quitting.\n");
                m->setControl_pressed(true);
            }
            
            delete data[i]->t;
            delete data[i]->ct;
            delete data[i];
            delete workerThreads[i];
        }
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "Weighted", "createProcesses");
        exit(1);
    }
}
/**************************************************************************************************/

