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

EstOutput Weighted::getValues(Tree* t, int p) {
    try {
		processors = p;
				
        return (createProcesses(t));
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
void getRoot(MothurOut* m, Tree* t, int v, vector<string> grouping, set<int>& rootForGrouping) {
    try {
        //you are a leaf so get your parent
        int index = t->tree[v].getParent();
        
        //my parent is a potential root
        rootForGrouping.insert(index);
        
        //while you aren't at root
        while(t->tree[index].getParent() != -1){
            
            if (m->getControl_pressed()) {  return; }
            
            //am I the root for this grouping? if so I want to stop "early"
            //does my sibling have descendants from the users groups?
            //if so I am not the root
            int parent = t->tree[index].getParent();
            int lc = t->tree[parent].getLChild();
            int rc = t->tree[parent].getRChild();
            
            int sib = lc;
            if (lc == index) { sib = rc; }
            
            map<string, int>::iterator itGroup;
            int pcountSize = 0;
            for (int j = 0; j < grouping.size(); j++) {
                map<string, int>::iterator itGroup = t->tree[sib].pcount.find(grouping[j]);
                if (itGroup != t->tree[sib].pcount.end()) { pcountSize++; if (pcountSize > 1) { break; } }
            }
            
            //if yes, I am not the root
            if (pcountSize != 0) {
                rootForGrouping.clear();
                rootForGrouping.insert(parent);
            }
            
            index = parent;
        }
        
        //get all nodes above the root to add so we don't add their u values above
        index = *(rootForGrouping.begin());
        while(t->tree[index].getParent() != -1){
            int parent = t->tree[index].getParent();
            rootForGrouping.insert(parent);
            
            index = parent;
        }
        
        
        return;
    }
    catch(exception& e) {
        m->errorOut(e, "Weighted", "getRoot");
        exit(1);
    }
}
/**************************************************************************************************/
double getLengthToRoot(MothurOut* m, Tree* t, int v, set<int> roots, vector<double>& leafToTreeRoot) {
    try {
        double sum = 0.0;
        int index = v;
        Utils util;
        
        //find length to complete tree root and save
        if (util.isEqual(leafToTreeRoot[v], 0.0)) {
            //you are a leaf
            if(!util.isEqual(t->tree[index].getBranchLength(), -1)){    leafToTreeRoot[v] += abs(t->tree[index].getBranchLength());    }
            
            index = t->tree[index].getParent();
                            
            //while you aren't at root
            while(t->tree[index].getParent() != -1){
                
                if (m->getControl_pressed()) {  return sum; }
                
                int parent = t->tree[index].getParent();
                                    
                if (!util.isEqual(t->tree[index].getBranchLength(), -1)) { leafToTreeRoot[v] += abs(t->tree[index].getBranchLength()); }
                
                index = parent;
            }
        }
        
        index = v;
        
        sum = leafToTreeRoot[v];
        
        //subtract excess root for this grouping
        for (set<int>::iterator it = roots.begin(); it != roots.end(); it++) {
            if (!util.isEqual(t->tree[*it].getBranchLength(), -1)) {
                sum -= abs(t->tree[*it].getBranchLength());
            }
        }
        
        return sum;
    }
    catch(exception& e) {
        m->errorOut(e, "Weighted", "getLengthToRoot");
        exit(1);
    }
}
/**************************************************************************************************/
int findNodeBelongingToThisComparison(MothurOut* m, vector<string>& namesOfGroupCombos, map< string, vector<int> >& groupNodeInfo) {
    try {
     
     int nodeBelonging = -1;
     for (int g = 0; g < namesOfGroupCombos.size(); g++) {
         if (groupNodeInfo[namesOfGroupCombos[g]].size() != 0) { nodeBelonging = groupNodeInfo[namesOfGroupCombos[g]][0]; break; }
     }
     
     //sanity check
     if (nodeBelonging == -1) {
         m->mothurOut("[WARNING]: cannot find a nodes in the tree from grouping ");
         for (int g = 0; g < namesOfGroupCombos.size()-1; g++) { m->mothurOut(namesOfGroupCombos[g] + "-"); }
         m->mothurOut(namesOfGroupCombos[namesOfGroupCombos.size()-1]);
         m->mothurOut(", skipping.\n");
     }
     
     return nodeBelonging;
 }
 catch(exception& e) {
     m->errorOut(e, "Weighted", "findNodeBelongingToThisComparison");
     exit(1);
 }
}
/**************************************************************************************************/
double findNumerator(MothurOut* m, Tree* t, set<int>& rootBranches, string groupA, string groupB, int groupACount, int groupBCount) {
    try {
        Utils util;
        double WScore = 0.0;
        
        for(int i=0;i<t->getNumNodes();i++){
            
            if (m->getControl_pressed()) { break; }
            
            double u = 0.00;
            
            //does this node have descendants from groupA
            map<string, int>::iterator it = t->tree[i].pcount.find(groupA);
            //if it does u = # of its descendants with a certain group / total number in tree with a certain group
            if (it != t->tree[i].pcount.end()) { u = (double) it->second / (double) groupACount; }
            
            //does this node have descendants from group l
            it = t->tree[i].pcount.find(groupB);
            
            //if it does subtract their percentage from u
            if (it != t->tree[i].pcount.end()) { u -= (double) it->second / (double) groupBCount; }
            
            if (!util.isEqual(t->tree[i].getBranchLength(), -1)) {
                //if this is not the root then add it
                if (rootBranches.count(i) == 0) {
                    u = abs(u * t->tree[i].getBranchLength());
                    WScore += u;
                }
            }
        }
        
        return WScore;
    }
    catch(exception& e) {
        m->errorOut(e, "Weighted", "findNumerator");
        exit(1);
    }
}
/**************************************************************************************************/
double findWeightedSums(MothurOut* m, Tree* t, set<int>& rootBranches, vector<double>& nodeToRootLength, string groupA, int groupACount) {
    try {
        double D = 0.0;
        
        //adding the wieghted sums from groupA
        for (int j = 0; j < t->groupNodeInfo[groupA].size(); j++) { //the leaf nodes that have seqs from groupA
            
            map<string, int>::iterator it = t->tree[t->groupNodeInfo[groupA][j]].pcount.find(groupA);
            int numSeqsInGroupI = it->second;

            double sum = getLengthToRoot(m, t, t->groupNodeInfo[groupA][j], rootBranches, nodeToRootLength);
            double weightedSum = ((numSeqsInGroupI * sum) / (double) groupACount);
            
            D += weightedSum;
        }
        
        return D;
    }
    catch(exception& e) {
        m->errorOut(e, "Weighted", "findNumerator");
        exit(1);
    }
}
/**************************************************************************************************/
void driverWeighted(weightedData* params) {
 try {
        Utils util;
		params->count = 0;
        vector<double> nodeToRootLength;  //length from leaf to tree root, grouping root maybe smaller. Used as reference and excess root is deducted if neccasary
        nodeToRootLength.resize(params->t->getNumLeaves(), 0.0); //set all leaf nodes length to root to zero

		for (int h = params->start; h < (params->start+params->num); h++) {
		
            if (params->m->getControl_pressed()) { break; }
		
			//initialize weighted score
			string groupA = params->namesOfGroupCombos[h][0];
			string groupB = params->namesOfGroupCombos[h][1];
            
            int groupACount = params->ct->getGroupCount(groupA);
            int groupBCount = params->ct->getGroupCount(groupB);
            
            set<int> rootBranches; //if not including root this will hold branches that are "above" the root for this comparison
			
            double WScore = 0.0;
			double D = 0.0;
            
            //find a node that belongs to one of the groups in this combo
            int nodeBelonging = findNodeBelongingToThisComparison(params->m, params->namesOfGroupCombos[h], params->t->groupNodeInfo);
            
            if (nodeBelonging != -1) {
            
                //fills rootBranches to exclude, if including the root then rootBranches should be empty.
                if (!params->includeRoot) { getRoot(params->m, params->t, nodeBelonging, params->namesOfGroupCombos[h], rootBranches); }

                WScore = findNumerator(params->m, params->t, rootBranches, groupA, groupB, groupACount, groupBCount);
                
                D += findWeightedSums(params->m, params->t, rootBranches, nodeToRootLength, groupA, groupACount);
                D += findWeightedSums(params->m, params->t, rootBranches, nodeToRootLength, groupB, groupBCount);
                
                double result = (WScore / D);
                if (isnan(result) || isinf(result)) { result = 0; }
                params->results.push_back(result);
                
            }else { params->results.push_back(0.0); }
            
			params->count++;
		}
	}
	catch(exception& e) {
		params->m->errorOut(e, "Weighted", "driverWeighted");
		exit(1);
	}
}
 /**************************************************************************************************/
EstOutput Weighted::getValues(Tree* t, string groupA, string groupB) { 
    try {
        EstOutput data;
        CountTable* ct = t->getCountTable();
        vector<double> nodeToRootLength;  //length from leaf to tree root, grouping root maybe smaller. Used as reference and excess root is deducted if neccesary
        nodeToRootLength.resize(t->getNumLeaves(), 0.0); //set all leaf nodes length to root to zero
        
        vector<string> grouping; grouping.push_back(groupA); grouping.push_back(groupB);
        
        int groupACount = ct->getGroupCount(groupA);
        int groupBCount = ct->getGroupCount(groupB);
        
        set<int> rootBranches; //if not including root this will hold branches that are "above" the root for this comparison
        
        double WScore = 0.0;
        double D = 0.0;
        
        //find a node that belongs to one of the groups in this combo
        int nodeBelonging = findNodeBelongingToThisComparison(m, grouping, t->groupNodeInfo);
        
        if (nodeBelonging != -1) {
            
            if (!includeRoot) { getRoot(m, t, nodeBelonging, grouping, rootBranches); }
            
            WScore = findNumerator(m, t, rootBranches, groupA, groupB, groupACount, groupBCount);
            
            D += findWeightedSums(m, t, rootBranches, nodeToRootLength, groupA, groupACount);
            D += findWeightedSums(m, t, rootBranches, nodeToRootLength, groupB, groupBCount);
            
            double result = (WScore / D);
            if (isnan(result) || isinf(result)) { result = 0; }
            data.push_back(result);
            
        }else { data.push_back(0.0); }
        
        return data;
    }
    catch(exception& e) {
        m->errorOut(e, "Weighted", "getValues");
        exit(1);
    }
}
/**************************************************************************************************/

EstOutput Weighted::createProcesses(Tree* t) {
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
        vector<std::thread*> workerThreads;
        vector<weightedData*> data;
        vector<string> Treenames; Treenames = t->getTreeNames();
        CountTable* ct = t->getCountTable();
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            CountTable* copyCount = new CountTable();
            copyCount->copy(ct);
            Tree* copyTree = new Tree(copyCount, Treenames);
            copyTree->getCopy(t);
            
            weightedData* dataBundle = new weightedData(lines[i+1].start, lines[i+1].end, namesOfGroupCombos, copyTree, copyCount, includeRoot);
            data.push_back(dataBundle);
            
            workerThreads.push_back(new std::thread(driverWeighted, dataBundle));
        }
        CountTable* copyCount = new CountTable();
        copyCount->copy(ct);
        Tree* copyTree = new Tree(copyCount, Treenames);
        copyTree->getCopy(t);
        
        weightedData* dataBundle = new weightedData(lines[0].start, lines[0].end, namesOfGroupCombos, copyTree, copyCount, includeRoot);
        driverWeighted(dataBundle);
        EstOutput results = dataBundle->results;
        delete copyTree; delete copyCount;
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

