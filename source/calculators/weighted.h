#ifndef WEIGHTED_H
#define WEIGHTED_H


/*
 *  weighted.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "treecalculator.h"
#include "counttable.h"

/***********************************************************************/

class Weighted : public TreeCalculator  {
	
	public:
        Weighted( bool r) : includeRoot(r) {};
		~Weighted() {};
		
		EstOutput getValues(Tree*, string, string);
		EstOutput getValues(Tree*, int, string);
		
	private:
		struct linePair {
			int start;
			int num;
			linePair(int i, int j) : start(i), num(j) {}
		};
		vector<linePair> lines;

		EstOutput data;
		map<string, int>::iterator it;
		map<string, double> WScore; //a score for each group combination i.e. AB, AC, BC.
		int processors;
		string outputDir;
		map< vector<string>, set<int> > rootForGrouping;  //maps a grouping combo to the root for that combo
		bool includeRoot;
		
		EstOutput driver(Tree*, vector< vector<string> >, int, int, CountTable*); 
		EstOutput createProcesses(Tree*, vector< vector<string> >, CountTable*);
		double getLengthToRoot(Tree*, int, string, string);
};

/***********************************************************************/
struct weightedData {
    int start;
	int num;
	MothurOut* m;
    EstOutput results;
    vector< vector<string> > namesOfGroupCombos;
    Tree* t;
    CountTable* ct;
    bool includeRoot;
    
	
	weightedData(){}
	weightedData(MothurOut* mout, int st, int en, vector< vector<string> > ngc, Tree* tree, CountTable* count, bool ir) {
        m = mout;
		start = st;
		num = en;
        namesOfGroupCombos = ngc;
        t = tree;
        ct = count;
        includeRoot = ir;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyWeightedThreadFunction(LPVOID lpParam){
	weightedData* pDataArray;
	pDataArray = (weightedData*)lpParam;
	try {
        map<string, int>::iterator it;
		vector<double> D;
		int count = 0;
        map< vector<string>, set<int> > rootForGrouping;
        map<string, double> WScore;
        
		for (int h = pDataArray->start; h < (pDataArray->start+pDataArray->num); h++) {
            
			if (pDataArray->m->getControl_pressed()) { return 0; }
            
			//initialize weighted score
			string groupA = pDataArray->namesOfGroupCombos[h][0];
			string groupB = pDataArray->namesOfGroupCombos[h][1];
			
			set<int> validBranches;
			WScore[groupA+groupB] = 0.0;
			D.push_back(0.0000); //initialize a spot in D for each combination
			
			//adding the wieghted sums from group i
			for (int j = 0; j < pDataArray->t->groupNodeInfo[groupA].size(); j++) { //the leaf nodes that have seqs from group i
				map<string, int>::iterator it = pDataArray->t->tree[pDataArray->t->groupNodeInfo[groupA][j]].pcount.find(groupA);
				int numSeqsInGroupI = it->second;

				//double sum = getLengthToRoot(pDataArray->t, pDataArray->t->groupNodeInfo[groupA][j], groupA, groupB);
                /*************************************************************************************/
                double sum = 0.0;
                int index = pDataArray->t->groupNodeInfo[groupA][j];
                
                //you are a leaf
                if(pDataArray->t->tree[index].getBranchLength() != -1){	sum += abs(pDataArray->t->tree[index].getBranchLength());	}
                double tempTotal = 0.0;
                index = pDataArray->t->tree[index].getParent();
                
                vector<string> grouping; grouping.push_back(groupA); grouping.push_back(groupB);
                
                rootForGrouping[grouping].insert(index);
                
                //while you aren't at root
                while(pDataArray->t->tree[index].getParent() != -1){
                    
                    if (pDataArray->m->getControl_pressed()) {  return 0; }
                    
                    int parent = pDataArray->t->tree[index].getParent();
                    
                    if (pDataArray->includeRoot) { //add everyone
                        if(pDataArray->t->tree[index].getBranchLength() != -1){	sum += abs(pDataArray->t->tree[index].getBranchLength());	}
                    }else {
                        
                        //am I the root for this grouping? if so I want to stop "early"
                        //does my sibling have descendants from the users groups?
                        int lc = pDataArray->t->tree[parent].getLChild();
                        int rc = pDataArray->t->tree[parent].getRChild();
                        
                        int sib = lc;
                        if (lc == index) { sib = rc; }
                        
                        map<string, int>::iterator itGroup;
                        int pcountSize = 0;
                        itGroup = pDataArray->t->tree[sib].pcount.find(groupA);
                        if (itGroup != pDataArray->t->tree[sib].pcount.end()) { pcountSize++;  }
                        itGroup = pDataArray->t->tree[sib].pcount.find(groupB);
                        if (itGroup != pDataArray->t->tree[sib].pcount.end()) { pcountSize++;  }
                        
                        //if yes, I am not the root so add me
                        if (pcountSize != 0) {
                            if (pDataArray->t->tree[index].getBranchLength() != -1) {
                                sum += abs(pDataArray->t->tree[index].getBranchLength()) + tempTotal;
                                tempTotal = 0.0;
                            }else {
                                sum += tempTotal;
                                tempTotal = 0.0;
                            }
                            rootForGrouping[grouping].clear();
                            rootForGrouping[grouping].insert(parent);
                        }else { //if no, I may be the root so add my br to tempTotal until I am proven innocent
                            if (pDataArray->t->tree[index].getBranchLength() != -1) {
                                tempTotal += abs(pDataArray->t->tree[index].getBranchLength()); 
                            }
                        }
                    }
                    
                    index = parent;	
                }
                
                //get all nodes above the root to add so we don't add their u values above
                index = *(rootForGrouping[grouping].begin());
                
                while(pDataArray->t->tree[index].getParent() != -1){
                    int parent = pDataArray->t->tree[index].getParent();
                    rootForGrouping[grouping].insert(parent);
                    index = parent;
                }
                
                /*************************************************************************************/
				double weightedSum = ((numSeqsInGroupI * sum) / (double)pDataArray->ct->getGroupCount(groupA));
                
				D[count] += weightedSum;
			}
			
			//adding the wieghted sums from group l
			for (int j = 0; j < pDataArray->t->groupNodeInfo[groupB].size(); j++) { //the leaf nodes that have seqs from group l
				map<string, int>::iterator it = pDataArray->t->tree[pDataArray->t->groupNodeInfo[groupB][j]].pcount.find(groupB);
				int numSeqsInGroupL = it->second;
				
				//double sum = getLengthToRoot(pDataArray->t, pDataArray->t->groupNodeInfo[groupB][j], groupA, groupB);
                /*************************************************************************************/
                double sum = 0.0;
                int index = pDataArray->t->groupNodeInfo[groupB][j];
                
                //you are a leaf
                if(pDataArray->t->tree[index].getBranchLength() != -1){	sum += abs(pDataArray->t->tree[index].getBranchLength());	}
                double tempTotal = 0.0;
                index = pDataArray->t->tree[index].getParent();
                
                vector<string> grouping; grouping.push_back(groupA); grouping.push_back(groupB);
                
                rootForGrouping[grouping].insert(index);
                
                //while you aren't at root
                while(pDataArray->t->tree[index].getParent() != -1){
                    
                    if (pDataArray->m->getControl_pressed()) {  return 0; }
                    
                    int parent = pDataArray->t->tree[index].getParent();
                    
                    if (pDataArray->includeRoot) { //add everyone
                        if(pDataArray->t->tree[index].getBranchLength() != -1){	sum += abs(pDataArray->t->tree[index].getBranchLength());	}
                    }else {
                        
                        //am I the root for this grouping? if so I want to stop "early"
                        //does my sibling have descendants from the users groups?
                        int lc = pDataArray->t->tree[parent].getLChild();
                        int rc = pDataArray->t->tree[parent].getRChild();
                        
                        int sib = lc;
                        if (lc == index) { sib = rc; }
                        
                        map<string, int>::iterator itGroup;
                        int pcountSize = 0;
                        itGroup = pDataArray->t->tree[sib].pcount.find(groupA);
                        if (itGroup != pDataArray->t->tree[sib].pcount.end()) { pcountSize++;  }
                        itGroup = pDataArray->t->tree[sib].pcount.find(groupB);
                        if (itGroup != pDataArray->t->tree[sib].pcount.end()) { pcountSize++;  }
                        
                        //if yes, I am not the root so add me
                        if (pcountSize != 0) {
                            if (pDataArray->t->tree[index].getBranchLength() != -1) {
                                sum += abs(pDataArray->t->tree[index].getBranchLength()) + tempTotal;
                                tempTotal = 0.0;
                            }else {
                                sum += tempTotal;
                                tempTotal = 0.0;
                            }
                            rootForGrouping[grouping].clear();
                            rootForGrouping[grouping].insert(parent);
                        }else { //if no, I may be the root so add my br to tempTotal until I am proven innocent
                            if (pDataArray->t->tree[index].getBranchLength() != -1) {
                                tempTotal += abs(pDataArray->t->tree[index].getBranchLength());
                            }
                        }
                    }
                    
                    index = parent;
                }
                
                //get all nodes above the root to add so we don't add their u values above
                index = *(rootForGrouping[grouping].begin());
                
                while(pDataArray->t->tree[index].getParent() != -1){
                    int parent = pDataArray->t->tree[index].getParent();
                    rootForGrouping[grouping].insert(parent);
                    index = parent;
                }
                
                /*************************************************************************************/

				double weightedSum = ((numSeqsInGroupL * sum) / (double)pDataArray->ct->getGroupCount(groupB));
                
				D[count] += weightedSum;
			}
			count++;
		}
        
		//calculate u for the group comb
		for (int h = pDataArray->start; h < (pDataArray->start+pDataArray->num); h++) {
			//report progress
			//pDataArray->m->mothurOut("Processing combo: " + toString(h)); pDataArray->m->mothurOutEndLine();
            
			string groupA = pDataArray->namesOfGroupCombos[h][0];
			string groupB = pDataArray->namesOfGroupCombos[h][1];
			
			//calculate u for the group comb
			for(int i=0;i<pDataArray->t->getNumNodes();i++){
				
				if (pDataArray->m->getControl_pressed()) { return 0; }
				
				double u;
				//int pcountSize = 0;
				//does this node have descendants from groupA
				it = pDataArray->t->tree[i].pcount.find(groupA);
				//if it does u = # of its descendants with a certain group / total number in tree with a certain group
				if (it != pDataArray->t->tree[i].pcount.end()) {
					u = (double) pDataArray->t->tree[i].pcount[groupA] / (double) pDataArray->ct->getGroupCount(groupA);
				}else { u = 0.00; }
				
				
				//does this node have descendants from group l
				it = pDataArray->t->tree[i].pcount.find(groupB);
				
				//if it does subtract their percentage from u
				if (it != pDataArray->t->tree[i].pcount.end()) {
					u -= (double) pDataArray->t->tree[i].pcount[groupB] / (double) pDataArray->ct->getGroupCount(groupB);
				}
				
				if (pDataArray->includeRoot) {
					if (pDataArray->t->tree[i].getBranchLength() != -1) {
						u = abs(u * pDataArray->t->tree[i].getBranchLength());
						WScore[(groupA+groupB)] += u;
					}
				}else {
					//if this is not the root then add it
					if (rootForGrouping[pDataArray->namesOfGroupCombos[h]].count(i) == 0) {
						if (pDataArray->t->tree[i].getBranchLength() != -1) {
							u = abs(u * pDataArray->t->tree[i].getBranchLength());
							WScore[(groupA+groupB)] += u;
						}
					}
				}
			}
			
		}
		
		/********************************************************/
		//calculate weighted score for the group combination
		double UN;
		count = 0;
		for (int h = pDataArray->start; h < (pDataArray->start+pDataArray->num); h++) {
			UN = (WScore[pDataArray->namesOfGroupCombos[h][0]+pDataArray->namesOfGroupCombos[h][1]] / D[count]);
			if (isnan(UN) || isinf(UN)) { UN = 0; }
			pDataArray->results.push_back(UN);
			count++;
		}
        
		return 0;

    }
	catch(exception& e) {
		pDataArray->m->errorOut(e, "Weighted", "MyWeightedThreadFunction");
		exit(1);
	}
}
#endif


#endif
