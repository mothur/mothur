#ifndef UNWEIGHTED_H
#define UNWEIGHTED_H


/*
 *  unweighted.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "treecalculator.h"
#include "counttable.h"

/***********************************************************************/

class Unweighted : public TreeCalculator  {
	
	public:
        Unweighted(bool r) : includeRoot(r) {};
		~Unweighted() {};
		EstOutput getValues(Tree*, int, string);
		EstOutput getValues(Tree*, string, string, int, string);
		
	private:
		struct linePair {
			int start;
			int num;
			linePair(int i, int j) : start(i), num(j) {}
		};
		vector<linePair> lines;
		
		EstOutput data;
		int processors;
		string outputDir;
		map< vector<string>, set<int> > rootForGrouping;  //maps a grouping combo to the roots for that combo
		bool includeRoot;
		
		EstOutput driver(Tree*, vector< vector<string> >, int, int, CountTable*); 
		EstOutput createProcesses(Tree*, vector< vector<string> >, CountTable*);
		EstOutput driver(Tree*, vector< vector<string> >, int, int, bool, CountTable*); 
		EstOutput createProcesses(Tree*, vector< vector<string> >, bool, CountTable*);
		int getRoot(Tree*, int, vector<string>);
};

/***********************************************************************/
struct unweightedData {
    int start;
	int num;
	MothurOut* m;
    EstOutput results;
    vector< vector<string> > namesOfGroupCombos;
    Tree* t;
    CountTable* ct;
    bool includeRoot;
    
	unweightedData(){}
	unweightedData(MothurOut* mout, int st, int en, vector< vector<string> > ngc, Tree* tree, CountTable* count, bool ir) {
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
static DWORD WINAPI MyUnWeightedThreadFunction(LPVOID lpParam){
	unweightedData* pDataArray;
	pDataArray = (unweightedData*)lpParam;
	try {
        pDataArray->results.resize(pDataArray->num);
        map< vector<string>, set<int> > rootForGrouping;
		
		int count = 0;
		        
		for (int h = pDataArray->start; h < (pDataArray->start+pDataArray->num); h++) {
			if (pDataArray->m->control_pressed) { return 0; }
            
			double UniqueBL=0.0000;  //a branch length is unique if it's chidren are from the same group
			double totalBL = 0.00;	//all branch lengths
			double UW = 0.00;		//Unweighted Value = UniqueBL / totalBL;
            
			//find a node that belongs to one of the groups in this combo
			int nodeBelonging = -1;
			for (int g = 0; g < pDataArray->namesOfGroupCombos[h].size(); g++) {
				if (pDataArray->t->groupNodeInfo[pDataArray->namesOfGroupCombos[h][g]].size() != 0) { nodeBelonging = pDataArray->t->groupNodeInfo[pDataArray->namesOfGroupCombos[h][g]][0]; break; }
			}
			
			//sanity check
			if (nodeBelonging == -1) {
				pDataArray->m->mothurOut("[WARNING]: cannot find a nodes in the tree from grouping ");
				for (int g = 0; g < pDataArray->namesOfGroupCombos[h].size()-1; g++) { pDataArray->m->mothurOut(pDataArray->namesOfGroupCombos[h][g] + "-"); }
				pDataArray->m->mothurOut(pDataArray->namesOfGroupCombos[h][pDataArray->namesOfGroupCombos[h].size()-1]);
				pDataArray->m->mothurOut(", skipping."); pDataArray->m->mothurOutEndLine(); pDataArray->results[count] = UW;
			}else{
				
				//if including the root this clears rootForGrouping[namesOfGroupCombos[h]]
				//getRoot(t, nodeBelonging, namesOfGroupCombos[h]);
				/////////////////////////////////////////////////////////////////////////////
                //you are a leaf so get your parent
                vector<string> grouping = pDataArray->namesOfGroupCombos[h];
                int index = pDataArray->t->tree[nodeBelonging].getParent();
                
                if (pDataArray->includeRoot) {
                    rootForGrouping[grouping].clear();
                }else {
                    
                    //my parent is a potential root
                    rootForGrouping[grouping].insert(index);
                    
                    //while you aren't at root
                    while(pDataArray->t->tree[index].getParent() != -1){
                        //cout << index << endl;
                        if (pDataArray->m->control_pressed) {  return 0; }
                        
                        //am I the root for this grouping? if so I want to stop "early"
                        //does my sibling have descendants from the users groups?
                        //if so I am not the root
                        int parent = pDataArray->t->tree[index].getParent();
                        int lc = pDataArray->t->tree[parent].getLChild();
                        int rc = pDataArray->t->tree[parent].getRChild();
                        
                        int sib = lc;
                        if (lc == index) { sib = rc; }
                        
                        map<string, int>::iterator itGroup;
                        int pcountSize = 0;
                        for (int j = 0; j < grouping.size(); j++) {
                            map<string, int>::iterator itGroup = pDataArray->t->tree[sib].pcount.find(grouping[j]);
                            if (itGroup != pDataArray->t->tree[sib].pcount.end()) { pcountSize++; if (pcountSize > 1) { break; } }
                        }
                        
                        //if yes, I am not the root
                        if (pcountSize != 0) {
                            rootForGrouping[grouping].clear();
                            rootForGrouping[grouping].insert(parent);
                        }
                        
                        index = parent;	
                    }
                    
                    //get all nodes above the root to add so we don't add their u values above
                    index = *(rootForGrouping[grouping].begin());
                    while(pDataArray->t->tree[index].getParent() != -1){
                        int parent = pDataArray->t->tree[index].getParent();
                        rootForGrouping[grouping].insert(parent);
                        //cout << parent << " in root" << endl;
                        index = parent;
                    }
                }
                /////////////////////////////////////////////////////////////////////////////
                
				for(int i=0;i<pDataArray->t->getNumNodes();i++){
					
					if (pDataArray->m->control_pressed) {  return 0; }
					//cout << i << endl;
					//pcountSize = 0, they are from a branch that is entirely from a group the user doesn't want
					//pcountSize = 2, not unique to one group
					//pcountSize = 1, unique to one group
					
					int pcountSize = 0;
					for (int j = 0; j < pDataArray->namesOfGroupCombos[h].size(); j++) {
						map<string, int>::iterator itGroup = pDataArray->t->tree[i].pcount.find(pDataArray->namesOfGroupCombos[h][j]);
						if (itGroup != pDataArray->t->tree[i].pcount.end()) { pcountSize++; if (pcountSize > 1) { break; } }
					}
					
					
					//unique calc
					if (pcountSize == 0) { }
					else if ((pDataArray->t->tree[i].getBranchLength() != -1) && (pcountSize == 1) && (rootForGrouping[pDataArray->namesOfGroupCombos[h]].count(i) == 0)) { //you have a unique branch length and you are not the root
						UniqueBL += abs(pDataArray->t->tree[i].getBranchLength());
					}
                    
					//total calc
					if (pcountSize == 0) { }
					else if ((pDataArray->t->tree[i].getBranchLength() != -1) && (pcountSize != 0) && (rootForGrouping[pDataArray->namesOfGroupCombos[h]].count(i) == 0)) { //you have a branch length and you are not the root
						totalBL += abs(pDataArray->t->tree[i].getBranchLength());
					}
				}
                //cout << UniqueBL << '\t' << totalBL << endl;
				UW = (UniqueBL / totalBL);
                
				if (isnan(UW) || isinf(UW)) { UW = 0; }
                
				pDataArray->results[count] = UW;
			}
			count++;
		}
		
		return 0;
        
    }
	catch(exception& e) {
		pDataArray->m->errorOut(e, "UnWeighted", "MyUnWeightedThreadFunction");
		exit(1);
	}
}
/**************************************************************************************************/

static DWORD WINAPI MyUnWeightedRandomThreadFunction(LPVOID lpParam){
	unweightedData* pDataArray;
	pDataArray = (unweightedData*)lpParam;
	try {
        pDataArray->results.resize(pDataArray->num);
		
		int count = 0;
		
		Tree* copyTree = new Tree(pDataArray->ct);
		
		for (int h = pDataArray->start; h < (pDataArray->start+pDataArray->num); h++) {
            
			if (pDataArray->m->control_pressed) { return 0; }
            
            map< vector<string>, set<int> > rootForGrouping;
            
			//copy random tree passed in
			copyTree->getCopy(pDataArray->t);
            
			//swap labels in the groups you want to compare
			copyTree->assembleRandomUnifracTree(pDataArray->namesOfGroupCombos[h]);
			
			double UniqueBL=0.0000;  //a branch length is unique if it's chidren are from the same group
			double totalBL = 0.00;	//all branch lengths
			double UW = 0.00;		//Unweighted Value = UniqueBL / totalBL;
			//find a node that belongs to one of the groups in this combo
			int nodeBelonging = -1;
			for (int g = 0; g < pDataArray->namesOfGroupCombos[h].size(); g++) {
				if (copyTree->groupNodeInfo[pDataArray->namesOfGroupCombos[h][g]].size() != 0) { nodeBelonging = copyTree->groupNodeInfo[pDataArray->namesOfGroupCombos[h][g]][0]; break; }
			}
			
			//sanity check
			if (nodeBelonging == -1) {
				pDataArray->m->mothurOut("[WARNING]: cannot find a nodes in the tree from grouping ");
				for (int g = 0; g < pDataArray->namesOfGroupCombos[h].size()-1; g++) { pDataArray->m->mothurOut(pDataArray->namesOfGroupCombos[h][g] + "-"); }
				pDataArray->m->mothurOut(pDataArray->namesOfGroupCombos[h][pDataArray->namesOfGroupCombos[h].size()-1]);
				pDataArray->m->mothurOut(", skipping."); pDataArray->m->mothurOutEndLine(); pDataArray->results[count] = UW;
			}else{
				
				//if including the root this clears rootForGrouping[namesOfGroupCombos[h]]
				//getRoot(copyTree, nodeBelonging, namesOfGroupCombos[h]);
                /////////////////////////////////////////////////////////////////////////////
                //you are a leaf so get your parent
                vector<string> grouping = pDataArray->namesOfGroupCombos[h];
                int index = copyTree->tree[nodeBelonging].getParent();
                
                if (pDataArray->includeRoot) {
                    rootForGrouping[grouping].clear();
                }else {
                    
                    //my parent is a potential root
                    rootForGrouping[grouping].insert(index);
                    
                    //while you aren't at root
                    while(copyTree->tree[index].getParent() != -1){
                        //cout << index << endl;
                        if (pDataArray->m->control_pressed) {  return 0; }
                        
                        //am I the root for this grouping? if so I want to stop "early"
                        //does my sibling have descendants from the users groups?
                        //if so I am not the root
                        int parent = copyTree->tree[index].getParent();
                        int lc = copyTree->tree[parent].getLChild();
                        int rc = copyTree->tree[parent].getRChild();
                        
                        int sib = lc;
                        if (lc == index) { sib = rc; }
                        
                        map<string, int>::iterator itGroup;
                        int pcountSize = 0;
                        for (int j = 0; j < grouping.size(); j++) {
                            map<string, int>::iterator itGroup = copyTree->tree[sib].pcount.find(grouping[j]);
                            if (itGroup != copyTree->tree[sib].pcount.end()) { pcountSize++; if (pcountSize > 1) { break; } }
                        }
                        
                        //if yes, I am not the root
                        if (pcountSize != 0) {
                            rootForGrouping[grouping].clear();
                            rootForGrouping[grouping].insert(parent);
                        }
                        
                        index = parent;
                    }
                    
                    //get all nodes above the root to add so we don't add their u values above
                    index = *(rootForGrouping[grouping].begin());
                    while(copyTree->tree[index].getParent() != -1){
                        int parent = copyTree->tree[index].getParent();
                        rootForGrouping[grouping].insert(parent);
                        //cout << parent << " in root" << endl;
                        index = parent;
                    }
                }
                /////////////////////////////////////////////////////////////////////////////
				for(int i=0;i<copyTree->getNumNodes();i++){
					
					if (pDataArray->m->control_pressed) {  return 0; }
					
					//pcountSize = 0, they are from a branch that is entirely from a group the user doesn't want
					//pcountSize = 2, not unique to one group
					//pcountSize = 1, unique to one group
                    int pcountSize = 0;
					for (int j = 0; j < pDataArray->namesOfGroupCombos[h].size(); j++) {
						map<string, int>::iterator itGroup = copyTree->tree[i].pcount.find(pDataArray->namesOfGroupCombos[h][j]);
						if (itGroup != copyTree->tree[i].pcount.end()) { pcountSize++; if (pcountSize > 1) { break; } }
					}
					
					//unique calc
					if (pcountSize == 0) { }
					else if ((copyTree->tree[i].getBranchLength() != -1) && (pcountSize == 1) && (rootForGrouping[pDataArray->namesOfGroupCombos[h]].count(i) == 0)) { //you have a unique branch length and you are not the root
						UniqueBL += abs(copyTree->tree[i].getBranchLength());
					}
					
					//total calc
					if (pcountSize == 0) { }
					else if ((copyTree->tree[i].getBranchLength() != -1) && (pcountSize != 0) && (rootForGrouping[pDataArray->namesOfGroupCombos[h]].count(i) == 0)) { //you have a branch length and you are not the root
						totalBL += abs(copyTree->tree[i].getBranchLength());
					}
					
				}
				cout << h << '\t' << UniqueBL << '\t' << totalBL << endl;
				UW = (UniqueBL / totalBL);
				
				if (isnan(UW) || isinf(UW)) { UW = 0; }
				
				pDataArray->results[count] = UW;
                cout << h << '\t' << UW << endl;
			}
			count++;
			
		}
		
		delete copyTree;
		
		return 0;
    }
	catch(exception& e) {
		pDataArray->m->errorOut(e, "UnWeighted", "MyUnWeightedRandomThreadFunction");
		exit(1);
	}
}

#endif


#endif
