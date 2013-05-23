#ifndef PARSIMONY_H
#define PARSIMONY_H


/*
 *  parsimony.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/26/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "treecalculator.h"
#include "counttable.h"

/***********************************************************************/

class Parsimony : public TreeCalculator  {
	
	public:
		Parsimony() {};
		~Parsimony() {};
		EstOutput getValues(Tree*, int, string);
		
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
	
		EstOutput driver(Tree*, vector< vector<string> >, int, int, CountTable*); 
		EstOutput createProcesses(Tree*, vector< vector<string> >, CountTable*);
};
/***********************************************************************/
struct parsData {
    int start;
	int num;
	MothurOut* m;
    EstOutput results;
    vector< vector<string> > namesOfGroupCombos;
    Tree* t;
    CountTable* ct;
    
	parsData(){}
	parsData(MothurOut* mout, int st, int en, vector< vector<string> > ngc, Tree* tree, CountTable* count) {
        m = mout;
		start = st;
		num = en;
        namesOfGroupCombos = ngc;
        t = tree;
        ct = count;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyParsimonyThreadFunction(LPVOID lpParam){
	parsData* pDataArray;
	pDataArray = (parsData*)lpParam;
	try {
        
        pDataArray->results.resize(pDataArray->num);
		
		Tree* copyTree = new Tree(pDataArray->ct);
		int count = 0;
		
		for (int h = pDataArray->start; h < (pDataArray->start+pDataArray->num); h++) {
            
			if (pDataArray->m->control_pressed) { delete copyTree; return 0; }
            
			int score = 0;
			
			//groups in this combo
			vector<string> groups = pDataArray->namesOfGroupCombos[h];
			
			//copy users tree so that you can redo pgroups
			copyTree->getCopy(pDataArray->t);
			
			//create pgroups that reflect the groups the user want to use
			for(int i=copyTree->getNumLeaves();i<copyTree->getNumNodes();i++){
				copyTree->tree[i].pGroups = (copyTree->mergeUserGroups(i, groups));
			}
			
			for(int i=copyTree->getNumLeaves();i<copyTree->getNumNodes();i++){
				
				if (pDataArray->m->control_pressed) { return 0; }
				
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
			
			pDataArray->results[count] = score;
			count++;
		}
        
		delete copyTree;
        
        return 0;
        
    }
	catch(exception& e) {
		pDataArray->m->errorOut(e, "Parsimony", "MyParsimonyThreadFunction");
		exit(1);
	}
}
#endif

#endif
