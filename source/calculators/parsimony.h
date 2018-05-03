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
        Parsimony(vector<string> G);
		~Parsimony() {};
		EstOutput getValues(Tree*, int, string);
		
	private:
        vector<string> Groups, Treenames;
		int processors;
		string outputDir;
        Utils util;
        vector< vector<string> > namesOfGroupCombos;
	
		EstOutput createProcesses(Tree*, CountTable*);
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
    Utils util;
    vector<string> Treenames;
    
	parsData(){}
	parsData(int st, int en, vector< vector<string> > ngc, Tree* tree, CountTable* count) {
        m = MothurOut::getInstance();
		start = st;
		num = en;
        namesOfGroupCombos = ngc;
        t = tree;
        ct = count;
        Treenames = t->getTreeNames();
        results.resize(num);
	}
};

/**************************************************************************************************/
#if defined NON_WINDOWS
#else
static DWORD WINAPI MyParsimonyThreadFunction(LPVOID lpParam){
	parsData* pDataArray;
	pDataArray = (parsData*)lpParam;
	try {
        
        pDataArray->results.resize(pDataArray->num);
		
		Tree* copyTree = new Tree(pDataArray->ct, pDataArray->Treenames);
		int count = 0;
		
		for (int h = pDataArray->start; h < (pDataArray->start+pDataArray->num); h++) {
            
			if (pDataArray->m->getControl_pressed()) { delete copyTree; return 0; }
            
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
				
				if (pDataArray->m->getControl_pressed()) { return 0; }
				
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
