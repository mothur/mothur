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
Parsimony::Parsimony(vector<string> G) : Groups(G) {
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
        m->errorOut(e, "Parsimony", "Parsimony");
        exit(1);
    }
}
/**************************************************************************************************/

EstOutput Parsimony::getValues(Tree* t, int p, string o) {
	try {
		processors = p; outputDir = o;
        CountTable* ct = t->getCountTable();
        Treenames = t->getTreeNames();
		
		return (createProcesses(t, ct));
		
	}
	catch(exception& e) {
		m->errorOut(e, "Parsimony", "getValues");
		exit(1);
	}
}
/**************************************************************************************************/
void driverPars(parsData* params) {
    try {
        
        Tree copyTree(params->ct, params->Treenames);
        int count = 0;
        
        for (int h = params->start; h < (params->start+params->num); h++) {
            
            if (params->m->getControl_pressed()) { break; }
            
            int score = 0;
            
            //groups in this combo
            vector<string> groups = params->namesOfGroupCombos[h];
            
            //copy users tree so that you can redo pgroups
            copyTree.getCopy(params->t);
            
            //create pgroups that reflect the groups the user want to use
            for(int i=copyTree.getNumLeaves();i<copyTree.getNumNodes();i++){
                copyTree.tree[i].pGroups = (copyTree.mergeUserGroups(i, groups));
            }
            
            for(int i=copyTree.getNumLeaves();i<copyTree.getNumNodes();i++){
                
                if (params->m->getControl_pressed()) { break; }
                
                int lc = copyTree.tree[i].getLChild();
                int rc = copyTree.tree[i].getRChild();
                
                int iSize = copyTree.tree[i].pGroups.size();
                int rcSize = copyTree.tree[rc].pGroups.size();
                int lcSize = copyTree.tree[lc].pGroups.size();
                
                //if isize are 0 then that branch is to be ignored
                if (iSize == 0) { }
                else if ((rcSize == 0) || (lcSize == 0)) { }
                //if you have more groups than either of your kids then theres been a change.
                else if(iSize > rcSize || iSize > lcSize){
                    score++;
                }
            } 
            
            params->results[count] = score;
            count++;
        }
    }
    catch(exception& e) {
        params->m->errorOut(e, "Parsimony", "driver");
        exit(1);
    }
}
/**************************************************************************************************/

EstOutput Parsimony::createProcesses(Tree* t, CountTable* ct) {
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
        vector<parsData*> data;
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            CountTable* copyCount = new CountTable();
            copyCount->copy(ct);
            Tree* copyTree = new Tree(copyCount, Treenames);
            copyTree->getCopy(t);
            
            parsData* dataBundle = new parsData(lines[i+1].start, lines[i+1].end, namesOfGroupCombos, copyTree, copyCount);
            data.push_back(dataBundle);
            
            workerThreads.push_back(new thread(driverPars, dataBundle));
        }
        
        parsData* dataBundle = new parsData(lines[0].start, lines[0].end, namesOfGroupCombos, t, ct);
        driverPars(dataBundle);
        EstOutput results = dataBundle->results;
        delete dataBundle;
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            
            for (int j = 0; j < data[i]->results.size(); j++) {  results.push_back(data[i]->results[j]);  }
            
            delete data[i]->t;
            delete data[i]->ct;
            delete data[i];
            delete workerThreads[i];
        }
		
        return results;
	}
	catch(exception& e) {
		m->errorOut(e, "Parsimony", "createProcesses");
		exit(1);
	}
}
/**************************************************************************************************/

