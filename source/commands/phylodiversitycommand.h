#ifndef PHYLODIVERSITYCOMMAND_H
#define PHYLODIVERSITYCOMMAND_H

/*
 *  phylodiversitycommand.h
 *  Mothur
 *
 *  Created by westcott on 4/30/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "counttable.h"

#include "tree.h"


class PhyloDiversityCommand : public Command {
	
	public:
		PhyloDiversityCommand(string);
		PhyloDiversityCommand();
		~PhyloDiversityCommand(){}
	
		vector<string> setParameters();
		string getCommandName()			{ return "phylo.diversity";			}
		string getCommandCategory()		{ return "Hypothesis Testing";		}
		
	string getHelpString();	
    string getOutputPattern(string);	
		string getCitation() { return "Faith DP (1994). Phylogenetic pattern and the quantification of organismal biodiversity. Philos Trans R Soc Lond B Biol Sci 345: 45-58. \nhttp://www.mothur.org/wiki/Phylo.diversity"; }
		string getDescription()		{ return "phylo.diversity"; }

		int execute();
		void help() { m->mothurOut(getHelpString()); }
private:
		CountTable* ct;
		float freq;
		int iters, processors, numUniquesInName, subsampleSize;
		bool abort, rarefy, summary, collect, scale, subsample;
		string groups, outputDir, treefile, groupfile, namefile, countfile;
		vector<string> Groups, outputNames; //holds groups to be used, and outputFile names
		
        map<string, int> getRootForGroups(Tree* t);
		int readNamesFile();
		void printData(set<int>&, map< string, vector<float> >&, ofstream&, int);
		void printSumData(map< string, vector<float> >&, ofstream&, int);
        vector<float> calcBranchLength(Tree*, int, vector< map<string, bool> >&, map<string, int>);
		int driver(Tree*, map< string, vector<float> >&, map<string, vector<float> >&, int, int, vector<int>&, set<int>&, ofstream&, ofstream&, bool);
		int createProcesses(vector<int>&, Tree*, map< string, vector<float> >&, map<string, vector<float> >&, int, int, vector<int>&, set<int>&, ofstream&, ofstream&);

};

/***********************************************************************/
struct phylodivData {
    int numIters;
	MothurOut* m;
    map< string, vector<float> > div;
    map<string, vector<float> > sumDiv;
    map<string, int> rootForGroup;
    vector<int> randomLeaf;
    set<int> numSampledList;
    int increment, subsampleSize;
    Tree* t;
    CountTable* ct;
    bool includeRoot, subsample;
	
   
	phylodivData(){}
	phylodivData(MothurOut* mout, int ni,  map< string, vector<float> > cd, map< string, vector<float> > csd, Tree* tree, CountTable* count, int incre, vector<int> crl, set<int> nsl, map<string, int> rfg, bool su, int suS) {
        m = mout;
        t = tree;
        ct = count;
        div = cd;
        numIters = ni;
        sumDiv = csd;
        increment = incre;
        randomLeaf = crl;
        numSampledList = nsl;
        rootForGroup = rfg;
        subsample = su;
        subsampleSize = suS;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyPhyloDivThreadFunction(LPVOID lpParam){
	phylodivData* pDataArray;
	pDataArray = (phylodivData*)lpParam;
	try {
        int numLeafNodes = pDataArray->randomLeaf.size();
		vector<string> mGroups = pDataArray->m->getGroups();
        
		for (int l = 0; l < pDataArray->numIters; l++) {
            pDataArray->util.mothurRandomShuffle(pDataArray->randomLeaf);
            
            //initialize counts
            map<string, int> counts;
            vector< map<string, bool> > countedBranch;
            for (int i = 0; i < pDataArray->t->getNumNodes(); i++) {
                map<string, bool> temp;
                for (int j = 0; j < mGroups.size(); j++) { temp[mGroups[j]] = false; }
                countedBranch.push_back(temp);
            }
            
            for (int j = 0; j < mGroups.size(); j++) {  counts[mGroups[j]] = 0;   }
            
            map<string, int> metCount; bool allDone = false;
            for (int j = 0; j < mGroups.size(); j++) {  counts[mGroups[j]] = false;   }
            for(int k = 0; k < numLeafNodes; k++){
                
                if (pDataArray->m->getControl_pressed()) { return 0; }
                
                //calc branch length of randomLeaf k
                //vector<float> br = calcBranchLength(t, randomLeaf[k], countedBranch, rootForGroup);
                //(Tree* t, int leaf, vector< map<string, bool> >& counted, map<string, int> roots
                /////////////////////////////////////////////////////////////////////////////////////
                vector<float> br;
                int index = pDataArray->randomLeaf[k];
                
                vector<string> groups = pDataArray->t->tree[pDataArray->randomLeaf[k]].getGroup();
                br.resize(groups.size(), 0.0);

                //you are a leaf
                if(pDataArray->t->tree[index].getBranchLength() != -1){
                    for (int k = 0; k < groups.size(); k++) {
                        br[k] += abs(pDataArray->t->tree[index].getBranchLength());
                    }
                }

                index = pDataArray->t->tree[index].getParent();
                
                //while you aren't at root
                while(pDataArray->t->tree[index].getParent() != -1){
                    
                    if (pDataArray->m->getControl_pressed()) {  return 0; }
                    
                    for (int k = 0; k < groups.size(); k++) {
                        
                        if (index >= pDataArray->rootForGroup[groups[k]]) { countedBranch[index][groups[k]] = true; } //if you are at this groups "root", then say we are done
                        
                        if (!countedBranch[index][groups[k]]){ //if counted[index][groups[k] is true this groups has already added all br from here to root, so quit early
                            if (pDataArray->t->tree[index].getBranchLength() != -1) {
                                br[k] += abs(pDataArray->t->tree[index].getBranchLength());
                            }
                            countedBranch[index][groups[k]] = true;
                        }
                    }
                    index = pDataArray->t->tree[index].getParent();	
                }
                /////////////////////////////////////////////////////////////////////////////////////
                
                //for each group in the groups update the total branch length accounting for the names file
                groups = pDataArray->t->tree[pDataArray->randomLeaf[k]].getGroup();
                
                for (int j = 0; j < groups.size(); j++) {
                    
                    if (pDataArray->util.inUsersGroups(groups[j], mGroups)) {
                        int numSeqsInGroupJ = 0;
                        map<string, int>::iterator it;
                        it = pDataArray->t->tree[pDataArray->randomLeaf[k]].pcount.find(groups[j]);
                        if (it != pDataArray->t->tree[pDataArray->randomLeaf[k]].pcount.end()) { //this leaf node contains seqs from group j
                            numSeqsInGroupJ = it->second;
                        }
                        
                        if (numSeqsInGroupJ != 0) {	pDataArray->div[groups[j]][(counts[groups[j]]+1)] = pDataArray->div[groups[j]][counts[groups[j]]] + br[j];  }
                        
                        for (int s = (counts[groups[j]]+2); s <= (counts[groups[j]]+numSeqsInGroupJ); s++) {
                            pDataArray->div[groups[j]][s] = pDataArray->div[groups[j]][s-1];  //update counts, but don't add in redundant branch lengths
                        }
                        counts[groups[j]] += numSeqsInGroupJ;
                        if (pDataArray->subsample) {
                            if (counts[groups[j]] >= pDataArray->subsampleSize) { metCount[groups[j]] = true; }
                            bool allTrue = true;
                            for (int h = 0; h < mGroups.size(); h++) {
                                if (!metCount[mGroups[h]]) { allTrue = false; }
                            }
                            if (allTrue) { allDone = true; }
                        }
                        if (allDone) { j+=groups.size(); k+=numLeafNodes; }
                    }
                }
            }
            
            
            //add this diversity to the sum
            for (int j = 0; j < mGroups.size(); j++) {
                for (int g = 0; g < pDataArray->div[mGroups[j]].size(); g++) {
                    pDataArray->sumDiv[mGroups[j]][g] += pDataArray->div[mGroups[j]][g];
                }
            }
            
        }

        return 0;
    }
	catch(exception& e) {
		pDataArray->m->errorOut(e, "PhyloDiversityCommand", "MyPhyloDivThreadFunction");
		exit(1);
	}
}
#endif

#endif

