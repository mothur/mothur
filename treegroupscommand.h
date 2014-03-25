#ifndef TREEGROUPCOMMAND_H
#define TREEGROUPCOMMAND_H

/*
 *  treegroupscommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
 
#include "command.hpp"
#include "inputdata.h"
#include "groupmap.h"
#include "validcalculator.h"
#include "tree.h"
#include "counttable.h"
#include "readmatrix.hpp"
#include "readcolumn.h"
#include "readphylip.h"
#include "sharedsobscollectsummary.h"
#include "sharedchao1.h"
#include "sharedace.h"
#include "sharednseqs.h"
#include "sharedjabund.h"
#include "sharedsorabund.h"
#include "sharedjclass.h"
#include "sharedsorclass.h"
#include "sharedjest.h"
#include "sharedsorest.h"
#include "sharedthetayc.h"
#include "sharedthetan.h"
#include "sharedkstest.h"
#include "whittaker.h"
#include "sharedochiai.h"
#include "sharedanderbergs.h"
#include "sharedkulczynski.h"
#include "sharedkulczynskicody.h"
#include "sharedlennon.h"
#include "sharedmorisitahorn.h"
#include "sharedbraycurtis.h"
#include "sharedjackknife.h"
#include "whittaker.h"
#include "odum.h"
#include "canberra.h"
#include "structeuclidean.h"
#include "structchord.h"
#include "hellinger.h"
#include "manhattan.h"
#include "structpearson.h"
#include "soergel.h"
#include "spearman.h"
#include "structkulczynski.h"
#include "structchi2.h"
#include "speciesprofile.h"
#include "hamming.h"
#include "gower.h"
#include "memchi2.h"
#include "memchord.h"
#include "memeuclidean.h"
#include "mempearson.h"
#include "sharedrjsd.h"
#include "sharedjsd.h"



/* This command create a tree file for each similarity calculator at distance level, using various calculators to find the similiarity between groups. 
	The user can select the lines or labels they wish to use as well as the groups they would like included.
	They can also use as many or as few calculators as they wish. */
	

class TreeGroupCommand : public Command {
	
public:
	TreeGroupCommand(string);	
	TreeGroupCommand();
	~TreeGroupCommand();
	
	vector<string> setParameters();
	string getCommandName()			{ return "tree.shared";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Tree.shared"; }
	string getDescription()		{ return "generate a tree file that describes the dissimilarity among groups"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
    
    struct linePair {
		int start;
		int end;
	};
	vector<linePair> lines;
    
	Tree* createTree(vector< vector<double> >&);
	void printSims(ostream&, vector< vector<double> >&);
	int makeSimsShared();
	vector< vector<double> > makeSimsDist(SparseDistanceMatrix*);
    int writeTree(string, Tree*);
    int driver(vector<SharedRAbundVector*>, int, int, vector< vector<seqDist> >&);
	
	NameAssignment* nameMap;
	ListVector* list;
	CountTable* ct;
	Tree* t;
    InputData* input;
	vector<Calculator*> treeCalculators;
	vector<SharedRAbundVector*> lookup;
	string lastLabel;
	string format, groupNames, filename, sharedfile, countfile, inputfile;
	int numGroups, subsampleSize, iters, processors;
	ofstream out;
	float precision, cutoff;

	bool abort, allLines, subsample;
	set<string> labels; //holds labels to be used
	string phylipfile, columnfile, namefile, calc, groups, label, outputDir;
	vector<string>  Estimators, Groups, outputNames; //holds estimators to be used
	
	//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
	int process(vector<SharedRAbundVector*>);
	
	

};

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct treeSharedData {
    vector<SharedRAbundVector*> thisLookup;
    vector< vector<seqDist> > calcDists;
    vector<string>  Estimators;
	unsigned long long start;
	unsigned long long end;
	MothurOut* m;
    int count;
	
	treeSharedData(){}
	treeSharedData(MothurOut* mout, unsigned long long st, unsigned long long en, vector<string> est, vector<SharedRAbundVector*> lu) {
		m = mout;
		start = st;
		end = en;
        Estimators = est;
        thisLookup = lu;
        count=0;
	}
};
/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyTreeSharedThreadFunction(LPVOID lpParam){ 
	treeSharedData* pDataArray;
	pDataArray = (treeSharedData*)lpParam;
	
	try {
        
        vector<Calculator*> treeCalculators;
        ValidCalculators validCalculator;
        for (int i=0; i<pDataArray->Estimators.size(); i++) {
            if (validCalculator.isValidCalculator("matrix", pDataArray->Estimators[i]) == true) { 
                if (pDataArray->Estimators[i] == "sharedsobs") { 
                    treeCalculators.push_back(new SharedSobsCS());
                }else if (pDataArray->Estimators[i] == "sharedchao") { 
                    treeCalculators.push_back(new SharedChao1());
                }else if (pDataArray->Estimators[i] == "sharedace") { 
                    treeCalculators.push_back(new SharedAce());
                }else if (pDataArray->Estimators[i] == "jabund") { 	
                    treeCalculators.push_back(new JAbund());
                }else if (pDataArray->Estimators[i] == "sorabund") { 
                    treeCalculators.push_back(new SorAbund());
                }else if (pDataArray->Estimators[i] == "jclass") { 
                    treeCalculators.push_back(new Jclass());
                }else if (pDataArray->Estimators[i] == "sorclass") { 
                    treeCalculators.push_back(new SorClass());
                }else if (pDataArray->Estimators[i] == "jest") { 
                    treeCalculators.push_back(new Jest());
                }else if (pDataArray->Estimators[i] == "sorest") { 
                    treeCalculators.push_back(new SorEst());
                }else if (pDataArray->Estimators[i] == "thetayc") { 
                    treeCalculators.push_back(new ThetaYC());
                }else if (pDataArray->Estimators[i] == "thetan") { 
                    treeCalculators.push_back(new ThetaN());
                }else if (pDataArray->Estimators[i] == "kstest") { 
                    treeCalculators.push_back(new KSTest());
                }else if (pDataArray->Estimators[i] == "sharednseqs") { 
                    treeCalculators.push_back(new SharedNSeqs());
                }else if (pDataArray->Estimators[i] == "ochiai") { 
                    treeCalculators.push_back(new Ochiai());
                }else if (pDataArray->Estimators[i] == "anderberg") { 
                    treeCalculators.push_back(new Anderberg());
                }else if (pDataArray->Estimators[i] == "kulczynski") { 
                    treeCalculators.push_back(new Kulczynski());
                }else if (pDataArray->Estimators[i] == "kulczynskicody") { 
                    treeCalculators.push_back(new KulczynskiCody());
                }else if (pDataArray->Estimators[i] == "lennon") { 
                    treeCalculators.push_back(new Lennon());
                }else if (pDataArray->Estimators[i] == "morisitahorn") { 
                    treeCalculators.push_back(new MorHorn());
                }else if (pDataArray->Estimators[i] == "braycurtis") { 
                    treeCalculators.push_back(new BrayCurtis());
                }else if (pDataArray->Estimators[i] == "whittaker") { 
                    treeCalculators.push_back(new Whittaker());
                }else if (pDataArray->Estimators[i] == "odum") { 
                    treeCalculators.push_back(new Odum());
                }else if (pDataArray->Estimators[i] == "canberra") { 
                    treeCalculators.push_back(new Canberra());
                }else if (pDataArray->Estimators[i] == "structeuclidean") { 
                    treeCalculators.push_back(new StructEuclidean());
                }else if (pDataArray->Estimators[i] == "structchord") { 
                    treeCalculators.push_back(new StructChord());
                }else if (pDataArray->Estimators[i] == "hellinger") { 
                    treeCalculators.push_back(new Hellinger());
                }else if (pDataArray->Estimators[i] == "manhattan") { 
                    treeCalculators.push_back(new Manhattan());
                }else if (pDataArray->Estimators[i] == "structpearson") { 
                    treeCalculators.push_back(new StructPearson());
                }else if (pDataArray->Estimators[i] == "soergel") { 
                    treeCalculators.push_back(new Soergel());
                }else if (pDataArray->Estimators[i] == "spearman") { 
                    treeCalculators.push_back(new Spearman());
                }else if (pDataArray->Estimators[i] == "structkulczynski") { 
                    treeCalculators.push_back(new StructKulczynski());
                }else if (pDataArray->Estimators[i] == "speciesprofile") { 
                    treeCalculators.push_back(new SpeciesProfile());
                }else if (pDataArray->Estimators[i] == "hamming") { 
                    treeCalculators.push_back(new Hamming());
                }else if (pDataArray->Estimators[i] == "structchi2") { 
                    treeCalculators.push_back(new StructChi2());
                }else if (pDataArray->Estimators[i] == "gower") { 
                    treeCalculators.push_back(new Gower());
                }else if (pDataArray->Estimators[i] == "memchi2") { 
                    treeCalculators.push_back(new MemChi2());
                }else if (pDataArray->Estimators[i] == "memchord") { 
                    treeCalculators.push_back(new MemChord());
                }else if (pDataArray->Estimators[i] == "memeuclidean") { 
                    treeCalculators.push_back(new MemEuclidean());
                }else if (pDataArray->Estimators[i] == "mempearson") { 
                    treeCalculators.push_back(new MemPearson());
                }
            }
        }
        
        pDataArray->calcDists.resize(treeCalculators.size());
        
		vector<SharedRAbundVector*> subset;
		for (int k = pDataArray->start; k < pDataArray->end; k++) { // pass cdd each set of groups to compare
			
            pDataArray->count++;
            
			for (int l = 0; l < k; l++) {
				
				if (k != l) { //we dont need to similiarity of a groups to itself
					subset.clear(); //clear out old pair of sharedrabunds
					//add new pair of sharedrabunds
					subset.push_back(pDataArray->thisLookup[k]); subset.push_back(pDataArray->thisLookup[l]); 
					
					for(int i=0;i<treeCalculators.size();i++) {
						
						//if this calc needs all groups to calculate the pair load all groups
						if (treeCalculators[i]->getNeedsAll()) { 
							//load subset with rest of lookup for those calcs that need everyone to calc for a pair
							for (int w = 0; w < pDataArray->thisLookup.size(); w++) {
								if ((w != k) && (w != l)) { subset.push_back(pDataArray->thisLookup[w]); }
							}
						}
						
						vector<double> tempdata = treeCalculators[i]->getValues(subset); //saves the calculator outputs
						
						if (pDataArray->m->control_pressed) { return 1; }
						
						seqDist temp(l, k, -(tempdata[0]-1.0));
						pDataArray->calcDists[i].push_back(temp);
					}
				}
			}
		}
        
        for(int i=0;i<treeCalculators.size();i++){  delete treeCalculators[i]; }
		
		return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "TreeGroupsCommand", "MyTreeSharedThreadFunction");
		exit(1);
	}
} 
#endif


	
#endif


