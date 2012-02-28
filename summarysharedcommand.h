#ifndef SUMMARYSHAREDCOMMAND_H
#define SUMMARYSHAREDCOMMAND_H
/*
 *  summarysharedcommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "command.hpp"
#include "sharedrabundvector.h"
#include "inputdata.h"
#include "calculator.h"
#include "validcalculator.h"
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

class SummarySharedCommand : public Command {

public:
	SummarySharedCommand(string);
	SummarySharedCommand();
	~SummarySharedCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "summary.shared";			}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Summary.shared"; }
	string getDescription()		{ return "generate a summary file containing calculator values for each line in the OTU data and for all possible comparisons between groups"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	struct linePair {
		int start;
		int end;
	};
	vector<linePair> lines;
	vector<Calculator*> sumCalculators;	
	InputData* input;
	
	bool abort, allLines, mult, all, createPhylip;
	set<string> labels; //holds labels to be used
	string label, calc, groups, sharedfile;
	vector<string>  Estimators, Groups, outputNames;
	vector<SharedRAbundVector*> lookup;
	string format, outputDir;
	int numGroups, processors;
	int process(vector<SharedRAbundVector*>, string, string);
	int driver(vector<SharedRAbundVector*>, int, int, string, string, vector< vector<seqDist> >&);

};

/**************************************************************************************************/
//custom data structure for threads to use.
//main process handling the calcs that can do more than 2 groups
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct summarySharedData {
    vector<SharedRAbundVector*> thisLookup;
    vector< vector<seqDist> > calcDists;
    vector<string>  Estimators;
	unsigned long long start;
	unsigned long long end;
	MothurOut* m;
	string sumFile;
	
	summarySharedData(){}
	summarySharedData(string sf, MothurOut* mout, unsigned long long st, unsigned long long en, vector<string> est, vector<SharedRAbundVector*> lu) {
		sumFile = sf;
		m = mout;
		start = st;
		end = en;
        Estimators = est;
        thisLookup = lu;
	}
};
/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
#else
static DWORD WINAPI MySummarySharedThreadFunction(LPVOID lpParam){ 
	summarySharedData* pDataArray;
	pDataArray = (summarySharedData*)lpParam;
	
	try {
        
        vector<Calculator*> sumCalculators;
        ValidCalculators validCalculator;
        for (int i=0; i<pDataArray->Estimators.size(); i++) {
            if (validCalculator.isValidCalculator("sharedsummary", pDataArray->Estimators[i]) == true) { 
                if (pDataArray->Estimators[i] == "sharedsobs") { 
                    sumCalculators.push_back(new SharedSobsCS());
                }else if (pDataArray->Estimators[i] == "sharedchao") { 
                    sumCalculators.push_back(new SharedChao1());
                }else if (pDataArray->Estimators[i] == "sharedace") { 
                    sumCalculators.push_back(new SharedAce());
                }else if (pDataArray->Estimators[i] == "jabund") { 	
                    sumCalculators.push_back(new JAbund());
                }else if (pDataArray->Estimators[i] == "sorabund") { 
                    sumCalculators.push_back(new SorAbund());
                }else if (pDataArray->Estimators[i] == "jclass") { 
                    sumCalculators.push_back(new Jclass());
                }else if (pDataArray->Estimators[i] == "sorclass") { 
                    sumCalculators.push_back(new SorClass());
                }else if (pDataArray->Estimators[i] == "jest") { 
                    sumCalculators.push_back(new Jest());
                }else if (pDataArray->Estimators[i] == "sorest") { 
                    sumCalculators.push_back(new SorEst());
                }else if (pDataArray->Estimators[i] == "thetayc") { 
                    sumCalculators.push_back(new ThetaYC());
                }else if (pDataArray->Estimators[i] == "thetan") { 
                    sumCalculators.push_back(new ThetaN());
                }else if (pDataArray->Estimators[i] == "kstest") { 
                    sumCalculators.push_back(new KSTest());
                }else if (pDataArray->Estimators[i] == "sharednseqs") { 
                    sumCalculators.push_back(new SharedNSeqs());
                }else if (pDataArray->Estimators[i] == "ochiai") { 
                    sumCalculators.push_back(new Ochiai());
                }else if (pDataArray->Estimators[i] == "anderberg") { 
                    sumCalculators.push_back(new Anderberg());
                }else if (pDataArray->Estimators[i] == "kulczynski") { 
                    sumCalculators.push_back(new Kulczynski());
                }else if (pDataArray->Estimators[i] == "kulczynskicody") { 
                    sumCalculators.push_back(new KulczynskiCody());
                }else if (pDataArray->Estimators[i] == "lennon") { 
                    sumCalculators.push_back(new Lennon());
                }else if (pDataArray->Estimators[i] == "morisitahorn") { 
                    sumCalculators.push_back(new MorHorn());
                }else if (pDataArray->Estimators[i] == "braycurtis") { 
                    sumCalculators.push_back(new BrayCurtis());
                }else if (pDataArray->Estimators[i] == "whittaker") { 
                    sumCalculators.push_back(new Whittaker());
                }else if (pDataArray->Estimators[i] == "odum") { 
                    sumCalculators.push_back(new Odum());
                }else if (pDataArray->Estimators[i] == "canberra") { 
                    sumCalculators.push_back(new Canberra());
                }else if (pDataArray->Estimators[i] == "structeuclidean") { 
                    sumCalculators.push_back(new StructEuclidean());
                }else if (pDataArray->Estimators[i] == "structchord") { 
                    sumCalculators.push_back(new StructChord());
                }else if (pDataArray->Estimators[i] == "hellinger") { 
                    sumCalculators.push_back(new Hellinger());
                }else if (pDataArray->Estimators[i] == "manhattan") { 
                    sumCalculators.push_back(new Manhattan());
                }else if (pDataArray->Estimators[i] == "structpearson") { 
                    sumCalculators.push_back(new StructPearson());
                }else if (pDataArray->Estimators[i] == "soergel") { 
                    sumCalculators.push_back(new Soergel());
                }else if (pDataArray->Estimators[i] == "spearman") { 
                    sumCalculators.push_back(new Spearman());
                }else if (pDataArray->Estimators[i] == "structkulczynski") { 
                    sumCalculators.push_back(new StructKulczynski());
                }else if (pDataArray->Estimators[i] == "speciesprofile") { 
                    sumCalculators.push_back(new SpeciesProfile());
                }else if (pDataArray->Estimators[i] == "hamming") { 
                    sumCalculators.push_back(new Hamming());
                }else if (pDataArray->Estimators[i] == "structchi2") { 
                    sumCalculators.push_back(new StructChi2());
                }else if (pDataArray->Estimators[i] == "gower") { 
                    sumCalculators.push_back(new Gower());
                }else if (pDataArray->Estimators[i] == "memchi2") { 
                    sumCalculators.push_back(new MemChi2());
                }else if (pDataArray->Estimators[i] == "memchord") { 
                    sumCalculators.push_back(new MemChord());
                }else if (pDataArray->Estimators[i] == "memeuclidean") { 
                    sumCalculators.push_back(new MemEuclidean());
                }else if (pDataArray->Estimators[i] == "mempearson") { 
                    sumCalculators.push_back(new MemPearson());
                }
            }
        }
        
        pDataArray->calcDists.resize(sumCalculators.size());
        
		ofstream outputFileHandle;
		pDataArray->m->openOutputFile(pDataArray->sumFile, outputFileHandle);
		
		vector<SharedRAbundVector*> subset;
		for (int k = pDataArray->start; k < pDataArray->end; k++) { // pass cdd each set of groups to compare
            
			for (int l = 0; l < k; l++) {
				
				outputFileHandle << pDataArray->thisLookup[0]->getLabel() << '\t';
				
				subset.clear(); //clear out old pair of sharedrabunds
				//add new pair of sharedrabunds
				subset.push_back(pDataArray->thisLookup[k]); subset.push_back(pDataArray->thisLookup[l]); 
				
				//sort groups to be alphanumeric
				if (pDataArray->thisLookup[k]->getGroup() > pDataArray->thisLookup[l]->getGroup()) {
					outputFileHandle << (pDataArray->thisLookup[l]->getGroup() +'\t' + pDataArray->thisLookup[k]->getGroup()) << '\t'; //print out groups
				}else{
					outputFileHandle << (pDataArray->thisLookup[k]->getGroup() +'\t' + pDataArray->thisLookup[l]->getGroup()) << '\t'; //print out groups
				}
				
				for(int i=0;i<sumCalculators.size();i++) {
					
					//if this calc needs all groups to calculate the pair load all groups
					if (sumCalculators[i]->getNeedsAll()) { 
						//load subset with rest of lookup for those calcs that need everyone to calc for a pair
						for (int w = 0; w < pDataArray->thisLookup.size(); w++) {
							if ((w != k) && (w != l)) { subset.push_back(pDataArray->thisLookup[w]); }
						}
					}
					
					vector<double> tempdata = sumCalculators[i]->getValues(subset); //saves the calculator outputs
					
					if (pDataArray->m->control_pressed) { for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; } outputFileHandle.close(); return 1; }
					
					outputFileHandle << '\t';
					sumCalculators[i]->print(outputFileHandle);
					
					seqDist temp(l, k, tempdata[0]);
					pDataArray->calcDists[i].push_back(temp);
				}
				outputFileHandle << endl;
			}
		}
		
		outputFileHandle.close();
        for(int i=0;i<sumCalculators.size();i++){  delete sumCalculators[i]; }
		
		return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "SummarySharedCommand", "MySummarySharedThreadFunction");
		exit(1);
	}
} 
#endif


#endif
