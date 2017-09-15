#ifndef MATRIXOUTPUTCOMMAND_H
#define MATRIXOUTPUTCOMMAND_H

/*
 *  matrixoutputcommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/20/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */ 
#include "command.hpp"
#include "inputdata.h"
#include "groupmap.h"
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
#include "sharedjsd.h"
#include "sharedrjsd.h"


// aka. dist.shared()

/* This command create a tree file for each similarity calculator at distance level, using various calculators to find the similiarity between groups. 
	The user can select the labels they wish to use as well as the groups they would like included.
	They can also use as many or as few calculators as they wish. */
	

class MatrixOutputCommand : public Command {
	
public:
	MatrixOutputCommand(string);
	MatrixOutputCommand();	
	~MatrixOutputCommand();
	
	vector<string> setParameters();
	string getCommandName()			{ return "dist.shared";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Dist.shared"; }
	string getDescription()		{ return "generate a distance matrix that describes the dissimilarity among multiple groups"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	vector<linePair> lines;
	
	void printSims(ostream&, vector< vector<double> >&, vector<string>);
	int process(SharedRAbundVectors*&);
	
	vector<Calculator*> matrixCalculators;
	//SharedRAbundVectors* lookup;
	string exportFileName, output, sharedfile;
	int numGroups, processors, iters, subsampleSize;
	ofstream out;

	bool abort, allLines, subsample;
	set<string> labels; //holds labels to be used
	string outputFile, calc, groups, label, outputDir, mode;
	vector<string>  Estimators, Groups, outputNames; //holds estimators to be used
	int process(SharedRAbundVectors*&, string, string);
	int driver(vector<SharedRAbundVector*>&, int, int, vector< vector<seqDist> >&);

};
	
/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct distSharedData {
    vector<SharedRAbundVector*> thisLookup;
    vector< vector<seqDist> > calcDists;
    vector<string>  Estimators;
	unsigned long long start;
	unsigned long long end;
	MothurOut* m;
    int count;
	
	distSharedData(){}
	distSharedData(MothurOut* mout, unsigned long long st, unsigned long long en, vector<string> est, vector<SharedRAbundVector*>& lu) {
		m = mout;
		start = st;
		end = en;
        Estimators = est;
        thisLookup = lu;
        count = 0;
	}
};
/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyDistSharedThreadFunction(LPVOID lpParam){ 
	distSharedData* pDataArray;
	pDataArray = (distSharedData*)lpParam;
	
	try {
        
        vector<Calculator*> matrixCalculators;
        ValidCalculators validCalculator;
        for (int i=0; i<pDataArray->Estimators.size(); i++) {
            if (validCalculator.isValidCalculator("matrix", pDataArray->Estimators[i]) == true) { 
                if (pDataArray->Estimators[i] == "sharedsobs") { 
                    matrixCalculators.push_back(new SharedSobsCS());
                }else if (pDataArray->Estimators[i] == "sharedchao") { 
                    matrixCalculators.push_back(new SharedChao1());
                }else if (pDataArray->Estimators[i] == "sharedace") { 
                    matrixCalculators.push_back(new SharedAce());
                }else if (pDataArray->Estimators[i] == "jabund") { 	
                    matrixCalculators.push_back(new JAbund());
                }else if (pDataArray->Estimators[i] == "sorabund") { 
                    matrixCalculators.push_back(new SorAbund());
                }else if (pDataArray->Estimators[i] == "jclass") { 
                    matrixCalculators.push_back(new Jclass());
                }else if (pDataArray->Estimators[i] == "sorclass") { 
                    matrixCalculators.push_back(new SorClass());
                }else if (pDataArray->Estimators[i] == "jest") { 
                    matrixCalculators.push_back(new Jest());
                }else if (pDataArray->Estimators[i] == "sorest") { 
                    matrixCalculators.push_back(new SorEst());
                }else if (pDataArray->Estimators[i] == "thetayc") { 
                    matrixCalculators.push_back(new ThetaYC());
                }else if (pDataArray->Estimators[i] == "thetan") { 
                    matrixCalculators.push_back(new ThetaN());
                }else if (pDataArray->Estimators[i] == "kstest") { 
                    matrixCalculators.push_back(new KSTest());
                }else if (pDataArray->Estimators[i] == "sharednseqs") { 
                    matrixCalculators.push_back(new SharedNSeqs());
                }else if (pDataArray->Estimators[i] == "ochiai") { 
                    matrixCalculators.push_back(new Ochiai());
                }else if (pDataArray->Estimators[i] == "anderberg") { 
                    matrixCalculators.push_back(new Anderberg());
                }else if (pDataArray->Estimators[i] == "kulczynski") { 
                    matrixCalculators.push_back(new Kulczynski());
                }else if (pDataArray->Estimators[i] == "kulczynskicody") { 
                    matrixCalculators.push_back(new KulczynskiCody());
                }else if (pDataArray->Estimators[i] == "lennon") { 
                    matrixCalculators.push_back(new Lennon());
                }else if (pDataArray->Estimators[i] == "morisitahorn") { 
                    matrixCalculators.push_back(new MorHorn());
                }else if (pDataArray->Estimators[i] == "braycurtis") { 
                    matrixCalculators.push_back(new BrayCurtis());
                }else if (pDataArray->Estimators[i] == "whittaker") { 
                    matrixCalculators.push_back(new Whittaker());
                }else if (pDataArray->Estimators[i] == "odum") { 
                    matrixCalculators.push_back(new Odum());
                }else if (pDataArray->Estimators[i] == "canberra") { 
                    matrixCalculators.push_back(new Canberra());
                }else if (pDataArray->Estimators[i] == "structeuclidean") { 
                    matrixCalculators.push_back(new StructEuclidean());
                }else if (pDataArray->Estimators[i] == "structchord") { 
                    matrixCalculators.push_back(new StructChord());
                }else if (pDataArray->Estimators[i] == "hellinger") { 
                    matrixCalculators.push_back(new Hellinger());
                }else if (pDataArray->Estimators[i] == "manhattan") { 
                    matrixCalculators.push_back(new Manhattan());
                }else if (pDataArray->Estimators[i] == "structpearson") { 
                    matrixCalculators.push_back(new StructPearson());
                }else if (pDataArray->Estimators[i] == "soergel") { 
                    matrixCalculators.push_back(new Soergel());
                }else if (pDataArray->Estimators[i] == "spearman") { 
                    matrixCalculators.push_back(new Spearman());
                }else if (pDataArray->Estimators[i] == "structkulczynski") { 
                    matrixCalculators.push_back(new StructKulczynski());
                }else if (pDataArray->Estimators[i] == "speciesprofile") { 
                    matrixCalculators.push_back(new SpeciesProfile());
                }else if (pDataArray->Estimators[i] == "hamming") { 
                    matrixCalculators.push_back(new Hamming());
                }else if (pDataArray->Estimators[i] == "structchi2") { 
                    matrixCalculators.push_back(new StructChi2());
                }else if (pDataArray->Estimators[i] == "gower") { 
                    matrixCalculators.push_back(new Gower());
                }else if (pDataArray->Estimators[i] == "memchi2") { 
                    matrixCalculators.push_back(new MemChi2());
                }else if (pDataArray->Estimators[i] == "memchord") { 
                    matrixCalculators.push_back(new MemChord());
                }else if (pDataArray->Estimators[i] == "memeuclidean") { 
                    matrixCalculators.push_back(new MemEuclidean());
                }else if (pDataArray->Estimators[i] == "mempearson") { 
                    matrixCalculators.push_back(new MemPearson());
                }else if (pDataArray->Estimators[i] == "jsd") {
                    matrixCalculators.push_back(new JSD());
                }else if (pDataArray->Estimators[i] == "rjsd") {
                    matrixCalculators.push_back(new RJSD());
                }

            }
        }
        
        pDataArray->calcDists.resize(matrixCalculators.size());
        		
		vector<SharedRAbundVector*> subset;
		for (int k = pDataArray->start; k < pDataArray->end; k++) { // pass cdd each set of groups to compare
			pDataArray->count++;
			for (int l = 0; l < k; l++) {
				
				if (k != l) { //we dont need to similiarity of a groups to itself
					subset.clear(); //clear out old pair of sharedrabunds
					//add new pair of sharedrabunds
					subset.push_back(pDataArray->thisLookup[k]); subset.push_back(pDataArray->thisLookup[l]); 
					
					for(int i=0;i<matrixCalculators.size();i++) {
						
						//if this calc needs all groups to calculate the pair load all groups
						if (matrixCalculators[i]->getNeedsAll()) { 
							//load subset with rest of lookup for those calcs that need everyone to calc for a pair
							for (int w = 0; w < pDataArray->thisLookup.size(); w++) {
								if ((w != k) && (w != l)) { subset.push_back(pDataArray->thisLookup[w]); }
							}
						}
						
						vector<double> tempdata = matrixCalculators[i]->getValues(subset); //saves the calculator outputs
						
						if (pDataArray->m->getControl_pressed()) { return 1; }
						
						seqDist temp(l, k, tempdata[0]);
						pDataArray->calcDists[i].push_back(temp);
					}
				}
			}
		}
        
        for(int i=0;i<matrixCalculators.size();i++){  delete matrixCalculators[i]; }
        for(int i=0;i<pDataArray->thisLookup.size();i++){  delete pDataArray->thisLookup[i]; }
		
		return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "MatrixOutputCommand", "MyDistSharedThreadFunction");
		exit(1);
	}
} 
#endif
	
#endif

