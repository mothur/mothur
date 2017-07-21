//
//  sparcccommand.h
//  Mothur
//
//  Created by SarahsWork on 5/10/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef Mothur_sparcccommand_h
#define Mothur_sparcccommand_h

#include "command.hpp"
#include "inputdata.h"
#include "calcsparcc.h"

/**************************************************************************************************/

class SparccCommand : public Command {
public:
    SparccCommand(string);
    SparccCommand();
    ~SparccCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "sparcc";			}
    string getCommandCategory()		{ return "OTU-Based Approaches";		}
    
    string getOutputPattern(string);
    //commmand category choices: Sequence Processing, OTU-Based Approaches, Hypothesis Testing, Phylotype Analysis, General, Clustering and Hidden
	string getHelpString();
    string getCitation() { return "Friedman J, Alm EJ (2012) Inferring Correlation Networks from Genomic Survey Data. PLoS Comput Biol 8(9): e1002687. doi:10.1371/journal.pcbi.1002687 http://www.mothur.org/wiki/Sparcc"; }
    string getDescription()		{ return "Calculates correlations between OTUs using a method that is insensitive to the use of relative abundance data"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    bool abort, allLines;
    string outputDir, sharedfile, normalizeMethod;
    int numSamplings, maxIterations, numPermutations, processors;
    set<string> labels;
    vector<string> Groups;
    vector<string> outputNames;
    
    int process(SharedRAbundVectors*);
    vector<vector<float> > createProcesses(vector<vector<float> >&, vector<vector<float> >&);
    vector<vector<float> > driver(vector<vector<float> >&, vector<vector<float> >&, int);
    vector<vector<float> > shuffleSharedVector(vector<vector<float> >&);
};

/**************************************************************************************************/

struct sparccData {
   	MothurOut* m;
    int numPerms;
    vector< vector<float> > sharedVector;
    vector< vector<float> > origCorrMatrix;
    vector<vector<float> > pValues;
    int numSamplings, maxIterations, numPermutations;
    string normalizeMethod;
	
	sparccData(){}
	sparccData(MothurOut* mout, int it, vector< vector<float> > cs, vector< vector<float> > co, int ns, int mi, int np, string nm) {
		m = mout;
        numPerms = it;
        sharedVector = cs;
        origCorrMatrix = co;
        numSamplings = ns;
        maxIterations = mi;
        numPermutations = np;
        normalizeMethod = nm;
    }
};
/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MySparccThreadFunction(LPVOID lpParam){
	sparccData* pDataArray;
	pDataArray = (sparccData*)lpParam;
	
	try {
        
        int numOTUs = pDataArray->sharedVector[0].size();
        vector<vector<float> > sharedShuffled = pDataArray->sharedVector;
        pDataArray->pValues.resize(numOTUs);
        for(int i=0;i<numOTUs;i++){ pDataArray->pValues[i].assign(numOTUs, 0);  }
        
        for(int i=0;i<pDataArray->numPerms;i++){
            if (pDataArray->m->control_pressed) { return 0; }
            
            //sharedShuffled = shuffleSharedVector(sharedVector);
            //////////////////////////////////////////////////////////
            int numGroups = (int)pDataArray->sharedVector.size();
            sharedShuffled = pDataArray->sharedVector;
            
            for(int k=0;k<numGroups;k++){
                for(int j=0;j<numOTUs;j++){
                    sharedShuffled[k][j] = pDataArray->sharedVector[pDataArray->m->getRandomIndex(numGroups-1)][j];
                }
            }
            /////////////////////////////////////////////////////////
            
            CalcSparcc permutedData(sharedShuffled, pDataArray->maxIterations, pDataArray->numSamplings, pDataArray->normalizeMethod);
            vector<vector<float> > permuteCorrMatrix = permutedData.getRho();
            
            for(int j=0;j<numOTUs;j++){
                for(int k=0;k<j;k++){
                    double randValue = permuteCorrMatrix[j][k];
                    double observedValue = pDataArray->origCorrMatrix[j][k];
                    if(observedValue >= 0 &&  randValue > observedValue)   { pDataArray->pValues[j][k]++; }//this method seems to deflate the
                    else if(observedValue < 0 && randValue < observedValue){ pDataArray->pValues[j][k]++; }//pvalues of small rho values
                }
            }
            if((i+1) % (int)(pDataArray->numPermutations * 0.05) == 0){ cout << i+1 << endl;  }
        }
        
        return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "SparccCommand", "MySparccThreadFunction");
		exit(1);
	}
}
#endif


/**************************************************************************************************/




#endif
