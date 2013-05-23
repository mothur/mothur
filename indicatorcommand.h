#ifndef INDICATORCOMMAND_H
#define INDICATORCOMMAND_H

/*
 *  indicatorcommand.h
 *  Mothur
 *
 *  Created by westcott on 11/12/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "readtree.h"
#include "counttable.h"
#include "sharedrabundvector.h"
#include "sharedrabundfloatvector.h"
#include "inputdata.h"

class IndicatorCommand : public Command {
public:
	IndicatorCommand(string);
	IndicatorCommand();
	~IndicatorCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "indicator";				}
	string getCommandCategory()		{ return "Hypothesis Testing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Dufrene M, Legendre P (1997). Species assemblages and indicator species: The need for a flexible asymmetrical approach. Ecol Monogr 67: 345-66.\n McCune B, Grace JB, Urban DL (2002). Analysis of ecological communities. MjM Software Design: Gleneden Beach, OR. \nLegendre P, Legendre L (1998). Numerical Ecology. Elsevier: New York. \nhttp://www.mothur.org/wiki/Indicator"; }
	string getDescription()		{ return "calculate the indicator value for each OTU"; }

	int execute();
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	ReadTree* read;
	CountTable* ct;
	GroupMap* designMap;
	string treefile, sharedfile, relabundfile, groups, label, inputFileName, outputDir, designfile;
	bool abort;
	int iters, processors;
	vector<string> outputNames, Groups;
	vector<SharedRAbundVector*> lookup;
	vector<SharedRAbundFloatVector*> lookupFloat;
	
	int getShared();
	int getSharedFloat();
	int GetIndicatorSpecies(Tree*&);
	int GetIndicatorSpecies();
	set<string> getDescendantList(Tree*&, int, map<int, set<string> >, map<int, set<int> >&);
	vector<float> getValues(vector< vector<SharedRAbundVector*> >&, vector<string>&, map< vector<int>, vector<int> >);
	vector<float> getValues(vector< vector<SharedRAbundFloatVector*> >&, vector<string>&, map< vector<int>, vector<int> >);
    
	map<int, float> getDistToRoot(Tree*&);
	map< vector<int>, vector<int> > randomizeGroupings(vector< vector<SharedRAbundVector*> >&, int);
	map< vector<int>, vector<int> > randomizeGroupings(vector< vector<SharedRAbundFloatVector*> >&, int);
    
	vector<float> driver(vector< vector<SharedRAbundFloatVector*> >&, int, vector<float>, int);
	vector<float> driver(vector< vector<SharedRAbundVector*> >&, int, vector<float>, int);
    
	vector<float> getPValues(vector< vector<SharedRAbundFloatVector*> >&, int, vector<float>);
	vector<float> getPValues(vector< vector<SharedRAbundVector*> >&, int, vector<float>);

	
};

/**************************************************************************************************/

struct indicatorData {
    vector< vector<SharedRAbundFloatVector*> > groupings;
   	MothurOut* m;
    int iters, num;
    vector<float> indicatorValues;
    vector<float> pvalues;
	
	indicatorData(){}
	indicatorData(MothurOut* mout, int it, vector< vector<SharedRAbundFloatVector*> > ng, int n, vector<float> iv) {
		m = mout;
        iters = it;
        groupings = ng;
        indicatorValues = iv;
        num = n;
    }
};
/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyIndicatorThreadFunction(LPVOID lpParam){
	indicatorData* pDataArray;
	pDataArray = (indicatorData*)lpParam;
	
	try {
        
		pDataArray->pvalues.resize(pDataArray->indicatorValues.size(), 0);
		
		for(int i=0;i<pDataArray->iters;i++){
			if (pDataArray->m->control_pressed) { break; }
            
			//groupingsMap = randomizeGroupings(groupings, num);
            ///////////////////////////////////////////////////////////////////////
            map< vector<int>, vector<int> > randomGroupings;
            
            for (int j = 0; j < pDataArray->num; j++) {
                
                //get random groups to swap to switch with
                //generate random int between 0 and groupings.size()-1
                int z = pDataArray->m->getRandomIndex(pDataArray->groupings.size()-1);
                int x = pDataArray->m->getRandomIndex(pDataArray->groupings.size()-1);
                int a = pDataArray->m->getRandomIndex(pDataArray->groupings[z].size()-1);
                int b = pDataArray->m->getRandomIndex(pDataArray->groupings[x].size()-1);
                //cout << i << '\t' << z << '\t' << x << '\t' << a << '\t' << b << endl;
                
                vector<int> from;
                vector<int> to;
                
                from.push_back(z); from.push_back(a);
                to.push_back(x); to.push_back(b);
                
                randomGroupings[from] = to;
            }
            ///////////////////////////////////////////////////////////////////////
            
			//vector<float> randomIndicatorValues = getValues(groupings, notUsedGroupings, randomGroupings);
            ///////////////////////////////////////////////////////////////////////
            vector<float> randomIndicatorValues;
            map< vector<int>, vector<int> >::iterator it;
            
            //for each otu
            for (int i = 0; i < pDataArray->groupings[0][0]->getNumBins(); i++) {
                
                if (pDataArray->m->control_pressed) { return 0; }
                
                vector<float> terms;
                float AijDenominator = 0.0;
                vector<float> Bij;
                
                //get overall abundance of each grouping
                for (int j = 0; j < pDataArray->groupings.size(); j++) {
                    
                    float totalAbund = 0;
                    int numNotZero = 0;
                    for (int k = 0; k < pDataArray->groupings[j].size(); k++) {
                        vector<int> temp; temp.push_back(j); temp.push_back(k);
                        it = randomGroupings.find(temp);
                        
                        if (it == randomGroupings.end()) { //this one didnt get moved
                            totalAbund += pDataArray->groupings[j][k]->getAbundance(i);
                            if (pDataArray->groupings[j][k]->getAbundance(i) != 0.0) { numNotZero++; }
                        }else {
                            totalAbund += pDataArray->groupings[(it->second)[0]][(it->second)[1]]->getAbundance(i);
                            if (pDataArray->groupings[(it->second)[0]][(it->second)[1]]->getAbundance(i) != 0.0) { numNotZero++; }
                        }
                        
                    }
                    
                    //mean abundance
                    float Aij = (totalAbund / (float) pDataArray->groupings[j].size());
                    terms.push_back(Aij);
                    
                    //percentage of sites represented
                    Bij.push_back(numNotZero / (float) pDataArray->groupings[j].size());
                    
                    AijDenominator += Aij;
                }
                
                float maxIndVal = 0.0;
                for (int j = 0; j < terms.size(); j++) { 
                    float thisAij = (terms[j] / AijDenominator); //relative abundance
                    float thisValue = thisAij * Bij[j] * 100.0;
                    
                    //save largest
                    if (thisValue > maxIndVal) { maxIndVal = thisValue; }
                }
                
                randomIndicatorValues.push_back(maxIndVal);
            }

            ///////////////////////////////////////////////////////////////////////
			
			for (int j = 0; j < pDataArray->indicatorValues.size(); j++) {
				if (randomIndicatorValues[j] >= pDataArray->indicatorValues[j]) { pDataArray->pvalues[j]++; }
			}
		}

        return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "IndicatorCommand", "MyIndicatorThreadFunction");
		exit(1);
	}
}
#endif



#endif

