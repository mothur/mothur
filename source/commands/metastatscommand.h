#ifndef METASTATSCOMMAND_H
#define METASTATSCOMMAND_H

/*
 *  metastatscommand.h
 *  Mothur
 *
 *  Created by westcott on 9/16/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */
 
#include "command.hpp"
#include "inputdata.h"
#include "sharedrabundvectors.hpp"
#include "mothurmetastats.h"
#include "designmap.h"

class MetaStatsCommand : public Command {

public:
	MetaStatsCommand(string);
	MetaStatsCommand();
	~MetaStatsCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "metastats";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "White JR, Nagarajan N, Pop M (2009). Statistical methods for detecting differentially abundant features in clinical metagenomic samples. PLoS Comput Biol 5: e1000352. \nhttp://www.mothur.org/wiki/Metastats"; }
	string getDescription()		{ return "detects differentially abundant features in clinical metagenomic samples"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	vector<linePair> lines;
	
	DesignMap* designMap;
	SharedRAbundVectors* lookup;
		
	bool abort, allLines, pickedGroups;
	set<string> labels; //holds labels to be used
	string groups, label, outputDir, inputDir, designfile, sets, sharedfile;
	vector<string> Groups, outputNames, Sets;
	vector< vector<string> > namesOfGroupCombos;
	int iters, processors;
	float threshold;
	
	int process(SharedRAbundVectors*);
	int driver(unsigned long long, unsigned long long, SharedRAbundVectors*);
    int convertToShared(string filename);
    int convertToInput(vector<RAbundVector*>&, vector<string>, string);
    bool convertSharedToInput;
};

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct metastatsData {
    SharedRAbundVectors* thisLookUp;
    vector< vector<string> > namesOfGroupCombos;
    vector<string> designMapGroups;
    vector<string> outputNames;
	int start;
	int num, iters, count;
	float threshold;
	MothurOut* m;
	string sharedfile;
    string outputDir;
	
	metastatsData(){}
	metastatsData(string sf, string oDir, MothurOut* mout, int st, int en, vector< vector<string> > ns, SharedRAbundVectors* lu, vector<string> dg, int i, float thr) {
		sharedfile = sf;
        outputDir = oDir;
		m = mout;
		start = st;
		num = en;
        namesOfGroupCombos = ns;
        thisLookUp = lu;
        designMapGroups = dg;
        iters = i;
        threshold = thr;
        count=0;
	}
};
/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyMetastatsThreadFunction(LPVOID lpParam){ 
	metastatsData* pDataArray;
	pDataArray = (metastatsData*)lpParam;
	
	try {
		
        vector<string> thisLookupNames = pDataArray->thisLookUp->getNamesGroups();
        vector<RAbundVector*> thisLookupRabunds = pDataArray->thisLookUp->getSharedRAbundVectors();
        
        //for each combo
		for (int c = pDataArray->start; c < (pDataArray->start+pDataArray->num); c++) {
			pDataArray->count++;
			//get set names
			string setA = pDataArray->namesOfGroupCombos[c][0]; 
			string setB = pDataArray->namesOfGroupCombos[c][1];
            
			//get filename
			string outputFileName = pDataArray->outputDir +  pDataArray->m->getRootName(pDataArray->m->getSimpleName(pDataArray->sharedfile)) + pDataArray->thisLookUp[0]->getLabel() + "." + setA + "-" + setB + ".metastats";
			pDataArray->outputNames.push_back(outputFileName); 
			
			vector< vector<double> > data2; data2.resize(pDataArray->thisLookUp[0]->getNumBins());
			
			vector<RAbundVector*> subset;
			int setACount = 0;
			int setBCount = 0;
			for (int i = 0; i < pDataArray->thisLookUp.size(); i++) {
				//is this group for a set we want to compare??
				//sorting the sets by putting setB at the back and setA in the front
				if (pDataArray->designMapGroups[i] == setB) {  
                    subset.push_back(thisLookupRabunds[i]);
					setBCount++;
				}else if (pDataArray->designMapGroups[i] == setA) {
                    subset.insert(subset.begin()+setACount, thisLookupRabunds[i]);
					setACount++;
				}
			}
            
			if ((setACount == 0) || (setBCount == 0))  { 
				pDataArray->m->mothurOut("Missing shared info for " + setA + " or " + setB + ". Skipping comparison."); pDataArray->m->mothurOutEndLine(); 
				pDataArray->outputNames.pop_back();
			}else {
				//fill data
				for (int j = 0; j < pDataArray->thisLookUp[0]->getNumBins(); j++) {
					data2[j].resize(subset.size(), 0.0);
					for (int i = 0; i < subset.size(); i++) {
						data2[j][i] = (subset[i]->getAbundance(j));
					}
				}
				
				pDataArray->m->mothurOut("Comparing " + setA + " and " + setB + "..."); pDataArray->m->mothurOutEndLine(); 
				
				pDataArray->m->mothurOutEndLine();
				MothurMetastats mothurMeta(pDataArray->threshold, pDataArray->iters);
				mothurMeta.runMetastats(outputFileName, data2, setACount);
				pDataArray->m->mothurOutEndLine();
				pDataArray->m->mothurOutEndLine(); 
			}
        }
		
        for(int i = 0; i < thisLookupRabunds.size(); i++)  {  delete thisLookupRabunds[i];  }
		return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "MetaStatsCommand", "MyMetastatsThreadFunction");
		exit(1);
	}
} 
#endif



#endif

