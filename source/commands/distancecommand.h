#ifndef DISTANCECOMMAND_H
#define DISTANCECOMMAND_H

/*
 *  distancecommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "validcalculator.h"
#include "dist.h"
#include "sequencedb.h"
#include "ignoregaps.h"
#include "eachgapdist.h"
#include "eachgapignore.h"
#include "onegapdist.h"
#include "onegapignore.h"

//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct distanceData {
	int startLine;
	int endLine;
	string dFileName;
	float cutoff;
	SequenceDB alignDB;
	vector<string> Estimators;
	MothurOut* m;
	string output;
	int numNewFasta, count;
	string countends;
	
	distanceData(){}
	distanceData(int s, int e, string dbname, float c, SequenceDB db, vector<string> Est, MothurOut* mout, string o, int num, string count) {
		startLine = s;
		endLine = e;
		dFileName = dbname;
		cutoff = c;
		alignDB = db;
		Estimators = Est;
		m = mout;
		output = o;
		numNewFasta = num;
		countends = count;
		
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyDistThreadFunction(LPVOID lpParam){ 
	distanceData* pDataArray;
	pDataArray = (distanceData*)lpParam;
	
	try {
		ValidCalculators validCalculator;
		Dist* distCalculator;
		if (pDataArray->m->isTrue(pDataArray->countends) == true) {
			for (int i=0; i<pDataArray->Estimators.size(); i++) {
				if (validCalculator.isValidCalculator("distance", pDataArray->Estimators[i]) == true) { 
					if (pDataArray->Estimators[i] == "nogaps")			{	distCalculator = new ignoreGaps();	}
					else if (pDataArray->Estimators[i] == "eachgap")	{	distCalculator = new eachGapDist();	}
					else if (pDataArray->Estimators[i] == "onegap")		{	distCalculator = new oneGapDist();	}
				}
			}
		}else {
			for (int i=0; i<pDataArray->Estimators.size(); i++) {
				if (validCalculator.isValidCalculator("distance", pDataArray->Estimators[i]) == true) { 
					if (pDataArray->Estimators[i] == "nogaps")		{	distCalculator = new ignoreGaps();					}
					else if (pDataArray->Estimators[i] == "eachgap"){	distCalculator = new eachGapIgnoreTermGapDist();	}
					else if (pDataArray->Estimators[i] == "onegap")	{	distCalculator = new oneGapIgnoreTermGapDist();		}
				}
			}
		}
		
		int startTime = time(NULL);
		
		//column file
		ofstream outFile(pDataArray->dFileName.c_str(), ios::trunc);
		outFile.setf(ios::fixed, ios::showpoint);
		outFile << setprecision(4);
		pDataArray->count = 0;
		
		if (pDataArray->output != "square") { 
			if((pDataArray->output == "lt") && (pDataArray->startLine == 0)){	outFile << pDataArray->alignDB.getNumSeqs() << endl;	}
			
			for(int i=pDataArray->startLine;i<pDataArray->endLine;i++){
				if(pDataArray->output == "lt")	{	
					string name = pDataArray->alignDB.get(i).getName();
					if (name.length() < 10) { //pad with spaces to make compatible
						while (name.length() < 10) {  name += " ";  }
					}
					outFile << name;
				}
				for(int j=0;j<i;j++){
					
					if (pDataArray->m->control_pressed) { delete distCalculator; outFile.close(); return 0;  }
					
					//if there was a column file given and we are appending, we don't want to calculate the distances that are already in the column file
					//the alignDB contains the new sequences and then the old, so if i an oldsequence and j is an old sequence then break out of this loop
					if ((i >= pDataArray->numNewFasta) && (j >= pDataArray->numNewFasta)) { break; }
					
					distCalculator->calcDist(pDataArray->alignDB.get(i), pDataArray->alignDB.get(j));
					double dist = distCalculator->getDist();
					
					if(dist <= pDataArray->cutoff){
						if (pDataArray->output == "column") { outFile << pDataArray->alignDB.get(i).getName() << ' ' << pDataArray->alignDB.get(j).getName() << ' ' << dist << endl; }
					}
					if (pDataArray->output == "lt") {  outFile  << '\t' << dist; }
				}
				
				if (pDataArray->output == "lt") { outFile << endl; }
				
				if(i % 100 == 0){
					pDataArray->m->mothurOutJustToScreen(toString(i) + "\t" + toString(time(NULL) - startTime)+"\n"); 				}
				pDataArray->count++;
			}
			pDataArray->m->mothurOutJustToScreen(toString(pDataArray->count) + "\t" + toString(time(NULL) - startTime)+"\n");
		}else{
			if(pDataArray->startLine == 0){	outFile << pDataArray->alignDB.getNumSeqs() << endl;	}
			
			for(int i=pDataArray->startLine;i<pDataArray->endLine;i++){
				
				string name = pDataArray->alignDB.get(i).getName();
				//pad with spaces to make compatible
				if (name.length() < 10) { while (name.length() < 10) {  name += " ";  } }
				
				outFile << name;
				
				for(int j=0;j<pDataArray->alignDB.getNumSeqs();j++){
					
					if (pDataArray->m->control_pressed) { delete distCalculator; outFile.close(); return 0;  }
					
					distCalculator->calcDist(pDataArray->alignDB.get(i), pDataArray->alignDB.get(j));
					double dist = distCalculator->getDist();
					
					outFile  << '\t' << dist;
				}
				
				outFile << endl; 
				
				if(i % 100 == 0){
					pDataArray->m->mothurOutJustToScreen(toString(i) + "\t" + toString(time(NULL) - startTime)+"\n"); 
				}
				pDataArray->count++;
			}
			pDataArray->m->mothurOutJustToScreen(toString(pDataArray->count) + "\t" + toString(time(NULL) - startTime)+"\n"); 
		}
		
		outFile.close();
		delete distCalculator;
		
		return 0; 
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "DistanceCommand", "MyDistThreadFunction");
		exit(1);
	}
} 
#endif

/**************************************************************************************************/
class DistanceCommand : public Command {

public:
	DistanceCommand(string);
	DistanceCommand();
	~DistanceCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "dist.seqs";			}
	string getCommandCategory()		{ return "Sequence Processing";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Schloss PD (2010). The effects of alignment quality, distance calculation method, sequence filtering, and region on the analysis of 16S rRNA gene-based studies. PLoS Comput Biol 6: e1000844. \nhttp://www.mothur.org/wiki/Dist.seqs"; }
	string getDescription()		{ return "calculate the pairwaise distances between aligned sequences"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	struct distlinePair {
		int start;
		int end;
	};

	SequenceDB alignDB;
	string countends, output, fastafile, calc, outputDir, oldfastafile, column, compress;
	int processors, numNewFasta;
	float cutoff;
	vector<int> processIDS;   //end line, processid
	vector<distlinePair> lines;
	
	bool abort;
	vector<string>  Estimators, outputNames; //holds estimators to be used
	
	void createProcesses(string, int);
	int driver(/*Dist*, SequenceDB, */int, int, string, float);
	int driver(int, int, string, string);
	bool sanityCheck();
};

#endif

/**************************************************************************************************/



