//
//  getmetacommunitycommand.h
//  Mothur
//
//  Created by SarahsWork on 4/9/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef Mothur_getmetacommunitycommand_h
#define Mothur_getmetacommunitycommand_h

#include "command.hpp"
#include "inputdata.h"
#include "qFinderDMM.h"

/**************************************************************************************************/

class GetMetaCommunityCommand : public Command {
public:
    GetMetaCommunityCommand(string);
    GetMetaCommunityCommand();
    ~GetMetaCommunityCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "get.metacommunity";		}
    string getCommandCategory()		{ return "OTU-Based Approaches";         }
    
    string getOutputPattern(string);
    
	string getHelpString();
    string getCitation() { return "Holmes I, Harris K, Quince C (2012) Dirichlet Multinomial Mixtures: Generative Models for Microbial Metagenomics. PLoS ONE 7(2): e30126. doi:10.1371/journal.pone.0030126 http://www.mothur.org/wiki/Get.metacommunity"; }
    string getDescription()		{ return "Assigns samples to bins using a Dirichlet multinomial mixture model"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    struct linePair {
		unsigned long long start;
		unsigned long long end;
		linePair(unsigned long long i, unsigned long long j) : start(i), end(j) {}
	};
    bool abort, allLines;
    string outputDir;
    vector<string> outputNames;
    string sharedfile;
    int minpartitions, maxpartitions, optimizegap, processors;
    vector<string> Groups;
    set<string> labels;
    
    int processDriver(vector<SharedRAbundVector*>&, vector<int>&, string, vector<string>, vector<string>, vector<string>, int);
    int createProcesses(vector<SharedRAbundVector*>&);
    vector<double> generateDesignFile(int, map<string,string>);
    int generateSummaryFile(int, map<string,string>);

};

/**************************************************************************************************/
struct summaryData {
    
    string name;
    double refMean, difference;
    vector<double> partMean, partLCI, partUCI;
    
};
/**************************************************************************************************/

struct metaCommunityData {
    vector<SharedRAbundVector*> thislookup;
	MothurOut* m;
	string outputFileName;
    vector<string> relabunds, matrix, outputNames;
    int minpartitions, maxpartitions, optimizegap;
    vector<int> parts;
    int minPartition;
	
	metaCommunityData(){}
	metaCommunityData(vector<SharedRAbundVector*> lu, MothurOut* mout, vector<int> dp, string fit, vector<string> rels, vector<string> mat, int minp, int maxp, int opg) {
		m = mout;
        thislookup = lu;
        parts = dp;
        outputFileName = fit;
        relabunds = rels;
        matrix = mat;
        minpartitions = minp;
        maxpartitions = maxp;
        optimizegap = opg;
        minPartition = 0;
	}
};
/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyMetaCommunityThreadFunction(LPVOID lpParam){
	metaCommunityData* pDataArray;
	pDataArray = (metaCommunityData*)lpParam;
	
	try {
        
        double minLaplace = 1e10;
        
		ofstream fitData;
		pDataArray->m->openOutputFile(pDataArray->outputFileName, fitData);
        fitData.setf(ios::fixed, ios::floatfield);
        fitData.setf(ios::showpoint);
        cout.setf(ios::fixed, ios::floatfield);
        cout.setf(ios::showpoint);
        
        vector< vector<int> > sharedMatrix;
        for (int i = 0; i < pDataArray->thislookup.size(); i++) { sharedMatrix.push_back(pDataArray->thislookup[i]->getAbundances()); }
        
        pDataArray->m->mothurOut("K\tNLE\t\tlogDet\tBIC\t\tAIC\t\tLaplace\n");
        fitData << "K\tNLE\tlogDet\tBIC\tAIC\tLaplace" << endl;
        
        for(int i=0;i<pDataArray->parts.size();i++){
            
            int numPartitions = pDataArray->parts[i];
            
            if (pDataArray->m->debug) { pDataArray->m->mothurOut("[DEBUG]: running partition " + toString(numPartitions) + "\n"); }
            
            if (pDataArray->m->control_pressed) { break; }
            
            qFinderDMM* findQ = new qFinderDMM(sharedMatrix, numPartitions);
            
            if (pDataArray->m->debug) { pDataArray->m->mothurOut("[DEBUG]: done finding Q " + toString(numPartitions) + "\n"); }
            
            double laplace = findQ->getLaplace();
            pDataArray->m->mothurOut(toString(numPartitions) + '\t');
            cout << setprecision (2) << findQ->getNLL() << '\t' << findQ->getLogDet() << '\t';
            pDataArray->m->mothurOutJustToLog(toString(findQ->getNLL()) + '\t' + toString(findQ->getLogDet()) + '\t');
            cout << findQ->getBIC() << '\t' << findQ->getAIC() << '\t' << laplace;
            pDataArray->m->mothurOutJustToLog(toString(findQ->getBIC()) + '\t' + toString(findQ->getAIC()) + '\t' + toString(laplace));
            
            fitData << numPartitions << '\t';
            fitData << setprecision (2) << findQ->getNLL() << '\t' << findQ->getLogDet() << '\t';
            fitData << findQ->getBIC() << '\t' << findQ->getAIC() << '\t' << laplace << endl;
            
            if(laplace < minLaplace){
                pDataArray->minPartition = numPartitions;
                minLaplace = laplace;
                pDataArray->m->mothurOut("***");
            }
            pDataArray->m->mothurOutEndLine();
            
            pDataArray->outputNames.push_back(pDataArray->relabunds[i]);
            pDataArray->outputNames.push_back(pDataArray->matrix[i]);
            
            findQ->printZMatrix(pDataArray->matrix[i], pDataArray->m->getGroups());
            findQ->printRelAbund(pDataArray->relabunds[i], pDataArray->m->currentBinLabels);
            
            if(pDataArray->optimizegap != -1 && (numPartitions - pDataArray->minPartition) >= pDataArray->optimizegap && numPartitions >= pDataArray->minpartitions){ break; }
            
            delete findQ;
        }
        fitData.close();
        
        //minPartition = 4;
        
        if (pDataArray->m->control_pressed) { return 0; }
        
        return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "GetMetaCommunityCommand", "MyMetaCommunityThreadFunction");
		exit(1);
	}
}
#endif



#endif
