#ifndef CLUSTERSPLITCOMMAND_H
#define CLUSTERSPLITCOMMAND_H

/*
 *  clustersplitcommand.h
 *  Mothur
 *
 *  Created by westcott on 5/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */
 
#include "command.hpp"
#include "rabundvector.hpp"
#include "sabundvector.hpp"
#include "listvector.hpp"
#include "cluster.hpp"
#include "sparsedistancematrix.h"
#include "readcluster.h"
#include "splitmatrix.h"
#include "readphylip.h"
#include "readcolumn.h"
#include "readmatrix.hpp"
#include "inputdata.h"
#include "clustercommand.h"
#include "clusterclassic.h"

class ClusterSplitCommand : public Command {
	
public:
	ClusterSplitCommand(string);
	ClusterSplitCommand();
	~ClusterSplitCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "cluster.split";		}
	string getCommandCategory()		{ return "Clustering";			}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Schloss PD, Westcott SL (2011). Assessing and improving methods used in OTU-based approaches for 16S rRNA gene sequence analysis. Appl Environ Microbiol 77:3219. \nhttp://www.mothur.org/wiki/Cluster.split"; }
	string getDescription()		{ return "splits your sequences by distance or taxonomy then clusters into OTUs"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	

private:
	vector<int> processIDS;   //processid
	vector<string> outputNames;
	
	string file, method, fileroot, tag, outputDir, phylipfile, columnfile, namefile, countfile, distfile, format, showabund, timing, splitmethod, taxFile, fastafile, inputDir, vsearchLocation, metric;
	double cutoff, splitcutoff, stableMetric;
	int precision, length, processors, taxLevelCutoff, maxIters;
	bool print_start, abort, hard, large, classic, runCluster, deleteFiles, isList, cutoffNotSet;
	time_t start;
	ofstream outList, outRabund, outSabund;
	
	void printData(ListVector*);
	vector<string> createProcesses(vector< map<string, string> >, set<string>&);
	vector<string> cluster(vector< map<string, string> >, set<string>&);
    string clusterFile(string, string, set<string>&, double&);
    string clusterClassicFile(string, string, set<string>&, double&);
	int mergeLists(vector<string>, map<float, int>, ListVector*);
	map<float, int> completeListFile(vector<string>, string, set<string>&, ListVector*&);
	int createMergedDistanceFile(vector< map<string, string> >);
    int createRabund(CountTable*& ct, ListVector*& list, RAbundVector*& rabund);
    string readFile(vector< map<string, string> >&);
    string printFile(string, vector< map<string, string> >&);
    int getLabels(string, set<string>& listLabels);
    bool findVsearch();
    int vsearchDriver(string, string, string, double);
    string runVsearchCluster(string, string, set<string>&, double&);
    string runOptiCluster(string, string, set<string>&, double&);
};

/////////////////not working for Windows////////////////////////////////////////////////////////////
// getting an access violation error.  This is most likely caused by the 
// threads stepping on eachother's structures, as I can run the thread function and the cluster fuction 
// in separately without errors occuring.  I suspect it may be in the use of the
// static class mothurOut, but I can't pinpoint the problem.  All other objects are made new
// within the thread.  MothurOut is used by almost all the classes in mothur, so if this was 
// really the cause I would expect to see all the windows threaded commands to have issues, but not 
// all do. So far, shhh.flows and trim.flows have similiar problems. Other thoughts, could it have 
// anything to do with mothur's use of copy constructors in many of our data structures. ie. listvector 
// is copied by nameassignment and passed to read which passes to the thread?  -westcott 2-8-12
////////////////////////////////////////////////////////////////////////////////////////////////////
/**************************************************************************************************
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct clusterData {
	set<string> labels;
	vector < map<string, string> > distNames; 
	string method; 
    MothurOut* m;
	double cutoff, precision;
    string tag, outputDir;
    vector<string> listFiles;
    bool hard;
    int length, threadID;
	
	
	clusterData(){}
	clusterData(vector < map<string, string> > dv, MothurOut* mout, double cu, string me, string ou, bool hd, double pre, int len, int th) {
		distNames = dv;
		m = mout;
		cutoff = cu;
        method = me;
		outputDir = ou;
        hard = hd;
        precision = pre;
        length = len;
        threadID = th;
	}
};

/**************************************************************************************************
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyClusterThreadFunction(LPVOID lpParam){ 
	clusterData* pDataArray;
	pDataArray = (clusterData*)lpParam;
	
	try {
		cout << "starting " << endl;		
		
		double smallestCutoff = pDataArray->cutoff;
		
		//cluster each distance file
		for (int i = 0; i < pDataArray->distNames.size(); i++) {
            
            Cluster* mycluster = NULL;
            SparseMatrix* mymatrix = NULL;
            ListVector* mylist = NULL;
            ListVector myoldList;
            RAbundVector* myrabund = NULL;
                        
			if (pDataArray->m->control_pressed) { break; }
			
			string thisNamefile = pDataArray->distNames[i].begin()->second;
			string thisDistFile = pDataArray->distNames[i].begin()->first;
            cout << thisNamefile << '\t' << thisDistFile << endl;	
			pDataArray->m->mothurOutEndLine(); pDataArray->m->mothurOut("Reading " + thisDistFile); pDataArray->m->mothurOutEndLine();
			
			ReadMatrix* myread = new ReadColumnMatrix(thisDistFile); 	
			myread->setCutoff(pDataArray->cutoff);
			NameAssignment* mynameMap = new NameAssignment(thisNamefile);
			mynameMap->readMap();
            cout << "done reading " << thisNamefile << endl;  
			myread->read(mynameMap);
			cout << "done reading " << thisDistFile << endl;  
			if (pDataArray->m->control_pressed) {  delete myread; delete mynameMap; break; }
            
			mylist = myread->getListVector();
			myoldList = *mylist;
			mymatrix = myread->getMatrix();
            cout << "here" << endl;	
			delete myread; myread = NULL;
			delete mynameMap; mynameMap = NULL;
			
            pDataArray->m->mothurOutEndLine(); pDataArray->m->mothurOut("Clustering " + thisDistFile); pDataArray->m->mothurOutEndLine();
            
			myrabund = new RAbundVector(mylist->getRAbundVector());
			 cout << "here" << endl;	
			//create cluster
			if (pDataArray->method == "furthest")	{	mycluster = new CompleteLinkage(myrabund, mylist, mymatrix, pDataArray->cutoff, pDataArray->method); }
			else if(pDataArray->method == "nearest"){	mycluster = new SingleLinkage(myrabund, mylist, mymatrix, pDataArray->cutoff, pDataArray->method); }
			else if(pDataArray->method == "average"){	mycluster = new AverageLinkage(myrabund, mylist, mymatrix, pDataArray->cutoff, pDataArray->method);	}
			pDataArray->tag = mycluster->getTag();
             cout << "here" << endl;	
			if (pDataArray->outputDir == "") { pDataArray->outputDir += pDataArray->m->hasPath(thisDistFile); }
			string fileroot = pDataArray->outputDir + pDataArray->m->getRootName(pDataArray->m->getSimpleName(thisDistFile));
			 cout << "here" << endl;	
			ofstream listFile;
			pDataArray->m->openOutputFile(fileroot+ pDataArray->tag + ".list",	listFile);
             cout << "here" << endl;	
			pDataArray->listFiles.push_back(fileroot+ pDataArray->tag + ".list");
            
			float previousDist = 0.00000;
			float rndPreviousDist = 0.00000;
			
			myoldList = *mylist;
        
			bool print_start = true;
			int start = time(NULL);
			double saveCutoff = pDataArray->cutoff;
            
			while (mymatrix->getSmallDist() < pDataArray->cutoff && mymatrix->getNNodes() > 0){
                
				if (pDataArray->m->control_pressed) { //clean up
					delete mymatrix; delete mylist;	delete mycluster; delete myrabund;
					listFile.close();
					for (int i = 0; i < pDataArray->listFiles.size(); i++) {	pDataArray->m->mothurRemove(pDataArray->listFiles[i]); 	}
					pDataArray->listFiles.clear(); break;
				}
                
				mycluster->update(saveCutoff);
                
				float dist = mymatrix->getSmallDist();
				float rndDist;
				if (pDataArray->hard) {
					rndDist = pDataArray->m->ceilDist(dist, pDataArray->precision); 
				}else{
					rndDist = pDataArray->m->roundDist(dist, pDataArray->precision); 
				}
                
				if(previousDist <= 0.0000 && dist != previousDist){
					myoldList.setLabel("unique");
					myoldList.print(listFile);
					if (pDataArray->labels.count("unique") == 0) {  pDataArray->labels.insert("unique");  }
				}
				else if(rndDist != rndPreviousDist){
					myoldList.setLabel(toString(rndPreviousDist,  pDataArray->length-1));
					myoldList.print(listFile);
					if (pDataArray->labels.count(toString(rndPreviousDist,  pDataArray->length-1)) == 0) { pDataArray->labels.insert(toString(rndPreviousDist,  pDataArray->length-1)); }
				}
               	
				previousDist = dist;
				rndPreviousDist = rndDist;
				myoldList = *mylist;
			}
            
             cout << "here2" << endl;	
			if(previousDist <= 0.0000){
				myoldList.setLabel("unique");
				myoldList.print(listFile);
				if (pDataArray->labels.count("unique") == 0) { pDataArray->labels.insert("unique"); }
			}
			else if(rndPreviousDist<pDataArray->cutoff){
				myoldList.setLabel(toString(rndPreviousDist,  pDataArray->length-1));
				myoldList.print(listFile);
				if (pDataArray->labels.count(toString(rndPreviousDist,  pDataArray->length-1)) == 0) { pDataArray->labels.insert(toString(rndPreviousDist,  pDataArray->length-1)); }
			}
            
			delete mymatrix; delete mylist;	delete mycluster; delete myrabund; 
            mymatrix = NULL; mylist = NULL; mycluster = NULL; myrabund = NULL;
			listFile.close();
			
			if (pDataArray->m->control_pressed) { //clean up
				for (int i = 0; i < pDataArray->listFiles.size(); i++) {	pDataArray->m->mothurRemove(pDataArray->listFiles[i]); 	}
				pDataArray->listFiles.clear(); break;
			}
			 cout << "here3" << endl;	
			pDataArray->m->mothurRemove(thisDistFile);
			pDataArray->m->mothurRemove(thisNamefile);
			 cout << "here4" << endl;	
			if (saveCutoff != pDataArray->cutoff) { 
				if (pDataArray->hard)	{  saveCutoff = pDataArray->m->ceilDist(saveCutoff, pDataArray->precision);	}
				else		{	saveCutoff = pDataArray->m->roundDist(saveCutoff, pDataArray->precision);  }
                
				pDataArray->m->mothurOut("Cutoff was " + toString(pDataArray->cutoff) + " changed cutoff to " + toString(saveCutoff)); pDataArray->m->mothurOutEndLine();  
			}
			 cout << "here5" << endl;	
			if (saveCutoff < smallestCutoff) { smallestCutoff = saveCutoff;  }
		}
		
		pDataArray->cutoff = smallestCutoff;
		
		return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "ClusterSplitCommand", "MyClusterThreadFunction");
		exit(1);
	}
} 
#endif

*/


#endif

