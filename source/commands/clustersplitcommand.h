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
#include "vsearchfileparser.h"
#include "opticluster.h"
#include "calculator.h"

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
    string getCommonQuestions();
	string getCitation() { return "Schloss PD, Westcott SL (2011). Assessing and improving methods used in OTU-based approaches for 16S rRNA gene sequence analysis. Appl Environ Microbiol 77:3219. \nhttp://www.mothur.org/wiki/Cluster.split"; }
	string getDescription()		{ return "splits your sequences by distance or taxonomy then clusters into OTUs"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	

private:
	vector<string> outputNames;
	string file, method, fileroot, tag, outputDir, phylipfile, columnfile, namefile, countfile, distfile, format, timing, splitmethod, taxFile, fastafile, inputDir, vsearchLocation, metricName, initialize;
	double cutoff, splitcutoff, stableMetric;
	int precision, length, processors, taxLevelCutoff, maxIters, numSingletons;
	bool  abort, large, classic, runCluster, deleteFiles, isList, cutoffNotSet, makeDist, runsensSpec, showabund; 
	
	void printData(ListVector*);
	vector<string> createProcesses(vector< map<string, string> >, set<string>&);
	int mergeLists(vector<string>, map<double, int>, ListVector*);
	map<double, int> completeListFile(vector<string>, string, set<string>&, ListVector*&);
	int createMergedDistanceFile(vector< map<string, string> >);
    string readFile(vector< map<string, string> >&);
    string printFile(string, vector< map<string, string> >&);
    int getLabels(string, set<string>& listLabels);
    bool findVsearch();
    int runSensSpec();
};

#endif

