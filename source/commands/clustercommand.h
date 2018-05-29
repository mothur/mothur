#ifndef CLUSTERCOMMAND_H
#define CLUSTERCOMMAND_H
/*
 *  clustercommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "rabundvector.hpp"
#include "sabundvector.hpp"
#include "listvector.hpp"
#include "cluster.hpp"
#include "counttable.h"
#include "vsearchfileparser.h"
#include "clusterdoturcommand.h"
#include "opticluster.h"
#include "optimatrix.h"
#include "calculator.h"


/* The cluster() command:
	The cluster command outputs a .list , .rabund and .sabund files.  
	The cluster command parameter options are method, cuttoff and precision. No parameters are required.  
	The cluster command should be in the following format: cluster(method=yourMethod, cutoff=yourCutoff, precision=yourPrecision).  
	The acceptable methods are furthest, nearest and average.  If you do not provide a method the default algorithm is furthest neighbor.  
	The cluster() command outputs three files *.list, *.rabund, and *.sabund.   */


class ClusterCommand : public Command {
	
public:
	ClusterCommand(string);
	ClusterCommand();
	~ClusterCommand();
	
	vector<string> setParameters();
	string getCommandName()			{ return "cluster";		}
	string getCommandCategory()		{ return "Clustering";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Schloss PD, Westcott SL (2011). Assessing and improving methods used in OTU-based approaches for 16S rRNA gene sequence analysis. Appl Environ Microbiol 77:3219.\nSchloss PD, Handelsman J (2005). Introducing DOTUR, a computer program for defining operational taxonomic units and estimating species richness. Appl Environ Microbiol 71: 1501-6.\nhttp://www.mothur.org/wiki/Cluster"; }
	string getDescription()		{ return "cluster your sequences into OTUs using a distance matrix"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	Cluster* cluster;
	SparseDistanceMatrix* matrix;
	ListVector* list;
	RAbundVector* rabund;
	RAbundVector oldRAbund;
	ListVector oldList;

	bool abort, sim, cutOffSet;

	string method, fileroot, tag, outputDir, phylipfile, columnfile, namefile, format, distfile, countfile, fastafile, inputDir, vsearchLocation, metric, initialize;
	double cutoff, stableMetric;
    float adjust;
	string showabund, timing, metricName;
	int precision, length, maxIters, processors;
	ofstream sabundFile, rabundFile, listFile;

	bool print_start;
	time_t start;
	unsigned long loops;
	
	void printData(string label, map<string, int>&, bool&);
	vector<string> outputNames;
    
    int createRabund(CountTable*&, ListVector*&, RAbundVector*&);
    int vsearchDriver(string, string, string);
    int runVsearchCluster();
    int runOptiCluster();
    int runMothurCluster();
    int runUniqueCluster();
};

#endif
