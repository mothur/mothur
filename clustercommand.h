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
#include "sparsematrix.hpp"

/* The cluster() command:
	The cluster command outputs a .list , .rabund and .sabund files.  
	The cluster command parameter options are method, cuttoff and precision. No parameters are required.  
	The cluster command should be in the following format: cluster(method=yourMethod, cutoff=yourCutoff, precision=yourPrecision).  
	The acceptable methods are furthest, nearest and average.  If you do not provide a method the default algorythm is furthest neighbor.  
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
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	Cluster* cluster;
	SparseMatrix* matrix;
	ListVector* list;
	RAbundVector* rabund;
	RAbundVector oldRAbund;
	ListVector oldList;

	bool abort, hard, sim;

	string method, fileroot, tag, outputDir, phylipfile, columnfile, namefile, format, distfile;
	double cutoff;
	string showabund, timing;
	int precision, length;
	ofstream sabundFile, rabundFile, listFile;

	bool print_start;
	time_t start;
	unsigned long loops;
	
	void printData(string label);
	vector<string> outputNames;
};

#endif
