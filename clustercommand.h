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

#include <Carbon/Carbon.h>
#include <iostream>
#include <fstream>
#include "command.hpp"
#include "rabundvector.hpp"
#include "listvector.hpp"
#include "cluster.hpp"
#include "sparsematrix.hpp"

/* The cluster() command:
	The cluster command can only be executed after a successful read.phylip or read.column command.   
	The cluster command outputs a .list , .rabund and .sabund files.  
	The cluster command parameter options are method, cuttoff and precision. No parameters are required.  
	The cluster command should be in the following format: cluster(method=yourMethod, cutoff=yourCutoff, precision=yourPrecision).  
	The acceptable methods are furthest, nearest and average.  If you do not provide a method the default algorythm is furthest neighbor.  
	The cluster() command outputs three files *.list, *.rabund, and *.sabund.   */


class GlobalData;

class ClusterCommand : public Command {
	
public:
	ClusterCommand();	
	~ClusterCommand();
	int execute();	
	
private:
	GlobalData* globaldata;
	Cluster* cluster;
	SparseMatrix* matrix;
	ListVector* list;
	RAbundVector* rabund;
	RAbundVector oldRAbund;
	ListVector oldList;
	
	string method, fileroot, tag;
	double cutoff;
	int precision, length;
	ofstream sabundFile, rabundFile, listFile;
	
	void printData(string label);
};

#endif