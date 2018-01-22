#ifndef READCLUSTER_H
#define READCLUSTER_H
/*
 *  readcluster.h
 *  Mothur
 *
 *  Created by westcott on 10/28/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */


#include "mothur.h"
#include "nameassignment.hpp"
#include "listvector.hpp"
#include "counttable.h"
#include "utils.hpp"


/******************************************************/

class ReadCluster {
	
public:
	ReadCluster(string, float, string, bool);
	~ReadCluster();
	int read(NameAssignment*&);
    int read(CountTable*&);
	string getOutputFile() { return OutPutFile; }
	void setFormat(string f) { format = f;	}
	ListVector* getListVector()		{	return list;	}
	
private:
	string distFile, outputDir;
	string OutPutFile, format;
	ListVector* list;
	float cutoff;
	MothurOut* m;
	bool sortWanted;
    Utils util;
	
	int convertPhylip2Column(NameAssignment*&);
    int convertPhylip2Column(CountTable*&);
};

/******************************************************/

#endif

