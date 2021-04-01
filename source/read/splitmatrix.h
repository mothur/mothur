#ifndef SPLITMATRIX_H
#define SPLITMATRIX_H
/*
 *  splitmatrix.h
 *  Mothur
 *
 *  Created by westcott on 5/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "mothur.h"
#include "mothurout.h"
#include "utils.hpp"
#include "counttable.h"

/******************************************************/

class SplitMatrix  {
	
	public:

		SplitMatrix(string, string, string, string, float, string, bool); //column formatted distance file, namesfile, countfile, cutoff, method, large
		SplitMatrix(string, string, string, string, float, float, string, int, bool, string, string); //fastafile, namefile, countfile, taxFile, taxcutoff, cutoff, method, processors, classic, outputDir, ("fasta" or "distance")
		
		~SplitMatrix();
    
		int split();
		vector< map<string, string> > getDistanceFiles();  //returns map of distance files -> namefile sorted by distance file size
		string getSingletonNames() { return singleton; } //returns namesfile containing singletons
        long long getNumSingleton() { return numSingleton; } //returns namesfile containing singletons
	
	private:
		MothurOut* m;
        Utils util;

		string distFile, namefile, singleton, method, taxFile, fastafile, outputDir, countfile, outputType;
		vector< map< string, string> > dists;
		float cutoff, distCutoff;
		bool large, classic;
        int processors;
        long long numSingleton;
				
		int splitDistance();
		int splitClassify();
		int splitDistanceLarge();
		int splitDistanceRAM();
		int splitNames(map<string, int>& groups, int, vector<string>&);
        int splitNamesVsearch(map<string, int>& groups, int, vector<string>&);
		int splitDistanceFileByTax(map<string, int>&, int);
		int createDistanceFilesFromTax(map<string, int>&, int);
};

/******************************************************/

#endif

