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

		SplitMatrix(string, string, string, string, float, float, int, bool, string, bool); //fastafile, namefile, countfile, taxFile, taxcutoff, cutoff, processors, classic, outputDir, usingVsearchToCLuster
		
        ~SplitMatrix() {}
    
		vector< map<string, string> > getDistanceFiles();  //returns map of distance files -> namefile sorted by distance file size
		string getSingletonNames() { return singleton; } //returns namesfile or countfile containing singletons
        //long long getNumSingleton() { return numSingleton; } //returns namesfile containing singletons
	
	private:
		MothurOut* m;
        Utils util;

		string distFile, namefile, singleton,  taxFile, fastafile, outputDir, countfile;
		vector< map< string, string> > dists;
		float cutoff, distCutoff;
		bool classic, usingVsearchToCLuster; 
        int processors;

		void splitClassify();
		int createDistanceFilesFromTax(vector<vector<string> >&, vector<string>);
        int createFastaFilesFromTax(vector<vector<string> >&, vector<string>);
};

/******************************************************/

#endif

