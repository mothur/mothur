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

/******************************************************/

class SplitMatrix  {
	
	public:

		SplitMatrix(string, string, string, string, float, string, bool); //column formatted distance file, namesfile, countfile, cutoff, method, large
		SplitMatrix(string, string, string, string, float, float, string, int, bool, string); //fastafile, namefile, countfile, taxFile, taxcutoff, cutoff, method, processors, classic, outputDir
		
		~SplitMatrix();
		int split();
		vector< map<string, string> > getDistanceFiles();  //returns map of distance files -> namefile sorted by distance file size
		string getSingletonNames() { return singleton; } //returns namesfile containing singletons
	
	private:
		MothurOut* m;

		string distFile, namefile, singleton, method, taxFile, fastafile, outputDir, countfile;
		vector< map< string, string> > dists;
		float cutoff, distCutoff;
		bool large, classic;
		int processors;
				
		int splitDistance();
		int splitClassify();
		int splitDistanceLarge();
		int splitDistanceRAM();
		int splitNames(map<string, int>& groups, int, vector<string>&);
		int splitDistanceFileByTax(map<string, int>&, int);
		int createDistanceFilesFromTax(map<string, int>&, int);
};

/******************************************************/

#endif

