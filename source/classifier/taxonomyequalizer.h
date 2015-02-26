#ifndef TAXONOMYEQUALIZER_H
#define TAXONOMYEQUALIZER_H


/*
 *  taxonomyequalizer.h
 *  Mothur
 *
 *  Created by westcott on 11/20/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "mothurout.h"

//reads in taxonomy file and makes all the taxonomies the same length 
//by appending the last taxon to a given taxonomy as many times as needed to 
//make it as long as the longest taxonomy in the file 

/**************************************************************************************************/

class TaxEqualizer  {
	
public:
	TaxEqualizer(string, int, string);
	~TaxEqualizer() {};
	
	string getEqualizedTaxFile()	{	return equalizedFile;	}
	int getHighestLevel()			{	return highestLevel;	}
	
	
private:
	string equalizedFile, testTax, outputDir;
	bool containsConfidence;
	int cutoff, highestLevel;
	map<string, int> seqLevels;  //maps name to level of taxonomy
	
	int getHighestLevel(ifstream&);  //scans taxonomy file to find taxonomy with highest level
	void extendTaxonomy(string, string&, int);  //name, taxonomy, desired level
	void truncateTaxonomy(string, string&, int);  //name, taxonomy, desired level
	MothurOut* m;
	
};

/**************************************************************************************************/


#endif


