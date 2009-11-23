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

//reads in taxonomy file and makes all the taxonomies the same length 
//by appending the last taxon to a given taxonomy as many times as needed to 
//make it as long as the longest taxonomy in the file 

/**************************************************************************************************/

class TaxEqualizer  {
	
public:
	TaxEqualizer(string, int);
	~TaxEqualizer() {};
	
	string getEqualizedTaxFile()	{  return equalizedFile;	}
	
	
private:
	string equalizedFile, testTax;
	bool containsConfidence;
	int cutoff;
	map<string, int> seqLevels;  //maps name to level of taxonomy
	
	int getHighestLevel(ifstream&);  //scans taxonomy file to find taxonomy with highest level
	void extendTaxonomy(string, string&, int);  //name, taxonomy, desired level
	void truncateTaxonomy(string, string&, int);  //name, taxonomy, desired level
	void removeConfidences(string&);  //removes the confidence limits on the taxon 

	
};

/**************************************************************************************************/


#endif


