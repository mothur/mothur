/*
 *  phylotype.cpp
 *  Mothur
 *
 *  Created by westcott on 11/3/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "phylotype.h"

/**************************************************************************************************/
PhyloType::PhyloType(string tfile, string tempFile, string method, int kmerSize, int gapOpen, int gapExtend, int match, int misMatch) 
: Classify(tfile, tempFile, method, kmerSize, gapOpen, gapExtend, match, misMatch)  {}
/**************************************************************************************************/
string PhyloType::getTaxonomy(Sequence* seq) {
	try {
		string tax;
		
		//use database to find closest seq
		vector<int> closest = database->findClosestSequences(seq, 1);
		
		//find that sequences taxonomy in map
		it = taxonomy.find(names[closest[0]]);
		
		//is this sequence in the taxonomy file
		if (it == taxonomy.end()) { //error not in file
			mothurOut("Error: sequence " + names[closest[0]] + " is not in the taxonomy file.  It is the closest match to sequence " + seq->getName() + ". " + seq->getName() + " will be disregarded."); mothurOutEndLine();
			tax = "bad seq";
		}else {
			tax = it->second;
		}
		
		return tax;	
	}
	catch(exception& e) {
		errorOut(e, "PhyloType", "getTaxonomy");
		exit(1);
	}
}
/**************************************************************************************************/

