/*
 *  bayesian.cpp
 *  Mothur
 *
 *  Created by westcott on 11/3/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "bayesian.h"

/**************************************************************************************************/
Bayesian::Bayesian(string tfile, string tempFile, string method, int kmerSize, int gapOpen, int gapExtend, int match, int misMatch) : 
Classify(tfile, tempFile, method, kmerSize, gapOpen, gapExtend, match, misMatch) {}
/**************************************************************************************************/
string Bayesian::getTaxonomy(Sequence* seq) {
	try {
		string tax;
		
				
		return tax;	
	}
	catch(exception& e) {
		errorOut(e, "Bayesian", "getTaxonomy");
		exit(1);
	}
}
/**************************************************************************************************/

