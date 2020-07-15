/*
 *  reportfile.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 12/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "mothurout.h"
#include "seqerrorreportfile.h"

/**************************************************************************************************/

SeqErrorReportFile::SeqErrorReportFile() : Report() {}

/**************************************************************************************************/

int SeqErrorReportFile::read(ifstream& repFile){
	try {
        
		repFile >> queryName;
		repFile >> queryLength;
		repFile >> templateName;
		repFile >> templateLength;
		repFile >> searchMethod;
		repFile >> dummySearchScore;
		repFile >> alignmentMethod;
		repFile >> queryStart;
		repFile >> queryEnd;
		repFile >> templateStart;
		repFile >> templateEnd;
		repFile >> pairwiseAlignmentLength;
		repFile >> gapsInQuery;
		repFile >> gapsInTemplate;
		repFile >> longestInsert;
		repFile >> simBtwnQueryAndTemplate;

		if(dummySearchScore != "nan"){
			istringstream stream(dummySearchScore);
			stream >> searchScore;
		}
		else{
			searchScore = 0;
		}
		
		util.gobble(repFile);
        
        return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SeqErrorReportFile", "read");
		exit(1);
	}							
	
}

/**************************************************************************************************/
