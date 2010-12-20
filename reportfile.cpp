/*
 *  reportfile.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 12/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "reportfile.h"

/**************************************************************************************************/

ReportFile::ReportFile(){
	try {
		m = MothurOut::getInstance();
	}
	catch(exception& e) {
		m->errorOut(e, "ReportFile", "ReportFile");
		exit(1);
	}							
}

/**************************************************************************************************/

ReportFile::ReportFile(ifstream& repFile, string repFileName){
	try {		
		m = MothurOut::getInstance();
		
		m->openInputFile(repFileName, repFile);
		m->getline(repFile);
	}
	catch(exception& e) {
		m->errorOut(e, "ReportFile", "ReportFile");
		exit(1);
	}							
}


/**************************************************************************************************/

ReportFile::ReportFile(ifstream& repFile){
	try {
		
		m = MothurOut::getInstance();
		
		repFile >> queryName;
		repFile >> queryLength;
		repFile >> templateName;
		repFile >> templateLength;
		repFile >> searchMethod;
		repFile >> searchScore;
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

		m->gobble(repFile);		
	}
	catch(exception& e) {
		m->errorOut(e, "ReportFile", "ReportFile");
		exit(1);
	}							
	
}

/**************************************************************************************************/
