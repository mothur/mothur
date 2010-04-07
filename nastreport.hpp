#ifndef NASTREPORT_HPP
#define NASTREPORT_HPP


/*
 *  nastreport.hpp
 *  
 *
 *  Created by Pat Schloss on 12/19/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mothur.h"

/******************************************************************************************************************/

class NastReport {

public:
	NastReport(string);
	NastReport();
	~NastReport();
	void setCandidate(Sequence*);
	void setTemplate(Sequence*);
	void setSearchParameters(string, float);
	void setAlignmentParameters(string, Alignment*);
	void setNastParameters(Nast);
	void print();
	string getReport();
	string getHeaders();
	
private:
	string queryName;
	string output;
	int queryLength;
	string templateName;
	int templateLength;
	string searchMethod;
	float searchScore;
	string alignmentMethod;
	int candidateStartPosition, candidateEndPosition;
	int templateStartPosition, templateEndPosition;

	int pairwiseAlignmentLength;
	int longestInsert;
	int totalGapsInQuery, totalGapsInTemplate;
	float similarityToTemplate;
	ofstream candidateReportFile;
};

/******************************************************************************************************************/

#endif
