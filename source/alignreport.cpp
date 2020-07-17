/*
 *  nastreport.cpp
 *  
 *
 *  Created by Pat Schloss on 12/19/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

#include "sequence.hpp"
#include "nast.hpp"
#include "alignment.hpp"
#include "alignreport.hpp"

/******************************************************************************************************************/

AlignReport::AlignReport() : Report() {
	try {
        fillHeaders();
	}
	catch(exception& e) {
		m->errorOut(e, "AlignReport", "AlignReport");
		exit(1);
	}
}
/**************************************************************************************************/

void AlignReport::read(ifstream& repFile){
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
        util.gobble(repFile);

        searchScore = 0;
        if(dummySearchScore != "nan"){ util.mothurConvert(dummySearchScore, searchScore); }
    }
    catch(exception& e) {
        m->errorOut(e, "AlignReport", "read");
        exit(1);
    }
    
}
/******************************************************************************************************************/
void AlignReport::fillHeaders() {
	try {
        reportHeaders.push_back("QueryName"); reportHeaders.push_back("QueryLength");
        reportHeaders.push_back("TemplateName"); reportHeaders.push_back("TemplateLength");
        reportHeaders.push_back("SearchMethod"); reportHeaders.push_back("SearchScore");
        
        reportHeaders.push_back("AlignmentMethod");
        reportHeaders.push_back("QueryStart"); reportHeaders.push_back("QueryEnd");
        reportHeaders.push_back("TemplateStart"); reportHeaders.push_back("TemplateEnd");
        reportHeaders.push_back("SearchMethod"); reportHeaders.push_back("SearchScore");

        reportHeaders.push_back("PairwiseAlignmentLength");
        reportHeaders.push_back("GapsInQuery"); reportHeaders.push_back("GapsInTemplate");
        reportHeaders.push_back("LongestInsert"); reportHeaders.push_back("SimBtwnQuery&Template");
		
	}
	catch(exception& e) {
		m->errorOut(e, "AlignReport", "fillHeaders");
		exit(1);
	}
}
/******************************************************************************************************************/

void AlignReport::print(ofstream& candidateReportFile){
	try {
		candidateReportFile << queryName << '\t' << queryLength << '\t' << templateName << '\t' << templateLength << '\t';
		candidateReportFile << searchMethod << '\t' << setprecision(2) << fixed << searchScore << '\t';

		candidateReportFile << alignmentMethod << '\t' << queryStart << "\t" << queryEnd << '\t';
		candidateReportFile << templateStart << "\t" << templateEnd << '\t';
		candidateReportFile << pairwiseAlignmentLength << '\t' << gapsInQuery << '\t' << gapsInTemplate << '\t';
		candidateReportFile << longestInsert << '\t';
		candidateReportFile << setprecision(2) << simBtwnQueryAndTemplate;
		
		candidateReportFile << endl;
		candidateReportFile.flush();
	}
	catch(exception& e) {
		m->errorOut(e, "AlignReport", "print");
		exit(1);
	}
}
/******************************************************************************************************************/

string AlignReport::getSeqReport(){
	try {
		string output = "";
		
		output += queryName + '\t' + toString(queryLength) + '\t' + templateName + '\t' + toString(templateLength) + '\t';
		
		string temp = toString(searchScore);
		int pos = temp.find_last_of('.');  //find deicmal point if their is one
		
		//if there is a decimal
		if (pos != -1) { temp = temp.substr(0, pos+3); } //set precision to 2 places
		else{	temp += ".00";	}
		
		output += searchMethod + '\t' + temp + '\t';
		output += alignmentMethod + '\t' + toString(queryStart) + "\t" + toString(queryEnd) + '\t';
		output += toString(templateStart) + "\t" + toString(templateEnd) + '\t';
		output += toString(pairwiseAlignmentLength) + '\t' + toString(gapsInQuery) + '\t' + toString(gapsInTemplate) + '\t';
		output += toString(longestInsert) + '\t';
		
		temp = toString(simBtwnQueryAndTemplate);
		pos = temp.find_last_of('.');  //find deicmal point if their is one
		
		//if there is a decimal
		if (pos != -1) { temp = temp.substr(0, pos+3); } //set precision to 2 places
		else{	temp += ".00";	}
		
		output += temp + '\n';
		
		return output;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignReport", "getSeqReport");
		exit(1);
	}
}

/******************************************************************************************************************/

void AlignReport::setCandidate(Sequence* candSeq){
	try {
		queryName = candSeq->getName();
		queryLength = candSeq->getNumBases();
	}
	catch(exception& e) {
		m->errorOut(e, "AlignReport", "setCandidate");
		exit(1);
	}
}

/******************************************************************************************************************/

void AlignReport::setTemplate(Sequence* tempSeq){
	try {
		templateName = tempSeq->getName();
		templateLength = tempSeq->getNumBases();
	}
	catch(exception& e) {
		m->errorOut(e, "AlignReport", "setTemplate");
		exit(1);
	}
}

/******************************************************************************************************************/

void AlignReport::setSearchParameters(string method, float score){
	try {
		searchMethod = method;
		searchScore = score;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignReport", "setSearchParameters");
		exit(1);
	}
}

/******************************************************************************************************************/

void AlignReport::setAlignmentParameters(string method, Alignment* align){
	try {
		alignmentMethod = method;
		
		queryStart = align->getCandidateStartPos();
		queryEnd = align->getCandidateEndPos();
		templateStart = align->getTemplateStartPos();
		templateEnd = align->getTemplateEndPos();
		pairwiseAlignmentLength = align->getPairwiseLength();

		gapsInQuery = pairwiseAlignmentLength - (queryEnd - queryStart + 1);
		gapsInTemplate = pairwiseAlignmentLength - (templateEnd - templateStart + 1);
	}
	catch(exception& e) {
		m->errorOut(e, "AlignReport", "setAlignmentParameters");
		exit(1);
	}
}
/******************************************************************************************************************/

void AlignReport::setNastParameters(Nast nast){
	try {

		longestInsert = nast.getMaxInsertLength();
		simBtwnQueryAndTemplate = nast.getSimilarityScore();
	}
	catch(exception& e) {
		m->errorOut(e, "AlignReport", "setNastParameters");
		exit(1);
	}
}

/******************************************************************************************************************/
