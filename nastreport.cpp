/*
 *  nastreport.cpp
 *  
 *
 *  Created by Pat Schloss on 12/19/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

using namespace std;

#include "sequence.hpp"
#include "nast.hpp"
#include "alignment.hpp"
#include "nastreport.hpp"

/******************************************************************************************************************/

NastReport::NastReport(string candidateReportFName) {
	openOutputFile(candidateReportFName, candidateReportFile);
	
	candidateReportFile << "QueryName\tQueryLength\tTemplateName\tTemplateLength\t";
	candidateReportFile << "SearchMethod\tSearchScore\t";
	candidateReportFile << "AlignmentMethod\tQueryStart\tQueryEnd\tTemplateStart\tTemplateEnd\t";
	candidateReportFile << "PairwiseAlignmentLength\tGapsInQuery\tGapsInTemplate\t";
	candidateReportFile << "LongestInsert\t";
	candidateReportFile << "SimBtwnQuery&Template" << endl;
}

/******************************************************************************************************************/

void NastReport::print(){
	
	candidateReportFile << queryName << '\t' << queryLength << '\t' << templateName << '\t' << templateLength << '\t';
	candidateReportFile << searchMethod << '\t' << setprecision(2) << fixed << searchScore << '\t';

	candidateReportFile << alignmentMethod << '\t' << candidateStartPosition << "\t" << candidateEndPosition << '\t';
	candidateReportFile << templateStartPosition << "\t" << templateEndPosition << '\t';
	candidateReportFile << pairwiseAlignmentLength << '\t' << totalGapsInQuery << '\t' << totalGapsInTemplate << '\t';
	candidateReportFile << longestInsert << '\t';;
	candidateReportFile << setprecision(2) << similarityToTemplate;
	
	candidateReportFile << endl;
	candidateReportFile.flush();
}

/******************************************************************************************************************/

void NastReport::setCandidate(Sequence* candSeq){ 
	queryName = candSeq->getName();
	queryLength = candSeq->getNumBases();
}

/******************************************************************************************************************/

void NastReport::setTemplate(Sequence* tempSeq){ 
	templateName = tempSeq->getName();
	templateLength = tempSeq->getNumBases();
}

/******************************************************************************************************************/

void NastReport::setSearchParameters(string method, float score){
	searchMethod = method;
	searchScore = score;
}

/******************************************************************************************************************/

void NastReport::setAlignmentParameters(string method, Alignment* align){
	alignmentMethod = method;
	
	candidateStartPosition = align->getCandidateStartPos();
	candidateEndPosition = align->getCandidateEndPos();
	templateStartPosition = align->getTemplateStartPos();
	templateEndPosition = align->getTemplateEndPos();
	pairwiseAlignmentLength = align->getPairwiseLength();

	totalGapsInQuery = pairwiseAlignmentLength - (candidateEndPosition - candidateStartPosition + 1);
	totalGapsInTemplate = pairwiseAlignmentLength - (templateEndPosition - templateStartPosition + 1);
}

/******************************************************************************************************************/

void NastReport::setNastParameters(Nast nast){

	longestInsert = nast.getMaxInsertLength();
	similarityToTemplate = nast.getSimilarityScore();
	
}

/******************************************************************************************************************/
