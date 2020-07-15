#ifndef REPORTFILE
#define REPORTFILE

/*
 *  reportfile.h
 *  Mothur
 *
 *  Created by Pat Schloss on 7/12/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "utils.hpp"
#include "report.hpp"

/**************************************************************************************************/

class SeqErrorReportFile : public Report {
    
public:
    SeqErrorReportFile();
    ~SeqErrorReportFile() {}
    
	int read(ifstream&);
	
	string getQueryName()				{	return queryName;				}
	string getTemplateName()			{	return templateName;			}
	string getSearchMethod()			{	return searchMethod;			}
	string getAlignmentMethod()			{	return alignmentMethod;			}
	
	int getQueryLength()				{	return queryLength;				}
	int getTemplateLength()				{	return templateLength;			}
	int getQueryStart()					{	return queryStart;				}
	int getQueryEnd()					{	return queryEnd;				}
	int getTemplateStart()				{	return templateStart;			}
	int getTemplateEnd()				{	return templateEnd;				}
	int getPairwiseAlignmentLength()	{	return pairwiseAlignmentLength;	}
	int getGapsInQuery()				{	return gapsInQuery;				}
	int getGapsInTemplate()				{	return gapsInTemplate;			}
	int getLongestInsert()				{	return longestInsert;			}
	
	float getSearchScore()				{	return searchScore;				}
	float getSimBtwnQueryAndTemplate()	{	return simBtwnQueryAndTemplate;	}
	

private:
    
	string queryName, templateName, searchMethod, alignmentMethod, dummySearchScore;
	int queryLength, templateLength, queryStart, queryEnd, templateStart, templateEnd, pairwiseAlignmentLength, gapsInQuery, gapsInTemplate, longestInsert;
	float searchScore, simBtwnQueryAndTemplate;
	
};

/**************************************************************************************************/

#endif
