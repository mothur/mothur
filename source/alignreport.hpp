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

#include "report.hpp"
#include "nast.hpp"
#include "alignment.hpp"

/******************************************************************************************************************/

class AlignReport : public Report {

public:
	
	AlignReport();
    ~AlignReport() = default;
    
    //io functions, note - printHeaders / readHeaders / getHeaders in Report parent class
    void read(ifstream&); //read line in report file
    void print(ofstream&); //print line in report file
    string getSeqReport(); //return string containing line from report file
    
    //set values from objects
	void setCandidate(Sequence*); //sets query name and length
	void setTemplate(Sequence*); //sets template name and length
	void setSearchParameters(string, float); //sets searchMethod, searchScore
	void setAlignmentParameters(string, Alignment*); //sets queryStart, queryEnd, templateStart, templateEnd, gapsInQuery, gapsInTemplate, alignmentMethod, pairwiseAlignmentLength
	void setNastParameters(Nast); //sets longestInsert and simBtwnQueryAndTemplate
    
    //set values
    void setQueryName(string n)         {    queryName = n;              }
    void setTemplateName(string n)      {    templateName = n;           }
    void setSearchMethod(string n)      {    searchMethod = n;           }
    void setAlignmentMethod(string n)   {   alignmentMethod = n;         }
    
    void setQueryLength(int n)               {    queryLength = n;            }
    void setTemplateLength(int n)            {    templateLength = n;         }
    void setQueryStart(int n)                {    queryStart = n;             }
    void setQueryEnd(int n)                  {    queryEnd = n;               }
    void setTemplateStart(int n)             {    templateStart = n;          }
    void setTemplateEnd(int n)               {    templateEnd = n;             }
    void setPairwiseAlignmentLength(int n)    {    pairwiseAlignmentLength = n;  }
    void setGapsInQuery(int n)               {    gapsInQuery = n;              }
    void setGapsInTemplate(int n)            {    gapsInTemplate = n;           }
    void setLongestInsert(int i)            { longestInsert = i;                }
    
    void setSearchScore(float i)             {    searchScore = i;               }
    void setSimBtwnQueryAndTemplate(float i)  { simBtwnQueryAndTemplate = i;    }
    
    //get values
    string getQueryName()              {    return queryName;              }
    string getTemplateName()           {    return templateName;           }
    string getSearchMethod()           {    return searchMethod;           }
    string getAlignmentMethod()        {   return alignmentMethod;         }
    
    int getQueryLength()               {    return queryLength;            }
    int getTemplateLength()            {    return templateLength;         }
    int getQueryStart()                {    return queryStart;             }
    int getQueryEnd()                  {    return queryEnd;               }
    int getTemplateStart()             {    return templateStart;          }
    int getTemplateEnd()               {    return templateEnd;             }
    int getPairwiseAlignmentLength()    {    return pairwiseAlignmentLength;  }
    int getGapsInQuery()               {    return gapsInQuery;              }
    int getGapsInTemplate()            {    return gapsInTemplate;           }
    int getLongestInsert()             {    return longestInsert;            }
    
    float getSearchScore()             {    return searchScore;               }
    float getSimBtwnQueryAndTemplate()  {    return simBtwnQueryAndTemplate;    }
    
private:
    
	void fillHeaders();
    
    string queryName, templateName, searchMethod, alignmentMethod, dummySearchScore;
    int queryLength, templateLength, queryStart, queryEnd, templateStart, templateEnd, pairwiseAlignmentLength, gapsInQuery, gapsInTemplate, longestInsert;
    float searchScore, simBtwnQueryAndTemplate;
    
};

/******************************************************************************************************************/

#endif
