//
//  contigsreport.hpp
//  Mothur
//
//  Created by Sarah Westcott on 7/17/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#ifndef contigsreport_hpp
#define contigsreport_hpp

#include "report.hpp"

/*
 This class stores information for one line in the contigs report file
 
 
 */
/******************************************************************************************************************/

class ContigsReport : public Report {

public:
    
    ContigsReport();
    ~ContigsReport() {}
    
    //io functions, note - printHeaders / readHeaders / getHeaders in Report parent class
    void read(ifstream&); //read line in report file
    void print(ofstream&); //print line in report file
    string getSeqReport(); //return string containing line from report file
    
    //set values
    void setName(string n)          {    name = n;              }
    void setLength(int n)           {    length = n;            }
    void setOverlapLength(int n)    {    overlapLength = n;     }
    void setOverlapStart(int n)     {    overlapStart = n;      }
    void setOverlapEnd(int n)       {    overlapEnd = n;        }
    void setMisMatches(int n)       {    misMatches = n;        }
    void setNumNs(int n)            {    numsNs = n;           }
    void setExpectedErrors(float i)  {    expectedErrors = i;   }
    
    //get values
    string getName()                {    return name;           }
    int getLength()                 {    return length;         }
    int getOverlapLength()          {    return overlapLength;  }
    int getOverlapStart()           {    return overlapStart;   }
    int getOverlapEnd()             {    return overlapEnd;     }
    int getMisMatches()             {    return misMatches;     }
    int getNumNs()                  {    return numsNs;         }
    float getExpectedErrors()       {    return expectedErrors; }
    
private:
    
    void fillHeaders();
   
    string name;
    int length, overlapLength, overlapStart, overlapEnd, misMatches, numsNs;
    float expectedErrors;
    
};

/******************************************************************************************************************/

#endif /* contigsreport_hpp */
