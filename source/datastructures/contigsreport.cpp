//
//  contigsreport.cpp
//  Mothur
//
//  Created by Sarah Westcott on 7/17/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "contigsreport.hpp"


//Name    Length    Overlap_Length    Overlap_Start    Overlap_End    MisMatches    Num_Ns    Expected_Errors

/******************************************************************************************************************/

ContigsReport::ContigsReport() : Report() {
    try {
        fillHeaders();
    }
    catch(exception& e) {
        m->errorOut(e, "ContigsReport", "ContigsReport");
        exit(1);
    }
}
/**************************************************************************************************/
void ContigsReport::read(ifstream& repFile){
    try {
        
        repFile >> name;
        repFile >> length;
        repFile >> overlapLength;
        repFile >> overlapStart;
        repFile >> overlapEnd;
        repFile >> misMatches;
        repFile >> numsNs;
        repFile >> expectedErrors;
        
        util.gobble(repFile);
    }
    catch(exception& e) {
        m->errorOut(e, "ContigsReport", "read");
        exit(1);
    }
    
}
/******************************************************************************************************************/
void ContigsReport::print(ofstream& reportFile){
    try {
        reportFile << name << '\t' << length << '\t' << overlapLength << '\t';
        reportFile << overlapStart << '\t' << overlapEnd << '\t';

        reportFile << misMatches << '\t' << numsNs << '\t';
        reportFile << setprecision(6) << expectedErrors << endl;
    }
    catch(exception& e) {
        m->errorOut(e, "ContigsReport", "print");
        exit(1);
    }
}
/******************************************************************************************************************/

string ContigsReport::getSeqReport(){
    try {
        string output = "";
        
        output += name + '\t' + toString(length) + '\t' + toString(overlapLength) + '\t';
        output +=  toString(overlapStart) + '\t' + toString(overlapEnd) + '\t';
        output +=  toString(misMatches) + '\t' + toString(numsNs) + '\t';
        
        string temp = toString(expectedErrors);
        int pos = temp.find_last_of('.');  //find deicmal point if their is one
        
        //if there is a decimal
        if (pos != -1) { temp = temp.substr(0, pos+6); } //set precision to 5 places
        else{    temp += ".00000";    }
        
        output +=  temp + '\n';
           
        return output;
    }
    catch(exception& e) {
        m->errorOut(e, "ContigsReport", "getSeqReport");
        exit(1);
    }
}
/******************************************************************************************************************/
void ContigsReport::fillHeaders() {
    try {
        reportHeaders.push_back("Name"); reportHeaders.push_back("Length");
        reportHeaders.push_back("Overlap_Length");
        reportHeaders.push_back("Overlap_Start"); reportHeaders.push_back("Overlap_End");
        reportHeaders.push_back("MisMatches");
        
        reportHeaders.push_back("Num_Ns"); reportHeaders.push_back("Expected_Errors");
    }
    catch(exception& e) {
        m->errorOut(e, "ContigsReport", "fillHeaders");
        exit(1);
    }
}
/******************************************************************************************************************/
