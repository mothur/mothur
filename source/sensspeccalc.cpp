//
//  sensspeccalc.cpp
//  Mothur
//
//  Created by Sarah Westcott on 1/22/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "sensspeccalc.hpp"


//***************************************************************************************************************
//removes anyone with no valid dists and changes name to matrix short names
SensSpecCalc::SensSpecCalc(OptiData& matrix, ListVector* list){
    try {
        m = MothurOut::getInstance();
        map<string, int> nameIndex = matrix.getNameIndexMap();
        
        if (list != NULL) {
            //for each bin
            for (int i = 0; i < list->getNumBins(); i++) {
                
                string binnames = list->get(i);
                vector<string> bnames;
                util.splitAtComma(binnames, bnames);
                
                vector<int> newNames;
                for (int j = 0; j < bnames.size(); j++) {
                    string name = bnames[j];
                    map<string, int>::iterator itSeq1 = nameIndex.find(name);
                    int seq1Index = -1;
                    if (itSeq1 != nameIndex.end()) { seq1Index = itSeq1->second; } //you have distances in the matrix
                    
                    newNames.push_back(seq1Index); 
                }
                
                //if there are names in this bin add to new list
                if (newNames.size() != 0) { otus.push_back(newNames); }
            }
        }
    }
    catch(exception& e) {
        m->errorOut(e, "SensSpecCalc", "SensSpecCalc");
        exit(1);
    }
}
//***************************************************************************************************************
void SensSpecCalc::getResults(OptiData& matrix, long long& tp, long long& tn, long long& fp, long long& fn){
    try {
        tp = 0; tn = 0; fp = 0; fn = 0;
        
        for(int otu=0;otu<otus.size();otu++){
            if (m->getControl_pressed()) { break; }
            
            for(int i=0;i<otus[otu].size();i++){
                for(int j=0;j<i;j++){
                    if (matrix.isClose(otus[otu][i], otus[otu][j])) { tp++; }
                    else { fp++; }
                }
            }
        }
        long long numSeqs = matrix.getNumSeqs() + matrix.getNumSingletons();
        long long numDists = matrix.getNumDists(); //square matrix
        
        fn = (numDists/2) - tp;
        tn = numSeqs * (numSeqs-1)/2  - (fp + fn + tp);
        
        //cout << numSeqs << '\t' << numDists << '\t' << tp << '\t' << tn << '\t' << fp << '\t' << fn << endl;
    }
    catch(exception& e) {
        m->errorOut(e, "SensSpecCalc", "getResults");
        exit(1);
    }
}

//***************************************************************************************************************

