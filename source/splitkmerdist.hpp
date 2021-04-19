//
//  splitkmerdist.hpp
//  Mothur
//
//  Created by Sarah Westcott on 3/30/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#ifndef splitkmerdist_hpp
#define splitkmerdist_hpp

#include "mothurout.h"
#include "utils.hpp"

/******************************************************/

class SplitKmerDistance  {
    
    public:
    
    SplitKmerDistance(string, string, double, int); //fastafile, outputDir, cutoff, kmerSize
    ~SplitKmerDistance() {}
    
    int getNumGroups() { return int(parsedFiles.size()); } //return number of split fasta files
    vector<string> getFastaFileNames() { return parsedFiles; } //return names fasta files for each group
    
    
    private:
        MothurOut* m;
        Utils util;
    
        double cutoff;
        int kmerSize;
        long long numSeqs;
        string fastafile, outputDir;
        vector<string> parsedFiles;
    
        vector< vector< bool > > kmerDB; //kmerDB[0] = vector<bool> maxKmers long
        vector<int> lengths;
    
    
        vector<int> getUniqueKmers(int i);
        vector<vector<long long> > split();
    
        int findRootGroup(vector<int>&, int);
        //void printGroup(SequenceDB*&, vector<long long>, string);
        int fillKmerDB();
        
};

/******************************************************/


#endif /* splitkmerdist_hpp */
