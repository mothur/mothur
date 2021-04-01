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
#include "sequencedb.h"

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
        string fastafile, outputDir;
        vector<string> parsedFiles;
    
        vector<vector<long long> > split(SequenceDB*&);
    
        void printGroup(SequenceDB*&, vector<long long>, string);
};

/******************************************************/


#endif /* splitkmerdist_hpp */
