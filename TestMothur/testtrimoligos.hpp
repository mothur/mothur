//
//  testtrimoligos.hpp
//  Mothur
//
//  Created by Sarah Westcott on 7/14/16.
//  Copyright © 2016 Schloss Lab. All rights reserved.
//

#ifndef testtrimoligos_hpp
#define testtrimoligos_hpp

#include "trimoligos.h"
#include "sequence.hpp"

class TestTrimOligos : public TrimOligos {
    
    
public:
    
    TestTrimOligos();
    ~TestTrimOligos();
    
    MothurOut* m;
    vector<Sequence> seqs;
    vector<Sequence> pairedSeqs;
    
    map<string, int> barcodes;
    map<string, int> primers;
    map<int, oligosPair> pairedPrimers;
    map<int, oligosPair> pairedBarcodes;
    
    //using TrimOligos::compareDNASeq(string, string);
    //using TrimOligos::countDiffs(string, string);
    
};


#endif /* testtrimoligos_hpp */
