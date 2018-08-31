//
//  fakeoligos.h
//  Mothur
//
//  Created by Sarah Westcott on 5/1/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef fakeoligos_h
#define fakeoligos_h

#include "mothurout.h"


class FakeOligos {
    
public:
    
    FakeOligos() {
        bdiffs=1; pdiffs=2; rdiffs=0; ldiffs=0; sdiffs=0;
        
        //single
        primers["CCGTCAATTCMTTTRAGT"] = 0;
        barcodes["AATGGTAC"] = 0; //F003D000
        barcodes["AACCTGGC"] = 1; //F003D002
        barcodes["TTCGTGGC"] = 2; //F003D004
        barcodes["TTCTTGAC"] = 3; //F003D006
        barcodes["TTCGCGAC"] = 4; //F003D008
        barcodes["TCCAGAAC"] = 5; //F003D142
        barcodes["AAGGCCTC"] = 6; //F003D144
        barcodes["TGACCGTC"] = 7; //F003D146
        barcodes["AGGTTGTC"] = 8; //F003D148
        barcodes["TGGTGAAC"] = 9; //F003D150
        barcodes["AACCGTGTC"] = 10; //MOCK.GQY1XT001
        
        
    }
    ~FakeOligos() {}
    
    int bdiffs, pdiffs, rdiffs, ldiffs, sdiffs;
    map<string, int> barcodes;
    map<string, int> primers;
    
    vector<string> revPrimer;
    vector<string> linker;
    vector<string> spacer;
    map<int, oligosPair> ipbarcodes;
    map<int, oligosPair> ipprimers;

};



#endif /* fakeoligos_h */
