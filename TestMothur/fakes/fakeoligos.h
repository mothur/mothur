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
    
    FakeOligos() { bdiffs=1; pdiffs=2; rdiffs=0; ldiffs=0; sdiffs=0;  }
    ~FakeOligos() {}
    
    //single
    void loadSingle() {
        primers.clear();
        barcodes.clear();
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
        primerNameVector.clear();
        barcodeNameVector.clear();
        barcodeNameVector.push_back("F003D000"); barcodeNameVector.push_back("F003D002");
        barcodeNameVector.push_back("F003D004"); barcodeNameVector.push_back("F003D006");
        barcodeNameVector.push_back("F003D008"); barcodeNameVector.push_back("F003D142");
        barcodeNameVector.push_back("F003D144"); barcodeNameVector.push_back("F003D146");
        barcodeNameVector.push_back("F003D148"); barcodeNameVector.push_back("F003D150"); barcodeNameVector.push_back("MOCK.GQY1XT001");
    }
    
    void loadPaired() {
        primerNameVector.clear();
        barcodeNameVector.clear();
        oligosPair F01R2A; F01R2A.forward = "ccaac";	F01R2A.reverse = "cactg";	barcodeNameVector.push_back("F01R2A"); ipbarcodes[0] = F01R2A;
        oligosPair F01R2B; F01R2B.forward = "ccaac";	F01R2B.reverse = "aacca";	barcodeNameVector.push_back("F01R2B"); ipbarcodes[1] = F01R2B;
        oligosPair F01R2C; F01R2C.forward = "ccaac";	F01R2C.reverse = "tgtca";	barcodeNameVector.push_back("F01R2C"); ipbarcodes[2] = F01R2C;
        oligosPair F01R2D; F01R2D.forward = "ccaac";	F01R2D.reverse = "aaacc";	barcodeNameVector.push_back("F01R2D"); ipbarcodes[3] = F01R2D;
        oligosPair F02R2A; F02R2A.forward = "ggttg";	F02R2A.reverse = "cactg";	barcodeNameVector.push_back("F02R2A"); ipbarcodes[4] = F02R2A;
        oligosPair F02R2B; F02R2B.forward = "ggttg";	F02R2B.reverse = "aacca";	barcodeNameVector.push_back("F02R2B"); ipbarcodes[5] = F02R2B;
        oligosPair F02R2C; F02R2C.forward = "ggttg";	F02R2C.reverse = "tgtca";	barcodeNameVector.push_back("F02R2C"); ipbarcodes[6] = F02R2C;
        oligosPair F02R2D; F02R2D.forward = "ggttg";	F02R2D.reverse = "aaacc";	barcodeNameVector.push_back("F02R2D"); ipbarcodes[7] = F02R2D;
    }

    int bdiffs, pdiffs, rdiffs, ldiffs, sdiffs;
    map<string, int> barcodes;
    map<string, int> primers;
    
    vector<string> revPrimer;
    vector<string> linker;
    vector<string> spacer;
    map<int, oligosPair> ipbarcodes;
    map<int, oligosPair> ipprimers;
    vector<string> primerNameVector;
    vector<string> barcodeNameVector;

};



#endif /* fakeoligos_h */
