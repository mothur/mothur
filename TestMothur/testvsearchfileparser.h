//
//  testvsearchfileparser.h
//  Mothur
//
//  Created by Sarah Westcott on 3/24/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__testvsearchfileparser__
#define __Mothur__testvsearchfileparser__


#include "vsearchfileparser.h"
#include "dataset.h"

class TestVsearchFileParser : public VsearchFileParser {
    
    
public:
    
    TestVsearchFileParser();
    ~TestVsearchFileParser();
    
    MothurOut* m;
    TestDataSet data;
    vector<string> filenames;
    CountTable* ct;
    
    using VsearchFileParser::removeAbundances;
    using VsearchFileParser::createListFile;
    
};


#endif /* defined(__Mothur__testvsearchfileparser__) */
