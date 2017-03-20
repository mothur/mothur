//
//  fastqdataset.h
//  Mothur
//
//  Created by Sarah Westcott on 3/31/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__fastqdataset__
#define __Mothur__fastqdataset__

#include "fastqread.h"

class TestFastqDataSet  {
    
public:
    
    TestFastqDataSet();
    ~TestFastqDataSet() {}
    vector<FastqRead> getForwardFastq() { fillForwardFastq(); return ffastqReads; }
    vector<FastqRead> getReverseFastq() { fillReverseFastq(); return rfastqReads; }
    vector<string> getSubsetFRFastq(int);
    
private:
    MothurOut* m;
    vector<FastqRead> ffastqReads;
    vector<FastqRead> rfastqReads; 
    void fillForwardFastq();
    void fillReverseFastq();
    
};


#endif /* defined(__Mothur__fastqdataset__) */
