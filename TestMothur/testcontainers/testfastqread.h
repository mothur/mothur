//
//  testfastqread.h
//  Mothur
//
//  Created by Sarah Westcott on 3/29/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__testfastqread__
#define __Mothur__testfastqread__

#include "fastqread.h"

class TestFastqRead : public FastqRead {
    
public:
    
    using FastqRead::convertQual;
    
};


#endif /* defined(__Mothur__testfastqread__) */
