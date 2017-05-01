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


class TestTrimOligos : public TrimOligos {
    
    
public:
    
    TestTrimOligos();
    ~TestTrimOligos();
    
    MothurOut* m;
    
    Sequence* fseq;
    //using TrimOligos::compareDNASeq(string, string);
    //using TrimOligos::countDiffs(string, string);
    
};


#endif /* testtrimoligos_hpp */
