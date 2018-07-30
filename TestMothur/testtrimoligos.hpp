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
    
    FakeOligos oligos;
    
    
};


#endif /* testtrimoligos_hpp */
