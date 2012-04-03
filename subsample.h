#ifndef Mothur_subsample_h
#define Mothur_subsample_h

//
//  subsample.h
//  Mothur
//
//  Created by Sarah Westcott on 4/2/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "mothurout.h"
#include "sharedrabundvector.h"

//subsampling overwrites the sharedRabunds.  If you need to reuse the original use the getSamplePreserve function.

class SubSample {
	
    public:
    
        SubSample() { m = MothurOut::getInstance(); }
        ~SubSample() {}
    
        vector<string> getSample(vector<SharedRAbundVector*>&, int); //returns the bin labels for the subsample, mothurOuts binlabels are preserved so you can run this multiple times. Overwrites original vector passed in, if you need to preserve it deep copy first.
    
    
    private:
    
        MothurOut* m;
        int eliminateZeroOTUS(vector<SharedRAbundVector*>&);

};

#endif
