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
#include "treemap.h"
#include "tree.h"
#include "counttable.h"


//subsampling overwrites the sharedRabunds.  If you need to reuse the original use the getSamplePreserve function.

class SubSample {
	
    public:
    
        SubSample() { m = MothurOut::getInstance(); }
        ~SubSample() {}
    
        vector<string> getSample(vector<SharedRAbundVector*>&, int); //returns the bin labels for the subsample, mothurOuts binlabels are preserved so you can run this multiple times. Overwrites original vector passed in, if you need to preserve it deep copy first.
        Tree* getSample(Tree*, CountTable*, CountTable*, int); //creates new subsampled tree. Uses first counttable to fill new counttable with sabsampled seqs. Sets groups of seqs not in subsample to "doNotIncludeMe".
        int getSample(SAbundVector*&, int); //destroys sabundvector passed in, so copy it if you need it
        CountTable getSample(CountTable&, int, vector<string>); //subsample a countTable bygroup(same number sampled from each group, returns subsampled countTable 
        CountTable getSample(CountTable&, int, vector<string>, bool); //subsample a countTable. If you want to only sample from specific groups, pass in groups in the vector and set bool=true, otherwise set bool=false.   
    
    private:
    
        MothurOut* m;
        int eliminateZeroOTUS(vector<SharedRAbundVector*>&);
         map<string, string> deconvolute(map<string, string> wholeSet, vector<string>& subsampleWanted); //returns new nameMap containing only subsampled names, and removes redundants from subsampled wanted because it makes the new nameMap.


};

#endif
