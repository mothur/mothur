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

//subsampling overwrites the sharedRabunds.  If you need to reuse the original use the getSamplePreserve function.

class SubSample {
	
    public:
    
        SubSample() { m = MothurOut::getInstance(); }
        ~SubSample() {}
    
        vector<string> getSample(vector<SharedRAbundVector*>&, int); //returns the bin labels for the subsample, mothurOuts binlabels are preserved so you can run this multiple times. Overwrites original vector passed in, if you need to preserve it deep copy first.
        
        //Tree* getSample(Tree*, TreeMap*, map<string, string>, int); //creates new subsampled tree, destroys treemap so copy if needed.
        Tree* getSample(Tree*, TreeMap*, TreeMap*, int, map<string, string>); //creates new subsampled tree. Uses first treemap to fill new treemap with sabsampled seqs. Sets groups of seqs not in subsample to "doNotIncludeMe".
        int getSample(SAbundVector*&, int); //destroys sabundvector passed in, so copy it if you need it
    
    private:
    
        MothurOut* m;
        int eliminateZeroOTUS(vector<SharedRAbundVector*>&);
    
        vector<string> getSample(TreeMap*, vector<string>);
        vector<string> getSample(TreeMap*, int); //names of seqs to include in sample tree 
        vector<string> getSample(TreeMap* tMap, int size, map<string, vector<string> >& sample); //sample maps group -> seqs in group. seqs not in sample are in doNotIncludeMe group
        map<string, string> deconvolute(map<string, string> wholeSet, vector<string>& subsampleWanted); //returns new nameMap containing only subsampled names, and removes redundants from subsampled wanted because it makes the new nameMap.


};

#endif
