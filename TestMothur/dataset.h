//
//  dataset.h
//  Mothur
//
//  Created by Sarah Westcott on 3/24/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__dataset__
#define __Mothur__dataset__

#include <vector>
#include "sequence.hpp"
#include "counttable.h"
#include "groupmap.h"

class TestDataSet  {
    
public:
    
    TestDataSet();
    vector<Sequence> getSeqs() { fillSeqs(); return seqs; }
    map<string, string> getNameMap() { fillNames(); return nameMap; }
    GroupMap* getGroupMap() {  fillGroup(); return gMap; }
    CountTable* getCountTable() { createCountable(); return ct; }
    
private:
    MothurOut* m;
    vector<Sequence> seqs;
    map<string, string> nameMap;
    CountTable ct;
    GroupMap* gMap;
    vector<SharedRAbundVector*> lookup;
    void fillNames();
    void fillSeqs();
    void fillGroup();
    void createCountTable();
    
};


#endif /* defined(__Mothur__dataset__) */
