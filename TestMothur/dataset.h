//
//  dataset.h
//  Mothur
//
//  Created by Sarah Westcott on 3/24/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__dataset__
#define __Mothur__dataset__

#include "sequence.hpp"
#include "counttable.h"
#include "groupmap.h"
#include "sharedrabundvector.h"
#include "fastqdataset.h"

class TestDataSet  {
    
public:
    
    TestDataSet();
    vector<Sequence> getSeqs()                  { fillSeqs(); return seqs;              }
    map<string, string> getNameMap()            { fillNames(); return nameMap;          }
    GroupMap* getGroupMap()                     { fillGroup(); return gMap;             }
    CountTable* getCountTable()                 { createCountTable(); return ct;        }
    vector<SharedRAbundVector*> getLookup()     { fillLookup(); return lookup;          }
    vector<string> getSubsetFNGFiles();  //Fasta, name, group returned - containing 100 seqs
    string getSubsetFNGDistFile();
    
private:
    MothurOut* m;
    TestFastqDataSet fastqData;
    vector<Sequence> seqs;
    map<string, string> nameMap;
    CountTable* ct;
    GroupMap* gMap;
    vector<SharedRAbundVector*> lookup;
    void fillNames();
    void fillSeqs();
    void fillGroup();
    void createCountTable();
    void fillLookup();
    
};


#endif /* defined(__Mothur__dataset__) */
