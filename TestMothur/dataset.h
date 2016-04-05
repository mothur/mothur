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
    vector<FastqRead> getForwardFastq()         { return fastqData.getForwardFastq();   }
    vector<FastqRead> getReverseFastq()         { return fastqData.getReverseFastq();   }
    vector<string> getSubsetFRFastq(int n)      { return fastqData.getSubsetFRFastq(n); }
    
    vector<string> getSubsetFNGFiles(int);  //number of uniques, Fasta, name, group returned
    
private:
    MothurOut* m;
    TestFastqDataSet fastqData;
    vector<Sequence> seqs;
    map<string, string> nameMap;
    CountTable* ct;
    GroupMap* gMap;
    vector<SharedRAbundVector*> lookup;
    vector<FastqRead> ffastqReads; //F8D0 Sample
    vector<FastqRead> rfastqReads; //F8D0 Sample
    void fillNames();
    void fillSeqs();
    void fillGroup();
    void createCountTable();
    void fillLookup();
    void fillForwardFastq();
    void fillReverseFastq();
    
};


#endif /* defined(__Mothur__dataset__) */
