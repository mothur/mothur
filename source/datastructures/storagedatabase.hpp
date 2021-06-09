//
//  storagedatabase.hpp
//  Mothur
//
//  Created by Sarah Westcott on 6/3/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#ifndef storagedatabase_hpp
#define storagedatabase_hpp

#include "mothurout.h"
#include "sequence.hpp"
#include "protein.hpp"

class StorageDatabase {
    
public:
    
    StorageDatabase() {  m = MothurOut::getInstance();  length = 0; samelength = true; }
    virtual ~StorageDatabase() {}             //loops through data and delete each sequence

    virtual int getNumSeqs() = 0;
    virtual bool sameLength() { return samelength; }
       
    virtual Sequence getSeq(int) { Sequence s; return s; }
    virtual Protein getProt(int) { Protein p; return p;  }
    virtual void push_back(Sequence) {}  //adds sequence
    virtual void push_back(Protein) {}  //adds protein
    
        
protected:
   
    MothurOut* m;
    Utils util;
    bool samelength;
    int length;

};

#endif /* storagedatabase_hpp */
