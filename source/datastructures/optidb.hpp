//
//  optidb.hpp
//  Mothur
//
//  Created by Sarah Westcott on 3/26/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#ifndef optidb_hpp
#define optidb_hpp

#include "sequence.hpp"
#include "database.hpp"
#include "calculator.h"

class OptiDB : public Database {

public:
    
    OptiDB(double);
    ~OptiDB() {}
    
    void addSequence(Sequence); //add otu with single seq
    void addSequences(vector<Sequence>); //add otu with multiple seqs
    
    void generateDB();
    
private:
    
    bool aligned;
    int alignedLength;
    double cutoff;
    
    vector<classifierOTU> reference;
    
    DistCalc* calc;
    
};


#endif /* optidb_hpp */
