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
#include "optimatrix.h"

class OptiDB : public Database {

public:
    
    OptiDB(string, string); //reference file name for shortcut file name generation, version
    ~OptiDB() {}
    
    void addSequence(Sequence);
    void generateDB();
    void readDB(ifstream&);
    
    vector< vector<int> > get(int i, char& allSame); //A,T,G,C,- returns vector[5][numSeqsWithBase] -> vector[0] = vector of indexes of reference with A in location i, vector[1] = vector of indexes of reference with T in location i,ect. If allSame!='x', all characters are the same in this column, and will return blank vector. ie if allSame='A', every reference in this location is an A
    
    vector<int> findClosestSequences(Sequence*, int, vector<float>&) const { return nullIntVector; }
   
    
    
private:
    

    string optiDBName, version;

    int alignedLength;
    classifierOTU reference;
    map<char, int> baseMap;
    
    void convertSequences();
    
    //only used when generating db, not when reading shortcut files
    set<int> lengths;
    vector<Sequence> refs;
    
    
};


#endif /* optidb_hpp */

