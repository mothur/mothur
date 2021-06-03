//
//  proteindb.hpp
//  Mothur
//
//  Created by Sarah Westcott on 6/3/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#ifndef proteindb_hpp
#define proteindb_hpp

#include "protein.hpp"

class ProteinDB {
    
public:
    
    ProteinDB();
    ProteinDB(int);           //makes data that size
    ProteinDB(ifstream&);       //reads file to fill data
    ProteinDB(const ProteinDB& sdb) : data(sdb.data) {};
    ~ProteinDB();             //loops through data and delete each protein sequence

    int getNumSeqs();
    Protein get(int);         //returns sequence name at that location
    void push_back(Protein);  //adds unaligned sequence
    bool sameLength() { return samelength; }
        
private:
    
    vector<Protein> data;
    MothurOut* m;
    bool samelength;
    int length;
    Utils util;

};

#endif /* proteindb_hpp */
