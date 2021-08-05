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
#include "storagedatabase.hpp"

class ProteinDB : public StorageDatabase {
    
public:
    
    ProteinDB();
    ProteinDB(int);           //makes data that size
    ProteinDB(ifstream&);       //reads file to fill data
    ProteinDB(const ProteinDB& sdb) : data(sdb.data) {};
    ~ProteinDB();             //loops through data and delete each protein sequence

    Protein getProt(int);         //returns sequence name at that location
    void push_back(Protein);  //adds unaligned sequence
    void print(string);  //prints fasta file containing sequences in this db
    int getNumSeqs();
    
private:
    
    vector<Protein> data;
    

};

#endif /* proteindb_hpp */
