//
//  proteindb.cpp
//  Mothur
//
//  Created by Sarah Westcott on 6/3/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#include "proteindb.hpp"

/***********************************************************************/

ProteinDB::ProteinDB() : StorageDatabase() { }
/***********************************************************************/
//the clear function free's the memory
ProteinDB::~ProteinDB() { data.clear(); }

/***********************************************************************/

ProteinDB::ProteinDB(int newSize) : StorageDatabase() {  data.resize(newSize, Protein()); }

/***********************************************************************/

ProteinDB::ProteinDB(ifstream& filehandle) : StorageDatabase() {
    try{
        
        //read through file
        while (!filehandle.eof()) {
            
            Protein newProteinSequence(filehandle);  gobble(filehandle);
            
            if (newProteinSequence.getName() != "") {
                if (length == 0) { length = newProteinSequence.getAligned().size(); }
                if (length != newProteinSequence.getAligned().size()) { samelength = false;  }
                data.push_back(newProteinSequence);
            }
        }
        
        filehandle.close();
    }
    catch(exception& e) {
        m->errorOut(e, "ProteinDB", "ProteinDB");
        exit(1);
    }
}
/***********************************************************************/

int ProteinDB::getNumSeqs() { return data.size(); }

/***********************************************************************/
Protein ProteinDB::getProt(int index) {
    if ((index >= 0) && (index < data.size()) ) { return data[index]; }
    else { m->mothurOut("[ERROR]: invalid database index, please correct.\n"); m->setControl_pressed(true); Protein p; return p; }
}
/***********************************************************************/

void ProteinDB::push_back(Protein newProteinSequence) {
    try {
        if (length == 0) { length = newProteinSequence.getAligned().size(); }
        if (length != newProteinSequence.getAligned().size()) { samelength = false; }

        data.push_back(newProteinSequence);
    }
    catch(exception& e) {
        m->errorOut(e, "ProteinDB", "push_back");
        exit(1);
    }
}
/***********************************************************************/

void ProteinDB::print(string outputFileName) {
    try {
        ofstream out; util.openOutputFile(outputFileName, out);
        
        for (int i = 0; i < data.size(); i++) {
            data[i].printProtein(out);
        }
        out.close();
    }
    catch(exception& e) {
        m->errorOut(e, "ProteinDB", "print");
        exit(1);
    }
}
/***********************************************************************/

