//
//  proteindb.cpp
//  Mothur
//
//  Created by Sarah Westcott on 6/3/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#include "proteindb.hpp"

/***********************************************************************/

ProteinDB::ProteinDB() {  m = MothurOut::getInstance();  length = 0; samelength = true; }
/***********************************************************************/
//the clear function free's the memory
ProteinDB::~ProteinDB() { data.clear(); }

/***********************************************************************/

ProteinDB::ProteinDB(int newSize) {
    m = MothurOut::getInstance();  length = 0; samelength = true;
    data.resize(newSize, Protein());
}

/***********************************************************************/

ProteinDB::ProteinDB(ifstream& filehandle) {
    try{
        m = MothurOut::getInstance(); length = 0; samelength = true;
        
        //read through file
        while (!filehandle.eof()) {
            
            Protein newProteinSequence(filehandle);  util.gobble(filehandle);
            
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
Protein ProteinDB::get(int index) {
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

