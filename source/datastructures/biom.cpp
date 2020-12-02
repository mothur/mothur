//
//  biom.cpp
//  Mothur
//
//  Created by Sarah Westcott on 10/26/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "biom.hpp"

/**************************************************************************************************/
Biom::Biom() {
    try {
        m = MothurOut::getInstance();
        
        formatURL = "http://biom-format.org";
        label = ""; version = "";
        
        tableID = "No Table ID";
        mothurVersion = ""; sharedFileName = "";
       
        shared = NULL;  sharedFloat = NULL;
    }
    catch(exception& e) {
        m->errorOut(e, "Biom", "Biom");
        exit(1);
    }
}
/**************************************************************************************************/
Biom::Biom(string v) : version(v) {
    try {
        m = MothurOut::getInstance();
        
        formatURL = "http://biom-format.org";
        label = "";
        
        tableID = "No Table ID";
        mothurVersion = ""; sharedFileName = "";
        
        shared = NULL; sharedFloat = NULL;
    }
    catch(exception& e) {
        m->errorOut(e, "Biom", "Biom");
        exit(1);
    }
}
/**************************************************************************************************/
Biom::~Biom() {
    if (shared != NULL) { delete shared; }
    if (sharedFloat != NULL) { delete sharedFloat; }
}
/**************************************************************************************************/
void Biom::load(SharedRAbundVectors* s, vector<Taxonomy> c){
    try {
        shared = new SharedRAbundVectors(*s);
        consTax = c;
        matrixElementType = "int";
        label = s->getLabel();
    }
    catch(exception& e) {
        m->errorOut(e, "Biom", "load-shared");
        exit(1);
    }
}
/**************************************************************************************************/
void Biom::load(SharedRAbundFloatVectors* s, vector<Taxonomy> c){
    try {
        sharedFloat = new SharedRAbundFloatVectors(*s);
        consTax = c;
        matrixElementType = "float";
        label = s->getLabel();
        
        vector<SharedRAbundVector*> sharedRabunds = s->getSharedRAbundVectors();
        shared = new SharedRAbundVectors();
        
        for (int i = 0; i < sharedRabunds.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            shared->push_back(sharedRabunds[i]);
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "Biom", "load-float");
        exit(1);
    }
}
/**************************************************************************************************/
