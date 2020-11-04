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
       
        shared = NULL; 
    }
    catch(exception& e) {
        m->errorOut(e, "Biom", "Biom");
        exit(1);
    }
}
/**************************************************************************************************/
Biom::Biom(string v) : version(v){
    try {
        m = MothurOut::getInstance();
        
        formatURL = "http://biom-format.org";
        label = "";
        
        shared = NULL;
    }
    catch(exception& e) {
        m->errorOut(e, "Biom", "Biom");
        exit(1);
    }
}
/**************************************************************************************************/
Biom::~Biom() {
    if (shared != NULL) { delete shared; }
    
}
/**************************************************************************************************/
