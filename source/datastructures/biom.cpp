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
        label = "";
        
        version = ""; basis = "otu"; printLevel = 0; relabund = false;
        
        shared = NULL; taxSum = NULL; consTaxSum = NULL;
        
    }
    catch(exception& e) {
        m->errorOut(e, "Biom", "Biom");
        exit(1);
    }
}
/**************************************************************************************************/
Biom::Biom(string v, string b, int pl, bool rel) : version(v), basis(b), printLevel(pl), relabund(rel){
    try {
        m = MothurOut::getInstance();
        
        formatURL = "http://biom-format.org";
        label = "";
        
        shared = NULL; taxSum = NULL; consTaxSum = NULL;
        
    }
    catch(exception& e) {
        m->errorOut(e, "Biom", "Biom");
        exit(1);
    }
}
/**************************************************************************************************/
Biom::~Biom() {
    if (shared != NULL) { delete shared; }
    if (taxSum != NULL) { delete taxSum; }
    if (consTaxSum != NULL) { delete consTaxSum; }
}
/**************************************************************************************************/
