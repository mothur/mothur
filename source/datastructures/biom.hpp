//
//  biom.hpp
//  Mothur
//
//  Created by Sarah Westcott on 10/26/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#ifndef biom_hpp
#define biom_hpp

#include "utils.hpp"
#include "mothurout.h"
#include "sharedrabundfloatvectors.hpp"
#include "sharedrabundvectors.hpp"
#include "phylosummary.h"
#include "taxonomy.hpp"

//http://biom-format.org
//http://biom-format.org/documentation/format_versions/biom-1.0.html
//http://biom-format.org/documentation/format_versions/biom-2.1.html

/**************************************************************************************************/
class Biom {

public:
    Biom();
    Biom(string); //version
    
    virtual ~Biom();
    
    virtual void read(string) = 0;
    virtual string getVersion() { return version; }
    virtual string getMatrixElementType() { return "int"; }
    
    virtual SharedRAbundVectors* getSharedRAbundVectors() { return shared; }
    
    //otu taxonomies
    virtual vector<Taxonomy> getConsTaxonomies() { return consTax; }
    
    //sample taxonomies
    virtual map<string, string> getGroupTaxonomies() {  return groupTaxonomies; }
    
protected:
    
    MothurOut* m;
    Utils util;
    string version, formatURL, label; //version = simple or hdf5, set by child
    int maxLevel;
   
    SharedRAbundVectors* shared;
    vector<Taxonomy> consTax;
    map<string, string> groupTaxonomies;
    
};
/**************************************************************************************************/
#endif /* biom_hpp */
