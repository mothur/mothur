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
#include "picrust.hpp"

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
    virtual void load(SharedRAbundVectors* s, vector<Taxonomy> c);
    virtual void load(SharedRAbundFloatVectors* s, vector<Taxonomy> c);
    virtual void print(ofstream&, vector<string>, Picrust*) = 0;
    
    virtual void printHeading(ofstream&, string, string) {}
    
    virtual string getVersion() { return version; }
    virtual string getMatrixElementType() { return matrixElementType; }
    
    virtual SharedRAbundVectors* getSharedRAbundVectors() { return shared; }
    virtual SharedRAbundFloatVectors* getSharedRAbundFloatVectors() { return sharedFloat; }
    
    //otu taxonomies
    virtual vector<Taxonomy> getConsTaxonomies() { return consTax; }
    
    //sample taxonomies
    virtual map<string, string> getGroupTaxonomies() {  return groupTaxonomies; }
    
protected:
    
    MothurOut* m;
    Utils util;
    string version, formatURL, label, matrixElementType; //version = simple or hdf5, set by child. matrixElementType = "int" or "float"
    int maxLevel;
   
    SharedRAbundVectors* shared; //always created with read
    SharedRAbundFloatVectors* sharedFloat; //only created if the matrixElementType is float
    vector<Taxonomy> consTax;
    map<string, string> groupTaxonomies;
    
};
/**************************************************************************************************/
#endif /* biom_hpp */
