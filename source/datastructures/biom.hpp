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

//http://biom-format.org
//http://biom-format.org/documentation/format_versions/biom-1.0.html
//http://biom-format.org/documentation/format_versions/biom-2.1.html

/**************************************************************************************************/
class Biom {

public:
    Biom(){ version = ""; formatURL = "http://biom-format.org"; m = MothurOut::getInstance(); shared = NULL; taxSum = NULL; }
    
    virtual ~Biom(){};
    virtual void read(string) = 0;
    virtual SharedRAbundVectors* getSharedRAbundVectors() { return shared; }
    virtual SharedRAbundFloatVectors* getSharedRAbundFloatVectors();
    virtual PhyloSummary* getTaxSummary() { return taxSum; }
    virtual vector<Taxonomy> getConsTaxonomies() { return consTax; }
    virtual map<string, string> getTaxonomies() {  return groupTaxonomies; }
    //virtual void addSequence(Sequence) = 0;  //add sequence to search engine
   // virtual void addSequences(vector<Sequence> seqs) { for (int i = 0; i < seqs.size(); i++) { addSequence(seqs[i]); } }

   // virtual void setNumSeqs(int i) {    numSeqs = i;     }
    
   // virtual vector<int> findClosestSequences(Sequence*, int, vector<float>&) const = 0;  // returns indexes of n closest sequences to query
    
    virtual string getVersion() { return version; }
   // virtual vector<int> getIndicatorColumns() { return nullIntVector; }
    //virtual map<int, int> getFilteredIndicatorColumns(string, vector<int>&) { return nullIntMap; }
    
protected:
    
    MothurOut* m;
    Utils util;
    string version, formatURL; //version = simple or hdf5, set by child
    int printLevel;
    bool relabund;
    
    SharedRAbundVectors* shared;
    PhyloSummary* taxSum;
    vector<Taxonomy> consTax;
    map<string, string> groupTaxonomies;
    
};
/**************************************************************************************************/
#endif /* biom_hpp */
