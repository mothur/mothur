//
//  opticlassifier.hpp
//  Mothur
//
//  Created by Sarah Westcott on 3/12/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#ifndef opticlassifier_hpp
#define opticlassifier_hpp

#include "mothur.h"
#include "classify.h"

/**************************************************************************************************/

class OptiClassifier : public Classify {
    
public:
    OptiClassifier(string reffasta, string reftax, string mothurVersion);
    ~OptiClassifier() {}
    
    string getTaxonomy(Sequence*, string&, bool&) { return "not done yet"; }
    
private:
    
    vector< vector< vector<float> > > charGenusProb;    //[locationInAlignment][base][genus]
                                        //charGenusProb[0][3][392] = probability that a sequence within genus that's index in the tree is 392 would contain 'C' at alignment position 0; charGenusProb[23][0][392] = probability that a sequence within genus that's index in the tree is 392 would contain 'A' at alignment position 23;
    
    vector<int> genusTotals; //number of sequence at each genus
    vector<int> genusNodes;  //indexes in phyloTree where genus' are located
    
    int alignmentLength;
    void readProbFile(ifstream&, ifstream&);

    
};

/**************************************************************************************************/


#endif /* opticlassifier_hpp */
