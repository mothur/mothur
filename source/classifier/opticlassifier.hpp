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
    
    OptiClassifier(string reftax, string reffasta, int cutoff, int iters, bool shortcuts, string mothurVersion, string filter);
    ~OptiClassifier() {}
    
    string getTaxonomy(Sequence*, string&, bool&); 
    
private:
    
    map<char, int> baseMap; map<char, int>::iterator itBase;
    map<int, float> baseProbs; map<int, float>::iterator itBaseProbs; //char -> prob where char is converted to int
    vector<int> allCols;

    int numBases, numFilteredColumns;
    string filter;
    vector< vector< vector<float> > > charGenusProb;    //[locationInAlignment][base][genus]
                                        //charGenusProb[0][3][392] = probability that a sequence within genus that's index in the tree is 392 would contain 'C' at alignment position 0; charGenusProb[23][0][392] = probability that a sequence within genus that's index in the tree is 392 would contain 'A' at alignment position 23;
    
    vector<int> genusTotals; //number of sequence at each genus
    vector<int> genusNodes;  //indexes in phyloTree where genus' are located
    
    vector< map< int, float> > reversedProbs; //reversedProbs[0][0] = Probability of 'A' being at alignment location 0 in the reference. reversedProbs[0][1] = Probability of 'T' being at alignment location 0 in the reference.
    map<int, int> reverse;
    
    int confidenceThreshold, iters;
    void readProbFile(ifstream&, ifstream&);
    bool isReversed(vector<int>&);
    string bootstrapResults(vector<int>&, int, string&);
    int getMostProbableTaxonomy(vector<int>&, vector<int>&);


    
};

/**************************************************************************************************/


#endif /* opticlassifier_hpp */
