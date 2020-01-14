//
//  constaxonomy.hpp
//  Mothur
//
//  Created by Sarah Westcott on 1/13/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#ifndef constaxonomy_hpp
#define constaxonomy_hpp

#include "mothurout.h"
#include "utils.hpp"
#include "writer.h"

/***********************************************************************/
struct Taxon {
    string name;
    float confidence;

    Taxon(string n, float conf) : name(n), confidence(conf) {}
    ~Taxon(){}
};

/**************************************************************************************************/

class ConsTaxonommy {
    
public:
    
    ConsTaxonommy();
    ConsTaxonommy(string, string, int);
    ConsTaxonommy(ifstream&);
    ~ConsTaxonommy() {}
    
    void setName(string n)          { name = n;         }
    void setNumReps(int n)          { numReps = n;      }
    string getName()                { return name;      }
    vector<Taxon> getTaxons()       { return taxonomy;  }
    void setTaxons(vector<Taxon> t) { taxonomy = t;     }
    int getNumSeqs()                { return numReps;   }
    void setTaxons(string);
    
    string getInlineConsTaxonomy();
    string getConsTaxString(bool); //pass in true to include confidences

    void printConsTax(ostream&);
    void printConsTax(OutputWriter*);
    void printConsTaxNoConfidence(ostream&);
    
protected:
    
    MothurOut* m;
    string name;
    int numReps;
    vector<Taxon> taxonomy;
    Utils util;
    
    vector<Taxon> parseTax(string);
};

/**************************************************************************************************/

#endif /* constaxonomy_hpp */
