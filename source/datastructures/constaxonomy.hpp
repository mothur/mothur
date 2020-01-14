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

class OTUTaxonomy {
    
public:
    
    OTUTaxonomy();
    OTUTaxonomy(string, string, int);
    OTUTaxonomy(ifstream&);
    ~OTUTaxonomy() {}
    
    void setName(string n)          { name = n;         }
    void setNumReps(int n)          { numReps = n;      }
    string getName()                { return name;      }
    vector<Taxon> getTaxons()       { return taxonomy;  }
    void setTaxons(vector<Taxon> t) { taxonomy = t;     }
    int getNumSeqs()                { return numReps;   }
    void setTaxons(string);
    
    string getInlineConsTaxonomy();
    string getConsTaxString(bool includeConfidence=true); //pass in true to include confidences

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
