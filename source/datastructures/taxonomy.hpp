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


/**************************************************************************************************/

class Taxonomy {
    
public:
    
    Taxonomy();
    Taxonomy(string, string, int); //name, tax, abund
    Taxonomy(string, string);
    Taxonomy(ifstream&);
    ~Taxonomy() {}
    
    void setName(string n)          { name = n;         }
    void setNumSeqs(int n)          { numReps = n;      }
    string getName()                { return name;      }
    vector<Taxon> getTaxons()       { return taxonomy;  }
    vector<string> getSimpleTaxons (bool includeConfidence=false);
    void setTaxons(vector<Taxon> t) { taxonomy = t;     }
    int getNumSeqs()                { return numReps;   }
    int getNumLevels()              { return taxonomy.size();   }
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
    bool containsConfidence;
    vector<Taxon> taxonomy;
    Utils util;
    
    vector<Taxon> parseTax(string);
};

/**************************************************************************************************/

#endif /* constaxonomy_hpp */
