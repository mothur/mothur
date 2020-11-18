//
//  biomsimple.hpp
//  Mothur
//
//  Created by Sarah Westcott on 10/26/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#ifndef biomsimple_hpp
#define biomsimple_hpp

//biom version 0.9.1

#include "biom.hpp"
#include "picrust.hpp"

class BiomSimple : public Biom {
    
public:
    
    BiomSimple();
    BiomSimple(string, string); // filename, label
    ~BiomSimple() {  }
    
    void read(string);
    void printHeading(ofstream&, string, string); //fileHandle, mothurVersion, sharedfileName
    void print(ofstream&, vector<string>, Picrust*); //filehandle, metatdata, picrust

    
private:
   
    //examples: tableType = "OTU table", matrixFormat = "sparse" or "dense", 
    string matrixFormat, tableType;
    
    string getTag(string&);
    void getDims(string, int&, int&);

    SharedRAbundVectors* extractOTUData(string, vector<string>&, int);
    vector< vector<string> > extractTaxonomyData(string, int&, bool&);
    vector<string> getNamesAndTaxonomies(string);
    string getName(string);
    string getTaxonomy(string, string);
    vector<string> getMetaData(Picrust*);
    
};


#endif /* biomsimple_hpp */
