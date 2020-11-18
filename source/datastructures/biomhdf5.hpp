//
//  biomhdf5.hpp
//  Mothur
//
//  Created by Sarah Westcott on 10/26/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//


//#ifdef USE_HDF5
//#endif

#ifndef biomhdf5_hpp
#define biomhdf5_hpp

#include "biom.hpp"

//http://biom-format.org/documentation/format_versions/biom-2.1.html

class BiomHDF5 : public Biom {
    
public:
    
    BiomHDF5();
    BiomHDF5(string, string);
    ~BiomHDF5() {  }
    
    void read(string);
    void print(ofstream&, vector<string>, Picrust*) {}
    
    
    
private:
    int nnz;
   
    set<string> requiredTopLevelAttrib;
    map<string, vector<string> > requiredOTUDatasets;
    map<string, vector<string> > requiredSampleDatasets;
    map<string, vector<string> > optionalDatasets;
    
    vector<string> otuNames, sampleNames, taxonomy, otuTaxonomies;
    vector<int> indices, indptr;
    vector<float> otudata;

#ifdef USE_HDF5
    void processAttributes(H5::Group&, set<string>&);
    void checkGroups(H5::H5File&, map<string, vector<string> >&);
#endif
    
};

#endif /* biomhdf5_hpp */
