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

/* Required Groups
observation/               : The HDF5 group that contains observation specific information and an observation oriented view of the data
observation/matrix         : The HDF5 group that contains matrix data oriented for observation-wise operations (e.g., in compressed sparse row format)
observation/metadata       : The HDF5 group that contains observation specific metadata information
observation/group-metadata : The HDF5 group that contains observation specific group metadata information (e.g., phylogenetic tree)
sample/                    : The HDF5 group that contains sample specific information and a sample oriented data oriented view of the data
sample/matrix              : The HDF5 group that contains matrix data oriented for sample-wise operations (e.g., in compressed sparse column format)
sample/metadata            : The HDF5 group that contains sample specific metadata information
sample/group-metadata      : The HDF5 group that contains sample specific group metadata information (e.g., relationships between samples)
 */

/* Required Datasets
observation/ids            : <string> or <variable length string> A (N,) dataset of the observation IDs, where N is the total number of IDs
observation/matrix/data    : <float64> A (nnz,) dataset containing the actual matrix data
observation/matrix/indices : <int32> A (nnz,) dataset containing the column indices (e.g., maps into samples/ids)
observation/matrix/indptr  : <int32> A (M+1,) dataset containing the compressed row offsets
sample/ids                 : <string> or <variable length string> A (M,) dataset of the sample IDs, where M is the total number of IDs
sample/matrix/data         : <float64> A (nnz,) dataset containing the actual matrix data
sample/matrix/indices      : <int32> A (nnz,) dataset containing the row indices (e.g., maps into observation/ids)
sample/matrix/indptr       : <int32> A (N+1,) dataset containing the compressed column offsets
 */

/* Required Attributes
id                   : <string or null> a field that can be used to id a table (or null)
type                 : <string> Table type (a controlled vocabulary)
                       Acceptable values:
                        "OTU table"
                        "Pathway table"
                        "Function table"
                        "Ortholog table"
                        "Gene table"
                        "Metabolite table"
                        "Taxon table"
format-url           : <url> A string with a static URL providing format details
format-version       : <tuple> The version of the current biom format, major and minor
generated-by         : <string> Package and revision that built the table
creation-date        : <datetime> Date the table was built (ISO 8601 format)
shape                : <list of ints>, the number of OTUs (rows) and number of Samples (cols) in data
nnz                  : <int> The number of non-zero elements in the table
*/

class BiomHDF5 : public Biom {
    
public:
    
    BiomHDF5();
    BiomHDF5(string, string);
    ~BiomHDF5() {  }
    
    void read(string);
    void print(string, vector<string>, Picrust*);
    
private:
    int nnz, numOTUs, numSamples;
    vector<string> otuNames, sampleNames, taxonomy, otuTaxonomies;
    vector<int> indices, indptr;
    vector<float> otudata;
    
    void printShared(string, vector<string>, Picrust*);
    void printFloat(string, vector<string>, Picrust*) {}
    
    vector< vector<string> > getMetaData(Picrust*, bool useRelabund=false);

#ifdef USE_HDF5
    void readAttributes(H5::H5File& file);
    void readStringAttributes(H5::Group& fileAttributes, string);
    void readIntAttributes(H5::Group& fileAttributes, string);
    //void checkGroups(H5::H5File&, map<string, vector<string> >&);
    
    void printRequiredFileAttributes(H5::H5File& file, int, int);
    void printOTULabels(H5::Group& group, vector<string>);
    void printOTUAbundances(H5::Group& group, int, int, bool);
    
    void readNames( H5::H5File& file, H5std_string groupName);
    void readTaxonomy( H5::H5File& file);
    void readOTUAbundances( H5::H5File& file);
    
    
#endif
    
};

#endif /* biomhdf5_hpp */
