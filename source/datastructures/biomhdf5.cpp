//
//  biomhdf5.cpp
//  Mothur
//
//  Created by Sarah Westcott on 10/26/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "biomhdf5.hpp"

/**************************************************************************************************/
BiomHDF5::BiomHDF5(string fname, string l) : Biom("Biological Observation Matrix 2.1.0"){
    try {
        label = l;
        numOTUs = 0; numSamples = 0;
        
        read(fname);
    }
    catch(exception& e) {
        m->errorOut(e, "BiomHDF5", "BiomHDF5");
        exit(1);
    }
}
/**************************************************************************************************/
BiomHDF5::BiomHDF5() : Biom("Biological Observation Matrix 2.1.0"){
    try {
        numOTUs = 0; numSamples = 0;
       
        
    }
    catch(exception& e) {
        m->errorOut(e, "BiomHDF5", "BiomHDF5");
        exit(1);
    }
}
/**************************************************************************************************/
void BiomHDF5::read(string fname){
    try {
        nnz = 0; maxLevel = 0;
        otuNames.clear(); sampleNames.clear(); taxonomy.clear(); otuTaxonomies.clear();
        
#ifdef USE_HDF5
        
        Picrust* picrust; vector<string> metadata;
        printShared("/Users/sarahwestcott/desktop/release/temp.biom", metadata, picrust );
        
        H5::H5File file( fname.c_str(), H5F_ACC_RDONLY );
        //H5::Group     what(file.openGroup( "/" ));

        readAttributes(file); if (m->getControl_pressed()) { return; }

        try {
            //read otu names
            readNames(file, "observation/"); if (m->getControl_pressed()) { return; }
    
        }catch(H5::Exception& e){ //do nothing taxonomy info does not exist
            m->mothurOut("[ERROR]: Missing /""observation/ids/"" needed for OTU names.\n");
            m->setControl_pressed(true);
        }
        
        try {
            //read group names
            readNames(file, "sample/"); if (m->getControl_pressed()) { return; }
    
        }catch(H5::Exception& e){ //do nothing taxonomy info does not exist
            m->mothurOut("[ERROR]: Missing /""sample/ids/"" needed for group names.\n");
            m->setControl_pressed(true);
        }
        
        bool hasConsTaxonomy = false;
        try {
            //read otu taxonomies
            readTaxonomy(file); if (m->getControl_pressed()) { return; }
            
            if (otuTaxonomies.size() == otuNames.size()) { hasConsTaxonomy = true; }
    
        }catch(H5::Exception& e){ //do nothing taxonomy info does not exist
            m->mothurOut("\nIgnore HDF5 errors, mothur was checking for OTU taxonomies and your file does not contain them. They are not required, continuing.\n");
            hasConsTaxonomy = false;
        }
        
        try {
            //read otu taxonomies
            readOTUAbundances(file); if (m->getControl_pressed()) { return; }
    
        }catch(H5::Exception& e){ //do nothing taxonomy info does not exist
            m->mothurOut("[ERROR]: Missing /""observation/matrix/"" needed for OTU abundances.\n");
        }
        
        bool error = false;
        if (nnz != otudata.size()) { error = true; }

        //create shared file
        sort(sampleNames.begin(), sampleNames.end());
        
        if (shared != NULL) { delete shared; }
        
        shared = new SharedRAbundVectors();

        //create empty sharedrabundvectors so we can add otus below
        for (int i = 0; i < sampleNames.size(); i++) {
            SharedRAbundVector* temp = new SharedRAbundVector();
            temp->setGroup(sampleNames[i]);
            shared->push_back(temp);
        }
        
        if (matrixElementType == "float") {
            
            if (sharedFloat != NULL) { delete sharedFloat; }
            sharedFloat = new SharedRAbundFloatVectors();
        
            //creates new sharedRAbunds
            for (int i = 0; i < sampleNames.size(); i++) {
                SharedRAbundFloatVector* temp = new SharedRAbundFloatVector(shared->getNumBins()); //sets all abunds to 0
                temp->setLabel(label);
                temp->setGroup(sampleNames[i]);
                sharedFloat->push_back(temp);
            }
            sharedFloat->setLabels(label);
        }

        shared->setLabels(label);
       
        //for each otu
        int count = 0;
        for (int h = 0; h < indptr.size()-1; h++) {
            int otuStart = indptr[h];
            int otuEnd = indptr[h+1];
            
            vector<int> otuAbunds; otuAbunds.resize(sampleNames.size(), 0); //initialze otus sample abundances to 0 - only non zero abunds are recorded
            
            vector<float> otuFloatAbunds; otuFloatAbunds.resize(sampleNames.size(), 0); //initialze otus sample abundances to 0 - only non zero abunds are recorded

            
            for (int i = otuStart; i < otuEnd; i++) {
                otuAbunds[indices[i]] = (int)otudata[count];
                otuFloatAbunds[indices[i]] = otudata[count];
                count++;
            }
            
            shared->push_back(otuAbunds, otuNames[h]);
            if (matrixElementType == "float") { sharedFloat->push_back(otuFloatAbunds, otuNames[h]); }
        }

        
        
        if (hasConsTaxonomy) {
            
            if (shared->getNumBins() == otuTaxonomies.size()) {
                for (int i = 0; i < otuTaxonomies.size(); i++) {
                    if (m->getControl_pressed()) { break; }
                    string thisOTUsTax = otuTaxonomies[i];
                    string newTax = util.addUnclassifieds(thisOTUsTax, maxLevel, false);
                    Taxonomy thisOTUsTaxonomy(otuNames[i], newTax, shared->getOTUTotal(i));
                    consTax.push_back(thisOTUsTaxonomy);
                }
            }
        }
        
        if (sampleNames.size() == taxonomy.size()) {
            for (int i = 0; sampleNames.size(); i++) {
                if (m->getControl_pressed()) { break; }
                groupTaxonomies[sampleNames[i]] = taxonomy[i];
            }
        }
        
        file.close();

#endif
        
    }
    catch(exception& e) {
        m->errorOut(e, "BiomHDF5", "read");
        exit(1);
    }
}
#ifdef USE_HDF5
//**********************************************************************************************************************
//Group = "observation/" or "sample/", datasetName = "ids"
void BiomHDF5::readNames( H5::H5File& file, H5std_string groupName) {
    try {
        
        H5std_string datasetName = "ids";
        H5::Group group(file.openGroup(groupName));
        
        hsize_t numObjects = group.getNumObjs();
        
        if (numObjects != 0) { //we have this group
            
            H5::DataSet dataset = group.openDataSet(datasetName);
            H5::DataType  dataType(dataset.getDataType());
            H5::DataSpace dataSpace = dataset.getSpace();
            
            int rank = dataSpace.getSimpleExtentNdims(); //number of dimensions, should be 1
            hsize_t dims[rank]; dataSpace.getSimpleExtentDims(dims); //size of each dimension
            
            char **data = new char*[dims[0]];
            H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
            dataset.read((void*)data, str_type);
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: " + datasetName + " = "); }
            
            for (int i = 0; i < dims[0]; i++) {
                
                if (m->getDebug()) { m->mothurOut(toString(data[i]) + "\t"); }
                
                if (groupName == "observation/") { otuNames.push_back(data[i]); }
                else if (groupName == "sample/") { sampleNames.push_back(data[i]); }
               
                delete[] data[i];
            }  if (m->getDebug()) { m->mothurOutEndLine(); }
            delete[] data;
            
            dataset.close();
        }
        
        group.close();
    }
    catch(exception& e) {
        m->errorOut(e, "BiomHDF5", "readNames");
        exit(1);
    }
}
//**********************************************************************************************************************
//Group = "observation/metadata", datasetName = "taxonomy"
void BiomHDF5::readTaxonomy( H5::H5File& file) {
    try {
        
        H5std_string groupName = "observation/metadata";
        H5std_string datasetName = "taxonomy";
        H5::Group group(file.openGroup(groupName));
        
        hsize_t numObjects = group.getNumObjs();
        
        if (numObjects != 0) { //we have this group
            
            H5::DataSet dataset = group.openDataSet(datasetName);
            H5::DataType  dataType(dataset.getDataType());
            H5::DataSpace dataSpace = dataset.getSpace();
            
            int rank = dataSpace.getSimpleExtentNdims(); //number of dimensions
            hsize_t dims[rank]; dataSpace.getSimpleExtentDims(dims); //size of each dimension
            
            if (dataset.getTypeClass() == H5T_STRING) {
                if (rank == 1) {
                    char **data = new char*[dims[0]];
                    H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
                    dataset.read((void*)data, str_type);
                    
                    if (m->getDebug()) { m->mothurOut("[DEBUG]: " + datasetName + " = "); }
                    
                    for (int i = 0; i < dims[0]; i++) {
                        if (m->getDebug()) { m->mothurOut(toString(data[i]) + "\t"); }
                        
                        taxonomy.push_back(data[i]);
                        
                        delete[] data[i];
                    }  if (m->getDebug()) { m->mothurOutEndLine(); }
                    delete[] data;
                }else if (rank == 2) {
                    
                    char **data = new char*[dims[0]*dims[1]];

                    H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
                    dataset.read((void*)data, str_type);
                    
                    string otuTaxonomy = ""; int count = 0;
                    for (int i = 0; i < dims[0]*dims[1]; i++) {
                        
                        otuTaxonomy += data[i];  otuTaxonomy += ";";  count++;
                        
                        if (count == dims[1]) {
                            if (m->getDebug()) { m->mothurOut("[DEBUG]: " + toString(otuTaxonomy) + "\n"); }
                            otuTaxonomies.push_back(otuTaxonomy);
                            if (count > maxLevel) { maxLevel = count; }
                            count = 0; otuTaxonomy = "";
                        }
                    }
                }
            }
            
            dataset.close();
        }
        
        group.close();
    }
    catch(exception& e) {
        m->errorOut(e, "BiomHDF5", "readTaxonomy");
        exit(1);
    }
}
//**********************************************************************************************************************
//Group = "observation/matrix", datasetName = "taxonomy"
void BiomHDF5::readOTUAbundances( H5::H5File& file) {
    try {
        
        H5std_string groupName = "observation/matrix/";
        
        vector<string> datasets;
        datasets.push_back("data"); //otu abundances for each non zero abundnace entry - 1, 5, 1, 2, 3, 1, 1, 4, 2, 2, 1, 1, 1, 1, 1
        datasets.push_back("indices"); //index of group - maps into samples/ids  2, 0, 1, 3, 4, 5, 2, 3, 5, 0, 1, 2, 5, 1, 2
        datasets.push_back("indptr"); //maps non zero abundance to OTU - 0, 1, 6, 9, 13, 15 - 0 start of OTU1s indexes, 1 start of OTU2s indexes, ... 15 start of OTU5s indexes
        
        /*
         label group numOtus GG_OTU_1  GG_OTU_2  GG_OTU_3  GG_OTU_4  GG_OTU_5
         userLabel  Sample1     0           5       0           2       0
         userLabel  Sample2     0           1       0           2       1
         userLabel  Sample3     1           0       1           1       1
         userLabel  Sample4     0           2       4           0       0
         userLabel  Sample5     0           3       0           0       0
         userLabel  Sample6     0           1       2           1       0
         */

        H5::Group group(file.openGroup(groupName));
        
        for (int h = 0; h < datasets.size(); h++) {
            
            H5std_string datasetName = datasets[h];
            hsize_t numObjects = group.getNumObjs();
        
            if (numObjects != 0) { //we have this dataset
            
                H5::DataSet dataset = group.openDataSet(datasetName);
                H5::DataType  dataType(dataset.getDataType());
                H5::DataSpace dataSpace = dataset.getSpace();
            
                int rank = dataSpace.getSimpleExtentNdims(); //number of dimensions, should be 1
                hsize_t dims[rank]; dataSpace.getSimpleExtentDims(dims); //size of each dimension
            
                if (dataset.getTypeClass() == H5T_INTEGER) {
                    matrixElementType = "int";
                    if (rank == 1) {
                        int* data = new int[dims[0]];
                        H5::DataSpace data_mspace(rank, dims);
                        dataset.read(data, H5::PredType::NATIVE_INT, data_mspace, dataSpace);
                    
                        if (m->getDebug()) { m->mothurOut("[DEBUG]: " + datasetName + " = "); }
                        for (int i = 0; i < dims[0]; i++) {
                            if (m->getDebug()) { m->mothurOut(toString(data[i]) + "\t"); }
                        
                        
                            if (datasetName == "data") { otudata.push_back(data[i]); }
                            else if (datasetName == "indices") { indices.push_back(data[i]); }
                            else if (datasetName == "indptr") { indptr.push_back(data[i]); }
                        
                        } if (m->getDebug()) { m->mothurOutEndLine(); }
                        delete[] data;
                    }
                }else if (dataset.getTypeClass() == H5T_FLOAT) {
                    matrixElementType = "float";
                    if (rank == 1) {
                        float* data = new float[dims[0]];
                        H5::DataSpace data_mspace(rank, dims);
                        dataset.read(data, H5::PredType::NATIVE_FLOAT, data_mspace, dataSpace);
                        
                        if (m->getDebug()) { m->mothurOut("[DEBUG]: " + datasetName + " = "); }
                        for (int i = 0; i < dims[0]; i++) {
                            if (m->getDebug()) { m->mothurOut(toString(data[i]) + "\t"); }
                            
                            if (groupName == "observation/matrix/") {
                                if (datasetName == "data") { otudata.push_back(data[i]); }
                                else if (datasetName == "indices") { indices.push_back((int)data[i]); }
                                else if (datasetName == "indptr") { indptr.push_back((int)data[i]); }
                            }
                        }  if (m->getDebug()) { m->mothurOutEndLine(); }
                        delete[] data;
                    }
                }
            
                dataset.close();
            }
        }
        
        group.close();
    }
    catch(exception& e) {
        m->errorOut(e, "BiomHDF5", "readOTUAbundances");
        exit(1);
    }
}
//**********************************************************************************************************************
//read attribute of group
void BiomHDF5::readStringAttributes(H5::Group& fileAttributes, string name) {
    try {
        H5::Attribute attribute(fileAttributes.openAttribute(name));
        H5std_string attributeName; attribute.getName(attributeName);
        H5::DataType  attributeType(attribute.getDataType());
        
        // Read the Attribute Data. Depends on the kind of data
        if (attributeType.getClass() == H5T_STRING) {
            
            H5std_string value; attribute.read(attributeType, value);
            
            if (name == "id")                   { tableID = value;  }
            else if (name == "type")            { tableType = value; }
            else if (name == "format-url")      { formatURL = value; }
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: " + attributeName + " = " + value + "\n");  }
        }
        
        attribute.close();
    }
    catch(exception& e) {
        m->errorOut(e, "BiomHDF5", "readStringAttributes");
        exit(1);
    }
}
//**********************************************************************************************************************
//read attribute of group
void BiomHDF5::readIntAttributes(H5::Group& fileAttributes, string name) {
    try {
        H5::Attribute attribute(fileAttributes.openAttribute(name));
        H5std_string attributeName; attribute.getName(attributeName);
        H5::DataType  attributeType(attribute.getDataType());
        H5::DataSpace attDataSpace = attribute.getSpace();
        
        int rank = attDataSpace.getSimpleExtentNdims(); //number of dimensions
        hsize_t dims[rank]; attDataSpace.getSimpleExtentDims(dims); //size of each dimension
        
        // Read the Attribute Data. Depends on the kind of data
        if (attributeType.getClass() == H5T_INTEGER) {

            if (attDataSpace.isSimple()) {
                
                if (rank == 0) {
                    hsize_t data = 0;
                    attribute.read(attributeType, &data);
                    if (m->getDebug()) { m->mothurOut("[DEBUG]: " + attributeName + " = " + toString(data) + "\n");  }
                    if (attributeName == "nnz") { nnz = data; }
                }else if (rank == 1) {
                    hsize_t data[dims[0]];
                    attribute.read(attributeType, data);
                    if (m->getDebug()) { m->mothurOut("[DEBUG]: " + attributeName + " = "); }
                    for (int i = 0; i < dims[0]; i++) {
                        if (m->getDebug()) { m->mothurOut(toString(data[i]) + "\t"); }
                        if (attributeName == "nnz") { nnz = data[i]; }
                    }  if (m->getDebug()) { m->mothurOutEndLine(); }
                    
                    if (attributeName == "shape") {
                        if (dims[0] == 2) { numOTUs = data[0]; numSamples = data[1]; }
                    }
                }
            }
            
        }
        
        attribute.close();
    }
    catch(exception& e) {
        m->errorOut(e, "BiomHDF5", "readIntAttributes");
        exit(1);
    }
}
//**********************************************************************************************************************
//Process attribute of group or dataset
void BiomHDF5::readAttributes(H5::H5File& file) {
    try {
        
        H5::Group  fileAttributes(file.openGroup( "/" ));
        
        //read table id
        readStringAttributes(fileAttributes, "id");
        
        //read table type
        readStringAttributes(fileAttributes, "type");
        
        //read format-url
        readStringAttributes(fileAttributes, "format-url");
        
        //read shape
        readIntAttributes(fileAttributes, "shape");
        
        //read number non zero
        readIntAttributes(fileAttributes, "nnz");
        
        fileAttributes.close();
    }
    catch(exception& e) {
        m->errorOut(e, "BiomHDF5", "readAttributes");
        exit(1);
    }
}
#endif
//**********************************************************************************************************************
void BiomHDF5::printShared(string outputFileName, vector<string> sampleMetadata, Picrust* picrust) {
    try {
        //set required fields
        set<string> requiredTopLevelAttrib;
        requiredTopLevelAttrib.insert("id"); requiredTopLevelAttrib.insert("type"); requiredTopLevelAttrib.insert("format-url");
        requiredTopLevelAttrib.insert("format-version"); requiredTopLevelAttrib.insert("generated-by"); requiredTopLevelAttrib.insert("creation-date");
        requiredTopLevelAttrib.insert("shape"); requiredTopLevelAttrib.insert("nnz");
        
        //set required datasets - groupname -> datasetname
        //"observation/ids" -> otuLabels - "GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5
        
        //"observation/matrix/data" -> otu abundances for each non zero abundnace entry - 1, 5, 1, 2, 3, 1, 1, 4, 2, 2, 1, 1, 1, 1, 1
        //"observation/matrix/indices" -> index of group - maps into samples/ids  2, 0, 1, 3, 4, 5, 2, 3, 5, 0, 1, 2, 5, 1, 2
        //"observation/matrix/indptr" -> maps non zero abundance to OTU - 0, 1, 6, 9, 13, 15 - 0 start of OTU1s indexes, 1 start of OTU2s indexes, ... 15 start of OTU5s indexes
        
        /*
         label group numOtus GG_OTU_1  GG_OTU_2  GG_OTU_3  GG_OTU_4  GG_OTU_5
         userLabel  Sample1     0           5       0           2       0
         userLabel  Sample2     0           1       0           2       1
         userLabel  Sample3     1           0       1           1       1
         userLabel  Sample4     0           2       4           0       0
         userLabel  Sample5     0           3       0           0       0
         userLabel  Sample6     0           1       2           1       0
         */
        
        //"sample/ids" -> group names - "Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"
        //"observation/metadata/taxonomy" -> taxonomy info - otu classifications
        //"sample/metadata/" -> group metadata (optional)
        
        int   RANK = 2;
        
        #ifdef USE_HDF5
           
        
            H5::H5File file(outputFileName.c_str(), H5F_ACC_TRUNC );
        
        H5::Group group = file.createGroup("/");
        
        // Create new dataspace for attribute
        H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
        
        H5::Attribute idAttribute = group.createAttribute("id", H5T_STRING, attr_dataspace);
        

            // Create new string datatype for attribute
            H5::StrType strdatatype(H5::PredType::C_S1, 256); // of length 256 characters

            // Set up write buffer for attribute
            const H5std_string strwritebuf ("a field that can be used to id a table (or null)");

            // Create attribute and write to it
        //H5::Attribute myatt_in = dataset.createAttribute(ATTR_NAME, strdatatype, attr_dataspace);
        idAttribute.write(strdatatype, strwritebuf);

            // Set up read buffer for attribute
            H5std_string strreadbuf ("");

            // Open attribute and read its contents
        H5::Group     what(file.openGroup( "/" ));
        H5::Attribute myatt_out = what.openAttribute("id");
            myatt_out.read(strdatatype, strreadbuf);

            // Display attribute contents
            cout << "Attribute contents: " << strreadbuf << endl;

        
            //int** data = new int*[numOtus*numGroups];
            
              
              // Define size of otudata
             // hsize_t     dimsf[2];              // dataset dimensions
             // dimsf[0] = numGroups;
             // dimsf[1] = numOtus;
        
             // H5::DataSpace dataspace( RANK, dimsf );
        
              //Initialize space for otu data
              
        
        //

        
        #endif
    }
    catch(exception& e) {
        m->errorOut(e, "BiomHDF5", "checkGroups");
        exit(1);
    }
}
//**********************************************************************************************************************
void BiomHDF5::print(string outputFileName, vector<string> sampleMetadata, Picrust* picrust) {
    try {
        
        if (matrixElementType == "int") {
            printShared(outputFileName, sampleMetadata, picrust);
        }else {
            printFloat(outputFileName, sampleMetadata, picrust);
        }
    }
    catch(exception& e) {
        m->errorOut(e, "BiomHDF5", "print");
        exit(1);
    }
}
/**************************************************************************************************/
