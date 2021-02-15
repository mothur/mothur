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
        printShared("/Users/swestcott/desktop/release/temp.biom", metadata, picrust );
        
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
//**********************************************************************************************************************
//print required dataset attributes
void BiomHDF5::printRequiredFileAttributes(H5::H5File& file, int numBins, int numSamples) {
    try {
        H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR); // Create new dataspace for attribute
        H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE);
    
        H5std_string idValue(tableID);
        H5::Attribute idAttribute = file.createAttribute("id", strdatatype, attr_dataspace);
        idAttribute.write(strdatatype, idValue);
    
        H5std_string typeValue(tableType);
        H5::Attribute typeAttribute = file.createAttribute("type", strdatatype, attr_dataspace);
        typeAttribute.write(strdatatype, typeValue);
    
        H5std_string formatUrl(formatURL);
        H5::Attribute urlAttribute = file.createAttribute("format-url", strdatatype, attr_dataspace);
        urlAttribute.write(strdatatype, formatUrl);
    
        H5std_string generatedByValue(mothurVersion);
        H5::Attribute generatedByAttribute = file.createAttribute("generated-by", strdatatype, attr_dataspace);
        generatedByAttribute.write(strdatatype, generatedByValue);
    
        time_t rawtime; struct tm * timeinfo;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        string dateString = asctime (timeinfo);
        int pos = dateString.find('\n');
        if (pos != string::npos) { dateString = dateString.substr(0, pos);}
    
        H5std_string dataValue(dateString);
        H5::Attribute dateAttribute = file.createAttribute("creation-date", strdatatype, attr_dataspace);
        dateAttribute.write(strdatatype, dataValue);

        const hsize_t dims=2;
        H5::DataType datatype = H5::ArrayType(H5::PredType::NATIVE_INT, 1, &dims); // Create new int datatype
        hsize_t data[dims]; data[0] = 2; data[1] = 1;
        H5::DataSpace intSpace(1, &dims);
    
        H5::Attribute formatVersionAttribute = file.createAttribute("format-version", datatype, intSpace);
        formatVersionAttribute.write(datatype, data);

        data[0] = numBins; data[1] = numSamples;
        H5::Attribute shapeAttribute = file.createAttribute("shape", datatype, intSpace);
        shapeAttribute.write(datatype, data);
    
        hsize_t nnzValue = nnz;
        H5::DataType nnzDatatype = H5::PredType::NATIVE_INT;
        H5::Attribute nnzAttribute = file.createAttribute("nnz", nnzDatatype, attr_dataspace);
        nnzAttribute.write(nnzDatatype, &nnzValue);
 
    }
    catch(exception& e) {
        m->errorOut(e, "BiomHDF5", "printRequiredAttributes");
        exit(1);
    }
}
//**********************************************************************************************************************
//print otuNames
//"observation/ids" -> otuLabels - "GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5
void BiomHDF5::printOTULabels(H5::Group& group, vector<string> otuNames) {
    try {
        hsize_t     dimsf[1]; dimsf[0] = otuNames.size();
        H5::DataSpace dataspace( 1, dimsf );
        H5::StrType datatype(H5::PredType::C_S1, H5T_VARIABLE); 
        
        //fill data with names
        char* data[dimsf[0]];
        for (int i = 0; i < otuNames.size(); i++) { data[i] = (char*) otuNames[i].c_str(); }
             
        const H5std_string  DATASET_NAME( "ids" );
        H5::DataSet dataset = group.createDataSet( DATASET_NAME, datatype, dataspace );
        
        dataset.write( data, datatype );
        dataset.close();
    }
    catch(exception& e) {
        m->errorOut(e, "BiomHDF5", "printOTULabels");
        exit(1);
    }
}
//**************************************************************************************************************
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

//group = "observation/matrix/";
void BiomHDF5::printOTUAbundances(H5::Group& group, int numBins, int numSamples, bool useRelabund=false) {
    try {
        
        int otuStartIndex = 0;
        vector<int> indptr, indices, abunds;
        vector<float> abundsFloat;
        
        //find numOTUs
        if (useRelabund)    { numBins = sharedFloat->getNumBins(); numSamples = sharedFloat->size();  }
        else                { numBins = shared->getNumBins(); numSamples = shared->size();            }
            
        //fill indices, indptr and data vectors
        for (int i = 0; i < numBins; i++) {
            
            if (m->getControl_pressed()) { return; }
            
            vector<int> thisOtusAbundances; vector<float> thisOtusFloatAbundances; float zero = 0.0;
            if (useRelabund) { thisOtusFloatAbundances = sharedFloat->getOTU(i);    }
            else             { thisOtusAbundances = shared->getOTU(i);              }
            
            indptr.push_back(otuStartIndex);
            for (int j = 0; j < numSamples; j++) {
                
                if (useRelabund) {
                    if (util.isEqual(thisOtusFloatAbundances[j], zero)) {} //skip zero values
                    else {
                        otuStartIndex++; //update number of non zero values for this OTU - use to create indptr values
                        indices.push_back(j); //index to sample providing this abund
                        abundsFloat.push_back(thisOtusFloatAbundances[j]); //save this samples OTU abundance
                    }
                }else {
                    if (thisOtusAbundances[j] == 0) {} //skip zero values
                    else {
                        otuStartIndex++; //update number of non zero values for this OTU - use to create indptr values
                        indices.push_back(j); //index to sample providing this abund
                        abunds.push_back(thisOtusAbundances[j]); //save this samples OTU abundance
                    }
                }
            }
        }
        
        // dataset dimensions
        hsize_t     dimsf[1]; dimsf[0] = otuStartIndex;
        H5::DataSpace dataspace( 1, dimsf ); //dataspace 1 x nnz
        H5::IntType datatype( H5::PredType::NATIVE_INT );
        
        int data[otuStartIndex];
        for (int i = 0; i < otuStartIndex; i++) { data[i] = indices[i]; } //fill data with indices
             
        const H5std_string  DATASET_NAME( "indices" );
        H5::DataSet dataset = group.createDataSet( DATASET_NAME, datatype, dataspace );
        
        dataset.write( data, H5::PredType::NATIVE_INT );
        dataset.close();
        
        //create data dataset - type depends on whether or not we are using the relabund values
        if (useRelabund) { //print float
            H5::FloatType datatypeFloat( H5::PredType::NATIVE_FLOAT );
            float dataFloat[otuStartIndex];
            
            for (int i = 0; i < otuStartIndex; i++) { dataFloat[i] = abundsFloat[i]; } //fill data with abunds
            const H5std_string  DATASET_NAME( "data" );
            H5::DataSet dataset = group.createDataSet( DATASET_NAME, datatypeFloat, dataspace );
            
            dataset.write( dataFloat, H5::PredType::NATIVE_FLOAT );
            dataset.close();
        }else { //print shared
            for (int i = 0; i < otuStartIndex; i++) { data[i] = abunds[i]; } //fill data with abunds
            const H5std_string  DATASET_NAME( "data" );
            H5::DataSet dataset = group.createDataSet( DATASET_NAME, datatype, dataspace );
            
            dataset.write( data, H5::PredType::NATIVE_INT );
            dataset.close();
        }
        
        //create indptr dataset
        dimsf[0] = numBins;
        H5::DataSpace dataspaceIndptr(1, dimsf); //dataspace 1 x numBins
        
        int dataIndptr[numBins];
        for (int i = 0; i < numBins; i++) { dataIndptr[i] = indptr[i]; } //fill data with indptr
             
        const H5std_string  DATASET_NAME_INDPTR( "indptr" );
        H5::DataSet datasetIndptr = group.createDataSet( DATASET_NAME_INDPTR, datatype, dataspaceIndptr );
        
        datasetIndptr.write(dataIndptr, H5::PredType::NATIVE_INT);
        datasetIndptr.close();
    }
    catch(exception& e) {
        m->errorOut(e, "BiomHDF5", "printOTUAbundances");
        exit(1);
    }
}
#endif
//**********************************************************************************************************************
void BiomHDF5::printShared(string outputFileName, vector<string> sampleMetadata, Picrust* picrust) {
    try {
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
        
        //"observation/metadata/taxonomy" -> taxonomy info - otu classifications
        
        //"sample/ids" -> group names - "Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"
        //"sample/metadata/" -> group metadata (optional)
        
        //run this first because if picrust alters the shared vector we will need to use the updated info
        //taxMetadata[0] = taxonomy for otu0
        vector< vector<string> > taxMetadata = getMetaData(picrust);
        
        //find number of non zero otus
        nnz = 0;
        for (int j = 0; j < shared->getNumBins(); j++) {
            vector<int> thisOTU = shared->getOTU(j);
            for (int i = 0; i < thisOTU.size(); i++) {
                if (thisOTU[i] != 0) { nnz++; }
            }
        }
        
    #ifdef USE_HDF5
           
        H5::H5File file(outputFileName.c_str(), H5F_ACC_TRUNC );
        
        //print required file attributes
        printRequiredFileAttributes(file, shared->getNumBins(), shared->size()); //id, type, format-url, format-version, generated-by, creation-date, shape, nnz
        
        H5::Group sampleGroup( file.createGroup( "sample/" ));
        
        //print otuLabels called "observation/ids" in biom file
        H5::Group observationGroup( file.createGroup( "observation/" ));
        printOTULabels(observationGroup, shared->getOTUNames());
        
        //print otuAbundances called "observation/matrix/" (data, indicies, indptr) in biom file
        H5::Group matrixGroup( file.createGroup( "observation/matrix/" ));
        printOTUAbundances(matrixGroup, shared->getNumBins(), shared->size());
        
       
    #endif
        
    }
    catch(exception& e) {
        m->errorOut(e, "BiomHDF5", "printShared");
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
//**********************************************************************************************************************
vector< vector<string> > BiomHDF5::getMetaData(Picrust* picrust, bool useRelabund){
    try {
        vector< vector<string> > metadata;
        
        if (consTax.size() == 0) {
            if (!useRelabund) {
                for (int i = 0; i < shared->getNumBins(); i++) {  vector<string> temp; temp.push_back("null"); metadata.push_back(temp);  } }
            else {
                for (int i = 0; i < sharedFloat->getNumBins(); i++) { vector<string> temp; temp.push_back("null"); metadata.push_back(temp);  }
            }
        }
        else {
            
            if (!useRelabund) {
                if (shared == NULL) { m->setControl_pressed(true); return metadata; }
            }else {
                if (sharedFloat == NULL) { m->setControl_pressed(true); return metadata; }
            }
            
            //should the labels be Otu001 or PhyloType001
            vector<string> otuNames;
            if (!useRelabund) { otuNames = shared->getOTUNames(); }
            else { otuNames = sharedFloat->getOTUNames(); }
            
            string firstBin = otuNames[0];
            string binTag = "Otu";
            if ((firstBin.find("Otu")) == string::npos) { binTag = "PhyloType";  }
            
            map<string, string> labelTaxMap;
            string snumBins = toString(otuNames.size());
            for (int i = 0; i < consTax.size(); i++) {
                
                if (m->getControl_pressed()) { return metadata; }
                
                string thisOtuLabel = consTax[i].getName();
                
                //if there is a bin label use it otherwise make one
                if (util.isContainingOnlyDigits(thisOtuLabel)) {
                    string binLabel = binTag;
                    string sbinNumber = thisOtuLabel;
                    if (sbinNumber.length() < snumBins.length()) {
                        int diff = snumBins.length() - sbinNumber.length();
                        for (int h = 0; h < diff; h++) { binLabel += "0"; }
                    }
                    binLabel += sbinNumber;
                    binLabel = util.getSimpleLabel(binLabel);
                    labelTaxMap[binLabel] = consTax[i].getConsTaxString();
                }else {
                    map<string, string>::iterator it = labelTaxMap.find(util.getSimpleLabel(thisOtuLabel));
                    if (it == labelTaxMap.end()) {
                        labelTaxMap[util.getSimpleLabel(thisOtuLabel)] = consTax[i].getConsTaxString();
                    }else {
                        m->mothurOut("[ERROR]: Cannot add OTULabel " +  thisOtuLabel + " because it's simple label " + util.getSimpleLabel(consTax[i].getName()) + " has already been added and will result in downstream errors. Have you mixed mothur labels and non mothur labels? To make the files work well together and backwards compatible mothur treats 1, OTU01, OTU001, OTU0001 all the same. We do this by removing any non numeric characters and leading zeros. For eaxample: Otu000018 and OtuMY18 both map to 18.\n"); m->setControl_pressed(true);
                    }
                }
            }
            
            //sanity check for file issues - do you have the same number of bins in the shared and constaxonomy file
            if (!useRelabund) {
                if (shared->getNumBins() != labelTaxMap.size()) {
                    m->mothurOut("[ERROR]: Your constaxonomy file contains " + toString(labelTaxMap.size()) + " otus and your shared file contain " + toString(shared->getNumBins()) + " otus, cannot continue.\n"); m->setControl_pressed(true); return metadata;
                }
            }else {
                if (sharedFloat->getNumBins() != labelTaxMap.size()) {
                    m->mothurOut("[ERROR]: Your constaxonomy file contains " + toString(labelTaxMap.size()) + " otus and your shared file contain " + toString(sharedFloat->getNumBins()) + " otus, cannot continue.\n"); m->setControl_pressed(true); return metadata;
                }
            }
            
            //merges OTUs classified to same gg otuid, sets otulabels to gg otuids, averages confidence scores of merged otus.  overwritting of otulabels is fine because constaxonomy only allows for one label to be processed.  If this assumption changes, could cause bug.
            if (picrust != NULL) {
                if (!useRelabund)   { picrust->setGGOTUIDs(labelTaxMap, shared);        }
                else                { picrust->setGGOTUIDs(labelTaxMap, sharedFloat);   }
            }
            
            //{"taxonomy":["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriaceae", "g__Escherichia", "s__"]}
            
            //traverse the binLabels forming the metadata strings and saving them
            //make sure to sanity check
            map<string, string>::iterator it;
            vector<string> currentLabels; int numBins = 0;
            if (!useRelabund)   { currentLabels = shared->getOTUNames();    numBins = shared->getNumBins();      }
            else                { currentLabels = sharedFloat->getOTUNames();   numBins = sharedFloat->getNumBins();
            }
            
            for (int i = 0; i < numBins; i++) {
                
                if (m->getControl_pressed()) { return metadata; }
                
                it = labelTaxMap.find(util.getSimpleLabel(currentLabels[i]));
                
                if (it == labelTaxMap.end()) { m->mothurOut("[ERROR]: can't find taxonomy information for " + currentLabels[i] + ".\n"); m->setControl_pressed(true); }
                else {
                    vector<string> scores;
                    vector<string> taxonomies = util.parseTax(it->second, scores);
                    
                    metadata.push_back(taxonomies);
                }
            }
        }
        
        return metadata;
    }
    catch(exception& e) {
        m->errorOut(e, "BiomHDF5", "getMetadata");
        exit(1);
    }

}
/**************************************************************************************************/
