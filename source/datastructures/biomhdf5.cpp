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
        requiredTopLevelAttrib.clear(); requiredOTUDatasets.clear(); requiredSampleDatasets.clear(); optionalDatasets.clear();
        otuNames.clear(); sampleNames.clear(); taxonomy.clear(); otuTaxonomies.clear();
        
        //set required fields
        requiredTopLevelAttrib.insert("id"); requiredTopLevelAttrib.insert("type"); requiredTopLevelAttrib.insert("format-url");
        requiredTopLevelAttrib.insert("format-version"); requiredTopLevelAttrib.insert("generated-by"); requiredTopLevelAttrib.insert("creation-date");
        requiredTopLevelAttrib.insert("shape"); requiredTopLevelAttrib.insert("nnz");
        
        //set required datasets - groupname -> datasetname
        vector<string> datasets; datasets.push_back("ids");
        requiredOTUDatasets["observation/"] = datasets; //otuLabels - "GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5"
        datasets.clear();
        datasets.push_back("data"); //otu abundances for each non zero abundnace entry - 1, 5, 1, 2, 3, 1, 1, 4, 2, 2, 1, 1, 1, 1, 1
        datasets.push_back("indices"); //index of group - maps into samples/ids  2, 0, 1, 3, 4, 5, 2, 3, 5, 0, 1, 2, 5, 1, 2
        datasets.push_back("indptr"); //maps non zero abundance to OTU - 0, 1, 6, 9, 13, 15 - 0 start of OTU1s indexes, 1 start of OTU2s indexes, ... 15 start of OTU5s indexes
        requiredOTUDatasets["observation/matrix/"] = datasets;
        
        datasets.clear(); datasets.push_back("ids"); //group names - "Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"
        requiredSampleDatasets["sample/"] = datasets;
        datasets.clear(); datasets.push_back("taxonomy"); //optional datasets to look for - taxonomy info - otu classifications
        optionalDatasets["observation/metadata"] = datasets;

        /*
         label group numOtus GG_OTU_1  GG_OTU_2  GG_OTU_3  GG_OTU_4  GG_OTU_5
         userLabel  Sample1     0           5       0           2       0
         userLabel  Sample2     0           1       0           2       1
         userLabel  Sample3     1           0       1           1       1
         userLabel  Sample4     0           2       4           0       0
         userLabel  Sample5     0           3       0           0       0
         userLabel  Sample6     0           1       2           1       0
         */
        
        
        
#ifdef USE_HDF5
        H5::H5File file( fname.c_str(), H5F_ACC_RDONLY );
        H5::Group     what(file.openGroup( "/" ));

        processAttributes(what, requiredTopLevelAttrib); if (m->getControl_pressed()) { return; }

        try {
    
            checkGroups(file, requiredOTUDatasets); if (m->getControl_pressed()) { return; }
    
        }catch(H5::Exception& e){ //do nothing taxonomy info does not exist
            m->mothurOut("[ERROR]: Missing required groups or datasets, aborting. Required datasets include: \n");
            for (map<string, vector< string> >::iterator it = requiredOTUDatasets.begin(); it != requiredOTUDatasets.end(); it++) {
                for (int i = 0; i < (it->second).size(); i++) { m->mothurOut(it->first + '\t' + it->second[i] + '\n'); }
            }
            m->mothurOutEndLine(); m->setControl_pressed(true);
        }

        try {
    
            checkGroups(file, requiredSampleDatasets); if (m->getControl_pressed()) { return; }
    
        }catch(H5::Exception& e){ //do nothing taxonomy info does not exist
            m->mothurOut("[ERROR]: Missing required groups or datasets, aborting. Required datasets include: \n");
            for (map<string, vector< string> >::iterator it = requiredOTUDatasets.begin(); it != requiredOTUDatasets.end(); it++) {
                for (int i = 0; i < (it->second).size(); i++) { m->mothurOut(it->first + '\t' + it->second[i] + '\n'); }
            }
            m->mothurOutEndLine(); m->setControl_pressed(true);
        }

        bool hasConsTaxonomy = false;
        try {
    
            checkGroups(file, optionalDatasets); if (m->getControl_pressed()) { return; }
    
            if (otuTaxonomies.size() == otuNames.size()) { hasConsTaxonomy = true; }
    
        }catch(H5::Exception& e){ //do nothing taxonomy info does not exist
            m->mothurOut("\nIgnore HDF5 errors, mothur was checking for OTU taxonomies and your file does not contain them. They are not required, continuing.\n");
            hasConsTaxonomy = false;
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

#endif
        
    }
    catch(exception& e) {
        m->errorOut(e, "BiomHDF5", "read");
        exit(1);
    }
}
#ifdef USE_HDF5
//**********************************************************************************************************************
//Process attribute of group or dataset
void BiomHDF5::processAttributes(H5::Group& groupID, set<string>& requiredAttributes) {
    try {
        
        for (set<string>::iterator it = requiredAttributes.begin(); it != requiredAttributes.end(); it++) {
            
            
            H5::Attribute attribute(groupID.openAttribute(*it));
            H5std_string attributeName; attribute.getName(attributeName);
            H5::DataType  attributeType(attribute.getDataType());
            H5::DataSpace attDataSpace = attribute.getSpace();
            
            int rank = attDataSpace.getSimpleExtentNdims(); //number of dimensions
            hsize_t dims[rank]; attDataSpace.getSimpleExtentDims(dims); //size of each dimension
            
            // Read the Attribute Data. Depends on the kind of data
            if (attributeType.getClass() == H5T_STRING) {
                if (rank == 0) {
                    H5std_string biomTableId;
                    attribute.read(attributeType, biomTableId);
                    if (m->getDebug()) { m->mothurOut("[DEBUG]: " + attributeName + " = " + biomTableId + "\n");  }
                }
            }else if (attributeType.getClass() == H5T_INTEGER) {
                matrixElementType == "int";
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
                    }
                }
                
            }else if (attributeType.getClass() == H5T_FLOAT) {
                matrixElementType == "float";
                //m->mothurOut("[WARNING]: the shared file only uses integers, any float values will be rounded down to the nearest integer.\n");
                if (rank == 0) {
                    hsize_t data = 0;
                    attribute.read(attributeType, &data);
                    if (m->getDebug()) { m->mothurOut("[DEBUG]: " + attributeName + " = " + toString(data) + "\n");  }
                }else if (rank == 1) {
                    hsize_t data[dims[0]];
                    attribute.read(attributeType, data);
                    if (m->getDebug()) { m->mothurOut("[DEBUG]: " + attributeName + " = "); }
                    for (int i = 0; i < dims[0]; i++) {
                        if (m->getDebug()) { m->mothurOut(toString(data[i]) + "\t"); }
                    } if (m->getDebug()) { m->mothurOutEndLine(); }
                }
   
            }else { m->mothurOut("[ERROR]: Unexpected datatype class, quitting.\n"); m->setControl_pressed(true);  }
            
            attribute.close();
        }
    }
    catch(exception& e) {
        m->errorOut(e, "BiomHDF5", "processAttributes");
        exit(1);
    }
}
//**********************************************************************************************************************
//check to make sure required groups are present
void BiomHDF5::checkGroups( H5::H5File& file, map<string, vector<string> >& requiredGroups) {
    try {
       
        for (map<string, vector<string> >::iterator it = requiredGroups.begin(); it != requiredGroups.end(); it++) {
            
            H5std_string groupName = it->first;
            vector<string> datasetNames = it->second;
            H5::Group group(file.openGroup(groupName));
            
            for (int h = 0; h < datasetNames.size(); h++) {
                string datasetName = datasetNames[h];
                int numObjects = group.getNumObjs();
                
                if (numObjects != 0) {
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
                                
                                if (groupName == "observation/") {
                                    otuNames.push_back(data[i]);
                                }else if (groupName == "sample/") {
                                    sampleNames.push_back(data[i]);
                                }else if (groupName == "observation/metadata") {
                                    if (datasetName == "taxonomy") { taxonomy.push_back(data[i]); }
                                }
                                delete[] data[i];
                            }  if (m->getDebug()) { m->mothurOutEndLine(); }
                            delete[] data;
                        }else if (rank == 2) {
                            
                            char **data = new char*[dims[0]*dims[1]];

                            H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
                            dataset.read((void*)data, str_type);
                            
                            string otuTaxonomy = ""; int count = 0;
                            for (int i = 0; i < dims[0]*dims[1]; i++) {
                                
                                if (groupName == "observation/metadata") {
                                    if (datasetName == "taxonomy") {  otuTaxonomy += data[i];  otuTaxonomy += ";";  count++; }
                                }
                                
                                if (count == dims[1]) {
                                    if (m->getDebug()) { m->mothurOut("[DEBUG]: " + toString(otuTaxonomy) + "\n"); }
                                    otuTaxonomies.push_back(otuTaxonomy);
                                    if (count > maxLevel) { maxLevel = count; }
                                    count = 0; otuTaxonomy = "";
                                }
                            }
                        }
                    }else if (dataset.getTypeClass() == H5T_INTEGER) {
                        
                        if (rank == 1) {
                            int* data = new int[dims[0]];
                            H5::DataSpace data_mspace(rank, dims);
                            dataset.read(data, H5::PredType::NATIVE_INT, data_mspace, dataSpace);
                            
                            if (m->getDebug()) { m->mothurOut("[DEBUG]: " + datasetName + " = "); }
                            for (int i = 0; i < dims[0]; i++) {
                                if (m->getDebug()) { m->mothurOut(toString(data[i]) + "\t"); }
                                
                                if (groupName == "observation/matrix/") {
                                    if (datasetName == "data") { otudata.push_back(data[i]); }
                                    else if (datasetName == "indices") { indices.push_back(data[i]); }
                                    else if (datasetName == "indptr") { indptr.push_back(data[i]); }
                                }
                            } if (m->getDebug()) { m->mothurOutEndLine(); }
                            delete[] data;
                        }
                        
                    }else if (dataset.getTypeClass() == H5T_FLOAT) {
                        
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
                    }else { m->mothurOut("[ERROR]: Unexpected datatype class, quitting.\n"); m->setControl_pressed(true);  }
                    
                    dataset.close();
                }
            }
            group.close();
        }
    }
    catch(exception& e) {
        m->errorOut(e, "BiomHDF5", "checkGroups");
        exit(1);
    }
}
#endif
/**************************************************************************************************/
