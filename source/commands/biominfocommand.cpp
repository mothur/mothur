//
//  biominfocommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 8/5/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#include "biominfocommand.h"



//**********************************************************************************************************************
vector<string> BiomInfoCommand::setParameters(){
    try {
        CommandParameter pbiom("biom", "InputTypes", "", "", "", "", "","",false,true, true); parameters.push_back(pbiom);
        CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
        CommandParameter prelabund("relabund", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(prelabund);
        CommandParameter pbasis("basis", "Multiple", "otu-sequence", "otu", "", "", "","",false,false); parameters.push_back(pbasis);
        CommandParameter pformat("format", "Multiple", "hdf5-simple", "hdf5", "", "", "","",false,false, true); parameters.push_back(pformat);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter poutput("output", "Multiple", "simple-detail", "detail", "", "", "","",false,false, true); parameters.push_back(poutput);
        CommandParameter pprintlevel("printlevel", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pprintlevel);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false; maxLevel = 0;
        
        vector<string> tempOutNames;
        outputTypes["taxonomy"] = tempOutNames;
        outputTypes["shared"] = tempOutNames;
        outputTypes["constaxonomy"] = tempOutNames;
        outputTypes["taxsummary"] = tempOutNames;
        
        vector<string> myArray;
        for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
        return myArray;
    }
    catch(exception& e) {
        m->errorOut(e, "BiomInfoCommand", "setParameters");
        exit(1);
    }
}
//**********************************************************************************************************************
string BiomInfoCommand::getHelpString(){
    try {
        string helpString = "";
        helpString += "The biom.info command reads a biom file creates a shared file. If your biom file contains metadata mothur will also create taxonomy or constaxonomy along with tax.summary files.\n";
        helpString += "The biom.info command parameters are " + getCommandParameters() + ". The biom parameter is required.\n";
        helpString += "The format parameter allows you indicate type of biom file you have. Options hdf5 or classic. Default is hdf5.\n";
        helpString += "The label parameter allows you to enter a distance label to be used in the shared file created from your biom file.\n";
        helpString += "The relabund parameter allows you to indicate you want the tax.summary file values to be relative abundances rather than raw abundances. Default=F. \n";
        helpString += "The basis parameter allows you indicate what you want the summary file to represent, options are otu and sequence. Default is otu.\n";
        helpString += "The output parameter allows you to specify format of your summary file. Options are simple and detail. The default is detail.\n";
        helpString += "The printlevel parameter allows you to specify taxlevel of your summary file to print to. Options are 1 to the maz level in the file.  The default is -1, meaning max level.  If you select a level greater than the level your sequences classify to, mothur will print to the level your max level. \n";
        helpString += "For example consider the following basis=sequence could give Clostridiales 3 105, where 105 is the total number of sequences whose OTU classified to Clostridiales. ";
        helpString += "Now for basis=otu could give Clostridiales 3 7, where 7 is the number of OTUs that classified to Clostridiales.\n";
        helpString += "The biom.info command should be in the following format: biom.info(biom=test.biom, label=0.03).\n";
        
        getCommonQuestions();
        
        return helpString;
    }
    catch(exception& e) {
        m->errorOut(e, "BiomInfoCommand", "getHelpString");
        exit(1);
    }
}
//**********************************************************************************************************************
string BiomInfoCommand::getCommonQuestions(){
    try {
        vector<string> questions, issues, qanswers, ianswers, howtos, hanswers;
        
        string issue = "Cannot convert error. What do I do?"; issues.push_back(issue);
        string ianswer = "\tThis issue is caused by a matrix_element_type mismatch. The biom file contains a field called 'matrix_element_type'. This field tells mothur what form your observation data is in: int or float. Mothur expects 'int' (an interger value) because the shared file contains interger value abundance counts. If your file contains float values mothur will round down to the nearest integer value. But if your matrix_element_type=int and yet the file contains integer counts in float form, (ie. 31.0 instead of 31) you will get this error. You can resolve this issue by setting matrix_element_type=float in the biom file.\n"; ianswers.push_back(ianswer);
        
        issue = "Mothur can't read my biom file. What does this mean?"; issues.push_back(issue);
        ianswer = "\tThis likely caused by failure to set format=hdf5. Mothur allows for 2 formats: classic (http://biom-format.org/documentation/format_versions/biom-1.0.html) and hdf5 (http://biom-format.org/documentation/format_versions/biom-2.0.html). By default mothur assumes your files are in classic form. If your file is in hdf5 format, then set format=hdf5. NOTE: you can only process hdf5 files if you are using our pre-built version or have built your version of mothur with USEHDF5=yes.\n"; ianswers.push_back(ianswer);
        
        string commonQuestions = util.getFormattedHelp(questions, qanswers, issues, ianswers, howtos, hanswers);
        
        return commonQuestions;
    }
    catch(exception& e) {
        m->errorOut(e, "BiomInfoCommand", "getCommonQuestions");
        exit(1);
    }
}

//**********************************************************************************************************************
string BiomInfoCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "shared") {  pattern = "[filename],[tag],shared"; }
        else if (type == "constaxonomy") {  pattern = "[filename],[tag],cons.taxonomy"; }
        else if (type == "taxonomy") {  pattern = "[filename],[tag],taxonomy"; }
        else if (type == "taxsummary") {  pattern = "[filename],[tag],[tag2],tax.summary"; } //tag2 = "" for taxonomy tag2 = cons for constaxonomy
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "BiomInfoCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
BiomInfoCommand::BiomInfoCommand(string option)  {
    try {
        maxLevel = 0;
        
        //allow user to run help
        if(option == "help") { help(); abort = true; calledHelp = true; }
        else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
        
        else {
            OptionParser parser(option, setParameters());
            map<string, string> parameters = parser.getParameters();
            
            ValidParameters validParameter;
            //check for required parameters
            biomfile = validParameter.validFile(parameters, "biom");
            if (biomfile == "not open") { biomfile = ""; abort = true; }
            else if (biomfile == "not found") { biomfile = ""; m->mothurOut("[ERROR]: You must provide a biom file, please correct.\n");  abort = true;}
            else { current->setBiomFile(biomfile); }
            
            label = validParameter.valid(parameters, "label");
            if (label == "not found") { label = "userLabel"; }
            
            output = validParameter.valid(parameters, "output");		if(output == "not found"){	output = "detail"; }
            if ((output != "simple") && (output != "detail")) { m->mothurOut(output + " is not a valid output form. Options are simple and detail. I will use detail.\n");  output = "detail"; }
            
            string temp = validParameter.valid(parameters, "relabund");
            if (temp == "not found"){	temp = "false";			}
            else { temp = util.getSimpleName(temp); }
            relabund = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "printlevel");		if (temp == "not found"){	temp = "-1";		}
            util.mothurConvert(temp, printlevel);
            
            basis = validParameter.valid(parameters, "basis");
            if (basis == "not found") { basis = "otu"; }
            
            if ((basis != "otu") && (basis != "sequence")) { m->mothurOut("Invalid option for basis. basis options are otu and sequence, using otu.\n"); }
            
            format = validParameter.valid(parameters, "format");
            if (format == "not found") { format = "classic"; }
            
            if ((format != "hdf5") && (format != "classic")) { m->mothurOut("Invalid option for format. format options are hdf5 and classic, using hdf5.\n"); }
            
            if (format == "hdf5") {
                #ifdef USE_HDF5
                //do nothing we have the api
                #else
                    m->mothurOut("[ERROR]: To read HDF5 biom files, you must have the API installed, quitting.\n"); abort=true;
                #endif
            }
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "BiomInfoCommand", "BiomInfoCommand");
        exit(1);
    }
}
//**********************************************************************************************************************

int BiomInfoCommand::execute(){
    try {
        
        if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        long start = time(NULL);
        
        if (format == "hdf5")   { extractFilesFromHDF5();   }
        else                    { createFilesFromBiomSimple();    }
        
        m->mothurOut("\nIt took " + toString(time(NULL) - start) + " create mothur files from your biom file.\n\n");
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); } }
        
        string currentName = "";
        itTypes = outputTypes.find("shared");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setSharedFile(currentName); }
        }
        
        //set taxonomy file as new current taxonomyfile
        itTypes = outputTypes.find("taxonomy");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setTaxonomyFile(currentName); }
        }
        
        //set constaxonomy file as new current constaxonomyfile
        itTypes = outputTypes.find("constaxonomy");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setConsTaxonomyFile(currentName); }
        }
        
        m->mothurOutEndLine();
        m->mothurOut("Output File Names: \n"); 
        for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
        m->mothurOutEndLine();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "BiomInfoCommand", "execute");
        exit(1);
    }
}
#ifdef USE_HDF5
//**********************************************************************************************************************
//Process attribute of group or dataset
void BiomInfoCommand::processAttributes(H5::Group& groupID, set<string>& requiredAttributes) {
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
                m->mothurOut("[WARNING]: the shared file only uses integers, any float values will be rounded down to the nearest integer.\n");
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
        m->errorOut(e, "BiomInfoCommand", "processAttributes");
        exit(1);
    }
}
//**********************************************************************************************************************
//check to make sure required groups are present
void BiomInfoCommand::checkGroups( H5::H5File& file, map<string, vector<string> >& requiredGroups) {
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
                                    taxonomy.push_back(otuTaxonomy);
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
                                    if (datasetName == "data") { otudata.push_back((int)data[i]); }
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
        m->errorOut(e, "BiomInfoCommand", "checkGroups");
        exit(1);
    }
}
#endif
//**********************************************************************************************************************
int BiomInfoCommand::extractFilesFromHDF5() {
    try {
        //getting output filename
        string filename = biomfile;
        if (outputdir == "") { outputdir += util.hasPath(filename); }
        
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(filename));
        variables["[tag]"] = label;
        string sharedFilename = getOutputFileName("shared",variables);
        string taxFilename = getOutputFileName("constaxonomy",variables);
        variables["[tag2]"] = "cons";
        string taxSumFilename = getOutputFileName("taxsummary",variables);
        
        nnz = 0;
        set<string> requiredTopLevelAttrib;
        map<string, vector<string> > requiredOTUDatasets;
        map<string, vector<string> > requiredSampleDatasets;
        map<string, vector<string> > optionalDatasets;

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
        H5::H5File file( filename.c_str(), H5F_ACC_RDONLY );
        H5::Group     what(file.openGroup( "/" ));
      
        processAttributes(what, requiredTopLevelAttrib); if (m->getControl_pressed()) { return 0; }

        try {
            
            checkGroups(file, requiredOTUDatasets); if (m->getControl_pressed()) { return 0; }
            
        }catch(H5::Exception& e){ //do nothing taxonomy info does not exist
            m->mothurOut("[ERROR]: Missing required groups or datasets, aborting. Required datasets include: \n");
            for (map<string, vector< string> >::iterator it = requiredOTUDatasets.begin(); it != requiredOTUDatasets.end(); it++) {
                for (int i = 0; i < (it->second).size(); i++) { m->mothurOut(it->first + '\t' + it->second[i] + '\n'); }
            }
            m->mothurOutEndLine(); m->setControl_pressed(true);
        }
        
        try {
            
            checkGroups(file, requiredSampleDatasets); if (m->getControl_pressed()) { return 0; }
            
        }catch(H5::Exception& e){ //do nothing taxonomy info does not exist
            m->mothurOut("[ERROR]: Missing required groups or datasets, aborting. Required datasets include: \n");
            for (map<string, vector< string> >::iterator it = requiredOTUDatasets.begin(); it != requiredOTUDatasets.end(); it++) {
                for (int i = 0; i < (it->second).size(); i++) { m->mothurOut(it->first + '\t' + it->second[i] + '\n'); }
            }
            m->mothurOutEndLine(); m->setControl_pressed(true);
        }
        
        bool hasTaxonomy = false;
        try {
            
            checkGroups(file, optionalDatasets); if (m->getControl_pressed()) { return 0; }
            
            if (taxonomy.size() == otuNames.size()) { hasTaxonomy = true; }
            
        }catch(H5::Exception& e){ //do nothing taxonomy info does not exist
            m->mothurOut("\nIgnore HDF5 errors, mothur was checking for OTU taxonomies and your file does not contain them. They are not required, continuing.\n");
            hasTaxonomy = false;
        }

        
        bool error = false;
        if (nnz != otudata.size()) { error = true; }
        
        //create shared file
        sort(sampleNames.begin(), sampleNames.end());
        
        //create empty sharedrabundvectors so we can add otus below
        SharedRAbundVectors lookup;
        for (int i = 0; i < sampleNames.size(); i++) {
            SharedRAbundVector* temp = new SharedRAbundVector();
            temp->setGroup(sampleNames[i]);
            lookup.push_back(temp);
        }
        
        lookup.setLabels(label);
        
        ofstream outTax;
        if (hasTaxonomy) {
            outputNames.push_back(taxFilename); outputTypes["constaxonomy"].push_back(taxFilename);
            util.openOutputFile(taxFilename, outTax);
            outTax << "OTU\tSize\tTaxonomy\n";
        }
        
        //for each otu
        int count = 0;
        for (int h = 0; h < indptr.size()-1; h++) {
            int otuStart = indptr[h];
            int otuEnd = indptr[h+1];
            
            vector<int> otuAbunds; otuAbunds.resize(sampleNames.size(), 0); //initialze otus sample abundances to 0 - only non zero abunds are recorded
            
            for (int i = otuStart; i < otuEnd; i++) {
                otuAbunds[indices[i]] = otudata[count]; count++;
            }
            
            lookup.push_back(otuAbunds, otuNames[h]);
        }
        
        ofstream out; util.openOutputFile(sharedFilename, out);
        m->mothurOut("\n"+lookup.getLabel()+"\n"); bool printHeaders = true;
        lookup.print(out, printHeaders);
        out.close();
        outputNames.push_back(sharedFilename); outputTypes["shared"].push_back(sharedFilename);
        
        if (hasTaxonomy) {
            
            CountTable ct; for (int j = 0; j < sampleNames.size(); j++) {  ct.addGroup(sampleNames[j]); }
            int numBins = lookup.getNumBins();
            
            for (int i = 0; i < numBins; i++) {
                vector<int> abunds;
                for (int j = 0; j < lookup.size(); j++) {
                    if (m->getControl_pressed()) { break; }
                    int abund = lookup.get(i, sampleNames[j]);
                    if (basis == "otu") { if (abund > 0) { abund = 1;  } } //count presence in otu
                    abunds.push_back(abund);
                }
                ct.push_back(otuNames[i], abunds);
            }
            
            PhyloSummary taxaSum(&ct, relabund, printlevel);
            
            for (int i = 0; i < lookup.getNumBins(); i++) {
                if (m->getControl_pressed()) { break; }
                
                int total = 0;
                map<string, bool> containsGroup;
                for (int j = 0; j < lookup.size(); j++) {
                    int abund = lookup.get(i, sampleNames[j]);
                    total += abund;
                    containsGroup[sampleNames[j]] = abund;
                }
                
                string newTax = util.addUnclassifieds(taxonomy[i], maxLevel, false);
                outTax << otuNames[i] << '\t' << total << '\t' << newTax << endl;
                
                if (basis == "sequence") {
                    taxaSum.addSeqToTree(otuNames[i], newTax);
                }else {
                    taxaSum.addSeqToTree(newTax, containsGroup); //add otu
                }
            }
            outTax.close();
            
            //write taxonomy summary file
            outputNames.push_back(taxSumFilename); outputTypes["taxsummary"].push_back(taxSumFilename);
            ofstream outTaxSum;
            util.openOutputFile(taxSumFilename, outTaxSum);
            taxaSum.print(outTaxSum, output);
            outTaxSum.close();
        }

        #endif
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "BiomInfoCommand", "extractFilesFromHDF5");
        exit(1);
    }
}
//**********************************************************************************************************************
//biom version 0.9.1
void BiomInfoCommand::createFilesFromBiomSimple() {
    try {
        //pass filename, printlevel
        BiomSimple biom(biomfile, label);
        
        //getting output filename
        string filename = biomfile;
        if (outputdir == "") { outputdir += util.hasPath(filename); }
        fileroot = outputdir + util.getRootName(util.getSimpleName(biomfile));
        
        SharedRAbundVectors* shared = biom.getSharedRAbundVectors();
        
        CountTable ct;
        vector< map<string, bool> > otuContainsGroups;
        
        if (shared != NULL) {
            map<string, string> variables;
            variables["[filename]"] = fileroot;
            variables["[tag]"] = label;
            string sharedFilename = getOutputFileName("shared",variables);
            outputNames.push_back(sharedFilename); outputTypes["shared"].push_back(sharedFilename);
            bool printHeaders = true;
            
            ofstream out; util.openOutputFile(sharedFilename, out);
            shared->print(out, printHeaders);
            out.close();
            
            vector<string> groupNames = shared->getNamesGroups();
            for (int j = 0; j < groupNames.size(); j++) {  ct.addGroup(groupNames[j]); }
                
            int numBins = shared->getNumBins();
            
            for (int i = 0; i < numBins; i++) {
                int total = 0;
                map<string, bool> containsGroup;
                vector<int> abunds;
                for (int j = 0; j < shared->size(); j++) {
                    if (m->getControl_pressed()) { break; }
                    int abund = shared->get(i, groupNames[j]);
                    
                    total += abund;
                    containsGroup[groupNames[j]] = abund;
                    
                    if (basis == "otu") { if (abund > 0) { abund = 1;  } } //count presence in otu
                    abunds.push_back(abund);
                }
                ct.push_back(shared->getOTUName(i), abunds);
                otuContainsGroups.push_back(containsGroup);
            }
        }
        
        //print group taxonomies if given
        map<string, string> groupTaxonomies = biom.getGroupTaxonomies();
        
        if (groupTaxonomies.size() != 0) {
            //write taxonomy file
            map<string, string> variables;
            variables["[filename]"] = fileroot;
            variables["[tag]"] = label;
            string taxFilename = getOutputFileName("taxonomy",variables);
            outputNames.push_back(taxFilename); outputTypes["taxonomy"].push_back(taxFilename);
            ofstream outTax; util.openOutputFile(taxFilename, outTax);
            
            GroupMap* g = NULL; PhyloSummary taxSum(g, relabund, printlevel);
            
            //print group taxonomy if given
            for (map<string, string>::iterator it = groupTaxonomies.begin(); it!= groupTaxonomies.end(); it++) {
                outTax << it->first << '\t' << it->second << endl;
                taxSum.addSeqToTree(it->first, it->second);
            }
            outTax.close();
            
            //write taxonomy file
            variables["[tag2]"] = "";
            string taxSumFilename = getOutputFileName("taxsummary",variables);
            outputNames.push_back(taxSumFilename); outputTypes["taxsummary"].push_back(taxSumFilename);
            ofstream outTaxSum; util.openOutputFile(taxSumFilename, outTaxSum);
                
            //write tax.summary
            if (relabund)   {   taxSum.print(outTaxSum, relabund);     }
            else            {   taxSum.print(outTaxSum, output);       }
                
            outTaxSum.close();
        }
        
        
        //print consTaxonomy if given
        vector<Taxonomy> consTax = biom.getConsTaxonomies();
        if (consTax.size() != 0) {

            //write taxonomy file
            map<string, string> variables;
            variables["[filename]"] = fileroot;
            variables["[tag]"] = label;
            string taxFilename = getOutputFileName("constaxonomy",variables);
            outputNames.push_back(taxFilename); outputTypes["constaxonomy"].push_back(taxFilename);
            ofstream outTax; util.openOutputFile(taxFilename, outTax);
            
            outTax << "OTU\tSize\tTaxonomy\n";
            
            PhyloSummary consTaxSum(&ct, relabund, printlevel);
            
            for (int i = 0; i < consTax.size(); i++) {
                
                consTax[i].printConsTax(outTax);
                
                if (basis == "sequence")    {  consTaxSum.addSeqToTree(consTax[i].getName(), consTax[i].getConsTaxString());  }
                else                       { consTaxSum.addSeqToTree(consTax[i].getConsTaxString(), otuContainsGroups[i]); } //add otu
            }
            outTax.close();
            
            variables["[tag2]"] = "cons";
            string taxSumFilename = getOutputFileName("taxsummary",variables);
            outputNames.push_back(taxSumFilename); outputTypes["taxsummary"].push_back(taxSumFilename);
            ofstream outTaxSum; util.openOutputFile(taxSumFilename, outTaxSum);
        
            //write tax.summary
            if (relabund)   {   consTaxSum.print(outTaxSum, relabund);     }
            else            {   consTaxSum.print(outTaxSum, output);       }
            
            outTaxSum.close();
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "BiomInfoCommand", "createFilesFromBiom");
        exit(1);
    }
}
//**********************************************************************************************************************/


