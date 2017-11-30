//
//  makebiomcommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/16/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "makebiomcommand.h"

#include "inputdata.h"

#include "phylotree.h"

//taken from http://biom-format.org/documentation/biom_format.html
/* Minimal Sparse 
 {
 "id":null,
 "format": "Biological Observation Matrix 0.9.1",
 "format_url": "http://biom-format.org",
 "type": "OTU table",
 "generated_by": "QIIME revision 1.4.0-dev",
 "date": "2011-12-19T19:00:00",
 "rows":[
 {"id":"GG_OTU_1", "metadata":null},
 {"id":"GG_OTU_2", "metadata":null},
 {"id":"GG_OTU_3", "metadata":null},
 {"id":"GG_OTU_4", "metadata":null},
 {"id":"GG_OTU_5", "metadata":null}
 ],
 "columns": [
 {"id":"Sample1", "metadata":null},
 {"id":"Sample2", "metadata":null},
 {"id":"Sample3", "metadata":null},
 {"id":"Sample4", "metadata":null},
 {"id":"Sample5", "metadata":null},
 {"id":"Sample6", "metadata":null}
 ],
 "matrix_type": "sparse",
 "matrix_element_type": "int",
 "shape": [5, 6],
 "data":[[0,2,1],
 [1,0,5],
 [1,1,1],
 [1,3,2],
 [1,4,3],
 [1,5,1],
 [2,2,1],
 [2,3,4],
 [2,4,2],
 [3,0,2],
 [3,1,1],
 [3,2,1],
 [3,5,1],
 [4,1,1],
 [4,2,1]
 ]
 }
 */
/* Minimal dense
 {
 "id":null,
 "format": "Biological Observation Matrix 0.9.1",
 "format_url": "http://biom-format.org",
 "type": "OTU table",
 "generated_by": "QIIME revision 1.4.0-dev",
 "date": "2011-12-19T19:00:00",
 "rows":[
 {"id":"GG_OTU_1", "metadata":null},
 {"id":"GG_OTU_2", "metadata":null},
 {"id":"GG_OTU_3", "metadata":null},
 {"id":"GG_OTU_4", "metadata":null},
 {"id":"GG_OTU_5", "metadata":null}
 ],
 "columns": [
 {"id":"Sample1", "metadata":null},
 {"id":"Sample2", "metadata":null},
 {"id":"Sample3", "metadata":null},
 {"id":"Sample4", "metadata":null},
 {"id":"Sample5", "metadata":null},
 {"id":"Sample6", "metadata":null}
 ],
 "matrix_type": "dense",
 "matrix_element_type": "int",
 "shape": [5,6],
 "data":  [[0,0,1,0,0,0],
 [5,1,0,2,3,1],
 [0,0,1,4,2,0],
 [2,1,1,0,0,1],
 [0,1,1,0,0,0]]
 }
 */
//**********************************************************************************************************************
vector<string> MakeBiomCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "SharedRel", "SharedRel", "none","biom",false,false,true); parameters.push_back(pshared);
        CommandParameter prelabund("relabund", "InputTypes", "", "", "SharedRel", "SharedRel", "none","biom",false,false,true); parameters.push_back(prelabund);
        CommandParameter pcontaxonomy("constaxonomy", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pcontaxonomy);
        CommandParameter preference("reftaxonomy", "InputTypes", "", "", "none", "none", "refPi","",false,false); parameters.push_back(preference);
        CommandParameter pmetadata("metadata", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pmetadata);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
        CommandParameter ppicrust("picrust", "InputTypes", "", "", "none", "none", "refPi","shared",false,false); parameters.push_back(ppicrust);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        CommandParameter pmatrixtype("matrixtype", "Multiple", "sparse-dense", "sparse", "", "", "","",false,false); parameters.push_back(pmatrixtype);

		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeBiomCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string MakeBiomCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The make.biom command parameters are shared, relabund, contaxonomy, metadata, groups, matrixtype, picrust, reftaxonomy and label.  shared or relabund are required, unless you have a valid current file.\n"; //
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like included. The group names are separated by dashes.\n";
		helpString += "The label parameter allows you to select what distance levels you would like, and are also separated by dashes.\n";
		helpString += "The matrixtype parameter allows you to select what type you would like to make. Choices are sparse and dense, default is sparse.\n";
        helpString += "The contaxonomy file is the taxonomy file outputted by classify.otu(list=yourListfile, taxonomy=yourTaxonomyFile). Be SURE that the you are the constaxonomy file distance matches the shared file distance.  ie, for *.0.03.cons.taxonomy set label=0.03. Mothur is smart enough to handle shared files that have been subsampled. It is used to assign taxonomy information to the metadata of rows.\n";
        helpString += "The metadata parameter is used to provide experimental parameters to the columns.  Things like 'sample1 gut human_gut'. \n";
        helpString += "The picrust parameter is used to provide the greengenes OTU IDs map table.  NOTE: Picrust requires a greengenes taxonomy. \n";
        helpString += "The referencetax parameter is used with the picrust parameter.  Picrust requires the greengenes OTU IDs to be in the biom file. \n";
		helpString += "The make.biom command should be in the following format: make.biom(shared=yourShared, groups=yourGroups, label=yourLabels).\n";
		helpString += "Example make.biom(shared=abrecovery.an.shared, groups=A-B-C).\n";
		helpString += "The default value for groups is all the groups in your groupfile, and all labels in your inputfile will be used.\n";
		helpString += "The make.biom command outputs a .biom file.\n";
		helpString += "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeBiomCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string MakeBiomCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "biom") {  pattern = "[filename],[distance],biom"; }
        else if (type == "shared") {  pattern = "[filename],[distance],biom_shared"; }
        else if (type == "relabund") {  pattern = "[filename],[distance],biom_relabund"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeBiomCommand", "getOutputPattern");
        exit(1);
    }
}

//**********************************************************************************************************************
MakeBiomCommand::MakeBiomCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["biom"] = tempOutNames;
        outputTypes["shared"] = tempOutNames;
        outputTypes["relabund"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeBiomCommand", "MakeBiomCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
MakeBiomCommand::MakeBiomCommand(string option) {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
        
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			map<string,string>::iterator it;
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}

			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["biom"] = tempOutNames;
            outputTypes["shared"] = tempOutNames;
            outputTypes["relabund"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
                
                it = parameters.find("relabund");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["relabund"] = inputDir + it->second;		}
                }
                
                it = parameters.find("constaxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["constaxonomy"] = inputDir + it->second;		}
				}
                
                it = parameters.find("reftaxonomy");
				//user has given a template file
				if(it != parameters.end()){
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["reftaxonomy"] = inputDir + it->second;		}
				}
                
                it = parameters.find("picrust");
				//user has given a template file
				if(it != parameters.end()){
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["picrust"] = inputDir + it->second;		}
				}
                
                it = parameters.find("metadata");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["metadata"] = inputDir + it->second;		}
				}
			}
            
            relabundfile = validParameter.validFile(parameters, "relabund");
            if (relabundfile == "not open") { abort = true; }
            else if (relabundfile == "not found") { relabundfile = ""; }
            else { inputFileName = relabundfile; fileFormat = "relabund"; current->setRelAbundFile(relabundfile); }
            
            sharedfile = validParameter.validFile(parameters, "shared");
            if (sharedfile == "not open") { abort = true; }
            else if (sharedfile == "not found") { sharedfile = ""; }
            else { inputFileName = sharedfile; fileFormat = "sharedfile"; current->setSharedFile(sharedfile); }
            
            
            if ((relabundfile == "") && (sharedfile == "")) {
                //is there are current file available for either of these?
                //give priority to shared, then relabund
                sharedfile = current->getSharedFile();
                if (sharedfile != "") {  inputFileName = sharedfile; fileFormat="sharedfile"; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
                else {
                    relabundfile = current->getRelAbundFile();
                    if (relabundfile != "") {  inputFileName = relabundfile; fileFormat="relabund"; m->mothurOut("Using " + relabundfile + " as input file for the relabund parameter."); m->mothurOutEndLine(); }
                    else {
                        m->mothurOut("No valid current files. You must provide a shared or relabund."); m->mothurOutEndLine(); abort = true;
                    }
                }
            }
            
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = util.hasPath(inputFileName);		}
            
            contaxonomyfile = validParameter.validFile(parameters, "constaxonomy");
			if (contaxonomyfile == "not found") {  contaxonomyfile = "";  }
			else if (contaxonomyfile == "not open") { contaxonomyfile = ""; abort = true; }
            
            referenceTax = validParameter.validFile(parameters, "reftaxonomy");
			if (referenceTax == "not found") {  referenceTax = "";  }
			else if (referenceTax == "not open") { referenceTax = ""; abort = true; }
            
            picrustOtuFile = validParameter.validFile(parameters, "picrust");
			if (picrustOtuFile == "not found") {  picrustOtuFile = "";  }
			else if (picrustOtuFile == "not open") { picrustOtuFile = ""; abort = true; }

            metadatafile = validParameter.validFile(parameters, "metadata");
			if (metadatafile == "not found") {  metadatafile = "";  }
			else if (metadatafile == "not open") { metadatafile = ""; abort = true; }
            
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = ""; }
			else { 
				util.splitAtDash(groups, Groups);
                if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
			}
			
            if (picrustOtuFile != "") {
                picrust=true;
                if (contaxonomyfile == "") {  m->mothurOut("[ERROR]: the constaxonomy parameter is required with the picrust parameter, aborting."); m->mothurOutEndLine(); abort = true;  }
                if (referenceTax == "") {  m->mothurOut("[ERROR]: the reftaxonomy parameter is required with the picrust parameter, aborting."); m->mothurOutEndLine(); abort = true;  }
            }else { picrust=false; }
            
            if ((contaxonomyfile != "") && (labels.size() > 1)) { m->mothurOut("[ERROR]: the contaxonomy parameter cannot be used with multiple labels."); m->mothurOutEndLine(); abort = true; }
            
			format = validParameter.valid(parameters, "matrixtype");				if (format == "not found") { format = "sparse"; }
			
			if ((format != "sparse") && (format != "dense")) {
				m->mothurOut(format + " is not a valid option for the matrixtype parameter. Options are sparse and dense."); m->mothurOutEndLine(); abort = true; 
			}
		}
        
	}
	catch(exception& e) {
		m->errorOut(e, "MakeBiomCommand", "MakeBiomCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int MakeBiomCommand::execute(){
	try {
        
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        SharedRAbundVectors* lookup = NULL;
        SharedRAbundFloatVectors* lookupRel = NULL;
        string lastLabel;
        
		InputData input(inputFileName, fileFormat, Groups);
        if (fileFormat == "sharedfile") {
            lookup = input.getSharedRAbundVectors();
            Groups = lookup->getNamesGroups();
            lastLabel = lookup->getLabel();
            getSampleMetaData(lookup);
        }else                        {
            lookupRel = input.getSharedRAbundFloatVectors();
            Groups = lookupRel->getNamesGroups();
            lastLabel = lookupRel->getLabel();
            getSampleMetaData(lookupRel);
        }
        
        //if user did not specify a label, then use first one
        if ((contaxonomyfile != "") && (labels.size() == 0)) {
            allLines = 0;
            labels.insert(lastLabel);
        }
        
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
        
        if (fileFormat == "sharedfile") {
            
            //as long as you are not at the end of the file or done wih the lines you want
            while((lookup != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
                
                if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } delete lookup;  return 0; }
                
                if(allLines == 1 || labels.count(lookup->getLabel()) == 1){
                    
                    m->mothurOut(lookup->getLabel()); m->mothurOutEndLine();
                    getBiom(lookup);
                    
                    processedLabels.insert(lookup->getLabel());
                    userLabels.erase(lookup->getLabel());
                }
                
                if ((util.anyLabelsToProcess(lookup->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
                    string saveLabel = lookup->getLabel();
                    
                    delete lookup;
                    lookup = input.getSharedRAbundVectors(lastLabel);
                    m->mothurOut(lookup->getLabel()); m->mothurOutEndLine();
                    
                    getBiom(lookup);
                    
                    processedLabels.insert(lookup->getLabel());
                    userLabels.erase(lookup->getLabel());
                    
                    //restore real lastlabel to save below
                    lookup->setLabels(saveLabel);
                }
                
                lastLabel = lookup->getLabel();
                
                //prevent memory leak and get next set
                delete lookup;
                lookup = input.getSharedRAbundVectors();				
            }
        }else {
            
            //as long as you are not at the end of the file or done wih the lines you want
            while((lookupRel != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
               
                if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } delete lookupRel; return 0; }
                
                if(allLines == 1 || labels.count(lookupRel->getLabel()) == 1){
                    
                    m->mothurOut(lookupRel->getLabel()); m->mothurOutEndLine();
                    getBiom(lookupRel);
                    
                    processedLabels.insert(lookupRel->getLabel());
                    userLabels.erase(lookupRel->getLabel());
                }
                
                if ((util.anyLabelsToProcess(lookupRel->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
                    string saveLabel = lookupRel->getLabel();
                    
                    delete lookupRel;
                    lookupRel = input.getSharedRAbundFloatVectors(lastLabel);
                    m->mothurOut(lookupRel->getLabel()); m->mothurOutEndLine();
                    
                    getBiom(lookupRel);
                    
                    processedLabels.insert(lookupRel->getLabel());
                    userLabels.erase(lookupRel->getLabel());
                    
                    //restore real lastlabel to save below
                    lookupRel->setLabels(saveLabel);
                }
                
                lastLabel = lookupRel->getLabel();
                
                //prevent memory leak and get next set
                delete lookupRel;
                lookupRel = input.getSharedRAbundFloatVectors();
            }
        }
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); }  return 0; }     
        
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			m->mothurOut("Your file does not include the label " + *it); 
			if (processedLabels.count(lastLabel) != 1) {
				m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
				needToRun = true;
			}else {
				m->mothurOut(". Please refer to " + lastLabel + "."); m->mothurOutEndLine();
			}
		}
        
		//run last label if you need to
		if (needToRun )  {
            if (fileFormat == "sharedfile") {
                delete lookup;
                lookup = input.getSharedRAbundVectors(lastLabel);
                
                m->mothurOut(lookup->getLabel()); m->mothurOutEndLine();
                getBiom(lookup);
                
                delete lookup;
            }else {
                delete lookupRel;
                lookupRel = input.getSharedRAbundFloatVectors(lastLabel);
                
                m->mothurOut(lookupRel->getLabel()); m->mothurOutEndLine();
                getBiom(lookupRel);
                
                delete lookupRel;
            }
		}
		
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); }  return 0; }     
		
        //set sabund file as new current sabundfile
        string currentName = "";
		itTypes = outputTypes.find("biom");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setBiomFile(currentName); }
		}

        
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeBiomCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int MakeBiomCommand::getBiom(SharedRAbundVectors*& lookup){
	try {
        map<string, string> variables; 
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(sharedfile));
        variables["[distance]"] = lookup->getLabel();
        string outputFileName = getOutputFileName("biom",variables);
		ofstream out;
		util.openOutputFile(outputFileName, out);
		outputNames.push_back(outputFileName); outputTypes["biom"].push_back(outputFileName);

        string mothurString = "mothur" + toString(current->getVersion());
        time_t rawtime;
        struct tm * timeinfo;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        string dateString = asctime (timeinfo);
        int pos = dateString.find('\n');
        if (pos != string::npos) { dateString = dateString.substr(0, pos);}
        string spaces = "      ";
        
        //standard 
        out << "{\n" + spaces + "\"id\":\"" + util.getSimpleName(sharedfile) + "-" + lookup->getLabel() + "\",\n" + spaces + "\"format\": \"Biological Observation Matrix 0.9.1\",\n" + spaces + "\"format_url\": \"http://biom-format.org\",\n";
        out << spaces + "\"type\": \"OTU table\",\n" + spaces + "\"generated_by\": \"" << mothurString << "\",\n" + spaces + "\"date\": \"" << dateString << "\",\n";
        
        
        vector<string> metadata = getMetaData(lookup);
        int numBins = lookup->getNumBins();
        
        if (m->getControl_pressed()) {  out.close(); return 0; }
        
        //get row info
        /*"rows":[
                {"id":"GG_OTU_1", "metadata":null},
                {"id":"GG_OTU_2", "metadata":null},
                {"id":"GG_OTU_3", "metadata":null},
                {"id":"GG_OTU_4", "metadata":null},
                {"id":"GG_OTU_5", "metadata":null}
                ],*/
        out << spaces + "\"rows\":[\n";
        string rowFront = spaces + spaces + "{\"id\":\"";
        string rowBack = "\", \"metadata\":";
        vector<string> currentLabels = lookup->getOTUNames();
        for (int i = 0; i < numBins-1; i++) {
            if (m->getControl_pressed()) {  out.close(); return 0; }
            out << rowFront << currentLabels[i] << rowBack << metadata[i] << "},\n";
        }
        out << rowFront << currentLabels[(numBins-1)] << rowBack << metadata[(numBins-1)] << "}\n" + spaces + "],\n";
       
        //get column info
        /*"columns": [
                    {"id":"Sample1", "metadata":null},
                    {"id":"Sample2", "metadata":null},
                    {"id":"Sample3", "metadata":null},
                    {"id":"Sample4", "metadata":null},
                    {"id":"Sample5", "metadata":null},
                    {"id":"Sample6", "metadata":null}
                    ],*/
        
        string colBack = "\", \"metadata\":";
        out << spaces + "\"columns\":[\n";
        vector<string> namesOfGroups = lookup->getNamesGroups();
        for (int i = 0; i < namesOfGroups.size()-1; i++) {
            if (m->getControl_pressed()) {  out.close(); return 0; }
            out << rowFront << namesOfGroups[i] << colBack << sampleMetadata[i] << "},\n";
        }
        out << rowFront << namesOfGroups[(namesOfGroups.size()-1)] << colBack << sampleMetadata[lookup->size()-1] << "}\n" + spaces + "],\n";
        
        out << spaces + "\"matrix_type\": \"" << format << "\",\n" + spaces + "\"matrix_element_type\": \"int\",\n";
        out <<  spaces + "\"shape\": [" << numBins << "," << lookup->size() << "],\n";
        out << spaces + "\"data\":  [";
        
        vector<string> dataRows;
        if (format == "sparse") {
            /*"data":[[0,2,1],
             [1,0,5],
             [1,1,1],
             [1,3,2],
             [1,4,3],
             [1,5,1],
             [2,2,1],
             [2,3,4],
             [2,4,2],
             [3,0,2],
             [3,1,1],
             [3,2,1],
             [3,5,1],
             [4,1,1],
             [4,2,1]
             ]*/
            string output = "";
            for (int i = 0; i < lookup->getNumBins(); i++) {
                
                if (m->getControl_pressed()) { out.close(); return 0; }
                vector<int> binAbunds = lookup->getOTU(i);
                
                for (int j = 0; j < binAbunds.size(); j++) {
                    int abund = binAbunds[j];
                    string binInfo = "[" + toString(i) + "," + toString(j) + "," + toString(abund) + "]";
                    //only print non zero values
                    if (abund != 0) { dataRows.push_back(binInfo); }
                }
            }
        }else {
            
            /* "matrix_type": "dense",
             "matrix_element_type": "int",
             "shape": [5,6],
             "data":  [[0,0,1,0,0,0],
             [5,1,0,2,3,1],
             [0,0,1,4,2,0],
             [2,1,1,0,0,1],
             [0,1,1,0,0,0]]*/
            
            for (int i = 0; i < lookup->getNumBins(); i++) {
                
                if (m->getControl_pressed()) { out.close(); return 0; }
                
                string binInfo = "[";
                vector<int> binAbund = lookup->getOTU(i);
                for (int j = 0; j < binAbund.size()-1; j++) {  binInfo += toString(binAbund[j]) + ","; }
                binInfo += toString(binAbund[binAbund.size()-1]) + "]";
                dataRows.push_back(binInfo);
            }
        }
        
        for (int i = 0; i < dataRows.size()-1; i++) {
            out << dataRows[i] << ",\n" + spaces  + spaces;
        }
        out << dataRows[dataRows.size()-1] << "]\n";
        
        out << "}\n";
        out.close();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "MakeBiomCommand", "getBiom");
		exit(1);
	}
}
//**********************************************************************************************************************
int MakeBiomCommand::getBiom(SharedRAbundFloatVectors*& lookup){
    try {
        map<string, string> variables;
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(inputFileName));
        variables["[distance]"] = lookup->getLabel();
        string outputFileName = getOutputFileName("biom",variables);
        ofstream out;
        util.openOutputFile(outputFileName, out);
        outputNames.push_back(outputFileName); outputTypes["biom"].push_back(outputFileName);
        
        string mothurString = "mothur" + toString(current->getVersion());
        time_t rawtime;
        struct tm * timeinfo;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        string dateString = asctime (timeinfo);
        int pos = dateString.find('\n');
        if (pos != string::npos) { dateString = dateString.substr(0, pos);}
        string spaces = "      ";
        
        //standard
        out << "{\n" + spaces + "\"id\":\"" + sharedfile + "-" + lookup->getLabel() + "\",\n" + spaces + "\"format\": \"Biological Observation Matrix 0.9.1\",\n" + spaces + "\"format_url\": \"http://biom-format.org\",\n";
        out << spaces + "\"type\": \"OTU table\",\n" + spaces + "\"generated_by\": \"" << mothurString << "\",\n" + spaces + "\"date\": \"" << dateString << "\",\n";
        
        
        vector<string> metadata = getMetaData(lookup);
        int numBins = lookup->getNumBins();
        
        if (m->getControl_pressed()) {  out.close(); return 0; }
        
        //get row info
        /*"rows":[
         {"id":"GG_OTU_1", "metadata":null},
         {"id":"GG_OTU_2", "metadata":null},
         {"id":"GG_OTU_3", "metadata":null},
         {"id":"GG_OTU_4", "metadata":null},
         {"id":"GG_OTU_5", "metadata":null}
         ],*/
        out << spaces + "\"rows\":[\n";
        string rowFront = spaces + spaces + "{\"id\":\"";
        string rowBack = "\", \"metadata\":";
        vector<string> currentLabels = lookup->getOTUNames();
        for (int i = 0; i < numBins-1; i++) {
            if (m->getControl_pressed()) {  out.close(); return 0; }
            out << rowFront << currentLabels[i] << rowBack << metadata[i] << "},\n";
        }
        out << rowFront << currentLabels[(numBins-1)] << rowBack << metadata[(numBins-1)] << "}\n" + spaces + "],\n";
        
        //get column info
        /*"columns": [
         {"id":"Sample1", "metadata":null},
         {"id":"Sample2", "metadata":null},
         {"id":"Sample3", "metadata":null},
         {"id":"Sample4", "metadata":null},
         {"id":"Sample5", "metadata":null},
         {"id":"Sample6", "metadata":null}
         ],*/
        
        string colBack = "\", \"metadata\":";
        out << spaces + "\"columns\":[\n";
        vector<string> namesOfGroups = lookup->getNamesGroups();
        for (int i = 0; i < namesOfGroups.size()-1; i++) {
            if (m->getControl_pressed()) {  out.close(); return 0; }
            out << rowFront << namesOfGroups[i] << colBack << sampleMetadata[i] << "},\n";
        }
        out << rowFront << namesOfGroups[(namesOfGroups.size()-1)] << colBack << sampleMetadata[lookup->size()-1] << "}\n" + spaces + "],\n";
        
        out << spaces + "\"matrix_type\": \"" << format << "\",\n" + spaces + "\"matrix_element_type\": \"float\",\n";
        out <<  spaces + "\"shape\": [" << numBins << "," << lookup->size() << "],\n";
        out << spaces + "\"data\":  [";
        
        vector<string> dataRows;
        if (format == "sparse") {
            /*"data":[[0,2,1],
             [1,0,5],
             [1,1,1],
             [1,3,2],
             [1,4,3],
             [1,5,1],
             [2,2,1],
             [2,3,4],
             [2,4,2],
             [3,0,2],
             [3,1,1],
             [3,2,1],
             [3,5,1],
             [4,1,1],
             [4,2,1]
             ]*/
            string output = "";
            for (int i = 0; i < lookup->getNumBins(); i++) {
                
                if (m->getControl_pressed()) { out.close(); return 0; }
                
                vector<float> binAbund = lookup->getOTU(i);
                for (int j = 0; j < binAbund.size(); j++) {
                    float abund = binAbund[j];
                    string binInfo = "[" + toString(i) + "," + toString(j) + "," + toString(abund) + "]";
                    //only print non zero values
                    if (abund != 0) { dataRows.push_back(binInfo); }
                }
            }
        }else {
            
            /* "matrix_type": "dense",
             "matrix_element_type": "int",
             "shape": [5,6],
             "data":  [[0,0,1,0,0,0],
             [5,1,0,2,3,1],
             [0,0,1,4,2,0],
             [2,1,1,0,0,1],
             [0,1,1,0,0,0]]*/
            
            for (int i = 0; i < lookup->getNumBins(); i++) {
                
                if (m->getControl_pressed()) { out.close(); return 0; }
                
                string binInfo = "[";
                vector<float> binAbund = lookup->getOTU(i);
                for (int j = 0; j < binAbund.size()-1; j++) {  binInfo += toString(binAbund[j]) + ","; }
                binInfo += toString(binAbund[binAbund.size()-1]) + "]";
                dataRows.push_back(binInfo);
            }

        }
        
        for (int i = 0; i < dataRows.size()-1; i++) {
            out << dataRows[i] << ",\n" + spaces  + spaces;
        }
        out << dataRows[dataRows.size()-1] << "]\n";
        
        out << "}\n";
        out.close();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeBiomCommand", "getBiom");
        exit(1);
    }
}
//**********************************************************************************************************************
vector<string> MakeBiomCommand::getMetaData(SharedRAbundVectors*& lookup){
	try {
        vector<string> metadata;
        
        if (contaxonomyfile == "") { for (int i = 0; i < lookup->getNumBins(); i++) {  metadata.push_back("null");  } }
        else {
            
            //read constaxonomy file storing in a map, otulabel -> taxonomy
            //constaxonomy file will most likely contain more labels than the shared file, because sharedfile could have been subsampled.
            ifstream in;
            util.openInputFile(contaxonomyfile, in);
            
            //grab headers
            util.getline(in); util.gobble(in);
            
            string otuLabel, tax;
            int size;
            vector<string> otuLabels;
            vector<string> taxs;
            while (!in.eof()) {
                
                if (m->getControl_pressed()) { in.close(); return metadata; }
                
                in >> otuLabel; util.gobble(in);
                in >> size; util.gobble(in);
                tax = util.getline(in); util.gobble(in);
                
                otuLabels.push_back(otuLabel);
                taxs.push_back(tax);
            }
            in.close();
            
            //should the labels be Otu001 or PhyloType001
            string firstBin = lookup->getOTUNames()[0];
            string binTag = "Otu";
            if ((firstBin.find("Otu")) == string::npos) { binTag = "PhyloType";  }
            
            //convert list file bin labels to shared file bin labels
            //parse tax strings
            //save in map
            map<string, string> labelTaxMap;
            string snumBins = toString(otuLabels.size());
            for (int i = 0; i < otuLabels.size(); i++) {  
                
                if (m->getControl_pressed()) { return metadata; }
                
                //if there is a bin label use it otherwise make one
                if (util.isContainingOnlyDigits(otuLabels[i])) {
                    string binLabel = binTag;
                    string sbinNumber = otuLabels[i];
                    if (sbinNumber.length() < snumBins.length()) { 
                        int diff = snumBins.length() - sbinNumber.length();
                        for (int h = 0; h < diff; h++) { binLabel += "0"; }
                    }
                    binLabel += sbinNumber;
                    binLabel = util.getSimpleLabel(binLabel);
                    labelTaxMap[binLabel] = taxs[i];
                }else {
                    map<string, string>::iterator it = labelTaxMap.find(util.getSimpleLabel(otuLabels[i]));
                    if (it == labelTaxMap.end()) {
                        labelTaxMap[util.getSimpleLabel(otuLabels[i])] = taxs[i];
                    }else {
                        m->mothurOut("[ERROR]: Cannot add OTULabel " +  otuLabels[i] + " because it's simple label " + util.getSimpleLabel(otuLabels[i]) + " has already been added and will result in downstream errors. Have you mixed mothur labels and non mothur labels? To make the files work well together and backwards compatible mothur treats 1, OTU01, OTU001, OTU0001 all the same. We do this by removing any non numeric characters and leading zeros. For eaxample: Otu000018 and OtuMY18 both map to 18.\n"); m->setControl_pressed(true);
                    }
                }
            }
            
            //merges OTUs classified to same gg otuid, sets otulabels to gg otuids, averages confidence scores of merged otus.  overwritting of otulabels is fine because constaxonomy only allows for one label to be processed.  If this assumption changes, could cause bug.
            if (picrust) {  getGreenGenesOTUIDs(lookup, labelTaxMap);  }
            
            //{"taxonomy":["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriaceae", "g__Escherichia", "s__"]}
            
            //traverse the binLabels forming the metadata strings and saving them
            //make sure to sanity check
            map<string, string>::iterator it;
            vector<string> currentLabels = lookup->getOTUNames();
            for (int i = 0; i < lookup->getNumBins(); i++) {
                
                if (m->getControl_pressed()) { return metadata; }
                
                it = labelTaxMap.find(util.getSimpleLabel(currentLabels[i]));
                
                if (it == labelTaxMap.end()) { m->mothurOut("[ERROR]: can't find taxonomy information for " + currentLabels[i] + ".\n"); m->setControl_pressed(true); }
                else {
                    vector<string> bootstrapValues;
                    string data = "{\"taxonomy\":[";
            
                    vector<string> scores;
                    vector<string> taxonomies = parseTax(it->second, scores);
                    
                    for (int j = 0; j < taxonomies.size()-1; j ++) { data += "\"" + taxonomies[j] + "\", "; }
                    data += "\"" + taxonomies[taxonomies.size()-1] + "\"]";
                    
                    //add bootstrap values if available
                    if (scores[0] != "null") {
                        data += ", \"bootstrap\":[";
                        
                        for (int j = 0; j < scores.size()-1; j ++) { data += scores[j] + ", "; }
                        data += scores[scores.size()-1] + "]";

                    }
                    data += "}";
                    
                    metadata.push_back(data);
                }
            }
        }
        
        return metadata;
        
    }
	catch(exception& e) {
		m->errorOut(e, "MakeBiomCommand", "getMetadata");
		exit(1);
	}

}
//**********************************************************************************************************************
vector<string> MakeBiomCommand::getMetaData(SharedRAbundFloatVectors*& lookup){
    try {
        vector<string> metadata;
        
        if (contaxonomyfile == "") { for (int i = 0; i < lookup->getNumBins(); i++) {  metadata.push_back("null");  } }
        else {
            
            //read constaxonomy file storing in a map, otulabel -> taxonomy
            //constaxonomy file will most likely contain more labels than the shared file, because sharedfile could have been subsampled.
            ifstream in;
            util.openInputFile(contaxonomyfile, in);
            
            //grab headers
            util.getline(in); util.gobble(in);
            
            string otuLabel, tax;
            int size;
            vector<string> otuLabels;
            vector<string> taxs;
            while (!in.eof()) {
                
                if (m->getControl_pressed()) { in.close(); return metadata; }
                
                in >> otuLabel; util.gobble(in);
                in >> size; util.gobble(in);
                tax = util.getline(in); util.gobble(in);
                
                otuLabels.push_back(otuLabel);
                taxs.push_back(tax);
            }
            in.close();
            
            //should the labels be Otu001 or PhyloType001
            string firstBin = lookup->getOTUNames()[0];
            string binTag = "Otu";
            if ((firstBin.find("Otu")) == string::npos) { binTag = "PhyloType";  }
            
            //convert list file bin labels to shared file bin labels
            //parse tax strings
            //save in map
            map<string, string> labelTaxMap;
            string snumBins = toString(otuLabels.size());
            for (int i = 0; i < otuLabels.size(); i++) {
                
                if (m->getControl_pressed()) { return metadata; }
                
                //if there is a bin label use it otherwise make one
                if (util.isContainingOnlyDigits(otuLabels[i])) {
                    string binLabel = binTag;
                    string sbinNumber = otuLabels[i];
                    if (sbinNumber.length() < snumBins.length()) {
                        int diff = snumBins.length() - sbinNumber.length();
                        for (int h = 0; h < diff; h++) { binLabel += "0"; }
                    }
                    binLabel += sbinNumber;
                    binLabel = util.getSimpleLabel(binLabel);
                    labelTaxMap[binLabel] = taxs[i];
                }else {  labelTaxMap[util.getSimpleLabel(otuLabels[i])] = taxs[i]; }
            }
            
            //merges OTUs classified to same gg otuid, sets otulabels to gg otuids, averages confidence scores of merged otus.  overwritting of otulabels is fine because constaxonomy only allows for one label to be processed.  If this assumption changes, could cause bug.
            if (picrust) {  getGreenGenesOTUIDs(lookup, labelTaxMap);  }
            
            //{"taxonomy":["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriaceae", "g__Escherichia", "s__"]}
            
            //traverse the binLabels forming the metadata strings and saving them
            //make sure to sanity check
            map<string, string>::iterator it;
            vector<string> currentLabels = lookup->getOTUNames();
            for (int i = 0; i < lookup->getNumBins(); i++) {
                
                if (m->getControl_pressed()) { return metadata; }
                
                it = labelTaxMap.find(util.getSimpleLabel(currentLabels[i]));
                
                if (it == labelTaxMap.end()) { m->mothurOut("[ERROR]: can't find taxonomy information for " + currentLabels[i] + ".\n"); m->setControl_pressed(true); }
                else {
                    vector<string> bootstrapValues;
                    string data = "{\"taxonomy\":[";
                    
                    vector<string> scores;
                    vector<string> taxonomies = parseTax(it->second, scores);
                    
                    for (int j = 0; j < taxonomies.size()-1; j ++) { data += "\"" + taxonomies[j] + "\", "; }
                    data += "\"" + taxonomies[taxonomies.size()-1] + "\"]";
                    
                    //add bootstrap values if available
                    if (scores[0] != "null") {
                        data += ", \"bootstrap\":[";
                        
                        for (int j = 0; j < scores.size()-1; j ++) { data += scores[j] + ", "; }
                        data += scores[scores.size()-1] + "]";
                        
                    }
                    data += "}";
                    
                    metadata.push_back(data);
                }
            }
        }
        
        return metadata;
        
    }
    catch(exception& e) {
        m->errorOut(e, "MakeBiomCommand", "getMetadata");
        exit(1);
    }
    
}
//**********************************************************************************************************************
int MakeBiomCommand::getGreenGenesOTUIDs(SharedRAbundVectors*& lookup, map<string, string>& labelTaxMap){
	try {
        //read reftaxonomy
        PhyloTree phylo(referenceTax);
        
        //read otu map file
        map<string, string> otuMap = readGGOtuMap(); //maps reference ID -> OTU ID
        
        if (m->getControl_pressed()) { return 0; }
        
        map<string, vector<string> > ggOTUIDs;
        //loop through otu taxonomies
        for (map<string, string>::iterator it = labelTaxMap.begin(); it != labelTaxMap.end(); it++) { //maps label -> consensus taxonomy
            if (m->getControl_pressed()) { break; }
            
            string OTUTaxonomy = it->second;
            
            //remove confidences
            util.removeConfidences(OTUTaxonomy);
            
            //remove unclassifieds to match template
            int thisPos = OTUTaxonomy.find("unclassified;"); //"Porphyromonadaceae"_unclassified;
            if (thisPos != string::npos) {
                OTUTaxonomy = OTUTaxonomy.substr(0, thisPos);
                thisPos = OTUTaxonomy.find_last_of(";"); //remove rest of parent taxon
                if (thisPos != string::npos) { OTUTaxonomy = OTUTaxonomy.substr(0, thisPos); }
                OTUTaxonomy += ";";
            }
           
            //get list of reference ids that map to this taxonomy
            vector<string> referenceIds = phylo.getSeqs(OTUTaxonomy);
            
            if (m->getControl_pressed()) { break; }
            
            //look for each one in otu map to find match
            string otuID = "not found";
            string referenceString = "";
            for (int i = 0; i < referenceIds.size(); i++) {
                referenceString += referenceIds[i] + " ";
                map<string, string>::iterator itMap = otuMap.find(referenceIds[i]);
                if (itMap != otuMap.end()) { //found it
                    otuID = itMap->second;
                    i += referenceIds.size(); //stop looking
                }
            }
            
            //if found, add otu to ggOTUID list
            if (otuID != "not found") {
                map<string, vector<string> >::iterator itGG = ggOTUIDs.find(otuID);
                if (itGG == ggOTUIDs.end()) {
                    vector<string> temp; temp.push_back(it->first); //save mothur OTU label
                    ggOTUIDs[otuID] = temp;
                }else { ggOTUIDs[otuID].push_back(it->first); } //add mothur OTU label to list
            }else {  m->mothurOut("[ERROR]: could not find OTUId for " + it->second + ". Its reference sequences are " + referenceString + ".\n"); m->setControl_pressed(true); }
            
        }
        
        vector<SharedRAbundVector*> newLookup;
        vector<string> namesOfGroups = lookup->getNamesGroups();
		for (int i = 0; i < namesOfGroups.size(); i++) {
			SharedRAbundVector* temp = new SharedRAbundVector();
			temp->setLabel(lookup->getLabel());
			temp->setGroup(namesOfGroups[i]);
			newLookup.push_back(temp);
		}
		
        map<string, int> labelIndex;
        vector<string> currentLabels = lookup->getOTUNames();
		for (int i = 0; i < currentLabels.size(); i++) {  labelIndex[util.getSimpleLabel(currentLabels[i])] = i; }
        
        vector<string> newBinLabels;
        map<string, string> newLabelTaxMap;
        //loop through ggOTUID list combining mothur otus and adjusting labels
        //ggOTUIDs = 16097 -> <OTU01, OTU10, OTU22>
        
        for (map<string, vector<string> >::iterator itMap = ggOTUIDs.begin(); itMap != ggOTUIDs.end(); itMap++) {
            if (m->getControl_pressed()) { delete lookup; return 0; }
            
            //set new gg otu id to taxonomy. OTU01 -> k__Bacteria becomes 16097 -> k__Bacteria
            //find taxonomy of this otu
            map<string, string>::iterator it = labelTaxMap.find(util.getSimpleLabel(itMap->second[0]));
            vector<string> scores;
            vector<string> taxonomies = parseTax(it->second, scores);
            
            //merge/set OTU abundances for new merged OTU - abunds[i] = new OTU total for group i
            vector<int> abunds; abunds.resize(lookup->size(), 0);
            string mergeString = "";
            vector<float> boots; boots.resize(scores.size(), 0);
            bool scoresNULL = false;
            for (int j = 0; j < itMap->second.size(); j++) { //<OTU01, OTU10, OTU22>
                
                if (scores[0] != "null") {
                    //merge bootstrap scores
                    vector<string> scores;
                    vector<string> taxonomies = parseTax(it->second, scores);
                    for (int i = 0; i < boots.size(); i++) {
                        if (scores[i] == "null") { scoresNULL = true; break; }
                        else {
                            float tempScore; util.mothurConvert(scores[i], tempScore);
                            boots[i] += tempScore;
                        }
                    }
                }else { scoresNULL = true; }
                
                //merge abunds
                mergeString += (itMap->second)[j] + " ";
                for (int i = 0; i < lookup->size(); i++) { abunds[i] += lookup->get(labelIndex[util.getSimpleLabel((itMap->second)[j])], namesOfGroups[i]); }
            }
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: merging " + mergeString + " for ggOTUid = " + itMap->first + ".\n");  }
            
            //average scores
            //add merged otu to new lookup
            //assemble new taxonomy
            string newTaxString = "";
            if (!scoresNULL) {
                for (int j = 0; j < boots.size(); j++) { boots[j] /= (float) itMap->second.size(); }
            
                for (int j = 0; j < boots.size(); j++) { newTaxString += taxonomies[j] + "(" + toString(boots[j]) + ");"; }
            }else {
                for (int j = 0; j < taxonomies.size(); j++) { newTaxString += taxonomies[j] + ";"; }
            }
            
            //set new gg otu id to taxonomy. OTU01 -> k__Bacteria becomes 16097 -> k__Bacteria
            //find taxonomy of this otu
            newLabelTaxMap[itMap->first] = newTaxString;
            
            //add merged otu to new lookup
            for (int j = 0; j < abunds.size(); j++) { newLookup[j]->push_back(abunds[j]); }
            
            //saved otu label
            newBinLabels.push_back(itMap->first);
        }
		
        lookup->clear();
        for (int i = 0; i < newLookup.size(); i++) { lookup->push_back(newLookup[i]);  }
        lookup->eliminateZeroOTUS();
		
		lookup->setOTUNames(newBinLabels);
        labelTaxMap = newLabelTaxMap;
        
        map<string, string> variables;
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(sharedfile));
        variables["[distance]"] = lookup->getLabel();
        string outputFileName = getOutputFileName("shared",variables);
		ofstream out;
		util.openOutputFile(outputFileName, out);
		outputNames.push_back(outputFileName); outputTypes["shared"].push_back(outputFileName);
        
        lookup->printHeaders(out);
        lookup->print(out);
        out.close();

        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "MakeBiomCommand", "getGreenGenesOTUIDs");
		exit(1);
	}
    
}
//**********************************************************************************************************************
int MakeBiomCommand::getGreenGenesOTUIDs(SharedRAbundFloatVectors*& lookup, map<string, string>& labelTaxMap){
    try {
        //read reftaxonomy
        PhyloTree phylo(referenceTax);
        
        //read otu map file
        map<string, string> otuMap = readGGOtuMap(); //maps reference ID -> OTU ID
        
        if (m->getControl_pressed()) { return 0; }
        
        map<string, vector<string> > ggOTUIDs;
        //loop through otu taxonomies
        for (map<string, string>::iterator it = labelTaxMap.begin(); it != labelTaxMap.end(); it++) { //maps label -> consensus taxonomy
            if (m->getControl_pressed()) { break; }
            
            string OTUTaxonomy = it->second;
            
            //remove confidences
            util.removeConfidences(OTUTaxonomy);
            
            //remove unclassifieds to match template
            int thisPos = OTUTaxonomy.find("unclassified;"); //"Porphyromonadaceae"_unclassified;
            if (thisPos != string::npos) {
                OTUTaxonomy = OTUTaxonomy.substr(0, thisPos);
                thisPos = OTUTaxonomy.find_last_of(";"); //remove rest of parent taxon
                if (thisPos != string::npos) {
                    OTUTaxonomy = OTUTaxonomy.substr(0, thisPos);
                }
            }
            
            //get list of reference ids that map to this taxonomy
            vector<string> referenceIds = phylo.getSeqs(OTUTaxonomy);
            
            if (m->getControl_pressed()) { break; }
            
            //look for each one in otu map to find match
            string otuID = "not found";
            string referenceString = "";
            for (int i = 0; i < referenceIds.size(); i++) {
                referenceString += referenceIds[i] + " ";
                map<string, string>::iterator itMap = otuMap.find(referenceIds[i]);
                if (itMap != otuMap.end()) { //found it
                    otuID = itMap->second;
                    i += referenceIds.size(); //stop looking
                }
            }
            
            //if found, add otu to ggOTUID list
            if (otuID != "not found") {
                map<string, vector<string> >::iterator itGG = ggOTUIDs.find(otuID);
                if (itGG == ggOTUIDs.end()) {
                    vector<string> temp; temp.push_back(it->first); //save mothur OTU label
                    ggOTUIDs[otuID] = temp;
                }else { ggOTUIDs[otuID].push_back(it->first); } //add mothur OTU label to list
            }else {  m->mothurOut("[ERROR]: could not find OTUId for " + it->second + ". Its reference sequences are " + referenceString + ".\n"); m->setControl_pressed(true); }
            
        }
        
        
        vector<SharedRAbundFloatVector*> newLookup;
        vector<string> namesOfGroups = lookup->getNamesGroups();
        for (int i = 0; i < namesOfGroups.size(); i++) {
            SharedRAbundFloatVector* temp = new SharedRAbundFloatVector();
            temp->setLabel(lookup->getLabel());
            temp->setGroup(namesOfGroups[i]);
            newLookup.push_back(temp);
        }
        
        map<string, int> labelIndex;
        vector<string> currentLabels = lookup->getOTUNames();
        for (int i = 0; i < currentLabels.size(); i++) {  labelIndex[util.getSimpleLabel(currentLabels[i])] = i; }
        
        vector<string> newBinLabels;
        map<string, string> newLabelTaxMap;
        //loop through ggOTUID list combining mothur otus and adjusting labels
        //ggOTUIDs = 16097 -> <OTU01, OTU10, OTU22>
        
        for (map<string, vector<string> >::iterator itMap = ggOTUIDs.begin(); itMap != ggOTUIDs.end(); itMap++) {
            if (m->getControl_pressed()) { for (int j = 0; j < newLookup.size(); j++) {  delete newLookup[j];  } return 0; }
            
            //set new gg otu id to taxonomy. OTU01 -> k__Bacteria becomes 16097 -> k__Bacteria
            //find taxonomy of this otu
            map<string, string>::iterator it = labelTaxMap.find(util.getSimpleLabel(itMap->second[0]));
            vector<string> scores;
            vector<string> taxonomies = parseTax(it->second, scores);
            
            //merge/set OTU abundances
            vector<float> abunds; abunds.resize(lookup->size(), 0.0);
            string mergeString = "";
            vector<float> boots; boots.resize(scores.size(), 0.0);
            bool scoresNULL = false;
            for (int j = 0; j < itMap->second.size(); j++) { //<OTU01, OTU10, OTU22>
                
                if (scores[0] != "null") {
                    //merge bootstrap scores
                    vector<string> scores;
                    vector<string> taxonomies = parseTax(it->second, scores);
                    for (int i = 0; i < boots.size(); i++) {
                        if (scores[i] == "null") { scoresNULL = true; break; }
                        else {
                            float tempScore; util.mothurConvert(scores[i], tempScore);
                            boots[i] += tempScore;
                        }
                    }
                }else { scoresNULL = true; }
                
                //merge abunds
                mergeString += (itMap->second)[j] + " ";
                for (int i = 0; i < lookup->size(); i++) { abunds[i] += lookup->get(labelIndex[util.getSimpleLabel((itMap->second)[j])], namesOfGroups[i]); }
            }
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: merging " + mergeString + " for ggOTUid = " + itMap->first + ".\n");  }
            
            //average scores
            //add merged otu to new lookup
            string newTaxString = "";
            if (!scoresNULL) {
                for (int j = 0; j < boots.size(); j++) { boots[j] /= (float) itMap->second.size(); }
                
                //assemble new taxomoy
                for (int j = 0; j < boots.size(); j++) {
                    newTaxString += taxonomies[j] + "(" + toString(boots[j]) + ");";
                }
            }else {
                //assemble new taxomoy
                for (int j = 0; j < taxonomies.size(); j++) {
                    newTaxString += taxonomies[j] + ";";
                }
            }
            
            //set new gg otu id to taxonomy. OTU01 -> k__Bacteria becomes 16097 -> k__Bacteria
            //find taxonomy of this otu
            newLabelTaxMap[itMap->first] = newTaxString;
            
            //add merged otu to new lookup
            for (int j = 0; j < abunds.size(); j++) { newLookup[j]->push_back(abunds[j]); }
            
            //saved otu label
            newBinLabels.push_back(itMap->first);
        }
        
        lookup->clear();
        for (int i = 0; i < newLookup.size(); i++) { lookup->push_back(newLookup[i]);  }
        lookup->eliminateZeroOTUS();
        
        lookup->setOTUNames(newBinLabels);
        labelTaxMap = newLabelTaxMap;
        
        map<string, string> variables;
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(inputFileName));
        variables["[distance]"] = lookup->getLabel();
        string outputFileName = getOutputFileName("relabund",variables);
        ofstream out;
    
        util.openOutputFile(outputFileName, out);
        
        outputNames.push_back(outputFileName); outputTypes["relabund"].push_back(outputFileName);
        
        lookup->printHeaders(out);
        lookup->print(out);
        out.close();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeBiomCommand", "getGreenGenesOTUIDs");
        exit(1);
    }
    
}

//**********************************************************************************************************************
map<string, string> MakeBiomCommand::readGGOtuMap(){
	try {
        map<string, string> otuMap;
        
        ifstream in;
        util.openInputFile(picrustOtuFile, in);
        
        //map referenceIDs -> otuIDs
        //lines look like:
        //16097	671376	616121	533566	683683	4332909	4434717	772666	611808	695209
        while(!in.eof()) {
            if (m->getControl_pressed()) { break; }
            
            string line = util.getline(in); util.gobble(in);
            vector<string> pieces = util.splitWhiteSpace(line);
            
            if (pieces.size() != 0) {
                string otuID = pieces[1];
                for (int i = 1; i < pieces.size(); i++) {  otuMap[pieces[i]] = otuID; }
            }
        }
        in.close();
        
        return otuMap;
    }
	catch(exception& e) {
		m->errorOut(e, "MakeBiomCommand", "readGGOtuMap");
		exit(1);
	}
    
}
//**********************************************************************************************************************
int MakeBiomCommand::getSampleMetaData(SharedRAbundVectors*& lookup){
	try {
        sampleMetadata.clear();
        if (metadatafile == "") {  for (int i = 0; i < lookup->size(); i++) {  sampleMetadata.push_back("null");  } }
        else {
            ifstream in;
            util.openInputFile(metadatafile, in);
            
            vector<string> groupNames, metadataLabels;
            map<string, vector<string> > lines;
            
            string headerLine = util.getline(in); util.gobble(in);
            vector<string> pieces = util.splitWhiteSpace(headerLine);
            
            //save names of columns you are reading
            for (int i = 1; i < pieces.size(); i++) {
                metadataLabels.push_back(pieces[i]);
            }
            int count = metadataLabels.size();
            
            //read rest of file
            while (!in.eof()) {
                
                if (m->getControl_pressed()) { in.close(); return 0; }
                
                string group = "";
                in >> group; util.gobble(in);
                groupNames.push_back(group);
                
                string line = util.getline(in); util.gobble(in);
                vector<string> thisPieces = util.splitWhiteSpaceWithQuotes(line);
                
                if (m->getDebug()) {  m->mothurOut("[DEBUG]: " + group + " " + util.getStringFromVector(thisPieces, ", ") + "\n"); }
        
                if (thisPieces.size() != count) { m->mothurOut("[ERROR]: expected " + toString(count) + " items of data for sample " + group + " read " + toString(thisPieces.size()) + ", quitting.\n"); }
                else {  if (util.inUsersGroups(group, Groups)) { lines[group] = thisPieces; } }
                
                util.gobble(in);
            }
            in.close();
            
            map<string, vector<string> >::iterator it;
            vector<string> namesOfGroups = lookup->getNamesGroups();
            for (int i = 0; i < namesOfGroups.size(); i++) {
                
                if (m->getControl_pressed()) { return 0; }
                
                it = lines.find(namesOfGroups[i]);
                
                if (it == lines.end()) { m->mothurOut("[ERROR]: can't find metadata information for " + namesOfGroups[i] + ", quitting.\n"); m->setControl_pressed(true); }
                else {
                    vector<string> values = it->second;
                    
                    string data = "{";
                    for (int j = 0; j < metadataLabels.size()-1; j++) { 
                        values[j] = util.removeQuotes(values[j]); 
                        data += "\"" + metadataLabels[j] + "\":\"" + values[j] + "\", "; 
                    }
                    values[metadataLabels.size()-1] = util.removeQuotes(values[metadataLabels.size()-1]);
                    data += "\"" + metadataLabels[metadataLabels.size()-1] + "\":\"" + values[metadataLabels.size()-1] + "\"}";
                    sampleMetadata.push_back(data);
                }
            }
        }
        
        return 0;
        
    }
	catch(exception& e) {
		m->errorOut(e, "MakeBiomCommand", "getSampleMetaData");
		exit(1);
	}
    
}
//**********************************************************************************************************************
int MakeBiomCommand::getSampleMetaData(SharedRAbundFloatVectors*& lookup){
    try {
        sampleMetadata.clear();
        if (metadatafile == "") {  for (int i = 0; i < lookup->size(); i++) {  sampleMetadata.push_back("null");  } }
        else {
            ifstream in;
            util.openInputFile(metadatafile, in);
            
            vector<string> groupNames, metadataLabels;
            map<string, vector<string> > lines;
            
            string headerLine = util.getline(in); util.gobble(in);
            vector<string> pieces = util.splitWhiteSpace(headerLine);
            
            //save names of columns you are reading
            for (int i = 1; i < pieces.size(); i++) {
                metadataLabels.push_back(pieces[i]);
            }
            int count = metadataLabels.size();
            
            //read rest of file
            while (!in.eof()) {
                
                if (m->getControl_pressed()) { in.close(); return 0; }
                
                string group = "";
                in >> group; util.gobble(in);
                groupNames.push_back(group);
                
                string line = util.getline(in); util.gobble(in);
                vector<string> thisPieces = util.splitWhiteSpaceWithQuotes(line);
                
                if (thisPieces.size() != count) { m->mothurOut("[ERROR]: expected " + toString(count) + " items of data for sample " + group + " read " + toString(thisPieces.size()) + ", quitting.\n"); }
                else {  if (util.inUsersGroups(group, Groups)) { lines[group] = thisPieces; } }
                
                util.gobble(in);
            }
            in.close();
            
            map<string, vector<string> >::iterator it;
            vector<string> namesOfGroups = lookup->getNamesGroups();
            for (int i = 0; i < namesOfGroups.size(); i++) {
                
                if (m->getControl_pressed()) { return 0; }
                
                it = lines.find(namesOfGroups[i]);
                
                if (it == lines.end()) { m->mothurOut("[ERROR]: can't find metadata information for " + namesOfGroups[i] + ", quitting.\n"); m->setControl_pressed(true); }
                else {
                    vector<string> values = it->second;
                    
                    string data = "{";
                    for (int j = 0; j < metadataLabels.size()-1; j++) {
                        values[j] = util.removeQuotes(values[j]);
                        data += "\"" + metadataLabels[j] + "\":\"" + values[j] + "\", ";
                    }
                    values[metadataLabels.size()-1] = util.removeQuotes(values[metadataLabels.size()-1]);
                    data += "\"" + metadataLabels[metadataLabels.size()-1] + "\":\"" + values[metadataLabels.size()-1] + "\"}";
                    sampleMetadata.push_back(data);
                }
            }
        }
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "MakeBiomCommand", "getSampleMetaData");
        exit(1);
    }
    
}
/**************************************************************************************************/
//returns {Bacteria, Bacteroidetes, ..} and scores is filled with {100, 98, ...} or {null, null, null}
vector<string> MakeBiomCommand::parseTax(string tax, vector<string>& scores) {
	try {
		
		string taxon;
        vector<string> taxs;
		
		while (tax.find_first_of(';') != -1) {
			
			if (m->getControl_pressed()) { return taxs; }
			
			//get taxon
			taxon = tax.substr(0,tax.find_first_of(';'));
            
			int pos = taxon.find_last_of('(');
			if (pos != -1) {
				//is it a number?
				int pos2 = taxon.find_last_of(')');
				if (pos2 != -1) {
					string confidenceScore = taxon.substr(pos+1, (pos2-(pos+1)));
					if (util.isNumeric1(confidenceScore)) {
						taxon = taxon.substr(0, pos); //rip off confidence 
                        scores.push_back(confidenceScore);
					}else{ scores.push_back("null"); }
				}
			}else{ scores.push_back("null"); }
			
            //strip "" if they are there
            pos = taxon.find("\"");
            if (pos != string::npos) {
                string newTax = "";
                for (int k = 0; k < taxon.length(); k++) {
                    if (taxon[k] != '\"') { newTax += taxon[k]; }
                }
                taxon = newTax;
            }
            
            //look for bootstrap value
			taxs.push_back(taxon);
            tax = tax.substr(tax.find_first_of(';')+1, tax.length());
		}
		
		return taxs;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeBiomCommand", "parseTax");
		exit(1);
	}
}

//**********************************************************************************************************************



