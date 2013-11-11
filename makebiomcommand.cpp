//
//  makebiomcommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/16/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "makebiomcommand.h"
#include "sharedrabundvector.h"
#include "inputdata.h"
#include "sharedutilities.h"
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
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none","biom",false,true,true); parameters.push_back(pshared);
        CommandParameter pcontaxonomy("constaxonomy", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pcontaxonomy);
        CommandParameter preference("reftaxonomy", "InputTypes", "", "", "none", "none", "refPi","",false,false); parameters.push_back(preference);
        CommandParameter pmetadata("metadata", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pmetadata);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
        CommandParameter ppicrust("picrust", "InputTypes", "", "", "none", "none", "refPi","shared",false,false); parameters.push_back(ppicrust);
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
		helpString += "The make.biom command parameters are shared, contaxonomy, metadata, groups, matrixtype, picrust, reftaxonomy and label.  shared is required, unless you have a valid current file.\n"; //
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
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
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
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
                
                it = parameters.find("constaxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["constaxonomy"] = inputDir + it->second;		}
				}
                
                it = parameters.find("reftaxonomy");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["reftaxonomy"] = inputDir + it->second;		}
				}
                
                it = parameters.find("picrust");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["picrust"] = inputDir + it->second;		}
				}
                
                it = parameters.find("metadata");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["metadata"] = inputDir + it->second;		}
				}
			}
            
			//get shared file
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }	
			else if (sharedfile == "not found") { 
				//if there is a current shared file, use it
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current sharedfile and the shared parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setSharedFile(sharedfile); }
			
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(sharedfile);		}
            
            contaxonomyfile = validParameter.validFile(parameters, "constaxonomy", true);
			if (contaxonomyfile == "not found") {  contaxonomyfile = "";  }
			else if (contaxonomyfile == "not open") { contaxonomyfile = ""; abort = true; }
            
            referenceTax = validParameter.validFile(parameters, "reftaxonomy", true);
			if (referenceTax == "not found") {  referenceTax = "";  }
			else if (referenceTax == "not open") { referenceTax = ""; abort = true; }
            
            picrustOtuFile = validParameter.validFile(parameters, "picrust", true);
			if (picrustOtuFile == "not found") {  picrustOtuFile = "";  }
			else if (picrustOtuFile == "not open") { picrustOtuFile = ""; abort = true; }

            metadatafile = validParameter.validFile(parameters, "metadata", true);
			if (metadatafile == "not found") {  metadatafile = "";  }
			else if (metadatafile == "not open") { metadatafile = ""; abort = true; }
            
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				m->splitAtDash(groups, Groups);
				m->setGroups(Groups);
			}
			
            if (picrustOtuFile != "") {
                picrust=true;
                if (contaxonomyfile == "") {  m->mothurOut("[ERROR]: the constaxonomy parameter is required with the picrust parameter, aborting."); m->mothurOutEndLine(); abort = true;  }
                if (referenceTax == "") {  m->mothurOut("[ERROR]: the reftaxonomy parameter is required with the picrust parameter, aborting."); m->mothurOutEndLine(); abort = true;  }
            }else { picrust=false; }
            
            if ((contaxonomyfile != "") && (labels.size() > 1)) { m->mothurOut("[ERROR]: the contaxonomy parameter cannot be used with multiple labels."); m->mothurOutEndLine(); abort = true; }
            
			format = validParameter.validFile(parameters, "matrixtype", false);				if (format == "not found") { format = "sparse"; }
			
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
        
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
            
		InputData input(sharedfile, "sharedfile");
		vector<SharedRAbundVector*> lookup = input.getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
        
        //if user did not specify a label, then use first one
        if ((contaxonomyfile != "") && (labels.size() == 0)) {
            allLines = 0;
            labels.insert(lastLabel);
        }
		
        getSampleMetaData(lookup);
        
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
        
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  return 0; }
            
			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			
                
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				getBiom(lookup);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}
			
			if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup[0]->getLabel();
                
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  
				lookup = input.getSharedRAbundVectors(lastLabel);
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				
				getBiom(lookup);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				
				//restore real lastlabel to save below
				lookup[0]->setLabel(saveLabel);
			}
			
			lastLabel = lookup[0]->getLabel();
            
			//prevent memory leak and get next set
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i]; lookup[i] = NULL; }
			lookup = input.getSharedRAbundVectors();				
		}
		
        if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); }  return 0; }     
        
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
		if (needToRun == true)  {
			for (int i = 0; i < lookup.size(); i++) { if (lookup[i] != NULL) { delete lookup[i]; } }  
			lookup = input.getSharedRAbundVectors(lastLabel);
			
			m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
            getBiom(lookup);
			
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
		}
		
        if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); }  return 0; }     
		
        //set sabund file as new current sabundfile
        string current = "";
		itTypes = outputTypes.find("biom");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setBiomFile(current); }
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
int MakeBiomCommand::getBiom(vector<SharedRAbundVector*>& lookup){
	try {
        map<string, string> variables; 
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(sharedfile));
        variables["[distance]"] = lookup[0]->getLabel();
        string outputFileName = getOutputFileName("biom",variables);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		outputNames.push_back(outputFileName); outputTypes["biom"].push_back(outputFileName);

        string mothurString = "mothur" + toString(m->getVersion());
        time_t rawtime;
        struct tm * timeinfo;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        string dateString = asctime (timeinfo);
        int pos = dateString.find('\n');
        if (pos != string::npos) { dateString = dateString.substr(0, pos);}
        string spaces = "      ";
        
        //standard 
        out << "{\n" + spaces + "\"id\":\"" + sharedfile + "-" + lookup[0]->getLabel() + "\",\n" + spaces + "\"format\": \"Biological Observation Matrix 0.9.1\",\n" + spaces + "\"format_url\": \"http://biom-format.org\",\n";
        out << spaces + "\"type\": \"OTU table\",\n" + spaces + "\"generated_by\": \"" << mothurString << "\",\n" + spaces + "\"date\": \"" << dateString << "\",\n";
        
        
        vector<string> metadata = getMetaData(lookup);
        int numBins = lookup[0]->getNumBins();
        
        if (m->control_pressed) {  out.close(); return 0; }
        
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
        for (int i = 0; i < numBins-1; i++) {
            if (m->control_pressed) {  out.close(); return 0; }
            out << rowFront << m->currentSharedBinLabels[i] << rowBack << metadata[i] << "},\n"; 
        }
        out << rowFront << m->currentSharedBinLabels[(numBins-1)] << rowBack << metadata[(numBins-1)] << "}\n" + spaces + "],\n"; 
       
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
        for (int i = 0; i < lookup.size()-1; i++) {
            if (m->control_pressed) {  out.close(); return 0; }
            out << rowFront << lookup[i]->getGroup() << colBack << sampleMetadata[i] << "},\n";
        }
        out << rowFront << lookup[(lookup.size()-1)]->getGroup() << colBack << sampleMetadata[lookup.size()-1] << "}\n" + spaces + "],\n";
        
        out << spaces + "\"matrix_type\": \"" << format << "\",\n" + spaces + "\"matrix_element_type\": \"int\",\n";
        out <<  spaces + "\"shape\": [" << numBins << "," << lookup.size() << "],\n";
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
            for (int i = 0; i < lookup[0]->getNumBins(); i++) {
                
                if (m->control_pressed) { out.close(); return 0; }
                
                for (int j = 0; j < lookup.size(); j++) {
                    string binInfo = "[" + toString(i) + "," + toString(j) + "," + toString(lookup[j]->getAbundance(i)) + "]";
                    //only print non zero values
                    if (lookup[j]->getAbundance(i) != 0) { dataRows.push_back(binInfo); }
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
            
            for (int i = 0; i < lookup[0]->getNumBins(); i++) {
                
                if (m->control_pressed) { out.close(); return 0; }
                
                string binInfo = "[";
                for (int j = 0; j < lookup.size()-1; j++) {
                    binInfo += toString(lookup[j]->getAbundance(i)) + ",";
                }
                binInfo += toString(lookup[lookup.size()-1]->getAbundance(i)) + "]";
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
vector<string> MakeBiomCommand::getMetaData(vector<SharedRAbundVector*>& lookup){
	try {
        vector<string> metadata;
        
        if (contaxonomyfile == "") { for (int i = 0; i < lookup[0]->getNumBins(); i++) {  metadata.push_back("null");  } }
        else {
            
            //read constaxonomy file storing in a map, otulabel -> taxonomy
            //constaxonomy file will most likely contain more labels than the shared file, because sharedfile could have been subsampled.
            ifstream in;
            m->openInputFile(contaxonomyfile, in);
            
            //grab headers
            m->getline(in); m->gobble(in);
            
            string otuLabel, tax;
            int size;
            vector<string> otuLabels;
            vector<string> taxs;
            while (!in.eof()) {
                
                if (m->control_pressed) { in.close(); return metadata; }
                
                in >> otuLabel >> size >> tax; m->gobble(in);
                
                otuLabels.push_back(otuLabel);
                taxs.push_back(tax);
            }
            in.close();
            
            //should the labels be Otu001 or PhyloType001
            string firstBin = m->currentSharedBinLabels[0];
            string binTag = "Otu";
            if ((firstBin.find("Otu")) == string::npos) { binTag = "PhyloType";  }
            
            //convert list file bin labels to shared file bin labels
            //parse tax strings
            //save in map
            map<string, string> labelTaxMap;
            string snumBins = toString(otuLabels.size());
            for (int i = 0; i < otuLabels.size(); i++) {  
                
                if (m->control_pressed) { return metadata; }
                
                //if there is a bin label use it otherwise make one
                if (m->isContainingOnlyDigits(otuLabels[i])) {
                    string binLabel = binTag;
                    string sbinNumber = otuLabels[i];
                    if (sbinNumber.length() < snumBins.length()) { 
                        int diff = snumBins.length() - sbinNumber.length();
                        for (int h = 0; h < diff; h++) { binLabel += "0"; }
                    }
                    binLabel += sbinNumber;
                    binLabel = m->getSimpleLabel(binLabel);
                    labelTaxMap[binLabel] = taxs[i];
                }else {  labelTaxMap[m->getSimpleLabel(otuLabels[i])] = taxs[i]; }
            }
            
            //merges OTUs classified to same gg otuid, sets otulabels to gg otuids, averages confidence scores of merged otus.  overwritting of otulabels is fine because constaxonomy only allows for one label to be processed.  If this assumption changes, could cause bug.
            if (picrust) {  getGreenGenesOTUIDs(lookup, labelTaxMap);  }
            
            //{"taxonomy":["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriaceae", "g__Escherichia", "s__"]}
            
            //traverse the binLabels forming the metadata strings and saving them
            //make sure to sanity check
            map<string, string>::iterator it;
            for (int i = 0; i < lookup[0]->getNumBins(); i++) {
                
                if (m->control_pressed) { return metadata; }
                
                it = labelTaxMap.find(m->getSimpleLabel(m->currentSharedBinLabels[i]));
                
                if (it == labelTaxMap.end()) { m->mothurOut("[ERROR]: can't find taxonomy information for " + m->currentSharedBinLabels[i] + ".\n"); m->control_pressed = true; }
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
int MakeBiomCommand::getGreenGenesOTUIDs(vector<SharedRAbundVector*>& lookup, map<string, string>& labelTaxMap){
	try {
        //read reftaxonomy
        PhyloTree phylo(referenceTax);
        
        //read otu map file
        map<string, string> otuMap = readGGOtuMap(); //maps reference ID -> OTU ID
        
        if (m->control_pressed) { return 0; }
        
        map<string, vector<string> > ggOTUIDs;
        //loop through otu taxonomies
        for (map<string, string>::iterator it = labelTaxMap.begin(); it != labelTaxMap.end(); it++) { //maps label -> consensus taxonomy
            if (m->control_pressed) { break; }
            
            //get list of reference ids that map to this taxonomy
            vector<string> referenceIds = phylo.getSeqs(it->second);
            
            if (m->control_pressed) { break; }
            
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
            }else {  m->mothurOut("[ERROR]: could not find OTUId for " + it->second + ". Its reference sequences are " + referenceString + ".\n"); m->control_pressed = true; }
            
        }
        
       
        vector<SharedRAbundVector*> newLookup;
		for (int i = 0; i < lookup.size(); i++) {
			SharedRAbundVector* temp = new SharedRAbundVector();
			temp->setLabel(lookup[i]->getLabel());
			temp->setGroup(lookup[i]->getGroup());
			newLookup.push_back(temp);
		}
		
        map<string, int> labelIndex;
		for (int i = 0; i < m->currentSharedBinLabels.size(); i++) {  labelIndex[m->getSimpleLabel(m->currentSharedBinLabels[i])] = i; }
        
        vector<string> newBinLabels;
        map<string, string> newLabelTaxMap;
        //loop through ggOTUID list combining mothur otus and adjusting labels
        //ggOTUIDs = 16097 -> <OTU01, OTU10, OTU22>
        for (map<string, vector<string> >::iterator itMap = ggOTUIDs.begin(); itMap != ggOTUIDs.end(); itMap++) {
			if (m->control_pressed) { for (int j = 0; j < newLookup.size(); j++) {  delete newLookup[j];  } return 0; }
            
            //set new gg otu id to taxonomy. OTU01 -> k__Bacteria becomes 16097 -> k__Bacteria
            //find taxonomy of this otu
            map<string, string>::iterator it = labelTaxMap.find(m->getSimpleLabel(itMap->second[0]));
            vector<string> scores;
            vector<string> taxonomies = parseTax(it->second, scores);
            
            //merge/set OTU abundances
            vector<int> abunds; abunds.resize(lookup.size(), 0);
            string mergeString = "";
            vector<float> boots; boots.resize(scores.size(), 0);
            for (int j = 0; j < itMap->second.size(); j++) { //<OTU01, OTU10, OTU22>
                //merge bootstrap scores
                vector<string> scores;
                vector<string> taxonomies = parseTax(it->second, scores);
                for (int i = 0; i < boots.size(); i++) {
                    float tempScore; m->mothurConvert(scores[i], tempScore);
                    boots[i] += tempScore;
                }
                
                //merge abunds
                mergeString += (itMap->second)[j] + " ";
                for (int i = 0; i < lookup.size(); i++) {
                    abunds[i] += lookup[i]->getAbundance(labelIndex[m->getSimpleLabel((itMap->second)[j])]);
                }
            }
            
            if (m->debug) { m->mothurOut("[DEBUG]: merging " + mergeString + " for ggOTUid = " + itMap->first + ".\n");  }
            
            //average scores
            //add merged otu to new lookup
            for (int j = 0; j < boots.size(); j++) { boots[j] /= (float) itMap->second.size(); }
            
            //assemble new taxomoy
            string newTaxString = "";
            for (int j = 0; j < boots.size(); j++) {
                newTaxString += taxonomies[j] + "(" + toString(boots[j]) + ");";
            }

            //set new gg otu id to taxonomy. OTU01 -> k__Bacteria becomes 16097 -> k__Bacteria
            //find taxonomy of this otu
            newLabelTaxMap[itMap->first] = newTaxString;
            
            //add merged otu to new lookup
            for (int j = 0; j < abunds.size(); j++) { newLookup[j]->push_back(abunds[j], newLookup[j]->getGroup()); }
            
            //saved otu label
            newBinLabels.push_back(itMap->first);
        }
		
		for (int j = 0; j < lookup.size(); j++) {  delete lookup[j];  }
		
		lookup = newLookup;
		m->currentSharedBinLabels = newBinLabels;
        labelTaxMap = newLabelTaxMap;
        
        map<string, string> variables;
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(sharedfile));
        variables["[distance]"] = lookup[0]->getLabel();
        string outputFileName = getOutputFileName("shared",variables);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		outputNames.push_back(outputFileName); outputTypes["shared"].push_back(outputFileName);
        
        lookup[0]->printHeaders(out);
		
		for (int i = 0; i < lookup.size(); i++) {
			out << lookup[i]->getLabel() << '\t' << lookup[i]->getGroup() << '\t';
			lookup[i]->print(out);
		}
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
        m->openInputFile(picrustOtuFile, in);
        
        //map referenceIDs -> otuIDs
        //lines look like:
        //16097	671376	616121	533566	683683	4332909	4434717	772666	611808	695209
        while(!in.eof()) {
            if (m->control_pressed) { break; }
            
            string line = m->getline(in); m->gobble(in);
            vector<string> pieces = m->splitWhiteSpace(line);
            
            if (pieces.size() != 0) {
                string otuID = pieces[0];
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
int MakeBiomCommand::getSampleMetaData(vector<SharedRAbundVector*>& lookup){
	try {
        sampleMetadata.clear();
        if (metadatafile == "") {  for (int i = 0; i < lookup.size(); i++) {  sampleMetadata.push_back("null");  } }
        else {
            ifstream in;
            m->openInputFile(metadatafile, in);
            
            vector<string> groupNames, metadataLabels;
            map<string, vector<string> > lines;
            
            string headerLine = m->getline(in); m->gobble(in);
            vector<string> pieces = m->splitWhiteSpace(headerLine);
            
            //save names of columns you are reading
            for (int i = 1; i < pieces.size(); i++) {
                metadataLabels.push_back(pieces[i]);
            }
            int count = metadataLabels.size();
			
            vector<string> groups = m->getGroups();
            
            //read rest of file
            while (!in.eof()) {
                
                if (m->control_pressed) { in.close(); return 0; }
                
                string group = "";
                in >> group; m->gobble(in);
                groupNames.push_back(group);
                
                string line = m->getline(in); m->gobble(in);
                vector<string> thisPieces = m->splitWhiteSpaceWithQuotes(line);
        
                if (thisPieces.size() != count) { m->mothurOut("[ERROR]: expected " + toString(count) + " items of data for sample " + group + " read " + toString(thisPieces.size()) + ", quitting.\n"); }
                else {  if (m->inUsersGroups(group, groups)) { lines[group] = thisPieces; } }
                
                m->gobble(in);
            }
            in.close();
            
            map<string, vector<string> >::iterator it;
            for (int i = 0; i < lookup.size(); i++) {
                
                if (m->control_pressed) { return 0; }
                
                it = lines.find(lookup[i]->getGroup());
                
                if (it == lines.end()) { m->mothurOut("[ERROR]: can't find metadata information for " + lookup[i]->getGroup() + ", quitting.\n"); m->control_pressed = true; }
                else {
                    vector<string> values = it->second;
                    
                    string data = "{";
                    for (int j = 0; j < metadataLabels.size()-1; j++) { 
                        values[j] = m->removeQuotes(values[j]); 
                        data += "\"" + metadataLabels[j] + "\":\"" + values[j] + "\", "; 
                    }
                    values[metadataLabels.size()-1] = m->removeQuotes(values[metadataLabels.size()-1]);
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
			
			if (m->control_pressed) { return taxs; }
			
			//get taxon
			taxon = tax.substr(0,tax.find_first_of(';'));
            
			int pos = taxon.find_last_of('(');
			if (pos != -1) {
				//is it a number?
				int pos2 = taxon.find_last_of(')');
				if (pos2 != -1) {
					string confidenceScore = taxon.substr(pos+1, (pos2-(pos+1)));
					if (m->isNumeric1(confidenceScore)) {
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



