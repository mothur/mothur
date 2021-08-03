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
        CommandParameter poutput("output", "Multiple", "simple-hdf5", "hdf5", "", "", "","",false,false, true); parameters.push_back(poutput);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        CommandParameter pmatrixtype("matrixtype", "Multiple", "sparse-dense", "sparse", "", "", "","",false,false); parameters.push_back(pmatrixtype);
        
        abort = false; calledHelp = false; allLines = true;
        
        //initialize outputTypes
        vector<string> tempOutNames;
        outputTypes["biom"] = tempOutNames;
        outputTypes["shared"] = tempOutNames;
        outputTypes["relabund"] = tempOutNames;

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
        helpString += "The output parameter allows you to specify format of your biom file. Options hdf5 or simple. Default is hdf5, unless you are running a version without HDF5 libraries.\n";
		helpString += "The make.biom command should be in the following format: make.biom(shared=yourShared, groups=yourGroups, label=yourLabels).\n";
		helpString += "Example make.biom(shared=abrecovery.an.shared, groups=A-B-C).\n";
		helpString += "The default value for groups is all the groups in your groupfile, and all labels in your inputfile will be used.\n";
		helpString += "The make.biom command outputs a .biom file.\n";
		
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
MakeBiomCommand::MakeBiomCommand(string option) : Command() {
	try {
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
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
                if (sharedfile != "") {  inputFileName = sharedfile; fileFormat="sharedfile"; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n");  }
                else {
                    relabundfile = current->getRelAbundFile();
                    if (relabundfile != "") {  inputFileName = relabundfile; fileFormat="relabund"; m->mothurOut("Using " + relabundfile + " as input file for the relabund parameter.\n");  }
                    else {
                        m->mothurOut("No valid current files. You must provide a shared or relabund.\n");  abort = true;
                    }
                }
            }
            
			 
					if (outputdir == ""){    outputdir = util.hasPath(inputFileName);		}
            
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
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				else { allLines = true;  }
			}
			
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = ""; }
			else { 
				util.splitAtDash(groups, Groups);
                if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
			}
			
            if (picrustOtuFile != "") {
                picrust=true;
                if (contaxonomyfile == "") {  m->mothurOut("[ERROR]: the constaxonomy parameter is required with the picrust parameter, aborting.\n");  abort = true;  }
                if (referenceTax == "") {  m->mothurOut("[ERROR]: the reftaxonomy parameter is required with the picrust parameter, aborting.\n");  abort = true;  }
            }else { picrust=false; }
            
            if ((contaxonomyfile != "") && (labels.size() > 1)) { m->mothurOut("[ERROR]: the contaxonomy parameter cannot be used with multiple labels.\n");  abort = true; }
            
			format = validParameter.valid(parameters, "matrixtype");				if (format == "not found") { format = "sparse"; }
			
			if ((format != "sparse") && (format != "dense")) {
				m->mothurOut(format + " is not a valid option for the matrixtype parameter. Options are sparse and dense.\n");  abort = true; 
			}
            
            output = validParameter.valid(parameters, "output");
            if (output == "not found") {
                #ifdef USE_HDF5
                    output = "hdf5";
                #else
                    output = "simple";
                #endif
            }
            
            if ((output != "hdf5") && (output != "simple")) {
                 m->mothurOut("Invalid option for output. output options are hdf5 and simple, quitting.\n"); abort = true;
            }
            
            if (output == "hdf5") {
                #ifdef USE_HDF5
                //do nothing we have the api
                #else
                    m->mothurOut("[ERROR]: To write HDF5 biom files, you must have the API installed, quitting.\n"); abort=true;
                #endif
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
        
        SharedRAbundVectors* lookup = NULL; SharedRAbundFloatVectors* lookupRel = NULL;
        
		InputData input(inputFileName, fileFormat, Groups);
        set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        vector<string> sampleMetadata;
        
        if (fileFormat == "sharedfile") {
            lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
            Groups = lookup->getNamesGroups();
            sampleMetadata = getSampleMetaData(lookup);
        }else                        {
            lookupRel = util.getNextRelabund(input, allLines, userLabels, processedLabels, lastLabel);
            Groups = lookupRel->getNamesGroups();
            sampleMetadata = getSampleMetaData(lookupRel);
        }
        
        //if user did not specify a label, then use first one
        if ((contaxonomyfile != "") && (labels.size() == 0)) { allLines = false; labels.insert(lastLabel); }
        
        Picrust* piCrust = NULL;
        if (picrust) { piCrust = new Picrust(referenceTax, picrustOtuFile); }
        
        vector<Taxonomy> consTax;
        if (contaxonomyfile != "") { util.readConsTax(contaxonomyfile, consTax); }
        
        if (fileFormat == "sharedfile") {
            while (lookup != NULL) {
                
                if (m->getControl_pressed()) { delete lookup; break; }
                
                getBiom(lookup, piCrust, consTax, sampleMetadata); delete lookup;
                
                lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
            }
        }else {
            while (lookupRel != NULL) {
                            
                if (m->getControl_pressed()) { delete lookupRel; break; }
                
                getBiom(lookupRel, piCrust, consTax, sampleMetadata); delete lookupRel;
                            
                lookupRel = util.getNextRelabund(input, allLines, userLabels, processedLabels, lastLabel);
            }
        }
        
        if (picrust) { delete piCrust; }
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); }  return 0; }     
        
        //set sabund file as new current sabundfile
        string currentName = "";
		itTypes = outputTypes.find("biom");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setBiomFile(currentName); }
		}

		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeBiomCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
void MakeBiomCommand::getBiom(SharedRAbundVectors*& lookup, Picrust* picrust, vector<Taxonomy> consTax, vector<string> sampleMetadata){
	try {
        map<string, string> variables; 
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(sharedfile));
        variables["[distance]"] = lookup->getLabel();
        string outputFileName = getOutputFileName("biom",variables);
		outputNames.push_back(outputFileName); outputTypes["biom"].push_back(outputFileName);
        
        string mothurString = "mothur_" + toString(current->getVersion());

        Biom* biom;
        if (output == "hdf5")   { biom = new BiomHDF5();   }
        else                   { biom = new BiomSimple(); }
        
        biom->load(lookup, consTax);
        biom->fillHeading(mothurString, sharedfile);
        
        biom->print(outputFileName, sampleMetadata, picrust);
        
        if (picrust) {
            string outputFileName2 = getOutputFileName("shared",variables);
            outputNames.push_back(outputFileName2); outputTypes["shared"].push_back(outputFileName2);
            ofstream out2; util.openOutputFile(outputFileName2, out2);  bool printHeaders = true;
        
            biom->getSharedRAbundVectors()->print(out2, printHeaders);
        
            out2.close();
        }
       
        delete biom;
    }
	catch(exception& e) {
		m->errorOut(e, "MakeBiomCommand", "getBiom");
		exit(1);
	}
}
//**********************************************************************************************************************
void MakeBiomCommand::getBiom(SharedRAbundFloatVectors*& lookup, Picrust* picrust, vector<Taxonomy> consTax, vector<string> sampleMetadata){
    try {
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(inputFileName));
        variables["[distance]"] = lookup->getLabel();
        string outputFileName = getOutputFileName("biom",variables);
        outputNames.push_back(outputFileName); outputTypes["biom"].push_back(outputFileName);
        
        string mothurString = "mothur_" + toString(current->getVersion());
        
        BiomSimple biom; biom.load(lookup, consTax);
       
        biom.fillHeading(mothurString, sharedfile);
        biom.print(outputFileName, sampleMetadata, picrust);
        
        if (picrust) {
            string outputFileName2 = getOutputFileName("relabund",variables);
            outputNames.push_back(outputFileName2); outputTypes["relabund"].push_back(outputFileName2);
            ofstream out2; util.openOutputFile(outputFileName2, out2);  bool printHeaders = true;
        
            biom.getSharedRAbundFloatVectors()->print(out2, printHeaders);
        
            out2.close();
        }
    }
    catch(exception& e) {
        m->errorOut(e, "MakeBiomCommand", "getBiom");
        exit(1);
    }
}
//**********************************************************************************************************************
vector<string> MakeBiomCommand::getSampleMetaData(SharedRAbundVectors*& lookup){
	try {
        vector<string> sampleMetadata;
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
                
                if (m->getControl_pressed()) { break; }
                
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
                
                if (m->getControl_pressed()) { break; }
                
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
        
        return sampleMetadata;
        
    }
	catch(exception& e) {
		m->errorOut(e, "MakeBiomCommand", "getSampleMetaData");
		exit(1);
	}
    
}
//**********************************************************************************************************************
vector<string> MakeBiomCommand::getSampleMetaData(SharedRAbundFloatVectors*& lookup){
    try {
        vector<string> sampleMetadata;
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
                
                if (m->getControl_pressed()) { break; }
                
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
                
                if (m->getControl_pressed()) { break; }
                
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
        
        return sampleMetadata;
        
    }
    catch(exception& e) {
        m->errorOut(e, "MakeBiomCommand", "getSampleMetaData");
        exit(1);
    }
    
}
/**************************************************************************************************/
