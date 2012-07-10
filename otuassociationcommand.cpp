/*
 *  otuassociationcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 1/19/12.
 *  Copyright 2012 Schloss Lab. All rights reserved.
 *
 */

#include "otuassociationcommand.h"
#include "linearalgebra.h"

//**********************************************************************************************************************
vector<string> OTUAssociationCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "SharedRelMeta", "SharedRelMeta", "none",false,false); parameters.push_back(pshared);
		CommandParameter prelabund("relabund", "InputTypes", "", "", "SharedRelMeta", "SharedRelMeta", "none",false,false); parameters.push_back(prelabund);
        CommandParameter pmetadata("metadata", "InputTypes", "", "", "SharedRelMeta", "SharedRelMeta", "none",false,false); parameters.push_back(pmetadata);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter pmethod("method", "Multiple", "pearson-spearman-kendall", "pearson", "", "", "",false,false); parameters.push_back(pmethod);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "OTUAssociationCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string OTUAssociationCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The otu.association command reads a shared or relabund file and calculates the correlation coefficients between otus.\n";
        helpString += "If you provide a metadata file, mothur will calculate te correlation bewteen the metadata and the otus.\n";
		helpString += "The otu.association command parameters are shared, relabund, metadata, groups, method and label.  The shared or relabund parameter is required.\n";
		helpString += "The groups parameter allows you to specify which of the groups you would like included. The group names are separated by dashes.\n";
		helpString += "The label parameter allows you to select what distances level you would like used, and are also separated by dashes.\n";
		helpString += "The method parameter allows you to select what method you would like to use. Options are pearson, spearman and kendall. Default=pearson.\n";
		helpString += "The otu.association command should be in the following format: otu.association(shared=yourSharedFile, method=yourMethod).\n";
		helpString += "Example otu.association(shared=genus.pool.shared, method=kendall).\n";
		helpString += "The otu.association command outputs a .otu.corr file.\n";
		helpString += "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "OTUAssociationCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string OTUAssociationCommand::getOutputFileNameTag(string type, string inputName=""){	
	try {
        string outputFileName = "";
		map<string, vector<string> >::iterator it;
        
        //is this a type this command creates
        it = outputTypes.find(type);
        if (it == outputTypes.end()) {  m->mothurOut("[ERROR]: this command doesn't create a " + type + " output file.\n"); }
        else {
            if (type == "otucorr") {  outputFileName =  "otu.corr"; }
            else { m->mothurOut("[ERROR]: No definition for type " + type + " output file tag.\n"); m->control_pressed = true;  }
        }
        return outputFileName;
	}
	catch(exception& e) {
		m->errorOut(e, "OTUAssociationCommand", "getOutputFileNameTag");
		exit(1);
	}
}
//**********************************************************************************************************************
OTUAssociationCommand::OTUAssociationCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["otucorr"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "OTUAssociationCommand", "OTUAssociationCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
OTUAssociationCommand::OTUAssociationCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			vector<string> tempOutNames;
			outputTypes["otucorr"] = tempOutNames;
			
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
				
				it = parameters.find("relabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["relabund"] = inputDir + it->second;		}
				}
                
                it = parameters.find("metadata");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["metadata"] = inputDir + it->second;		}
				}
			}
			
			
			//check for required parameters			
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { abort = true; }
			else if (sharedfile == "not found") { sharedfile = ""; }
			else { inputFileName = sharedfile; m->setSharedFile(sharedfile); }
			
			relabundfile = validParameter.validFile(parameters, "relabund", true);
			if (relabundfile == "not open") { abort = true; }
			else if (relabundfile == "not found") { relabundfile = ""; }
			else { inputFileName = relabundfile; m->setRelAbundFile(relabundfile); }
			
            metadatafile = validParameter.validFile(parameters, "metadata", true);
			if (metadatafile == "not open") { abort = true; metadatafile = ""; }
			else if (metadatafile == "not found") { metadatafile = ""; }
            
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = "";  pickedGroups = false;  }
			else { 
				pickedGroups = true;
				m->splitAtDash(groups, Groups);	
			}			
			m->setGroups(Groups);
			
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(inputFileName);	}
			
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			if ((relabundfile == "") && (sharedfile == "")) { 
				//is there are current file available for any of these?
				//give priority to shared, then relabund
				//if there is a current shared file, use it
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") { inputFileName = sharedfile; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 
					relabundfile = m->getRelAbundFile(); 
					if (relabundfile != "") { inputFileName = relabundfile;  m->mothurOut("Using " + relabundfile + " as input file for the relabund parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("You must provide either a shared or relabund file."); m->mothurOutEndLine(); abort = true; 
					}
				}
			}	
			
			
			if ((relabundfile != "") && (sharedfile != "")) { m->mothurOut("You may only use one of the following : shared or relabund file."); m->mothurOutEndLine(); abort = true;  }
			
			method = validParameter.validFile(parameters, "method", false);		if (method == "not found"){	method = "pearson";		}
			
			if ((method != "pearson") && (method != "spearman") && (method != "kendall")) { m->mothurOut(method + " is not a valid method. Valid methods are pearson, spearman, and kendall."); m->mothurOutEndLine(); abort = true; }
			
		}
	}
	catch(exception& e) {
		m->errorOut(e, "OTUAssociationCommand", "OTUAssociationCommand");		
		exit(1);
	}
}
//**********************************************************************************************************************

int OTUAssociationCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        if (metadatafile != "") {  readMetadata(); }
		
		//function are identical just different datatypes
		if (sharedfile != "")			{  processShared();		} 
		else if (relabundfile != "")	{  processRelabund();	}
				
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "OTUAssociationCommand", "execute");	
		exit(1);
	}
}
//**********************************************************************************************************************
int OTUAssociationCommand::processShared(){
	try {
		InputData* input = new InputData(sharedfile, "sharedfile");
		vector<SharedRAbundVector*> lookup = input->getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
        
        if (metadatafile != "") { 
            getMetadata();  
            bool error = false;
            if (metadata[0].size() != lookup.size()) { m->mothurOut("[ERROR]: You have selected to use " + toString(metadata[0].size()) + " data rows from the metadata file, but " + toString(lookup.size()) + " from the shared file.\n");  m->control_pressed = true; error=true;}
            if (error) {
                //maybe add extra info here?? compare groups in each file??
            }
        }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (m->control_pressed) {  delete input; return 0;  }
			
			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){	
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				process(lookup);
			}
			
			if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup[0]->getLabel();
				
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
				lookup = input->getSharedRAbundVectors(lastLabel);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				
				//restore real lastlabel to save below
				lookup[0]->setLabel(saveLabel);
				
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				process(lookup);
			}
			
			lastLabel = lookup[0]->getLabel();			
			
			//get next line to process
			//prevent memory leak
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
			lookup = input->getSharedRAbundVectors();
		}
		
		
		if (m->control_pressed) { delete input; return 0;  }
		
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
			for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) {	delete lookup[i];	} } 
			lookup = input->getSharedRAbundVectors(lastLabel);
			
			m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
			process(lookup);
		}	
		
		delete input;
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "OTUAssociationCommand", "processShared");	
		exit(1);
	}
}
//**********************************************************************************************************************
int OTUAssociationCommand::process(vector<SharedRAbundVector*>& lookup){
	try {
		
		string outputFileName = outputDir + m->getRootName(m->getSimpleName(inputFileName)) + lookup[0]->getLabel() + "." + method + "." + getOutputFileNameTag("otucorr");
		outputNames.push_back(outputFileName); outputTypes["otucorr"].push_back(outputFileName);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		//column headings
		if (metadatafile == "") { out << "OTUA\tOTUB\t" << method << "Coef\tSignificance\n"; }
        else { out << "OTUA\tMetadata\t" << method << "Coef\tSignificance\n";  }

		
		vector< vector<double> > xy; xy.resize(lookup[0]->getNumBins());
		for (int i = 0; i < lookup[0]->getNumBins(); i++) { for (int j = 0; j < lookup.size(); j++) { xy[i].push_back(lookup[j]->getAbundance(i)); } }
		
		LinearAlgebra linear;
        if (metadatafile == "") {//compare otus
            for (int i = 0; i < xy.size(); i++) {
                
                for (int k = 0; k < i; k++) {
                    
                    if (m->control_pressed) { out.close(); return 0; }
                    
                    double coef = 0.0;
                    double sig = 0.0;
                    if (method == "spearman")		{   coef = linear.calcSpearman(xy[i], xy[k], sig);	}
                    else if (method == "pearson")	{	coef = linear.calcPearson(xy[i], xy[k], sig);	}
                    else if (method == "kendall")	{	coef = linear.calcKendall(xy[i], xy[k], sig);	}                   
                    else { m->mothurOut("[ERROR]: invalid method, choices are spearman, pearson or kendall."); m->mothurOutEndLine(); m->control_pressed = true; }
                    
                    out << m->binLabelsInFile[i] << '\t' << m->binLabelsInFile[k] << '\t' << coef << '\t' << sig << endl;
                }
            }
		}else { //compare otus to metadata
            for (int i = 0; i < xy.size(); i++) {
                
                for (int k = 0; k < metadata.size(); k++) {
                    
                    if (m->control_pressed) { out.close(); return 0; }
                    
                    double coef = 0.0;
                    double sig = 0.0;
                    if (method == "spearman")		{   coef = linear.calcSpearman(xy[i], metadata[k], sig);	}
                    else if (method == "pearson")	{	coef = linear.calcPearson(xy[i], metadata[k], sig);	}
                    else if (method == "kendall")	{	coef = linear.calcKendall(xy[i], metadata[k], sig);	}                   
                    else { m->mothurOut("[ERROR]: invalid method, choices are spearman, pearson or kendall."); m->mothurOutEndLine(); m->control_pressed = true; }
                    
                    out << m->binLabelsInFile[i] << '\t' << metadataLabels[k] << '\t' << coef << '\t' << sig << endl;
                }
            }

        }
		out.close();
		
               
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "OTUAssociationCommand", "process");	
		exit(1);
	}
}
//**********************************************************************************************************************
int OTUAssociationCommand::processRelabund(){
	try {
		InputData* input = new InputData(relabundfile, "relabund");
		vector<SharedRAbundFloatVector*> lookup = input->getSharedRAbundFloatVectors();
		string lastLabel = lookup[0]->getLabel();
        
        if (metadatafile != "") { 
            getMetadata(); 
            bool error = false;
            if (metadata[0].size() != lookup.size()) { m->mothurOut("[ERROR]: You have selected to use " + toString(metadata[0].size()) + " data rows from the metadata file, but " + toString(lookup.size()) + " from the relabund file.\n");  m->control_pressed = true; error=true;}
            if (error) {
                //maybe add extra info here?? compare groups in each file??
            }
        }
        
        
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (m->control_pressed) {  delete input; return 0;  }
			
			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){	
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				process(lookup);
			}
			
			if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup[0]->getLabel();
				
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
				lookup = input->getSharedRAbundFloatVectors(lastLabel);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				
				//restore real lastlabel to save below
				lookup[0]->setLabel(saveLabel);
				
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				process(lookup);
			}
			
			lastLabel = lookup[0]->getLabel();			
			
			//get next line to process
			//prevent memory leak
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
			lookup = input->getSharedRAbundFloatVectors();
		}
		
		
		if (m->control_pressed) { delete input; return 0;  }
		
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
			for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) {	delete lookup[i];	} } 
			lookup = input->getSharedRAbundFloatVectors(lastLabel);
			
			m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
			process(lookup);
		}	
		
		delete input;
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "OTUAssociationCommand", "processRelabund");	
		exit(1);
	}
}
//**********************************************************************************************************************
int OTUAssociationCommand::process(vector<SharedRAbundFloatVector*>& lookup){
	try {
		
		string outputFileName = outputDir + m->getRootName(m->getSimpleName(inputFileName)) + lookup[0]->getLabel() + "." + method + "." + getOutputFileNameTag("otucorr");
		outputNames.push_back(outputFileName); outputTypes["otucorr"].push_back(outputFileName);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		//column headings
		if (metadatafile == "") { out << "OTUA\tOTUB\t" << method << "Coef\tSignificance\n"; }
        else { out << "OTUA\tMetadata\t" << method << "Coef\tSignificance\n";  }
		
		vector< vector<double> > xy; xy.resize(lookup[0]->getNumBins());
		for (int i = 0; i < lookup[0]->getNumBins(); i++) { for (int j = 0; j < lookup.size(); j++) { xy[i].push_back(lookup[j]->getAbundance(i)); } }
		
		LinearAlgebra linear;
        if (metadatafile == "") {//compare otus
            for (int i = 0; i < xy.size(); i++) {
                
                for (int k = 0; k < i; k++) {
                    
                    if (m->control_pressed) { out.close(); return 0; }
                    
                    double coef = 0.0;
                    double sig = 0.0;
                    if (method == "spearman")		{   coef = linear.calcSpearman(xy[i], xy[k], sig);	}
                    else if (method == "pearson")	{	coef = linear.calcPearson(xy[i], xy[k], sig);	}
                    else if (method == "kendall")	{	coef = linear.calcKendall(xy[i], xy[k], sig);	}                   
                    else { m->mothurOut("[ERROR]: invalid method, choices are spearman, pearson or kendall."); m->mothurOutEndLine(); m->control_pressed = true; }
                    
                    out << m->binLabelsInFile[i] << '\t' << m->binLabelsInFile[k] << '\t' << coef << '\t' << sig << endl;
                }
            }
		}else { //compare otus to metadata
            for (int i = 0; i < xy.size(); i++) {
                
                for (int k = 0; k < metadata.size(); k++) {
                    
                    if (m->control_pressed) { out.close(); return 0; }
                    
                    double coef = 0.0;
                    double sig = 0.0;
                    if (method == "spearman")		{   coef = linear.calcSpearman(xy[i], metadata[k], sig);	}
                    else if (method == "pearson")	{	coef = linear.calcPearson(xy[i], metadata[k], sig);	}
                    else if (method == "kendall")	{	coef = linear.calcKendall(xy[i], metadata[k], sig);	}                   
                    else { m->mothurOut("[ERROR]: invalid method, choices are spearman, pearson or kendall."); m->mothurOutEndLine(); m->control_pressed = true; }
                    
                    out << m->binLabelsInFile[i] << '\t' << metadataLabels[k] << '\t' << coef << '\t' << sig << endl;
                }
            }
            
        }
		
		out.close();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "OTUAssociationCommand", "process");	
		exit(1);
	}
}
/*****************************************************************/
int OTUAssociationCommand::readMetadata(){
	try {
		ifstream in;
		m->openInputFile(metadatafile, in);
		
		string headerLine = m->getline(in); m->gobble(in);
		istringstream iss (headerLine,istringstream::in);
		
		//read the first label, because it refers to the groups
		string columnLabel;
		iss >> columnLabel; m->gobble(iss); 
		
		//save names of columns you are reading
		while (!iss.eof()) {
			iss >> columnLabel; m->gobble(iss);
			metadataLabels.push_back(columnLabel);
		}
		int count = metadataLabels.size();
        
		//read rest of file
		while (!in.eof()) {
			
			if (m->control_pressed) { in.close(); return 0; }
			
			string group = "";
			in >> group; m->gobble(in);
            
			SharedRAbundFloatVector* tempLookup = new SharedRAbundFloatVector();
			tempLookup->setGroup(group);
			tempLookup->setLabel("1");
			
			for (int i = 0; i < count; i++) {
				float temp = 0.0;
				in >> temp; 
				tempLookup->push_back(temp, group);
			}
			
			metadataLookup.push_back(tempLookup);
			
			m->gobble(in);
		}
		in.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "OTUAssociationCommand", "readMetadata");	
		exit(1);
	}
}
/*****************************************************************/
//eliminate groups user did not pick, remove zeroed out otus, fill metadata vector.
int OTUAssociationCommand::getMetadata(){
	try {
        
		vector<string> mGroups = m->getGroups();
        
		bool remove = false;
		for (int i = 0; i < metadataLookup.size(); i++) {
			//if this sharedrabund is not from a group the user wants then delete it.
			if (!(m->inUsersGroups(metadataLookup[i]->getGroup(), mGroups))) { 
				delete metadataLookup[i]; metadataLookup[i] = NULL;
				metadataLookup.erase(metadataLookup.begin()+i); 
				i--; 
				remove = true;
			}
		}
        
        vector<SharedRAbundFloatVector*> newLookup;
		for (int i = 0; i < metadataLookup.size(); i++) {
			SharedRAbundFloatVector* temp = new SharedRAbundFloatVector();
			temp->setLabel(metadataLookup[i]->getLabel());
			temp->setGroup(metadataLookup[i]->getGroup());
			newLookup.push_back(temp);
		}
		
		//for each bin
        vector<string> newBinLabels;
		for (int i = 0; i < metadataLookup[0]->getNumBins(); i++) {
			if (m->control_pressed) { for (int j = 0; j < newLookup.size(); j++) {  delete newLookup[j];  } return 0; }
			
			//look at each sharedRabund and make sure they are not all zero
			bool allZero = true;
			for (int j = 0; j < metadataLookup.size(); j++) {
				if (metadataLookup[j]->getAbundance(i) != 0) { allZero = false;  break;  }
			}
			
			//if they are not all zero add this bin
			if (!allZero) {
				for (int j = 0; j < metadataLookup.size(); j++) {
					newLookup[j]->push_back(metadataLookup[j]->getAbundance(i), metadataLookup[j]->getGroup());
				}
                newBinLabels.push_back(metadataLabels[i]);
			}
		}
		
        metadataLabels = newBinLabels;
        
		for (int j = 0; j < metadataLookup.size(); j++) {  delete metadataLookup[j];  } 
        metadataLookup.clear();
		
        metadata.resize(newLookup[0]->getNumBins());
		for (int i = 0; i < newLookup[0]->getNumBins(); i++) { for (int j = 0; j < newLookup.size(); j++) { metadata[i].push_back(newLookup[j]->getAbundance(i)); } }
        
        for (int j = 0; j < newLookup.size(); j++) {  delete newLookup[j];  }
	        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "OTUAssociationCommand", "getMetadata");	
		exit(1);
	}
}
/*****************************************************************/








