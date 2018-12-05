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
		CommandParameter pshared("shared", "InputTypes", "", "", "SharedRelMeta", "SharedRelMeta", "none","otucorr",false,false,true); parameters.push_back(pshared);
		CommandParameter prelabund("relabund", "InputTypes", "", "", "SharedRelMeta", "SharedRelMeta", "none","otucorr",false,false); parameters.push_back(prelabund);
        CommandParameter pmetadata("metadata", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pmetadata);
        CommandParameter pcutoff("cutoff", "Number", "", "10", "", "", "","",false,false,true); parameters.push_back(pcutoff);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter pmethod("method", "Multiple", "pearson-spearman-kendall", "pearson", "", "", "","",false,false,true); parameters.push_back(pmethod);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
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
		helpString += "The otu.association command parameters are shared, relabund, metadata, groups, method, cutoff and label.  The shared or relabund parameter is required.\n";
		helpString += "The groups parameter allows you to specify which of the groups you would like included. The group names are separated by dashes.\n";
		helpString += "The label parameter allows you to select what distances level you would like used, and are also separated by dashes.\n";
        helpString += "The cutoff parameter allows you to set a pvalue at which the otu will be reported.\n";
		helpString += "The method parameter allows you to select what method you would like to use. Options are pearson, spearman and kendall. Default=pearson.\n";
		helpString += "The otu.association command should be in the following format: otu.association(shared=yourSharedFile, method=yourMethod).\n";
		helpString += "Example otu.association(shared=genus.pool.shared, method=kendall).\n";
		helpString += "The otu.association command outputs a .otu.corr file.\n";
		
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "OTUAssociationCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string OTUAssociationCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "otucorr") {  pattern = "[filename],[distance],[tag],otu.corr"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "OTUAssociationCommand", "getOutputPattern");
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
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
			vector<string> tempOutNames;
			outputTypes["otucorr"] = tempOutNames;
			
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
                
                it = parameters.find("metadata");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["metadata"] = inputDir + it->second;		}
				}
			}
			
			
			//check for required parameters			
			sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { abort = true; }
			else if (sharedfile == "not found") { sharedfile = ""; }
			else { inputFileName = sharedfile; current->setSharedFile(sharedfile); }
			
			relabundfile = validParameter.validFile(parameters, "relabund");
			if (relabundfile == "not open") { abort = true; }
			else if (relabundfile == "not found") { relabundfile = ""; }
			else { inputFileName = relabundfile; current->setRelAbundFile(relabundfile); }
			
            metadatafile = validParameter.valid(parameters, "metadata");
			if (metadatafile == "not open") { abort = true; metadatafile = ""; }
			else if (metadatafile == "not found") { metadatafile = ""; }
            
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = "";  pickedGroups = false;  }
			else { 
				pickedGroups = true;
				util.splitAtDash(groups, Groups);
                if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
			}
			
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = util.hasPath(inputFileName);	}
			
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			if ((relabundfile == "") && (sharedfile == "")) { 
				//is there are current file available for any of these?
				//give priority to shared, then relabund
				//if there is a current shared file, use it
				sharedfile = current->getSharedFile(); 
				if (sharedfile != "") { inputFileName = sharedfile; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 
					relabundfile = current->getRelAbundFile(); 
					if (relabundfile != "") { inputFileName = relabundfile;  m->mothurOut("Using " + relabundfile + " as input file for the relabund parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("You must provide either a shared or relabund file."); m->mothurOutEndLine(); abort = true; 
					}
				}
			}	
			
			
			if ((relabundfile != "") && (sharedfile != "")) { m->mothurOut("You may only use one of the following : shared or relabund file."); m->mothurOutEndLine(); abort = true;  }
			
			method = validParameter.valid(parameters, "method");		if (method == "not found"){	method = "pearson";		}
			
            string temp = validParameter.valid(parameters, "cutoff");
			if (temp == "not found") { temp = "10"; }
			util.mothurConvert(temp, cutoff); 
            
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
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        if (metadatafile != "") {  readMetadata(); }
		
		//function are identical just different datatypes
		if (sharedfile != "")			{  processShared();		} 
		else if (relabundfile != "")	{  processRelabund();	}
				
		if (m->getControl_pressed()) {  for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
		
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
		
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
		InputData* input = new InputData(sharedfile, "sharedfile", Groups);
		SharedRAbundVectors* lookup = input->getSharedRAbundVectors();
        Groups = lookup->getNamesGroups();
		string lastLabel = lookup->getLabel();
        
        if (metadatafile != "") {
            bool error = false;
            if (metadata[0].size() != lookup->size()) { m->mothurOut("[ERROR]: You have selected to use " + toString(metadata[0].size()) + " data rows from the metadata file, but " + toString(lookup->size()) + " from the shared file.\n");  m->setControl_pressed(true); error=true;}
            if (error) {
                //maybe add extra info here?? compare groups in each file??
            }
        }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (m->getControl_pressed()) {  delete input; return 0;  }
			
			if(allLines == 1 || labels.count(lookup->getLabel()) == 1){
				processedLabels.insert(lookup->getLabel()); userLabels.erase(lookup->getLabel());
				
				m->mothurOut(lookup->getLabel()+"\n"); 
				process(lookup);
			}
			
			if ((util.anyLabelsToProcess(lookup->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup->getLabel();
				
				delete lookup;
				lookup = input->getSharedRAbundVectors(lastLabel);
				
				processedLabels.insert(lookup->getLabel()); userLabels.erase(lookup->getLabel());
				
				//restore real lastlabel to save below
				lookup->setLabels(saveLabel);
				
				m->mothurOut(lookup->getLabel()+"\n"); 
				process(lookup);
			}
			
			lastLabel = lookup->getLabel();
			
			//get next line to process
			//prevent memory leak
			delete lookup;
			lookup = input->getSharedRAbundVectors();
		}
		
		
		if (m->getControl_pressed()) { delete input; return 0;  }
		
		//output error messages about any remaining user labels
		bool needToRun = false;
		for (set<string>::iterator it = userLabels.begin(); it != userLabels.end(); it++) {
			m->mothurOut("Your file does not include the label " + *it); 
            if (processedLabels.count(lastLabel) != 1)  { m->mothurOut(". I will use " + lastLabel + ".\n"); needToRun = true;  }
			else                                        { m->mothurOut(". Please refer to " + lastLabel + ".\n");               }
		}
		
		//run last label if you need to
		if (needToRun )  {
			delete lookup;
			lookup = input->getSharedRAbundVectors(lastLabel);
			
			m->mothurOut(lookup->getLabel()+"\n"); 
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
int OTUAssociationCommand::process(SharedRAbundVectors*& lookup){
	try {
		map<string, string> variables; 
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(inputFileName));
        variables["[distance]"] = lookup->getLabel();
        variables["[tag]"] = method;
		string outputFileName = getOutputFileName("otucorr",variables);
		outputNames.push_back(outputFileName); outputTypes["otucorr"].push_back(outputFileName);
		
		ofstream out;
		util.openOutputFile(outputFileName, out);
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		//column headings
		if (metadatafile == "") { out << "OTUA\tOTUB\t" << method << "Coef\tSignificance\n"; }
        else { out << "OTUA\tMetadata\t" << method << "Coef\tSignificance\n";  }

		
		vector< vector<double> > xy; xy.resize(lookup->getNumBins());
        vector<string> sampleNames = lookup->getNamesGroups();
		for (int i = 0; i < lookup->getNumBins(); i++) {
            vector<int> abunds = lookup->getOTU(i);
            for (int j = 0; j < abunds.size(); j++) { xy[i].push_back(abunds[j]); }
        }
		
		LinearAlgebra linear;
        vector<string> currentLabels = lookup->getOTUNames();
        if (metadatafile == "") {//compare otus
            for (int i = 0; i < xy.size(); i++) {
                
                for (int k = 0; k < i; k++) {
                    
                    if (m->getControl_pressed()) { out.close(); return 0; }
                    
                    double coef = 0.0;
                    double sig = 0.0;
                    if (method == "spearman")		{   coef = linear.calcSpearman(xy[i], xy[k], sig);	}
                    else if (method == "pearson")	{	coef = linear.calcPearson(xy[i], xy[k], sig);	}
                    else if (method == "kendall")	{	coef = linear.calcKendall(xy[i], xy[k], sig);	}                   
                    else { m->mothurOut("[ERROR]: invalid method, choices are spearman, pearson or kendall."); m->mothurOutEndLine(); m->setControl_pressed(true); }
                    
                    if (sig < cutoff) { out << currentLabels[i] << '\t' << currentLabels[k] << '\t' << coef << '\t' << sig << endl; }
                }
            }
		}else { //compare otus to metadata
            for (int i = 0; i < xy.size(); i++) {
                
                for (int k = 0; k < metadata.size(); k++) {
                    
                    if (m->getControl_pressed()) { out.close(); return 0; }
                    
                    double coef = 0.0;
                    double sig = 0.0;
                    if (method == "spearman")		{   coef = linear.calcSpearman(xy[i], metadata[k], sig);	}
                    else if (method == "pearson")	{	coef = linear.calcPearson(xy[i], metadata[k], sig);	}
                    else if (method == "kendall")	{	coef = linear.calcKendall(xy[i], metadata[k], sig);	}                   
                    else { m->mothurOut("[ERROR]: invalid method, choices are spearman, pearson or kendall."); m->mothurOutEndLine(); m->setControl_pressed(true); }
                    
                    if (sig < cutoff) { out << currentLabels[i] << '\t' << metadataLabels[k] << '\t' << coef << '\t' << sig << endl; }
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
		InputData* input = new InputData(relabundfile, "relabund", Groups);
		SharedRAbundFloatVectors* lookup = input->getSharedRAbundFloatVectors();
        Groups = lookup->getNamesGroups();
		string lastLabel = lookup->getLabel();
        
        if (metadatafile != "") {
            bool error = false;
            if (metadata[0].size() != lookup->size()) { m->mothurOut("[ERROR]: You have selected to use " + toString(metadata[0].size()) + " data rows from the metadata file, but " + toString(lookup->size()) + " from the relabund file.\n");  m->setControl_pressed(true); error=true;}
            if (error) {
                //maybe add extra info here?? compare groups in each file??
            }
        }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (m->getControl_pressed()) {  delete input; return 0;  }
			
			if(allLines == 1 || labels.count(lookup->getLabel()) == 1){
				processedLabels.insert(lookup->getLabel()); userLabels.erase(lookup->getLabel());
				
				m->mothurOut(lookup->getLabel()+"\n"); 
				process(lookup);
			}
			
			if ((util.anyLabelsToProcess(lookup->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup->getLabel();
				
				delete lookup;
				lookup = input->getSharedRAbundFloatVectors(lastLabel);
				
				processedLabels.insert(lookup->getLabel()); userLabels.erase(lookup->getLabel());
				
				//restore real lastlabel to save below
				lookup->setLabels(saveLabel);
				
				m->mothurOut(lookup->getLabel()+"\n"); 
				process(lookup);
			}
			
			lastLabel = lookup->getLabel();
			
			//get next line to process
			//prevent memory leak
			delete lookup;
			lookup = input->getSharedRAbundFloatVectors();
		}
		
		
		if (m->getControl_pressed()) { delete input; return 0;  }
		
		//output error messages about any remaining user labels
		bool needToRun = false;
		for (set<string>::iterator it = userLabels.begin(); it != userLabels.end(); it++) {
			m->mothurOut("Your file does not include the label " + *it); 
            if (processedLabels.count(lastLabel) != 1)  { m->mothurOut(". I will use " + lastLabel + ".\n"); needToRun = true;  }
			else                                        { m->mothurOut(". Please refer to " + lastLabel + ".\n");               }
		}
		
		//run last label if you need to
		if (needToRun )  {
			delete lookup;
			lookup = input->getSharedRAbundFloatVectors(lastLabel);
			
			m->mothurOut(lookup->getLabel()+"\n"); 
			process(lookup);
            delete lookup;
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
int OTUAssociationCommand::process(SharedRAbundFloatVectors*& lookup){
	try {
		
		map<string, string> variables; 
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(inputFileName));
        variables["[distance]"] = lookup->getLabel();
        variables["[tag]"] = method;
        string outputFileName = getOutputFileName("otucorr",variables);
		outputNames.push_back(outputFileName); outputTypes["otucorr"].push_back(outputFileName);
		
		ofstream out;
		util.openOutputFile(outputFileName, out);
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		//column headings
		if (metadatafile == "") { out << "OTUA\tOTUB\t" << method << "Coef\tSignificance\n"; }
        else { out << "OTUA\tMetadata\t" << method << "Coef\tSignificance\n";  }
		
		vector< vector<double> > xy; xy.resize(lookup->getNumBins());
        vector<string> sampleNames = lookup->getNamesGroups();
        for (int i = 0; i < lookup->getNumBins(); i++) {
            for (int j = 0; j < sampleNames.size(); j++) { xy[i].push_back(lookup->get(i, sampleNames[j])); }
        }
		
		LinearAlgebra linear;
        vector<string> currentLabels = lookup->getOTUNames();
        if (metadatafile == "") {//compare otus
            for (int i = 0; i < xy.size(); i++) {
                
                for (int k = 0; k < i; k++) {
                    
                    if (m->getControl_pressed()) { out.close(); return 0; }
                    
                    double coef = 0.0;
                    double sig = 0.0;
                    if (method == "spearman")		{   coef = linear.calcSpearman(xy[i], xy[k], sig);	}
                    else if (method == "pearson")	{	coef = linear.calcPearson(xy[i], xy[k], sig);	}
                    else if (method == "kendall")	{	coef = linear.calcKendall(xy[i], xy[k], sig);	}                   
                    else { m->mothurOut("[ERROR]: invalid method, choices are spearman, pearson or kendall."); m->mothurOutEndLine(); m->setControl_pressed(true); }
                    
                    if (sig < cutoff) { out << currentLabels[i] << '\t' << currentLabels[k] << '\t' << coef << '\t' << sig << endl; }
                }
            }
		}else { //compare otus to metadata
            for (int i = 0; i < xy.size(); i++) {
                
                for (int k = 0; k < metadata.size(); k++) {
                    
                    if (m->getControl_pressed()) { out.close(); return 0; }
                    
                    double coef = 0.0;
                    double sig = 0.0;
                    if (method == "spearman")		{   coef = linear.calcSpearman(xy[i], metadata[k], sig);	}
                    else if (method == "pearson")	{	coef = linear.calcPearson(xy[i], metadata[k], sig);	}
                    else if (method == "kendall")	{	coef = linear.calcKendall(xy[i], metadata[k], sig);	}                   
                    else { m->mothurOut("[ERROR]: invalid method, choices are spearman, pearson or kendall."); m->mothurOutEndLine(); m->setControl_pressed(true); }
                    
                    if (sig < cutoff) { out << currentLabels[i] << '\t' << metadataLabels[k] << '\t' << coef << '\t' << sig << endl; }
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
		util.openInputFile(metadatafile, in);
		
		string headerLine = util.getline(in); util.gobble(in);
		istringstream iss (headerLine,istringstream::in);
		
		//read the first label, because it refers to the groups
		string columnLabel;
		iss >> columnLabel; util.gobble(iss); 
		
		//save names of columns you are reading
		while (!iss.eof()) {
			iss >> columnLabel; util.gobble(iss);
            if (m->getDebug()) { m->mothurOut("[DEBUG]: metadata column Label = " + columnLabel + "\n"); }
			metadataLabels.push_back(columnLabel);
		}
        
		int count = metadataLabels.size();
        SharedRAbundFloatVectors* metadataLookup = new SharedRAbundFloatVectors();
        metadataLookup->setLabels("1");
        
		//read rest of file
		while (!in.eof()) {
			
			if (m->getControl_pressed()) { in.close(); return 0; }
			
			string group = "";
			in >> group; util.gobble(in);
            if (m->getDebug()) { m->mothurOut("[DEBUG]: metadata group = " + group + "\n"); }
            
            SharedRAbundFloatVector* tempLookup = new SharedRAbundFloatVector();
            tempLookup->setLabel("1");
            tempLookup->setGroup(group);
			
			for (int i = 0; i < count; i++) {
				float temp = 0.0;
				in >> temp;
                if (m->getDebug()) { m->mothurOut("[DEBUG]: metadata value = " + toString(temp) + "\n"); }
				tempLookup->push_back(temp);
			}
			
            if (util.inUsersGroups(group, Groups)) {  metadataLookup->push_back(tempLookup);  }
			
			util.gobble(in);
		}
		in.close();
        
        metadataLookup->setOTUNames(metadataLabels);
        metadataLookup->eliminateZeroOTUS();
        
        metadata.resize(metadataLookup->getNumBins());
        vector<string> sampleNames = metadataLookup->getNamesGroups();
        for (int i = 0; i < metadataLookup->getNumBins(); i++) {
            for (int j = 0; j < sampleNames.size(); j++) { metadata[i].push_back(metadataLookup->get(i, sampleNames[j])); }
        }
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "OTUAssociationCommand", "readMetadata");	
		exit(1);
	}
}
/*****************************************************************/









