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
        
        abort = false; calledHelp = false;    allLines = true;
		
        vector<string> tempOutNames;
        outputTypes["otucorr"] = tempOutNames;
        
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
OTUAssociationCommand::OTUAssociationCommand(string option) : Command()  {
	try {
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { abort = true; }
			else if (sharedfile == "not found") { sharedfile = ""; }
			else { inputFileName = sharedfile; current->setSharedFile(sharedfile); }
			
			relabundfile = validParameter.validFile(parameters, "relabund");
			if (relabundfile == "not open") { abort = true; }
			else if (relabundfile == "not found") { relabundfile = ""; }
			else { inputFileName = relabundfile; current->setRelAbundFile(relabundfile); }
			
            metadatafile = validParameter.validFile(parameters, "metadata");
			if (metadatafile == "not open") { abort = true; metadatafile = ""; }
			else if (metadatafile == "not found") { metadatafile = ""; }
            
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = "";  pickedGroups = false;  }
			else { 
				pickedGroups = true;
				util.splitAtDash(groups, Groups);
                if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
			}
			
					if (outputdir == ""){    outputdir = util.hasPath(inputFileName);	}
			
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				else { allLines = true;  }
			}
			
			if ((relabundfile == "") && (sharedfile == "")) { 
				//is there are current file available for any of these?
				//give priority to shared, then relabund
				//if there is a current shared file, use it
				sharedfile = current->getSharedFile(); 
				if (sharedfile != "") { inputFileName = sharedfile; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n");  }
				else { 
					relabundfile = current->getRelAbundFile(); 
					if (relabundfile != "") { inputFileName = relabundfile;  m->mothurOut("Using " + relabundfile + " as input file for the relabund parameter.\n");  }
					else { 
						m->mothurOut("You must provide either a shared or relabund file.\n");  abort = true; 
					}
				}
			}	
			
			
			if ((relabundfile != "") && (sharedfile != "")) { m->mothurOut("You may only use one of the following : shared or relabund file.\n");  abort = true;  }
			
			method = validParameter.valid(parameters, "method");		if (method == "not found"){	method = "pearson";		}
			
            string temp = validParameter.valid(parameters, "cutoff");
			if (temp == "not found") { temp = "10"; }
			util.mothurConvert(temp, cutoff); 
            
			if ((method != "pearson") && (method != "spearman") && (method != "kendall")) { m->mothurOut(method + " is not a valid method. Valid methods are pearson, spearman, and kendall.\n");  abort = true; }
			
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
void OTUAssociationCommand::processShared(){
	try {
		InputData input(sharedfile, "sharedfile", Groups);
        set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
        SharedRAbundVectors* lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
        Groups = lookup->getNamesGroups();
        
        if (metadatafile != "") {
            bool error = false;
            if (metadata[0].size() != lookup->size()) { m->mothurOut("[ERROR]: You have selected to use " + toString(metadata[0].size()) + " data rows from the metadata file, but " + toString(lookup->size()) + " from the shared file.\n");  m->setControl_pressed(true); error=true; }
        }
        
        while (lookup != nullptr) {
            
            if (m->getControl_pressed()) { delete lookup; break; }
            
            process(lookup); delete lookup;
            
            lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
        }
	}
	catch(exception& e) {
		m->errorOut(e, "OTUAssociationCommand", "processShared");	
		exit(1);
	}
}
//**********************************************************************************************************************
void OTUAssociationCommand::process(SharedRAbundVectors*& lookup){
	try {
		map<string, string> variables; 
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(inputFileName));
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
                    
                    if (m->getControl_pressed()) { out.close(); return; }
                    
                    double coef = 0.0;
                    double sig = 0.0;
                    if (method == "spearman")		{   coef = linear.calcSpearman(xy[i], xy[k], sig);	}
                    else if (method == "pearson")	{	coef = linear.calcPearson(xy[i], xy[k], sig);	}
                    else if (method == "kendall")	{	coef = linear.calcKendall(xy[i], xy[k], sig);	}                   
                    else { m->mothurOut("[ERROR]: invalid method, choices are spearman, pearson or kendall.\n");  m->setControl_pressed(true); }
                    
                    if (sig < cutoff) { out << currentLabels[i] << '\t' << currentLabels[k] << '\t' << coef << '\t' << sig << endl; }
                }
            }
		}else { //compare otus to metadata
            for (int i = 0; i < xy.size(); i++) {
                
                for (int k = 0; k < metadata.size(); k++) {
                    
                    if (m->getControl_pressed()) { out.close(); return; }
                    
                    double coef = 0.0;
                    double sig = 0.0;
                    if (method == "spearman")		{   coef = linear.calcSpearman(xy[i], metadata[k], sig);	}
                    else if (method == "pearson")	{	coef = linear.calcPearson(xy[i], metadata[k], sig);	}
                    else if (method == "kendall")	{	coef = linear.calcKendall(xy[i], metadata[k], sig);	}                   
                    else { m->mothurOut("[ERROR]: invalid method, choices are spearman, pearson or kendall.\n");  m->setControl_pressed(true); }
                    
                    if (sig < cutoff) { out << currentLabels[i] << '\t' << metadataLabels[k] << '\t' << coef << '\t' << sig << endl; }
                }
            }

        }
		out.close();
	}
	catch(exception& e) {
		m->errorOut(e, "OTUAssociationCommand", "process");	
		exit(1);
	}
}
//**********************************************************************************************************************
void OTUAssociationCommand::processRelabund(){
	try {
		InputData input(relabundfile, "relabund", Groups);
        set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
		SharedRAbundFloatVectors* lookup = util.getNextRelabund(input, allLines, userLabels, processedLabels, lastLabel);
        Groups = lookup->getNamesGroups();
        
        if (metadatafile != "") {
            bool error = false;
            if (metadata[0].size() != lookup->size()) { m->mothurOut("[ERROR]: You have selected to use " + toString(metadata[0].size()) + " data rows from the metadata file, but " + toString(lookup->size()) + " from the relabund file.\n");  m->setControl_pressed(true); error=true;}
        }
		
        while (lookup != nullptr) {
            
            if (m->getControl_pressed()) { delete lookup; break; }
            
            process(lookup); delete lookup;
            
            lookup = util.getNextRelabund(input, allLines, userLabels, processedLabels, lastLabel);
        }
	}
	catch(exception& e) {
		m->errorOut(e, "OTUAssociationCommand", "processRelabund");	
		exit(1);
	}
}
//**********************************************************************************************************************
void OTUAssociationCommand::process(SharedRAbundFloatVectors*& lookup){
	try {
		
		map<string, string> variables; 
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(inputFileName));
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
                    
                    if (m->getControl_pressed()) { out.close(); return; }
                    
                    double coef = 0.0;
                    double sig = 0.0;
                    if (method == "spearman")		{   coef = linear.calcSpearman(xy[i], xy[k], sig);	}
                    else if (method == "pearson")	{	coef = linear.calcPearson(xy[i], xy[k], sig);	}
                    else if (method == "kendall")	{	coef = linear.calcKendall(xy[i], xy[k], sig);	}                   
                    else { m->mothurOut("[ERROR]: invalid method, choices are spearman, pearson or kendall.\n");  m->setControl_pressed(true); }
                    
                    if (sig < cutoff) { out << currentLabels[i] << '\t' << currentLabels[k] << '\t' << coef << '\t' << sig << endl; }
                }
            }
		}else { //compare otus to metadata
            for (int i = 0; i < xy.size(); i++) {
                
                for (int k = 0; k < metadata.size(); k++) {
                    
                    if (m->getControl_pressed()) { out.close(); return; }
                    
                    double coef = 0.0;
                    double sig = 0.0;
                    if (method == "spearman")		{   coef = linear.calcSpearman(xy[i], metadata[k], sig);	}
                    else if (method == "pearson")	{	coef = linear.calcPearson(xy[i], metadata[k], sig);	}
                    else if (method == "kendall")	{	coef = linear.calcKendall(xy[i], metadata[k], sig);	}                   
                    else { m->mothurOut("[ERROR]: invalid method, choices are spearman, pearson or kendall.\n");  m->setControl_pressed(true); }
                    
                    if (sig < cutoff) { out << currentLabels[i] << '\t' << metadataLabels[k] << '\t' << coef << '\t' << sig << endl; }
                }
            }
            
        }
		out.close();
	}
	catch(exception& e) {
		m->errorOut(e, "OTUAssociationCommand", "process");	
		exit(1);
	}
}
/*****************************************************************/
void OTUAssociationCommand::readMetadata(){
	try {
		ifstream in; util.openInputFile(metadatafile, in);
		
		string headerLine = util.getline(in); gobble(in);
        metadataLabels = util.splitWhiteSpace(headerLine);
        metadataLabels.erase(metadataLabels.begin());
        
		int count = metadataLabels.size();
        SharedRAbundFloatVectors* metadataLookup = new SharedRAbundFloatVectors();
        metadataLookup->setLabels("1");
        
		//read rest of file
		while (!in.eof()) {
			
			if (m->getControl_pressed()) { in.close(); return; }
			
			string group = "";
			in >> group; gobble(in);
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
			
            if (Groups.size() == 0) { metadataLookup->push_back(tempLookup);  }
            else if (util.inUsersGroups(group, Groups)) {  metadataLookup->push_back(tempLookup);  }
			
			gobble(in);
		}
		in.close();
        
        metadataLookup->setOTUNames(metadataLabels);
        metadataLookup->eliminateZeroOTUS();
        
        metadata.resize(metadataLookup->getNumBins());
        vector<string> sampleNames = metadataLookup->getNamesGroups();
        for (int i = 0; i < metadataLookup->getNumBins(); i++) {
            for (int j = 0; j < sampleNames.size(); j++) { metadata[i].push_back(metadataLookup->get(i, sampleNames[j])); }
        }
        delete metadataLookup;
	}
	catch(exception& e) {
		m->errorOut(e, "OTUAssociationCommand", "readMetadata");	
		exit(1);
	}
}
/*****************************************************************/









