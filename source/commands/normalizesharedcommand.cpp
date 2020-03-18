/*
 *  normalizesharedcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 9/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "normalizesharedcommand.h"

//**********************************************************************************************************************
vector<string> NormalizeSharedCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "LRSS", "LRSS", "none","shared",false,false,true); parameters.push_back(pshared);	
		CommandParameter prelabund("relabund", "InputTypes", "", "", "LRSS", "LRSS", "none","shared",false,false,true); parameters.push_back(prelabund);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter pmethod("method", "Multiple", "totalgroup-zscore", "totalgroup", "", "", "","",false,false,true); parameters.push_back(pmethod);
		CommandParameter pnorm("norm", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pnorm);
		CommandParameter pmakerelabund("makerelabund", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pmakerelabund);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "NormalizeSharedCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string NormalizeSharedCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The normalize.shared command parameters are shared, relabund, groups, method, norm, makerelabund and label.  shared or relabund is required, unless you have a valid current file.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like included. The group names are separated by dashes.\n";
		helpString += "The label parameter allows you to select what distance levels you would like, and are also separated by dashes.\n";
		helpString += "The method parameter allows you to select what method you would like to use to normalize. The options are totalgroup and zscore. We hope to add more ways to normalize in the future, suggestions are welcome!\n";
		helpString += "The makerelabund parameter allows you to convert a shared file to a relabund file before you normalize. default=f.\n";
		helpString += "The norm parameter allows you to number you would like to normalize to. By default this is set to the number of sequences in your smallest group.\n";
		helpString += "The normalize.shared command should be in the following format: normalize.shared(groups=yourGroups, label=yourLabels).\n";
		helpString += "Example normalize.shared(groups=A-B-C, scale=totalgroup).\n";
		helpString += "The default value for groups is all the groups in your groupfile, and all labels in your inputfile will be used.\n";
		helpString += "The normalize.shared command outputs a .norm.shared file.\n";
		
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "NormalizeSharedCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string NormalizeSharedCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "shared") {  pattern = "[filename],[distance],norm.shared"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "NormalizeSharedCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
NormalizeSharedCommand::NormalizeSharedCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["shared"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "NormalizeSharedCommand", "NormalizeSharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

NormalizeSharedCommand::NormalizeSharedCommand(string option) {
	try {
		abort = false; calledHelp = false;   
		allLines = true;
		
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
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["shared"] = tempOutNames;
			
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
			}
			
			sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }	
			else if (sharedfile == "not found") { sharedfile = ""; }
			else {  format = "sharedfile"; inputfile = sharedfile; current->setSharedFile(sharedfile); }
			
			relabundfile = validParameter.validFile(parameters, "relabund");
			if (relabundfile == "not open") { relabundfile = ""; abort = true; }	
			else if (relabundfile == "not found") { relabundfile = ""; }
			else {  format = "relabund"; inputfile = relabundfile; current->setRelAbundFile(relabundfile); }
			
			
			if ((sharedfile == "") && (relabundfile == "")) { 
				//is there are current file available for any of these?
				//give priority to shared, then list, then rabund, then sabund
				//if there is a current shared file, use it
				sharedfile = current->getSharedFile(); 
				if (sharedfile != "") { inputfile = sharedfile; format = "sharedfile"; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 
					relabundfile = current->getRelAbundFile(); 
					if (relabundfile != "") { inputfile = relabundfile; format = "relabund"; m->mothurOut("Using " + relabundfile + " as input file for the relabund parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("No valid current files. You must provide a list, sabund, rabund, relabund or shared file."); m->mothurOutEndLine(); 
						abort = true;
					}
				}
			}
			
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = util.hasPath(inputfile);		}
			
			

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				else { allLines = true;  }
			}
			
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = ""; pickedGroups = false; }
			else { 
				pickedGroups = true;
				util.splitAtDash(groups, Groups);
                if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
			}
			
			method = validParameter.valid(parameters, "method");				if (method == "not found") { method = "totalgroup"; }
			if ((method != "totalgroup") && (method != "zscore")) {  m->mothurOut(method + " is not a valid scaling option for the normalize.shared command. The options are totalgroup and zscore. We hope to add more ways to normalize in the future, suggestions are welcome!"); m->mothurOutEndLine(); abort = true; }
		
			string temp = validParameter.valid(parameters, "norm");
			if (temp == "not found") {  
				norm = 0;  //once you have read, set norm to smallest group number
			}else { 
				util.mothurConvert(temp, norm);
				if (norm < 0) { m->mothurOut("norm must be positive."); m->mothurOutEndLine(); abort=true; }
			}
			
			temp = validParameter.valid(parameters, "makerelabund");	if (temp == "") { temp = "f"; }
			makeRelabund = util.isTrue(temp);
		}

	}
	catch(exception& e) {
		m->errorOut(e, "NormalizeSharedCommand", "NormalizeSharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int NormalizeSharedCommand::execute(){
	try {
	
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		InputData input(inputfile, format, Groups);
        set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
		//you are reading a sharedfile and you do not want to make relabund
		if ((format == "sharedfile") && (!makeRelabund)) {
			SharedRAbundVectors* lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
			
			//look for groups whose numseqs is below norm and remove them, warning the user
			if (norm != 0) { lookup->removeGroups(norm); Groups = lookup->getNamesGroups(); }
			
			if (method == "totalgroup") {
				//set norm to smallest group number
				if (norm == 0) { norm = lookup->getNumSeqsSmallestGroup(); }
				
				m->mothurOut("Normalizing to " + toString(norm) + ".\n");
			}
			
			bool printHeaders = true;
            while (lookup != NULL) {
                
                if (m->getControl_pressed()) { delete lookup; break; }
                
                normalize(lookup, printHeaders); delete lookup;
                
                lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
            }
			
		}else{ //relabund values
			SharedRAbundFloatVectors* lookupFloat = util.getNextRelabund(input, allLines, userLabels, processedLabels, lastLabel);
			
			//look for groups whose numseqs is below norm and remove them, warning the user
			if (norm != 0) {  lookupFloat->removeGroups(norm); Groups = lookupFloat->getNamesGroups(); }
			
			//set norm to smallest group number
			if (method == "totalgroup") {
                if (norm == 0) {
                    norm = lookupFloat->getNumSeqsSmallestGroup();
                    Groups = lookupFloat->getNamesGroups();
                }
				
				m->mothurOut("Normalizing to " + toString(norm) + ".\n");
			}
			
            bool printHeaders = true;
            while (lookupFloat != NULL) {
                    
                if (m->getControl_pressed()) { delete lookupFloat; break; }

                normalize(lookupFloat, printHeaders); delete lookupFloat;
                    
                lookupFloat = util.getNextRelabund(input, allLines, userLabels, processedLabels, lastLabel);
            }
        }
        
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);	} outputTypes.clear(); return 0;}
		
		m->mothurOut("\nOutput File Names: \n"); 
		//m->mothurOut(outputFileName); m->mothurOutEndLine(); outputNames.push_back(outputFileName); outputTypes["shared"].push_back(outputFileName);
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
		
		//set shared file as new current sharedfile
		string currentName = "";
		itTypes = outputTypes.find("shared");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setSharedFile(currentName); }
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "NormalizeSharedCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int NormalizeSharedCommand::normalize(SharedRAbundVectors*& thisLookUp, bool& printHeaders){
	try {
        vector<string> lookupGroups = thisLookUp->getNamesGroups();
        map<string, string> variables; 
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(inputfile));
        variables["[distance]"] = thisLookUp->getLabel();
		string outputFileName = getOutputFileName("shared",variables);
        
		ofstream out;
		util.openOutputFile(outputFileName, out);
		outputNames.push_back(outputFileName); outputTypes["shared"].push_back(outputFileName);
				
		if (method == "totalgroup") { 
			
			//save numSeqs since they will change as the data is normalized
			vector<int> sizes;
			for (int i = 0; i < lookupGroups.size(); i++) {  sizes.push_back(thisLookUp->getNumSeqs(lookupGroups[i])); }
					
            for (int j = 0; j < thisLookUp->getNumBins(); j++) {
                
                for (int i = 0; i < lookupGroups.size(); i++) {
                    
                    if (m->getControl_pressed()) { out.close(); return 0; }
                    
                    int abund = thisLookUp->get(j, lookupGroups[i]);
                    
                    float relabund = abund / (float) sizes[i];
                    float newNorm = relabund * norm;
                    
                    //round to nearest int
                    int finalNorm = (int) floor((newNorm + 0.5));
                    
                    thisLookUp->set(j, finalNorm, lookupGroups[i]);
                }
            }
					
		}else if (method == "zscore") {
			
			for (int j = 0; j < thisLookUp->getNumBins(); j++) {
				
				if (m->getControl_pressed()) { out.close(); return 0; }
				
				//calc mean
				float mean = 0.0;
				for (int i = 0; i < lookupGroups.size(); i++) {  mean += thisLookUp->get(j, lookupGroups[i]); }
				mean /= (float) lookupGroups.size();
					
				//calc standard deviation
				float sumSquared = 0.0;
				for (int i = 0; i < lookupGroups.size(); i++) { sumSquared += (((float)thisLookUp->get(j, lookupGroups[i]) - mean) * ((float)thisLookUp->get(j, lookupGroups[i]) - mean)); }
				sumSquared /= (float) lookupGroups.size();
				
				float standardDev = sqrt(sumSquared);
					
				for (int i = 0; i < lookupGroups.size(); i++) {
					int finalNorm = 0;
					if (!util.isEqual(standardDev, 0)) { // stop divide by zero
						float newNorm = ((float)thisLookUp->get(j, lookupGroups[i]) - mean) / standardDev;
						//round to nearest int
						finalNorm = (int) floor((newNorm + 0.5));
					}
					
					thisLookUp->set(j, finalNorm, lookupGroups[i]);
				}
			}
						
		}else{ m->mothurOut(method + " is not a valid scaling option."); m->mothurOutEndLine(); m->setControl_pressed(true); return 0; }
				
		thisLookUp->eliminateZeroOTUS();
        thisLookUp->print(out, printHeaders);
		out.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "NormalizeSharedCommand", "normalize");
		exit(1);
	}
}
//**********************************************************************************************************************

int NormalizeSharedCommand::normalize(SharedRAbundFloatVectors*& thisLookUp, bool& printHeaders){
	try {
        vector<string> lookupGroups = thisLookUp->getNamesGroups();
        map<string, string> variables; 
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(inputfile));
        variables["[distance]"] = thisLookUp->getLabel();
		string outputFileName = getOutputFileName("shared",variables);
		ofstream out;
		util.openOutputFile(outputFileName, out);
		outputNames.push_back(outputFileName); outputTypes["shared"].push_back(outputFileName);
		
		if (method == "totalgroup") { 
			
			//save numSeqs since they will change as the data is normalized
			vector<float> sizes;
			for (int i = 0; i < lookupGroups.size(); i++) {  sizes.push_back(thisLookUp->getNumSeqs(lookupGroups[i])); }
			
			for (int j = 0; j < thisLookUp->getNumBins(); j++) {
				
				for (int i = 0; i < lookupGroups.size(); i++) {
					
					if (m->getControl_pressed()) { out.close(); return 0; }
					
					float abund = thisLookUp->get(j, lookupGroups[i]);
					
					float relabund = abund / (float) sizes[i];
					float newNorm = relabund * norm;
					
					thisLookUp->set(j, newNorm, lookupGroups[i]);
				}
			}
			
		}else if (method == "zscore") {
			for (int j = 0; j < thisLookUp->getNumBins(); j++) {
				
				if (m->getControl_pressed()) { out.close(); return 0; }
				
				//calc mean
				float mean = 0.0;
				for (int i = 0; i < lookupGroups.size(); i++) {  mean += thisLookUp->get(j, lookupGroups[i]); }
				mean /= (float) lookupGroups.size();
				
				//calc standard deviation
				float sumSquared = 0.0;
				for (int i = 0; i < lookupGroups.size(); i++) { sumSquared += ((thisLookUp->get(j, lookupGroups[i]) - mean) * (thisLookUp->get(j, lookupGroups[i]) - mean)); }
				sumSquared /= (float) lookupGroups.size();
				
				float standardDev = sqrt(sumSquared);
				
				for (int i = 0; i < lookupGroups.size(); i++) {
					float newNorm = 0.0;
                    if (!util.isEqual(standardDev, 0)) { // stop divide by zero
                        newNorm = ((float)thisLookUp->get(j, lookupGroups[i]) - mean) / standardDev;
                    }
                    thisLookUp->set(j, newNorm, lookupGroups[i]);
                }
			}			
			
		}else{ m->mothurOut(method + " is not a valid scaling option."); m->mothurOutEndLine(); m->setControl_pressed(true); return 0; }
		
		
        thisLookUp->eliminateZeroOTUS();
        thisLookUp->print(out, printHeaders);
        out.close();
	
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "NormalizeSharedCommand", "normalize");
		exit(1);
	}
}
//**********************************************************************************************************************
