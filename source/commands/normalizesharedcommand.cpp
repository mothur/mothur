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
		helpString += "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n";
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
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
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
				
				it = parameters.find("relabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["relabund"] = inputDir + it->second;		}
				}
			}
			
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }	
			else if (sharedfile == "not found") { sharedfile = ""; }
			else {  format = "sharedfile"; inputfile = sharedfile; m->setSharedFile(sharedfile); }
			
			relabundfile = validParameter.validFile(parameters, "relabund", true);
			if (relabundfile == "not open") { relabundfile = ""; abort = true; }	
			else if (relabundfile == "not found") { relabundfile = ""; }
			else {  format = "relabund"; inputfile = relabundfile; m->setRelAbundFile(relabundfile); }
			
			
			if ((sharedfile == "") && (relabundfile == "")) { 
				//is there are current file available for any of these?
				//give priority to shared, then list, then rabund, then sabund
				//if there is a current shared file, use it
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") { inputfile = sharedfile; format = "sharedfile"; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 
					relabundfile = m->getRelAbundFile(); 
					if (relabundfile != "") { inputfile = relabundfile; format = "relabund"; m->mothurOut("Using " + relabundfile + " as input file for the relabund parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("No valid current files. You must provide a list, sabund, rabund, relabund or shared file."); m->mothurOutEndLine(); 
						abort = true;
					}
				}
			}
			
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(inputfile);		}
			
			

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; pickedGroups = false; }
			else { 
				pickedGroups = true;
				m->splitAtDash(groups, Groups);
				m->setGroups(Groups);
			}
			
			method = validParameter.validFile(parameters, "method", false);				if (method == "not found") { method = "totalgroup"; }
			if ((method != "totalgroup") && (method != "zscore")) {  m->mothurOut(method + " is not a valid scaling option for the normalize.shared command. The options are totalgroup and zscore. We hope to add more ways to normalize in the future, suggestions are welcome!"); m->mothurOutEndLine(); abort = true; }
		
			string temp = validParameter.validFile(parameters, "norm", false);				
			if (temp == "not found") {  
				norm = 0;  //once you have read, set norm to smallest group number
			}else { 
				m->mothurConvert(temp, norm);
				if (norm < 0) { m->mothurOut("norm must be positive."); m->mothurOutEndLine(); abort=true; }
			}
			
			temp = validParameter.validFile(parameters, "makerelabund", false);	if (temp == "") { temp = "f"; }
			makeRelabund = m->isTrue(temp);
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
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		InputData input(inputfile, format);
        
		//you are reading a sharedfile and you do not want to make relabund
		if ((format == "sharedfile") && (!makeRelabund)) {
			SharedRAbundVectors* lookup = input.getSharedRAbundVectors();
			string lastLabel = lookup->getLabel();
			
			//look for groups whose numseqs is below norm and remove them, warning the user
			if (norm != 0) { 
				m->clearGroups();
                lookup->removeGroups(norm);
                vector<string> mGroups = lookup->getNamesGroups();
				m->setGroups(mGroups);
			}
			
			//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
			set<string> processedLabels;
			set<string> userLabels = labels;
			
			if (method == "totalgroup") {
				//set norm to smallest group number
				if (norm == 0) { norm = lookup->getNumSeqsSmallestGroup(); }
				
				m->mothurOut("Normalizing to " + toString(norm) + "."); m->mothurOutEndLine();
			}
			
			
			//as long as you are not at the end of the file or done wih the lines you want
			while((lookup != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
				
                if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);	} outputTypes.clear();  delete lookup; m->clearGroups();   return 0; }
				
				if(allLines == 1 || labels.count(lookup->getLabel()) == 1){
					
					m->mothurOut(lookup->getLabel()); m->mothurOutEndLine();
					normalize(lookup);
					
					processedLabels.insert(lookup->getLabel());
					userLabels.erase(lookup->getLabel());
				}
				
				if ((m->anyLabelsToProcess(lookup->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = lookup->getLabel();
					
					delete lookup;
					lookup = input.getSharedRAbundVectors(lastLabel);
					m->mothurOut(lookup->getLabel()); m->mothurOutEndLine();
					
					normalize(lookup);
					
					processedLabels.insert(lookup->getLabel());
					userLabels.erase(lookup->getLabel());
					
					//restore real lastlabel to save below
					lookup->setLabels(saveLabel);
				}
				
				lastLabel = lookup->getLabel();
				//prevent memory leak
				delete lookup;
				
				if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);	} outputTypes.clear(); m->clearGroups();  return 0; }
				
				//get next line to process
				lookup = input.getSharedRAbundVectors();
			}
			
			if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);	} outputTypes.clear(); m->clearGroups();   return 0; }
			
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
				delete lookup;
				lookup = input.getSharedRAbundVectors(lastLabel);
				m->mothurOut(lookup->getLabel()); m->mothurOutEndLine();
				normalize(lookup);
				delete lookup;
			}
			
		}else{ //relabund values
			SharedRAbundFloatVectors* lookupFloat = input.getSharedRAbundFloatVectors();
			string lastLabel = lookupFloat->getLabel();
			
			//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
			set<string> processedLabels;
			set<string> userLabels = labels;
			
			//look for groups whose numseqs is below norm and remove them, warning the user
			if (norm != 0) { 
                m->clearGroups();
                vector<string> mGroups;
                vector<string> temp;
                vector<string> lookupGroups = lookupFloat->getNamesGroups();
                for (int i = 0; i < lookupGroups.size(); i++) {
                    if (lookupFloat->getNumSeqs(lookupGroups[i]) < norm) {
                        m->mothurOut(lookupGroups[i] + " contains " + toString(lookupFloat->getNumSeqs(lookupGroups[i])) + ". Eliminating."); m->mothurOutEndLine();
                        temp.push_back(lookupGroups[i]);
                    }else { Groups.push_back(lookupGroups[i]); }
                }
                lookupFloat->removeGroups(temp);
                m->setGroups(mGroups);
			}
			
			//set norm to smallest group number
			if (method == "totalgroup") {
                if (norm == 0) {
                    vector<string> lookupGroups = lookupFloat->getNamesGroups();
                    norm = lookupFloat->getNumSeqs(lookupGroups[0]);
                    for (int i = 1; i < lookupGroups.size(); i++) {
                        if (lookupFloat->getNumSeqs(lookupGroups[i]) < norm) { norm = lookupFloat->getNumSeqs(lookupGroups[i]);  }
                    }  
                }
				
				m->mothurOut("Normalizing to " + toString(norm) + "."); m->mothurOutEndLine();
			}
			
			//as long as you are not at the end of the file or done wih the lines you want
			while((lookupFloat != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
				
                if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);	} outputTypes.clear();  delete lookupFloat; m->clearGroups();  return 0; }
				
				if(allLines == 1 || labels.count(lookupFloat->getLabel()) == 1){
					
					m->mothurOut(lookupFloat->getLabel()); m->mothurOutEndLine();
					
					normalize(lookupFloat);
					
					processedLabels.insert(lookupFloat->getLabel());
					userLabels.erase(lookupFloat->getLabel());
				}
				
				if ((m->anyLabelsToProcess(lookupFloat->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = lookupFloat->getLabel();
					
					delete lookupFloat;
					lookupFloat = input.getSharedRAbundFloatVectors(lastLabel);
					
					m->mothurOut(lookupFloat->getLabel()); m->mothurOutEndLine();
		
					normalize(lookupFloat);
					
					processedLabels.insert(lookupFloat->getLabel());
					userLabels.erase(lookupFloat->getLabel());
					
					//restore real lastlabel to save below
					lookupFloat->setLabels(saveLabel);
				}
				
				lastLabel = lookupFloat->getLabel();
				//prevent memory leak
				delete lookupFloat;
				
				if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);	} outputTypes.clear(); m->clearGroups();   return 0; }
				
				//get next line to process
				lookupFloat = input.getSharedRAbundFloatVectors();
			}
			
			if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);	} outputTypes.clear(); m->clearGroups();   return 0; }
			
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
				delete lookupFloat;
				lookupFloat = input.getSharedRAbundFloatVectors(lastLabel);
				
				m->mothurOut(lookupFloat->getLabel()); m->mothurOutEndLine();
				
				normalize(lookupFloat);
				
				delete lookupFloat;
			}
			
		}
		//reset groups parameter
		m->clearGroups();
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);	} outputTypes.clear(); return 0;}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		//m->mothurOut(outputFileName); m->mothurOutEndLine(); outputNames.push_back(outputFileName); outputTypes["shared"].push_back(outputFileName);
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		//set shared file as new current sharedfile
		string current = "";
		itTypes = outputTypes.find("shared");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setSharedFile(current); }
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "NormalizeSharedCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int NormalizeSharedCommand::normalize(SharedRAbundVectors*& thisLookUp){
	try {
		//save mothurOut's binLabels to restore for next label
		vector<string> saveBinLabels = m->currentSharedBinLabels;
        vector<string> lookupGroups = thisLookUp->getNamesGroups();
		
        map<string, string> variables; 
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(inputfile));
        variables["[distance]"] = thisLookUp->getLabel();
		string outputFileName = getOutputFileName("shared",variables);
        
		ofstream out;
		m->openOutputFile(outputFileName, out);
		outputNames.push_back(outputFileName); outputTypes["shared"].push_back(outputFileName);
				
		if (method == "totalgroup") { 
			
			//save numSeqs since they will change as the data is normalized
			vector<int> sizes;
			for (int i = 0; i < lookupGroups.size(); i++) {  sizes.push_back(thisLookUp->getNumSeqs(lookupGroups[i])); }
					
            for (int j = 0; j < thisLookUp->getNumBins(); j++) {
                
                for (int i = 0; i < lookupGroups.size(); i++) {
                    
                    if (m->control_pressed) { out.close(); return 0; }
                    
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
				
				if (m->control_pressed) { out.close(); return 0; }
				
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
					if (standardDev != 0) { // stop divide by zero
						float newNorm = ((float)thisLookUp->get(j, lookupGroups[i]) - mean) / standardDev;
						//round to nearest int
						finalNorm = (int) floor((newNorm + 0.5));
					}
					
					thisLookUp->set(j, finalNorm, lookupGroups[i]);
				}
			}
						
		}else{ m->mothurOut(method + " is not a valid scaling option."); m->mothurOutEndLine(); m->control_pressed = true; return 0; }
				
		thisLookUp->eliminateZeroOTUS();
		thisLookUp->printHeaders(out);
        thisLookUp->print(out);
		out.close();
		
		m->currentSharedBinLabels = saveBinLabels;
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "NormalizeSharedCommand", "normalize");
		exit(1);
	}
}
//**********************************************************************************************************************

int NormalizeSharedCommand::normalize(SharedRAbundFloatVectors*& thisLookUp){
	try {
		
		//save mothurOut's binLabels to restore for next label
		vector<string> saveBinLabels = m->currentSharedBinLabels;
        vector<string> lookupGroups = thisLookUp->getNamesGroups();
		
        map<string, string> variables; 
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(inputfile));
        variables["[distance]"] = thisLookUp->getLabel();
		string outputFileName = getOutputFileName("shared",variables);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		outputNames.push_back(outputFileName); outputTypes["shared"].push_back(outputFileName);
		
		if (method == "totalgroup") { 
			
			//save numSeqs since they will change as the data is normalized
			vector<float> sizes;
			for (int i = 0; i < lookupGroups.size(); i++) {  sizes.push_back(thisLookUp->getNumSeqs(lookupGroups[i])); }
			
			for (int j = 0; j < thisLookUp->getNumBins(); j++) {
				
				for (int i = 0; i < lookupGroups.size(); i++) {
					
					if (m->control_pressed) { out.close(); return 0; }
					
					float abund = thisLookUp->get(j, lookupGroups[i]);
					
					float relabund = abund / (float) sizes[i];
					float newNorm = relabund * norm;
					
					thisLookUp->set(j, newNorm, lookupGroups[i]);
				}
			}
			
		}else if (method == "zscore") {
			for (int j = 0; j < thisLookUp->getNumBins(); j++) {
				
				if (m->control_pressed) { out.close(); return 0; }
				
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
                    if (standardDev != 0) { // stop divide by zero
                        newNorm = ((float)thisLookUp->get(j, lookupGroups[i]) - mean) / standardDev;
                    }
                    thisLookUp->set(j, newNorm, lookupGroups[i]);
                }
			}			
			
		}else{ m->mothurOut(method + " is not a valid scaling option."); m->mothurOutEndLine(); m->control_pressed = true; return 0; }
		
		
        thisLookUp->eliminateZeroOTUS();
        thisLookUp->printHeaders(out);
        thisLookUp->print(out);
        out.close();
		
		m->currentSharedBinLabels = saveBinLabels;
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "NormalizeSharedCommand", "normalize");
		exit(1);
	}
}
//**********************************************************************************************************************
