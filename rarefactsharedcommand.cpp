/*
 *  rarefactsharedcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/6/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "rarefactsharedcommand.h"
#include "sharedsobs.h"
#include "sharednseqs.h"
#include "sharedutilities.h"

//**********************************************************************************************************************
vector<string> RareFactSharedCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pshared);
        CommandParameter pdesign("design", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pdesign);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pfreq("freq", "Number", "", "100", "", "", "",false,false); parameters.push_back(pfreq);
		CommandParameter piters("iters", "Number", "", "1000", "", "", "",false,false); parameters.push_back(piters);
		CommandParameter pcalc("calc", "Multiple", "sharednseqs-sharedobserved", "sharedobserved", "", "", "",true,false); parameters.push_back(pcalc);
		CommandParameter pjumble("jumble", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(pjumble);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
        CommandParameter psets("sets", "String", "", "", "", "", "",false,false); parameters.push_back(psets);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactSharedCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string RareFactSharedCommand::getHelpString(){	
	try {
		string helpString = "";
		ValidCalculators validCalculator;
		helpString += "The collect.shared command parameters are shared, label, freq, calc and groups.  shared is required if there is no current sharedfile. \n";
		helpString += "The rarefaction.shared command parameters are shared, design, label, iters, groups, sets, jumble and calc.  shared is required if there is no current sharedfile. \n";
        helpString += "The design parameter allows you to assign your groups to sets. If provided mothur will run rarefaction.shared on a per set basis. \n";
        helpString += "The sets parameter allows you to specify which of the sets in your designfile you would like to analyze. The set names are separated by dashes. THe default is all sets in the designfile.\n";
		helpString += "The rarefaction command should be in the following format: \n";
		helpString += "rarefaction.shared(label=yourLabel, iters=yourIters, calc=yourEstimators, jumble=yourJumble, groups=yourGroups).\n";
		helpString += "The freq parameter is used indicate when to output your data, by default it is set to 100. But you can set it to a percentage of the number of sequence. For example freq=0.10, means 10%. \n";
		helpString += "Example rarefaction.shared(label=unique-0.01-0.03,  iters=10000, groups=B-C, jumble=T, calc=sharedobserved).\n";
		helpString += "The default values for iters is 1000, freq is 100, and calc is sharedobserved which calculates the shared rarefaction curve for the observed richness.\n";
		helpString += "The default value for groups is all the groups in your groupfile, and jumble is true.\n";
		helpString += validCalculator.printCalc("sharedrarefaction");
		helpString += "The label parameter is used to analyze specific labels in your input.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 2 valid groups.\n";
		helpString += "Note: No spaces between parameter labels (i.e. freq), '=' and parameters (i.e.yourFreq).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactSharedCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
RareFactSharedCommand::RareFactSharedCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["sharedrarefaction"] = tempOutNames;
		outputTypes["sharedr_nseqs"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactSharedCommand", "RareFactSharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

RareFactSharedCommand::RareFactSharedCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
				
		//allow user to run help
		if(option == "help") {  help(); abort = true; calledHelp = true; }
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
			outputTypes["sharedrarefaction"] = tempOutNames;
			outputTypes["sharedr_nseqs"] = tempOutNames;
			
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
                
                it = parameters.find("design");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["design"] = inputDir + it->second;		}
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
            
            designfile = validParameter.validFile(parameters, "design", true);
			if (designfile == "not open") { abort = true; designfile = ""; }
			else if (designfile == "not found") {  	designfile = "";	}
			else { m->setDesignFile(designfile); }
			
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(sharedfile);		}
			
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
				
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "sharedobserved";  }
			else { 
				 if (calc == "default")  {  calc = "sharedobserved";  }
			}
			m->splitAtDash(calc, Estimators);
			if (m->inUsersGroups("citation", Estimators)) { 
				ValidCalculators validCalc; validCalc.printCitations(Estimators); 
				//remove citation from list of calcs
				for (int i = 0; i < Estimators.size(); i++) { if (Estimators[i] == "citation") {  Estimators.erase(Estimators.begin()+i); break; } }
			}
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				m->splitAtDash(groups, Groups);
			}
			m->setGroups(Groups);
            
            string sets = validParameter.validFile(parameters, "sets", false);			
			if (sets == "not found") { sets = ""; }
			else { 
				m->splitAtDash(sets, Sets);
			}
			
			string temp;
			temp = validParameter.validFile(parameters, "freq", false);			if (temp == "not found") { temp = "100"; }
			m->mothurConvert(temp, freq); 
			
			temp = validParameter.validFile(parameters, "iters", false);			if (temp == "not found") { temp = "1000"; }
			m->mothurConvert(temp, nIters); 
			
			temp = validParameter.validFile(parameters, "jumble", false);			if (temp == "not found") { temp = "T"; }
			if (m->isTrue(temp)) { jumble = true; }
			else { jumble = false; }
			m->jumble = jumble;
				
		}

	}
	catch(exception& e) {
		m->errorOut(e, "RareFactSharedCommand", "RareFactSharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int RareFactSharedCommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        GroupMap designMap;
        if (designfile == "") { //fake out designMap to run with process
            process(designMap, "");
        }else {
            designMap.readDesignMap(designfile);
            
            //fill Sets - checks for "all" and for any typo groups
			SharedUtil util;
			vector<string> nameSets = designMap.getNamesOfGroups();
			util.setGroups(Sets, nameSets);
            
            for (int i = 0; i < Sets.size(); i++) {
                process(designMap, Sets[i]);
            }
        }
                    
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactSharedCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int RareFactSharedCommand::process(GroupMap& designMap, string thisSet){
	try {
        Rarefact* rCurve;
        vector<Display*> rDisplays;
        
        InputData input(sharedfile, "sharedfile");
		lookup = input.getSharedRAbundVectors();
        if (lookup.size() < 2) { 
			m->mothurOut("I cannot run the command without at least 2 valid groups."); 
            for (int i = 0; i < lookup.size(); i++) { delete lookup[i]; }
            return 0;
        }
        
        string fileNameRoot = outputDir + m->getRootName(m->getSimpleName(sharedfile));
        
        vector<string> newGroups = m->getGroups();
        if (thisSet != "") {  //make groups only filled with groups from this set so that's all inputdata will read
            vector<string> thisSets; thisSets.push_back(thisSet);
            newGroups = designMap.getNamesSeqs(thisSets);
            fileNameRoot += thisSet + ".";
        }
        
        vector<SharedRAbundVector*> subset;
        if (thisSet == "") {  subset.clear(); subset = lookup;  }
        else {//fill subset with this sets groups
            subset.clear();
            for (int i = 0; i < lookup.size(); i++) {
                if (m->inUsersGroups(lookup[i]->getGroup(), newGroups)) {
                    subset.push_back(lookup[i]);
                }
            }
        }
		
		ValidCalculators validCalculator;
		for (int i=0; i<Estimators.size(); i++) {
			if (validCalculator.isValidCalculator("sharedrarefaction", Estimators[i]) == true) { 
				if (Estimators[i] == "sharedobserved") { 
					rDisplays.push_back(new RareDisplay(new SharedSobs(), new SharedThreeColumnFile(fileNameRoot+"shared.rarefaction", "")));
					outputNames.push_back(fileNameRoot+"shared.rarefaction"); outputTypes["sharedrarefaction"].push_back(fileNameRoot+"shared.rarefaction");
				}else if (Estimators[i] == "sharednseqs") { 
					rDisplays.push_back(new RareDisplay(new SharedNSeqs(), new SharedThreeColumnFile(fileNameRoot+"shared.r_nseqs", "")));
					outputNames.push_back(fileNameRoot+"shared.r_nseqs"); outputTypes["sharedr_nseqs"].push_back(fileNameRoot+"shared.r_nseqs");
				}
			}
		}
		
		//if the users entered no valid calculators don't execute command
		if (rDisplays.size() == 0) { for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  return 0; }
		
		if (m->control_pressed) { 
			for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	}
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	}
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
			return 0;
		}
        
        
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
        string lastLabel = subset[0]->getLabel();
		set<string> processedLabels;
		set<string> userLabels = labels;
        
		//as long as you are not at the end of the file or done wih the lines you want
		while((subset[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			if (m->control_pressed) { 
				for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	}
				for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	}
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
				return 0;
			}
			
			if(allLines == 1 || labels.count(subset[0]->getLabel()) == 1){
				m->mothurOut(subset[0]->getLabel() + '\t' + thisSet); m->mothurOutEndLine();
				rCurve = new Rarefact(subset, rDisplays);
				rCurve->getSharedCurve(freq, nIters);
				delete rCurve;
                
				processedLabels.insert(subset[0]->getLabel());
				userLabels.erase(subset[0]->getLabel());
			}
			
			if ((m->anyLabelsToProcess(subset[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
                string saveLabel = subset[0]->getLabel();
                
                for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
                lookup = input.getSharedRAbundVectors(lastLabel);
                
                if (thisSet == "") {  subset.clear(); subset = lookup;  }
                else {//fill subset with this sets groups
                    subset.clear();
                    for (int i = 0; i < lookup.size(); i++) {
                        if (m->inUsersGroups(lookup[i]->getGroup(), newGroups)) {
                            subset.push_back(lookup[i]);
                        }
                    }
                }

                m->mothurOut(subset[0]->getLabel() + '\t' + thisSet); m->mothurOutEndLine();
                rCurve = new Rarefact(subset, rDisplays);
                rCurve->getSharedCurve(freq, nIters);
                delete rCurve;
                
                processedLabels.insert(subset[0]->getLabel());
                userLabels.erase(subset[0]->getLabel());
                
                //restore real lastlabel to save below
                subset[0]->setLabel(saveLabel);
			}
            
			
			lastLabel = subset[0]->getLabel();
			
			//get next line to process
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
			lookup = input.getSharedRAbundVectors();
            
            if (lookup[0] != NULL) {
                if (thisSet == "") {  subset.clear(); subset = lookup;  }
                else {//fill subset with this sets groups
                    subset.clear();
                    for (int i = 0; i < lookup.size(); i++) {
                        if (m->inUsersGroups(lookup[i]->getGroup(), newGroups)) {
                            subset.push_back(lookup[i]);
                        }
                    }
                }
            }else {  subset.clear(); subset.push_back(NULL); }

		}
		
		if (m->control_pressed) { 
			for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	}
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	}
			return 0;
		}
		
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
		
		if (m->control_pressed) { 
			for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	}
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	}
			return 0;
		}
		
		//run last label if you need to
		if (needToRun == true)  {
			for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) {	delete lookup[i]; }  } 
			lookup = input.getSharedRAbundVectors(lastLabel);
            
            if (thisSet == "") {  subset.clear(); subset = lookup;  }
            else {//fill subset with this sets groups
                subset.clear();
                for (int i = 0; i < lookup.size(); i++) {
                    if (m->inUsersGroups(lookup[i]->getGroup(), newGroups)) {
                        subset.push_back(lookup[i]);
                    }
                }
            }
            
			m->mothurOut(subset[0]->getLabel() + '\t' + thisSet); m->mothurOutEndLine();
			rCurve = new Rarefact(subset, rDisplays);
			rCurve->getSharedCurve(freq, nIters);
			delete rCurve;
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
		}
		
		for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	}	

        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "RareFactSharedCommand", "process");
		exit(1);
	}
}
//**********************************************************************************************************************
