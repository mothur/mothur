//
//  filtersharedcommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 1/4/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "filtersharedcommand.h"

//**********************************************************************************************************************
vector<string> FilterSharedCommand::setParameters(){	
	try {		
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none","shared",false,true,true); parameters.push_back(pshared);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter pminpercent("minpercent", "Number", "", "-1", "", "", "","",false,false,true); parameters.push_back(pminpercent);
        CommandParameter prarepercent("rarepercent", "Number", "", "-1", "", "", "","",false,false,true); parameters.push_back(prarepercent);
        CommandParameter pminabund("minabund", "Number", "", "-1", "", "", "","",false,false,true); parameters.push_back(pminabund);
        CommandParameter pmintotal("mintotal", "Number", "", "-1", "", "", "","",false,false,true); parameters.push_back(pmintotal);
        CommandParameter pminnumsamples("minnumsamples", "Number", "", "-1", "", "", "","",false,false,true); parameters.push_back(pminnumsamples);
        CommandParameter pminpercentsamples("minpercentsamples", "Number", "", "-1", "", "", "","",false,false,true); parameters.push_back(pminpercentsamples);
        CommandParameter pkeepties("keepties", "Boolean", "", "T", "", "", "","",false,false,true); parameters.push_back(pkeepties);
		CommandParameter pmakerare("makerare", "Boolean", "", "T", "", "", "","",false,false,true); parameters.push_back(pmakerare);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSharedCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string FilterSharedCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The filter.shared command is used to remove OTUs based on various critieria.\n";
		helpString += "The filter.shared command parameters are shared, minpercent, minabund, mintotal, minnumsamples, minpercentsamples, rarepercent, makerare, keepties, groups and label.  You must provide a shared file.\n";
		helpString += "The groups parameter allows you to specify which of the groups you would like included. The group names are separated by dashes.\n";
		helpString += "The label parameter allows you to select what distance levels you would like, and are also separated by dashes.\n";
        
		helpString += "The minabund parameter allows you indicate the minimum abundance required for each sample in a given OTU.  If any samples abundance falls below the minimum, the OTU is removed. Default=0\n";
        helpString += "The minpercent parameter allows you indicate the minimum relative abundance of an OTU. For example, if the OTUs total abundance across all samples is 8, and the total abundance across all OTUs is 1000, and minpercent=1. The OTU's relative abundance is 0.008, the minimum is 0.01, so the OTU will be removed. Default=0.\n";
        helpString += "The rarepercent parameter allows you indicate the percentage of otus to remove.  The OTUs chosen to be removed are the rarest.  For example if you have 1000 OTUs, rarepercent=20 would remove the 200 OTUs with the lowest abundance. Default=0.\n";
        helpString += "The keepties parameter is used with the rarepercent parameter.  It allows you indicate you want to keep the OTUs with the same abundance as the first 'not rare' OTU. For example if you have 10 OTUs, rarepercent=20 abundances of 20, 18, 15, 15, 10, 5, 3, 3, 3, 1. keepties=t, would remove the 10th OTU, but keep the 9th because its abundance ties the 8th OTU. keepties=f would remove OTUs 9 and 10.  Default=T\n";   
        helpString += "The minnumsamples parameter allows you indicate the minimum number of samples present in an OTU. If the number of samples present falls below the minimum, the OTU is removed. Default=0.\n";
        helpString += "The minpercentsamples parameter allows you indicate the minimum percent of sample present in an OTU. For example, if the total number of samples is 10, the number present is 3, and the minpercentsamples=50. The OTU's precent of samples is 0.333, the minimum is 0.50, so the OTU will be removed. Default=0.\n";
        helpString += "The mintotal parameter allows you indicate the minimum abundance required for a given OTU.  If abundance across all samples falls below the minimum, the OTU is removed. Default=0.\n";
        
		helpString += "The makerare parameter allows you indicate you want the abundances of any removed OTUs to be saved and a new \"rare\" OTU created with its abundances equal to the sum of the OTUs removed.  This will preserve the number of reads in your dataset. Default=T\n";
		helpString += "The filter.shared command should be in the following format: filter.shared(shared=yourSharedFile, minabund=yourMinAbund, groups=yourGroups, label=yourLabels).\n";
		helpString += "Example filter.shared(shared=final.an.shared, minabund=10).\n";
		helpString += "The default value for groups is all the groups in your sharedfile, and all labels in your inputfile will be used.\n";
		helpString += "The filter.shared command outputs a .filter.shared file.\n";
		helpString += "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSharedCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string FilterSharedCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "shared")      {   pattern = "[filename],[distance],filter,[extension]";    }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "FilterSharedCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
FilterSharedCommand::FilterSharedCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["shared"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSharedCommand", "GetRelAbundCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
FilterSharedCommand::FilterSharedCommand(string option) {
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
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			map<string,string>::iterator it;
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
            }
					
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }	
			else if (sharedfile == "not found") { 
				//if there is a current shared file, use it
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current sharedfile and the shared parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setSharedFile(sharedfile); }
			
            //if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(sharedfile);	}
            
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
			
            bool setSomething = false;
			string temp = validParameter.validFile(parameters, "minabund", false);		
            if (temp == "not found"){	temp = "-1";		}
            else { setSomething = true; }
			m->mothurConvert(temp, minAbund);  
            
            temp = validParameter.validFile(parameters, "mintotal", false);		
            if (temp == "not found"){	temp = "-1";		}
            else { setSomething = true; }
			m->mothurConvert(temp, minTotal);
            
            temp = validParameter.validFile(parameters, "minnumsamples", false);		
            if (temp == "not found"){	temp = "-1";		}
            else { setSomething = true; }
			m->mothurConvert(temp, minSamples);
            
            temp = validParameter.validFile(parameters, "minpercent", false);		
            if (temp == "not found"){	temp = "-1";		}
            else { setSomething = true; }
            
			m->mothurConvert(temp, minPercent);
            if (minPercent < 1) {} //already in percent form
            else {  minPercent = minPercent / 100.0; } //user gave us a whole number version so convert to %
            
            temp = validParameter.validFile(parameters, "rarepercent", false);
            if (temp == "not found"){	temp = "-1";		}
            else { setSomething = true; }
            
			m->mothurConvert(temp, rarePercent);
            if (rarePercent < 1) {} //already in percent form
            else {  rarePercent = rarePercent / 100.0; } //user gave us a whole number version so convert to %
            
            temp = validParameter.validFile(parameters, "minpercentsamples", false);		
            if (temp == "not found"){	temp = "-1";		}
            else { setSomething = true; }
            
			m->mothurConvert(temp, minPercentSamples);
            if (minPercentSamples < 1) {} //already in percent form
            else {  minPercentSamples = minPercentSamples / 100.0; } //user gave us a whole number version so convert to %
			
			temp = validParameter.validFile(parameters, "makerare", false);		if (temp == "not found"){	temp = "T";		}
			makeRare = m->isTrue(temp);
            
            temp = validParameter.validFile(parameters, "keepties", false);		if (temp == "not found"){	temp = "T";		}
			keepties = m->isTrue(temp);
            
            if (!setSomething) { m->mothurOut("\nYou did not set any parameters. I will filter using minabund=1.\n\n"); minAbund = 1; }
        }
        
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSharedCommand", "FilterSharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int FilterSharedCommand::execute(){
	try {
        
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
        InputData input(sharedfile, "sharedfile");
		vector<SharedRAbundVector*> lookup = input.getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			if (m->control_pressed) {   for (int i = 0; i < lookup.size(); i++) {  delete lookup[i]; lookup[i] = NULL; }  
                for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0;  }
			
			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			
				
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				
				processShared(lookup);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}
			
			if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup[0]->getLabel();
				
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  
				
				lookup = input.getSharedRAbundVectors(lastLabel);
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				
				processShared(lookup);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				
				//restore real lastlabel to save below
				lookup[0]->setLabel(saveLabel);
			}
			
			lastLabel = lookup[0]->getLabel();
			//prevent memory leak
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i]; lookup[i] = NULL; }
			
			//get next line to process
			lookup = input.getSharedRAbundVectors();				
		}
		
		
		if (m->control_pressed) {   return 0;  }
		
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
			
			processShared(lookup);
			
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
		}
        
		//set shared file as new current sharedfile
		string current = "";
        itTypes = outputTypes.find("shared");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setSharedFile(current); }
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSharedCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int FilterSharedCommand::processShared(vector<SharedRAbundVector*>& thislookup) {
	try {
		
		//save mothurOut's binLabels to restore for next label
		vector<string> saveBinLabels = m->currentSharedBinLabels;
		
        map<string, string> variables; 
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(sharedfile));
        variables["[extension]"] = m->getExtension(sharedfile);
        variables["[distance]"] = thislookup[0]->getLabel();
		string outputFileName = getOutputFileName("shared", variables);        
        
        if (m->control_pressed) {  return 0; }
        
        map<string, int> labelsForRare;
        vector<string> filteredLabels;
        vector<int> rareCounts; rareCounts.resize(m->getGroups().size(), 0);
        
        //create new "filtered" lookup
        vector<SharedRAbundVector*> filteredLookup;
        for (int i = 0; i < thislookup.size(); i++) {
            SharedRAbundVector* temp = new SharedRAbundVector();
			temp->setLabel(thislookup[i]->getLabel());
			temp->setGroup(thislookup[i]->getGroup());
			filteredLookup.push_back(temp);
        }
        
        //you want to remove a percentage of OTUs
        set<string> removeLabels;
        if (rarePercent != -0.01) {
            vector<spearmanRank> otus;
            //rank otus by abundance
            for (int i = 0; i < thislookup[0]->getNumBins(); i++) {
                float otuTotal = 0.0;
                for (int j = 0; j < thislookup.size(); j++) {
                    otuTotal += thislookup[j]->getAbundance(i);
                }
                spearmanRank temp(saveBinLabels[i], otuTotal);
                otus.push_back(temp);
            }
            
            //sort by abundance
            sort(otus.begin(), otus.end(), compareSpearman);
            
            //find index of cutoff
            int indexFirstNotRare = ceil(rarePercent * (float)thislookup[0]->getNumBins());
            
            //handle ties
            if (keepties) { //adjust indexFirstNotRare if needed
                if (indexFirstNotRare != 0) { //not out of bounds
                    if (otus[indexFirstNotRare].score == otus[indexFirstNotRare-1].score) { //you have a tie
                        bool tie = true;
                        for (int i = indexFirstNotRare-1; i >=0; i--) {
                            if (otus[indexFirstNotRare].score != otus[i].score) { //found value below tie
                                indexFirstNotRare = i+1; tie = false; break;
                            }
                        }
                        if (tie) { if (m->debug) { m->mothurOut("For distance " + thislookup[0]->getLabel() + " all rare OTUs abundance tie with first 'non rare' OTU, not removing any for rarepercent parameter.\n"); }indexFirstNotRare = 0; }
                    }
                }
            }
            
            //saved labels for OTUs above rarepercent
            for (int i = 0; i < indexFirstNotRare; i++) { removeLabels.insert(otus[i].name); }
        }
        
        bool filteredSomething = false;
        int numRemoved = 0;
        for (int i = 0; i < thislookup[0]->getNumBins(); i++) {
            
            if (m->control_pressed) { for (int j = 0; j < filteredLookup.size(); j++) { delete filteredLookup[j]; } return 0; }
            
            bool okay = true; //innocent until proven guilty
            if (minAbund != -1) {
                for (int j = 0; j < thislookup.size(); j++) { 
                    if (thislookup[j]->getAbundance(i) < minAbund) { okay = false; break; }
                }
            }
            
            if (okay && (minTotal != -1)) {
                int otuTotal = 0;
                for (int j = 0; j < thislookup.size(); j++) { 
                    otuTotal += thislookup[j]->getAbundance(i);
                }
                if (otuTotal < minTotal) { okay = false; }
            }
            
            if (okay && (minPercent != -0.01)) {
                double otuTotal = 0;
                double total = 0;
                for (int j = 0; j < thislookup.size(); j++) { 
                    otuTotal += thislookup[j]->getAbundance(i);
                    total += thislookup[j]->getNumSeqs();
                }
                double percent = otuTotal / total; 
                if (percent < minPercent) { okay = false; }
            }
            
            
            if (okay && (minSamples != -1)) {
                int samples = 0;
                for (int j = 0; j < thislookup.size(); j++) { 
                    if (thislookup[j]->getAbundance(i) != 0) { samples++; }
                }
                if (samples < minSamples) { okay = false; }
            }
            
            if (okay && (minPercentSamples != -0.01)) {
                double samples = 0;
                double total = thislookup.size();
                for (int j = 0; j < thislookup.size(); j++) { 
                    if (thislookup[j]->getAbundance(i) != 0) { samples++; }
                }
                double percent = samples / total; 
                if (percent < minPercentSamples) { okay = false; }
            }
            
            if (okay && (rarePercent != -0.01)) {
                if (removeLabels.count(saveBinLabels[i]) != 0) { //are we on the 'bad' list
                    okay = false;
                }
            }
            
            //did this OTU pass the filter criteria
            if (okay) {
                filteredLabels.push_back(saveBinLabels[i]);
                labelsForRare[m->getSimpleLabel(saveBinLabels[i])] = i;
                for (int j = 0; j < filteredLookup.size(); j++) { //add this OTU to the filtered lookup
                    filteredLookup[j]->push_back(thislookup[j]->getAbundance(i), thislookup[j]->getGroup());
                }
            }else { //if not, do we want to save the counts
                filteredSomething = true;
                if (makeRare) {
                    for (int j = 0; j < rareCounts.size(); j++) {  rareCounts[j] += thislookup[j]->getAbundance(i); }
                }
                numRemoved++;
            }
            
        }
        
        //if we are saving the counts add a "rare" OTU if anything was filtered
        if (makeRare) {
            if (filteredSomething) {
                for (int j = 0; j < rareCounts.size(); j++) { //add "rare" OTU to the filtered lookup
                    filteredLookup[j]->push_back(rareCounts[j], thislookup[j]->getGroup());
                }
                //create label for rare OTUs
                map<string, int>::iterator it;
                int otuNum = 0; bool notDone = true;
                
                //find label prefix
                string prefix = "Otu";
                if (filteredLabels[filteredLabels.size()-1][0] == 'P') { prefix = "PhyloType"; }
                
                string tempLabel = filteredLabels[filteredLabels.size()-1];
                string simpleLastLabel = m->getSimpleLabel(tempLabel);
                m->mothurConvert(simpleLastLabel, otuNum); otuNum++;
                while (notDone) {
                    if (m->control_pressed) { notDone = false; break; }
                    
                    string potentialLabel = toString(otuNum);
                    it = labelsForRare.find(potentialLabel);
                    if (it == labelsForRare.end()) {
                        potentialLabel = prefix + toString(otuNum);
                        it = labelsForRare.find(potentialLabel);
                        if (it == labelsForRare.end()) {
                            notDone = false; break;
                        }
                    }
                    otuNum++;
                }
                filteredLabels.push_back("rareOTUs" + toString(otuNum));
            }
        }
        
        ofstream out;
		m->openOutputFile(outputFileName, out);
		outputTypes["shared"].push_back(outputFileName);  outputNames.push_back(outputFileName);
		
        m->currentSharedBinLabels = filteredLabels;
        
		filteredLookup[0]->printHeaders(out);
		
		for (int i = 0; i < filteredLookup.size(); i++) {
			out << filteredLookup[i]->getLabel() << '\t' << filteredLookup[i]->getGroup() << '\t';
			filteredLookup[i]->print(out);
		}
		out.close();
        
        
        //save mothurOut's binLabels to restore for next label
		m->currentSharedBinLabels = saveBinLabels;
        
        for (int j = 0; j < filteredLookup.size(); j++) { delete filteredLookup[j]; }
		
        m->mothurOut("\nRemoved " + toString(numRemoved) + " OTUs.\n");
        
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSharedCommand", "processShared");
		exit(1);
	}
}			
//**********************************************************************************************************************


