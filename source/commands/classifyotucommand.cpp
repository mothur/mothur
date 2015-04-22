/*
 *  classifyotucommand.cpp
 *  Mothur
 *
 *  Created by westcott on 6/1/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "classifyotucommand.h"
#include "phylotree.h"
#include "phylosummary.h"
#include "sharedutilities.h"

//**********************************************************************************************************************
vector<string> ClassifyOtuCommand::setParameters(){	
	try {
		CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(plist);
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "none", "none","constaxonomy",false,true,true); parameters.push_back(ptaxonomy);
		CommandParameter preftaxonomy("reftaxonomy", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(preftaxonomy);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
        CommandParameter ppersample("persample", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(ppersample);
        CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pbasis("basis", "Multiple", "otu-sequence", "otu", "", "", "","",false,false); parameters.push_back(pbasis);
		CommandParameter pcutoff("cutoff", "Number", "", "51", "", "", "","",false,true); parameters.push_back(pcutoff);
		CommandParameter pprobs("probs", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pprobs);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyOtuCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ClassifyOtuCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The classify.otu command parameters are list, taxonomy, reftaxonomy, name, group, count, persample, cutoff, label, basis and probs.  The taxonomy and list parameters are required unless you have a valid current file.\n";
		helpString += "The reftaxonomy parameter allows you give the name of the reference taxonomy file used when you classified your sequences. Providing it will keep the rankIDs in the summary file static.\n";
		helpString += "The name parameter allows you add a names file with your taxonomy file.\n";
		helpString += "The group parameter allows you provide a group file to use in creating the summary file breakdown.\n";
		helpString += "The count parameter allows you add a count file associated with your list file. When using the count parameter mothur assumes your list file contains only uniques.\n";
        helpString += "The basis parameter allows you indicate what you want the summary file to represent, options are otu and sequence. Default is otu.\n";
		helpString += "For example consider the following basis=sequence could give Clostridiales	3	105	16	43	46, where 105 is the total number of sequences whose otu classified to Clostridiales.\n";
		helpString += "16 is the number of sequences in the otus from groupA, 43 is the number of sequences in the otus from groupB, and 46 is the number of sequences in the otus from groupC.\n";
		helpString += "Now for basis=otu could give Clostridiales	3	7	6	1	2, where 7 is the number of otus that classified to Clostridiales.\n";
		helpString += "6 is the number of otus containing sequences from groupA, 1 is the number of otus containing sequences from groupB, and 2 is the number of otus containing sequences from groupC.\n";
		helpString += "The label parameter allows you to select what distance levels you would like a output files created for, and is separated by dashes.\n";
        helpString += "The persample parameter allows you to find a consensus taxonomy for each group. Default=f\n";
		helpString += "The default value for label is all labels in your inputfile.\n";
		helpString += "The cutoff parameter allows you to specify a consensus confidence threshold for your taxonomy.  The default is 51, meaning 51%. Cutoff cannot be below 51.\n";
		helpString += "The probs parameter shuts off the outputting of the consensus confidence results. The default is true, meaning you want the confidence to be shown.\n";
		helpString += "The classify.otu command should be in the following format: classify.otu(taxonomy=yourTaxonomyFile, list=yourListFile, name=yourNamesFile, label=yourLabels).\n";
		helpString += "Example classify.otu(taxonomy=abrecovery.silva.full.taxonomy, list=abrecovery.fn.list, label=0.10).\n";
		helpString += "Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListFile).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyOtuCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ClassifyOtuCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "constaxonomy") {  pattern = "[filename],[distance],cons.taxonomy"; } 
        else if (type == "taxsummary") {  pattern = "[filename],[distance],cons.tax.summary"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ClassifyOtuCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ClassifyOtuCommand::ClassifyOtuCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["constaxonomy"] = tempOutNames;
		outputTypes["taxsummary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyOtuCommand", "ClassifyOtuCommand");
		exit(1);
	}
}

//**********************************************************************************************************************
ClassifyOtuCommand::ClassifyOtuCommand(string option)  {
	try{
		abort = false; calledHelp = false;   
		allLines = 1;
		labels.clear();
				
		//allow user to run help
		if (option == "help") { 
			help(); abort = true; calledHelp = true;
		}else if(option == "citation") { citation(); abort = true; calledHelp = true;} 
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
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["constaxonomy"] = tempOutNames;
			outputTypes["taxsummary"] = tempOutNames;
		
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
				
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
				
				it = parameters.find("reftaxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["reftaxonomy"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
			}

			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			//check for required parameters
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not found") {				
				//if there is a current list file, use it
				listfile = m->getListFile(); 
				if (listfile != "") {  m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current listfile and the list parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (listfile == "not open") { abort = true; }	
			else { m->setListFile(listfile); }
			
			taxfile = validParameter.validFile(parameters, "taxonomy", true);
			if (taxfile == "not found") {  //if there is a current list file, use it
				taxfile = m->getTaxonomyFile(); 
				if (taxfile != "") {  m->mothurOut("Using " + taxfile + " as input file for the taxonomy parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current taxonomy file and the taxonomy parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (taxfile == "not open") { abort = true; }
			else { m->setTaxonomyFile(taxfile); }
			
			refTaxonomy = validParameter.validFile(parameters, "reftaxonomy", true);
			if (refTaxonomy == "not found") { refTaxonomy = ""; m->mothurOut("reftaxonomy is not required, but if given will keep the rankIDs in the summary file static."); m->mothurOutEndLine(); }
			else if (refTaxonomy == "not open") { abort = true; }
	
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			else { m->setNameFile(namefile); }
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }	
			else if (groupfile == "not found") { groupfile = ""; }
			else { m->setGroupFile(groupfile); }
            
            countfile = validParameter.validFile(parameters, "count", true);
			if (countfile == "not open") { countfile = ""; abort = true; }
			else if (countfile == "not found") { countfile = "";  }	
			else { m->setCountTableFile(countfile); }
            
            if ((namefile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: name or count."); m->mothurOutEndLine(); abort = true;
            }
			
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count."); m->mothurOutEndLine(); abort=true;
            }

			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; allLines = 1;  }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
            
			basis = validParameter.validFile(parameters, "basis", false);
			if (basis == "not found") { basis = "otu"; }	
			
			if ((basis != "otu") && (basis != "sequence")) { m->mothurOut("Invalid option for basis. basis options are otu and sequence, using otu."); m->mothurOutEndLine(); }
			
			string temp = validParameter.validFile(parameters, "cutoff", false);			if (temp == "not found") { temp = "51"; }
			m->mothurConvert(temp, cutoff); 
			
			temp = validParameter.validFile(parameters, "probs", false);					if (temp == "not found"){	temp = "true";			}
			probs = m->isTrue(temp);
            
            temp = validParameter.validFile(parameters, "persample", false);		if (temp == "not found"){	temp = "f";		}
			persample = m->isTrue(temp);
			
            if ((groupfile == "") && (countfile == "")) { if (persample) { m->mothurOut("persample is only valid with a group file, or count file with group information. Setting persample=f.\n"); persample = false; } 
            }
            if (countfile != "") {
                CountTable cts;
                if (!cts.testGroups(countfile)) { 
                    if (persample) { m->mothurOut("persample is only valid with a group file, or count file with group information. Setting persample=f.\n"); persample = false; }
                }
            }
			
			if ((cutoff < 51) || (cutoff > 100)) { m->mothurOut("cutoff must be above 50, and no greater than 100."); m->mothurOutEndLine(); abort = true;  }
			
            if (countfile == "") {
                if (namefile == ""){
                    vector<string> files; files.push_back(taxfile);
                    parser.getNameFile(files);
                }
            }
			
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyOtuCommand", "ClassifyOtuCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int ClassifyOtuCommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//if user gave a namesfile then use it
		if (namefile != "")     {	m->readNames(namefile, nameMap, true);	}
        if (groupfile != "")    {   groupMap = new GroupMap(groupfile);  groupMap->readMap();  groups = groupMap->getNamesOfGroups(); }
        else { groupMap = NULL;  }
        if (countfile != "") {  ct = new CountTable(); ct->readTable(countfile, true, false);  if (ct->hasGroupInfo()) { groups = ct->getNamesOfGroups(); } }
        else {  ct = NULL;    }
        
		//read taxonomy file and save in map for easy access in building bin trees
		m->readTax(taxfile, taxMap);
		
		if (m->control_pressed) { return 0; }
		
		input = new InputData(listfile, "list");
		list = input->getListVector();
		string lastLabel = list->getLabel();

		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		if (m->control_pressed) { outputTypes.clear(); if (ct != NULL) { delete ct; } if (groupMap != NULL) { delete groupMap; } delete input; delete list; for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  }  return 0; }
	
		while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (allLines == 1 || labels.count(list->getLabel()) == 1){
			
					m->mothurOut(list->getLabel() + "\t" + toString(list->size())); m->mothurOutEndLine();
					process(list);
					if (m->control_pressed) { outputTypes.clear(); for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } if (ct != NULL) { delete ct; } if (groupMap != NULL) { delete groupMap; } delete input; delete list; return 0; }
										
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
			}
			
			if ((m->anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = list->getLabel();
					
					delete list;
					list = input->getListVector(lastLabel);
					m->mothurOut(list->getLabel() + "\t" + toString(list->size())); m->mothurOutEndLine();
					process(list);
				
					
					if (m->control_pressed) { outputTypes.clear(); if (ct != NULL) { delete ct; }  if (groupMap != NULL) { delete groupMap; } for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } delete input; delete list; return 0; }
										
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
					
					//restore real lastlabel to save below
					list->setLabel(saveLabel);
			}
			
			lastLabel = list->getLabel();
	
			delete list;
			list = input->getListVector();
		}
		
		//output error messages about any remaining user labels
		bool needToRun = false;
		for (set<string>::iterator it = userLabels.begin(); it != userLabels.end(); it++) {  
			m->mothurOut("Your file does not include the label " + (*it)); 
			if (processedLabels.count(lastLabel) != 1) {
				m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
				needToRun = true;
			}else {
				m->mothurOut(". Please refer to " + lastLabel + "."); m->mothurOutEndLine();
			}
		}
		
		//run last label if you need to
		if (needToRun == true)  {
			if (list != NULL) {	delete list;	}
			list = input->getListVector(lastLabel);
			m->mothurOut(list->getLabel() + "\t" + toString(list->size())); m->mothurOutEndLine();
			
			process(list);
			delete list;
			
			if (m->control_pressed) { outputTypes.clear();  if (ct != NULL) { delete ct; } if (groupMap != NULL) { delete groupMap; } for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } delete input; delete list; return 0; }
		}
		
		delete input;  
        if (groupMap != NULL) { delete groupMap; }
        if (ct != NULL) { delete ct; }
				
		if (m->control_pressed) { outputTypes.clear(); for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyOtuCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> ClassifyOtuCommand::findConsensusTaxonomy(vector<string> names, int& size, string& conTax) {
	try{
		conTax = "";
		vector<string> allNames;
		map<string, string>::iterator it;
		map<string, string>::iterator it2;

		//create a tree containing sequences from this bin
		PhyloTree* phylo = new PhyloTree();
		
		size = 0;
		for (int i = 0; i < names.size(); i++) {
	
			//if namesfile include the names
			if (namefile != "") {
	
				//is this sequence in the name file - namemap maps seqName -> repSeqName
				it2 = nameMap.find(names[i]);
				
				if (it2 == nameMap.end()) { //this name is not in name file, skip it
					m->mothurOut(names[i] + " is not in your name file.  I will not include it in the consensus."); m->mothurOutEndLine();
				}else{
					
					//is this sequence in the taxonomy file - look for repSeqName since we are assuming the taxonomy file is unique
					it = taxMap.find(it2->second);
			
					if (it == taxMap.end()) { //this name is not in taxonomy file, skip it
					
						if (names[i] != it2->second) { m->mothurOut(names[i] + " is represented by " +  it2->second + " and is not in your taxonomy file.  I will not include it in the consensus."); m->mothurOutEndLine(); }
						else {  m->mothurOut(names[i] + " is not in your taxonomy file.  I will not include it in the consensus."); m->mothurOutEndLine(); }
					}else{
				
						//add seq to tree
						phylo->addSeqToTree(names[i], it->second);
						size++;
						allNames.push_back(names[i]);
					}
				}
				
			}else{
				//is this sequence in the taxonomy file - look for repSeqName since we are assuming the taxonomy file is unique
				it = taxMap.find(names[i]);
		
				if (it == taxMap.end()) { //this name is not in taxonomy file, skip it
					m->mothurOut(names[i] + " is not in your taxonomy file.  I will not include it in the consensus."); m->mothurOutEndLine();
				}else{
                    if (countfile != "") {
                        int numDups = ct->getNumSeqs(names[i]); 
                        for (int j = 0; j < numDups; j++) {  phylo->addSeqToTree(names[i], it->second);  }
                        size += numDups;
                    }else{
					//add seq to tree
                        phylo->addSeqToTree(names[i], it->second);
                        size++;  
                    }
                    allNames.push_back(names[i]);
				}
			}

			
			if (m->control_pressed) { delete phylo; return allNames; }
			
		}
		
		//build tree
		phylo->assignHeirarchyIDs(0);
		
		TaxNode currentNode = phylo->get(0);
		int myLevel = 0; 	
		//at each level
		while (currentNode.children.size() != 0) { //you still have more to explore
		
			TaxNode bestChild;
			int bestChildSize = 0;
			
			//go through children
			for (map<string, int>::iterator itChild = currentNode.children.begin(); itChild != currentNode.children.end(); itChild++) {
				
				TaxNode temp = phylo->get(itChild->second);
				
				//select child with largest accesions - most seqs assigned to it
				if (temp.accessions.size() > bestChildSize) {
					bestChild = phylo->get(itChild->second);
					bestChildSize = temp.accessions.size();
				}
				
			}
            
            //phylotree adds an extra unknown so we want to remove that
            if (bestChild.name == "unknown") { bestChildSize--; }
				
			//is this taxonomy above cutoff
			int consensusConfidence = ceil((bestChildSize / (float) size) * 100);
			
			if (consensusConfidence >= cutoff) { //if yes, add it
				if (probs) {
					conTax += bestChild.name + "(" + toString(consensusConfidence) + ");";
				}else{
					conTax += bestChild.name + ";";
				}
				myLevel++;
			}else{ //if no, quit
				break;
			}
			
			//move down a level
			currentNode = bestChild;
		}
		
		if (myLevel != phylo->getMaxLevel()) {
			while (myLevel != phylo->getMaxLevel()) {
                if (probs) {
                    conTax += "unclassified(100);";
                }else{
                    conTax += "unclassified;";
                }
				myLevel++;
			}
		}		
		if (conTax == "") {  conTax = "no_consensus;";  }
		
		delete phylo;	
		
		return allNames;
			
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyOtuCommand", "findConsensusTaxonomy");
		exit(1);
	}
}

//**********************************************************************************************************************
int ClassifyOtuCommand::process(ListVector* processList) {
	try{
		string conTax;
		int size;
		
		//create output file
		if (outputDir == "") { outputDir += m->hasPath(listfile); }
				
		ofstream out;
        map<string, string> variables; 
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(listfile));
        variables["[distance]"] = processList->getLabel();
		string outputFile = getOutputFileName("constaxonomy", variables);
		m->openOutputFile(outputFile, out);
		outputNames.push_back(outputFile); outputTypes["constaxonomy"].push_back(outputFile);
		
		ofstream outSum;
		string outputSumFile = getOutputFileName("taxsummary", variables);
		m->openOutputFile(outputSumFile, outSum);
		outputNames.push_back(outputSumFile); outputTypes["taxsummary"].push_back(outputSumFile);
		
		out << "OTU\tSize\tTaxonomy" << endl;
		
		PhyloSummary* taxaSum;
        if (countfile != "") {
            if (refTaxonomy != "") { taxaSum = new PhyloSummary(refTaxonomy, ct,false);  }
            else {  taxaSum = new PhyloSummary(ct,false); }
		}else {
            if (refTaxonomy != "") { taxaSum = new PhyloSummary(refTaxonomy, groupMap,false);  }
            else {  taxaSum = new PhyloSummary(groupMap,false); }
        }
        
        vector<ofstream*> outSums;
        vector<ofstream*> outs;
        vector<PhyloSummary*> taxaSums;
        map<string, int> groupIndex;
        if (persample) {
            for (int i = 0; i < groups.size(); i++) {
                groupIndex[groups[i]] = i;
                ofstream* temp = new ofstream();
                variables["[distance]"] = processList->getLabel() + "." + groups[i];
                string outputFile = getOutputFileName("constaxonomy", variables);
                m->openOutputFile(outputFile, *temp);
                (*temp) << "OTU\tSize\tTaxonomy" << endl;
                outs.push_back(temp);
                outputNames.push_back(outputFile); outputTypes["constaxonomy"].push_back(outputFile);
                
                ofstream* tempSum = new ofstream();
                string outputSumFile = getOutputFileName("taxsummary", variables);
                m->openOutputFile(outputSumFile, *tempSum);
                outSums.push_back(tempSum);
                outputNames.push_back(outputSumFile); outputTypes["taxsummary"].push_back(outputSumFile);
                
                PhyloSummary* taxaSumt;
                if (countfile != "") {
                    if (refTaxonomy != "") { taxaSumt = new PhyloSummary(refTaxonomy, ct, false);  }
                    else {  taxaSumt = new PhyloSummary(ct, false); }
                }else {
                    if (refTaxonomy != "") { taxaSumt = new PhyloSummary(refTaxonomy, groupMap,false);  }
                    else {  taxaSumt = new PhyloSummary(groupMap,false); }
                }
                taxaSums.push_back(taxaSumt);
            }
        }
        
		//for each bin in the list vector
        string snumBins = toString(processList->getNumBins());
        vector<string> binLabels = processList->getLabels();
		for (int i = 0; i < processList->getNumBins(); i++) {
			
			if (m->control_pressed) { break; }
			
			vector<string> names;
            string binnames = processList->get(i);
            vector<string> thisNames;
            m->splitAtComma(binnames, thisNames);
            
			names = findConsensusTaxonomy(thisNames, size, conTax);
		
			if (m->control_pressed) { break; }

			out << binLabels[i] << '\t' << size << '\t' << conTax << endl;
			
			string noConfidenceConTax = conTax;
			m->removeConfidences(noConfidenceConTax);
			
			//add this bins taxonomy to summary
			if (basis == "sequence") {
				for(int j = 0; j < names.size(); j++) {  
                    //int numReps = 1;
                    //if (countfile != "") {  numReps = ct->getNumSeqs(names[j]); }
                    //for(int k = 0; k < numReps; k++) {  taxaSum->addSeqToTree(names[j], noConfidenceConTax);  }
                    taxaSum->addSeqToTree(names[j], noConfidenceConTax);
                }
			}else { //otu
                map<string, bool> containsGroup; 
                if (countfile != "") {
                    if (ct->hasGroupInfo()) {
                        vector<string> mGroups = ct->getNamesOfGroups();
                        for (int k = 0; k < names.size(); k++) {
                            vector<int> counts = ct->getGroupCounts(names[k]);
                            for (int h = 0; h < counts.size(); h++) {  
                                if (counts[h] != 0) {  containsGroup[mGroups[h]] = true; }
                            }
                        }
                    }
                }else {
                    if (groupfile != "") {
                        vector<string> mGroups = groupMap->getNamesOfGroups();
                        for (int j = 0; j < mGroups.size(); j++) { containsGroup[mGroups[j]] = false; }
                        
                        for (int k = 0; k < names.size(); k++) {
                            //find out the sequences group
                            string group = groupMap->getGroup(names[k]);
                            
                            if (group == "not found") {  m->mothurOut("[WARNING]: " + names[k] + " is not in your groupfile, and will be included in the overall total, but not any group total."); m->mothurOutEndLine();  }
                            else {
                                containsGroup[group] = true;
                            }
                        }
                    }
                }
				taxaSum->addSeqToTree(noConfidenceConTax, containsGroup);
			}
            
            
            if (persample) {
                //divide names by group
                map<string, vector<string> > parsedNames;
                map<string, vector<string> >::iterator itParsed;
                
                //parse names by group
                for (int j = 0; j < names.size(); j++) {
                    if (groupfile != "") { 
                        string group = groupMap->getGroup(names[j]); 
                        itParsed = parsedNames.find(group);
                        
                        if (itParsed != parsedNames.end()) { itParsed->second.push_back(names[j]); }
                        else { vector<string> tempNames; tempNames.push_back(names[j]); parsedNames[group] = tempNames; }
                    }else { //count file was used
                        vector<string> thisSeqsGroups = ct->getGroups(names[j]);
                        for (int k = 0; k < thisSeqsGroups.size(); k++) {
                            string group = thisSeqsGroups[k]; 
                            itParsed = parsedNames.find(group);
                            
                            if (itParsed != parsedNames.end()) { itParsed->second.push_back(names[j]); }
                            else { vector<string> tempNames; tempNames.push_back(names[j]); parsedNames[group] = tempNames; }
                        }
                    }
                }
                
                for (itParsed = parsedNames.begin(); itParsed != parsedNames.end(); itParsed++) {
                    vector<string> theseNames = findConsensusTaxonomy(itParsed->second, size, conTax);
                    
                    if (m->control_pressed) { break; }
                    
                    
                    (*outs[groupIndex[itParsed->first]]) << binLabels[i] << '\t' << size << '\t' << conTax << endl;
                    
                    string noConfidenceConTax = conTax;
                    m->removeConfidences(noConfidenceConTax);
                    
                    //add this bins taxonomy to summary
                    if (basis == "sequence") {
                        for(int j = 0; j < theseNames.size(); j++) {  
                            int numReps = 1;
                            if (countfile != "") {  numReps = ct->getGroupCount(theseNames[j], itParsed->first); } //get num seqs for this seq from this group
                            for(int k = 0; k < numReps; k++) {  (taxaSums[groupIndex[itParsed->first]])->addSeqToTree(theseNames[j], noConfidenceConTax);  }
                        }
                    }else { //otu
                        map<string, bool> containsGroup; 
                        containsGroup[itParsed->first] = true;
                        (taxaSums[groupIndex[itParsed->first]])->addSeqToTree(noConfidenceConTax, containsGroup);
                    }
                }
            }
		}

		out.close();
		
		//print summary file
		taxaSum->print(outSum);
		outSum.close();
        
        if (persample) {
            for (int i = 0; i < groups.size(); i++) {
                (*outs[i]).close();
                taxaSums[i]->print(*outSums[i]);
                (*outSums[i]).close();
                delete outs[i];
                delete outSums[i];
                delete taxaSums[i];
            }
        }
		
		delete taxaSum;
		
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyOtuCommand", "process");
		exit(1);
	}
}
/**************************************************************************************************/
string ClassifyOtuCommand::addUnclassifieds(string tax, int maxlevel) {
	try{
		string newTax, taxon;
		int level = 0;
		
		//keep what you have counting the levels
		while (tax.find_first_of(';') != -1) {
			//get taxon
			taxon = tax.substr(0,tax.find_first_of(';'))+';';
			tax = tax.substr(tax.find_first_of(';')+1, tax.length());
			newTax += taxon;
			level++;
		}
		
		//add "unclassified" until you reach maxLevel
		while (level < maxlevel) {
			newTax += "unclassified;";
			level++;
		}
		
		return newTax;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyOtuCommand", "addUnclassifieds");
		exit(1);
	}
}
//**********************************************************************************************************************


