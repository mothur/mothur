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


//**********************************************************************************************************************
vector<string> ClassifyOtuCommand::setParameters(){	
	try {
		CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(plist);
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "none", "none","constaxonomy",false,true,true); parameters.push_back(ptaxonomy);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","",false,false,true); parameters.push_back(pcount);
        CommandParameter poutput("output", "Multiple", "plain-detail", "detail", "", "", "","",false,false, true); parameters.push_back(poutput);
        CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
        CommandParameter prelabund("relabund", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(prelabund);
        CommandParameter pprintlevel("printlevel", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pprintlevel);
        CommandParameter ppersample("persample", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(ppersample);
        CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pbasis("basis", "Multiple", "otu-sequence", "otu", "", "", "","",false,false); parameters.push_back(pbasis);
		CommandParameter pcutoff("cutoff", "Number", "", "51", "", "", "","",false,true); parameters.push_back(pcutoff);
        CommandParameter pthreshold("threshold", "Number", "", "0", "", "", "","",false,true); parameters.push_back(pthreshold);
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
		helpString += "The classify.otu command parameters are list, taxonomy, name, group, count, persample, cutoff, label, basis, relabund and probs.  The taxonomy and list parameters are required unless you have a valid current file.\n";
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
        helpString += "The relabund parameter allows you to indicate you want the summary file values to be relative abundances rather than raw abundances. Default=F. \n";
		helpString += "The default value for label is all labels in your inputfile.\n";
        helpString += "The output parameter allows you to specify format of your summary file. Options are simple and detail. The default is detail.\n";
        helpString += "The printlevel parameter allows you to specify taxlevel of your summary file to print to. Options are 1 to the maz level in the file.  The default is -1, meaning max level.  If you select a level greater than the level your sequences classify to, mothur will print to the level your max level. \n";
		helpString += "The cutoff parameter allows you to specify a consensus confidence threshold for your otu taxonomy output.  The default is 51, meaning 51%. Cutoff cannot be below 51.\n";
		helpString += "The probs parameter shuts off the outputting of the consensus confidence results. The default is true, meaning you want the confidence to be shown.\n";
        helpString += "The threshold parameter allows you to specify a cutoff for the taxonomy file that is being inputted. Once the classification falls below the threshold the mothur will refer to it as unclassified when calculating the concensus.  This feature is similar to adjusting the cutoff in classify.seqs. Default=0.\n";
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
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
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
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
				
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
			}

			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";		}
			
			//check for required parameters
			listfile = validParameter.validFile(parameters, "list");
			if (listfile == "not found") {				
				//if there is a current list file, use it
				listfile = current->getListFile(); 
				if (listfile != "") {  m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current listfile and the list parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (listfile == "not open") { abort = true; }	
			else { current->setListFile(listfile); }
			
			taxfile = validParameter.validFile(parameters, "taxonomy");
			if (taxfile == "not found") {  //if there is a current list file, use it
				taxfile = current->getTaxonomyFile(); 
				if (taxfile != "") {  m->mothurOut("Using " + taxfile + " as input file for the taxonomy parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current taxonomy file and the taxonomy parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (taxfile == "not open") { abort = true; }
			else { current->setTaxonomyFile(taxfile); }
	
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { namefile = ""; abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			else { current->setNameFile(namefile); }
			
			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not open") { abort = true; }	
			else if (groupfile == "not found") { groupfile = ""; }
			else { current->setGroupFile(groupfile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { countfile = ""; abort = true; }
			else if (countfile == "not found") { countfile = "";  }	
			else { current->setCountFile(countfile); }
            
            if ((namefile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: name or count."); m->mothurOutEndLine(); abort = true;
            }
			
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count."); m->mothurOutEndLine(); abort=true;
            }

			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; allLines = 1;  }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
            
			basis = validParameter.valid(parameters, "basis");
			if (basis == "not found") { basis = "otu"; }	
			
			if ((basis != "otu") && (basis != "sequence")) { m->mothurOut("Invalid option for basis. basis options are otu and sequence, using otu."); m->mothurOutEndLine(); }
			
			string temp = validParameter.valid(parameters, "cutoff");			if (temp == "not found") { temp = "51"; }
			util.mothurConvert(temp, cutoff);
            
            temp = validParameter.valid(parameters, "threshold");			if (temp == "not found") { temp = "0"; }
            util.mothurConvert(temp, threshold);
			
			temp = validParameter.valid(parameters, "probs");					if (temp == "not found"){	temp = "true";			}
			probs = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "persample");		if (temp == "not found"){	temp = "f";		}
			persample = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "relabund");		if (temp == "not found"){	temp = "false";			}
            relabund = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "printlevel");		if (temp == "not found"){	temp = "-1";		}
            util.mothurConvert(temp, printlevel);
            
            output = validParameter.valid(parameters, "output");		if(output == "not found"){	output = "detail"; }
            if ((output != "simple") && (output != "detail")) { m->mothurOut(output + " is not a valid output form. Options are simple and detail. I will use detail."); m->mothurOutEndLine(); output = "detail"; }
			
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
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
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
	
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		//if user gave a namesfile then use it
		if (namefile != "")     {	util.readNames(namefile, nameMap, true);	}
        if (groupfile != "")    {   groupMap = new GroupMap(groupfile);  groupMap->readMap();  groups = groupMap->getNamesOfGroups(); }
        else { groupMap = NULL;  }
        if (countfile != "") {  ct = new CountTable(); ct->readTable(countfile, true, false);  if (ct->hasGroupInfo()) { groups = ct->getNamesOfGroups(); } }
        else {  ct = NULL;    }
        
		//read taxonomy file and save in map for easy access in building bin trees
        bool removeConfidences = false;
        if (threshold == 0) { removeConfidences = true; }
		util.readTax(taxfile, taxMap, removeConfidences);
        
        if (threshold != 0) {  processTaxMap();  }
		
		if (m->getControl_pressed()) { return 0; }
		
		input = new InputData(listfile, "list", nullVector);
		list = input->getListVector();
		string lastLabel = list->getLabel();

		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		if (m->getControl_pressed()) { outputTypes.clear(); if (ct != NULL) { delete ct; } if (groupMap != NULL) { delete groupMap; } delete input; delete list; for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);  }  return 0; }
	
		while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (allLines == 1 || labels.count(list->getLabel()) == 1){
                    string output = toString(list->size());  if (basis == "sequence") { output = toString(list->getNumSeqs()); }
                
					m->mothurOut(list->getLabel() + "\t" + toString(output)); m->mothurOutEndLine();
					process(list);
					if (m->getControl_pressed()) { outputTypes.clear(); for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);  } if (ct != NULL) { delete ct; } if (groupMap != NULL) { delete groupMap; } delete input; delete list; return 0; }
										
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
			}
			
			if ((util.anyLabelsToProcess(list->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = list->getLabel();
					
					delete list;
					list = input->getListVector(lastLabel);
                
                    string output = toString(list->size());  if (basis == "sequence") { output = toString(list->getNumSeqs()); }
					m->mothurOut(list->getLabel() + "\t" + toString(output)); m->mothurOutEndLine();
					process(list);
				
					
					if (m->getControl_pressed()) { outputTypes.clear(); if (ct != NULL) { delete ct; }  if (groupMap != NULL) { delete groupMap; } for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);  } delete input; delete list; return 0; }
										
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
		if (needToRun )  {
			if (list != NULL) {	delete list;	}
			list = input->getListVector(lastLabel);
            string output = toString(list->size());  if (basis == "sequence") { output = toString(list->getNumSeqs()); }
			m->mothurOut(list->getLabel() + "\t" + toString(output)); m->mothurOutEndLine();
			
			process(list);
			delete list;
			
			if (m->getControl_pressed()) { outputTypes.clear();  if (ct != NULL) { delete ct; } if (groupMap != NULL) { delete groupMap; } for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);  } delete input; delete list; return 0; }
		}
		
		delete input;  
        if (groupMap != NULL) { delete groupMap; }
        if (ct != NULL) { delete ct; }
				
		if (m->getControl_pressed()) { outputTypes.clear(); for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);  } return 0; }
		
        //set constaxonomy file as new current constaxonomyfile
        string currentName = "";
        itTypes = outputTypes.find("constaxonomy");
        if (itTypes != outputTypes.end()) { if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setConsTaxonomyFile(currentName); } }
        
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifyOtuCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> ClassifyOtuCommand::findConsensusTaxonomy(vector<string> names, int& size, string& conTax, string group) {
	try{
		conTax = "";
		vector<string> allNames;
		map<string, string>::iterator it;
		map<string, string>::iterator it2;

		//create a tree containing sequences from this bin
		PhyloTree* phylo = new PhyloTree();
		
		size = 0;
		for (int i = 0; i < names.size(); i++) {
            
            if (group != "") { //no need to check for name file, names already added in previous step
                //is this sequence in the taxonomy file - look for repSeqName since we are assuming the taxonomy file is unique
                it = taxMap.find(names[i]);
                
                if (it == taxMap.end()) { //this name is not in taxonomy file, skip it
                    m->mothurOut("[WARNING]: " + names[i] + " is not in your taxonomy file.  I will not include it in the consensus."); m->mothurOutEndLine();
                }else{
                    if (countfile != "") {
                        int numDups = ct->getGroupCount(names[i], group);
                        for (int j = 0; j < numDups; j++) {  phylo->addSeqToTree(names[i], it->second);  }
                        size += numDups;
                    }else{
                        //add seq to tree
                        phylo->addSeqToTree(names[i], it->second);
                        size++;
                    }
                    allNames.push_back(names[i]);
                }

            }else {
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
                        m->mothurOut("[WARNING]: " + names[i] + " is not in your taxonomy file.  I will not include it in the consensus."); m->mothurOutEndLine();
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
            }
			
			if (m->getControl_pressed()) { delete phylo; return allNames; }
			
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
        
        if (conTax == "") {  conTax = "unknown;";  }
        
		if (myLevel != phylo->getMaxLevel()) {  conTax = util.addUnclassifieds(conTax, phylo->getMaxLevel(), probs);  }
		
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
		if (outputDir == "") { outputDir += util.hasPath(listfile); }
				
		ofstream out;
        map<string, string> variables; 
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(listfile));
        variables["[distance]"] = processList->getLabel();
		string outputFile = getOutputFileName("constaxonomy", variables);
		util.openOutputFile(outputFile, out);
		outputNames.push_back(outputFile); outputTypes["constaxonomy"].push_back(outputFile);
		
		ofstream outSum;
		string outputSumFile = getOutputFileName("taxsummary", variables);
		util.openOutputFile(outputSumFile, outSum);
		outputNames.push_back(outputSumFile); outputTypes["taxsummary"].push_back(outputSumFile);
		
		out << "OTU\tSize\tTaxonomy" << endl;
		
		PhyloSummary* taxaSum;
        if (countfile != "") { taxaSum = new PhyloSummary(ct,relabund, printlevel); }
        else { taxaSum = new PhyloSummary(groupMap,relabund, printlevel); }
        
        vector<string> outs;
        vector<PhyloSummary*> taxaSums;
        map<string, int> groupIndex;
        if (persample) {
            for (int i = 0; i < groups.size(); i++) {
                groupIndex[groups[i]] = i;
                variables["[distance]"] = processList->getLabel() + "." + groups[i];
                string outputFile = getOutputFileName("constaxonomy", variables);
                ofstream temp;
                util.openOutputFile(outputFile, temp);
                outs.push_back(outputFile);
                temp << "OTU\tSize\tTaxonomy" << endl;
                outputNames.push_back(outputFile); outputTypes["constaxonomy"].push_back(outputFile);
                
                PhyloSummary* taxaSumt;
                if (countfile != "") { taxaSumt = new PhyloSummary(ct, relabund, printlevel);
                }else { taxaSumt = new PhyloSummary(groupMap,relabund, printlevel); }
                taxaSums.push_back(taxaSumt);
            }
        }
        
		//for each bin in the list vector
        string snumBins = toString(processList->getNumBins());
        vector<string> binLabels = processList->getLabels();
		for (int i = 0; i < processList->getNumBins(); i++) {
			
			if (m->getControl_pressed()) { break; }
			
			vector<string> names;
            string binnames = processList->get(i);
            vector<string> thisNames;
            util.splitAtComma(binnames, thisNames);
            
			names = findConsensusTaxonomy(thisNames, size, conTax, "");
		
			if (m->getControl_pressed()) { break; }

			out << binLabels[i] << '\t' << size << '\t' << conTax << endl;
			
			string noConfidenceConTax = conTax;
			util.removeConfidences(noConfidenceConTax);
			
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
                    vector<string> theseNames = findConsensusTaxonomy(itParsed->second, size, conTax, itParsed->first);
                    
                    if (m->getControl_pressed()) { break; }
                    
                    ofstream out; util.openOutputFileAppend(outs[groupIndex[itParsed->first]], out);
                    out << binLabels[i] << '\t' << size << '\t' << conTax << endl;
                    out.close();
                    
                    string noConfidenceConTax = conTax;
                    util.removeConfidences(noConfidenceConTax);
                    
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
		taxaSum->print(outSum, output);
		outSum.close();
        
        if (persample) {
            for (int i = 0; i < groups.size(); i++) {
                ofstream outSums;
                variables["[distance]"] = processList->getLabel() + "." + groups[i];
                string outputSumFile = getOutputFileName("taxsummary", variables);
                util.openOutputFile(outputSumFile, outSums);
                outputNames.push_back(outputSumFile); outputTypes["taxsummary"].push_back(outputSumFile);
                taxaSums[i]->print(outSums, output);
                outSums.close();
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
int ClassifyOtuCommand::processTaxMap() {
    try{
        
        for (map<string, string>::iterator it = taxMap.begin(); it != taxMap.end(); it++) {
            
            if (m->getControl_pressed()) { break; }
            
            vector<string> taxons;
            string tax = it->second;
            int taxLength = tax.length();
            string taxon = "";
            int spot = 0;
            
            for(int i=0;i<taxLength;i++){
                
                
                if(tax[i] == ';'){
                    
                    int openParen = taxon.find_last_of('(');
                    int closeParen = taxon.find_last_of(')');
                    
                    string newtaxon, confidence;
                    if ((openParen != string::npos) && (closeParen != string::npos)) {
                        string confidenceScore = taxon.substr(openParen+1, (closeParen-(openParen+1)));
                        if (util.isNumeric1(confidenceScore)) {  //its a confidence
                            newtaxon = taxon.substr(0, openParen); //rip off confidence
                            confidence = taxon.substr((openParen+1), (closeParen-openParen-1));
                        }else { //its part of the taxon
                            newtaxon = taxon;
                            confidence = "0";
                        }
                    }else{
                        newtaxon = taxon;
                        confidence = "-1";
                    }
                    float con = 0;
                    convert(confidence, con);
                    
                    if (con == -1) { i += taxLength; } //not a confidence score, no confidence scores on this taxonomy
                    else if ( con < threshold)  { spot = i; break; } //below threshold, set all to unclassified
                    else {} //acceptable, move on
                    taxons.push_back(taxon);
                    
                    taxon = "";
                }
                else{
                    taxon += tax[i];
                }
                
            }
            
            if (spot != 0) {
                string newTax = "";
                for (int i = 0; i < taxons.size(); i++) {  newTax += taxons[i] + ";";  }
                //for (int i = spot; i < taxLength; i++) {
                    //if(tax[i] == ';'){   newTax += "unclassified;"; }
                    util.removeConfidences(newTax);
                    it->second = newTax;
                //}
            }else { util.removeConfidences(tax); it->second = tax; } //leave tax alone
        }
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "ClassifyOtuCommand", "processTaxMap");
        exit(1);
    }
}
//**********************************************************************************************************************


