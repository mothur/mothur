/*
 *  indicatorcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 11/12/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "indicatorcommand.h"
#include "sharedutilities.h"


//**********************************************************************************************************************
vector<string> IndicatorCommand::setParameters(){	
	try {
		CommandParameter piters("iters", "Number", "", "1000", "", "", "","",false,false); parameters.push_back(piters);
		CommandParameter pdesign("design", "InputTypes", "", "", "TreeDesign", "TreeDesign", "none","summary",false,false,true); parameters.push_back(pdesign);
		CommandParameter pshared("shared", "InputTypes", "", "", "SharedRel", "SharedRel", "none","summary",false,false,true); parameters.push_back(pshared);	
		CommandParameter prelabund("relabund", "InputTypes", "", "", "SharedRel", "SharedRel", "none","summary",false,false); parameters.push_back(prelabund);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter ptree("tree", "InputTypes", "", "", "TreeDesign", "TreeDesign", "none","tree-summary",false,false,true); parameters.push_back(ptree);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false); parameters.push_back(pprocessors);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string IndicatorCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The indicator command can be run in 3 ways: with a shared or relabund file and a design file, or with a shared or relabund file and a tree file, or with a shared or relabund file, tree file and design file. \n";
		helpString += "The indicator command outputs a .indicator.summary file and a .indicator.tre if a tree is given. \n";
		helpString += "The new tree contains labels at each internal node.  The label is the node number so you can relate the tree to the summary file.\n";
		helpString += "The summary file lists the indicator value for each OTU for each node.\n";
		helpString += "The indicator command parameters are tree, groups, shared, relabund, design and label. \n";
		helpString += "The design parameter allows you to relate the tree to the shared or relabund file, if your tree contains the grouping names, or if no tree is provided to group your groups into groupings.\n";			
		helpString += "The groups parameter allows you to specify which of the groups in your shared or relabund you would like analyzed, or if you provide a design file the groups in your design file.  The groups may be entered separated by dashes.\n";
		helpString += "The label parameter indicates at what distance your tree relates to the shared or relabund.\n";
		helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n";
		helpString += "The iters parameter allows you to set number of randomization for the P value.  The default is 1000.";
		helpString += "The indicator command should be used in the following format: indicator(tree=test.tre, shared=test.shared, label=0.03)\n";
		helpString += "Note: No spaces between parameter labels (i.e. tree), '=' and parameters (i.e.yourTreefile).\n"; 
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string IndicatorCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "tree") {  pattern = "[filename],indicator.tre"; } 
        else if (type == "summary") {  pattern = "[filename],indicator.summary"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "IndicatorCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
IndicatorCommand::IndicatorCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["tree"] = tempOutNames;
		outputTypes["summary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "IndicatorCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
IndicatorCommand::IndicatorCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
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
			
			m->runParse = true;
			m->clearGroups();
			m->clearAllGroups();
			m->Treenames.clear();
			
			vector<string> tempOutNames;
			outputTypes["tree"] = tempOutNames;
			outputTypes["summary"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("tree");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["tree"] = inputDir + it->second;		}
				}
				
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
				
				it = parameters.find("design");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["design"] = inputDir + it->second;		}
				}
			}
			
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}

			//check for required parameters
			treefile = validParameter.validFile(parameters, "tree", true);
			if (treefile == "not open") { treefile = ""; abort = true; }
			else if (treefile == "not found") { treefile = ""; 	}		
			else { m->setTreeFile(treefile); }	
			
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { abort = true; }
			else if (sharedfile == "not found") { sharedfile = ""; }
			else { inputFileName = sharedfile; m->setSharedFile(sharedfile); }
			
			relabundfile = validParameter.validFile(parameters, "relabund", true);
			if (relabundfile == "not open") { abort = true; }
			else if (relabundfile == "not found") { relabundfile = ""; }
			else { inputFileName = relabundfile; m->setRelAbundFile(relabundfile); }
			
			designfile = validParameter.validFile(parameters, "design", true);
			if (designfile == "not open") { designfile = ""; abort = true; }
			else if (designfile == "not found") { designfile = ""; }
			else { m->setDesignFile(designfile); }
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = "";  Groups.push_back("all"); }
			else { m->splitAtDash(groups, Groups);	}			
			m->setGroups(Groups);
			
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; m->mothurOut("You did not provide a label, I will use the first label in your inputfile."); m->mothurOutEndLine(); label=""; }	
			
			string temp = validParameter.validFile(parameters, "iters", false);		if (temp == "not found") { temp = "1000"; }
			m->mothurConvert(temp, iters); 
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors); 
			
			if ((relabundfile == "") && (sharedfile == "")) { 
				//is there are current file available for either of these?
				//give priority to shared, then relabund
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") {  inputFileName = sharedfile; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 
					relabundfile = m->getRelAbundFile(); 
					if (relabundfile != "") {  inputFileName = relabundfile; m->mothurOut("Using " + relabundfile + " as input file for the relabund parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("No valid current files. You must provide a shared or relabund."); m->mothurOutEndLine(); 
						abort = true;
					}
				}
			}
			
			
			if ((designfile == "") && (treefile == "")) { 
				treefile = m->getTreeFile(); 
				if (treefile != "") {  m->mothurOut("Using " + treefile + " as input file for the tree parameter."); m->mothurOutEndLine(); }
				else { 
					designfile = m->getDesignFile(); 
					if (designfile != "") {  m->mothurOut("Using " + designfile + " as input file for the design parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("[ERROR]: You must provide either a tree or design file."); m->mothurOutEndLine(); abort = true; 
					}
				}
			}
			
			if ((relabundfile != "") && (sharedfile != "")) { m->mothurOut("[ERROR]: You may not use both a shared and relabund file."); m->mothurOutEndLine(); abort = true;  }
			
			
		}
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "IndicatorCommand");		
		exit(1);
	}
}
//**********************************************************************************************************************

int IndicatorCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		cout.setf(ios::fixed, ios::floatfield); cout.setf(ios::showpoint);
		
		int start = time(NULL);
	
		//read designfile if given and set up groups for read of sharedfiles
		if (designfile != "") {
			designMap = new DesignMap(designfile);
			
			//fill Groups - checks for "all" and for any typo groups
			SharedUtil util;
			vector<string> nameGroups = designMap->getCategory();
			util.setGroups(Groups, nameGroups);
			
			vector<string> namesSeqs = designMap->getNamesGroups(Groups);
			m->setGroups(namesSeqs);
		}
	
		/***************************************************/
		// use smart distancing to get right sharedRabund  //
		/***************************************************/
		if (sharedfile != "") {  
			getShared(); 
            if (m->control_pressed) { if (designfile != "") { delete designMap; } delete lookup;   return 0; }
			if (lookup == NULL) { m->mothurOut("[ERROR] reading shared file."); m->mothurOutEndLine(); return 0; }
		}else { 
			getSharedFloat(); 
			if (m->control_pressed) {  if (designfile != "") { delete designMap; } delete lookupFloat; return 0; }
			if (lookupFloat == NULL) { m->mothurOut("[ERROR] reading relabund file."); m->mothurOutEndLine(); return 0; }
		}
		
		//reset groups if needed
		if (designfile != "") { m->setGroups(Groups); }
			
		/***************************************************/
		//    reading tree info							   //
		/***************************************************/
		if (treefile != "") {
			string groupfile = ""; 
			m->setTreeFile(treefile);
			Tree* tree = new Tree(treefile); delete tree;  //extracts names from tree to make faked out groupmap
			ct = new CountTable();
			bool mismatch = false;
            
            set<string> nameMap;
            map<string, string> groupMap;
            set<string> gps;
            for (int i = 0; i < m->Treenames.size(); i++) { 
                nameMap.insert(m->Treenames[i]); 
                //sanity check - is this a group that is not in the sharedfile?
                if (i == 0) { gps.insert("Group1"); }
				if (designfile == "") {
					if (!(m->inUsersGroups(m->Treenames[i], m->getAllGroups()))) {
						m->mothurOut("[ERROR]: " + m->Treenames[i] + " is not a group in your shared or relabund file."); m->mothurOutEndLine();
						mismatch = true;
					}
					groupMap[m->Treenames[i]] = "Group1"; 
				}else{
					vector<string> myGroups; myGroups.push_back(m->Treenames[i]);
					vector<string> myNames = designMap->getNamesGroups(myGroups);
					
					for(int k = 0; k < myNames.size(); k++) {
						if (!(m->inUsersGroups(myNames[k], m->getAllGroups()))) {
							m->mothurOut("[ERROR]: " + myNames[k] + " is not a group in your shared or relabund file."); m->mothurOutEndLine();
							mismatch = true;
						}
					}
					groupMap[m->Treenames[i]] = "Group1";
				}
            }
            ct->createTable(nameMap, groupMap, gps);
			
			if ((designfile != "") && (m->Treenames.size() != Groups.size())) { cout << Groups.size() << '\t' << m->Treenames.size() << endl; m->mothurOut("[ERROR]: You design file does not match your tree, aborting."); m->mothurOutEndLine(); mismatch = true; }
					
			if (mismatch) { //cleanup and exit
				if (designfile != "") { delete designMap; }
				if (sharedfile != "") {  delete lookup; }
				else { delete lookupFloat; }
				delete ct;
				return 0;
			}
		 
			read = new ReadNewickTree(treefile);
			int readOk = read->read(ct); 
			
			if (readOk != 0) { m->mothurOut("Read Terminated."); m->mothurOutEndLine(); delete ct; delete read; return 0; }
			
			vector<Tree*> T = read->getTrees();
			
			delete read;
			
			if (m->control_pressed) { 
				if (designfile != "") { delete designMap; }
				if (sharedfile != "") {  delete lookup; }
				else { delete lookupFloat; }
				for (int i = 0; i < T.size(); i++) {  delete T[i];  }  delete ct; return 0; 
			}
            
			T[0]->assembleTree();
					
			/***************************************************/
			//    create ouptut tree - respecting pickedGroups //
			/***************************************************/
			Tree* outputTree = new Tree(m->getNumGroups(), ct); 
			
			outputTree->getSubTree(T[0], m->getGroups());
			outputTree->assembleTree();
				
			//no longer need original tree, we have output tree to use and label
			for (int i = 0; i < T.size(); i++) {  delete T[i];  } 
			
			if (m->control_pressed) { 
				if (designfile != "") { delete designMap; }
                if (sharedfile != "") {  delete lookup; }
                else { delete lookupFloat; }
				delete outputTree; delete ct;  return 0; 
			}
			
			/***************************************************/
			//		get indicator species values			   //
			/***************************************************/
			GetIndicatorSpecies(outputTree);
			delete outputTree; delete ct;
			
		}else { //run with design file only
			//get indicator species
			GetIndicatorSpecies();
		}
		
		if (designfile != "") { delete designMap; }
        if (sharedfile != "") {  delete lookup; }
        else { delete lookupFloat; }
        
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);	} return 0; }
		
		//set tree file as new current treefile
		if (treefile != "") {
			string current = "";
			itTypes = outputTypes.find("tree");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setTreeFile(current); }
			}
		}
		
		m->mothurOutEndLine(); m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to find the indicator species."); m->mothurOutEndLine();
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
				
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "execute");	
		exit(1);
	}
}
//**********************************************************************************************************************
//divide shared or relabund file by groupings in the design file
//report all otu values to file
int IndicatorCommand::GetIndicatorSpecies(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(inputFileName);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(inputFileName));
		string outputFileName = getOutputFileName("summary", variables);
		outputNames.push_back(outputFileName); outputTypes["summary"].push_back(outputFileName);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		m->mothurOutEndLine(); m->mothurOut("Species\tIndicator_Groups\tIndicatorValue\tpValue\n");
		
		int numBins = 0;
		if (sharedfile != "") { numBins = lookup->getNumBins(); }
		else { numBins = lookupFloat->getNumBins(); }
		
		if (m->control_pressed) { out.close(); return 0; }
			
		/*****************************************************/
		//create vectors containing rabund info				 //
		/*****************************************************/
			
		vector<float> indicatorValues; //size of numBins
		vector<float> pValues;
        vector<string> indicatorGroups;
		map< vector<int>, vector<int> > randomGroupingsMap; //maps location in groupings to location in groupings, ie, [0][0] -> [1][2]. This is so we don't have to actually move the sharedRabundVectors.
			
		if (sharedfile != "") {
			vector< vector<SharedRAbundVector*> > groupings;
            vector< vector<string> > groupingNames;
			set<string> groupsAlreadyAdded;
			vector<SharedRAbundVector*> subset;
            vector<string> subsetNames;
            
            vector<SharedRAbundVector*> data = lookup->getSharedRAbundVectors();
            vector<string> dataGroupNames = lookup->getNamesGroups();
			
			//for each grouping
			for (int i = 0; i < (designMap->getCategory()).size(); i++) {
				
				for (int k = 0; k < data.size(); k++) {
					//are you from this grouping?
					if (designMap->get(dataGroupNames[k]) == (designMap->getCategory())[i]) {
						subset.push_back(data[k]);
                        subsetNames.push_back(dataGroupNames[k]);
						groupsAlreadyAdded.insert(dataGroupNames[k]);
					}
				}
				if (subset.size() != 0) { groupings.push_back(subset); groupingNames.push_back(subsetNames); }
				subset.clear();
			}
				
			if (groupsAlreadyAdded.size() != data.size()) {  m->mothurOut("[ERROR]: could not make proper groupings."); m->mothurOutEndLine(); }
				
			indicatorValues = getValues(groupings, groupingNames, indicatorGroups, randomGroupingsMap);
			
			pValues = getPValues(groupings, groupingNames, lookup->size(), indicatorValues);
		}else {
			vector< vector<SharedRAbundFloatVector*> > groupings;
            vector< vector<string> > groupingNames;
			set<string> groupsAlreadyAdded;
			vector<SharedRAbundFloatVector*> subset;
            vector<string> subsetNames;
            
            vector<SharedRAbundFloatVector*> data = lookupFloat->getSharedRAbundFloatVectors();
            vector<string> dataGroupNames = lookupFloat->getNamesGroups();
			
			//for each grouping
			for (int i = 0; i < (designMap->getCategory()).size(); i++) {
				for (int k = 0; k < data.size(); k++) {
					//are you from this grouping?
					if (designMap->get(dataGroupNames[k]) == (designMap->getCategory())[i]) {
						subset.push_back(data[k]);
                        subsetNames.push_back(dataGroupNames[k]);
						groupsAlreadyAdded.insert(dataGroupNames[k]);
					}
				}
				if (subset.size() != 0) { groupings.push_back(subset); groupingNames.push_back(subsetNames); }
				subset.clear();
			}
			
			if (groupsAlreadyAdded.size() != data.size()) {  m->mothurOut("[ERROR]: could not make proper groupings."); m->mothurOutEndLine(); }
			
			indicatorValues = getValues(groupings, groupingNames, indicatorGroups, randomGroupingsMap);
			
			pValues = getPValues(groupings, groupingNames, lookupFloat->size(), indicatorValues);
		}
			
		if (m->control_pressed) { out.close(); return 0; }
			
			
		/******************************************************/
		//output indicator values to table form               //
		/*****************************************************/
		out << "OTU\tIndicator_Groups\tIndicator_Value\tpValue" << endl;
		for (int j = 0; j < indicatorValues.size(); j++) {
				
			if (m->control_pressed) { out.close(); return 0; }
			
			out << m->currentSharedBinLabels[j] << '\t' << indicatorGroups[j] << '\t' << indicatorValues[j] << '\t';
			
			if (pValues[j] > (1/(float)iters)) { out << pValues[j] << endl; } 
			else { out << "<" << (1/(float)iters) << endl; }
			
			if (pValues[j] <= 0.05) {
				cout << m->currentSharedBinLabels[j] << '\t' << indicatorGroups[j] << '\t' << indicatorValues[j]  << '\t';
				string pValueString = "<" + toString((1/(float)iters)); 
				if (pValues[j] > (1/(float)iters)) { pValueString = toString(pValues[j]); cout << pValues[j];} 
				else { cout << "<" << (1/(float)iters); }
				m->mothurOutJustToLog(m->currentSharedBinLabels[j] + "\t" + indicatorGroups[j] + "\t" + toString(indicatorValues[j]) + "\t" + pValueString);
				m->mothurOutEndLine(); 
			}
		}
		
		out.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "GetIndicatorSpecies");	
		exit(1);
	}
}
//**********************************************************************************************************************
//traverse tree finding indicator species values for each otu at each node
//label node with otu number that has highest indicator value
//report all otu values to file
int IndicatorCommand::GetIndicatorSpecies(Tree*& T){
	try {
		
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(inputFileName);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(inputFileName));
		string outputFileName = getOutputFileName("summary",variables);
		outputNames.push_back(outputFileName); outputTypes["summary"].push_back(outputFileName);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		int numBins = 0;
		if (sharedfile != "") { numBins = lookup->getNumBins(); }
		else { numBins = lookupFloat->getNumBins(); }
		
		//print headings
		out << "TreeNode\t";
		for (int i = 0; i < numBins; i++) { out << m->currentSharedBinLabels[i] << "_IndGroups" << '\t' << m->currentSharedBinLabels[i] << "_IndValue" << '\t' << "pValue" << '\t'; }
		out << endl;
		
		m->mothurOutEndLine(); m->mothurOut("Node\tSpecies\tIndicator_Groups\tIndicatorValue\tpValue\n");
		
		string treeOutputDir = outputDir;
		if (outputDir == "") {  treeOutputDir += m->hasPath(treefile);  }
        variables["[filename]"] = treeOutputDir + m->getRootName(m->getSimpleName(treefile));
		string outputTreeFileName = getOutputFileName("tree", variables);
		
		
		//create a map from tree node index to names of descendants, save time later to know which sharedRabund you need
		map<int, set<string> > nodeToDescendants;
		map<int, set<int> > descendantNodes;
		for (int i = 0; i < T->getNumNodes(); i++) { if (m->control_pressed) { return 0; }
			nodeToDescendants[i] = getDescendantList(T, i, nodeToDescendants, descendantNodes);
		}
		
		//you need the distances to leaf to decide grouping below
		//this will also set branch lengths if the tree does not include them
		map<int, float> distToRoot = getDistToRoot(T);
			
		//for each node
		for (int i = T->getNumLeaves(); i < T->getNumNodes(); i++) {
			
			if (m->control_pressed) { out.close(); return 0; }
			
			/*****************************************************/
			//create vectors containing rabund info				 //
			/*****************************************************/
			
			vector<float> indicatorValues; //size of numBins
			vector<float> pValues;
            vector<string> indicatorGroups;
			map< vector<int>, vector<int> > randomGroupingsMap; //maps location in groupings to location in groupings, ie, [0][0] -> [1][2]. This is so we don't have to actually move the sharedRabundVectors.
			
			if (sharedfile != "") {
				vector< vector<SharedRAbundVector*> > groupings;
                vector< vector<string> > groupingNames;
				
				//get nodes that will be a valid grouping
				//you are valid if you are not one of my descendants
				//AND your distToRoot is >= mine
				//AND you were not added as part of a larger grouping. Largest nodes are added first.
				
				set<string> groupsAlreadyAdded;
				//create a grouping with my grouping
				vector<SharedRAbundVector*> subset;
                vector<string> subsetNames;
				int count = 0;
				int doneCount = nodeToDescendants[i].size();
                vector<SharedRAbundVector*> data = lookupFloat->getSharedRAbundVectors();
                vector<string> dataGroupNames = lookupFloat->getNamesGroups();

				for (int k = 0; k < data.size(); k++) {
					//is this descendant of i
					if ((nodeToDescendants[i].count(dataGroupNames[k]) != 0)) {
						subset.push_back(data[k]);
                        subsetNames.push_back(dataGroupNames[k]);
						groupsAlreadyAdded.insert(dataGroupNames[k]);
						count++;
					}
					if (count == doneCount) { break; } //quit once you get the rabunds for this grouping
				}
				if (subset.size() != 0) { groupings.push_back(subset); groupingNames.push_back(subsetNames); }
				
				for (int j = (T->getNumNodes()-1); j >= 0; j--) {
	
					if ((descendantNodes[i].count(j) == 0) && (distToRoot[j] >= distToRoot[i])) { 
						vector<SharedRAbundVector*> subset;
                        vector<string> subsetNames;
						int count = 0;
						int doneCount = nodeToDescendants[j].size();
						for (int k = 0; k < data.size(); k++) {
							//is this descendant of j, and we didn't already add this as part of a larger grouping
							if ((nodeToDescendants[j].count(dataGroupNames[k]) != 0) && (groupsAlreadyAdded.count(dataGroupNames[k]) == 0)) {
								subset.push_back(data[k]);
                                subsetNames.push_back(dataGroupNames[k]);
								groupsAlreadyAdded.insert(dataGroupNames[k]);
								count++;
							}
							if (count == doneCount) { break; } //quit once you get the rabunds for this grouping
						}
						
						//if subset.size == 0 then the node was added as part of a larger grouping
                        if (subset.size() != 0) { groupings.push_back(subset); groupingNames.push_back(subsetNames); }
					}
				}
				
				if (groupsAlreadyAdded.size() != data.size()) {  m->mothurOut("[ERROR]: could not make proper groupings."); m->mothurOutEndLine(); }
								
				indicatorValues = getValues(groupings, groupingNames, indicatorGroups, randomGroupingsMap);
				
				pValues = getPValues(groupings, groupingNames, data.size(), indicatorValues);
			}else {
				vector< vector<SharedRAbundFloatVector*> > groupings;
                vector< vector<string> > groupingNames;
				
				//get nodes that will be a valid grouping
				//you are valid if you are not one of my descendants
				//AND your distToRoot is >= mine
				//AND you were not added as part of a larger grouping. Largest nodes are added first.
				
				set<string> groupsAlreadyAdded;
				//create a grouping with my grouping
				vector<SharedRAbundFloatVector*> subset;
                vector<string> subsetNames;
				int count = 0;
				int doneCount = nodeToDescendants[i].size();
                
                vector<SharedRAbundFloatVector*> data = lookupFloat->getSharedRAbundFloatVectors();
                vector<string> dataGroupNames = lookupFloat->getNamesGroups();
                
				for (int k = 0; k < data.size(); k++) {
					//is this descendant of i
					if ((nodeToDescendants[i].count(dataGroupNames[k]) != 0)) {
						subset.push_back(data[k]);
                        subsetNames.push_back(dataGroupNames[k]);
						groupsAlreadyAdded.insert(dataGroupNames[k]);
						count++;
					}
					if (count == doneCount) { break; } //quit once you get the rabunds for this grouping
				}
				if (subset.size() != 0) { groupings.push_back(subset); groupingNames.push_back(subsetNames); }
				
				for (int j = (T->getNumNodes()-1); j >= 0; j--) {
					if ((descendantNodes[i].count(j) == 0) && (distToRoot[j] >= distToRoot[i])) {
						vector<SharedRAbundFloatVector*> subset;
                        vector<string> subsetNames;
						int count = 0;
						int doneCount = nodeToDescendants[j].size();
						for (int k = 0; k < data.size(); k++) {
							//is this descendant of j, and we didn't already add this as part of a larger grouping
							if ((nodeToDescendants[j].count(dataGroupNames[k]) != 0) && (groupsAlreadyAdded.count(dataGroupNames[k]) == 0)) {
								subset.push_back(data[k]);
                                subsetNames.push_back(dataGroupNames[k]);
								groupsAlreadyAdded.insert(dataGroupNames[k]);
								count++;
							}
							if (count == doneCount) { break; } //quit once you get the rabunds for this grouping
						}
						
						//if subset.size == 0 then the node was added as part of a larger grouping
						if (subset.size() != 0) { groupings.push_back(subset); groupingNames.push_back(subsetNames); }
					}
				}
				
				if (groupsAlreadyAdded.size() != data.size()) { m->mothurOut("[ERROR]: could not make proper groupings."); m->mothurOutEndLine(); }
				
				indicatorValues = getValues(groupings, groupingNames, indicatorGroups, randomGroupingsMap);
				
				pValues = getPValues(groupings, groupingNames, data.size(), indicatorValues);
			}
			
			if (m->control_pressed) { out.close(); return 0; }
			
			
			/******************************************************/
			//output indicator values to table form + label tree  //
			/*****************************************************/
			out << (i+1);
			for (int j = 0; j < indicatorValues.size(); j++) {
				
				if (m->control_pressed) { out.close(); return 0; }
				
				if (pValues[j] < (1/(float)iters)) {
					out  << '\t' << indicatorGroups[j] << '\t' << indicatorValues[j] << '\t' << '<' << (1/(float)iters);
				}else {
					out  << '\t' << indicatorGroups[j] << '\t' << indicatorValues[j] << '\t' << pValues[j];
				}
				
				if (pValues[j] <= 0.05) {
					cout << i+1 << '\t' << m->currentSharedBinLabels[j] << '\t' << indicatorGroups[j] << '\t' << indicatorValues[j];
					string pValueString = "\t<" + toString((1/(float)iters));
					if (pValues[j] > (1/(float)iters)) { pValueString = toString('\t' + pValues[j]); cout << '\t' << pValues[j];}
					else { cout << "\t<" << (1/(float)iters); }
					m->mothurOutJustToLog(toString(i) + "\t" + m->currentSharedBinLabels[j] + "\t" + indicatorGroups[j] + "\t" + toString(indicatorValues[j]) + "\t" + pValueString);
					m->mothurOutEndLine(); 
				}
			}
			out << endl;
			
			T->tree[i].setLabel(toString(i+1));
		}
		out.close();
	
		ofstream outTree;
		m->openOutputFile(outputTreeFileName, outTree);
		outputNames.push_back(outputTreeFileName); outputTypes["tree"].push_back(outputTreeFileName);
		
		T->print(outTree, "both");
		outTree.close();
	
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "GetIndicatorSpecies");	
		exit(1);
	}
}
//**********************************************************************************************************************
vector<float> IndicatorCommand::getValues(vector< vector<SharedRAbundFloatVector*> >& groupings, vector< vector<string> >& groupingNames, vector<string>& indicatorGroupings, map< vector<int>, vector<int> > groupingsMap){
	try {
		vector<float> values;
		map< vector<int>, vector<int> >::iterator it;
        
        indicatorGroupings.clear();
        
        //create grouping strings
        vector<string> groupingsGroups;
        for (int j = 0; j < groupings.size(); j++) {
            string tempGrouping = "";
            for (int k = 0; k < groupings[j].size()-1; k++) { 
                tempGrouping += groupingNames[j][k] + "-";
            }
            tempGrouping += groupingNames[j][groupingNames[j].size()-1];
            groupingsGroups.push_back(tempGrouping);
        }
        
        
		//for each otu
		for (int i = 0; i < groupings[0][0]->getNumBins(); i++) {
			
			if (m->control_pressed) { return values; }
			
			vector<float> terms; 
			float AijDenominator = 0.0;
			vector<float> Bij;
			
			//get overall abundance of each grouping
			for (int j = 0; j < groupings.size(); j++) {
				
				float totalAbund = 0;
				int numNotZero = 0;
				for (int k = 0; k < groupings[j].size(); k++) { 
					vector<int> temp; temp.push_back(j); temp.push_back(k);
					it = groupingsMap.find(temp);
					
					if (it == groupingsMap.end()) { //this one didnt get moved
						totalAbund += groupings[j][k]->get(i);
						if (groupings[j][k]->get(i) != 0.0) { numNotZero++; }
					}else {
						totalAbund += groupings[(it->second)[0]][(it->second)[1]]->get(i);
						if (groupings[(it->second)[0]][(it->second)[1]]->get(i) != 0.0) { numNotZero++; }
					}
					
				}
				
				//mean abundance
				float Aij = (totalAbund / (float) groupings[j].size());
				terms.push_back(Aij);
				
				//percentage of sites represented
				Bij.push_back(numNotZero / (float) groupings[j].size());
				
				AijDenominator += Aij;
			}
			
			float maxIndVal = 0.0;
            string maxGrouping = "";
			for (int j = 0; j < terms.size(); j++) { 
				float thisAij = (terms[j] / AijDenominator); //relative abundance
				float thisValue = thisAij * Bij[j] * 100.0;
				
				//save largest
				if (thisValue > maxIndVal) { maxIndVal = thisValue;  maxGrouping = groupingsGroups[j]; }
			}
			
			values.push_back(maxIndVal);
            indicatorGroupings.push_back(maxGrouping);
		}
		
		return values;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "getValues");	
		exit(1);
	}
}
//**********************************************************************************************************************
//same as above, just data type difference
vector<float> IndicatorCommand::getValues(vector< vector<SharedRAbundVector*> >& groupings, vector< vector<string> >& groupingNames,vector<string>& indicatorGroupings, map< vector<int>, vector<int> > groupingsMap){
	try {
		vector<float> values;
		map< vector<int>, vector<int> >::iterator it;
        
        indicatorGroupings.clear();
        
        //create grouping strings
        vector<string> groupingsGroups;
        for (int j = 0; j < groupings.size(); j++) {
            string tempGrouping = "";
            for (int k = 0; k < groupings[j].size()-1; k++) { 
                tempGrouping += groupingNames[j][k] + "-";
            }
            tempGrouping += groupingNames[j][groupingNames[j].size()-1];
            groupingsGroups.push_back(tempGrouping);
        }

		
		//for each otu
		for (int i = 0; i < groupings[0][0]->getNumBins(); i++) {
			vector<float> terms; 
			float AijDenominator = 0.0;
			vector<float> Bij;
			//get overall abundance of each grouping
			for (int j = 0; j < groupings.size(); j++) {
	
				int totalAbund = 0.0;
				int numNotZero = 0;
				for (int k = 0; k < groupings[j].size(); k++) { 
					vector<int> temp; temp.push_back(j); temp.push_back(k);
					it = groupingsMap.find(temp);
					
					if (it == groupingsMap.end()) { //this one didnt get moved
						totalAbund += groupings[j][k]->get(i);
						if (groupings[j][k]->get(i) != 0.0) { numNotZero++; }
					}else {
						totalAbund += groupings[(it->second)[0]][(it->second)[1]]->get(i);
						if (groupings[(it->second)[0]][(it->second)[1]]->get(i) != 0.0) { numNotZero++; }
					}
					
				}
				
				//mean abundance	
				float Aij = (totalAbund / (float) groupings[j].size());
				terms.push_back(Aij);
				
				//percentage of sites represented
				Bij.push_back(numNotZero / (float) groupings[j].size());
				
				AijDenominator += Aij;
			}
			
			float maxIndVal = 0.0;
            string maxGrouping = "";
			for (int j = 0; j < terms.size(); j++) { 
				float thisAij = (terms[j] / AijDenominator); //relative abundance
				float thisValue = thisAij * Bij[j] * 100.0;
					
				//save largest
				if (thisValue > maxIndVal) { maxIndVal = thisValue;  maxGrouping = groupingsGroups[j]; }
			}
			
			values.push_back(maxIndVal);
            indicatorGroupings.push_back(maxGrouping);
		}
		
		return values;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "getValues");	
		exit(1);
	}
}
//**********************************************************************************************************************
//you need the distances to root to decide groupings
//this will also set branch lengths if the tree does not include them
map<int, float> IndicatorCommand::getDistToRoot(Tree*& T){
	try {
		map<int, float> dists;
		
		bool hasBranchLengths = false;
		for (int i = 0; i < T->getNumNodes(); i++) { 
			if (T->tree[i].getBranchLength() > 0.0) {  hasBranchLengths = true; break; }
		}
		
		//set branchlengths if needed
		if (!hasBranchLengths) { 
			for (int i = 0; i < T->getNumNodes(); i++) {
				
				int lc = T->tree[i].getLChild();
				int rc = T->tree[i].getRChild();
				
				if (lc == -1) { // you are a leaf
					//if you are a leaf set you priliminary length to 1.0, this may adjust later
					T->tree[i].setBranchLength(1.0);
					dists[i] = 1.0;
				}else{ // you are an internal node
					//look at your children's length to leaf
					float ldist = dists[lc];
					float rdist = dists[rc];
					
					float greater = ldist;
					if (rdist > greater) { greater = rdist; dists[i] = ldist + 1.0;}
					else { dists[i] = rdist + 1.0; }
					
					
					//branch length = difference + 1
					T->tree[lc].setBranchLength((abs(ldist-greater) + 1.0));
					T->tree[rc].setBranchLength((abs(rdist-greater) + 1.0));
				}
			}
		}
			
		dists.clear();
		
		for (int i = 0; i < T->getNumNodes(); i++) {
			
			double sum = 0.0;
			int index = i;
			
			while(T->tree[index].getParent() != -1){
				if (T->tree[index].getBranchLength() != -1) {
					sum += abs(T->tree[index].getBranchLength()); 
				}
				index = T->tree[index].getParent();
			}
			
			dists[i] = sum;
		}
		
		return dists;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "getLengthToLeaf");	
		exit(1);
	}
}
//**********************************************************************************************************************
set<string> IndicatorCommand::getDescendantList(Tree*& T, int i, map<int, set<string> > descendants, map<int, set<int> >& nodes){
	try {
		set<string> names;
		
		set<string>::iterator it;
		
		int lc = T->tree[i].getLChild();
		int rc = T->tree[i].getRChild();
		
		if (lc == -1) { //you are a leaf your only descendant is yourself
			set<int> temp; temp.insert(i);
			nodes[i] = temp;
			
			if (designfile == "") {
				names.insert(T->tree[i].getName());
			}else {
				vector<string> myGroup; myGroup.push_back(T->tree[i].getName());
				vector<string> myReps = designMap->getNamesGroups(myGroup);
				for (int k = 0; k < myReps.size(); k++) {
					names.insert(myReps[k]);
				}
			}
			
		}else{ //your descedants are the combination of your childrens descendants
			names = descendants[lc];
			nodes[i] = nodes[lc];
			for (it = descendants[rc].begin(); it != descendants[rc].end(); it++) {
				names.insert(*it);
			}
			for (set<int>::iterator itNum = nodes[rc].begin(); itNum != nodes[rc].end(); itNum++) {
				nodes[i].insert(*itNum);
			}
			//you are your own descendant
			nodes[i].insert(i);
		}
		
		return names;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "getDescendantList");	
		exit(1);
	}
}
//**********************************************************************************************************************
int IndicatorCommand::getShared(){
	try {
		InputData input(sharedfile, "sharedfile");
		lookup = input.getSharedRAbundVectors();
		string lastLabel = lookup->getLabel();
		
		if (label == "") { label = lastLabel;  return 0; }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> labels; labels.insert(label);
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup != NULL) && (userLabels.size() != 0)) {
			if (m->control_pressed) {  return 0;  }
			
			if(labels.count(lookup->getLabel()) == 1){
				processedLabels.insert(lookup->getLabel());
				userLabels.erase(lookup->getLabel());
				break;
			}
			
			if ((m->anyLabelsToProcess(lookup->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup->getLabel();
				
                delete lookup;
				lookup = input.getSharedRAbundVectors(lastLabel);
				
				processedLabels.insert(lookup->getLabel());
				userLabels.erase(lookup->getLabel());
				
				//restore real lastlabel to save below
				lookup->setLabels(saveLabel);
				break;
			}
			
			lastLabel = lookup->getLabel();
			
			//get next line to process
			//prevent memory leak
			delete lookup;
			lookup = input.getSharedRAbundVectors();
		}
		
		
		if (m->control_pressed) {  return 0;  }
		
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
		}	
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "getShared");	
		exit(1);
	}
}
//**********************************************************************************************************************
int IndicatorCommand::getSharedFloat(){
	try {
		InputData input(relabundfile, "relabund");
		lookupFloat = input.getSharedRAbundFloatVectors();
		string lastLabel = lookupFloat->getLabel();
		
		if (label == "") { label = lastLabel;  return 0; }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> labels; labels.insert(label);
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookupFloat != NULL) && (userLabels.size() != 0)) {
			
			if (m->control_pressed) {   return 0;  }
			
			if(labels.count(lookupFloat->getLabel()) == 1){
				processedLabels.insert(lookupFloat->getLabel());
				userLabels.erase(lookupFloat->getLabel());
				break;
			}
			
			if ((m->anyLabelsToProcess(lookupFloat->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookupFloat->getLabel();
				
				delete lookupFloat;
				lookupFloat = input.getSharedRAbundFloatVectors(lastLabel);
				
				processedLabels.insert(lookupFloat->getLabel());
				userLabels.erase(lookupFloat->getLabel());
				
				//restore real lastlabel to save below
				lookupFloat->setLabels(saveLabel);
				break;
			}
			
			lastLabel = lookupFloat->getLabel();
			
			//get next line to process
			//prevent memory leak
			delete lookupFloat;
			lookupFloat = input.getSharedRAbundFloatVectors();
		}
		
		
		if (m->control_pressed) {  return 0;  }
		
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
		}	
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "getShared");	
	exit(1);
	}
}
//**********************************************************************************************************************
vector<float> IndicatorCommand::driver(vector< vector<SharedRAbundFloatVector*> >& groupings, vector< vector<string> >& groupingNames, int num, vector<float> indicatorValues, int numIters){
	try {
		vector<float> pvalues;
		pvalues.resize(indicatorValues.size(), 0);
        vector<string> notUsedGroupings;  //we dont care about the grouping for the pvalues since they are randomized, but we need to pass the function something to make it work.
		
		for(int i=0;i<numIters;i++){
			if (m->control_pressed) { break; }
			map< vector<int>, vector<int> > groupingsMap = randomizeGroupings(groupings, num);
			vector<float> randomIndicatorValues = getValues(groupings,groupingNames, notUsedGroupings, groupingsMap);
			
			for (int j = 0; j < indicatorValues.size(); j++) {
				if (randomIndicatorValues[j] >= indicatorValues[j]) { pvalues[j]++; }
			}
		}
		
		return pvalues;
		
	}catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "driver");	
		exit(1);
	}
}
//**********************************************************************************************************************
vector<float> IndicatorCommand::getPValues(vector< vector<SharedRAbundFloatVector*> >& groupings, vector< vector<string> >& groupingNames, int num, vector<float> indicatorValues){
	try {
		vector<float> pvalues;
        bool recalc = false;

		if(processors == 1){
			pvalues = driver(groupings, groupingNames, num, indicatorValues, iters);
            for (int i = 0; i < pvalues.size(); i++) { pvalues[i] /= (double)iters; }
		}else{
            //divide iters between processors
			vector<int> procIters;
            int numItersPerProcessor = iters / processors;
            
            //divide iters between processes
            for (int h = 0; h < processors; h++) {
                if(h == processors - 1){ numItersPerProcessor = iters - h * numItersPerProcessor; }
                procIters.push_back(numItersPerProcessor);
            }
            
            vector<int> processIDS;
            int process = 1;
			
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
						
			//loop through and create all the processes you want
			while (process != processors) {
				pid_t pid = fork();
				
				if (pid > 0) {
					processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
					process++;
				}else if (pid == 0){
					pvalues = driver(groupings, groupingNames, num, indicatorValues, procIters[process]);
					
					//pass pvalues to parent
					ofstream out;
					string tempFile = m->mothurGetpid(process) + ".pvalues.temp";
					m->openOutputFile(tempFile, out);
					
                    //pass values
                    for (int i = 0; i < pvalues.size(); i++) {  out << pvalues[i] << '\t'; } out << endl;
					
					out.close();
					
					exit(0);
				}else { 
                    m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(process) + "\n"); processors = process;
                    for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                    //wait to die
                    for (int i=0;i<processIDS.size();i++) {
                        int temp = processIDS[i];
                        wait(&temp);
                    }
                    m->control_pressed = false;
                    for (int i=0;i<processIDS.size();i++) {
                        m->mothurRemove((toString(processIDS[i]) + ".pvalues.temp"));
                    }
                    recalc = true;
                    break;
				}
			}
			
            if (recalc) {
                //test line, also set recalc to true.
                //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->control_pressed = false;  for (int i=0;i<processIDS.size();i++) {m->mothurRemove((toString(processIDS[i]) + ".pvalues.temp"));}processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
                
                //divide iters between processors
                processIDS.resize(0);
                process = 1;
                procIters.clear();
                int numItersPerProcessor = iters / processors;
                
                //divide iters between processes
                for (int h = 0; h < processors; h++) {
                    if(h == processors - 1){ numItersPerProcessor = iters - h * numItersPerProcessor; }
                    procIters.push_back(numItersPerProcessor);
                }
                
                //loop through and create all the processes you want
                while (process != processors) {
                    pid_t pid = fork();
                    
                    if (pid > 0) {
                        processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                        process++;
                    }else if (pid == 0){
                        pvalues = driver(groupings, groupingNames, num, indicatorValues, procIters[process]);
                        
                        //pass pvalues to parent
                        ofstream out;
                        string tempFile = m->mothurGetpid(process) + ".pvalues.temp";
                        m->openOutputFile(tempFile, out);
                        
                        //pass values
                        for (int i = 0; i < pvalues.size(); i++) {  out << pvalues[i] << '\t'; } out << endl;
                        
                        out.close();
                        
                        exit(0);
                    }else { 
                        m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
                        for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                        exit(0);
                    }
                }
            }

			//do my part
			pvalues = driver(groupings, groupingNames, num, indicatorValues, procIters[0]);
			
			//force parent to wait until all the processes are done
			for (int i=0;i<processIDS.size();i++) { 
				int temp = processIDS[i];
				wait(&temp);
			}
			
			//combine results			
			for (int i = 0; i < processIDS.size(); i++) {
				ifstream in;
				string tempFile =  toString(processIDS[i]) + ".pvalues.temp";
				m->openInputFile(tempFile, in);
				
				////// to do ///////////
				int numTemp; numTemp = 0; 
				for (int j = 0; j < pvalues.size(); j++) {
					in >> numTemp; m->gobble(in);
					pvalues[j] += numTemp;
				}
				in.close(); m->mothurRemove(tempFile);
			}
			for (int i = 0; i < pvalues.size(); i++) { pvalues[i] /= (double)iters; }
#else
                       
            //fill in functions
            vector<indicatorData*> pDataArray;
            DWORD   dwThreadIdArray[processors-1];
            HANDLE  hThreadArray[processors-1];
            
            //Create processor worker threads.
            for( int i=1; i<processors; i++ ){
                
                //make copy of lookup so we don't get access violations
                vector< vector<SharedRAbundFloatVector*> > newGroupings;

                for (int k = 0; k < groupings.size(); k++) {
                    vector<SharedRAbundFloatVector*> newLookup;
                    for (int l = 0; l < groupings[k].size(); l++) {
                        SharedRAbundFloatVector* temp = new SharedRAbundFloatVector();
                        temp->setLabel(groupings[k][l]->getLabel());
                        temp->setGroup(groupings[k][l]->getGroup());
                        newLookup.push_back(temp);
                    }
                    newGroupings.push_back(newLookup);
                }
                
                //for each bin
                for (int l = 0; l < groupings.size(); l++) {
                    for (int k = 0; k < groupings[l][0]->getNumBins(); k++) {
                        if (m->control_pressed) { for (int j = 0; j < newGroupings.size(); j++) { for (int u = 0; u < newGroupings[j].size(); u++) { delete newGroupings[j][u];  } } return pvalues; }
                        
                        for (int j = 0; j < groupings[l].size(); j++) { newGroupings[l][j]->push_back(groupings[l][j]->get(k)); }
                    }
                }
        
                vector<float> copyIValues = indicatorValues;

                indicatorData* temp = new indicatorData(m, procIters[i], newGroupings, num, copyIValues);
                pDataArray.push_back(temp);
                processIDS.push_back(i);
                
                hThreadArray[i-1] = CreateThread(NULL, 0, MyIndicatorThreadFunction, pDataArray[i-1], 0, &dwThreadIdArray[i-1]);
            }
            
            //do my part
			pvalues = driver(groupings, groupingNames, num, indicatorValues, procIters[0]);
           
            //Wait until all threads have terminated.
            WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
            
            //Close all thread handles and free memory allocations.
            for(int i=0; i < pDataArray.size(); i++){
                for (int j = 0; j < pDataArray[i]->pvalues.size(); j++) { pvalues[j] += pDataArray[i]->pvalues[j];  }
                
                for (int l = 0; l < pDataArray[i]->groupings.size(); l++) {
                    for (int j = 0; j < pDataArray[i]->groupings[l].size(); j++) { delete pDataArray[i]->groupings[l][j]; }
                }
                
                CloseHandle(hThreadArray[i]);
                delete pDataArray[i];
            }
            
            for (int i = 0; i < pvalues.size(); i++) { pvalues[i] /= (double)iters; }
#endif
		}

		
		return pvalues;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "getPValues");	
		exit(1);
	}
}

//**********************************************************************************************************************
//same as above, just data type difference
vector<float> IndicatorCommand::driver(vector< vector<SharedRAbundVector*> >& groupings, vector< vector<string> >& groupingNames, int num, vector<float> indicatorValues, int numIters){
	try {
		vector<float> pvalues;
		pvalues.resize(indicatorValues.size(), 0);
        vector<string> notUsedGroupings;  //we dont care about the grouping for the pvalues since they are randomized, but we need to pass the function something to make it work.
		
		for(int i=0;i<numIters;i++){
			if (m->control_pressed) { break; }
			map< vector<int>, vector<int> > groupingsMap = randomizeGroupings(groupings, num);
			vector<float> randomIndicatorValues = getValues(groupings, groupingNames, notUsedGroupings, groupingsMap);
			
			for (int j = 0; j < indicatorValues.size(); j++) {
				if (randomIndicatorValues[j] >= indicatorValues[j]) { pvalues[j]++; }
			}
		}
		
		return pvalues;
		
	}catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "driver");	
		exit(1);
	}
}
//**********************************************************************************************************************
//same as above, just data type difference
vector<float> IndicatorCommand::getPValues(vector< vector<SharedRAbundVector*> >& groupings, vector< vector<string> >& groupingNames, int num, vector<float> indicatorValues){
	try {
		vector<float> pvalues;
        bool recalc = false;
        
		if(processors == 1){
			pvalues = driver(groupings, groupingNames, num, indicatorValues, iters);
            for (int i = 0; i < pvalues.size(); i++) { pvalues[i] /= (double)iters; }
		}else{
            //divide iters between processors
			vector<int> procIters;
            int numItersPerProcessor = iters / processors;
            
            //divide iters between processes
            for (int h = 0; h < processors; h++) {
                if(h == processors - 1){ numItersPerProcessor = iters - h * numItersPerProcessor; }
                procIters.push_back(numItersPerProcessor);
            }
            
            vector<int> processIDS;
            int process = 1;
			
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
            
			//loop through and create all the processes you want
			while (process != processors) {
				pid_t pid = fork();
				
				if (pid > 0) {
					processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
					process++;
				}else if (pid == 0){
					pvalues = driver(groupings, groupingNames, num, indicatorValues, procIters[process]);
					
					//pass pvalues to parent
					ofstream out;
					string tempFile = m->mothurGetpid(process) + ".pvalues.temp";
					m->openOutputFile(tempFile, out);
					
					//pass values
					for (int i = 0; i < pvalues.size(); i++) {
						out << pvalues[i] << '\t';
					}
					out << endl;
					
					out.close();
					
					exit(0);
				}else {
                    m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(process) + "\n"); processors = process;
                    for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                    //wait to die
                    for (int i=0;i<processIDS.size();i++) {
                        int temp = processIDS[i];
                        wait(&temp);
                    }
                    m->control_pressed = false;
                    for (int i=0;i<processIDS.size();i++) {
                        m->mothurRemove((toString(processIDS[i]) + ".pvalues.temp"));
                    }
                    recalc = true;
                    break;
				}
			}
			
            if (recalc) {
                //test line, also set recalc to true.
                //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->control_pressed = false;  for (int i=0;i<processIDS.size();i++) {m->mothurRemove((toString(processIDS[i]) + ".pvalues.temp"));}processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
                
                //divide iters between processors
                processIDS.resize(0);
                process = 1;
                procIters.clear();
                int numItersPerProcessor = iters / processors;
                
                //divide iters between processes
                for (int h = 0; h < processors; h++) {
                    if(h == processors - 1){ numItersPerProcessor = iters - h * numItersPerProcessor; }
                    procIters.push_back(numItersPerProcessor);
                }
                
                //loop through and create all the processes you want
                while (process != processors) {
                    pid_t pid = fork();
                    
                    if (pid > 0) {
                        processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                        process++;
                    }else if (pid == 0){
                        pvalues = driver(groupings, groupingNames, num, indicatorValues, procIters[process]);
                        
                        //pass pvalues to parent
                        ofstream out;
                        string tempFile = m->mothurGetpid(process) + ".pvalues.temp";
                        m->openOutputFile(tempFile, out);
                        
                        //pass values
                        for (int i = 0; i < pvalues.size(); i++) {  out << pvalues[i] << '\t'; } out << endl;
                        
                        out.close();
                        
                        exit(0);
                    }else {
                        m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine();
                        for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                        exit(0);
                    }
                }
            }

			//do my part
			pvalues = driver(groupings, groupingNames, num, indicatorValues, procIters[0]);
			
			//force parent to wait until all the processes are done
			for (int i=0;i<processIDS.size();i++) {
				int temp = processIDS[i];
				wait(&temp);
			}
			
			//combine results
			for (int i = 0; i < processIDS.size(); i++) {
				ifstream in;
				string tempFile =  toString(processIDS[i]) + ".pvalues.temp";
				m->openInputFile(tempFile, in);
				
				////// to do ///////////
				int numTemp; numTemp = 0;
				for (int j = 0; j < pvalues.size(); j++) {
					in >> numTemp; m->gobble(in);
					pvalues[j] += numTemp;
				}
				in.close(); m->mothurRemove(tempFile);
			}
			for (int i = 0; i < pvalues.size(); i++) { pvalues[i] /= (double)iters; }
#else
            
            //fill in functions
            vector<indicatorData*> pDataArray;
            DWORD   dwThreadIdArray[processors-1];
            HANDLE  hThreadArray[processors-1];
            
            //Create processor worker threads.
            for( int i=1; i<processors; i++ ){
                
                //make copy of lookup so we don't get access violations
                vector< vector<SharedRAbundFloatVector*> > newGroupings;
                
                for (int k = 0; k < groupings.size(); k++) {
                    vector<SharedRAbundFloatVector*> newLookup;
                    for (int l = 0; l < groupings[k].size(); l++) {
                        SharedRAbundFloatVector* temp = new SharedRAbundFloatVector();
                        temp->setLabel(groupings[k][l]->getLabel());
                        temp->setGroup(groupings[k][l]->getGroup());
                        newLookup.push_back(temp);
                    }
                    newGroupings.push_back(newLookup);
                }
                
                //for each bin
                for (int l = 0; l < groupings.size(); l++) {
                    for (int k = 0; k < groupings[l][0]->getNumBins(); k++) {
                        if (m->control_pressed) { for (int j = 0; j < newGroupings.size(); j++) { for (int u = 0; u < newGroupings[j].size(); u++) { delete newGroupings[j][u];  } } return pvalues; }
                        
                        for (int j = 0; j < groupings[l].size(); j++) { newGroupings[l][j]->push_back((float)(groupings[l][j]->get(k))); }
                    }
                }
                
                vector<float> copyIValues = indicatorValues;
                
                indicatorData* temp = new indicatorData(m, procIters[i], newGroupings, num, copyIValues);
                pDataArray.push_back(temp);
                processIDS.push_back(i);
                
                hThreadArray[i-1] = CreateThread(NULL, 0, MyIndicatorThreadFunction, pDataArray[i-1], 0, &dwThreadIdArray[i-1]);
            }
            
            //do my part
			pvalues = driver(groupings, num, indicatorValues, procIters[0]);
            
            //Wait until all threads have terminated.
            WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
            
            //Close all thread handles and free memory allocations.
            for(int i=0; i < pDataArray.size(); i++){
                for (int j = 0; j < pDataArray[i]->pvalues.size(); j++) { pvalues[j] += pDataArray[i]->pvalues[j];  }
                
                for (int l = 0; l < pDataArray[i]->groupings.size(); l++) {
                    for (int j = 0; j < pDataArray[i]->groupings[l].size(); j++) { delete pDataArray[i]->groupings[l][j]; }
                }
                
                CloseHandle(hThreadArray[i]);
                delete pDataArray[i];
            }
            
            for (int i = 0; i < pvalues.size(); i++) { pvalues[i] /= (double)iters; }
#endif
		}
        
		
		return pvalues;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "getPValues");
		exit(1);
	}
}
//**********************************************************************************************************************
//swap groups between groupings, in essence randomizing the second column of the design file
map< vector<int>, vector<int> > IndicatorCommand::randomizeGroupings(vector< vector<SharedRAbundVector*> >& groupings, int numLookupGroups){
	try {
		
		map< vector<int>, vector<int> > randomGroupings; 
		
		for (int i = 0; i < numLookupGroups; i++) {
			if (m->control_pressed) {break;}
			
			//get random groups to swap to switch with
			//generate random int between 0 and groupings.size()-1
			int z = m->getRandomIndex(groupings.size()-1);
			int x = m->getRandomIndex(groupings.size()-1);
			int a = m->getRandomIndex(groupings[z].size()-1);
			int b = m->getRandomIndex(groupings[x].size()-1);
			//cout << i << '\t' << z << '\t' << x << '\t' << a << '\t' << b << endl;	
			//if ((z < 0) || (z > 1) || x<0 || x>1 || a <0 || a>groupings[z].size()-1 || b<0 || b>groupings[x].size()-1) { cout << "probelm" << i << '\t' << z << '\t' << x << '\t' << a << '\t' << b << endl;	}
			
			vector<int> from;
			vector<int> to;
			
			from.push_back(z); from.push_back(a);
			to.push_back(x); to.push_back(b);
			
			randomGroupings[from] = to;
		}
		//cout << "done" << endl;	
		return randomGroupings;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "randomizeGroupings");	
		exit(1);
	}
}	
//**********************************************************************************************************************
//swap groups between groupings, in essence randomizing the second column of the design file
map< vector<int>, vector<int> > IndicatorCommand::randomizeGroupings(vector< vector<SharedRAbundFloatVector*> >& groupings, int numLookupGroups){
	try {
		
		map< vector<int>, vector<int> > randomGroupings; 
		
		for (int i = 0; i < numLookupGroups; i++) {
			
			//get random groups to swap to switch with
			//generate random int between 0 and groupings.size()-1
			int z = m->getRandomIndex(groupings.size()-1);
			int x = m->getRandomIndex(groupings.size()-1);
			int a = m->getRandomIndex(groupings[z].size()-1);
			int b = m->getRandomIndex(groupings[x].size()-1);
			//cout << i << '\t' << z << '\t' << x << '\t' << a << '\t' << b << endl;		
			
			vector<int> from;
			vector<int> to;
			
			from.push_back(z); from.push_back(a);
			to.push_back(x); to.push_back(b);
			
			randomGroupings[from] = to;
		}
		
		return randomGroupings;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "randomizeGroupings");	
		exit(1);
	}
}			
								
/*****************************************************************/


