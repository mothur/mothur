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
		CommandParameter pdesign("design", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pdesign);
		CommandParameter pshared("shared", "InputTypes", "", "", "SharedRel", "SharedRel", "none",false,false); parameters.push_back(pshared);	
		CommandParameter prelabund("relabund", "InputTypes", "", "", "SharedRel", "SharedRel", "none",false,false); parameters.push_back(prelabund);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter ptree("tree", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(ptree);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
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
		helpString += "The indicator command reads a shared or relabund file and a tree file, and outputs a .indicator.tre and .indicator.summary file. \n";
		helpString += "The new tree contains labels at each internal node.  The label is the node number so you can relate the tree to the summary file.\n";
		helpString += "The summary file lists the indicator value for each OTU for each node.\n";
		helpString += "The indicator command parameters are tree, groups, shared, relabund, design and label. The tree parameter is required as well as either shared or relabund.\n";
		helpString += "The design parameter allows you to provide a design file to relate the tree to the shared or relabund file.\n";		
		helpString += "The groups parameter allows you to specify which of the groups in your shared or relabund you would like analyzed, or if you provide a design file the groups in your design file.  The groups may be entered separated by dashes.\n";
		helpString += "The label parameter indicates at what distance your tree relates to the shared or relabund.\n";
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
			m->Groups.clear();
			m->namesOfGroups.clear();
			m->Treenames.clear();
			m->names.clear();
			
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
			if (treefile == "not open") { abort = true; }
			else if (treefile == "not found") { 				//if there is a current design file, use it
				treefile = m->getTreeFile(); 
				if (treefile != "") { m->mothurOut("Using " + treefile + " as input file for the tree parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current tree file and the tree parameter is required."); m->mothurOutEndLine(); abort = true; }								
			}	
			
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { abort = true; }
			else if (sharedfile == "not found") { sharedfile = ""; }
			else { inputFileName = sharedfile; }
			
			relabundfile = validParameter.validFile(parameters, "relabund", true);
			if (relabundfile == "not open") { abort = true; }
			else if (relabundfile == "not found") { relabundfile = ""; }
			else { inputFileName = relabundfile; }
			
			designfile = validParameter.validFile(parameters, "design", true);
			if (designfile == "not open") { abort = true; }
			else if (designfile == "not found") { designfile = ""; }
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = "";  Groups.push_back("all"); }
			else { m->splitAtDash(groups, Groups);	}			
			m->Groups = Groups;
			
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; m->mothurOut("You did not provide a label, I will use the first label in your inputfile."); m->mothurOutEndLine(); label=""; }	
			
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
			
			if ((relabundfile != "") && (sharedfile != "")) { m->mothurOut("You may not use both a shared and relabund file."); m->mothurOutEndLine(); abort = true;  }
			
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
	
		//read designfile if given and set up globaldatas groups for read of sharedfiles
		if (designfile != "") {
			designMap = new GroupMap(designfile);
			designMap->readDesignMap();
			
			//fill Groups - checks for "all" and for any typo groups
			SharedUtil* util = new SharedUtil();
			util->setGroups(Groups, designMap->namesOfGroups);
			delete util;
			
			//loop through the Groups and fill Globaldata's Groups with the design file info
			m->Groups = designMap->getNamesSeqs(Groups);
		}
	
		/***************************************************/
		// use smart distancing to get right sharedRabund  //
		/***************************************************/
		if (sharedfile != "") {  
			getShared(); 
			if (m->control_pressed) { if (designfile != "") { delete designMap; } for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } return 0; }
			if (lookup[0] == NULL) { m->mothurOut("[ERROR] reading shared file."); m->mothurOutEndLine(); return 0; }
		}else { 
			getSharedFloat(); 
			if (m->control_pressed) {  if (designfile != "") { delete designMap; } for (int i = 0; i < lookupFloat.size(); i++) {  delete lookupFloat[i];  } return 0; }
			if (lookupFloat[0] == NULL) { m->mothurOut("[ERROR] reading relabund file."); m->mothurOutEndLine(); return 0; }
		}
		
		//reset Globaldatas groups if needed
		if (designfile != "") { m->Groups = Groups; }
			
		/***************************************************/
		//    reading tree info							   //
		/***************************************************/
		string groupfile = ""; 
		m->setTreeFile(treefile);
		Tree* tree = new Tree(treefile); delete tree;  //extracts names from tree to make faked out groupmap
		treeMap = new TreeMap();
		bool mismatch = false;
			
		for (int i = 0; i < m->Treenames.size(); i++) { 
			//sanity check - is this a group that is not in the sharedfile?
			if (designfile == "") {
				if (!(m->inUsersGroups(m->Treenames[i], m->namesOfGroups))) {
					m->mothurOut("[ERROR]: " + m->Treenames[i] + " is not a group in your shared or relabund file."); m->mothurOutEndLine();
					mismatch = true;
				}
				treeMap->addSeq(m->Treenames[i], "Group1"); 
			}else{
				vector<string> myGroups; myGroups.push_back(m->Treenames[i]);
				vector<string> myNames = designMap->getNamesSeqs(myGroups);
				
				for(int k = 0; k < myNames.size(); k++) {
					if (!(m->inUsersGroups(myNames[k], m->namesOfGroups))) {
						m->mothurOut("[ERROR]: " + myNames[k] + " is not a group in your shared or relabund file."); m->mothurOutEndLine();
						mismatch = true;
					}
				}
				treeMap->addSeq(m->Treenames[i], "Group1");
			}
		}
		
		if ((designfile != "") && (m->Treenames.size() != Groups.size())) { m->mothurOut("[ERROR]: You design file does not match your tree, aborting."); m->mothurOutEndLine(); mismatch = true; }
				
		if (mismatch) { //cleanup and exit
			if (designfile != "") { delete designMap; }
			if (sharedfile != "") {  for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } }
			else { for (int i = 0; i < lookupFloat.size(); i++) {  delete lookupFloat[i];  } }
			delete treeMap;
			return 0;
		}
	 
		read = new ReadNewickTree(treefile);
		int readOk = read->read(treeMap); 
		
		if (readOk != 0) { m->mothurOut("Read Terminated."); m->mothurOutEndLine(); delete treeMap; delete read; return 0; }
		
		vector<Tree*> T = read->getTrees();
		
		delete read;
		
		if (m->control_pressed) { 
			if (designfile != "") { delete designMap; }
			if (sharedfile != "") {  for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } }
			else { for (int i = 0; i < lookupFloat.size(); i++) {  delete lookupFloat[i];  } }
			for (int i = 0; i < T.size(); i++) {  delete T[i];  }  delete treeMap; return 0; 
		}
			
		T[0]->assembleTree();
				
		/***************************************************/
		//    create ouptut tree - respecting pickedGroups //
		/***************************************************/
		Tree* outputTree = new Tree(m->Groups.size(), treeMap); 
		
		outputTree->getSubTree(T[0], m->Groups);
		outputTree->assembleTree();
			
		//no longer need original tree, we have output tree to use and label
		for (int i = 0; i < T.size(); i++) {  delete T[i];  } 
		
				
		if (m->control_pressed) { 
			if (designfile != "") { delete designMap; }
			if (sharedfile != "") {  for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } }
			else { for (int i = 0; i < lookupFloat.size(); i++) {  delete lookupFloat[i];  } }
			delete outputTree; delete treeMap;  return 0; 
		}
		
		/***************************************************/
		//		get indicator species values			   //
		/***************************************************/
		GetIndicatorSpecies(outputTree);
		
		if (designfile != "") { delete designMap; }
		if (sharedfile != "") {  for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } }
		else { for (int i = 0; i < lookupFloat.size(); i++) {  delete lookupFloat[i];  } }
		delete outputTree; delete treeMap;
		
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());	} return 0; }
		
		//set tree file as new current treefile
		string current = "";
		itTypes = outputTypes.find("tree");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setTreeFile(current); }
		}
		
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
//traverse tree finding indicator species values for each otu at each node
//label node with otu number that has highest indicator value
//report all otu values to file
int IndicatorCommand::GetIndicatorSpecies(Tree*& T){
	try {
		
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(inputFileName);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(inputFileName)) + "indicator.summary";
		outputNames.push_back(outputFileName); outputTypes["summary"].push_back(outputFileName);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		int numBins = 0;
		if (sharedfile != "") { numBins = lookup[0]->getNumBins(); }
		else { numBins = lookupFloat[0]->getNumBins(); }
		
		//print headings
		out << "TreeNode\t";
		for (int i = 0; i < numBins; i++) { out << "OTU-" << (i+1) << '\t'; }
		out << endl;
		
		string treeOutputDir = outputDir;
		if (outputDir == "") {  treeOutputDir += m->hasPath(treefile);  }
		string outputTreeFileName = treeOutputDir + m->getRootName(m->getSimpleName(treefile)) + "indicator.tre";
		
		
		//create a map from tree node index to names of descendants, save time later to know which sharedRabund you need
		map<int, set<string> > nodeToDescendants;
		map<int, set<int> > descendantNodes;
		for (int i = 0; i < T->getNumNodes(); i++) {
			if (m->control_pressed) { return 0; }
			
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
			
			if (sharedfile != "") {
				vector< vector<SharedRAbundVector*> > groupings;
				
				//get nodes that will be a valid grouping
				//you are valid if you are not one of my descendants
				//AND your distToRoot is >= mine
				//AND you were not added as part of a larger grouping. Largest nodes are added first.
				
				set<string> groupsAlreadyAdded;
				//create a grouping with my grouping
				vector<SharedRAbundVector*> subset;
				int count = 0;
				int doneCount = nodeToDescendants[i].size();
				for (int k = 0; k < lookup.size(); k++) {
					//is this descendant of i
					if ((nodeToDescendants[i].count(lookup[k]->getGroup()) != 0)) { 
						subset.push_back(lookup[k]);
						groupsAlreadyAdded.insert(lookup[k]->getGroup());
						count++;
					}
					if (count == doneCount) { break; } //quit once you get the rabunds for this grouping
				}
				if (subset.size() != 0) { groupings.push_back(subset); }
				
				
				for (int j = (T->getNumNodes()-1); j >= 0; j--) {
	
					
					if ((descendantNodes[i].count(j) == 0) && (distToRoot[j] >= distToRoot[i])) { 
						vector<SharedRAbundVector*> subset;
						int count = 0;
						int doneCount = nodeToDescendants[j].size();
						for (int k = 0; k < lookup.size(); k++) {
							//is this descendant of j, and we didn't already add this as part of a larger grouping
							if ((nodeToDescendants[j].count(lookup[k]->getGroup()) != 0) && (groupsAlreadyAdded.count(lookup[k]->getGroup()) == 0)) { 
								subset.push_back(lookup[k]);
								groupsAlreadyAdded.insert(lookup[k]->getGroup());
								count++;
							}
							if (count == doneCount) { break; } //quit once you get the rabunds for this grouping
						}
						
						//if subset.size == 0 then the node was added as part of a larger grouping
						if (subset.size() != 0) { groupings.push_back(subset); }
					}
				}
				
				if (groupsAlreadyAdded.size() != lookup.size()) {  m->mothurOut("[ERROR]: could not make proper groupings."); m->mothurOutEndLine(); }
								
				indicatorValues = getValues(groupings);
				
			}else {
				vector< vector<SharedRAbundFloatVector*> > groupings;
				
				//get nodes that will be a valid grouping
				//you are valid if you are not one of my descendants
				//AND your distToRoot is >= mine
				//AND you were not added as part of a larger grouping. Largest nodes are added first.
				
				set<string> groupsAlreadyAdded;
				//create a grouping with my grouping
				vector<SharedRAbundFloatVector*> subset;
				int count = 0;
				int doneCount = nodeToDescendants[i].size();
				for (int k = 0; k < lookupFloat.size(); k++) {
					//is this descendant of i
					if ((nodeToDescendants[i].count(lookupFloat[k]->getGroup()) != 0)) { 
						subset.push_back(lookupFloat[k]);
						groupsAlreadyAdded.insert(lookupFloat[k]->getGroup());
						count++;
					}
					if (count == doneCount) { break; } //quit once you get the rabunds for this grouping
				}
				if (subset.size() != 0) { groupings.push_back(subset); }
				
				for (int j = (T->getNumNodes()-1); j >= 0; j--) {
					if ((descendantNodes[i].count(j) == 0) && (distToRoot[j] >= distToRoot[i])) {
						vector<SharedRAbundFloatVector*> subset;
						int count = 0;
						int doneCount = nodeToDescendants[j].size();
						for (int k = 0; k < lookupFloat.size(); k++) {
							//is this descendant of j, and we didn't already add this as part of a larger grouping
							if ((nodeToDescendants[j].count(lookupFloat[k]->getGroup()) != 0) && (groupsAlreadyAdded.count(lookupFloat[k]->getGroup()) == 0)) { 
								subset.push_back(lookupFloat[k]);
								groupsAlreadyAdded.insert(lookupFloat[k]->getGroup());
								count++;
							}
							if (count == doneCount) { break; } //quit once you get the rabunds for this grouping
						}
						
						//if subset.size == 0 then the node was added as part of a larger grouping
						if (subset.size() != 0) { groupings.push_back(subset); }
					}
				}
				
				if (groupsAlreadyAdded.size() != lookupFloat.size()) { m->mothurOut("[ERROR]: could not make proper groupings."); m->mothurOutEndLine(); }
				
				indicatorValues = getValues(groupings);
			}
			
			if (m->control_pressed) { out.close(); return 0; }
			
			
			/******************************************************/
			//output indicator values to table form + label tree  //
			/*****************************************************/
			out << (i+1) << '\t';
			for (int j = 0; j < indicatorValues.size(); j++) {
				
				if (m->control_pressed) { out.close(); return 0; }
				
				out << indicatorValues[j] << '\t';
			}
			out << endl;
			
			T->tree[i].setLabel((i+1));
			
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
vector<float> IndicatorCommand::getValues(vector< vector<SharedRAbundFloatVector*> >& groupings){
	try {
		vector<float> values;
		
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
					totalAbund += groupings[j][k]->getAbundance(i); 
					if (groupings[j][k]->getAbundance(i) != 0) { numNotZero++; }
				}
				
				float Aij = (totalAbund / (float) groupings[j].size());
				terms.push_back(Aij);
				
				//percentage of sites represented
				Bij.push_back(numNotZero / (float) groupings[j].size());
				
				AijDenominator += Aij;
			}
			
			float maxIndVal = 0.0;
			for (int j = 0; j < terms.size(); j++) { 
				float thisAij = (terms[j] / AijDenominator); 
				float thisValue = thisAij * Bij[j] * 100.0;
				
				//save largest
				if (thisValue > maxIndVal) { maxIndVal = thisValue; }
			}
			
			values.push_back(maxIndVal);
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
vector<float> IndicatorCommand::getValues(vector< vector<SharedRAbundVector*> >& groupings){
	try {
		vector<float> values;
		
	/*for (int j = 0; j < groupings.size(); j++) {		
		cout << "grouping " << j << endl;
		for (int k = 0; k < groupings[j].size(); k++) { 
			cout << groupings[j][k]->getGroup() << endl;
		}
	}*/
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
					totalAbund += groupings[j][k]->getAbundance(i); 
					if (groupings[j][k]->getAbundance(i) != 0.0) { numNotZero++; }
					
				}
					
				float Aij = (totalAbund / (float) groupings[j].size());
				terms.push_back(Aij);
				
				//percentage of sites represented
				Bij.push_back(numNotZero / (float) groupings[j].size());
				
				AijDenominator += Aij;
			}
			
			float maxIndVal = 0.0;
			for (int j = 0; j < terms.size(); j++) { 
				float thisAij = (terms[j] / AijDenominator); 
				float thisValue = thisAij * Bij[j] * 100.0;
					
				//save largest
				if (thisValue > maxIndVal) { maxIndVal = thisValue; }
			}
			
			values.push_back(maxIndVal);
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
				vector<string> myReps = designMap->getNamesSeqs(myGroup);
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
		InputData* input = new InputData(sharedfile, "sharedfile");
		lookup = input->getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
		
		if (label == "") { label = lastLabel; delete input; return 0; }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> labels; labels.insert(label);
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && (userLabels.size() != 0)) {
			if (m->control_pressed) {  delete input; return 0;  }
			
			if(labels.count(lookup[0]->getLabel()) == 1){
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				break;
			}
			
			if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup[0]->getLabel();
				
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
				lookup = input->getSharedRAbundVectors(lastLabel);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				
				//restore real lastlabel to save below
				lookup[0]->setLabel(saveLabel);
				break;
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
		}	
		
		delete input;
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
		InputData* input = new InputData(relabundfile, "relabund");
		lookupFloat = input->getSharedRAbundFloatVectors();
		string lastLabel = lookupFloat[0]->getLabel();
		
		if (label == "") { label = lastLabel; delete input; return 0; }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> labels; labels.insert(label);
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookupFloat[0] != NULL) && (userLabels.size() != 0)) {
			
			if (m->control_pressed) {  delete input; return 0;  }
			
			if(labels.count(lookupFloat[0]->getLabel()) == 1){
				processedLabels.insert(lookupFloat[0]->getLabel());
				userLabels.erase(lookupFloat[0]->getLabel());
				break;
			}
			
			if ((m->anyLabelsToProcess(lookupFloat[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookupFloat[0]->getLabel();
				
				for (int i = 0; i < lookupFloat.size(); i++) {  delete lookupFloat[i];  } 
				lookupFloat = input->getSharedRAbundFloatVectors(lastLabel);
				
				processedLabels.insert(lookupFloat[0]->getLabel());
				userLabels.erase(lookupFloat[0]->getLabel());
				
				//restore real lastlabel to save below
				lookupFloat[0]->setLabel(saveLabel);
				break;
			}
			
			lastLabel = lookupFloat[0]->getLabel();			
			
			//get next line to process
			//prevent memory leak
			for (int i = 0; i < lookupFloat.size(); i++) {  delete lookupFloat[i];  } 
			lookupFloat = input->getSharedRAbundFloatVectors();
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
			for (int i = 0; i < lookupFloat.size(); i++) {  if (lookupFloat[i] != NULL) {	delete lookupFloat[i];	} } 
			lookupFloat = input->getSharedRAbundFloatVectors(lastLabel);
		}	
		
		delete input;
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "getShared");	
	exit(1);
	}
}		
/*****************************************************************/


