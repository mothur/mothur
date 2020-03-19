/*
 *  indicatorcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 11/12/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "indicatorcommand.h"



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
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
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
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
			vector<string> tempOutNames;
			outputTypes["tree"] = tempOutNames;
			outputTypes["summary"] = tempOutNames;
			
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";	}

			//check for required parameters
			treefile = validParameter.validFile(parameters, "tree");
			if (treefile == "not open") { treefile = ""; abort = true; }
			else if (treefile == "not found") { treefile = ""; 	}		
			else { current->setTreeFile(treefile); }	
			
			sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { abort = true; }
			else if (sharedfile == "not found") { sharedfile = ""; }
			else { inputFileName = sharedfile; current->setSharedFile(sharedfile); }
			
			relabundfile = validParameter.validFile(parameters, "relabund");
			if (relabundfile == "not open") { abort = true; }
			else if (relabundfile == "not found") { relabundfile = ""; }
			else { inputFileName = relabundfile; current->setRelAbundFile(relabundfile); }
			
			designfile = validParameter.validFile(parameters, "design");
			if (designfile == "not open") { designfile = ""; abort = true; }
			else if (designfile == "not found") { designfile = ""; }
			else { current->setDesignFile(designfile); }
			
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = "";   }
			else { util.splitAtDash(groups, Groups); if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }	}
			
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; m->mothurOut("You did not provide a label, I will use the first label in your inputfile."); m->mothurOutEndLine(); label=""; }	
			
			string temp = validParameter.valid(parameters, "iters");		if (temp == "not found") { temp = "1000"; }
			util.mothurConvert(temp, iters); 
			
			temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
            
			if ((relabundfile == "") && (sharedfile == "")) { 
				//is there are current file available for either of these?
				//give priority to shared, then relabund
				sharedfile = current->getSharedFile(); 
				if (sharedfile != "") {  inputFileName = sharedfile; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 
					relabundfile = current->getRelAbundFile(); 
					if (relabundfile != "") {  inputFileName = relabundfile; m->mothurOut("Using " + relabundfile + " as input file for the relabund parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("No valid current files. You must provide a shared or relabund."); m->mothurOutEndLine(); 
						abort = true;
					}
				}
			}
			
			
			if ((designfile == "") && (treefile == "")) { 
				treefile = current->getTreeFile(); 
				if (treefile != "") {  m->mothurOut("Using " + treefile + " as input file for the tree parameter."); m->mothurOutEndLine(); }
				else { 
					designfile = current->getDesignFile(); 
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
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		cout.setf(ios::fixed, ios::floatfield); cout.setf(ios::showpoint);
		long start = time(NULL);

		//read designfile if given and set up groups for read of sharedfiles
        vector<string> allGroups;
		if (designfile != "") {
			designMap = new DesignMap(designfile);
			
            if (Groups.size() == 0) { Groups = designMap->getCategory(); }
            allGroups = designMap->getCategory();
						
			namesSeqs = designMap->getNamesGroups(Groups);
        }else {
            InputData* input = NULL;
            if (sharedfile != "") { input = new InputData(sharedfile, "sharedfile", Groups);  }
            else { input = new InputData(sharedfile, "relabundfile", Groups); }
            
            SharedRAbundVectors* lookup = input->getSharedRAbundVectors();
            Groups = lookup->getNamesGroups();
            namesSeqs = Groups;
            delete lookup; delete input;
        }
	
		if (treefile != "") {
			string groupfile = ""; 
			current->setTreeFile(treefile);
            vector<string> Treenames = util.parseTreeFile(treefile);
			CountTable ct;
			bool mismatch = false;
            set<string> nameMap; map<string, string> groupMap; set<string> gps;
            
            for (int i = 0; i < Treenames.size(); i++) {
                nameMap.insert(Treenames[i]);
                //sanity check - is this a group that is not in the sharedfile?
                if (i == 0) { gps.insert("Group1"); }
				if (designfile == "") {
					if (!(util.inUsersGroups(Treenames[i], namesSeqs))) { m->mothurOut("[ERROR]: " + Treenames[i] + " is not a group in your shared or relabund file.\n"); mismatch = true; }
					groupMap[Treenames[i]] = "Group1";
				}else{
					vector<string> myGroups; myGroups.push_back(Treenames[i]);
					vector<string> myNames = designMap->getNamesGroups(myGroups);
					
					for(int k = 0; k < myNames.size(); k++) {
						if (!(util.inUsersGroups(myNames[k], allGroups))) { m->mothurOut("[ERROR]: " + myNames[k] + " is not a group in your shared or relabund file.\n"); mismatch = true; }
					}
					groupMap[Treenames[i]] = designMap->get(Treenames[i]);
				}
            }
            ct.createTable(nameMap, groupMap, gps);
			
			if ((designfile != "") && (Treenames.size() != namesSeqs.size())) {  m->mothurOut("[ERROR]: You design file does not match your tree, aborting.\n"); mismatch = true; }
					
			if (mismatch) { //cleanup and exit
				if (designfile != "") { delete designMap; }  return 0;
			}
		 
			ReadTree* read = new ReadNewickTree(treefile, Treenames);
			int readOk = read->read(&ct);   if (readOk != 0) { m->mothurOut("Read Terminated.\n");  delete read; return 0; }
			vector<Tree*> T = read->getTrees();  delete read;
			
			if (m->getControl_pressed()) {  if (designfile != "") { delete designMap; } for (int i = 0; i < T.size(); i++) {  delete T[i];  }   return 0; }
            
			T[0]->assembleTree();
					
			Tree* outputTree = new Tree(namesSeqs.size(), &ct, Treenames);  //create ouptut tree - respecting pickedGroups
			outputTree->getSubTree(T[0], namesSeqs);
			outputTree->assembleTree();
				
			//no longer need original tree, we have output tree to use and label
			for (int i = 0; i < T.size(); i++) {  delete T[i];  } 
			
			if (m->getControl_pressed()) {  if (designfile != "") { delete designMap; } delete outputTree;   return 0; }
			
			GetIndicatorSpecies(outputTree);  //get indicator species values
			delete outputTree;
            
		}else { GetIndicatorSpecies(); }  //run with design file only
		
		if (designfile != "") { delete designMap; }
        
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);	} return 0; }
		
		//set tree file as new current treefile
		if (treefile != "") {
			string currentName = "";
			itTypes = outputTypes.find("tree");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setTreeFile(currentName); }
			}
		}
		
		m->mothurOut("\n\nIt took " + toString(time(NULL) - start) + " secs to find the indicator species.\n");
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
				
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "execute");	
		exit(1);
	}
}
//**********************************************************************************************************************
//indicatorValues = getValues(groupings, groupingNames, indicatorGroups, randomGroupingsMap, m);
vector<float> getValues(vector< vector<SharedRAbundFloatVector*> >& groupings, vector< vector<string> >& groupingNames, vector<string>& indicatorGroupings, map<sharedIndexes, sharedIndexes> groupingsMap, MothurOut* m){
    try {
        vector<float> values;
        map<sharedIndexes, sharedIndexes>::iterator it;
        
        indicatorGroupings.clear();
        
        //create grouping strings
        vector<string> groupingsGroups;
        for (int j = 0; j < groupings.size(); j++) {
            string tempGrouping = "";
            for (int k = 0; k < groupings[j].size()-1; k++) {  tempGrouping += groupingNames[j][k] + "-"; }
            
            tempGrouping += groupingNames[j][groupingNames[j].size()-1];
            groupingsGroups.push_back(tempGrouping);
        }
        
        Utils util;
        //for each otu
        for (int i = 0; i < groupings[0][0]->getNumBins(); i++) {
            
            if (m->getControl_pressed()) { return values; }
            
            vector<float> terms;
            float AijDenominator = 0.0;
            vector<float> Bij;
            
            //get overall abundance of each grouping
            for (int j = 0; j < groupings.size(); j++) {
                
                float totalAbund = 0;
                int numNotZero = 0;
                for (int k = 0; k < groupings[j].size(); k++) {
                    sharedIndexes temp(j,k);
                    it = groupingsMap.find(temp);
                    
                    if (it == groupingsMap.end()) { //this one didnt get moved, or initial calc
                        totalAbund += groupings[j][k]->get(i);
                        if (!util.isEqual(groupings[j][k]->get(i), 0)) { numNotZero++; }
                    }else {
                        float thisAbund = groupings[(it->second).groupIndex][(it->second).otuIndex]->get(i);
                        totalAbund += thisAbund;
                        if (!util.isEqual(thisAbund, 0)) { numNotZero++; }
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
//indicatorValues = getValues(groupings, groupingNames, indicatorGroups, randomGroupingsMap, m);
vector<float> getValues(vector< vector<SharedRAbundVector*> >& groupings, vector< vector<string> >& groupingNames, vector<string>& indicatorGroupings, map<sharedIndexes, sharedIndexes> groupingsMap, MothurOut* m){
    try {
        vector<float> values;
        map<sharedIndexes, sharedIndexes>::iterator it;
        
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
                    sharedIndexes temp(j,k);
                    it = groupingsMap.find(temp);
                    
                    if (it == groupingsMap.end()) { //this one didnt get moved
                        totalAbund += groupings[j][k]->get(i);
                        if (groupings[j][k]->get(i) != 0) { numNotZero++; }
                    }else {
                        totalAbund += groupings[(it->second).groupIndex][(it->second).otuIndex]->get(i);
                        if (groupings[(it->second).groupIndex][(it->second).otuIndex]->get(i) != 0) { numNotZero++; }
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
//divide shared or relabund file by groupings in the design file
//report all otu values to file
int IndicatorCommand::GetIndicatorSpecies(){
	try {
        
        SharedRAbundVectors* lookup = NULL; SharedRAbundFloatVectors* lookupFloat = NULL;
        if (sharedfile != "") {
            lookup = getShared();
            if (m->getControl_pressed()) { if (designfile != "") { delete designMap; } delete lookup;   return 0; }
            if (lookup == NULL) { m->mothurOut("[ERROR] reading shared file."); m->mothurOutEndLine(); return 0; }
        }else {
            lookupFloat = getSharedFloat();
            if (m->getControl_pressed()) {  if (designfile != "") { delete designMap; } delete lookupFloat; return 0; }
            if (lookupFloat == NULL) { m->mothurOut("[ERROR] reading relabund file."); m->mothurOutEndLine(); return 0; }
        }

        
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(inputFileName);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(inputFileName));
		string outputFileName = getOutputFileName("summary", variables);
		outputNames.push_back(outputFileName); outputTypes["summary"].push_back(outputFileName);
		
		ofstream out;
		util.openOutputFile(outputFileName, out);
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		m->mothurOut("\nSpecies\tIndicator_Groups\tIndicatorValue\tpValue\n");
		
		int numBins = 0;
        vector<string> currentLabels;
        if (sharedfile != "") { numBins = lookup->getNumBins(); currentLabels = lookup->getOTUNames(); }
		else { numBins = lookupFloat->getNumBins(); currentLabels = lookupFloat->getOTUNames(); }
		
        if (m->getControl_pressed()) { out.close(); if (sharedfile != "") {  delete lookup; }
            else { delete lookupFloat; } return 0; }
			
		/*****************************************************/
		//create vectors containing rabund info				 //
		/*****************************************************/
			
		vector<float> indicatorValues; //size of numBins
		vector<float> pValues;
        vector<string> indicatorGroups;
		map<sharedIndexes, sharedIndexes> randomGroupingsMap; //maps location in groupings to location in groupings, ie, [0][0] -> [1][2]. This is so we don't have to actually move the sharedRabundVectors.
			
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
				
			if (groupsAlreadyAdded.size() != data.size()) {  m->mothurOut("[ERROR]: could not make proper groupings.\n"); }
				
			indicatorValues = getValues(groupings, groupingNames, indicatorGroups, randomGroupingsMap, m);
			
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
			
			indicatorValues = getValues(groupings, groupingNames, indicatorGroups, randomGroupingsMap, m);
			
			pValues = getPValues(groupings, groupingNames, lookupFloat->size(), indicatorValues);
		}
			
		if (m->getControl_pressed()) { out.close(); return 0; }
			
			
		/******************************************************/
		//output indicator values to table form               //
		/*****************************************************/
		out << "OTU\tIndicator_Groups\tIndicator_Value\tpValue" << endl;
		for (int j = 0; j < indicatorValues.size(); j++) {
				
			if (m->getControl_pressed()) { out.close(); return 0; }
			
			out << currentLabels[j] << '\t' << indicatorGroups[j] << '\t' << indicatorValues[j] << '\t';
			
			if (pValues[j] > (1/(float)iters)) { out << pValues[j] << endl; } 
			else { out << "<" << (1/(float)iters) << endl; }
			
			if (pValues[j] <= 0.05) {
				string pValueString = "<" + toString((1/(float)iters)); 
				if (pValues[j] > (1/(float)iters)) { pValueString = toString(pValues[j]); }
				m->mothurOut(currentLabels[j] + "\t" + indicatorGroups[j] + "\t" + toString(indicatorValues[j]) + "\t" + pValueString + "\n");
			}
		}
		
		out.close();
        if (sharedfile != "") {  delete lookup; }
        else { delete lookupFloat; }
		
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
        SharedRAbundVectors* lookup = NULL; SharedRAbundFloatVectors* lookupFloat = NULL;
        if (sharedfile != "") {
            lookup = getShared();
            if (m->getControl_pressed()) { if (designfile != "") { delete designMap; } delete lookup;   return 0; }
            if (lookup == NULL) { m->mothurOut("[ERROR] reading shared file."); m->mothurOutEndLine(); return 0; }
        }else {
            lookupFloat = getSharedFloat();
            if (m->getControl_pressed()) {  if (designfile != "") { delete designMap; } delete lookupFloat; return 0; }
            if (lookupFloat == NULL) { m->mothurOut("[ERROR] reading relabund file."); m->mothurOutEndLine(); return 0; }
        }
        
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(inputFileName);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(inputFileName));
		string outputFileName = getOutputFileName("summary",variables);
		outputNames.push_back(outputFileName); outputTypes["summary"].push_back(outputFileName);
		
		ofstream out;
		util.openOutputFile(outputFileName, out);
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		int numBins = 0;
        vector<string> currentLabels;
        if (sharedfile != "") { numBins = lookup->getNumBins(); currentLabels = lookup->getOTUNames(); }
        else { numBins = lookupFloat->getNumBins(); currentLabels = lookupFloat->getOTUNames(); }
		
		//print headings
		out << "TreeNode\t";
		for (int i = 0; i < numBins; i++) { out << currentLabels[i] << "_IndGroups" << '\t' << currentLabels[i] << "_IndValue" << '\t' << "pValue" << '\t'; } out << endl;
		
		m->mothurOut("\nNode\tSpecies\tIndicator_Groups\tIndicatorValue\tpValue\n");
		
		string treeOutputDir = outputDir;
		if (outputDir == "") {  treeOutputDir += util.hasPath(treefile);  }
        variables["[filename]"] = treeOutputDir + util.getRootName(util.getSimpleName(treefile));
		string outputTreeFileName = getOutputFileName("tree", variables);
		
		//create a map from tree node index to names of descendants, save time later to know which sharedRabund you need
		map<int, set<string> > nodeToDescendants;
		map<int, set<int> > descendantNodes;
		for (int i = 0; i < T->getNumNodes(); i++) { if (m->getControl_pressed()) { return 0; }
			nodeToDescendants[i] = getDescendantList(T, i, nodeToDescendants, descendantNodes);
		}
		
		//you need the distances to leaf to decide grouping below
		//this will also set branch lengths if the tree does not include them
		map<int, float> distToRoot = getDistToRoot(T);
			
		//for each node
		for (int i = T->getNumLeaves(); i < T->getNumNodes(); i++) {
			
			if (m->getControl_pressed()) { out.close(); return 0; }
			
			/*****************************************************/
			//create vectors containing rabund info				 //
			/*****************************************************/
			
			vector<float> indicatorValues; //size of numBins
			vector<float> pValues;
            vector<string> indicatorGroups;
			vector< map<sharedIndexes, sharedIndexes> > randomGroupingsMap; //maps location in groupings to location in groupings, ie, [0][0] -> [1][2]. This is so we don't have to actually move the sharedRabundVectors.
			
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
                vector<SharedRAbundVector*> data = lookup->getSharedRAbundVectors();
                vector<string> dataGroupNames = lookup->getNamesGroups();

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
				
				if (groupsAlreadyAdded.size() != data.size()) {  m->mothurOut("[ERROR]: could not make proper groupings.\n"); }
							
                map<sharedIndexes, sharedIndexes> placeHolder; //don't need randomization for initial calc
				indicatorValues = getValues(groupings, groupingNames, indicatorGroups, placeHolder, m);
				
				pValues = getPValues(groupings, groupingNames, lookup->getNumGroups(), indicatorValues);
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
				
                map<sharedIndexes, sharedIndexes> placeHolder; //don't need randomization for initial calc
				indicatorValues = getValues(groupings, groupingNames, indicatorGroups, placeHolder, m);
				
				pValues = getPValues(groupings, groupingNames, lookupFloat->getNumGroups(),  indicatorValues);
			}
			
			if (m->getControl_pressed()) { out.close(); return 0; }
			
			
			/******************************************************/
			//output indicator values to table form + label tree  //
			/*****************************************************/
			out << (i+1);
			for (int j = 0; j < indicatorValues.size(); j++) {
				
                if (m->getControl_pressed()) { out.close(); if (sharedfile != "") {  delete lookup; }
                    else { delete lookupFloat; } return 0; }
				
				if (pValues[j] < (1/(float)iters)) {
					out  << '\t' << indicatorGroups[j] << '\t' << indicatorValues[j] << '\t' << '<' << (1/(float)iters);
				}else {
					out  << '\t' << indicatorGroups[j] << '\t' << indicatorValues[j] << '\t' << pValues[j];
				}
				
				if (pValues[j] <= 0.05) {
					string pValueString = "\t<" + toString((1/(float)iters));
					if (pValues[j] > (1/(float)iters)) { pValueString = toString('\t' + pValues[j]); }
					m->mothurOut(toString(i+1) + "\t" + currentLabels[j] + "\t" + indicatorGroups[j] + "\t" + toString(indicatorValues[j]) + "\t" + pValueString + "\n");
				}
			}
			out << endl;
			
			T->tree[i].setLabel(toString(i+1));
		}
		out.close();
	
		ofstream outTree;
		util.openOutputFile(outputTreeFileName, outTree);
		outputNames.push_back(outputTreeFileName); outputTypes["tree"].push_back(outputTreeFileName);
		
		T->print(outTree, "both");
		outTree.close();
        if (sharedfile != "") {  delete lookup; }
        else { delete lookupFloat; }
	
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "GetIndicatorSpecies");	
		exit(1);
	}
}
/**************************************************************************************************/

struct indicatorFloatData {
    vector< vector<SharedRAbundFloatVector*> > groupings;
    vector< vector<string> > groupingNames;
    vector< map<sharedIndexes, sharedIndexes> > randomGroupings;
   	MothurOut* m;
    int iters, numGroups;
    vector<float> indicatorValues, pvalues;
    
    indicatorFloatData(){}
    indicatorFloatData(int it, vector< map<sharedIndexes, sharedIndexes> > ran, vector< vector<SharedRAbundFloatVector*> > ng, vector< vector<string> > gn, int n, vector<float> iv) {
        m = MothurOut::getInstance();
        iters = it;
        groupings = ng;
        groupingNames = gn;
        randomGroupings = ran;
        indicatorValues = iv;
        pvalues.resize(indicatorValues.size(), 0);
        numGroups = n;
    }
};

struct indicatorData {
    vector< vector<SharedRAbundVector*> > groupings;
    vector< vector<string> > groupingNames;
    vector< map<sharedIndexes, sharedIndexes> > randomGroupings;
   	MothurOut* m;
    int iters, numGroups;
    vector<float> indicatorValues, pvalues;
    
    indicatorData(){}
    indicatorData(int it, vector< map<sharedIndexes, sharedIndexes> > ran, vector< vector<SharedRAbundVector*> > ng, vector< vector<string> > gn, int n, vector<float> iv) {
        m = MothurOut::getInstance();
        iters = it;
        groupings = ng;
        groupingNames = gn;
        randomGroupings = ran;
        indicatorValues = iv;
        pvalues.resize(indicatorValues.size(), 0);
        numGroups = n;
    }
};
//**********************************************************************************************************************
void driverValues(indicatorData* params){
    try {
        vector<string> notUsedGroupings;  //we dont care about the grouping for the pvalues since they are randomized, but we need to pass the function something
        
        for(int i=0;i<params->iters;i++){
            if (params->m->getControl_pressed()) { break; }
            map<sharedIndexes, sharedIndexes> groupingsMap = params->randomGroupings[i];
            
            vector<float> randomIndicatorValues = getValues(params->groupings, params->groupingNames, notUsedGroupings, groupingsMap, params->m);
            
            for (int j = 0; j < params->indicatorValues.size(); j++) {  if (randomIndicatorValues[j] >= params->indicatorValues[j]) { params->pvalues[j]++; } }
        }
        
    }catch(exception& e) {
        params->m->errorOut(e, "IndicatorCommand", "driver");
        exit(1);
    }
}
//**********************************************************************************************************************
void driverValuesFloat(indicatorFloatData* params){
	try {
        vector<string> notUsedGroupings;  //we dont care about the grouping for the pvalues since they are randomized, but we need to pass the function something
		
		for(int i=0;i<params->iters;i++){
			if (params->m->getControl_pressed()) { break; }
			map<sharedIndexes, sharedIndexes> groupingsMap = params->randomGroupings[i];
			vector<float> randomIndicatorValues = getValues(params->groupings, params->groupingNames, notUsedGroupings, groupingsMap, params->m);
			
			for (int j = 0; j < params->indicatorValues.size(); j++) {
				if (randomIndicatorValues[j] >= params->indicatorValues[j]) { params->pvalues[j]++; }
            }
		}
		
	}catch(exception& e) {
		params->m->errorOut(e, "IndicatorCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<float> IndicatorCommand::getPValues(vector< vector<SharedRAbundFloatVector*> >& groupings, vector< vector<string> >& groupingNames, int num, vector<float> indicatorValues){
    try {
        vector<float> pvalues;
        vector<int> groupingsSize;
        for (int i = 0; i < groupings.size(); i++) {  groupingsSize.push_back(groupings[i].size());  }
        vector< map<sharedIndexes, sharedIndexes> > randomize = randomizeGroupings(groupingsSize, groupings.size());
        
        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<indicatorFloatData*> data;
        
        //divide iters between processors
        vector<int> procIters;
        int numItersPerProcessor = iters / processors;
        
        //divide iters between processes
        for (int h = 0; h < processors; h++) {
            if(h == processors - 1){ numItersPerProcessor = iters - h * numItersPerProcessor; }
            procIters.push_back(numItersPerProcessor);
        }
        
        //Lauch worker threads
        int start = 0;
        for (int i = 0; i < processors-1; i++) {
            
            //make copy of lookup so we don't get access violations
            vector< vector<SharedRAbundFloatVector*> > newGroupings;
            for (int k = 0; k < groupings.size(); k++) {
                vector<SharedRAbundFloatVector*> newLookup;
                for (int l = 0; l < groupings[k].size(); l++) {
                    SharedRAbundFloatVector* temp = new SharedRAbundFloatVector(*groupings[k][l]);
                    newLookup.push_back(temp);
                }
                newGroupings.push_back(newLookup);
            }
            
            vector< map<sharedIndexes, sharedIndexes> > thisProcessorsRandom;
            thisProcessorsRandom.insert(thisProcessorsRandom.begin(), randomize.begin()+start, randomize.begin()+start+procIters[i]); start += procIters[i];
            indicatorFloatData* dataBundle = new indicatorFloatData(procIters[i], thisProcessorsRandom, newGroupings, groupingNames, num, indicatorValues);
            
            data.push_back(dataBundle);
            
            workerThreads.push_back(new std::thread(driverValuesFloat, dataBundle));
        }
        
        //make copy of lookup so we don't get access violations
        vector< vector<SharedRAbundFloatVector*> > newGroupings;
        for (int k = 0; k < groupings.size(); k++) {
            vector<SharedRAbundFloatVector*> newLookup;
            for (int l = 0; l < groupings[k].size(); l++) {
                SharedRAbundFloatVector* temp = new SharedRAbundFloatVector(*groupings[k][l]);
                newLookup.push_back(temp);
            }
            newGroupings.push_back(newLookup);
        }
        
        vector< map<sharedIndexes, sharedIndexes> > thisProcessorsRandom;
        thisProcessorsRandom.insert(thisProcessorsRandom.begin(), randomize.begin()+start, randomize.begin()+start+procIters[processors-1]);
        indicatorFloatData* dataBundle = new indicatorFloatData(procIters[processors-1], thisProcessorsRandom, newGroupings, groupingNames, num, indicatorValues);
        
        driverValuesFloat(dataBundle);
        pvalues = dataBundle->pvalues;
        for (int l = 0; l < newGroupings.size(); l++) { for (int j = 0; j < newGroupings[l].size(); j++) { delete newGroupings[l][j]; } }
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            
            for (int j = 0; j < data[i]->pvalues.size(); j++) { pvalues[j] += data[i]->pvalues[j];  }
            for (int l = 0; l < data[i]->groupings.size(); l++) { for (int j = 0; j < data[i]->groupings[l].size(); j++) { delete data[i]->groupings[l][j]; } }
            
            delete data[i];
            delete workerThreads[i];
        }
        delete dataBundle;
        
        
        for (int i = 0; i < pvalues.size(); i++) { pvalues[i] /= (double)iters; }
        
        return pvalues;
    }
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "getPValues");
		exit(1);
	}
}
//**********************************************************************************************************************
//same as above, just data type difference
vector<float> IndicatorCommand::getPValues(vector< vector<SharedRAbundVector*> >& groupings, vector< vector<string> >& groupingNames, int num, vector<float> indicatorValues){
	try {
		vector<float> pvalues;
        vector<int> groupingsSize;
        for (int i = 0; i < groupings.size(); i++) {  groupingsSize.push_back(groupings[i].size());  }
        vector< map<sharedIndexes, sharedIndexes> > randomize = randomizeGroupings(groupingsSize, groupings.size());
        
        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<indicatorData*> data;
        
        //divide iters between processors
        vector<int> procIters;
        int numItersPerProcessor = iters / processors;
        
        //divide iters between processes
        for (int h = 0; h < processors; h++) {
            if(h == processors - 1){ numItersPerProcessor = iters - h * numItersPerProcessor; }
            procIters.push_back(numItersPerProcessor);
        }
        
        //Lauch worker threads
        int start = 0;
        for (int i = 0; i < processors-1; i++) {
            
            //make copy of lookup so we don't get access violations
            vector< vector<SharedRAbundVector*> > newGroupings;
            for (int k = 0; k < groupings.size(); k++) {
                vector<SharedRAbundVector*> newLookup;
                for (int l = 0; l < groupings[k].size(); l++) {
                    SharedRAbundVector* temp = new SharedRAbundVector(*groupings[k][l]);
                    newLookup.push_back(temp);
                }
                newGroupings.push_back(newLookup);
            }

            vector< map<sharedIndexes, sharedIndexes> > thisProcessorsRandom;
            thisProcessorsRandom.insert(thisProcessorsRandom.begin(), randomize.begin()+start, randomize.begin()+start+procIters[i]); start += procIters[i];
            indicatorData* dataBundle = new indicatorData(procIters[i], thisProcessorsRandom, newGroupings, groupingNames, num, indicatorValues);
           
            data.push_back(dataBundle);
            
            workerThreads.push_back(new std::thread(driverValues, dataBundle));
        }
        
        //make copy of lookup so we don't get access violations
        vector< vector<SharedRAbundVector*> > newGroupings;
        for (int k = 0; k < groupings.size(); k++) {
            vector<SharedRAbundVector*> newLookup;
            for (int l = 0; l < groupings[k].size(); l++) {
                SharedRAbundVector* temp = new SharedRAbundVector(*groupings[k][l]);
                newLookup.push_back(temp);
            }
            newGroupings.push_back(newLookup);
        }

        vector< map<sharedIndexes, sharedIndexes> > thisProcessorsRandom;
        thisProcessorsRandom.insert(thisProcessorsRandom.begin(), randomize.begin()+start, randomize.begin()+start+procIters[processors-1]);
        indicatorData* dataBundle = new indicatorData(procIters[processors-1], thisProcessorsRandom, newGroupings, groupingNames, num, indicatorValues);

        driverValues(dataBundle);
        pvalues = dataBundle->pvalues;
        for (int l = 0; l < newGroupings.size(); l++) { for (int j = 0; j < newGroupings[l].size(); j++) { delete newGroupings[l][j]; } }
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            
            for (int j = 0; j < data[i]->pvalues.size(); j++) { pvalues[j] += data[i]->pvalues[j];  }
            for (int l = 0; l < data[i]->groupings.size(); l++) { for (int j = 0; j < data[i]->groupings[l].size(); j++) { delete data[i]->groupings[l][j]; } }
            
            delete data[i];
            delete workerThreads[i];
        }
        delete dataBundle;
 
        
        for (int i = 0; i < pvalues.size(); i++) { pvalues[i] /= (double)iters; }
        
        return pvalues;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "getPValues");
		exit(1);
	}
}

//**********************************************************************************************************************
//swap groups between groupings, in essence randomizing the second column of the design file
vector< map<sharedIndexes, sharedIndexes> > IndicatorCommand::randomizeGroupings(vector<int> sizesOfEachGrouping, int numLookupGroups){
	try {
		
		vector< map<sharedIndexes, sharedIndexes> > randomGroupings;
        
        for (int j = 0; j < iters; j++) {
            
            if (m->getControl_pressed()) {break;}
            map<sharedIndexes, sharedIndexes> thisRandomization;
            
            for (int i = 0; i < numLookupGroups; i++) {
                
                //select random group and random OTU to swap
                int z = util.getRandomIndex(numLookupGroups-1);
                int x = util.getRandomIndex(numLookupGroups-1);
                int a = util.getRandomIndex(sizesOfEachGrouping[z]-1);
                int b = util.getRandomIndex(sizesOfEachGrouping[x]-1);
                
                sharedIndexes from(z, a);
                sharedIndexes to(x, b);
                thisRandomization[from] = to;
            }
            randomGroupings.push_back(thisRandomization);
        }
		return randomGroupings;
	}
	catch(exception& e) {
		m->errorOut(e, "IndicatorCommand", "randomizeGroupings");	
		exit(1);
	}
}
//**********************************************************************************************************************
SharedRAbundVectors* IndicatorCommand::getShared(){
    try {
        InputData input(sharedfile, "sharedfile", namesSeqs);
        SharedRAbundVectors* lookup = input.getSharedRAbundVectors();
        Groups = lookup->getNamesGroups();
        string lastLabel = lookup->getLabel();
        
        if (label == "") { label = lastLabel;  return lookup; }
        
        //if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
        set<string> labels; labels.insert(label);
        set<string> processedLabels;
        set<string> userLabels = labels;
        
        //as long as you are not at the end of the file or done wih the lines you want
        while((lookup != NULL) && (userLabels.size() != 0)) {
            if (m->getControl_pressed()) {  return lookup;  }
            
            if(labels.count(lookup->getLabel()) == 1){
                processedLabels.insert(lookup->getLabel()); userLabels.erase(lookup->getLabel());
                break;
            }
            
            if ((util.anyLabelsToProcess(lookup->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
                string saveLabel = lookup->getLabel();
                
                delete lookup;
                lookup = input.getSharedRAbundVectors(lastLabel);
                
                processedLabels.insert(lookup->getLabel()); userLabels.erase(lookup->getLabel());
                
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
        
        
        if (m->getControl_pressed()) {  return 0;  }
        
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
            lookup = input.getSharedRAbundVectors(lastLabel);
        }
        
        return lookup;
    }
    catch(exception& e) {
        m->errorOut(e, "IndicatorCommand", "getShared");
        exit(1);
    }
}
//**********************************************************************************************************************
SharedRAbundFloatVectors* IndicatorCommand::getSharedFloat(){
    try {
        InputData input(relabundfile, "relabund", namesSeqs);
        SharedRAbundFloatVectors* lookupFloat = input.getSharedRAbundFloatVectors();
        Groups = lookupFloat->getNamesGroups();
        string lastLabel = lookupFloat->getLabel();
        
        if (label == "") { label = lastLabel;  return lookupFloat; }
        
        //if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
        set<string> labels; labels.insert(label);
        set<string> processedLabels;
        set<string> userLabels = labels;
        
        //as long as you are not at the end of the file or done wih the lines you want
        while((lookupFloat != NULL) && (userLabels.size() != 0)) {
            
            if (m->getControl_pressed()) {   return lookupFloat;  }
            
            if(labels.count(lookupFloat->getLabel()) == 1){
                processedLabels.insert(lookupFloat->getLabel());
                userLabels.erase(lookupFloat->getLabel());
                break;
            }
            
            if ((util.anyLabelsToProcess(lookupFloat->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
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
        
        
        if (m->getControl_pressed()) {  return 0;  }
        
        //output error messages about any remaining user labels
        bool needToRun = false;
        for (set<string>::iterator it = userLabels.begin(); it != userLabels.end(); it++) {
            m->mothurOut("Your file does not include the label " + *it); 
            if (processedLabels.count(lastLabel) != 1)  { m->mothurOut(". I will use " + lastLabel + ".\n"); needToRun = true;  }
            else                                        { m->mothurOut(". Please refer to " + lastLabel + ".\n");               }
        }
        
        //run last label if you need to
        if (needToRun )  {
            delete lookupFloat;
            lookupFloat = input.getSharedRAbundFloatVectors(lastLabel);
        }	
        
        return lookupFloat;
    }
    catch(exception& e) {
        m->errorOut(e, "IndicatorCommand", "getShared");	
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
                if (!util.isEqual(T->tree[index].getBranchLength(), -1)) {
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
            
            if (designfile == "") { names.insert(T->tree[i].getName());  }
            else {
                //string myRep = designMap->get(T->tree[i].getName());
                names.insert(T->tree[i].getName());
            }
            
        }else{ //your descedants are the combination of your childrens descendants
            names = descendants[lc];
            nodes[i] = nodes[lc];
            for (it = descendants[rc].begin(); it != descendants[rc].end(); it++) {  names.insert(*it); }
            
            for (set<int>::iterator itNum = nodes[rc].begin(); itNum != nodes[rc].end(); itNum++) { nodes[i].insert(*itNum); }
            
            nodes[i].insert(i); //you are your own descendant
        }
        
        return names;
    }
    catch(exception& e) {
        m->errorOut(e, "IndicatorCommand", "getDescendantList");	
        exit(1);
    }
}
/*****************************************************************/


