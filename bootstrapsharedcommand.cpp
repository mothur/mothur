/*
 *  bootstrapsharedcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/16/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "bootstrapsharedcommand.h"
#include "sharedjabund.h"
#include "sharedsorabund.h"
#include "sharedjclass.h"
#include "sharedsorclass.h"
#include "sharedjest.h"
#include "sharedsorest.h"
#include "sharedthetayc.h"
#include "sharedthetan.h"
#include "sharedmorisitahorn.h"
#include "sharedbraycurtis.h"


//**********************************************************************************************************************
vector<string> BootSharedCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pshared);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter piters("iters", "Number", "", "1000", "", "", "",false,false); parameters.push_back(piters);
		CommandParameter pcalc("calc", "Multiple", "jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan-morisitahorn-braycurtis", "jclass-thetayc", "", "", "",true,false); parameters.push_back(pcalc);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);

		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "BootSharedCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string BootSharedCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The bootstrap.shared command parameters are shared, groups, calc, iters and label. shared is required.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like included used.\n";
		helpString += "The group names are separated by dashes. The label parameter allows you to select what distance levels you would like trees created for, and is also separated by dashes.\n";
		helpString += "The bootstrap.shared command should be in the following format: bootstrap.shared(groups=yourGroups, calc=yourCalcs, label=yourLabels, iters=yourIters).\n";
		helpString += "Example bootstrap.shared(groups=A-B-C, calc=jabund-sorabund, iters=100).\n";
		helpString += "The default value for groups is all the groups in your groupfile.\n";
		helpString += "The default value for calc is jclass-thetayc. The default for iters is 1000.\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "BootSharedCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
BootSharedCommand::BootSharedCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["tree"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "BootSharedCommand", "BootSharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
BootSharedCommand::BootSharedCommand(string option) {
	try {
		abort = false; calledHelp = false;   
		
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
			outputTypes["tree"] = tempOutNames;
			
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
			if (sharedfile == "not found") { m->mothurOut("shared is a required parameter.");  m->mothurOutEndLine(); sharedfile = ""; abort = true; }
			else if (sharedfile == "not open") { sharedfile = ""; abort = true; }
		
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(sharedfile); //if user entered a file with a path then preserve it	
			}
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
							
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				m->splitAtDash(groups, Groups);
				m->Groups = Groups;
			}
				
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "jclass-thetayc";  }
			else { 
				 if (calc == "default")  {  calc = "jclass-thetayc";  }
			}
			m->splitAtDash(calc, Estimators);

			string temp;
			temp = validParameter.validFile(parameters, "iters", false);  if (temp == "not found") { temp = "1000"; }
			convert(temp, iters); 
				
			if (abort == false) {
			
				//used in tree constructor 
				m->runParse = false;
			
				validCalculator = new ValidCalculators();
				
				
				//NOTE: if you add a calc to this if statement you must add it to the setParameters function
				//or it will not be visible in the gui
				int i;
				for (i=0; i<Estimators.size(); i++) {
					if (validCalculator->isValidCalculator("boot", Estimators[i]) == true) { 
						if (Estimators[i] == "jabund") { 	
							treeCalculators.push_back(new JAbund());
						}else if (Estimators[i] == "sorabund") { 
							treeCalculators.push_back(new SorAbund());
						}else if (Estimators[i] == "jclass") { 
							treeCalculators.push_back(new Jclass());
						}else if (Estimators[i] == "sorclass") { 
							treeCalculators.push_back(new SorClass());
						}else if (Estimators[i] == "jest") { 
							treeCalculators.push_back(new Jest());
						}else if (Estimators[i] == "sorest") { 
							treeCalculators.push_back(new SorEst());
						}else if (Estimators[i] == "thetayc") { 
							treeCalculators.push_back(new ThetaYC());
						}else if (Estimators[i] == "thetan") { 
							treeCalculators.push_back(new ThetaN());
						}else if (Estimators[i] == "morisitahorn") { 
							treeCalculators.push_back(new MorHorn());
						}else if (Estimators[i] == "braycurtis") { 
							treeCalculators.push_back(new BrayCurtis());
						}
					}
				}
				
				delete validCalculator;
				
				ofstream* tempo;
				for (int i=0; i < treeCalculators.size(); i++) {
					tempo = new ofstream;
					out.push_back(tempo);
				}
				
				//make a vector of tree* for each calculator
				trees.resize(treeCalculators.size());
			}
		}

	}
	catch(exception& e) {
		m->errorOut(e, "BootSharedCommand", "BootSharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
BootSharedCommand::~BootSharedCommand(){}
//**********************************************************************************************************************

int BootSharedCommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		m->mothurOut("bootstrap.shared command is no longer available."); m->mothurOutEndLine();
	/*
		//read first line
		input = new InputData(sharedfile, "sharedfile");
		order = input->getSharedOrderVector();
		string lastLabel = order->getLabel();
		
		//if the users entered no valid calculators don't execute command
		if (treeCalculators.size() == 0) { return 0; }

		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
				
		//set users groups
		util = new SharedUtil();
		util->setGroups(m->Groups, m->namesOfGroups, "treegroup");
		
		numGroups = m->Groups.size();
		
		//clear globaldatas old tree names if any
		globaldata->Treenames.clear();
		
		//fills globaldatas tree names
		globaldata->Treenames = m->Groups;
		
		//create treemap class from groupmap for tree class to use
		tmap = new TreeMap();
		tmap->makeSim(globaldata->gGroupmap);
		globaldata->gTreemap = tmap;
			
		while((order != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } globaldata->Groups.clear(); delete input;delete util; return 0;	} 
	
			if(allLines == 1 || labels.count(order->getLabel()) == 1){			
				
				m->mothurOut(order->getLabel()); m->mothurOutEndLine();
				int error = process(order);
				if (error == 1) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } globaldata->Groups.clear(); delete input;delete util; return 0;	} 
				
				processedLabels.insert(order->getLabel());
				userLabels.erase(order->getLabel());
			}
			
			//you have a label the user want that is smaller than this label and the last label has not already been processed
			if ((m->anyLabelsToProcess(order->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = order->getLabel();
				
				delete order;
				order = input->getSharedOrderVector(lastLabel);													
				m->mothurOut(order->getLabel()); m->mothurOutEndLine();
				int error = process(order);
				if (error == 1) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } globaldata->Groups.clear(); delete input;delete util; return 0;	} 

				processedLabels.insert(order->getLabel());
				userLabels.erase(order->getLabel());
				
				//restore real lastlabel to save below
				order->setLabel(saveLabel);
			}
			
			
			lastLabel = order->getLabel();			

			//get next line to process
			delete order;
			order = input->getSharedOrderVector();
		}
		
		
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } globaldata->Groups.clear(); delete input; delete util;  return 0;	} 

		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			m->mothurOut("Your file does not include the label " + *it); 
			if (processedLabels.count(lastLabel) != 1) {
				m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
				needToRun = true;
			}else {
				m->mothurOut(". Please refer to " + lastLabel + ".");  m->mothurOutEndLine();
			}
		}
		
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } globaldata->Groups.clear(); delete input; delete util; return 0;	} 

		//run last line if you need to
		if (needToRun == true)  {
				if (order != NULL) {	delete order;	}
				order = input->getSharedOrderVector(lastLabel);													
				m->mothurOut(order->getLabel()); m->mothurOutEndLine();
				int error = process(order);
				if (error == 1) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } globaldata->Groups.clear(); delete input; delete util; return 0;	} 
				
				delete order;

		}
		
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } globaldata->Groups.clear();delete input; delete util; return 0;	} 

		//reset groups parameter
		globaldata->Groups.clear();  
		
		//set first tree file as new current treefile
		string currentTree = "";
		itTypes = outputTypes.find("tree");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentTree = (itTypes->second)[0]; m->setTreeFile(currentTree); }
		}
		
		delete input;
		delete util;
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
*/

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "BootSharedCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int BootSharedCommand::createTree(ostream* out, Tree* t){
	try {
		/*
		//do merges and create tree structure by setting parents and children
		//there are numGroups - 1 merges to do
		for (int i = 0; i < (numGroups - 1); i++) {
		
			if (m->control_pressed) {  return 1; }
		
			float largest = -1000.0;
			int row, column;
			//find largest value in sims matrix by searching lower triangle
			for (int j = 1; j < simMatrix.size(); j++) {
				for (int k = 0; k < j; k++) {
					if (simMatrix[j][k] > largest) {  largest = simMatrix[j][k]; row = j; column = k;  }
				}
			}

			//set non-leaf node info and update leaves to know their parents
			//non-leaf
			t->tree[numGroups + i].setChildren(index[row], index[column]);
		
			//parents
			t->tree[index[row]].setParent(numGroups + i);
			t->tree[index[column]].setParent(numGroups + i);
			
			//blength = distance / 2;
			float blength = ((1.0 - largest) / 2);
			
			//branchlengths
			t->tree[index[row]].setBranchLength(blength - t->tree[index[row]].getLengthToLeaves());
			t->tree[index[column]].setBranchLength(blength - t->tree[index[column]].getLengthToLeaves());
		
			//set your length to leaves to your childs length plus branchlength
			t->tree[numGroups + i].setLengthToLeaves(t->tree[index[row]].getLengthToLeaves() + t->tree[index[row]].getBranchLength());
			
		
			//update index 
			index[row] = numGroups+i;
			index[column] = numGroups+i;
			
			//zero out highest value that caused the merge.
			simMatrix[row][column] = -1000.0;
			simMatrix[column][row] = -1000.0;
		
			//merge values in simsMatrix
			for (int n = 0; n < simMatrix.size(); n++)	{
				//row becomes merge of 2 groups
				simMatrix[row][n] = (simMatrix[row][n] + simMatrix[column][n]) / 2;
				simMatrix[n][row] = simMatrix[row][n];
				//delete column
				simMatrix[column][n] = -1000.0;
				simMatrix[n][column] = -1000.0;
			}
		}
		
		//adjust tree to make sure root to tip length is .5
		int root = t->findRoot();
		t->tree[root].setBranchLength((0.5 - t->tree[root].getLengthToLeaves()));

		//assemble tree
		t->assembleTree();
		
		if (m->control_pressed) { return 1; }
	
		//print newick file
		t->print(*out);*/
		
		return 0;
	
	}
	catch(exception& e) {
		m->errorOut(e, "BootSharedCommand", "createTree");
		exit(1);
	}
}
/***********************************************************/
void BootSharedCommand::printSims() {
	try {
		m->mothurOut("simsMatrix"); m->mothurOutEndLine(); 
		for (int k = 0; k < simMatrix.size(); k++)	{
			for (int n = 0; n < simMatrix.size(); n++)	{
				m->mothurOut(toString(simMatrix[k][n]));  m->mothurOut("\t"); 
			}
			m->mothurOutEndLine(); 
		}

	}
	catch(exception& e) {
		m->errorOut(e, "BootSharedCommand", "printSims");
		exit(1);
	}
}
/***********************************************************/
int BootSharedCommand::process(SharedOrderVector* order) {
	try{
			/*	EstOutput data;
				vector<SharedRAbundVector*> subset;
								
				//open an ostream for each calc to print to
				for (int z = 0; z < treeCalculators.size(); z++) {
					//create a new filename
					outputFile = outputDir + m->getRootName(m->getSimpleName(sharedfile)) + treeCalculators[z]->getName() + ".boot" + order->getLabel() + ".tre";
					m->openOutputFile(outputFile, *(out[z]));
					outputNames.push_back(outputFile); outputTypes["tree"].push_back(outputFile);
				}
				
				m->mothurOut("Generating bootstrap trees..."); cout.flush();
				
				//create a file for each calculator with the 1000 trees in it.
				for (int p = 0; p < iters; p++) {
					
					if (m->control_pressed) {  return 1; }
					
					util->getSharedVectorswithReplacement(m->Groups, lookup, order);  //fills group vectors from order vector.

				
					//for each calculator												
					for(int i = 0 ; i < treeCalculators.size(); i++) {
						
						if (m->control_pressed) {  return 1; }
						
						//initialize simMatrix
						simMatrix.clear();
						simMatrix.resize(numGroups);
						for (int o = 0; o < simMatrix.size(); o++)	{
							for (int j = 0; j < simMatrix.size(); j++)	{
								simMatrix[o].push_back(0.0);
							}
						}
				
						//initialize index
						index.clear();
						for (int g = 0; g < numGroups; g++) {	index[g] = g;	}
							
						for (int k = 0; k < lookup.size(); k++) { // pass cdd each set of groups to commpare
							for (int l = k; l < lookup.size(); l++) {
								if (k != l) { //we dont need to similiarity of a groups to itself
									subset.clear(); //clear out old pair of sharedrabunds
									//add new pair of sharedrabunds
									subset.push_back(lookup[k]); subset.push_back(lookup[l]); 
									
									//get estimated similarity between 2 groups
									data = treeCalculators[i]->getValues(subset); //saves the calculator outputs
									//save values in similarity matrix
									simMatrix[k][l] = data[0];
									simMatrix[l][k] = data[0];
								}
							}
						}
						
						tempTree = new Tree();
						
						if (m->control_pressed) {   delete tempTree; return 1; }
						
						//creates tree from similarity matrix and write out file
						createTree(out[i], tempTree);
						
						if (m->control_pressed) {   delete tempTree; return 1; }
						
						//save trees for consensus command.
						trees[i].push_back(tempTree);
					}
				}
				
				m->mothurOut("\tDone."); m->mothurOutEndLine();
								
				
				//create consensus trees for each bootstrapped tree set
				for (int k = 0; k < trees.size(); k++) {
					
					m->mothurOut("Generating consensus tree for " + treeCalculators[k]->getName()); m->mothurOutEndLine();
					
					if (m->control_pressed) {  return 1; }
					
					//set global data to calc trees
					globaldata->gTree = trees[k];
					
					string filename = outputDir + m->getRootName(m->getSimpleName(sharedfile)) + treeCalculators[k]->getName() + ".boot" + order->getLabel();
					consensus = new ConcensusCommand(filename);
					consensus->execute();
					delete consensus;
					
					outputNames.push_back(filename + ".cons.pairs");
					outputNames.push_back(filename + ".cons.tre");
					
				}
				
				
					
				//close ostream for each calc
				for (int z = 0; z < treeCalculators.size(); z++) { out[z]->close(); }
				*/
				return 0;
	
	}
	catch(exception& e) {
		m->errorOut(e, "BootSharedCommand", "process");
		exit(1);
	}
}
/***********************************************************/



