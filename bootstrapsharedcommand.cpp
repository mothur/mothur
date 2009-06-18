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

BootSharedCommand::BootSharedCommand(string option){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		lines.clear();
		labels.clear();
		Groups.clear();
		Estimators.clear();
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"line","label","calc","groups","iters"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//make sure the user has already run the read.otu command
			if (globaldata->getSharedFile() == "") {
				if (globaldata->getListFile() == "") { cout << "You must read a list and a group, or a shared before you can use the bootstrap.shared command." << endl; abort = true; }
				else if (globaldata->getGroupFile() == "") { cout << "You must read a list and a group, or a shared before you can use the bootstrap.shared command." << endl; abort = true; }
			}
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			line = validParameter.validFile(parameters, "line", false);				
			if (line == "not found") { line = "";  }
			else { 
				if(line != "all") {  splitAtDash(line, lines);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			//make sure user did not use both the line and label parameters
			if ((line != "") && (label != "")) { cout << "You cannot use both the line and label parameters at the same time. " << endl; abort = true; }
			//if the user has not specified any line or labels use the ones from read.otu
			else if((line == "") && (label == "")) {  
				allLines = globaldata->allLines; 
				labels = globaldata->labels; 
				lines = globaldata->lines;
			}
				
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
			}
				
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "jclass-thetayc";  }
			else { 
				 if (calc == "default")  {  calc = "jclass-thetayc";  }
			}
			splitAtDash(calc, Estimators);

			string temp;
			temp = validParameter.validFile(parameters, "iters", false);  if (temp == "not found") { temp = "1000"; }
			convert(temp, iters); 
				
			if (abort == false) {
			
				validCalculator = new ValidCalculators();
				
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
			}
		}

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the BootSharedCommand class Function BootSharedCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the BootSharedCommand class function BootSharedCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************

void BootSharedCommand::help(){
	try {
		cout << "The bootstrap.shared command can only be executed after a successful read.otu command." << "\n";
		cout << "The bootstrap.shared command parameters are groups, calc, iters, line and label.  You may not use line and label at the same time." << "\n";
		cout << "The groups parameter allows you to specify which of the groups in your groupfile you would like included used." << "\n";
		cout << "The group names are separated by dashes. The line and label allow you to select what distance levels you would like trees created for, and are also separated by dashes." << "\n";
		cout << "The bootstrap.shared command should be in the following format: bootstrap.shared(groups=yourGroups, calc=yourCalcs, line=yourLines, label=yourLabels, iters=yourIters)." << "\n";
		cout << "Example bootstrap.shared(groups=A-B-C, line=1-3-5, calc=jabund-sorabund, iters=100)." << "\n";
		cout << "The default value for groups is all the groups in your groupfile." << "\n";
		cout << "The default value for calc is jclass-thetayc. The default for iters is 1000." << "\n";
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the BootSharedCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the BootSharedCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************

BootSharedCommand::~BootSharedCommand(){
	//made new in execute
	if (abort == false) {
		delete input; globaldata->ginput = NULL;
		delete read;
		delete util;
		globaldata->gorder = NULL;
	}
}

//**********************************************************************************************************************

int BootSharedCommand::execute(){
	try {
	
		if (abort == true) {	return 0;	}
	
		int count = 1;
		util = new SharedUtil();	
	
		//read first line
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
		input = globaldata->ginput;
		order = input->getSharedOrderVector();
		SharedOrderVector* lastOrder = order;
		
		//if the users entered no valid calculators don't execute command
		if (treeCalculators.size() == 0) { return 0; }

		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		set<int> userLines = lines;
				
		//set users groups
		util->setGroups(globaldata->Groups, globaldata->gGroupmap->namesOfGroups, "treegroup");
		numGroups = globaldata->Groups.size();
		
		//clear globaldatas old tree names if any
		globaldata->Treenames.clear();
		
		//fills globaldatas tree names
		globaldata->Treenames = globaldata->Groups;
		
		//create treemap class from groupmap for tree class to use
		tmap = new TreeMap();
		tmap->makeSim(globaldata->gGroupmap);
		globaldata->gTreemap = tmap;
			
		while((order != NULL) && ((allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
		
			if(allLines == 1 || lines.count(count) == 1 || labels.count(order->getLabel()) == 1){			
				
				cout << order->getLabel() << '\t' << count << endl;
				process(order);
				
				processedLabels.insert(order->getLabel());
				userLabels.erase(order->getLabel());
				userLines.erase(count);
			}
			
			//you have a label the user want that is smaller than this line and the last line has not already been processed
			if ((anyLabelsToProcess(order->getLabel(), userLabels, "") == true) && (processedLabels.count(lastOrder->getLabel()) != 1)) {
											
				cout << lastOrder->getLabel() << '\t' << count << endl;
				process(lastOrder);

				processedLabels.insert(lastOrder->getLabel());
				userLabels.erase(lastOrder->getLabel());
			}
			
			if (count != 1) { delete lastOrder; }
			lastOrder = order;			

			//get next line to process
			order = input->getSharedOrderVector();
			count++;
		}
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			cout << "Your file does not include the label "<< *it; 
			if (processedLabels.count(lastOrder->getLabel()) != 1) {
				cout << ". I will use " << lastOrder->getLabel() << "." << endl;
				needToRun = true;
			}else {
				cout << ". Please refer to " << lastOrder->getLabel() << "." << endl;
			}
		}
		
		//run last line if you need to
		if (needToRun == true)  {
			process(lastOrder);			
			cout << lastOrder->getLabel() << '\t' << count << endl;
		}
		
		delete lastOrder;

		//reset groups parameter
		globaldata->Groups.clear();  

		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the BootSharedCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the BootSharedCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
//**********************************************************************************************************************

void BootSharedCommand::createTree(ostream* out){
	try {
		//create tree
		t = new Tree();
		
		//do merges and create tree structure by setting parents and children
		//there are numGroups - 1 merges to do
		for (int i = 0; i < (numGroups - 1); i++) {
		
			float largest = -1.0;
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
			simMatrix[row][column] = -1.0;
			simMatrix[column][row] = -1.0;
		
			//merge values in simsMatrix
			for (int n = 0; n < simMatrix.size(); n++)	{
				//row becomes merge of 2 groups
				simMatrix[row][n] = (simMatrix[row][n] + simMatrix[column][n]) / 2;
				simMatrix[n][row] = simMatrix[row][n];
				//delete column
				simMatrix[column][n] = -1.0;
				simMatrix[n][column] = -1.0;
			}
		}

		//assemble tree
		t->assembleTree();
	
		//print newick file
		t->print(*out);
	
		//delete tree
		delete t;
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the BootSharedCommand class Function createTree. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the BootSharedCommand class function createTree. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/***********************************************************/
void BootSharedCommand::printSims() {
	try {
		cout << "simsMatrix" << endl;
		for (int m = 0; m < simMatrix.size(); m++)	{
			for (int n = 0; n < simMatrix.size(); n++)	{
				cout << simMatrix[m][n] << '\t'; 
			}
			cout << endl;
		}

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the BootSharedCommand class Function printSims. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the BootSharedCommand class function printSims. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
/***********************************************************/
void BootSharedCommand::process(SharedOrderVector* order) {
	try{
				EstOutput data;
				vector<SharedRAbundVector*> subset;
				
				//open an ostream for each calc to print to
				for (int z = 0; z < treeCalculators.size(); z++) {
					//create a new filename
					outputFile = getRootName(globaldata->inputFileName) + treeCalculators[z]->getName() + ".boot" + order->getLabel() + ".tre";
					openOutputFile(outputFile, *(out[z]));
				}
				
				//create a file for each calculator with the 1000 trees in it.
				for (int p = 0; p < iters; p++) {
					
					util->getSharedVectorswithReplacement(Groups, lookup, order);  //fills group vectors from order vector.
				
					//for each calculator												
					for(int i = 0 ; i < treeCalculators.size(); i++) {
					
						//initialize simMatrix
						simMatrix.clear();
						simMatrix.resize(numGroups);
						for (int m = 0; m < simMatrix.size(); m++)	{
							for (int j = 0; j < simMatrix.size(); j++)	{
								simMatrix[m].push_back(0.0);
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
				
						//creates tree from similarity matrix and write out file
						createTree(out[i]);
					}
				}
				//close ostream for each calc
				for (int z = 0; z < treeCalculators.size(); z++) { out[z]->close(); }

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the BootSharedCommand class Function process. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the BootSharedCommand class function process. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
/***********************************************************/



