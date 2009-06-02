/*
 *  treegroupscommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "treegroupscommand.h"
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

TreeGroupCommand::TreeGroupCommand(){
	try {
		globaldata = GlobalData::getInstance();
		format = globaldata->getFormat();
		validCalculator = new ValidCalculators();
		
		if (format == "sharedfile") {
			int i;
			for (i=0; i<globaldata->Estimators.size(); i++) {
				if (validCalculator->isValidCalculator("treegroup", globaldata->Estimators[i]) == true) { 
					if (globaldata->Estimators[i] == "jabund") { 	
						treeCalculators.push_back(new JAbund());
					}else if (globaldata->Estimators[i] == "sorabund") { 
						treeCalculators.push_back(new SorAbund());
					}else if (globaldata->Estimators[i] == "jclass") { 
						treeCalculators.push_back(new Jclass());
					}else if (globaldata->Estimators[i] == "sorclass") { 
						treeCalculators.push_back(new SorClass());
					}else if (globaldata->Estimators[i] == "jest") { 
						treeCalculators.push_back(new Jest());
					}else if (globaldata->Estimators[i] == "sorest") { 
						treeCalculators.push_back(new SorEst());
					}else if (globaldata->Estimators[i] == "thetayc") { 
						treeCalculators.push_back(new ThetaYC());
					}else if (globaldata->Estimators[i] == "thetan") { 
						treeCalculators.push_back(new ThetaN());
					}else if (globaldata->Estimators[i] == "morisitahorn") { 
						treeCalculators.push_back(new MorHorn());
					}else if (globaldata->Estimators[i] == "braycurtis") { 
						treeCalculators.push_back(new BrayCurtis());
					}
				}
			}
		}
		
		//reset calc for next command
		globaldata->setCalc("");

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the TreeGroupCommand class Function TreeGroupCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the TreeGroupCommand class function TreeGroupCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************

TreeGroupCommand::~TreeGroupCommand(){
	delete input;
	if (format == "sharedfile") {delete read;}
	else { delete readMatrix;  delete matrix; delete list; }
	delete tmap;
	
}

//**********************************************************************************************************************

int TreeGroupCommand::execute(){
	try {
		if (format == "sharedfile") {
			//if the users entered no valid calculators don't execute command
			if (treeCalculators.size() == 0) { cout << "You have given no valid calculators." << endl; return 0; }

			//you have groups
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
			
			input = globaldata->ginput;
			lookup = input->getSharedRAbundVectors();
			lastLookup = lookup;
			
			if (lookup.size() < 2) { cout << "You have not provided enough valid groups.  I cannot run the command." << endl; return 0; }
		
			//create tree file
			makeSimsShared();
		}else{
			//read in dist file
			filename = globaldata->inputFileName;
		
			if (format == "column") { readMatrix = new ReadColumnMatrix(filename); }	
			else if (format == "phylip") { readMatrix = new ReadPhylipMatrix(filename); }
				
			if(globaldata->getPrecision() != ""){
				convert(globaldata->getPrecision(), precision);	
			}
		
			if(globaldata->getCutOff() != ""){
				convert(globaldata->getCutOff(), cutoff);	
				cutoff += (5 / (precision * 10.0));
			}
			readMatrix->setCutoff(cutoff);
	
			if(globaldata->getNameFile() != ""){	
				nameMap = new NameAssignment(globaldata->getNameFile());
				nameMap->readMap(1,2);
			}
			else{
				nameMap = NULL;
			}
	
			readMatrix->read(nameMap);
			list = readMatrix->getListVector();
			matrix = readMatrix->getMatrix();

			//make treemap
			tmap = new TreeMap();
			tmap->makeSim(list);
			globaldata->gTreemap = tmap;
			
			globaldata->Groups = tmap->namesOfGroups;
		
			//clear globaldatas old tree names if any
			globaldata->Treenames.clear();
		
			//fills globaldatas tree names
			globaldata->Treenames = globaldata->Groups;

			makeSimsDist();

			//create a new filename
			outputFile = getRootName(globaldata->inputFileName) + "tre";	
				
			createTree();
		}
				
		//reset groups parameter
		globaldata->Groups.clear();  globaldata->setGroups("");

		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the TreeGroupCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the TreeGroupCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
//**********************************************************************************************************************

void TreeGroupCommand::createTree(){
	try {
		//create tree
		t = new Tree();
		
		//do merges and create tree structure by setting parents and children
		//there are numGroups - 1 merges to do
		for (int i = 0; i < (numGroups - 1); i++) {
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
			
			//remove highest value that caused the merge.
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
		
		//print newick file
		t->createNewickFile(outputFile);
		
		//delete tree
		delete t;
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the TreeGroupCommand class Function createTree. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the TreeGroupCommand class function createTree. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/***********************************************************/
void TreeGroupCommand::printSims(ostream& out) {
	try {
		
		//output column headers
		//out << '\t';
		//for (int i = 0; i < lookup.size(); i++) {	out << lookup[i]->getGroup() << '\t';		}
		//out << endl;
		
		
		for (int m = 0; m < simMatrix.size(); m++)	{
			//out << lookup[m]->getGroup() << '\t';
			for (int n = 0; n < simMatrix.size(); n++)	{
				out << simMatrix[m][n] << '\t'; 
			}
			out << endl;
		}

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the TreeGroupCommand class Function printSims. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the TreeGroupCommand class function printSims. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
/***********************************************************/
void TreeGroupCommand::makeSimsDist() {
	try {
		numGroups = list->size();
		
		//initialize index
		index.clear();
		for (int g = 0; g < numGroups; g++) {	index[g] = g;	}
		
		//initialize simMatrix
		simMatrix.clear();
		simMatrix.resize(numGroups);
		for (int m = 0; m < simMatrix.size(); m++)	{
			for (int j = 0; j < simMatrix.size(); j++)	{
				simMatrix[m].push_back(0.0);
			}
		}
		
		//go through sparse matrix and fill sims
		//go through each cell in the sparsematrix
		for(MatData currentCell = matrix->begin(); currentCell != matrix->end(); currentCell++){
			//similairity = -(distance-1)
			simMatrix[currentCell->row][currentCell->column] = -(currentCell->dist -1.0);	
			simMatrix[currentCell->column][currentCell->row] = -(currentCell->dist -1.0);				
		}


	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the TreeGroupCommand class Function makeSimsDist. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the TreeGroupCommand class function makeSimsDist. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

/***********************************************************/
void TreeGroupCommand::makeSimsShared() {
	try {
		int count = 1;	
	
		//clear globaldatas old tree names if any
		globaldata->Treenames.clear();
		
		//fills globaldatas tree names
		globaldata->Treenames = globaldata->Groups;
		
		//create treemap class from groupmap for tree class to use
		tmap = new TreeMap();
		tmap->makeSim(globaldata->gGroupmap);
		globaldata->gTreemap = tmap;
		
		set<string> processedLabels;
		set<string> userLabels = globaldata->labels;

		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((globaldata->allLines == 1) || (userLabels.size() != 0))) {
		
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(lookup[0]->getLabel()) == 1){			
				cout << lookup[0]->getLabel() << '\t' << count << endl;
				process(lookup);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}
			
			if ((anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLookup[0]->getLabel()) != 1)) {
				cout << lastLookup[0]->getLabel() << '\t' << count << endl;
				process(lastLookup);
					
				processedLabels.insert(lastLookup[0]->getLabel());
				userLabels.erase(lastLookup[0]->getLabel());
			}

			//prevent memory leak
			if (count != 1) { for (int i = 0; i < lastLookup.size(); i++) {  delete lastLookup[i];  } }
			lastLookup = lookup;			
			
			//get next line to process
			lookup = input->getSharedRAbundVectors();
			count++;
		}
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			cout << "Your file does not include the label "<< *it; 
			if (processedLabels.count(lastLookup[0]->getLabel()) != 1) {
				cout << ". I will use " << lastLookup[0]->getLabel() << "." << endl;
				needToRun = true;
			}else {
				cout << ". Please refer to " << lastLookup[0]->getLabel() << "." << endl;
			}
		}
		
		//run last line if you need to
		if (needToRun == true)  {
			cout << lastLookup[0]->getLabel() << '\t' << count << endl;
			process(lastLookup);
		}
		
		for (int i = 0; i < lastLookup.size(); i++) {  delete lastLookup[i];  }
		for(int i = 0 ; i < treeCalculators.size(); i++) {  delete treeCalculators[i]; }
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the TreeGroupCommand class Function makeSimsShared. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the TreeGroupCommand class function makeSimsShared. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

/***********************************************************/
void TreeGroupCommand::process(vector<SharedRAbundVector*> thisLookup) {
	try{
				EstOutput data;
				vector<SharedRAbundVector*> subset;
				numGroups = globaldata->Groups.size();
				
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
		
					//create a new filename
					outputFile = getRootName(globaldata->inputFileName) + treeCalculators[i]->getName() + "." + thisLookup[0]->getLabel() + ".tre";				
												
					for (int k = 0; k < thisLookup.size(); k++) { 
						for (int l = k; l < thisLookup.size(); l++) {
							if (k != l) { //we dont need to similiarity of a groups to itself
								//get estimated similarity between 2 groups
								
								subset.clear(); //clear out old pair of sharedrabunds
								//add new pair of sharedrabunds
								subset.push_back(thisLookup[k]); subset.push_back(thisLookup[l]); 
								
								data = treeCalculators[i]->getValues(subset); //saves the calculator outputs
								//save values in similarity matrix
								simMatrix[k][l] = data[0];
								simMatrix[l][k] = data[0];
							}
						}
					}
					
					//creates tree from similarity matrix and write out file
					createTree();
				}

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the TreeGroupCommand class Function process. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the TreeGroupCommand class function process. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
/***********************************************************/

	

