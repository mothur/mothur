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


//**********************************************************************************************************************

TreeGroupCommand::TreeGroupCommand(){
	try {
		globaldata = GlobalData::getInstance();
		format = globaldata->getFormat();
		validCalculator = new ValidCalculators();
		util = new SharedUtil();
		
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
	delete read;
	delete util;
}

//**********************************************************************************************************************

int TreeGroupCommand::execute(){
	try {
		int count = 1;	
		EstOutput data;
	
		//if the users entered no valid calculators don't execute command
		if (treeCalculators.size() == 0) { return 0; }

		if (format == "sharedfile") {
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
			
			input = globaldata->ginput;
			order = input->getSharedOrderVector();
		}else {
			//you are using a list and a groupfile
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
		
			input = globaldata->ginput;
			SharedList = globaldata->gSharedList;
			order = SharedList->getSharedOrderVector();
		}
		
		//set users groups
		util->setGroups(globaldata->Groups, globaldata->gGroupmap->namesOfGroups, "treegroup");
		numGroups = globaldata->Groups.size();
		groupNames = "";
		for (int i = 0; i < numGroups; i++) { groupNames += globaldata->Groups[i]; }
		
		//clear globaldatas old tree names if any
		globaldata->Treenames.clear();
		
		//fills globaldatas tree names
		globaldata->Treenames = globaldata->Groups;
		
		//create treemap class from groupmap for tree class to use
		tmap = new TreeMap();
		tmap->makeSim(globaldata->gGroupmap);
		globaldata->gTreemap = tmap;
			
		while(order != NULL){
		
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(order->getLabel()) == 1){			
				
				cout << order->getLabel() << '\t' << count << endl;
				util->getSharedVectors(globaldata->Groups, lookup, order);  //fills group vectors from order vector.
				
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
					outputFile = getRootName(globaldata->inputFileName) + treeCalculators[i]->getName() + "." + order->getLabel() + ".tre";
							
					for (int k = 0; k < lookup.size(); k++) { 
						for (int l = k; l < lookup.size(); l++) {
							if (k != l) { //we dont need to similiarity of a groups to itself
								//get estimated similarity between 2 groups
								data = treeCalculators[i]->getValues(lookup[k], lookup[l]); //saves the calculator outputs
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
		
			//get next line to process
			if (format == "sharedfile") {
				order = input->getSharedOrderVector();
			}else {
				//you are using a list and a groupfile
				SharedList = input->getSharedListVector(); //get new list vector to process
				if (SharedList != NULL) {
					order = SharedList->getSharedOrderVector(); //gets new order vector with group info.
				}else {
					break;
				}
			}
			count++;
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
			
			float largest = 0.0;
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
			simMatrix[row][column] = 0.0;
			simMatrix[column][row] = 0.0;
			
			//merge values in simsMatrix
			for (int n = 0; n < simMatrix.size(); n++)	{
				//row becomes merge of 2 groups
				simMatrix[row][n] = (simMatrix[row][n] + simMatrix[column][n]) / 2;
				simMatrix[n][row] = simMatrix[row][n];
				//delete column
				simMatrix[column][n] = 0.0;
				simMatrix[n][column] = 0.0;
			}
		}
	
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
void TreeGroupCommand::printSims() {
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
		cout << "Standard Error: " << e.what() << " has occurred in the TreeGroupCommand class Function printSims. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the TreeGroupCommand class function printSims. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
/***********************************************************/


