/*
 *  clustercommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "clustercommand.h"

//**********************************************************************************************************************
//This function checks to make sure the cluster command has no errors and then clusters based on the method chosen.
ClusterCommand::ClusterCommand(string option){
	try{
		globaldata = GlobalData::getInstance();
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"cutoff","precision","method"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			parser = new OptionParser();
			parser->parse(option, parameters);  delete parser;
			
			ValidParameters* validParameter = new ValidParameters();
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter->isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//error checking to make sure they read a distance file
			if ((globaldata->gSparseMatrix == NULL) || (globaldata->gListVector == NULL)) {
				cout << "Before you use the cluster command, you first need to read in a distance matrix." << endl;  abort = true;
			} 
		
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			//get user cutoff and precision or use defaults
			string temp;
			temp = validParameter->validFile(parameters, "precision", false);		if (temp == "not found") { temp = "100"; }
			//saves precision legnth for formatting below
			length = temp.length();
			convert(temp, precision); 
			
			temp = validParameter->validFile(parameters, "cutoff", false);			if (temp == "not found") { temp = "10"; }
			convert(temp, cutoff); 
			cutoff += (5 / (precision * 10.0));
			
			method = validParameter->validFile(parameters, "method", false);			if (method == "not found") { method = "furthest"; }

			delete validParameter;
			
			if ((method == "furthest") || (method == "nearest") || (method == "average")) { }
			else {cout << "Not a valid clustering method.  Valid clustering algorithms are furthest, nearest or average." << endl; abort = true; }

			
			if (abort == false) {
			
				//get matrix, list and rabund for execute
				if(globaldata->gSparseMatrix != NULL)	{	matrix = new SparseMatrix(*globaldata->gSparseMatrix);		}
			
				if(globaldata->gListVector != NULL){
					list = new ListVector(*globaldata->gListVector);
					rabund = new RAbundVector(list->getRAbundVector());
				}
				
				//create cluster
				if(method == "furthest")	{	cluster = new CompleteLinkage(rabund, list, matrix);	tag = "fn";	}
				else if(method == "nearest"){	cluster = new SingleLinkage(rabund, list, matrix);		tag = "nn";	}
				else if(method == "average"){	cluster = new AverageLinkage(rabund, list, matrix);		tag = "an";	}
				else						{	cout << "error - not recognized method" << endl;	abort = true;	}	
				
				fileroot = getRootName(globaldata->inputFileName);
			
				openOutputFile(fileroot+ tag + ".sabund",	sabundFile);
				openOutputFile(fileroot+ tag + ".rabund",	rabundFile);
				openOutputFile(fileroot+ tag + ".list",		listFile);
				
				
			}

		}
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ClusterCommand class Function ClusterCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ClusterCommand class function ClusterCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

//**********************************************************************************************************************

void ClusterCommand::help(){
	try {
		cout << "The cluster command can only be executed after a successful read.dist command." << "\n";
		cout << "The cluster command parameter options are method, cuttoff and precision. No parameters are required." << "\n";
		cout << "The cluster command should be in the following format: " << "\n";
		cout << "cluster(method=yourMethod, cutoff=yourCutoff, precision=yourPrecision) " << "\n";
		cout << "The acceptable cluster methods are furthest, nearest and average.  If no method is provided then furthest is assumed." << "\n" << "\n";	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ClusterCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ClusterCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************

ClusterCommand::~ClusterCommand(){
	delete cluster;
	delete matrix;
	delete list;
	delete rabund;
}

//**********************************************************************************************************************

int ClusterCommand::execute(){
	try {
	
		if (abort == true) {	return 0;	}
		
		float previousDist = 0.00000;
		float rndPreviousDist = 0.00000;
		oldRAbund = *rabund;
		oldList = *list;
		
		float x;
		x=0.1;
		toString(x, 2);
	
		while(matrix->getSmallDist() < cutoff && matrix->getNNodes() > 0){
			cluster->update();
			float dist = matrix->getSmallDist();
			float rndDist = roundDist(dist, precision);

			if(previousDist <= 0.0000 && dist != previousDist){
				printData("unique");
			}
			else if(rndDist != rndPreviousDist){
				printData(toString(rndPreviousDist,  length-1));
			}
		
			previousDist = dist;
			rndPreviousDist = rndDist;
			oldRAbund = *rabund;
			oldList = *list;
		}
	
		if(previousDist <= 0.0000){
			printData("unique");
		}
		else if(rndPreviousDist<cutoff){
			printData(toString(rndPreviousDist, length-1));
		}
		
		//delete globaldata's copy of the sparsematrix and listvector to free up memory
		delete globaldata->gSparseMatrix;  globaldata->gSparseMatrix = NULL;
		delete globaldata->gListVector;	 globaldata->gListVector = NULL;
		
		//saves .list file so you can do the collect, rarefaction and summary commands without doing a read.list
		if (globaldata->getFormat() == "phylip") { globaldata->setPhylipFile(""); }
		else if (globaldata->getFormat() == "column") { globaldata->setColumnFile(""); }
		
		globaldata->setListFile(fileroot+ tag + ".list");
		globaldata->setNameFile("");
		globaldata->setFormat("list");
		
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ClusterCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ClusterCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

//**********************************************************************************************************************

void ClusterCommand::printData(string label){
	try {
		oldRAbund.setLabel(label);
		oldRAbund.getSAbundVector().print(cout);
		oldRAbund.print(rabundFile);
		oldRAbund.getSAbundVector().print(sabundFile);
	
		oldList.setLabel(label);
		oldList.print(listFile);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ClusterCommand class Function printData. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ClusterCommand class function printData. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
//**********************************************************************************************************************
