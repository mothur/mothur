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
ClusterCommand::ClusterCommand(){
	try{
		globaldata = GlobalData::getInstance();
	
		if(globaldata->getSparseMatrix() != NULL)	{	matrix = new SparseMatrix(*globaldata->getSparseMatrix());		}
	//  Not sure if we want the address or an entire new memory allocation.  Might be nice to have new memory so data
	//  doesn't need to be re-read, but then again, it could suck up a ton of memory.  Dunno.
	//	if(globaldata->getSparseMatrix() != NULL)	{	matrix = globaldata->getSparseMatrix();		}
	
		if(globaldata->getListVector() != NULL){
			list = new ListVector(*globaldata->getListVector());
			rabund = new RAbundVector(list->getRAbundVector());
			//rabund->print(cout);
		}
	
		if(globaldata->getMethod() != "")	{	method = globaldata->getMethod();		}		
		//if no method given use furthest, initialized in globaldata
		if(method == "furthest")	{	cluster = new CompleteLinkage(rabund, list, matrix);	tag = "fn";	}
		else if(method == "nearest"){	cluster = new SingleLinkage(rabund, list, matrix);		tag = "nn";	}
		else if(method == "average"){	cluster = new AverageLinkage(rabund, list, matrix);		tag = "an";	}
		else						{	cout << "error - not recognized method" << endl;											}
	
		if(globaldata->getPrecision() != ""){
			convert(globaldata->getPrecision(), precision);	
		}
		
		//saves precision legnth for formatting below
		length = globaldata->getPrecision().length();
		
		if(globaldata->getCutOff() != ""){
			convert(globaldata->getCutOff(), cutoff);	
			cutoff += (5 / (precision * 10.0));
		}
	
		fileroot = getRootName(globaldata->getFileRoot());
		
		openOutputFile(fileroot+ tag + ".sabund",	sabundFile);
		openOutputFile(fileroot+ tag + ".rabund",	rabundFile);
		openOutputFile(fileroot+ tag + ".list",		listFile);
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

ClusterCommand::~ClusterCommand(){
	delete cluster;
	delete matrix;
	delete list;
	delete rabund;
}

//**********************************************************************************************************************

int ClusterCommand::execute(){
	try {
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
		SparseMatrix* clearM = NULL;
		globaldata->setSparseMatrix(clearM);
		ListVector* clearL = NULL;
		globaldata->setListVector(clearL);

		
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
