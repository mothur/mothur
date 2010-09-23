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
ClusterCommand::ClusterCommand(string option)  {
	try{
		globaldata = GlobalData::getInstance();
		
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"cutoff","precision","method","showabund","timing","hard","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {
					abort = true;
				}
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			//error checking to make sure they read a distance file
			if ((globaldata->gSparseMatrix == NULL) || (globaldata->gListVector == NULL)) {
				m->mothurOut("Before you use the cluster command, you first need to read in a distance matrix."); m->mothurOutEndLine();
				abort = true;
			} 
		
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			//get user cutoff and precision or use defaults
			string temp;
			temp = validParameter.validFile(parameters, "precision", false);
			if (temp == "not found") { temp = "100"; }
			//saves precision legnth for formatting below
			length = temp.length();
			convert(temp, precision); 
			
			temp = validParameter.validFile(parameters, "hard", false);			if (temp == "not found") { temp = "F"; }
			hard = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "cutoff", false);
			if (temp == "not found") { temp = "10"; }
			convert(temp, cutoff); 
			cutoff += (5 / (precision * 10.0));  
			
			method = validParameter.validFile(parameters, "method", false);
			if (method == "not found") { method = "furthest"; }
			
			if ((method == "furthest") || (method == "nearest") || (method == "average") || (method == "weighted")) { }
			else { m->mothurOut("Not a valid clustering method.  Valid clustering algorithms are furthest, nearest, average, and weighted."); m->mothurOutEndLine(); abort = true; }

			showabund = validParameter.validFile(parameters, "showabund", false);
			if (showabund == "not found") { showabund = "T"; }

			timing = validParameter.validFile(parameters, "timing", false);
			if (timing == "not found") { timing = "F"; }
			
			if (abort == false) {
			
	
							//get matrix, list and rabund for execute
				if(globaldata->gSparseMatrix != NULL)	{	matrix = globaldata->gSparseMatrix;		}
			
				if(globaldata->gListVector != NULL){
					list = globaldata->gListVector;
					rabund = new RAbundVector(list->getRAbundVector());
				}
				
				//create cluster
				if (method == "furthest")	{	cluster = new CompleteLinkage(rabund, list, matrix, cutoff, method); }
				else if(method == "nearest"){	cluster = new SingleLinkage(rabund, list, matrix, cutoff, method); }
				else if(method == "average"){	cluster = new AverageLinkage(rabund, list, matrix, cutoff, method);	}
				else if(method == "weighted"){	cluster = new WeightedLinkage(rabund, list, matrix, cutoff, method);	}
				tag = cluster->getTag();
				
				if (outputDir == "") { outputDir += m->hasPath(globaldata->inputFileName); }
				fileroot = outputDir + m->getRootName(m->getSimpleName(globaldata->inputFileName));
			
				m->openOutputFile(fileroot+ tag + ".sabund",	sabundFile);
				m->openOutputFile(fileroot+ tag + ".rabund",	rabundFile);
				m->openOutputFile(fileroot+ tag + ".list",		listFile);
				
				outputNames.push_back(fileroot+ tag + ".sabund");
				outputNames.push_back(fileroot+ tag + ".rabund");
				outputNames.push_back(fileroot+ tag + ".list");
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterCommand", "ClusterCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void ClusterCommand::help(){
	try {
		m->mothurOut("The cluster command can only be executed after a successful read.dist command.\n");
		m->mothurOut("The cluster command parameter options are method, cuttoff, hard, precision, showabund and timing. No parameters are required.\n");
		m->mothurOut("The cluster command should be in the following format: \n");
		m->mothurOut("cluster(method=yourMethod, cutoff=yourCutoff, precision=yourPrecision) \n");
		m->mothurOut("The acceptable cluster methods are furthest, nearest and average.  If no method is provided then furthest is assumed.\n\n");	
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

ClusterCommand::~ClusterCommand(){
	if (abort == false) {
		delete cluster;
		delete rabund;
	}
}

//**********************************************************************************************************************

int ClusterCommand::execute(){
	try {
	
		if (abort == true) {	return 0;	}
		
		time_t estart = time(NULL);
		//int ndist = matrix->getNNodes();
		float previousDist = 0.00000;
		float rndPreviousDist = 0.00000;
		oldRAbund = *rabund;
		oldList = *list;

		print_start = true;
		start = time(NULL);
		loops = 0;
		double saveCutoff = cutoff;
		
		while (matrix->getSmallDist() < cutoff && matrix->getNNodes() > 0){
	cout << matrix->getSmallDist() << '\t' << cutoff << '\t' << matrix->getNNodes() << endl;	
			if (m->control_pressed) { //clean up
				delete globaldata->gSparseMatrix;  globaldata->gSparseMatrix = NULL;
				delete globaldata->gListVector;	 globaldata->gListVector = NULL;
				if (globaldata->getFormat() == "phylip") { globaldata->setPhylipFile(""); }
				else if (globaldata->getFormat() == "column") { globaldata->setColumnFile(""); }
				sabundFile.close();rabundFile.close();listFile.close();
				for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); 	}
				return 0;
			}
		
			if (print_start && m->isTrue(timing)) {
				m->mothurOut("Clustering (" + tag + ") dist " + toString(matrix->getSmallDist()) + "/" 
					+ toString(m->roundDist(matrix->getSmallDist(), precision)) 
					+ "\t(precision: " + toString(precision) + ", Nodes: " + toString(matrix->getNNodes()) + ")");
				cout.flush();
				print_start = false;
			}

			loops++;

			cluster->update(cutoff);
	
			float dist = matrix->getSmallDist();
			float rndDist;
			if (hard) {
				rndDist = m->ceilDist(dist, precision); 
			}else{
				rndDist = m->roundDist(dist, precision); 
			}

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

		if (print_start && m->isTrue(timing)) {
			m->mothurOut("Clustering (" + tag + ") for distance " + toString(previousDist) + "/" + toString(rndPreviousDist) 
					 + "\t(precision: " + toString(precision) + ", Nodes: " + toString(matrix->getNNodes()) + ")");
			cout.flush();
	 		print_start = false;
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
		
		sabundFile.close();
		rabundFile.close();
		listFile.close();
	
		if (saveCutoff != cutoff) { 
			if (hard)	{  saveCutoff = m->ceilDist(saveCutoff, precision);	}
			else		{	saveCutoff = m->roundDist(saveCutoff, precision);  }

			m->mothurOut("changed cutoff to " + toString(cutoff)); m->mothurOutEndLine(); 
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

		
		//if (m->isTrue(timing)) {
			m->mothurOut("It took " + toString(time(NULL) - estart) + " seconds to cluster"); m->mothurOutEndLine();
		//}
		
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************

void ClusterCommand::printData(string label){
	try {
		if (m->isTrue(timing)) {
			m->mothurOut("\tTime: " + toString(time(NULL) - start) + "\tsecs for " + toString(oldRAbund.getNumBins()) 
		     + "\tclusters. Updates: " + toString(loops)); m->mothurOutEndLine();
		}
		print_start = true;
		loops = 0;
		start = time(NULL);

		oldRAbund.setLabel(label);
		if (m->isTrue(showabund)) {
			oldRAbund.getSAbundVector().print(cout);
		}
		oldRAbund.print(rabundFile);
		oldRAbund.getSAbundVector().print(sabundFile);
	
		oldList.setLabel(label);
		oldList.print(listFile);
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterCommand", "printData");
		exit(1);
	}


}
//**********************************************************************************************************************
