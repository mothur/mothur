/*
 *  hclustercommand.cpp
 *  Mothur
 *
 *  Created by westcott on 10/13/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "hclustercommand.h"

//**********************************************************************************************************************
//This function checks to make sure the cluster command has no errors and then clusters based on the method chosen.
HClusterCommand::HClusterCommand(string option){
	try{
		globaldata = GlobalData::getInstance();
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"cutoff","precision","method","phylip","column","name","sorted","showabund","timing","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {
					abort = true;
				}
			}
			
			globaldata->newRead();
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("phylip");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["phylip"] = inputDir + it->second;		}
				}
				
				it = parameters.find("column");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["column"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
			}

			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			//check for required parameters
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not open") { abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; }	
			else {  distfile = phylipfile;  format = "phylip"; 	}
			
			columnfile = validParameter.validFile(parameters, "column", true);
			if (columnfile == "not open") { abort = true; }	
			else if (columnfile == "not found") { columnfile = ""; }
			else {  distfile = columnfile; format = "column";	}
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			
			if ((phylipfile == "") && (columnfile == "")) { mothurOut("When executing a hcluster command you must enter a phylip or a column."); mothurOutEndLine(); abort = true; }
			else if ((phylipfile != "") && (columnfile != "")) { mothurOut("When executing a hcluster command you must enter ONLY ONE of the following: phylip or column."); mothurOutEndLine(); abort = true; }
		
			if (columnfile != "") {
				if (namefile == "") {  cout << "You need to provide a namefile if you are going to use the column format." << endl; abort = true; }
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
			
			temp = validParameter.validFile(parameters, "cutoff", false);
			if (temp == "not found") { temp = "10"; }
			convert(temp, cutoff); 
			cutoff += (5 / (precision * 10.0));
			
			method = validParameter.validFile(parameters, "method", false);
			if (method == "not found") { method = "furthest"; }
			
			if ((method == "furthest") || (method == "nearest") || (method == "average")) { }
			else { mothurOut("Not a valid clustering method.  Valid clustering algorithms are furthest, nearest or average."); mothurOutEndLine(); abort = true; }

			showabund = validParameter.validFile(parameters, "showabund", false);
			if (showabund == "not found") { showabund = "T"; }
			
			sort = validParameter.validFile(parameters, "sorted", false);
			if (sort == "not found") { sort = "F"; }
			sorted = isTrue(sort);

			timing = validParameter.validFile(parameters, "timing", false);
			if (timing == "not found") { timing = "F"; }
			
				
			if (abort == false) {
				
				if (outputDir == "") {  outputDir += hasPath(distfile); }
				fileroot = outputDir + getRootName(getSimpleName(distfile));
				
				if (method == "furthest")		{ tag = "fn";  }
				else if (method == "nearest")	{ tag = "nn";  }
				else							{ tag = "an";  }
			
				openOutputFile(fileroot+ tag + ".sabund",	sabundFile);
				openOutputFile(fileroot+ tag + ".rabund",	rabundFile);
				openOutputFile(fileroot+ tag + ".list",		listFile);
			}
		}
	}
	catch(exception& e) {
		errorOut(e, "HClusterCommand", "HClusterCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void HClusterCommand::help(){
	try {
		mothurOut("The hcluster command parameter options are cutoff, precision, method, phylip, column, name, showabund, timing and sorted. Phylip or column and name are required.\n");
		mothurOut("The phylip and column parameter allow you to enter your distance file, and sorted indicates whether your column distance file is already sorted. \n");
		mothurOut("The name parameter allows you to enter your name file and is required if your distance file is in column format. \n");
		mothurOut("The hcluster command should be in the following format: \n");
		mothurOut("hcluster(column=youDistanceFile, name=yourNameFile, method=yourMethod, cutoff=yourCutoff, precision=yourPrecision) \n");
		mothurOut("The acceptable hcluster methods are furthest and nearest, but we hope to add average in the future.\n\n");	
	}
	catch(exception& e) {
		errorOut(e, "HClusterCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

HClusterCommand::~HClusterCommand(){}

//**********************************************************************************************************************

int HClusterCommand::execute(){
	try {
	
		if (abort == true) {	return 0;	}
		
		if(namefile != ""){	
			globaldata->nameMap = new NameAssignment(namefile);
			globaldata->nameMap->readMap();
		}else{
			globaldata->nameMap = NULL;
		}
		
		time_t estart = time(NULL);
		
		if (!sorted) {
			read = new ReadCluster(distfile, cutoff); 	
			read->setFormat(format);
			read->read(globaldata->nameMap);
			distfile = read->getOutputFile();
		
			list = read->getListVector();
			delete read;
		}else {
			list = new ListVector(globaldata->nameMap->getListVector());
		}
	
		mothurOut("It took " + toString(time(NULL) - estart) + " seconds to sort. "); mothurOutEndLine();
		estart = time(NULL);
	
		//list vector made by read contains all sequence names
		if(list != NULL){
			rabund = new RAbundVector(list->getRAbundVector());
		}else{
			mothurOut("Error: no list vector!"); mothurOutEndLine(); return 0;
		}
		
		float previousDist = 0.00000;
		float rndPreviousDist = 0.00000;
		oldRAbund = *rabund;
		oldList = *list;
		
		print_start = true;
		start = time(NULL);
				
		cluster = new HCluster(rabund, list, method, distfile, globaldata->nameMap, cutoff);
		vector<seqDist> seqs; seqs.resize(1); // to start loop
		
		while (seqs.size() != 0){
		
			seqs = cluster->getSeqs();
				
			for (int i = 0; i < seqs.size(); i++) {  //-1 means skip me

				if (seqs[i].seq1 != seqs[i].seq2) {
					cluster->update(seqs[i].seq1, seqs[i].seq2, seqs[i].dist);
					
					float rndDist = roundDist(seqs[i].dist, precision);
					
					if((previousDist <= 0.0000) && (seqs[i].dist != previousDist)){
						printData("unique");
					}
					else if((rndDist != rndPreviousDist)){
						printData(toString(rndPreviousDist,  length-1));
					}
				
					previousDist = seqs[i].dist;
					rndPreviousDist = rndDist;
					oldRAbund = *rabund;
					oldList = *list;
				}
			}
		}

		if(previousDist <= 0.0000){
			printData("unique");
		}
		else if(rndPreviousDist<cutoff){
			printData(toString(rndPreviousDist, length-1));
		}
		
		//delete globaldata's copy of the sparsematrix and listvector to free up memory
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
		delete cluster;
	
		mothurOut("It took " + toString(time(NULL) - estart) + " seconds to cluster. "); mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "HClusterCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************

void HClusterCommand::printData(string label){
	try {
		if (isTrue(timing)) {
			mothurOut("\tTime: " + toString(time(NULL) - start) + "\tsecs for " + toString(oldRAbund.getNumBins()) 
		     + "\tclusters. Updates: " + toString(loops)); mothurOutEndLine();
		}
		print_start = true;
		loops = 0;
		start = time(NULL);

		oldRAbund.setLabel(label);
		if (isTrue(showabund)) {
			oldRAbund.getSAbundVector().print(cout);
		}
		oldRAbund.print(rabundFile);
		oldRAbund.getSAbundVector().print(sabundFile);
	
		oldList.setLabel(label);
		oldList.print(listFile);
	}
	catch(exception& e) {
		errorOut(e, "HClusterCommand", "printData");
		exit(1);
	}


}
//**********************************************************************************************************************

