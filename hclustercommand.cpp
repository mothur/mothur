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
vector<string> HClusterCommand::getValidParameters(){	
	try {
		string Array[] =  {"cutoff","hard","precision","method","phylip","column","name","sorted","showabund","timing","outputdir","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "HClusterCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
HClusterCommand::HClusterCommand(){	
	try {
		abort = true;
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["list"] = tempOutNames;
		outputTypes["rabund"] = tempOutNames;
		outputTypes["sabund"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "HClusterCommand", "HClusterCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> HClusterCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"phylip","column","or"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "HClusterCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> HClusterCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "HClusterCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************
//This function checks to make sure the cluster command has no errors and then clusters based on the method chosen.
HClusterCommand::HClusterCommand(string option)  {
	try{
		globaldata = GlobalData::getInstance();
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"cutoff","hard","precision","method","phylip","column","name","sorted","showabund","timing","outputdir","inputdir"};
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
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["list"] = tempOutNames;
			outputTypes["rabund"] = tempOutNames;
			outputTypes["sabund"] = tempOutNames;
		
			globaldata->newRead();
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("phylip");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["phylip"] = inputDir + it->second;		}
				}
				
				it = parameters.find("column");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["column"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
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
			
			if ((phylipfile == "") && (columnfile == "")) { m->mothurOut("When executing a hcluster command you must enter a phylip or a column."); m->mothurOutEndLine(); abort = true; }
			else if ((phylipfile != "") && (columnfile != "")) { m->mothurOut("When executing a hcluster command you must enter ONLY ONE of the following: phylip or column."); m->mothurOutEndLine(); abort = true; }
		
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
			
			temp = validParameter.validFile(parameters, "hard", false);			if (temp == "not found") { temp = "F"; }
			hard = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "cutoff", false);
			if (temp == "not found") { temp = "10"; }
			convert(temp, cutoff); 
			cutoff += (5 / (precision * 10.0)); 
			
			method = validParameter.validFile(parameters, "method", false);
			if (method == "not found") { method = "furthest"; }
			
			if ((method == "furthest") || (method == "nearest") || (method == "average")) { }
			else { m->mothurOut("Not a valid clustering method.  Valid clustering algorithms are furthest, nearest or average."); m->mothurOutEndLine(); abort = true; }

			showabund = validParameter.validFile(parameters, "showabund", false);
			if (showabund == "not found") { showabund = "T"; }
			
			sort = validParameter.validFile(parameters, "sorted", false);
			if (sort == "not found") { sort = "F"; }
			sorted = m->isTrue(sort);

			timing = validParameter.validFile(parameters, "timing", false);
			if (timing == "not found") { timing = "F"; }
			
				
			if (abort == false) {
				
				if (outputDir == "") {  outputDir += m->hasPath(distfile); }
				fileroot = outputDir + m->getRootName(m->getSimpleName(distfile));
				
				if (method == "furthest")		{ tag = "fn";  }
				else if (method == "nearest")	{ tag = "nn";  }
				else							{ tag = "an";  }
			
				m->openOutputFile(fileroot+ tag + ".sabund",	sabundFile);
				m->openOutputFile(fileroot+ tag + ".rabund",	rabundFile);
				m->openOutputFile(fileroot+ tag + ".list",		listFile);
				
				outputNames.push_back(fileroot+ tag + ".sabund"); outputTypes["sabund"].push_back(fileroot+ tag + ".sabund");
				outputNames.push_back(fileroot+ tag + ".rabund"); outputTypes["rabund"].push_back(fileroot+ tag + ".rabund");
				outputNames.push_back(fileroot+ tag + ".list"); outputTypes["list"].push_back(fileroot+ tag + ".list");
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "HClusterCommand", "HClusterCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void HClusterCommand::help(){
	try {
		m->mothurOut("The hcluster command parameter options are cutoff, precision, method, phylip, column, name, showabund, timing and sorted. Phylip or column and name are required.\n");
		m->mothurOut("The phylip and column parameter allow you to enter your distance file, and sorted indicates whether your column distance file is already sorted. \n");
		m->mothurOut("The name parameter allows you to enter your name file and is required if your distance file is in column format. \n");
		m->mothurOut("The hcluster command should be in the following format: \n");
		m->mothurOut("hcluster(column=youDistanceFile, name=yourNameFile, method=yourMethod, cutoff=yourCutoff, precision=yourPrecision) \n");
		m->mothurOut("The acceptable hcluster methods are furthest, nearest and average.\n\n");	
	}
	catch(exception& e) {
		m->errorOut(e, "HClusterCommand", "help");
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
			read = new ReadCluster(distfile, cutoff, outputDir, true); 	
			read->setFormat(format);
			read->read(globaldata->nameMap);
			
			if (m->control_pressed) {  
				delete read; 
				sabundFile.close();
				rabundFile.close();
				listFile.close();
				for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } outputTypes.clear();
				return 0;  
			}
			
			distfile = read->getOutputFile();
		
			list = read->getListVector();
			delete read;
		}else {
			list = new ListVector(globaldata->nameMap->getListVector());
		}
		
		if (m->control_pressed) {  
			sabundFile.close();
			rabundFile.close();
			listFile.close();
			for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } outputTypes.clear();
			return 0;  
		}

		m->mothurOut("It took " + toString(time(NULL) - estart) + " seconds to sort. "); m->mothurOutEndLine();
		estart = time(NULL);
	
		//list vector made by read contains all sequence names
		if(list != NULL){
			rabund = new RAbundVector(list->getRAbundVector());
		}else{
			m->mothurOut("Error: no list vector!"); m->mothurOutEndLine(); return 0;
		}
		
		float previousDist = 0.00000;
		float rndPreviousDist = 0.00000;
		oldRAbund = *rabund;
		oldList = *list;
		
		print_start = true;
		start = time(NULL);
				
		cluster = new HCluster(rabund, list, method, distfile, globaldata->nameMap, cutoff);
		vector<seqDist> seqs; seqs.resize(1); // to start loop
		
		if (m->control_pressed) {  
				delete cluster;
				sabundFile.close();
				rabundFile.close();
				listFile.close();
				for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } outputTypes.clear();
				return 0;  
		}

		
		while (seqs.size() != 0){
		
			seqs = cluster->getSeqs();
			
			if (m->control_pressed) {  
				delete cluster;
				sabundFile.close();
				rabundFile.close();
				listFile.close();
				for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } outputTypes.clear();
				return 0;  
			}

			for (int i = 0; i < seqs.size(); i++) {  //-1 means skip me
				
				if (seqs[i].seq1 != seqs[i].seq2) {
					cluster->update(seqs[i].seq1, seqs[i].seq2, seqs[i].dist);
					
					if (m->control_pressed) {  
						delete cluster;
						sabundFile.close();
						rabundFile.close();
						listFile.close();
						for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } outputTypes.clear();
						return 0;  
					}

			
					float rndDist;
					if (hard) {
						rndDist = m->ceilDist(seqs[i].dist, precision); 
					}else{
						rndDist = m->roundDist(seqs[i].dist, precision); 
					}

					
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

		if (m->control_pressed) {  
			delete cluster;
			sabundFile.close();
			rabundFile.close();
			listFile.close();
			for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } outputTypes.clear();
			return 0;  
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
		
		if (m->control_pressed) {  
			for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } outputTypes.clear();
			return 0;  
		}

		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		m->mothurOut("It took " + toString(time(NULL) - estart) + " seconds to cluster. "); m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "HClusterCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************

void HClusterCommand::printData(string label){
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
		m->errorOut(e, "HClusterCommand", "printData");
		exit(1);
	}


}
//**********************************************************************************************************************

