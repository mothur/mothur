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
vector<string> HClusterCommand::setParameters(){	
	try {
		CommandParameter pphylip("phylip", "InputTypes", "", "", "PhylipColumn", "PhylipColumn", "none",false,false); parameters.push_back(pphylip);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "ColumnName",false,false); parameters.push_back(pname);
		CommandParameter pcolumn("column", "InputTypes", "", "", "PhylipColumn", "PhylipColumn", "ColumnName",false,false); parameters.push_back(pcolumn);
		CommandParameter pcutoff("cutoff", "Number", "", "10", "", "", "",false,false); parameters.push_back(pcutoff);
		CommandParameter pprecision("precision", "Number", "", "100", "", "", "",false,false); parameters.push_back(pprecision);
		CommandParameter pmethod("method", "Multiple", "furthest-nearest-average-weighted", "average", "", "", "",false,false); parameters.push_back(pmethod);
		CommandParameter phard("hard", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(phard);
		CommandParameter psorted("sorted", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(psorted);
		CommandParameter pshowabund("showabund", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(pshowabund);
		CommandParameter ptiming("timing", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(ptiming);		
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
			
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "HClusterCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string HClusterCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The hcluster command parameter options are cutoff, precision, method, phylip, column, name, showabund, timing and sorted. Phylip or column and name are required, unless you have valid current files.\n";
		helpString += "The phylip and column parameter allow you to enter your distance file, and sorted indicates whether your column distance file is already sorted. \n";
		helpString += "The name parameter allows you to enter your name file and is required if your distance file is in column format. \n";
		helpString += "The hcluster command should be in the following format: \n";
		helpString += "hcluster(column=youDistanceFile, name=yourNameFile, method=yourMethod, cutoff=yourCutoff, precision=yourPrecision) \n";
		helpString += "The acceptable hcluster methods are furthest, nearest, weighted and average.\n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "HClusterCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
HClusterCommand::HClusterCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
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
//This function checks to make sure the cluster command has no errors and then clusters based on the method chosen.
HClusterCommand::HClusterCommand(string option)  {
	try{
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
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
			else {  distfile = phylipfile;  format = "phylip"; 	m->setPhylipFile(phylipfile); }
			
			columnfile = validParameter.validFile(parameters, "column", true);
			if (columnfile == "not open") { abort = true; }	
			else if (columnfile == "not found") { columnfile = ""; }
			else {  distfile = columnfile; format = "column";	m->setColumnFile(columnfile); }
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			else { m->setNameFile(namefile); }
			
			if ((phylipfile == "") && (columnfile == "")) { 
				//is there are current file available for either of these?
				//give priority to column, then phylip
				columnfile = m->getColumnFile(); 
				if (columnfile != "") {  m->mothurOut("Using " + columnfile + " as input file for the column parameter."); m->mothurOutEndLine(); }
				else { 
					phylipfile = m->getPhylipFile(); 
					if (phylipfile != "") {  m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("No valid current files. You must provide a phylip or column file before you can use the hcluster command."); m->mothurOutEndLine(); 
						abort = true;
					}
				}
			}
			else if ((phylipfile != "") && (columnfile != "")) { m->mothurOut("When executing a hcluster command you must enter ONLY ONE of the following: phylip or column."); m->mothurOutEndLine(); abort = true; }
		
			if (columnfile != "") {
				if (namefile == "") { 
					namefile = m->getNameFile(); 
					if (namefile != "") {  m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("You need to provide a namefile if you are going to use the column format."); m->mothurOutEndLine(); 
						abort = true; 
					}	
				}
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
			
			temp = validParameter.validFile(parameters, "hard", false);			if (temp == "not found") { temp = "T"; }
			hard = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "cutoff", false);
			if (temp == "not found") { temp = "10"; }
			convert(temp, cutoff); 
			cutoff += (5 / (precision * 10.0)); 
			
			method = validParameter.validFile(parameters, "method", false);
			if (method == "not found") { method = "average"; }
			
			if ((method == "furthest") || (method == "nearest") || (method == "average") || (method == "weighted")) { }
			else { m->mothurOut("Not a valid clustering method.  Valid clustering algorithms are furthest, nearest, average or weighted."); m->mothurOutEndLine(); abort = true; }

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
				else if (method == "weighted")	{ tag = "wn";  }
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

int HClusterCommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		NameAssignment* nameMap = NULL;
		if(namefile != ""){	
			nameMap = new NameAssignment(namefile);
			nameMap->readMap();
		}		
		
		time_t estart = time(NULL);
		
		if (!sorted) {
			read = new ReadCluster(distfile, cutoff, outputDir, true); 	
			read->setFormat(format);
			read->read(nameMap);
			
			if (m->control_pressed) {  
				delete read; 
				sabundFile.close();
				rabundFile.close();
				listFile.close();
				for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } outputTypes.clear();
				return 0;  
			}
			
			distfile = read->getOutputFile();
		
			list = read->getListVector();
			delete read;
		}else {
			list = new ListVector(nameMap->getListVector());
		}
		
		if (m->control_pressed) {  
			sabundFile.close();
			rabundFile.close();
			listFile.close();
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } outputTypes.clear();
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
				
		cluster = new HCluster(rabund, list, method, distfile, nameMap, cutoff);
		vector<seqDist> seqs; seqs.resize(1); // to start loop
		
		if (m->control_pressed) {  
				delete cluster;
				sabundFile.close();
				rabundFile.close();
				listFile.close();
				for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } outputTypes.clear();
				return 0;  
		}

		float saveCutoff = cutoff;
		
		while (seqs.size() != 0){
		
			seqs = cluster->getSeqs();
			
			//to account for cutoff change in average neighbor
			if (seqs.size() != 0) {
				if (seqs[0].dist > cutoff) { break; }
			}
			
			if (m->control_pressed) {  
				delete cluster;
				sabundFile.close();
				rabundFile.close();
				listFile.close();
				for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } outputTypes.clear();
				return 0;  
			}

			for (int i = 0; i < seqs.size(); i++) {  //-1 means skip me
				
				if (seqs[i].seq1 != seqs[i].seq2) {
					cutoff = cluster->update(seqs[i].seq1, seqs[i].seq2, seqs[i].dist);
					
					if (m->control_pressed) {  
						delete cluster;
						sabundFile.close();
						rabundFile.close();
						listFile.close();
						for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } outputTypes.clear();
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
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } outputTypes.clear();
			return 0;  
		}
					
		if(previousDist <= 0.0000){
			printData("unique");
		}
		else if(rndPreviousDist<cutoff){
			printData(toString(rndPreviousDist, length-1));
		}
				
		sabundFile.close();
		rabundFile.close();
		listFile.close();
		delete cluster;
		
		if (m->control_pressed) {  
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } outputTypes.clear();
			return 0;  
		}
		
		
		if (saveCutoff != cutoff) { 
			if (hard)	{  saveCutoff = m->ceilDist(saveCutoff, precision);	}
			else		{	saveCutoff = m->roundDist(saveCutoff, precision);  }
			
			m->mothurOut("changed cutoff to " + toString(cutoff)); m->mothurOutEndLine(); 
		}
		
		//set list file as new current listfile
		string current = "";
		itTypes = outputTypes.find("list");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setListFile(current); }
		}
		
		//set rabund file as new current rabundfile
		itTypes = outputTypes.find("rabund");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setRabundFile(current); }
		}
		
		//set sabund file as new current sabundfile
		itTypes = outputTypes.find("sabund");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setSabundFile(current); }
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

