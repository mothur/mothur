/*
 *  clusterdoturcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 10/27/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "clusterdoturcommand.h"
#include "clusterclassic.h"

//**********************************************************************************************************************
vector<string> ClusterDoturCommand::setParameters(){	
	try {
		CommandParameter pphylip("phylip", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pphylip);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pname);
		CommandParameter pcutoff("cutoff", "Number", "", "10", "", "", "",false,false); parameters.push_back(pcutoff);
		CommandParameter pprecision("precision", "Number", "", "100", "", "", "",false,false); parameters.push_back(pprecision);
		CommandParameter pmethod("method", "Multiple", "furthest-nearest-average-weighted", "furthest", "", "", "",false,false); parameters.push_back(pmethod);
		CommandParameter phard("hard", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(phard);
		CommandParameter psim("sim", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(psim);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterDoturCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ClusterDoturCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The cluster.classic command clusters using the algorithm from dotur. \n";
		helpString += "The cluster.classic command parameter options are phylip, name, method, cuttoff, hard, sim, precision. Phylip is required, unless you have a valid current file.\n";
		helpString += "The cluster.classic command should be in the following format: \n";
		helpString += "cluster.classic(phylip=yourDistanceMatrix, method=yourMethod, cutoff=yourCutoff, precision=yourPrecision) \n";
		helpString += "The acceptable cluster methods are furthest, nearest, weighted and average.  If no method is provided then furthest is assumed.\n\n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterDoturCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
ClusterDoturCommand::ClusterDoturCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
		vector<string> tempOutNames;
		outputTypes["list"] = tempOutNames;
		outputTypes["rabund"] = tempOutNames;
		outputTypes["sabund"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterDoturCommand", "ClusterCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
//This function checks to make sure the cluster command has no errors and then clusters based on the method chosen.
ClusterDoturCommand::ClusterDoturCommand(string option)  {
	try{
		
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			map<string,string>::iterator it;
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {
					abort = true;
				}
			}
			
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
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}

			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["list"] = tempOutNames;
			outputTypes["rabund"] = tempOutNames;
			outputTypes["sabund"] = tempOutNames;
		
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			//check for required parameters
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not open") { abort = true; }
			else if (phylipfile == "not found") { 
				phylipfile = m->getPhylipFile(); 
				if (phylipfile != "") {  m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter."); m->mothurOutEndLine(); }
				else { 
					m->mothurOut("You need to provide a phylip file with the cluster.classic command."); m->mothurOutEndLine(); 
					abort = true; 
				}	
			}	

		
			//check for optional parameter and set defaults
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			
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
			
			temp = validParameter.validFile(parameters, "hard", false);			if (temp == "not found") { temp = "F"; }
			hard = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "sim", false);				if (temp == "not found") { temp = "F"; }
			sim = m->isTrue(temp); 
			
			method = validParameter.validFile(parameters, "method", false);
			if (method == "not found") { method = "furthest"; }
			
			if ((method == "furthest") || (method == "nearest") || (method == "average") || (method == "weighted")) { 
				if (method == "furthest") { tag = "fn"; }
				else if (method == "nearest") { tag = "nn"; }
				else if (method == "average") { tag = "an"; }
				else if (method == "weighted") { tag = "wn"; }
			}else { m->mothurOut("Not a valid clustering method.  Valid clustering algorithms are furthest, nearest, average, weighted."); m->mothurOutEndLine(); abort = true; }
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterDoturCommand", "ClusterCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int ClusterDoturCommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		if(namefile != ""){	
			nameMap = new NameAssignment(namefile);
			nameMap->readMap();
		}else{
			nameMap = NULL;
		}
		
		//reads phylip file storing data in 2D vector, also fills list and rabund
		ClusterClassic* cluster = new ClusterClassic(cutoff, method, sim);
		cluster->readPhylipFile(phylipfile, nameMap);
		
		if (m->control_pressed) { delete cluster; delete list; delete rabund; return 0; }
		
		list = cluster->getListVector();
		rabund = cluster->getRAbundVector();
						
		if (outputDir == "") { outputDir += m->hasPath(phylipfile); }
		fileroot = outputDir + m->getRootName(m->getSimpleName(phylipfile));
			
		m->openOutputFile(fileroot+ tag + ".sabund",	sabundFile);
		m->openOutputFile(fileroot+ tag + ".rabund",	rabundFile);
		m->openOutputFile(fileroot+ tag + ".list",		listFile);
				
		outputNames.push_back(fileroot+ tag + ".sabund"); outputTypes["sabund"].push_back(fileroot+ tag + ".sabund");
		outputNames.push_back(fileroot+ tag + ".rabund"); outputTypes["rabund"].push_back(fileroot+ tag + ".rabund");
		outputNames.push_back(fileroot+ tag + ".list"); outputTypes["list"].push_back(fileroot+ tag + ".list");
		
		float previousDist = 0.00000;
		float rndPreviousDist = 0.00000;
		oldRAbund = *rabund;
		oldList = *list;

		//double saveCutoff = cutoff;
		
		int estart = time(NULL);
	
		while ((cluster->getSmallDist() < cutoff) && (cluster->getNSeqs() > 1)){
			if (m->control_pressed) { delete cluster; delete list; delete rabund; sabundFile.close();rabundFile.close();listFile.close();  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); 	} outputTypes.clear();  return 0;  }
		
			cluster->update(cutoff);
	
			float dist = cluster->getSmallDist();
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
	
		if(previousDist <= 0.0000){
			printData("unique");
		}
		else if(rndPreviousDist<cutoff){
			printData(toString(rndPreviousDist, length-1));
		}
					
		sabundFile.close();
		rabundFile.close();
		listFile.close();
		
		delete cluster; delete nameMap; delete list; delete rabund;
	
		//if (saveCutoff != cutoff) { 
		//	if (hard)	{  saveCutoff = m->ceilDist(saveCutoff, precision);	}
		//	else		{	saveCutoff = m->roundDist(saveCutoff, precision);  }
		//	m->mothurOut("changed cutoff to " + toString(cutoff)); m->mothurOutEndLine(); 
		//}
		
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

		m->mothurOut("It took " + toString(time(NULL) - estart) + " seconds to cluster"); m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterDoturCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************

void ClusterDoturCommand::printData(string label){
	try {
	
		oldRAbund.setLabel(label);
		oldRAbund.print(rabundFile);
		oldRAbund.getSAbundVector().print(sabundFile);
		
		oldRAbund.getSAbundVector().print(cout);
		
		oldList.setLabel(label);
		oldList.print(listFile);
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterDoturCommand", "printData");
		exit(1);
	}
}
//**********************************************************************************************************************
