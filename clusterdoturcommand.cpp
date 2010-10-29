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
vector<string> ClusterDoturCommand::getValidParameters(){	
	try {
		string AlignArray[] =  {"phylip","cutoff","precision","method","outputdir","inputdir"};
		vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterDoturCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
ClusterDoturCommand::ClusterDoturCommand(){	
	try {
		abort = true;
		//initialize outputTypes
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
vector<string> ClusterDoturCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"phylip"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterDoturCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> ClusterDoturCommand::getRequiredFiles(){	
	try {
		vector<string> myArray; 
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterDoturCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************
//This function checks to make sure the cluster command has no errors and then clusters based on the method chosen.
ClusterDoturCommand::ClusterDoturCommand(string option)  {
	try{
		
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"phylip","cutoff","precision","method","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
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
			else if (phylipfile == "not found") { phylipfile = ""; m->mothurOut("When executing the cluster.dotur command you must enter a phylip file."); m->mothurOutEndLine(); abort = true; }	

		
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

void ClusterDoturCommand::help(){
	try {
		m->mothurOut("The cluster.classic command clusters using the algorithm from dotur. \n");
		m->mothurOut("The cluster.classic command parameter options are phylip, name, method, cuttoff, hard, precision. Phylip is required.\n");
		m->mothurOut("The cluster.classic command should be in the following format: \n");
		m->mothurOut("cluster.classic(phylip=yourDistanceMatrix, method=yourMethod, cutoff=yourCutoff, precision=yourPrecision) \n");
		m->mothurOut("The acceptable cluster methods are furthest, nearest, weighted and average.  If no method is provided then furthest is assumed.\n\n");	

	}
	catch(exception& e) {
		m->errorOut(e, "ClusterDoturCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

ClusterDoturCommand::~ClusterDoturCommand(){}

//**********************************************************************************************************************

int ClusterDoturCommand::execute(){
	try {
	
		if (abort == true) {	return 0;	}
		
		if(namefile != ""){	
			nameMap = new NameAssignment(namefile);
			nameMap->readMap();
		}else{
			nameMap = NULL;
		}
		
		//reads phylip file storing data in 2D vector, also fills list and rabund
		ClusterClassic* cluster = new ClusterClassic(cutoff, method);
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

		double saveCutoff = cutoff;
		
		int estart = time(NULL);
	
		while (cluster->getSmallDist() < cutoff && cluster->getNSeqs() > 0){
		
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
	
		if (saveCutoff != cutoff) { 
			if (hard)	{  saveCutoff = m->ceilDist(saveCutoff, precision);	}
			else		{	saveCutoff = m->roundDist(saveCutoff, precision);  }
			m->mothurOut("changed cutoff to " + toString(cutoff)); m->mothurOutEndLine(); 
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
