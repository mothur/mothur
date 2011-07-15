/*
 *  clustercommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "clustercommand.h"
#include "readphylip.h"
#include "readcolumn.h"
#include "readmatrix.hpp"

//**********************************************************************************************************************
vector<string> ClusterCommand::setParameters(){	
	try {
		CommandParameter pphylip("phylip", "InputTypes", "", "", "PhylipColumn", "PhylipColumn", "none",false,false); parameters.push_back(pphylip);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "ColumnName",false,false); parameters.push_back(pname);
		CommandParameter pcolumn("column", "InputTypes", "", "", "PhylipColumn", "PhylipColumn", "ColumnName",false,false); parameters.push_back(pcolumn);		
		CommandParameter pcutoff("cutoff", "Number", "", "10", "", "", "",false,false); parameters.push_back(pcutoff);
		CommandParameter pprecision("precision", "Number", "", "100", "", "", "",false,false); parameters.push_back(pprecision);
		CommandParameter pmethod("method", "Multiple", "furthest-nearest-average-weighted", "average", "", "", "",false,false); parameters.push_back(pmethod);
		CommandParameter pshowabund("showabund", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(pshowabund);
		CommandParameter ptiming("timing", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(ptiming);
		CommandParameter psim("sim", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(psim);
		CommandParameter phard("hard", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(phard);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ClusterCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The cluster command parameter options are phylip, column, name, method, cuttoff, hard, precision, sim, showabund and timing. Phylip or column and name are required, unless you have a valid current file.\n";
		helpString += "The cluster command should be in the following format: \n";
		helpString += "cluster(method=yourMethod, cutoff=yourCutoff, precision=yourPrecision) \n";
		helpString += "The acceptable cluster methods are furthest, nearest, average and weighted.  If no method is provided then average is assumed.\n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
ClusterCommand::ClusterCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["list"] = tempOutNames;
		outputTypes["rabund"] = tempOutNames;
		outputTypes["sabund"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterCommand", "ClusterCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
//This function checks to make sure the cluster command has no errors and then clusters based on the method chosen.
ClusterCommand::ClusterCommand(string option)  {
	try{
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			map<string,string>::iterator it;
			
			ValidParameters validParameter;
		
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
		
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
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
			
			//check for required parameters
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not open") { phylipfile = ""; abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; }	
			else {  distfile = phylipfile;  format = "phylip"; 	m->setPhylipFile(phylipfile); }
			
			columnfile = validParameter.validFile(parameters, "column", true);
			if (columnfile == "not open") { columnfile = ""; abort = true; }	
			else if (columnfile == "not found") { columnfile = ""; }
			else {  distfile = columnfile; format = "column"; m->setColumnFile(columnfile);	}
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			else { m->setNameFile(namefile); }
			
			if ((phylipfile == "") && (columnfile == "")) { 
				//is there are current file available for either of these?
				//give priority to column, then phylip
				columnfile = m->getColumnFile(); 
				if (columnfile != "") {  distfile = columnfile; format = "column"; m->mothurOut("Using " + columnfile + " as input file for the column parameter."); m->mothurOutEndLine(); }
				else { 
					phylipfile = m->getPhylipFile(); 
					if (phylipfile != "") { distfile = phylipfile;  format = "phylip"; m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("No valid current files. You must provide a phylip or column file before you can use the cluster command."); m->mothurOutEndLine(); 
						abort = true;
					}
				}
			}
			else if ((phylipfile != "") && (columnfile != "")) { m->mothurOut("When executing a cluster command you must enter ONLY ONE of the following: phylip or column."); m->mothurOutEndLine(); abort = true; }
			
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
			
			temp = validParameter.validFile(parameters, "sim", false);				if (temp == "not found") { temp = "F"; }
			sim = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "cutoff", false);
			if (temp == "not found") { temp = "10"; }
			convert(temp, cutoff); 
			cutoff += (5 / (precision * 10.0));  
			
			method = validParameter.validFile(parameters, "method", false);
			if (method == "not found") { method = "average"; }
			
			if ((method == "furthest") || (method == "nearest") || (method == "average") || (method == "weighted")) { }
			else { m->mothurOut("Not a valid clustering method.  Valid clustering algorithms are furthest, nearest, average, and weighted."); m->mothurOutEndLine(); abort = true; }

			showabund = validParameter.validFile(parameters, "showabund", false);
			if (showabund == "not found") { showabund = "T"; }

			timing = validParameter.validFile(parameters, "timing", false);
			if (timing == "not found") { timing = "F"; }
			
			}
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterCommand", "ClusterCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
ClusterCommand::~ClusterCommand(){}
//**********************************************************************************************************************

int ClusterCommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		ReadMatrix* read;
		if (format == "column") { read = new ReadColumnMatrix(columnfile, sim); }	//sim indicates whether its a similarity matrix
		else if (format == "phylip") { read = new ReadPhylipMatrix(phylipfile, sim); }
		
		read->setCutoff(cutoff);
		
		NameAssignment* nameMap = NULL;
		if(namefile != ""){	
			nameMap = new NameAssignment(namefile);
			nameMap->readMap();
		}
		
		read->read(nameMap);
		list = read->getListVector();
		matrix = read->getMatrix();
		rabund = new RAbundVector(list->getRAbundVector());
		delete read;
		
		if (m->control_pressed) { //clean up
			delete list; delete matrix; delete rabund;
			sabundFile.close();rabundFile.close();listFile.close();
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} outputTypes.clear();
			return 0;
		}
		
		//create cluster
		if (method == "furthest")	{	cluster = new CompleteLinkage(rabund, list, matrix, cutoff, method); }
		else if(method == "nearest"){	cluster = new SingleLinkage(rabund, list, matrix, cutoff, method); }
		else if(method == "average"){	cluster = new AverageLinkage(rabund, list, matrix, cutoff, method);	}
		else if(method == "weighted"){	cluster = new WeightedLinkage(rabund, list, matrix, cutoff, method);	}
		tag = cluster->getTag();
		
		if (outputDir == "") { outputDir += m->hasPath(distfile); }
		fileroot = outputDir + m->getRootName(m->getSimpleName(distfile));
		
		m->openOutputFile(fileroot+ tag + ".sabund",	sabundFile);
		m->openOutputFile(fileroot+ tag + ".rabund",	rabundFile);
		m->openOutputFile(fileroot+ tag + ".list",		listFile);
		
		outputNames.push_back(fileroot+ tag + ".sabund"); outputTypes["sabund"].push_back(fileroot+ tag + ".sabund");
		outputNames.push_back(fileroot+ tag + ".rabund"); outputTypes["rabund"].push_back(fileroot+ tag + ".rabund");
		outputNames.push_back(fileroot+ tag + ".list"); outputTypes["list"].push_back(fileroot+ tag + ".list");
		
		
		time_t estart = time(NULL);
		float previousDist = 0.00000;
		float rndPreviousDist = 0.00000;
		oldRAbund = *rabund;
		oldList = *list;

		print_start = true;
		start = time(NULL);
		loops = 0;
		double saveCutoff = cutoff;
		
		while (matrix->getSmallDist() < cutoff && matrix->getNNodes() > 0){
		
			if (m->control_pressed) { //clean up
				delete list; delete matrix; delete rabund; delete cluster;
				sabundFile.close();rabundFile.close();listFile.close();
				for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} outputTypes.clear();
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
		
		delete matrix;
		delete list;
		delete rabund;
		delete cluster;
		
		sabundFile.close();
		rabundFile.close();
		listFile.close();
	
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
