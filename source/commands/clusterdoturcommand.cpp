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
		CommandParameter pphylip("phylip", "InputTypes", "", "", "none", "none", "none","list",false,true,true); parameters.push_back(pphylip);
		CommandParameter pname("name", "InputTypes", "", "", "namecount", "none", "none","rabund-sabund",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "namecount", "none", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter pcutoff("cutoff", "Number", "", "10", "", "", "","",false,false,true); parameters.push_back(pcutoff);
		CommandParameter pprecision("precision", "Number", "", "100", "", "", "","",false,false); parameters.push_back(pprecision);
		CommandParameter pmethod("method", "Multiple", "furthest-nearest-average-weighted", "average", "", "", "","",false,false); parameters.push_back(pmethod);
		CommandParameter phard("hard", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(phard);
		CommandParameter psim("sim", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(psim);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
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
		helpString += "The cluster.classic command parameter options are phylip, name, count, method, cuttoff, hard, sim, precision. Phylip is required, unless you have a valid current file.\n";
		helpString += "The cluster.classic command should be in the following format: \n";
		helpString += "cluster.classic(phylip=yourDistanceMatrix, method=yourMethod, cutoff=yourCutoff, precision=yourPrecision) \n";
		helpString += "The acceptable cluster methods are furthest, nearest, weighted and average.  If no method is provided then average is assumed.\n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterDoturCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ClusterDoturCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "list") {  pattern = "[filename],[clustertag],list-[filename],[clustertag],[tag2],list"; } 
        else if (type == "rabund") {  pattern = "[filename],[clustertag],rabund"; } 
        else if (type == "sabund") {  pattern = "[filename],[clustertag],sabund"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterDoturCommand", "getOutputPattern");
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
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
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
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
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
			}else { m->setPhylipFile(phylipfile); }	

		
			//check for optional parameter and set defaults
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; namefile = ""; }	
			else if (namefile == "not found") { namefile = ""; }
			else { m->setNameFile(namefile); }
            
            countfile = validParameter.validFile(parameters, "count", true);
			if (countfile == "not open") { abort = true; countfile = ""; }	
			else if (countfile == "not found") { countfile = ""; }
			else { m->setCountTableFile(countfile); }
			
            if ((countfile != "") && (namefile != "")) { m->mothurOut("When executing a cluster.classic command you must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
            
			string temp;
			temp = validParameter.validFile(parameters, "precision", false);
			if (temp == "not found") { temp = "100"; }
			//saves precision legnth for formatting below
			length = temp.length();
			m->mothurConvert(temp, precision); 
			
			temp = validParameter.validFile(parameters, "cutoff", false);
			if (temp == "not found") { temp = "10"; }
			m->mothurConvert(temp, cutoff); 
			cutoff += (5 / (precision * 10.0));  
			
			temp = validParameter.validFile(parameters, "hard", false);			if (temp == "not found") { temp = "T"; }
			hard = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "sim", false);				if (temp == "not found") { temp = "F"; }
			sim = m->isTrue(temp); 
			
			method = validParameter.validFile(parameters, "method", false);
			if (method == "not found") { method = "average"; }
			
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
		
        
        ClusterClassic* cluster = new ClusterClassic(cutoff, method, sim);
        
        NameAssignment* nameMap = NULL;
        CountTable* ct = NULL;
        map<string, int> counts;
        if(namefile != "") {	
			nameMap = new NameAssignment(namefile);
			nameMap->readMap();
            cluster->readPhylipFile(phylipfile, nameMap);
            delete nameMap;
		}else if (countfile != "") {
            ct = new CountTable();
            ct->readTable(countfile, false, false);
            cluster->readPhylipFile(phylipfile, ct);
            counts = ct->getNameMap();
            delete ct;
        }else {
            cluster->readPhylipFile(phylipfile, nameMap);
        }
        tag = cluster->getTag();
        
		if (m->control_pressed) { delete cluster; return 0; }
		
		list = cluster->getListVector();
		rabund = cluster->getRAbundVector();
								
		if (outputDir == "") { outputDir += m->hasPath(phylipfile); }
		fileroot = outputDir + m->getRootName(m->getSimpleName(phylipfile));
			
        map<string, string> variables; 
        variables["[filename]"] = fileroot;
        variables["[clustertag]"] = tag;
        string sabundFileName = getOutputFileName("sabund", variables);
        string rabundFileName = getOutputFileName("rabund", variables);
        if (countfile != "") { variables["[tag2]"] = "unique_list"; }
        string listFileName = getOutputFileName("list", variables);
        
        if (countfile == "") {
            m->openOutputFile(sabundFileName,	sabundFile);
            m->openOutputFile(rabundFileName,	rabundFile);
            outputNames.push_back(sabundFileName); outputTypes["sabund"].push_back(sabundFileName);
            outputNames.push_back(rabundFileName); outputTypes["rabund"].push_back(rabundFileName);
            
        }
		m->openOutputFile(listFileName,	listFile);
        outputNames.push_back(listFileName); outputTypes["list"].push_back(listFileName);
        list->printHeaders(listFile);
		
		float previousDist = 0.00000;
		float rndPreviousDist = 0.00000;
		oldRAbund = *rabund;
		oldList = *list;

		//double saveCutoff = cutoff;
		
		int estart = time(NULL);
	
		while ((cluster->getSmallDist() < cutoff) && (cluster->getNSeqs() > 1)){
			if (m->control_pressed) { delete cluster; delete list; delete rabund; if(countfile == "") {rabundFile.close(); sabundFile.close();  m->mothurRemove((fileroot+ tag + ".rabund")); m->mothurRemove((fileroot+ tag + ".sabund")); }
                listFile.close(); m->mothurRemove((fileroot+ tag + ".list")); outputTypes.clear();  return 0;  }
		
			cluster->update(cutoff);
	
			float dist = cluster->getSmallDist();
			float rndDist;
			if (hard) {
				rndDist = m->ceilDist(dist, precision); 
			}else{
				rndDist = m->roundDist(dist, precision); 
			}

			if(previousDist <= 0.0000 && dist != previousDist){
				printData("unique", counts);
			}
			else if(rndDist != rndPreviousDist){
				printData(toString(rndPreviousDist,  length-1), counts);
			}
		
			previousDist = dist;
			rndPreviousDist = rndDist;
			oldRAbund = *rabund;
			oldList = *list;
		}
	
		if(previousDist <= 0.0000){
			printData("unique", counts);
		}
		else if(rndPreviousDist<cutoff){
			printData(toString(rndPreviousDist, length-1), counts);
		}
		
        if (countfile == "") {
            sabundFile.close();
            rabundFile.close();
        }
		listFile.close();
		
		delete cluster;  delete list; delete rabund;
		
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

void ClusterDoturCommand::printData(string label, map<string, int>& counts){
	try {
        oldRAbund.setLabel(label);
        if (countfile == "") {
            oldRAbund.print(rabundFile);
            oldRAbund.getSAbundVector().print(sabundFile);
        }

		oldRAbund.getSAbundVector().print(cout);
		
		oldList.setLabel(label);
        
        if(countfile != "") {
            oldList.print(listFile, counts);
        }else {
            oldList.print(listFile);
        }

	}
	catch(exception& e) {
		m->errorOut(e, "ClusterDoturCommand", "printData");
		exit(1);
	}
}
//**********************************************************************************************************************
