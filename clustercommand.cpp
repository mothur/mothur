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
#include "clusterdoturcommand.h"


//**********************************************************************************************************************
vector<string> ClusterCommand::setParameters(){	
	try {
		CommandParameter pphylip("phylip", "InputTypes", "", "", "PhylipColumn", "PhylipColumn", "none",false,false); parameters.push_back(pphylip);
		CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "ColumnName",false,false); parameters.push_back(pname);
		CommandParameter pcount("count", "InputTypes", "", "", "NameCount", "none", "none",false,false); parameters.push_back(pcount);
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
		helpString += "The cluster command parameter options are phylip, column, name, count, method, cuttoff, hard, precision, sim, showabund and timing. Phylip or column and name are required, unless you have a valid current file.\n";
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
string ClusterCommand::getOutputFileNameTag(string type, string inputName=""){	
	try {
        string outputFileName = "";
		map<string, vector<string> >::iterator it;
        
        //is this a type this command creates
        it = outputTypes.find(type);
        if (it == outputTypes.end()) {  m->mothurOut("[ERROR]: this command doesn't create a " + type + " output file.\n"); }
        else {
            if (type == "list") {  outputFileName =  "list"; }
            else if (type == "rabund") {  outputFileName =  "rabund"; }
            else if (type == "sabund") {  outputFileName =  "sabund"; }
            else { m->mothurOut("[ERROR]: No definition for type " + type + " output file tag.\n"); m->control_pressed = true;  }
        }
        return outputFileName;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterCommand", "getOutputFileNameTag");
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
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
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
            
            countfile = validParameter.validFile(parameters, "count", true);
			if (countfile == "not open") { abort = true; countfile = ""; }	
			else if (countfile == "not found") { countfile = ""; }
			else { m->setCountTableFile(countfile); }
			
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
				if ((namefile == "") && (countfile == "")){ 
					namefile = m->getNameFile(); 
					if (namefile != "") {  m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
					else { 
						countfile = m->getCountTableFile();
                        if (countfile != "") {  m->mothurOut("Using " + countfile + " as input file for the count parameter."); m->mothurOutEndLine(); }
                        else { 
                            m->mothurOut("You need to provide a namefile or countfile if you are going to use the column format."); m->mothurOutEndLine(); 
                            abort = true; 
                        }	
					}	
				}
			}
			
            if ((countfile != "") && (namefile != "")) { m->mothurOut("When executing a cluster command you must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
            
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			//get user cutoff and precision or use defaults
			string temp;
			temp = validParameter.validFile(parameters, "precision", false);
			if (temp == "not found") { temp = "100"; }
			//saves precision legnth for formatting below
			length = temp.length();
			m->mothurConvert(temp, precision); 
			
			temp = validParameter.validFile(parameters, "hard", false);			if (temp == "not found") { temp = "T"; }
			hard = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "sim", false);				if (temp == "not found") { temp = "F"; }
			sim = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "cutoff", false);
			if (temp == "not found") { temp = "10"; }
			m->mothurConvert(temp, cutoff); 
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
		
		//phylip file given and cutoff not given - use cluster.classic because it uses less memory and is faster
		if ((format == "phylip") && (cutoff > 10.0)) {
			m->mothurOutEndLine(); m->mothurOut("You are using a phylip file and no cutoff.  I will run cluster.classic to save memory and time."); m->mothurOutEndLine();
			
			//run unique.seqs for deconvolute results
			string inputString = "phylip=" + distfile;
			if (namefile != "") { inputString += ", name=" + namefile; }
            else if (countfile != "") { inputString += ", count=" + countfile; }
			inputString += ", precision=" + toString(precision);
			inputString += ", method=" + method;
			if (hard)	{ inputString += ", hard=T";	}
			else		{ inputString += ", hard=F";	}
			if (sim)	{ inputString += ", sim=T";		}
			else		{ inputString += ", sim=F";		}

			
			m->mothurOutEndLine(); 
			m->mothurOut("/------------------------------------------------------------/"); m->mothurOutEndLine(); 
			m->mothurOut("Running command: cluster.classic(" + inputString + ")"); m->mothurOutEndLine(); 
			
			Command* clusterClassicCommand = new ClusterDoturCommand(inputString);
			clusterClassicCommand->execute();
			delete clusterClassicCommand;
			
			m->mothurOut("/------------------------------------------------------------/"); m->mothurOutEndLine();  

			return 0;
		}
		
		ReadMatrix* read;
		if (format == "column") { read = new ReadColumnMatrix(columnfile, sim); }	//sim indicates whether its a similarity matrix
		else if (format == "phylip") { read = new ReadPhylipMatrix(phylipfile, sim); }
		
		read->setCutoff(cutoff);
		
		NameAssignment* nameMap = NULL;
        CountTable* ct = NULL;
		if(namefile != ""){	
			nameMap = new NameAssignment(namefile);
			nameMap->readMap();
            read->read(nameMap);
		}else if (countfile != "") {
            ct = new CountTable();
            ct->readTable(countfile);
            read->read(ct);
        }
		
		list = read->getListVector();
		matrix = read->getDMatrix();
        
		if(countfile != "") {
            rabund = new RAbundVector();
            createRabund(ct, list, rabund); //creates an rabund that includes the counts for the unique list
            delete ct;
        }else { rabund = new RAbundVector(list->getRAbundVector()); }
		delete read;
		
		if (m->control_pressed) { //clean up
			delete list; delete matrix; delete rabund; if(countfile == ""){rabundFile.close(); sabundFile.close();  m->mothurRemove((fileroot+ tag + ".rabund")); m->mothurRemove((fileroot+ tag + ".sabund")); }
			listFile.close(); m->mothurRemove((fileroot+ tag + ".list")); outputTypes.clear(); return 0;
		}
		
		//create cluster
		if (method == "furthest")	{	cluster = new CompleteLinkage(rabund, list, matrix, cutoff, method); }
		else if(method == "nearest"){	cluster = new SingleLinkage(rabund, list, matrix, cutoff, method); }
		else if(method == "average"){	cluster = new AverageLinkage(rabund, list, matrix, cutoff, method);	}
		else if(method == "weighted"){	cluster = new WeightedLinkage(rabund, list, matrix, cutoff, method);	}
		tag = cluster->getTag();
		
		if (outputDir == "") { outputDir += m->hasPath(distfile); }
		fileroot = outputDir + m->getRootName(m->getSimpleName(distfile));
		
        string sabundFileName = fileroot+ tag + "." + getOutputFileNameTag("sabund");
        string rabundFileName = fileroot+ tag + "." + getOutputFileNameTag("rabund");
        string listFileName = fileroot+ tag + ".";
        if (countfile != "") { listFileName += "unique_"; }
        listFileName += getOutputFileNameTag("list");
        
        if (countfile == "") {
            m->openOutputFile(sabundFileName,	sabundFile);
            m->openOutputFile(rabundFileName,	rabundFile);
            outputNames.push_back(sabundFileName); outputTypes["sabund"].push_back(sabundFileName);
            outputNames.push_back(rabundFileName); outputTypes["rabund"].push_back(rabundFileName);

        }
		m->openOutputFile(listFileName,	listFile);
        outputNames.push_back(listFileName); outputTypes["list"].push_back(listFileName);
		
		
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
				if(countfile == "") {rabundFile.close(); sabundFile.close();  m->mothurRemove((fileroot+ tag + ".rabund")); m->mothurRemove((fileroot+ tag + ".sabund")); }
                listFile.close(); m->mothurRemove((fileroot+ tag + ".list")); outputTypes.clear(); return 0;
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
        if (countfile == "") {
            sabundFile.close();
            rabundFile.close();
        }
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
        if (countfile == "") {
            oldRAbund.print(rabundFile);
            oldRAbund.getSAbundVector().print(sabundFile);
        }
       
        if (m->isTrue(showabund)) {
            oldRAbund.getSAbundVector().print(cout);
        }
        
		oldList.setLabel(label);
		oldList.print(listFile);
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterCommand", "printData");
		exit(1);
	}


}
//**********************************************************************************************************************

int ClusterCommand::createRabund(CountTable*& ct, ListVector*& list, RAbundVector*& rabund){
    try {
        rabund->setLabel(list->getLabel());        
        for(int i = 0; i < list->getNumBins(); i++) { 
            if (m->control_pressed) { break; }
            vector<string> binNames;
            string bin = list->get(i);
            m->splitAtComma(bin, binNames);
            int total = 0;
            for (int j = 0; j < binNames.size(); j++) { total += ct->getNumSeqs(binNames[j]);  }
            rabund->push_back(total);   
        }
        return 0;
    }
    catch(exception& e) {
		m->errorOut(e, "ClusterCommand", "createRabund");
		exit(1);
	}
    
}
//**********************************************************************************************************************
