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

#include "sequence.hpp"
#include "systemcommand.h"

//**********************************************************************************************************************
vector<string> ClusterCommand::setParameters(){	
	try {
        CommandParameter pphylip("phylip", "InputTypes", "", "", "PhylipColumnFasta", "PhylipColumnFasta", "none","list",false,false,true); parameters.push_back(pphylip);
        CommandParameter pfasta("fasta", "InputTypes", "", "", "PhylipColumnFasta", "PhylipColumnFasta", "FastaTaxName","list",false,false,true); parameters.push_back(pfasta);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "ColumnName-FastaTaxName","rabund-sabund",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount", "none", "","",false,false,true); parameters.push_back(pcount);
        CommandParameter pcolumn("column", "InputTypes", "", "", "PhylipColumnFasta", "PhylipColumnFasta", "ColumnName","list",false,false,true); parameters.push_back(pcolumn);
		CommandParameter pcutoff("cutoff", "Number", "", "10", "", "", "","",false,false,true); parameters.push_back(pcutoff);
		CommandParameter pprecision("precision", "Number", "", "100", "", "", "","",false,false); parameters.push_back(pprecision);
		CommandParameter pmethod("method", "Multiple", "furthest-nearest-average-weighted-agc-dgc-opti", "average", "", "", "","",false,false,true); parameters.push_back(pmethod);
        CommandParameter pmetric("metric", "Multiple", "mcc-sens-spec-tptn-fpfn-tp2tn", "mcc", "", "", "","",false,false,true); parameters.push_back(pmetric);
        CommandParameter pmetriccutoff("delta", "Number", "", "0.001", "", "", "","",false,false,true); parameters.push_back(pmetriccutoff);
        CommandParameter piters("iters", "Number", "", "10000", "", "", "","",false,false,true); parameters.push_back(piters);
		CommandParameter pshowabund("showabund", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pshowabund);
		CommandParameter ptiming("timing", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(ptiming);
		CommandParameter psim("sim", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(psim);
		CommandParameter phard("hard", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(phard);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        //CommandParameter padjust("adjust", "String", "", "F", "", "", "","",false,false); parameters.push_back(padjust);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
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
		helpString += "The cluster command parameter options are phylip, column, name, count, method, cuttoff, hard, precision, sim, showabund, metric, delta, iters and timing. Fasta or Phylip or column and name are required.\n";
		//helpString += "The adjust parameter is used to handle missing distances.  If you set a cutoff, adjust=f by default.  If not, adjust=t by default. Adjust=f, means ignore missing distances and adjust cutoff as needed with the average neighbor method.  Adjust=t, will treat missing distances as 1.0. You can also set the value the missing distances should be set to, adjust=0.5 would give missing distances a value of 0.5.\n";
        helpString += "The phylip and column parameter allow you to enter your distance file. \n";
        helpString += "The fasta parameter allows you to enter your fasta file for use with the agc or dgc methods. \n";
        helpString += "The name parameter allows you to enter your name file. \n";
        helpString += "The count parameter allows you to enter your count file. \n A count or name file is required if your distance file is in column format.\n";
        helpString += "The iters parameter allow you to set the maxiters for the opticluster method. \n";
        helpString += "The metric parameter allows to select the metric in the opticluster method. Options are Matthews correlation coefficient (mcc), sensitivity (sens), specificity (spec), true positives + true negatives (tptn), false positives + false negatives (fpfn), true positives + 2* true negatives (tp2tn). Default=mcc.\n";
        helpString += "The delta parameter allows to set the stable value for the metric in the opticluster method. \n";
        helpString += "The method parameter allows you to enter your clustering mothod. Options are furthest, nearest, average, weighted, agc, dgc and opti. Default=average.  The agc and dgc methods require a fasta file.";
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
string ClusterCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "list") {  pattern = "[filename],[clustertag],list-[filename],[clustertag],[tag2],list"; } 
        else if (type == "rabund") {  pattern = "[filename],[clustertag],rabund"; } 
        else if (type == "sabund") {  pattern = "[filename],[clustertag],sabund"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterCommand", "getOutputPattern");
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
			
			inputDir = validParameter.validFile(parameters, "inputdir", false);		
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
                
                it = parameters.find("fasta");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
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
			
            fastafile = validParameter.validFile(parameters, "fasta", true);
            if (fastafile == "not open") { abort = true; }
            else if (fastafile == "not found") { fastafile = ""; }
            else { distfile = fastafile;  format = "fasta"; m->setFastaFile(fastafile); }
            
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			else { m->setNameFile(namefile); }
            
            countfile = validParameter.validFile(parameters, "count", true);
			if (countfile == "not open") { abort = true; countfile = ""; }	
			else if (countfile == "not found") { countfile = ""; }
			else { m->setCountTableFile(countfile); }
			
			if ((phylipfile == "") && (columnfile == "") && (fastafile == "")) {
				//is there are current file available for either of these?
				//give priority to column, then phylip
				columnfile = m->getColumnFile(); 
				if (columnfile != "") {  distfile = columnfile; format = "column"; m->mothurOut("Using " + columnfile + " as input file for the column parameter."); m->mothurOutEndLine(); }
				else { 
					phylipfile = m->getPhylipFile(); 
					if (phylipfile != "") { distfile = phylipfile;  format = "phylip"; m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter."); m->mothurOutEndLine(); }
					else {
                        fastafile = m->getFastaFile();
                        if (fastafile != "") {  distfile = fastafile; format = "fasta"; m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
                        else {
                            m->mothurOut("No valid current files. You must provide a phylip, column or fasta file before you can use the cluster command."); m->mothurOutEndLine();
                            abort = true;
                        }
					}
				}
			}
			else if (((phylipfile != "") && (columnfile != "")) || ((phylipfile != "") && (fastafile != "")) || ((fastafile != "") && (columnfile != "")))  { m->mothurOut("When executing a cluster command you must enter ONLY ONE of the following: phylip, column or fasta."); m->mothurOutEndLine(); abort = true; }
			
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
			
            temp = validParameter.validFile(parameters, "delta", false);		if (temp == "not found")  { temp = "0.001"; }
            m->mothurConvert(temp, stableMetric);
            
            metric = validParameter.validFile(parameters, "metric", false);		if (metric == "not found") { metric = "mcc"; }
            
            if ((metric == "mcc") || (metric == "sens") || (metric == "spec") || (metric == "tptn") || (metric == "tp2tn") || (metric == "fpfn")){ }
            else { m->mothurOut("[ERROR]: Not a valid metric.  Valid metrics are mcc."); m->mothurOutEndLine(); abort = true; }

            temp = validParameter.validFile(parameters, "iters", false);		if (temp == "not found")  { temp = "1000"; }
            m->mothurConvert(temp, maxIters);
            
            //bool cutoffSet = false;
			temp = validParameter.validFile(parameters, "cutoff", false);
            if (temp == "not found") { temp = "1.0"; cutoffNotSet = true; }
            //else { cutoffSet = true; }
			m->mothurConvert(temp, cutoff); 
			cutoff += (5 / (precision * 10.0));
            
            //temp = validParameter.validFile(parameters, "adjust", false);				if (temp == "not found") { temp = "F"; }
            //if (m->isNumeric1(temp))    { m->mothurConvert(temp, adjust);   }
            //else if (m->isTrue(temp))   { adjust = 1.0;                     }
            //else                        { adjust = -1.0;                    }
            adjust=-1.0;
			
			method = validParameter.validFile(parameters, "method", false);
			if (method == "not found") { method = "average"; }
			
            if ((method == "furthest") || (method == "nearest") || (method == "average") || (method == "weighted") || (method == "agc") || (method == "dgc") || (method == "opti")) { }
            else { m->mothurOut("[ERROR]: Not a valid clustering method.  Valid clustering algorithms are furthest, nearest, average, weighted, agc, dgc and opti."); m->mothurOutEndLine(); abort = true; }
        
            
            if ((method == "agc") || (method == "dgc")) {
                if (fastafile == "") { m->mothurOut("[ERROR]: You must provide a fasta file when using the agc or dgc clustering methods, aborting\n."); abort = true;}
            }
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
            if ((method == "agc") || (method == "dgc")) { m->mothurOut("[ERROR]: The agc and dgc clustering methods are not available for Windows, aborting\n."); abort = true; }
#endif
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
		
        time_t estart = time(NULL);
        
        if (format == "fasta")      {   runVsearchCluster();    }
        else if (method == "opti")  {   runOptiCluster();       }
        else                        {   runMothurCluster();     }
        
		if (m->control_pressed) { 	for (int j = 0; j < outputNames.size(); j++) { m->mothurRemove(outputNames[j]); }  return 0; }
        
        m->mothurOut("It took " + toString(time(NULL) - estart) + " seconds to cluster"); m->mothurOutEndLine();
        
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

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int ClusterCommand::runVsearchCluster(){
    try {
        //look for uchime exe
        string path = m->argv;
        string tempPath = path;
        for (int i = 0; i < path.length(); i++) { tempPath[i] = tolower(path[i]); }
        path = path.substr(0, (tempPath.find_last_of('m')));
        
        string vsearchCommand;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        vsearchCommand = path + "vsearch";	//	format the database, -o option gives us the ability
        if (m->debug) {
            m->mothurOut("[DEBUG]: vsearch location using \"which vsearch\" = ");
            Command* newCommand = new SystemCommand("which vsearch"); m->mothurOutEndLine();
            newCommand->execute();
            delete newCommand;
            m->mothurOut("[DEBUG]: Mothur's location using \"which mothur\" = ");
            newCommand = new SystemCommand("which mothur"); m->mothurOutEndLine();
            newCommand->execute();
            delete newCommand;
        }
#else
        vsearchCommand = path + "vsearch.exe";
#endif
        
        //test to make sure uchime exists
        ifstream in;
        vsearchCommand = m->getFullPathName(vsearchCommand);
        int ableToOpen = m->openInputFile(vsearchCommand, in, "no error"); in.close();
        if(ableToOpen == 1) {
            m->mothurOut(vsearchCommand + " file does not exist. Checking path... \n");
            //check to see if uchime is in the path??
            
            string uLocation = m->findProgramPath("vsearch");
            
            
            ifstream in2;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
            ableToOpen = m->openInputFile(uLocation, in2, "no error"); in2.close();
#else
            ableToOpen = m->openInputFile((uLocation + ".exe"), in2, "no error"); in2.close();
#endif
            
            if(ableToOpen == 1) { m->mothurOut("[ERROR]: " + uLocation + " file does not exist. mothur requires the vsearch executable."); m->mothurOutEndLine(); abort = true; }
            else {  m->mothurOut("Found vsearch in your path, using " + uLocation + "\n");vsearchLocation = uLocation; }
        }else {  vsearchLocation = vsearchCommand; }
        
        vsearchLocation = m->getFullPathName(vsearchLocation);
        
        string vsearchFastafile = ""; VsearchFileParser* vParse;
        if ((namefile == "") && (countfile == ""))  { vParse = new VsearchFileParser(fastafile);                        }
        else if (namefile != "")                    { vParse = new VsearchFileParser(fastafile, namefile, "name");      }
        else if (countfile != "")                   { vParse = new VsearchFileParser(fastafile, countfile, "count");    }
        else                                        { m->mothurOut("[ERROR]: Opps, should never get here. ClusterCommand::runVsearchCluster() \n"); m->control_pressed = true; }
    
        if (m->control_pressed) {  return 0; }
        
        vsearchFastafile = vParse->getVsearchFile();
        
        if (cutoff > 1.0) {  m->mothurOut("You did not set a cutoff, using 0.03.\n"); cutoff = 0.03; }
        
        //Run vsearch
        string ucVsearchFile = m->getSimpleName(vsearchFastafile) + ".clustered.uc";
        string logfile = m->getSimpleName(vsearchFastafile) + ".clustered.log";
        vsearchDriver(vsearchFastafile, ucVsearchFile, logfile);
        
        if (m->control_pressed) { m->mothurRemove(ucVsearchFile); m->mothurRemove(logfile);  m->mothurRemove(vsearchFastafile); return 0; }
        
        if (outputDir == "") { outputDir += m->hasPath(distfile); }
        fileroot = outputDir + m->getRootName(m->getSimpleName(distfile));
        tag = method;
        
        map<string, string> variables;
        variables["[filename]"] = fileroot;
        variables["[clustertag]"] = tag;
        string sabundFileName = getOutputFileName("sabund", variables);
        string rabundFileName = getOutputFileName("rabund", variables);
        if (countfile != "") { variables["[tag2]"] = "unique_list"; }
        string listFileName = getOutputFileName("list", variables);
        outputNames.push_back(listFileName); outputTypes["list"].push_back(listFileName);
        if (countfile == "") {
            outputNames.push_back(sabundFileName); outputTypes["sabund"].push_back(sabundFileName);
            outputNames.push_back(rabundFileName); outputTypes["rabund"].push_back(rabundFileName);
        }
        
        //Convert outputted *.uc file into a list file
        vParse->createListFile(ucVsearchFile, listFileName, sabundFileName, rabundFileName, vParse->getNumBins(logfile), toString(1.0-cutoff));  delete vParse;
        
        //remove temp files
        m->mothurRemove(ucVsearchFile); m->mothurRemove(logfile);  m->mothurRemove(vsearchFastafile);
        
        if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); } return 0; }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterCommand", "runVsearchCluster");
        exit(1);
    }
}
//**********************************************************************************************************************

int ClusterCommand::vsearchDriver(string inputFile, string ucClusteredFile, string logfile){
    try {
        
        //vsearch --maxaccepts 16 --usersort --id 0.97 --minseqlength 30 --wordlength 8 --uc $ROOT.clustered.uc --cluster_smallmem $ROOT.sorted.fna --maxrejects 64 --strand both --log $ROOT.clustered.log --sizeorder

        
        //no sizeorder for dgc
        
        ucClusteredFile = m->getFullPathName(ucClusteredFile);
        inputFile = m->getFullPathName(inputFile);
        logfile = m->getFullPathName(logfile);
        
        //to allow for spaces in the path
        ucClusteredFile = "\"" + ucClusteredFile + "\"";
        inputFile = "\"" + inputFile + "\"";
        logfile = "\"" + logfile + "\"";
        
        vector<char*> cPara;
        
        string vsearchCommand = vsearchLocation;
        vsearchCommand = "\"" + vsearchCommand + "\" ";
        
        vector<char*> vsearchParameters;
        char* vsearchParameter = new char[vsearchCommand.length()+1];  vsearchParameter[0] = '\0'; strncat(vsearchParameter, vsearchCommand.c_str(), vsearchCommand.length());
        vsearchParameters.push_back(vsearchParameter);
        
        //--maxaccepts=16
        char* maxaccepts = new char[16];  maxaccepts[0] = '\0'; strncat(maxaccepts, "--maxaccepts=16", 15);
        vsearchParameters.push_back(maxaccepts);
        
        //--usersort
        char* usersort = new char[11];  usersort[0] = '\0'; strncat(usersort, "--usersort", 10);
        vsearchParameters.push_back(usersort);
        
        //--id=0.97
        cutoff = abs(1.0 - cutoff); string cutoffString = toString(cutoff);
        if (cutoffString.length() > 4) {  cutoffString = cutoffString.substr(0, 4);  }
        else if (cutoffString.length() < 4)  {  for (int i = cutoffString.length(); i < 4; i++)  { cutoffString += "0";  } }
        
        cutoffString = "--id=" +  cutoffString;
        char* cutoffParameter = new char[cutoffString.length()+1];  cutoffParameter[0] = '\0'; strncat(cutoffParameter, cutoffString.c_str(), cutoffString.length());
        vsearchParameters.push_back(cutoffParameter);
        
        //--minseqlength=30
        char* minseqlength = new char[18];  minseqlength[0] = '\0'; strncat(minseqlength, "--minseqlength=30", 17);
        vsearchParameters.push_back(minseqlength);
        
        //--wordlength=8
        char* wordlength = new char[15];  wordlength[0] = '\0'; strncat(wordlength, "--wordlength=8", 14);
        vsearchParameters.push_back(wordlength);

        //--uc=$ROOT.clustered.uc
        string tempIn = "--uc=" + ucClusteredFile;
        char* uc = new char[tempIn.length()+1];  uc[0] = '\0'; strncat(uc, tempIn.c_str(), tempIn.length());
        vsearchParameters.push_back(uc);

        //--cluster_smallmem $ROOT.sorted.fna
        string tempSorted = "--cluster_smallmem=" + inputFile;
        char* cluster_smallmen = new char[tempSorted.length()+1];  cluster_smallmen[0] = '\0'; strncat(cluster_smallmen, tempSorted.c_str(), tempSorted.length());
        vsearchParameters.push_back(cluster_smallmen);
        
        //--maxrejects=64
        char* maxrejects = new char[16];  maxrejects[0] = '\0'; strncat(maxrejects, "--maxrejects=64", 15);
        vsearchParameters.push_back(maxrejects);
        
        //--strand=both
        char* strand = new char[14];  strand[0] = '\0'; strncat(strand, "--strand=both", 13);
        vsearchParameters.push_back(strand);
        
        //--log=$ROOT.clustered.log
        string tempLog = "--log=" + logfile;
        char* log = new char[tempLog.length()+1];  log[0] = '\0'; strncat(log, tempLog.c_str(), tempLog.length());
        vsearchParameters.push_back(log);

        if (method == "agc") {
            //--sizeorder
            char* sizeorder = new char[12];  sizeorder[0] = '\0'; strncat(sizeorder, "--sizeorder", 11);
            vsearchParameters.push_back(sizeorder);
         }

        if (m->debug) {  for(int i = 0; i < vsearchParameters.size(); i++)  { cout << vsearchParameters[i]; } cout << endl;  }
        
        string commandString = "";
        for (int i = 0; i < vsearchParameters.size(); i++) {    commandString += toString(vsearchParameters[i]) + " "; }
 
        //cout << "commandString = " << commandString << endl;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
        commandString = "\"" + commandString + "\"";
#endif
        if (m->debug) { m->mothurOut("[DEBUG]: vsearch cluster command = " + commandString + ".\n"); }
        system(commandString.c_str());
 
        //free memory
        for(int i = 0; i < vsearchParameters.size(); i++)  {  delete vsearchParameters[i];  }
        
        //remove "" from filenames
        ucClusteredFile = ucClusteredFile.substr(1, ucClusteredFile.length()-2);
        inputFile = inputFile.substr(1, inputFile.length()-2);
        logfile = logfile.substr(1, logfile.length()-2);

        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterCommand", "vsearchDriver");
        exit(1);
    }
}
//**********************************************************************************************************************

int ClusterCommand::runMothurCluster(){
    try {
        
        ReadMatrix* read;
        if (format == "column") { read = new ReadColumnMatrix(columnfile, sim); }	//sim indicates whether its a similarity matrix
        else if (format == "phylip") { read = new ReadPhylipMatrix(phylipfile, sim); }
        
        read->setCutoff(cutoff);
        
        NameAssignment* nameMap = NULL;
        CountTable* ct = NULL;
        map<string, int> counts;
        if(namefile != ""){
            nameMap = new NameAssignment(namefile);
            nameMap->readMap();
            read->read(nameMap);
        }else if (countfile != "") {
            ct = new CountTable();
            ct->readTable(countfile, false, false);
            read->read(ct);
            counts = ct->getNameMap();
        }else { read->read(nameMap); }
        
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
        if (method == "furthest")	{	cluster = new CompleteLinkage(rabund, list, matrix, cutoff, method, adjust); }
        else if(method == "nearest"){	cluster = new SingleLinkage(rabund, list, matrix, cutoff, method, adjust); }
        else if(method == "average"){	cluster = new AverageLinkage(rabund, list, matrix, cutoff, method, adjust);	}
        else if(method == "weighted"){	cluster = new WeightedLinkage(rabund, list, matrix, cutoff, method, adjust);	}
        tag = cluster->getTag();
        
        if (outputDir == "") { outputDir += m->hasPath(distfile); }
        fileroot = outputDir + m->getRootName(m->getSimpleName(distfile));
        
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
        
        if (print_start && m->isTrue(timing)) {
            m->mothurOut("Clustering (" + tag + ") for distance " + toString(previousDist) + "/" + toString(rndPreviousDist)
                         + "\t(precision: " + toString(precision) + ", Nodes: " + toString(matrix->getNNodes()) + ")");
            cout.flush();
            print_start = false;
        }
        
        if(previousDist <= 0.0000){
            printData("unique", counts);
        }
        else if(rndPreviousDist<cutoff){
            printData(toString(rndPreviousDist, length-1), counts);
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

        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterCommand", "runMothurCluster");
        exit(1);
    }
}
//**********************************************************************************************************************

void ClusterCommand::printData(string label, map<string, int>& counts){
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
        if(countfile != "") {
            oldList.print(listFile, counts);
        }else {
            oldList.print(listFile);
        }
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

int ClusterCommand::runOptiCluster(){
    try {
        if (cutoffNotSet) {  m->mothurOut("\nYou did not set a cutoff, using 0.03.\n"); cutoff = 0.03; cutoff += (5 / (precision * 10.0)); }
        cutoff -= (5 / (precision * 10.0));
        
        string nameOrCount = "name";
        string thisNamefile = namefile;
        map<string, int> counts;
        if (countfile != "") { nameOrCount = "count"; thisNamefile = countfile; CountTable ct; ct.readTable(countfile, false, false); counts = ct.getNameMap(); }
        string distfile = columnfile;
        if (format == "phylip") { distfile = phylipfile; }
        
        OptiMatrix matrix(distfile, thisNamefile, nameOrCount, cutoff, false);
        
        OptiCluster cluster(&matrix, metric, stableMetric);
        tag = cluster.getTag();
        
        m->mothurOutEndLine(); m->mothurOut("Clustering " + distfile); m->mothurOutEndLine();
        
        if (outputDir == "") { outputDir += m->hasPath(distfile); }
        fileroot = outputDir + m->getRootName(m->getSimpleName(distfile));
        
        string listFileName = fileroot+ tag + ".list";
        
        ofstream listFile;
        m->openOutputFile(listFileName,	listFile);
        outputNames.push_back(listFileName); outputTypes["list"].push_back(listFileName);
        
        int iters = 0;
        double listVectorMetric = 0; //worst state
        double delta = 1;
        
        cluster.initialize(listVectorMetric);
    
        while ((delta > stableMetric) && (iters < maxIters)) {
            
            if (m->control_pressed) { break; }
            double oldMetric = listVectorMetric;
            
            cluster.update(listVectorMetric);

            delta = abs(oldMetric - listVectorMetric);
            iters++;
        }
        
        ListVector* list = cluster.getList();
        list->setLabel(toString(cutoff));
        if (!m->printedListHeaders) { m->listBinLabelsInFile.clear(); list->printHeaders(listFile); }
        if(countfile != "") { list->print(listFile, counts); }
        else { list->print(listFile); }
        listFile.close();
        
        map<string, string> variables;
        variables["[filename]"] = fileroot;
        variables["[clustertag]"] = tag;
        string sabundFileName = getOutputFileName("sabund", variables);
        string rabundFileName = getOutputFileName("rabund", variables);
        
        if (countfile == "") {
            m->openOutputFile(sabundFileName,	sabundFile);
            m->openOutputFile(rabundFileName,	rabundFile);
            outputNames.push_back(sabundFileName); outputTypes["sabund"].push_back(sabundFileName);
            outputNames.push_back(rabundFileName); outputTypes["rabund"].push_back(rabundFileName);
            
            SAbundVector sabund = list->getSAbundVector();
            sabund.print(sabundFile);
            sabundFile.close();
            
            RAbundVector rabund = list->getRAbundVector();
            rabund.print(rabundFile);
            rabundFile.close();
        }
        delete list;
        
        
        m->mothurOut("\ncutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");
        vector<double> results = cluster.getStats();
        m->mothurOut(toString(cutoff) + "\t");
        for (int i = 0; i < results.size(); i++) { m->mothurOut(toString(results[i]) + "\t");  }
        m->mothurOut("\n\n");
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterCommand", "runOptiCluster");
        exit(1);
    }
    
}
//**********************************************************************************************************************
