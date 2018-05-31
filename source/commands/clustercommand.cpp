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
#include "sensspeccommand.h"
#include "mcc.hpp"
#include "sensitivity.hpp"
#include "specificity.hpp"
#include "fdr.hpp"
#include "npv.hpp"
#include "ppv.hpp"
#include "f1score.hpp"
#include "tp.hpp"
#include "fp.hpp"
#include "fpfn.hpp"
#include "tptn.hpp"
#include "tn.hpp"
#include "fn.hpp"
#include "accuracy.hpp"



//**********************************************************************************************************************
vector<string> ClusterCommand::setParameters(){	
	try {
        CommandParameter pphylip("phylip", "InputTypes", "", "", "PhylipColumnFasta", "PhylipColumnFasta", "none","list",false,false,true); parameters.push_back(pphylip);
        CommandParameter pfasta("fasta", "InputTypes", "", "", "PhylipColumnFasta", "PhylipColumnFasta", "FastaTaxName","list",false,false,true); parameters.push_back(pfasta);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "ColumnName-FastaTaxName","rabund-sabund",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount", "none", "","",false,false,true); parameters.push_back(pcount);
        CommandParameter pcolumn("column", "InputTypes", "", "", "PhylipColumnFasta", "PhylipColumnFasta", "ColumnName","list",false,false,true); parameters.push_back(pcolumn);
		CommandParameter pcutoff("cutoff", "Number", "", "0.03", "", "", "","",false,false,true); parameters.push_back(pcutoff);
		CommandParameter pprecision("precision", "Number", "", "100", "", "", "","",false,false); parameters.push_back(pprecision);
		CommandParameter pmethod("method", "Multiple", "furthest-nearest-average-weighted-agc-dgc-opti-unique", "opti", "", "", "","",false,false,true); parameters.push_back(pmethod);
        CommandParameter pinitialize("initialize", "Multiple", "oneotu-singleton", "singleton", "", "", "","",false,false,true); parameters.push_back(pinitialize);
        CommandParameter pmetric("metric", "Multiple", "mcc-sens-spec-tptn-fpfn-tp-tn-fp-fn-f1score-accuracy-ppv-npv-fdr", "mcc", "", "", "","",false,false,true); parameters.push_back(pmetric);
        CommandParameter pmetriccutoff("delta", "Number", "", "0.0001", "", "", "","",false,false,true); parameters.push_back(pmetriccutoff);
        CommandParameter piters("iters", "Number", "", "100", "", "", "","",false,false,true); parameters.push_back(piters);
		CommandParameter pshowabund("showabund", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pshowabund);
		CommandParameter ptiming("timing", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(ptiming);
		CommandParameter psim("sim", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(psim);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        //CommandParameter padjust("adjust", "String", "", "F", "", "", "","",false,false); parameters.push_back(padjust);
        CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
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
		helpString += "The cluster command parameter options are phylip, column, name, count, method, cutoff, precision, sim, showabund, timing, metric, iters, initialize. Fasta or Phylip or column and name are required.\n";
        helpString += "The phylip and column parameter allow you to enter your distance file. \n";
        helpString += "The fasta parameter allows you to enter your fasta file for use with the agc or dgc methods. \n";
        helpString += "The name parameter allows you to enter your name file. \n";
        helpString += "The count parameter allows you to enter your count file. \n A count or name file is required if your distance file is in column format.\n";
        helpString += "The iters parameter allow you to set the maxiters for the opticluster method. \n";
        helpString += "The metric parameter allows to select the metric in the opticluster method. Options are Matthews correlation coefficient (mcc), sensitivity (sens), specificity (spec), true positives + true negatives (tptn), false positives + false negatives (fpfn), true positives (tp), true negative (tn), false positive (fp), false negative (fn), f1score (f1score), accuracy (accuracy), positive predictive value (ppv), negative predictive value (npv), false discovery rate (fdr). Default=mcc.\n";
        helpString += "The initialize parameter allows to select the initial randomization for the opticluster method. Options are singleton, meaning each sequence is randomly assigned to its own OTU, or oneotu meaning all sequences are assigned to one otu. Default=singleton.\n";
        helpString += "The delta parameter allows to set the stable value for the metric in the opticluster method (delta=0.0001). \n";
        helpString += "The method parameter allows you to enter your clustering mothod. Options are furthest, nearest, average, weighted, agc, dgc, unique and opti. Default=opti.  The agc and dgc methods require a fasta file.";
        helpString += "The processors parameter allows you to specify the number of processors to use. The default is 1.\n";
       helpString += "The cluster command should be in the following format: \n";
		helpString += "cluster(method=yourMethod, cutoff=yourCutoff, precision=yourPrecision) \n";
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
        else if (type == "sensspec") {  pattern = "[filename],[clustertag],sensspec"; }
        else if (type == "steps") {  pattern = "[filename],[clustertag],steps"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
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
        outputTypes["sensspec"] = tempOutNames;
		outputTypes["rabund"] = tempOutNames;
		outputTypes["sabund"] = tempOutNames;
        outputTypes["steps"] = tempOutNames;
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
            outputTypes["sensspec"] = tempOutNames;
			outputTypes["rabund"] = tempOutNames;
			outputTypes["sabund"] = tempOutNames;
            outputTypes["steps"] = tempOutNames;
		
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";		}
			
			inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("phylip");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["phylip"] = inputDir + it->second;		}
				}
				
				it = parameters.find("column");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["column"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
                
                it = parameters.find("fasta");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
                }
			}
			
			//check for required parameters
			phylipfile = validParameter.validFile(parameters, "phylip");
			if (phylipfile == "not open") { phylipfile = ""; abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; }	
			else {  distfile = phylipfile;  format = "phylip"; 	current->setPhylipFile(phylipfile); }
			
			columnfile = validParameter.validFile(parameters, "column");
			if (columnfile == "not open") { columnfile = ""; abort = true; }	
			else if (columnfile == "not found") { columnfile = ""; }
			else {  distfile = columnfile; format = "column"; current->setColumnFile(columnfile);	}
			
            fastafile = validParameter.validFile(parameters, "fasta");
            if (fastafile == "not open") { abort = true; }
            else if (fastafile == "not found") { fastafile = ""; }
            else { distfile = fastafile;  format = "fasta"; current->setFastaFile(fastafile); }
            
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			else { current->setNameFile(namefile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { abort = true; countfile = ""; }	
			else if (countfile == "not found") { countfile = ""; }
			else { current->setCountFile(countfile); }
			
            method = validParameter.valid(parameters, "method");
            if (method == "not found") {  method = "opti";}
            
            if ((method == "furthest") || (method == "nearest") || (method == "average") || (method == "weighted") || (method == "agc") || (method == "dgc") || (method == "opti") || (method == "unique")) { }
            else { m->mothurOut("[ERROR]: Not a valid clustering method.  Valid clustering algorithms are furthest, nearest, average, weighted, agc, dgc, unique and opti."); m->mothurOutEndLine(); abort = true; }
            
            if (method != "unique") {
                if ((phylipfile == "") && (columnfile == "") && (fastafile == "")) {
                    //is there are current file available for either of these?
                    //give priority to column, then phylip
                    columnfile = current->getColumnFile();
                    if (columnfile != "") {  distfile = columnfile; format = "column"; m->mothurOut("Using " + columnfile + " as input file for the column parameter."); m->mothurOutEndLine(); }
                    else {
                        phylipfile = current->getPhylipFile();
                        if (phylipfile != "") { distfile = phylipfile;  format = "phylip"; m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter."); m->mothurOutEndLine(); }
                        else {
                            fastafile = current->getFastaFile();
                            if (fastafile != "") {  distfile = fastafile; format = "fasta"; m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
                            else {
                                m->mothurOut("No valid current files. You must provide a phylip, column or fasta file before you can use the cluster command, unless using the unique method."); m->mothurOutEndLine();
                                abort = true;
                            }
                        }
                    }
                }
                else if (((phylipfile != "") && (columnfile != "")) || ((phylipfile != "") && (fastafile != "")) || ((fastafile != "") && (columnfile != "")))  { m->mothurOut("When executing a cluster command you must enter ONLY ONE of the following: phylip, column or fasta."); m->mothurOutEndLine(); abort = true; }
                
                if (columnfile != "") {
                    if ((namefile == "") && (countfile == "")){
                        namefile = current->getNameFile();
                        if (namefile != "") {  m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
                        else {
                            countfile = current->getCountFile();
                            if (countfile != "") {  m->mothurOut("Using " + countfile + " as input file for the count parameter."); m->mothurOutEndLine(); }
                            else { 
                                m->mothurOut("You need to provide a namefile or countfile if you are going to use the column format."); m->mothurOutEndLine(); 
                                abort = true; 
                            }	
                        }	
                    }
                }
                if ((method != "agc") && (method != "dgc")) {
                    if ((columnfile == "") && (phylipfile == "")) {
                        m->mothurOut("[ERROR]: You must provide a distance file unless you are using the agc, dgc or unique clustering methods, aborting\n."); abort = true;
                    }
                }
            }else {
                if ((countfile == "") && (namefile == "")) {
                    countfile = current->getCountFile();
                    if (countfile != "") { distfile = countfile;  format = "count"; m->mothurOut("Using " + countfile + " as input file for the count parameter."); m->mothurOutEndLine(); }
                    else {
                        namefile = current->getNameFile();
                        if (namefile != "") {  distfile = namefile; format = "name"; m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
                        else {
                            m->mothurOut("No valid current files. You must provide a count or name file before you can use the cluster command with the unique method."); m->mothurOutEndLine();
                            abort = true;
                        }
                    }

                }
                else if(countfile != "")    { format = "count"; }
                else if(namefile != "")     { format = "name";  }
            }
            if ((countfile != "") && (namefile != "")) { m->mothurOut("When executing a cluster command you must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
            
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			//get user cutoff and precision or use defaults
			string temp;
			temp = validParameter.valid(parameters, "precision");
			if (temp == "not found") { temp = "100"; }
			//saves precision legnth for formatting below
			length = temp.length();
			util.mothurConvert(temp, precision);
			
			temp = validParameter.valid(parameters, "sim");				if (temp == "not found") { temp = "F"; }
			sim = util.isTrue(temp); 
			
            temp = validParameter.valid(parameters, "delta");		if (temp == "not found")  { temp = "0.0001"; }
            util.mothurConvert(temp, stableMetric);
            
            metricName = validParameter.valid(parameters, "metric");		if (metricName == "not found") { metricName = "mcc"; }
            
            if ((metricName == "mcc") || (metricName == "sens") || (metricName == "spec") || (metricName == "tptn") || (metricName == "tp") || (metricName == "tn") || (metricName == "fp") || (metricName == "fn") || (metricName == "f1score") || (metricName == "accuracy") || (metricName == "ppv") || (metricName == "npv") || (metricName == "fdr") || (metricName == "fpfn") ){ }
            else { m->mothurOut("[ERROR]: Not a valid metric.  Valid metrics are mcc, sens, spec, tp, tn, fp, fn, tptn, fpfn, f1score, accuracy, ppv, npv, fdr."); m->mothurOutEndLine(); abort = true; }
            
            initialize = validParameter.valid(parameters, "initialize");		if (initialize == "not found") { initialize = "singleton"; }
            
            if ((initialize == "singleton") || (initialize == "oneotu")){ }
            else { m->mothurOut("[ERROR]: Not a valid initialization.  Valid initializations are singleton and oneotu."); m->mothurOutEndLine(); abort = true; }

            temp = validParameter.valid(parameters, "iters");		if (temp == "not found")  { temp = "100"; }
            util.mothurConvert(temp, maxIters);
            
            adjust=-1.0;
			
            bool setProcessors = true;
            temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){ setProcessors=false;	temp = current->getProcessors();	}
            processors = current->setProcessors(temp);
            
            if ((method == "agc") || (method == "dgc")) {
                if (fastafile == "") { m->mothurOut("[ERROR]: You must provide a fasta file when using the agc or dgc clustering methods, aborting\n."); abort = true;}
            }else if (setProcessors) {
                m->mothurOut("[WARNING]: You can only use the processors option when using the agc or dgc clustering methods. Using 1 processor.\n.");
            }
            
            cutOffSet = false;
            temp = validParameter.valid(parameters, "cutoff");
            if (temp == "not found") { if ((method == "opti") || (method == "agc") || (method == "dgc")) { temp = "0.03"; }else { temp = "0.15"; } }
            else { cutOffSet = true; }
            util.mothurConvert(temp, cutoff);
            
			showabund = validParameter.valid(parameters, "showabund");
			if (showabund == "not found") { showabund = "T"; }

			timing = validParameter.valid(parameters, "timing");
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
	
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		//phylip file given and cutoff not given - use cluster.classic because it uses less memory and is faster
		if ((format == "phylip") && (!cutOffSet) && (method != "opti")) {
			m->mothurOutEndLine(); m->mothurOut("You are using a phylip file and no cutoff.  I will run cluster.classic to save memory and time."); m->mothurOutEndLine();
			
			//run unique.seqs for deconvolute results
			string inputString = "phylip=" + distfile;
			if (namefile != "") { inputString += ", name=" + namefile; }
            else if (countfile != "") { inputString += ", count=" + countfile; }
			inputString += ", precision=" + toString(precision);
			inputString += ", method=" + method;
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
        
        if (format == "fasta")          {   runVsearchCluster();    }
        else if (method == "opti")      {   runOptiCluster();       }
        else if (method == "unique")    {   runUniqueCluster();     }
        else                            {   runMothurCluster();     }
        
		if (m->getControl_pressed()) { 	for (int j = 0; j < outputNames.size(); j++) { util.mothurRemove(outputNames[j]); }  return 0; }
        
        m->mothurOut("It took " + toString(time(NULL) - estart) + " seconds to cluster"); m->mothurOutEndLine();
        
		//set list file as new current listfile
		string currentName = "";
		itTypes = outputTypes.find("list");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setListFile(currentName); }
		}
		
		//set rabund file as new current rabundfile
		itTypes = outputTypes.find("rabund");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setRabundFile(currentName); }
		}
		
		//set sabund file as new current sabundfile
		itTypes = outputTypes.find("sabund");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setSabundFile(currentName); }
		}
		
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();

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
        //look for vsearch exe
        string path = current->getProgramPath();
        string vsearchCommand = path + PATH_SEPARATOR + "vsearch" + EXECUTABLE_EXT;
#if defined NON_WINDOWS
        if (m->getDebug()) {
            m->mothurOut("[DEBUG]: vsearch location using \"which vsearch\" = ");
            Command* newCommand = new SystemCommand("which vsearch"); m->mothurOutEndLine();
            newCommand->execute(); delete newCommand;
            m->mothurOut("[DEBUG]: Mothur's location using \"which mothur\" = ");
            newCommand = new SystemCommand("which mothur"); m->mothurOutEndLine();
            newCommand->execute(); delete newCommand;
        }
#endif
        
        //test to make sure vsearch exists
        ifstream in;
        vsearchCommand = util.getFullPathName(vsearchCommand);
        bool ableToOpen = util.openInputFile(vsearchCommand, in, "no error"); in.close();
        if(!ableToOpen) {
            m->mothurOut(vsearchCommand + " file does not exist. Checking path... \n");

            ifstream in2;
            string programName = "vsearch"; programName += EXECUTABLE_EXT;
            string uLocation = util.findProgramPath(programName);
            ableToOpen = util.openInputFile(uLocation+programName, in2, "no error"); in2.close();
            
            if(!ableToOpen) { m->mothurOut("[ERROR]: " + uLocation + " file does not exist. mothur requires the vsearch executable.\n");  m->setControl_pressed(true);  }
            else {  m->mothurOut("Found vsearch in your path, using " + uLocation + "\n");vsearchLocation = uLocation+programName; }
        }else {  vsearchLocation = vsearchCommand; }
        
        if (m->getControl_pressed()) {  return 0; }
        
        vsearchLocation = util.getFullPathName(vsearchLocation);
        
        string vsearchFastafile = ""; VsearchFileParser* vParse;
        if ((namefile == "") && (countfile == ""))  { vParse = new VsearchFileParser(fastafile);                        }
        else if (namefile != "")                    { vParse = new VsearchFileParser(fastafile, namefile, "name");      }
        else if (countfile != "")                   { vParse = new VsearchFileParser(fastafile, countfile, "count");    }
        else                                        { m->mothurOut("[ERROR]: Opps, should never get here. ClusterCommand::runVsearchCluster() \n"); m->setControl_pressed(true); }
    
        if (m->getControl_pressed()) {  return 0; }
        
        vsearchFastafile = vParse->getVsearchFile();
        
        if (cutoff > 1.0) {  m->mothurOut("You did not set a cutoff, using 0.03.\n"); cutoff = 0.03; }
        
        //Run vsearch
        string ucVsearchFile = util.getSimpleName(vsearchFastafile) + ".clustered.uc";
        string logfile = util.getSimpleName(vsearchFastafile) + ".clustered.log";
        vsearchDriver(vsearchFastafile, ucVsearchFile, logfile);
        
        if (m->getControl_pressed()) { util.mothurRemove(ucVsearchFile); util.mothurRemove(logfile);  util.mothurRemove(vsearchFastafile); return 0; }
        
        if (outputDir == "") { outputDir += util.hasPath(distfile); }
        fileroot = outputDir + util.getRootName(util.getSimpleName(distfile));
        tag = method;
        
        map<string, string> variables;
        variables["[filename]"] = fileroot;
        variables["[clustertag]"] = tag;
        string sabundFileName = getOutputFileName("sabund", variables);
        string rabundFileName = getOutputFileName("rabund", variables);
        //if (countfile != "") { variables["[tag2]"] = "unique_list"; }
        string listFileName = getOutputFileName("list", variables);
        outputNames.push_back(listFileName); outputTypes["list"].push_back(listFileName);
        if (countfile == "") {
            outputNames.push_back(sabundFileName); outputTypes["sabund"].push_back(sabundFileName);
            outputNames.push_back(rabundFileName); outputTypes["rabund"].push_back(rabundFileName);
        }
        
        //Convert outputted *.uc file into a list file
        vParse->createListFile(ucVsearchFile, listFileName, sabundFileName, rabundFileName, vParse->getNumBins(logfile), toString(1.0-cutoff));  delete vParse;
        
        //remove temp files
        util.mothurRemove(ucVsearchFile); util.mothurRemove(logfile);  util.mothurRemove(vsearchFastafile);
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); } return 0; }
        
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
        
        ucClusteredFile = util.getFullPathName(ucClusteredFile);
        inputFile = util.getFullPathName(inputFile);
        logfile = util.getFullPathName(logfile);
        
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
        
        //--threads=1
        string processorsString = "--threads=" + toString(processors);
        char* processorsParameter = new char[processorsString.length()+1];  processorsParameter[0] = '\0'; strncat(processorsParameter, processorsString.c_str(), processorsString.length());
        vsearchParameters.push_back(processorsParameter);
        
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

        if (m->getDebug()) {  for(int i = 0; i < vsearchParameters.size(); i++)  { cout << vsearchParameters[i]; } cout << endl;  }
        
        string commandString = "";
        for (int i = 0; i < vsearchParameters.size(); i++) {    commandString += toString(vsearchParameters[i]) + " "; }
 
        //cout << "commandString = " << commandString << endl;
        //exit(1);
        
#if defined NON_WINDOWS
#else
        commandString = "\"" + commandString + "\"";
#endif
        if (m->getDebug()) { m->mothurOut("[DEBUG]: vsearch cluster command = " + commandString + ".\n"); }
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
        
        if (m->getControl_pressed()) { //clean up
            delete list; delete matrix; delete rabund; if(countfile == ""){rabundFile.close(); sabundFile.close();  util.mothurRemove((fileroot+ tag + ".rabund")); util.mothurRemove((fileroot+ tag + ".sabund")); }
            listFile.close(); util.mothurRemove((fileroot+ tag + ".list")); outputTypes.clear(); return 0;
        }
        
        //create cluster
        if (method == "furthest")	{	cluster = new CompleteLinkage(rabund, list, matrix, cutoff, method, adjust); }
        else if(method == "nearest"){	cluster = new SingleLinkage(rabund, list, matrix, cutoff, method, adjust); }
        else if(method == "average"){	cluster = new AverageLinkage(rabund, list, matrix, cutoff, method, adjust);	}
        else if(method == "weighted"){	cluster = new WeightedLinkage(rabund, list, matrix, cutoff, method, adjust);	}
        tag = cluster->getTag();
        
        if (outputDir == "") { outputDir += util.hasPath(distfile); }
        fileroot = outputDir + util.getRootName(util.getSimpleName(distfile));
        
        map<string, string> variables;
        variables["[filename]"] = fileroot;
        variables["[clustertag]"] = tag;
        string sabundFileName = getOutputFileName("sabund", variables);
        string rabundFileName = getOutputFileName("rabund", variables);
        //if (countfile != "") { variables["[tag2]"] = "unique_list"; }
        string listFileName = getOutputFileName("list", variables);
        
        if (countfile == "") {
            util.openOutputFile(sabundFileName,	sabundFile);
            util.openOutputFile(rabundFileName,	rabundFile);
            outputNames.push_back(sabundFileName); outputTypes["sabund"].push_back(sabundFileName);
            outputNames.push_back(rabundFileName); outputTypes["rabund"].push_back(rabundFileName);
            
        }
        util.openOutputFile(listFileName,	listFile);
        outputNames.push_back(listFileName); outputTypes["list"].push_back(listFileName);
        
        float previousDist = 0.00000;
        float rndPreviousDist = 0.00000;
        oldRAbund = *rabund;
        oldList = *list;
        
        print_start = true;
        start = time(NULL);
        loops = 0;
        double saveCutoff = cutoff;
        bool printHeaders = true;
        
        while (matrix->getSmallDist() < cutoff && matrix->getNNodes() > 0){
            
            if (m->getControl_pressed()) { //clean up
                delete list; delete matrix; delete rabund; delete cluster;
                if(countfile == "") {rabundFile.close(); sabundFile.close();  util.mothurRemove((fileroot+ tag + ".rabund")); util.mothurRemove((fileroot+ tag + ".sabund")); }
                listFile.close(); util.mothurRemove((fileroot+ tag + ".list")); outputTypes.clear(); return 0;
            }
            
            if (print_start && util.isTrue(timing)) {
                m->mothurOut("Clustering (" + tag + ") dist " + toString(matrix->getSmallDist()) + "/"
                             + toString(util.roundDist(matrix->getSmallDist(), precision))
                             + "\t(precision: " + toString(precision) + ", Nodes: " + toString(matrix->getNNodes()) + ")");
                cout.flush();
                print_start = false;
            }
            
            loops++;
            
            cluster->update(cutoff);
            
            float dist = matrix->getSmallDist();
            float rndDist = util.ceilDist(dist, precision);
            
            if(previousDist <= 0.0000 && dist != previousDist){
                printData("unique", counts, printHeaders);
            }
            else if(rndDist != rndPreviousDist){
                printData(toString(rndPreviousDist,  length-1), counts, printHeaders);
            }
            
            previousDist = dist;
            rndPreviousDist = rndDist;
            oldRAbund = *rabund;
            oldList = *list;
        }
        
        if (print_start && util.isTrue(timing)) {
            m->mothurOut("Clustering (" + tag + ") for distance " + toString(previousDist) + "/" + toString(rndPreviousDist)
                         + "\t(precision: " + toString(precision) + ", Nodes: " + toString(matrix->getNNodes()) + ")");
            cout.flush();
            print_start = false;
        }
        
        if(previousDist <= 0.0000){
            printData("unique", counts, printHeaders);
        }
        else if(rndPreviousDist<cutoff){
            printData(toString(rndPreviousDist, length-1), counts, printHeaders);
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
            saveCutoff = util.ceilDist(saveCutoff, precision);
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

void ClusterCommand::printData(string label, map<string, int>& counts, bool& ph){
	try {
        oldList.setPrintedLabels(ph); ph = false;
        
		if (util.isTrue(timing)) {
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
       
        if (util.isTrue(showabund)) {
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
            if (m->getControl_pressed()) { break; }
            vector<string> binNames;
            string bin = list->get(i);
            util.splitAtComma(bin, binNames);
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
        if (!cutOffSet) {  m->mothurOut("\nYou did not set a cutoff, using 0.03.\n"); cutoff = 0.03;  }
        
        string nameOrCount = "";
        string thisNamefile = "";
        map<string, int> counts;
        if (countfile != "") { nameOrCount = "count"; thisNamefile = countfile; CountTable ct; ct.readTable(countfile, false, false); counts = ct.getNameMap(); }
        else if (namefile != "") { nameOrCount = "name"; thisNamefile = namefile; }
        
        string distfile = columnfile;
        if (format == "phylip") { distfile = phylipfile; }
        
        
        OptiMatrix matrix(distfile, thisNamefile, nameOrCount, format, cutoff, false);
    
        ClusterMetric* metric = NULL;
        if (metricName == "mcc")             { metric = new MCC();              }
        else if (metricName == "sens")       { metric = new Sensitivity();      }
        else if (metricName == "spec")       { metric = new Specificity();      }
        else if (metricName == "tptn")       { metric = new TPTN();             }
        else if (metricName == "tp")         { metric = new TP();               }
        else if (metricName == "tn")         { metric = new TN();               }
        else if (metricName == "fp")         { metric = new FP();               }
        else if (metricName == "fn")         { metric = new FN();               }
        else if (metricName == "f1score")    { metric = new F1Score();          }
        else if (metricName == "accuracy")   { metric = new Accuracy();         }
        else if (metricName == "ppv")        { metric = new PPV();              }
        else if (metricName == "npv")        { metric = new NPV();              }
        else if (metricName == "fdr")        { metric = new FDR();              }
        else if (metricName == "fpfn")       { metric = new FPFN();             }
        
        
        OptiCluster cluster(&matrix, metric, 0);
        tag = cluster.getTag();
        
        m->mothurOutEndLine(); m->mothurOut("Clustering " + distfile); m->mothurOutEndLine();
        
        if (outputDir == "") { outputDir += util.hasPath(distfile); }
        fileroot = outputDir + util.getRootName(util.getSimpleName(distfile));
        
        string listFileName = fileroot+ tag + ".list";
        
        ofstream listFile;
        util.openOutputFile(listFileName,	listFile);
        outputNames.push_back(listFileName); outputTypes["list"].push_back(listFileName);
        
        map<string, string> variables;
        variables["[filename]"] = fileroot;
        variables["[clustertag]"] = tag;
        string outputName = getOutputFileName("steps", variables);
        outputNames.push_back(outputName); outputTypes["steps"].push_back(outputName);
        ofstream outStep;
        util.openOutputFile(outputName, outStep);
        
        int iters = 0;
        double listVectorMetric = 0; //worst state
        double delta = 1;
        
        cluster.initialize(listVectorMetric, true, initialize);
        
        long long numBins = cluster.getNumBins();
        m->mothurOut("\n\niter\ttime\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");
        outStep << "iter\ttime\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";
        long long tp, tn, fp, fn;
        vector<double> results = cluster.getStats(tp, tn, fp, fn);
        m->mothurOut("0\t0\t" + toString(cutoff) + "\t" + toString(numBins) + "\t"+ toString(cutoff) + "\t" + toString(tp) + "\t" + toString(tn) + "\t" + toString(fp) + "\t" + toString(fn) + "\t");
        outStep << "0\t0\t" + toString(cutoff) + "\t" + toString(numBins) + "\t" + toString(cutoff) + "\t" << tp << '\t' << tn << '\t' << fp << '\t' << fn << '\t';
        for (int i = 0; i < results.size(); i++) { m->mothurOut(toString(results[i]) + "\t"); outStep << results[i] << "\t"; }
        m->mothurOutEndLine();
        outStep << endl;
        
        while ((delta > stableMetric) && (iters < maxIters)) {
            
            long start = time(NULL);
            
            if (m->getControl_pressed()) { break; }
            double oldMetric = listVectorMetric;
            
            cluster.update(listVectorMetric);

            delta = abs(oldMetric - listVectorMetric);
            iters++;
            
            results = cluster.getStats(tp, tn, fp, fn);
            numBins = cluster.getNumBins();
            
            m->mothurOut(toString(iters) + "\t" + toString(time(NULL) - start) + "\t" + toString(cutoff) + "\t" + toString(numBins) + "\t" + toString(cutoff) + "\t"+ toString(tp) + "\t" + toString(tn) + "\t" + toString(fp) + "\t" + toString(fn) + "\t");
            outStep << (toString(iters) + "\t" + toString(time(NULL) - start) + "\t" + toString(cutoff) + "\t" + toString(numBins) + "\t" + toString(cutoff) + "\t") << tp << '\t' << tn << '\t' << fp << '\t' << fn << '\t';
            for (int i = 0; i < results.size(); i++) { m->mothurOut(toString(results[i]) + "\t"); outStep << results[i] << "\t"; }
            m->mothurOutEndLine();
            outStep << endl;
        }
        m->mothurOutEndLine(); m->mothurOutEndLine();
        
        if (m->getControl_pressed()) { delete metric; metric = NULL; return 0; }
        
        ListVector* list = cluster.getList();
        list->setLabel(toString(cutoff));
        if(countfile != "") { list->print(listFile, counts); }
        else { list->print(listFile); }
        listFile.close();
        
        variables["[filename]"] = fileroot;
        variables["[clustertag]"] = tag;
        string sabundFileName = getOutputFileName("sabund", variables);
        string rabundFileName = getOutputFileName("rabund", variables);
        
        if (countfile == "") {
            util.openOutputFile(sabundFileName,	sabundFile);
            util.openOutputFile(rabundFileName,	rabundFile);
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
        
        string sensspecFilename = fileroot+ tag + ".sensspec";
        ofstream sensFile;
        util.openOutputFile(sensspecFilename,	sensFile);
        outputNames.push_back(sensspecFilename); outputTypes["sensspec"].push_back(sensspecFilename);
        
        
        sensFile << "label\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";
        
        results = cluster.getStats(tp, tn, fp, fn);
        
        sensFile << cutoff << '\t' << cutoff << '\t' << tp << '\t' << tn << '\t' << fp << '\t' << fn << '\t';
        for (int i = 0; i < results.size(); i++) {  sensFile << results[i] << '\t'; }
        sensFile << '\n';
        sensFile.close();
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterCommand", "runOptiCluster");
        exit(1);
    }
    
}
//**********************************************************************************************************************

int ClusterCommand::runUniqueCluster(){
    try {
        if (!cutOffSet) {  m->mothurOut("\nYou did not set a cutoff, using 0.03.\n"); cutoff = 0.03;  }
        
        if (countfile != "") {  distfile = countfile;   }
        else if (namefile != "") {  distfile = namefile;  }
        
        m->mothurOutEndLine(); m->mothurOut("Clustering " + distfile); m->mothurOutEndLine();
        
        ListVector list; list.setLabel(toString(cutoff));
        
        map<string, int> counts;
        if (countfile != "") {
            CountTable ct; ct.readTable(countfile, false, false); counts = ct.getNameMap();
            for (map<string, int>::iterator it = counts.begin(); it != counts.end(); it++) {
                if (m->getControl_pressed()) { return 0; }
                list.push_back(it->first);
            }
        }else {
            map<string, string> nameMap;
            util.readNames(namefile, nameMap);
            for (map<string, string>::iterator it = nameMap.begin(); it != nameMap.end(); it++) {
                if (m->getControl_pressed()) { return 0; }
                list.push_back(it->second);
            }
        }
        
        if (outputDir == "") { outputDir += util.hasPath(distfile); }
        fileroot = outputDir + util.getRootName(util.getSimpleName(distfile));
        
        tag = "unique";
        string listFileName = fileroot+ tag + ".list";
        
        ofstream listFile;
        util.openOutputFile(listFileName,	listFile);
        outputNames.push_back(listFileName); outputTypes["list"].push_back(listFileName);

        if(countfile != "") { list.print(listFile, counts); }
        else { list.print(listFile); }
        listFile.close();
        
        map<string, string> variables;
        variables["[filename]"] = fileroot;
        variables["[clustertag]"] = tag;
        string sabundFileName = getOutputFileName("sabund", variables);
        string rabundFileName = getOutputFileName("rabund", variables);
        
        if (countfile == "") {
            util.openOutputFile(sabundFileName,	sabundFile);
            util.openOutputFile(rabundFileName,	rabundFile);
            outputNames.push_back(sabundFileName); outputTypes["sabund"].push_back(sabundFileName);
            outputNames.push_back(rabundFileName); outputTypes["rabund"].push_back(rabundFileName);
            
            SAbundVector sabund = list.getSAbundVector();
            sabund.print(sabundFile);
            sabundFile.close();
            
            RAbundVector rabund = list.getRAbundVector();
            rabund.print(rabundFile);
            rabundFile.close();
        }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterCommand", "runUniqueCluster");
        exit(1);
    }
    
}
//**********************************************************************************************************************
