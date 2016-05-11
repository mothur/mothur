/*
 *  clustersplitcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 5/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "clustersplitcommand.h"
#include "vsearchfileparser.h"
#include "systemcommand.h"

//**********************************************************************************************************************
vector<string> ClusterSplitCommand::setParameters(){	
	try {
        CommandParameter pfile("file", "InputTypes", "", "", "PhylipColumnFasta", "PhylipColumnFasta", "none","",false,false,true); parameters.push_back(pfile);
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "none", "FastaTaxName","",false,false,true); parameters.push_back(ptaxonomy);
		CommandParameter pphylip("phylip", "InputTypes", "", "", "PhylipColumnFasta", "PhylipColumnFasta", "none","list",false,false,true); parameters.push_back(pphylip);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "PhylipColumnFasta", "PhylipColumnFasta", "FastaTaxName","list",false,false,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "ColumnName-FastaTaxName","rabund-sabund",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount", "none", "","",false,false,true); parameters.push_back(pcount);
		CommandParameter pcolumn("column", "InputTypes", "", "", "PhylipColumnFasta", "PhylipColumnFasta", "ColumnName","list",false,false,true); parameters.push_back(pcolumn);
		CommandParameter ptaxlevel("taxlevel", "Number", "", "3", "", "", "","",false,false,true); parameters.push_back(ptaxlevel);
		CommandParameter psplitmethod("splitmethod", "Multiple", "classify-fasta-distance", "distance", "", "", "","",false,false,true); parameters.push_back(psplitmethod);
		CommandParameter plarge("large", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(plarge);
		CommandParameter pshowabund("showabund", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pshowabund);
        CommandParameter pcluster("cluster", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pcluster);
		CommandParameter ptiming("timing", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(ptiming);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pcutoff("cutoff", "Number", "", "0.25", "", "", "","",false,false,true); parameters.push_back(pcutoff);
		CommandParameter pprecision("precision", "Number", "", "100", "", "", "","",false,false); parameters.push_back(pprecision);
        CommandParameter pmethod("method", "Multiple", "furthest-nearest-average-weighted-agc-dgc-opti", "average", "", "", "","",false,false,true); parameters.push_back(pmethod);
		CommandParameter phard("hard", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(phard);
        CommandParameter pislist("islist", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pislist);
        CommandParameter pclassic("classic", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pclassic);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
			
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ClusterSplitCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The cluster.split command parameter options are file, fasta, phylip, column, name, count, cutoff, precision, method, splitmethod, taxonomy, taxlevel, showabund, timing, hard, large, cluster, processors. Fasta or Phylip or column and name are required.\n";
		helpString += "The cluster.split command can split your files in 3 ways. Splitting by distance file, by classification, or by classification also using a fasta file. \n";
		helpString += "For the distance file method, you need only provide your distance file and mothur will split the file into distinct groups. \n";
		helpString += "For the classification method, you need to provide your distance file and taxonomy file, and set the splitmethod to classify.  \n";
		helpString += "You will also need to set the taxlevel you want to split by. mothur will split the sequences into distinct taxonomy groups, and split the distance file based on those groups. \n";
		helpString += "For the classification method using a fasta file, you need to provide your fasta file, names file and taxonomy file.  \n";
		helpString += "You will also need to set the taxlevel you want to split by. mothur will split the sequence into distinct taxonomy groups, and create distance files for each grouping. \n";
        helpString += "The file option allows you to enter your file containing your list of column and names/count files as well as the singleton file.  This file is mothur generated, when you run cluster.split() with the cluster=f parameter.  This can be helpful when you have a large dataset that you may be able to use all your processors for the splitting step, but have to reduce them for the cluster step due to RAM constraints. For example: cluster.split(fasta=yourFasta, taxonomy=yourTax, count=yourCount, taxlevel=3, cluster=f, processors=8) then cluster.split(file=yourFile, processors=4).  This allows your to maximize your processors during the splitting step.  Also, if you are unsure if the cluster step will have RAM issue with multiple processors, you can avoid running the first part of the command multiple times.\n";
		helpString += "The phylip and column parameter allow you to enter your distance file. \n";
		helpString += "The fasta parameter allows you to enter your aligned fasta file. \n";
		helpString += "The name parameter allows you to enter your name file. \n";
        helpString += "The count parameter allows you to enter your count file. \n A count or name file is required if your distance file is in column format";
        helpString += "The cluster parameter allows you to indicate whether you want to run the clustering or just split the distance matrix, default=t";
		helpString += "The cutoff parameter allow you to set the distance you want to cluster to, default is 0.25. \n";
		helpString += "The precision parameter allows you specify the precision of the precision of the distances outputted, default=100, meaning 2 decimal places. \n";
		helpString += "The method parameter allows you to enter your clustering mothod. Options are furthest, nearest, average, weighted, agc, dgc and opti. Default=average.  The agc and dgc methods require a fasta file.";
		helpString += "The splitmethod parameter allows you to specify how you want to split your distance file before you cluster, default=distance, options distance, classify or fasta. \n";
		helpString += "The taxonomy parameter allows you to enter the taxonomy file for your sequences, this is only valid if you are using splitmethod=classify. Be sure your taxonomy file does not include the probability scores. \n";
		helpString += "The taxlevel parameter allows you to specify the taxonomy level you want to use to split the distance file, default=3, meaning use the first taxon in each list. \n";
		helpString += "The large parameter allows you to indicate that your distance matrix is too large to fit in RAM.  The default value is false.\n";
        helpString += "The classic parameter allows you to indicate that you want to run your files with cluster.classic.  It is only valid with splitmethod=fasta. Default=f.\n";
		helpString += "The cluster.split command should be in the following format: \n";
		helpString += "cluster.split(column=youDistanceFile, name=yourNameFile, method=yourMethod, cutoff=yourCutoff, precision=yourPrecision, splitmethod=yourSplitmethod, taxonomy=yourTaxonomyfile, taxlevel=yourtaxlevel) \n";
		helpString += "Example: cluster.split(column=abrecovery.dist, name=abrecovery.names, method=furthest, cutoff=0.10, precision=1000, splitmethod=classify, taxonomy=abrecovery.silva.slv.taxonomy, taxlevel=5) \n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ClusterSplitCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "list") {  pattern = "[filename],[clustertag],list-[filename],[clustertag],[tag2],list"; } 
        else if (type == "rabund") {  pattern = "[filename],[clustertag],rabund"; } 
        else if (type == "sabund") {  pattern = "[filename],[clustertag],sabund"; }
        else if (type == "column") {  pattern = "[filename],dist"; }
        else if (type == "file")   {  pattern = "[filename],file"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterSplitCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ClusterSplitCommand::ClusterSplitCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["list"] = tempOutNames;
		outputTypes["rabund"] = tempOutNames;
		outputTypes["sabund"] = tempOutNames;
		outputTypes["column"] = tempOutNames;
        outputTypes["file"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "ClusterSplitCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
//This function checks to make sure the cluster command has no errors and then clusters based on the method chosen.
ClusterSplitCommand::ClusterSplitCommand(string option)  {
	try{
		abort = false; calledHelp = false;   
		format = "";
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter("cluster.split");
		
			//check to make sure all parameters are valid for command
			map<string,string>::iterator it;
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
			outputTypes["column"] = tempOutNames;
            outputTypes["file"] = tempOutNames;
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
				//if the user changes the input directory command factory will send this info to us in the output parameter 
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
				
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
				
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
                
                it = parameters.find("file");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["file"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			file = validParameter.validFile(parameters, "file", true);
			if (file == "not open") { file = ""; abort = true; }
			else if (file == "not found") { file = ""; }
            else { distfile = file; }
            
            phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not open") { abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; }
			else {  distfile = phylipfile;  format = "phylip"; 	m->setPhylipFile(phylipfile); }
			
			columnfile = validParameter.validFile(parameters, "column", true);
			if (columnfile == "not open") { abort = true; }	
			else if (columnfile == "not found") { columnfile = ""; }
			else {  distfile = columnfile; format = "column";	m->setColumnFile(columnfile); }
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; namefile = "";}	
			else if (namefile == "not found") { namefile = "";  }
			else { m->setNameFile(namefile); }
            
            countfile = validParameter.validFile(parameters, "count", true);
			if (countfile == "not open") { abort = true; countfile = "";}	
			else if (countfile == "not found") { countfile = "";  }
			else { m->setCountTableFile(countfile); }
			
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }	
			else if (fastafile == "not found") { fastafile = ""; }
			else { distfile = fastafile;  splitmethod = "fasta";  m->setFastaFile(fastafile); }
			
			taxFile = validParameter.validFile(parameters, "taxonomy", true);
			if (taxFile == "not open") { taxFile = ""; abort = true; }	
			else if (taxFile == "not found") { taxFile = ""; }
			else {  m->setTaxonomyFile(taxFile); if (splitmethod != "fasta") { splitmethod = "classify"; } }
			
			if ((phylipfile == "") && (columnfile == "") && (fastafile == "") && (file == "")) {
				//is there are current file available for either of these?
				//give priority to column, then phylip, then fasta
				columnfile = m->getColumnFile(); 
				if (columnfile != "") {  format = "column"; m->mothurOut("Using " + columnfile + " as input file for the column parameter."); m->mothurOutEndLine(); }
				else { 
					phylipfile = m->getPhylipFile(); 
					if (phylipfile != "") {  format = "phylip"; m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter."); m->mothurOutEndLine(); }
					else { 
						fastafile = m->getFastaFile(); 
						if (fastafile != "") {   m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
						else { 
							m->mothurOut("No valid current files. When executing a cluster.split command you must enter a file, phylip or a column or fastafile."); m->mothurOutEndLine();
							abort = true; 
						}
					}
				}
			}
			else if ((phylipfile != "") && (columnfile != "") && (fastafile != "") && (file != "")) { m->mothurOut("When executing a cluster.split command you must enter ONLY ONE of the following: file, fasta, phylip or column."); m->mothurOutEndLine(); abort = true; }
            
            if ((countfile != "") && (namefile != "")) { m->mothurOut("When executing a cluster.split command you must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
            
			if (columnfile != "") {
				if ((namefile == "") && (countfile == "")) { 
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
			
			if (fastafile != "") {
				if (taxFile == "") { 
					taxFile = m->getTaxonomyFile(); 
					if (taxFile != "") {  m->mothurOut("Using " + taxFile + " as input file for the taxonomy parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("You need to provide a taxonomy file if you are if you are using a fasta file to generate the split."); m->mothurOutEndLine(); 
						abort = true; 
					}	
				}
				
				if ((namefile == "") && (countfile == "")) { 
					namefile = m->getNameFile(); 
					if (namefile != "") {  m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
					else { 
						countfile = m->getCountTableFile();
                        if (countfile != "") {  m->mothurOut("Using " + countfile + " as input file for the count parameter."); m->mothurOutEndLine(); }
                        else { 
                            m->mothurOut("You need to provide a namefile or countfile if you are going to use the fasta file to generate the split."); m->mothurOutEndLine(); 
                            abort = true; 
                        }	
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
			m->mothurConvert(temp, precision); 
			
			temp = validParameter.validFile(parameters, "hard", false);			if (temp == "not found") { temp = "T"; }
			hard = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "large", false);			if (temp == "not found") { temp = "F"; }
			large = m->isTrue(temp);
            
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
			
			temp = validParameter.validFile(parameters, "splitmethod", false);	
			if ((splitmethod != "fasta") && (splitmethod != "classify")) {
				if (temp == "not found")  { splitmethod = "distance"; }
				else {  splitmethod = temp; }
			}
			
            temp = validParameter.validFile(parameters, "classic", false);			if (temp == "not found") { temp = "F"; }
			classic = m->isTrue(temp);
            
            //not using file option and don't have fasta method with classic
            if (((splitmethod != "fasta") && classic) && (file == "")) { m->mothurOut("[ERROR]: splitmethod must be fasta to use cluster.classic, or you must use the file option.\n"); abort=true; }
            
            cutoffNotSet = false;
            temp = validParameter.validFile(parameters, "cutoff", false);		if (temp == "not found")  { cutoffNotSet = true; temp = "0.25"; }
			m->mothurConvert(temp, cutoff); 
			cutoff += (5 / (precision * 10.0));  
			
			temp = validParameter.validFile(parameters, "taxlevel", false);		if (temp == "not found")  { temp = "3"; }
			m->mothurConvert(temp, taxLevelCutoff); 
			
			method = validParameter.validFile(parameters, "method", false);		if (method == "not found") { method = "average"; }
			
            if ((method == "furthest") || (method == "nearest") || (method == "average") || (method == "weighted") || (method == "agc") || (method == "dgc") || (method == "opti")) { }
            else { m->mothurOut("[ERROR]: Not a valid clustering method.  Valid clustering algorithms are furthest, nearest, average, weighted, agc, dgc and opti."); m->mothurOutEndLine(); abort = true; }
            
            if ((method == "agc") || (method == "dgc")) {
                if (fastafile == "") { m->mothurOut("[ERROR]: You must provide a fasta file when using the agc or dgc clustering methods, aborting\n."); abort = true;}
                if (classic) { m->mothurOut("[ERROR]: You cannot use cluster.classic with the agc or dgc clustering methods, aborting\n."); abort = true; }
            }
            #if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
            #else
            if ((method == "agc") || (method == "dgc")) { m->mothurOut("[ERROR]: The agc and dgc clustering methods are not available for Windows, aborting\n."); abort = true; }
            #endif
            
			if ((splitmethod == "distance") || (splitmethod == "classify") || (splitmethod == "fasta")) { }
			else { m->mothurOut("[ERROR]: " + splitmethod + " is not a valid splitting method.  Valid splitting algorithms are distance, classify or fasta."); m->mothurOutEndLine(); abort = true; }
			
			if ((splitmethod == "classify") && (taxFile == "")) {  m->mothurOut("[ERROR]: You need to provide a taxonomy file if you are going to use the classify splitmethod."); m->mothurOutEndLine(); abort = true;  }

			showabund = validParameter.validFile(parameters, "showabund", false);
			if (showabund == "not found") { showabund = "T"; }
            
            temp = validParameter.validFile(parameters, "cluster", false);  if (temp == "not found") { temp = "T"; }
            runCluster = m->isTrue(temp);
            
            temp = validParameter.validFile(parameters, "islist", false);  if (temp == "not found") { temp = "F"; }
            isList = m->isTrue(temp);

			timing = validParameter.validFile(parameters, "timing", false);
			if (timing == "not found") { timing = "F"; }
			
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "ClusterSplitCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

int ClusterSplitCommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        time_t estart;
        vector<string> listFileNames;
        vector< map<string, string> > distName;
        set<string> labels;
        string singletonName = "";
        
        double saveCutoff = cutoff;

        if (file != "") {
            deleteFiles = false; estart = time(NULL);
            singletonName = readFile(distName);

            if (isList) {
                
                //set list file as new current listfile
                string current = "";
                itTypes = outputTypes.find("list");
                if (itTypes != outputTypes.end()) {
                    if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setListFile(current); }
                }
                
                m->mothurOutEndLine();
                m->mothurOut("Output File Names: "); m->mothurOutEndLine();
                for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
                m->mothurOutEndLine();
                
                return 0;
            }
                    
        }else {
		          
            //****************** file prep work ******************************//
            
                //if user gave a phylip file convert to column file
                if (format == "phylip") {
                    estart = time(NULL);
                    m->mothurOut("Converting to column format..."); m->mothurOutEndLine();
                    
                    ReadCluster* convert = new ReadCluster(distfile, cutoff, outputDir, false);
                    
                    NameAssignment* nameMap = NULL;
                    convert->setFormat("phylip");
                    convert->read(nameMap);
                    
                    if (m->control_pressed) {  delete convert;  return 0;  }
                    
                    distfile = convert->getOutputFile();
                    
                    //if no names file given with phylip file, create it
                    ListVector* listToMakeNameFile =  convert->getListVector();
                    if ((namefile == "") && (countfile == "")) {  //you need to make a namefile for split matrix
                        ofstream out;
                        namefile = phylipfile + ".names";
                        m->openOutputFile(namefile, out);
                        for (int i = 0; i < listToMakeNameFile->getNumBins(); i++) {
                            string bin = listToMakeNameFile->get(i);
                            out << bin << '\t' << bin << endl;
                        }
                        out.close();
                    }
                    delete listToMakeNameFile;
                    delete convert;
                    
                    m->mothurOut("It took " + toString(time(NULL) - estart) + " seconds to convert the distance file."); m->mothurOutEndLine();
                }
                if (m->control_pressed) { return 0; }
                
                estart = time(NULL);
                m->mothurOut("Splitting the file..."); m->mothurOutEndLine();

                //split matrix into non-overlapping groups
                SplitMatrix* split;
                if (splitmethod == "distance")			{	split = new SplitMatrix(distfile, namefile, countfile, taxFile, cutoff, splitmethod, large);							}
                else if (splitmethod == "classify")		{	split = new SplitMatrix(distfile, namefile, countfile, taxFile, taxLevelCutoff, splitmethod, large);					}
                else if (splitmethod == "fasta")		{
                    if ((method == "agc") || (method == "dgc")) {
                        if (!findVsearch()) { m->mothurOut("[ERROR] cannot find vsearch executable, aborting.\n"); return 0; }
                        split = new SplitMatrix(fastafile, namefile, countfile, taxFile, taxLevelCutoff, cutoff, "vsearch", processors, classic, outputDir, "fasta");
                    }else{
                        split = new SplitMatrix(fastafile, namefile, countfile, taxFile, taxLevelCutoff, cutoff, splitmethod, processors, classic, outputDir, "distance");
                    }
                }
                else { m->mothurOut("Not a valid splitting method.  Valid splitting algorithms are distance, classify or fasta."); m->mothurOutEndLine(); return 0;		}
                split->split();

                if (m->control_pressed) { delete split; return 0; }
                
                singletonName = split->getSingletonNames();
                distName = split->getDistanceFiles();  //returns map of distance files -> namefile sorted by distance file size
                delete split;
                
                if (m->debug) { m->mothurOut("[DEBUG]: distName.size() = " + toString(distName.size()) + ".\n"); }
                
                //output a merged distance file
                //if (splitmethod == "fasta")		{ createMergedDistanceFile(distName); }
				
                if (m->control_pressed) { return 0; }
                
                m->mothurOut("It took " + toString(time(NULL) - estart) + " seconds to split the distance file."); m->mothurOutEndLine();
                estart = time(NULL);

                if (!runCluster) {
                    string filename = printFile(singletonName, distName);
                    
                    m->mothurOutEndLine();
                    m->mothurOut("Output File Names: "); m->mothurOutEndLine();
                    m->mothurOutEndLine(); m->mothurOut(filename); m->mothurOutEndLine();
                    for (int i = 0; i < distName.size(); i++) {	m->mothurOut(distName[i].begin()->first); m->mothurOutEndLine(); m->mothurOut(distName[i].begin()->second); m->mothurOutEndLine();	}
                    m->mothurOutEndLine();

                    return 0;
                }
                deleteFiles = true;

            }
		//****************** break up files between processes and cluster each file set ******************************//
		///////////////////// WINDOWS CAN ONLY USE 1 PROCESSORS ACCESS VIOLATION UNRESOLVED ///////////////////////
		//sanity check
		if (processors > distName.size()) { processors = distName.size(); }
		
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				if(processors == 1){
					listFileNames = cluster(distName, labels); //clusters individual files and returns names of list files
				}else{
					listFileNames = createProcesses(distName, labels);
                }
		#else
				listFileNames = cluster(distName, labels); //clusters individual files and returns names of list files
		#endif
		
		if (m->control_pressed) { for (int i = 0; i < listFileNames.size(); i++) { m->mothurRemove(listFileNames[i]); } return 0; }
		
		if (saveCutoff != cutoff) { m->mothurOut("\nCutoff was " + toString(saveCutoff) + " changed cutoff to " + toString(cutoff)); m->mothurOutEndLine();  }
		
		m->mothurOut("It took " + toString(time(NULL) - estart) + " seconds to cluster"); m->mothurOutEndLine();
		
		//****************** merge list file and create rabund and sabund files ******************************//
		estart = time(NULL);
		m->mothurOut("Merging the clustered files..."); m->mothurOutEndLine();

		ListVector* listSingle;
		map<float, int> labelBins = completeListFile(listFileNames, singletonName, labels, listSingle); //returns map of label to numBins
		
		if (m->control_pressed) { if (listSingle != NULL) { delete listSingle; } for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); } return 0; }
		
		mergeLists(listFileNames, labelBins, listSingle);

		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); } return 0; }
        
        //delete after all are complete incase a crash happens
        if (!deleteFiles) { for (int i = 0; i < distName.size(); i++) {	m->mothurRemove(distName[i].begin()->first); m->mothurRemove(distName[i].begin()->second); 	} }
		
		m->mothurOut("It took " + toString(time(NULL) - estart) + " seconds to merge."); m->mothurOutEndLine();
		
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
		m->errorOut(e, "ClusterSplitCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
map<float, int> ClusterSplitCommand::completeListFile(vector<string> listNames, string singleton, set<string>& userLabels, ListVector*& listSingle){
	try {
		map<float, int> labelBin;
		vector<float> orderFloat;
		int numSingleBins;
		
		//read in singletons
		if (singleton != "none") {
            
            ifstream in;
            m->openInputFile(singleton, in);
				
			string firstCol, secondCol;
			listSingle = new ListVector();
            
            if (countfile != "") { m->getline(in); m->gobble(in); }
            
			while (!in.eof()) {
				in >> firstCol >> secondCol; m->getline(in); m->gobble(in);
				if (countfile == "") { listSingle->push_back(secondCol); }
                else { listSingle->push_back(firstCol); }
			}
            
			in.close();
			m->mothurRemove(singleton);
			
			numSingleBins = listSingle->getNumBins();
		}else{  listSingle = NULL; numSingleBins = 0;  }
		
		//go through users set and make them floats so we can sort them 
		for(set<string>::iterator it = userLabels.begin(); it != userLabels.end(); ++it) {
			float temp = -10.0;

			if ((*it != "unique") && (convertTestFloat(*it, temp) == true))	{	convert(*it, temp);	}
			else if (*it == "unique")										{	temp = -1.0;		}
			
			if (temp <= cutoff) {
				orderFloat.push_back(temp);
				labelBin[temp] = numSingleBins; //initialize numbins 
			}
		}
	
		//sort order
		sort(orderFloat.begin(), orderFloat.end());
		userLabels.clear();
			
		//get the list info from each file
		for (int k = 0; k < listNames.size(); k++) {
            
			if (m->control_pressed) {  
				if (listSingle != NULL) { delete listSingle; listSingle = NULL; m->mothurRemove(singleton);  }
				for (int i = 0; i < listNames.size(); i++) {   m->mothurRemove(listNames[i]);  }
				return labelBin;
			}
			
			InputData* input = new InputData(listNames[k], "list");
			ListVector* list = input->getListVector();
			string lastLabel = list->getLabel();
            
			string filledInList = listNames[k] + "filledInTemp";
			ofstream outFilled;
			m->openOutputFile(filledInList, outFilled);
	
			//for each label needed
			for(int l = 0; l < orderFloat.size(); l++){
			
				string thisLabel;
				if (orderFloat[l] == -1) { thisLabel = "unique"; }
				else { thisLabel = toString(orderFloat[l],  length-1);  } 

				//this file has reached the end
				if (list == NULL) { 
					list = input->getListVector(lastLabel, true); 
				}else{	//do you have the distance, or do you need to fill in
						
					float labelFloat;
					if (list->getLabel() == "unique") {  labelFloat = -1.0;  }
					else { convert(list->getLabel(), labelFloat); }

					//check for missing labels
					if (labelFloat > orderFloat[l]) { //you are missing the label, get the next smallest one
						//if its bigger get last label, otherwise keep it
						delete list;
						list = input->getListVector(lastLabel, true);  //get last list vector to use, you actually want to move back in the file
					}
					lastLabel = list->getLabel();
				}
				
				//print to new file
				list->setLabel(thisLabel);
				list->print(outFilled, true);
		
				//update labelBin
				labelBin[orderFloat[l]] += list->getNumBins();
									
				delete list;
									
				list = input->getListVector();
			}
			
			if (list != NULL) { delete list; }
			delete input;
			
			outFilled.close();
			m->mothurRemove(listNames[k]);
			rename(filledInList.c_str(), listNames[k].c_str());
		}
		
		return labelBin;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "completeListFile");
		exit(1);
	}
}
//**********************************************************************************************************************
int ClusterSplitCommand::mergeLists(vector<string> listNames, map<float, int> userLabels, ListVector* listSingle){
	try {
		if (outputDir == "") { outputDir += m->hasPath(distfile); }
		fileroot = outputDir + m->getRootName(m->getSimpleName(distfile));
		
        map<string, string> variables; 
        variables["[filename]"] = fileroot;
        variables["[clustertag]"] = tag;
        string sabundFileName = getOutputFileName("sabund", variables);
        string rabundFileName = getOutputFileName("rabund", variables);
        if (countfile != "") { variables["[tag2]"] = "unique_list"; }
        string listFileName = getOutputFileName("list", variables);
        
        map<string, int> counts;
        if (countfile == "") {
            m->openOutputFile(sabundFileName,	outSabund);
            m->openOutputFile(rabundFileName,	outRabund);
            outputNames.push_back(sabundFileName); outputTypes["sabund"].push_back(sabundFileName);
            outputNames.push_back(rabundFileName); outputTypes["rabund"].push_back(rabundFileName);
            
        }else {
            if (file == "") {
                CountTable ct;
                ct.readTable(countfile, false, false);
                counts = ct.getNameMap();
            }
        }
        
		m->openOutputFile(listFileName,	outList);
        outputNames.push_back(listFileName); outputTypes["list"].push_back(listFileName);
		
		map<float, int>::iterator itLabel;
        
        //clears out junk for autocompleting of list files above.  Perhaps there is a beter way to handle this from within the data structure?
        m->printedListHeaders = false;

		//for each label needed
		for(itLabel = userLabels.begin(); itLabel != userLabels.end(); itLabel++) {
			
			string thisLabel;
			if (itLabel->first == -1) { thisLabel = "unique"; }
			else { thisLabel = toString(itLabel->first,  length-1);  } 
			
			//outList << thisLabel << '\t' << itLabel->second << '\t';
            
            RAbundVector* rabund = NULL;
            ListVector completeList;
            completeList.setLabel(thisLabel);
            
            if (countfile == "") {
                rabund = new RAbundVector();
                rabund->setLabel(thisLabel);
            }

			//add in singletons
			if (listSingle != NULL) {
				for (int j = 0; j < listSingle->getNumBins(); j++) {
					//outList << listSingle->get(j) << '\t';
                    completeList.push_back(listSingle->get(j));
					if (countfile == "") { rabund->push_back(m->getNumNames(listSingle->get(j))); }
				}
			}
			
			//get the list info from each file
			for (int k = 0; k < listNames.size(); k++) {
	
				if (m->control_pressed) {  if (listSingle != NULL) { delete listSingle;   } for (int i = 0; i < listNames.size(); i++) { m->mothurRemove(listNames[i]);  } if (rabund != NULL) { delete rabund; } return 0; }
				
				InputData* input = new InputData(listNames[k], "list");
				ListVector* list = input->getListVector(thisLabel);
				
				//this file has reached the end
				if (list == NULL) { m->mothurOut("Error merging listvectors in file " + listNames[k]); m->mothurOutEndLine();  }	
				else {		
					for (int j = 0; j < list->getNumBins(); j++) {
						//outList << list->get(j) << '\t';
                        completeList.push_back(list->get(j));
						if (countfile == "") { rabund->push_back(m->getNumNames(list->get(j))); }
					}
					delete list;
				}
				delete input;
			}
			
            if (countfile == "") {
                SAbundVector sabund = rabund->getSAbundVector();
                sabund.print(outSabund);
                rabund->print(outRabund);
            }
			//outList << endl;
            if (!m->printedListHeaders) { m->listBinLabelsInFile.clear(); completeList.printHeaders(outList); }
            if (countfile == "") { completeList.print(outList); }
            else if ((file == "") && (countfile != "")) { completeList.print(outList, counts);   }
            else { completeList.print(outList); }
			
			if (rabund != NULL) { delete rabund; }
		}
		
		outList.close();
        if (countfile == "") {
            outRabund.close();
            outSabund.close();
		}
		if (listSingle != NULL) { delete listSingle;  }
		
		for (int i = 0; i < listNames.size(); i++) {  m->mothurRemove(listNames[i]);  }
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "mergeLists");
		exit(1);
	}
}

//**********************************************************************************************************************

void ClusterSplitCommand::printData(ListVector* oldList){
	try {
		string label = oldList->getLabel();
		RAbundVector oldRAbund = oldList->getRAbundVector();
		
		oldRAbund.setLabel(label);
		if (m->isTrue(showabund)) {
			oldRAbund.getSAbundVector().print(cout);
		}
		oldRAbund.print(outRabund);
		oldRAbund.getSAbundVector().print(outSabund);
	
		oldList->print(outList, true);
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "printData");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string>  ClusterSplitCommand::createProcesses(vector< map<string, string> > distName, set<string>& labels){
	try {
        deleteFiles = false; //so if we need to recalc the processors the files are still there
        bool recalc = false;
        vector<string> listFiles;
        vector < vector < map<string, string> > > dividedNames; //distNames[1] = vector of filenames for process 1...
        dividedNames.resize(processors);
        
        //for each file group figure out which process will complete it
        //want to divide the load intelligently so the big files are spread between processes
        for (int i = 0; i < distName.size(); i++) { 
            //cout << i << endl;
            int processToAssign = (i+1) % processors; 
            if (processToAssign == 0) { processToAssign = processors; }
            
            dividedNames[(processToAssign-1)].push_back(distName[i]);
            if ((processToAssign-1) == 1) { m->mothurOut(distName[i].begin()->first + "\n"); }
        }
        
        //now lets reverse the order of ever other process, so we balance big files running with little ones
        for (int i = 0; i < processors; i++) {
            //cout << i << endl;
            int remainder = ((i+1) % processors);
            if (remainder) {  reverse(dividedNames[i].begin(), dividedNames[i].end());  }
        }
        
        if (m->control_pressed) { return listFiles; }
	
	#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		int process = 1;
		processIDS.clear();
		
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				set<string> labels;
				vector<string> listFileNames = cluster(dividedNames[process], labels);
				
				//write out names to file
				string filename = m->mothurGetpid(process) + ".temp";
				ofstream out;
				m->openOutputFile(filename, out);
				out << tag << endl;
				for (int j = 0; j < listFileNames.size(); j++) { out << listFileNames[j] << endl;  }
				out.close();
				
				//print out labels
				ofstream outLabels;
				filename = m->mothurGetpid(process) + ".temp.labels";
				m->openOutputFile(filename, outLabels);
				
				outLabels << cutoff << endl;
				for (set<string>::iterator it = labels.begin(); it != labels.end(); it++) {
					outLabels << (*it) << endl;
				}
				outLabels.close();

				exit(0);
            }else {
                m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(process) + "\n"); processors = process;
                for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                //wait to die
                for (int i=0;i<processIDS.size();i++) {
                    int temp = processIDS[i];
                    wait(&temp);
                }
                for (int i=0;i<processIDS.size();i++) {
                    m->mothurRemove((toString(processIDS[i]) + ".temp"));
                    m->mothurRemove((toString(processIDS[i]) + ".temp.labels"));
                }
                m->control_pressed = false;
                recalc = true;
                break;
            }

		}
        
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->control_pressed = false;  for (int i=0;i<processIDS.size();i++) {m->mothurRemove((toString(processIDS[i]) + ".temp"));m->mothurRemove((toString(processIDS[i]) + ".temp.labels"));} processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            
            listFiles.clear();
            dividedNames.clear(); //distNames[1] = vector of filenames for process 1...
            dividedNames.resize(processors);
            
            //for each file group figure out which process will complete it
            //want to divide the load intelligently so the big files are spread between processes
            for (int i = 0; i < distName.size(); i++) {
                //cout << i << endl;
                int processToAssign = (i+1) % processors;
                if (processToAssign == 0) { processToAssign = processors; }
                
                dividedNames[(processToAssign-1)].push_back(distName[i]);
                if ((processToAssign-1) == 1) { m->mothurOut(distName[i].begin()->first + "\n"); }
            }
            
            //now lets reverse the order of ever other process, so we balance big files running with little ones
            for (int i = 0; i < processors; i++) {
                //cout << i << endl;
                int remainder = ((i+1) % processors);
                if (remainder) {  reverse(dividedNames[i].begin(), dividedNames[i].end());  }
            }
            
            processIDS.resize(0);
            process = 1;
            
            while (process != processors) {
                pid_t pid = fork();
                
                if (pid > 0) {
                    processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                    process++;
                }else if (pid == 0){
                    set<string> labels;
                    vector<string> listFileNames = cluster(dividedNames[process], labels);
                    
                    //write out names to file
                    string filename = m->mothurGetpid(process) + ".temp";
                    ofstream out;
                    m->openOutputFile(filename, out);
                    out << tag << endl;
                    for (int j = 0; j < listFileNames.size(); j++) { out << listFileNames[j] << endl;  }
                    out.close();
                    
                    //print out labels
                    ofstream outLabels;
                    filename = m->mothurGetpid(process) + ".temp.labels";
                    m->openOutputFile(filename, outLabels);
                    
                    outLabels << cutoff << endl;
                    for (set<string>::iterator it = labels.begin(); it != labels.end(); it++) {
                        outLabels << (*it) << endl;
                    }
                    outLabels.close();
                    
                    exit(0);
                }else {
                    m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine();
                    for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                    exit(0);
                }
            }
        }

		
        //do your part
        listFiles = cluster(dividedNames[0], labels);
        
		//force parent to wait until all the processes are done
		for (int i=0;i< processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
        
        //get list of list file names from each process
        for(int i=0;i<processIDS.size();i++){
            string filename = toString(processIDS[i]) + ".temp";
            ifstream in;
            m->openInputFile(filename, in);
            
            in >> tag; m->gobble(in);
            
            while(!in.eof()) {
                string tempName;
                in >> tempName; m->gobble(in);
                listFiles.push_back(tempName);
            }
            in.close();
            m->mothurRemove((toString(processIDS[i]) + ".temp"));
            
            //get labels
            filename = toString(processIDS[i]) + ".temp.labels";
            ifstream in2;
            m->openInputFile(filename, in2);
            
            float tempCutoff;
            in2 >> tempCutoff; m->gobble(in2);
            if (tempCutoff < cutoff) { cutoff = tempCutoff; }
            
            while(!in2.eof()) {
                string tempName;
                in2 >> tempName; m->gobble(in2);
                if (labels.count(tempName) == 0) { labels.insert(tempName); }
            }
            in2.close();
            m->mothurRemove((toString(processIDS[i]) + ".temp.labels"));
        }
        
        deleteFiles = true;
        
        //delete the temp files now that we are done
        for (int i = 0; i < distName.size(); i++) {
            string thisNamefile = distName[i].begin()->second;
            string thisDistFile = distName[i].begin()->first;
            m->mothurRemove(thisNamefile);
            m->mothurRemove(thisDistFile);
        }

    #else
    #endif
        
        return listFiles;
	
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "createProcesses");
		exit(1);
	}
}
//**********************************************************************************************************************

vector<string> ClusterSplitCommand::cluster(vector< map<string, string> > distNames, set<string>& labels){
	try {
		vector<string> listFileNames;
		double smallestCutoff = cutoff;
		
		//cluster each distance file
		for (int i = 0; i < distNames.size(); i++) {
            
			string thisNamefile = distNames[i].begin()->second;
			string thisDistFile = distNames[i].begin()->first;
			
			string listFileName = "";
            if (classic)    {  listFileName = clusterClassicFile(thisDistFile, thisNamefile, labels, smallestCutoff);   }
            else            {  listFileName = clusterFile(thisDistFile, thisNamefile, labels, smallestCutoff);          }

			if (m->control_pressed) { //clean up
				for (int i = 0; i < listFileNames.size(); i++) {	m->mothurRemove(listFileNames[i]); 	}
				listFileNames.clear(); return listFileNames;
			}
            
            listFileNames.push_back(listFileName);
        }
		
		cutoff = smallestCutoff;
        
		return listFileNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "cluster");
		exit(1);
	}


}
//**********************************************************************************************************************
string ClusterSplitCommand::clusterClassicFile(string thisDistFile, string thisNamefile, set<string>& labels, double& smallestCutoff){
	try {
        string listFileName = "";
        
        ListVector* list = NULL;
        ListVector oldList;
        RAbundVector* rabund = NULL;

        m->mothurOutEndLine(); m->mothurOut("Reading " + thisDistFile); m->mothurOutEndLine();
        
        //reads phylip file storing data in 2D vector, also fills list and rabund
        bool sim = false;
		ClusterClassic* cluster = new ClusterClassic(cutoff, method, sim);
        
        NameAssignment* nameMap = NULL;
        CountTable* ct = NULL;
        if(namefile != ""){	
			nameMap = new NameAssignment(thisNamefile);
			nameMap->readMap();
            cluster->readPhylipFile(thisDistFile, nameMap);
		}else if (countfile != "") {
            ct = new CountTable();
            ct->readTable(thisNamefile, false, false);
            cluster->readPhylipFile(thisDistFile, ct);
        }
        tag = cluster->getTag();
        
		if (m->control_pressed) { if(namefile != ""){	delete nameMap; }
            else { delete ct; } delete cluster; return 0; }
		
		list = cluster->getListVector();
		rabund = cluster->getRAbundVector();
        
		if (outputDir == "") { outputDir += m->hasPath(thisDistFile); }
		fileroot = outputDir + m->getRootName(m->getSimpleName(thisDistFile));
        listFileName = fileroot+ tag + ".list";
        
        ofstream listFile;
		m->openOutputFile(fileroot+ tag + ".list",	listFile);
		
		float previousDist = 0.00000;
		float rndPreviousDist = 0.00000;
		oldList = *list;

        m->mothurOutEndLine(); m->mothurOut("Clustering " + thisDistFile); m->mothurOutEndLine();
        
		while ((cluster->getSmallDist() < cutoff) && (cluster->getNSeqs() > 1)){
			if (m->control_pressed) { delete cluster; delete list; delete rabund; listFile.close();  if(namefile != ""){	delete nameMap; }
                else { delete ct; } return listFileName;  }
            
			cluster->update(cutoff);
            
			float dist = cluster->getSmallDist();
			float rndDist;
			if (hard) {
				rndDist = m->ceilDist(dist, precision); 
			}else{
				rndDist = m->roundDist(dist, precision); 
			}
            
            if(previousDist <= 0.0000 && dist != previousDist){
                oldList.setLabel("unique");
                oldList.print(listFile);
                if (labels.count("unique") == 0) {  labels.insert("unique");  }
            }
            else if(rndDist != rndPreviousDist){
                oldList.setLabel(toString(rndPreviousDist,  length-1));
                oldList.print(listFile);
                if (labels.count(toString(rndPreviousDist,  length-1)) == 0) { labels.insert(toString(rndPreviousDist,  length-1)); }
            }

            
			previousDist = dist;
			rndPreviousDist = rndDist;
			oldList = *list;
		}
        
		if(previousDist <= 0.0000){
            oldList.setLabel("unique");
            oldList.print(listFile);
            if (labels.count("unique") == 0) { labels.insert("unique"); }
        }
        else if(rndPreviousDist<cutoff){
            oldList.setLabel(toString(rndPreviousDist,  length-1));
            oldList.print(listFile);
            if (labels.count(toString(rndPreviousDist,  length-1)) == 0) { labels.insert(toString(rndPreviousDist,  length-1)); }
        }

        
		listFile.close();
		
		delete cluster;  delete list; delete rabund;
        if(namefile != ""){	delete nameMap; }
        else { delete ct; }
        
        if (deleteFiles) {
            m->mothurRemove(thisDistFile);
            m->mothurRemove(thisNamefile);
        }
        return listFileName;
        
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "clusterClassicFile");
		exit(1);
	}
}

//**********************************************************************************************************************
string ClusterSplitCommand::clusterFile(string thisDistFile, string thisNamefile, set<string>& labels, double& smallestCutoff){
	try {
        string listFileName = "";
        
        if ((method == "agc") || (method == "dgc")) {  listFileName = runVsearchCluster(thisDistFile, thisNamefile, labels, smallestCutoff);  }
        else if (method == "opti")                  {  listFileName = runOptiCluster(thisDistFile, thisNamefile, labels, smallestCutoff);     }
        else {
            
            Cluster* cluster = NULL;
            SparseDistanceMatrix* matrix = NULL;
            ListVector* list = NULL;
            ListVector oldList;
            RAbundVector* rabund = NULL;
            
            if (m->control_pressed) { return listFileName; }
            
            m->mothurOutEndLine(); m->mothurOut("Reading " + thisDistFile); m->mothurOutEndLine();
            
            ReadMatrix* read = new ReadColumnMatrix(thisDistFile);
            read->setCutoff(cutoff);
            
            NameAssignment* nameMap = NULL;
            CountTable* ct = NULL;
            if(namefile != ""){
                nameMap = new NameAssignment(thisNamefile);
                nameMap->readMap();
                read->read(nameMap);
            }else if (countfile != "") {
                ct = new CountTable();
                ct->readTable(thisNamefile, false, false);
                read->read(ct);
            }else { read->read(nameMap); }
            
            list = read->getListVector();
            oldList = *list;
            matrix = read->getDMatrix();
            
            if(countfile != "") {
                rabund = new RAbundVector();
                createRabund(ct, list, rabund); //creates an rabund that includes the counts for the unique list
                delete ct;
            }else { rabund = new RAbundVector(list->getRAbundVector()); }
            
            delete read;  read = NULL;
            if (namefile != "") { delete nameMap; nameMap = NULL; }
            
            m->mothurOutEndLine(); m->mothurOut("Clustering " + thisDistFile); m->mothurOutEndLine();
            
            //create cluster
            float adjust = -1.0;
            if (method == "furthest")	{	cluster = new CompleteLinkage(rabund, list, matrix, cutoff, method, adjust); }
            else if(method == "nearest"){	cluster = new SingleLinkage(rabund, list, matrix, cutoff, method, adjust); }
            else if(method == "average"){	cluster = new AverageLinkage(rabund, list, matrix, cutoff, method, adjust);	}
            tag = cluster->getTag();
            
            if (outputDir == "") { outputDir += m->hasPath(thisDistFile); }
            fileroot = outputDir + m->getRootName(m->getSimpleName(thisDistFile));
            
            ofstream listFile;
            m->openOutputFile(fileroot+ tag + ".list",	listFile);
            
            listFileName = fileroot+ tag + ".list";
            
            float previousDist = 0.00000;
            float rndPreviousDist = 0.00000;
            
            oldList = *list;
            
            print_start = true;
            start = time(NULL);
            double saveCutoff = cutoff;
            
            while (matrix->getSmallDist() < cutoff && matrix->getNNodes() > 0){
                
                if (m->control_pressed) { //clean up
                    delete matrix; delete list;	delete cluster; delete rabund;
                    listFile.close();
                    m->mothurRemove(listFileName);
                    return listFileName;
                }
                
                cluster->update(saveCutoff);
                
                float dist = matrix->getSmallDist();
                float rndDist;
                if (hard) {
                    rndDist = m->ceilDist(dist, precision);
                }else{
                    rndDist = m->roundDist(dist, precision);
                }
                
                if(previousDist <= 0.0000 && dist != previousDist){
                    oldList.setLabel("unique");
                    oldList.print(listFile);
                    if (labels.count("unique") == 0) {  labels.insert("unique");  }
                }
                else if(rndDist != rndPreviousDist){
                    oldList.setLabel(toString(rndPreviousDist,  length-1));
                    oldList.print(listFile);
                    if (labels.count(toString(rndPreviousDist,  length-1)) == 0) { labels.insert(toString(rndPreviousDist,  length-1)); }
                }
                
                previousDist = dist;
                rndPreviousDist = rndDist;
                oldList = *list;
            }
            
            
            if(previousDist <= 0.0000){
                oldList.setLabel("unique");
                oldList.print(listFile);
                if (labels.count("unique") == 0) { labels.insert("unique"); }
            }
            else if(rndPreviousDist<cutoff){
                oldList.setLabel(toString(rndPreviousDist,  length-1));
                oldList.print(listFile);
                if (labels.count(toString(rndPreviousDist,  length-1)) == 0) { labels.insert(toString(rndPreviousDist,  length-1)); }
            }
            
            delete matrix; delete list;	delete cluster; delete rabund;
            matrix = NULL; list = NULL; cluster = NULL; rabund = NULL;
            listFile.close();
            
            if (m->control_pressed) { //clean up
                m->mothurRemove(listFileName);
                return listFileName;
            }
            
            if (deleteFiles) {
                m->mothurRemove(thisDistFile);
                m->mothurRemove(thisNamefile);
            }
            
            if (saveCutoff != cutoff) { 
                if (hard)	{  saveCutoff = m->ceilDist(saveCutoff, precision);	}
                else		{	saveCutoff = m->roundDist(saveCutoff, precision);  }
                
                m->mothurOut("Cutoff was " + toString(cutoff) + " changed cutoff to " + toString(saveCutoff)); m->mothurOutEndLine();  
            }
            
            if (saveCutoff < smallestCutoff) { smallestCutoff = saveCutoff;  }
        }
        return listFileName;
        
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "clusterFile");
		exit(1);
	}
}
//**********************************************************************************************************************
string ClusterSplitCommand::runOptiCluster(string thisDistFile, string thisNamefile, set<string>& labels, double& smallestCutoff){
    try {
        string listFileName = "";
        
        string nameOrCount = "count";
        if (namefile != "") { nameOrCount = "name"; }
        
        OptiMatrix optiMatrix(thisDistFile, thisNamefile, nameOrCount, cutoff, -1);
        
        
        int iters = 0;
        double listVectorMetric = -1; //worst state
        while ((listVectorMetric < stableMetric) && (iters < maxIters)) {
            
            if (m->control_pressed) { break; }
        }
    
        return listFileName;
        
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterSplitCommand", "runOptiCluster");
        exit(1);
    }
}
//**********************************************************************************************************************

string ClusterSplitCommand::runVsearchCluster(string thisDistFile, string thisNamefile, set<string>& labels, double& smallestCutoff){
    try {

        m->mothurOutEndLine(); m->mothurOut("Clustering " + thisDistFile); m->mothurOutEndLine();
        
        string vsearchFastafile = ""; VsearchFileParser* vParse;
        if (namefile != "")                    { vParse = new VsearchFileParser(thisDistFile, thisNamefile, "name");      }
        else if (countfile != "")              { vParse = new VsearchFileParser(thisDistFile, thisNamefile, "count");    }
        else                                   { m->mothurOut("[ERROR]: Opps, should never get here. ClusterSplitCommand::runVsearchCluster() \n"); m->control_pressed = true; }
        
        if (m->control_pressed) {  return 0; }
        
        vsearchFastafile = vParse->getVsearchFile();
        
        if (cutoff > 1.0) {  m->mothurOut("You did not set a cutoff, using 0.03.\n"); cutoff = 0.03; }
        
        //Run vsearch
        string ucVsearchFile = m->getSimpleName(vsearchFastafile) + ".clustered.uc";
        string logfile = m->getSimpleName(vsearchFastafile) + ".clustered.log";
        vsearchDriver(vsearchFastafile, ucVsearchFile, logfile, smallestCutoff);
        
        if (m->control_pressed) { m->mothurRemove(ucVsearchFile); m->mothurRemove(logfile);  m->mothurRemove(vsearchFastafile); return 0; }
        
        if (outputDir == "") { outputDir += m->hasPath(thisDistFile); }
        tag = method;
        string listFileName = outputDir + m->getRootName(m->getSimpleName(thisDistFile)) + tag + ".list";
        
        //Convert outputted *.uc file into a list file
        vParse->createListFile(ucVsearchFile, listFileName, "", "", vParse->getNumBins(logfile), toString(cutoff));  delete vParse;
        
        //remove temp files
        m->mothurRemove(ucVsearchFile); m->mothurRemove(logfile);  m->mothurRemove(vsearchFastafile);
        
        if (deleteFiles) {
            m->mothurRemove(thisDistFile);
            m->mothurRemove(thisNamefile);
        }
        
        if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); } return 0; }
        
        labels.insert(toString(cutoff));
        
        return listFileName;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterSplitCommand", "runVsearchCluster");
        exit(1);
    }
}
//**********************************************************************************************************************

int ClusterSplitCommand::vsearchDriver(string inputFile, string ucClusteredFile, string logfile, double cutoff){
    try {
        
        //vsearch --maxaccepts 16 --usersort --id 0.97 --minseqlength 30 --wordlength 8 --uc $ROOT.clustered.uc --cluster_smallmem $ROOT.sorted.fna --maxrejects 64 --strand both --log $ROOT.clustered.log --sizeorder
        
        
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
        if (m->debug) {  m->mothurOut("[DEBUG]: vsearch cluster command = " + commandString + ".\n"); }
        
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
        m->errorOut(e, "ClusterSplitCommand", "vsearchDriver");
        exit(1);
    }
}
//**********************************************************************************************************************

int ClusterSplitCommand::createMergedDistanceFile(vector< map<string, string> > distNames) {
	try{
		
		string thisOutputDir = outputDir;
		if (outputDir == "") { thisOutputDir = m->hasPath(fastafile); }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(fastafile));
		string outputFileName = getOutputFileName("column", variables);
		m->mothurRemove(outputFileName);
		
		
		for (int i = 0; i < distNames.size(); i++) {
			if (m->control_pressed) {  return 0; }
			
			string thisDistFile = distNames[i].begin()->first;
			
			m->appendFiles(thisDistFile, outputFileName);
		}	
			
		outputTypes["column"].push_back(outputFileName); outputNames.push_back(outputFileName);
			
		return 0;	
		
		
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "createMergedDistanceFile");
		exit(1);
	}
}
//**********************************************************************************************************************
int ClusterSplitCommand::createRabund(CountTable*& ct, ListVector*& list, RAbundVector*& rabund){
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
		m->errorOut(e, "ClusterSplitCommand", "createRabund");
		exit(1);
	}
    
}
//**********************************************************************************************************************
string ClusterSplitCommand::printFile(string singleton, vector< map<string, string> >& distName){
    try {
        ofstream out;
        map<string, string> variables;
        string thisOutputDir = outputDir;
		if (outputDir == "") { thisOutputDir = m->hasPath(distfile); }
        variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(distfile));
		string outputFileName = getOutputFileName("file", variables);
        m->openOutputFile(outputFileName, out);
        outputTypes["file"].push_back(outputFileName); outputNames.push_back(outputFileName);
        m->setFileFile(outputFileName);
        
        out << singleton << endl;
        if (namefile != "") { out << "name" << endl; }
        else if (countfile != "") { out << "count" << endl; }
        else { out << "unknown" << endl; }
        
        for (int i = 0; i < distName.size(); i++) {    out << distName[i].begin()->first << '\t' << distName[i].begin()->second << endl;	}
        out.close();
        
        return outputFileName;
    }
    catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "printFile");
		exit(1);
	}
    
}
//**********************************************************************************************************************
string ClusterSplitCommand::readFile(vector< map<string, string> >& distName){
    try {
        
        string singleton, thiscolumn, thisname, type;
        
        ifstream in;
        m->openInputFile(file, in);
        
        in >> singleton; m->gobble(in);
        
        string path = m->hasPath(singleton);
        if (path == "") {  singleton = inputDir + singleton;  }
        
        in >> type; m->gobble(in);
        
        if (type == "name") { namefile = "name"; }
        else if (type == "count") { countfile = "count"; }
        else {  m->mothurOut("[ERROR]: unknown file type. Are the files in column 2 of the file name files or count files? Please change unknown to name or count.\n"); m->control_pressed = true; }
        
        if (isList) {

            vector<string> listFileNames;
            string thisListFileName = "";
            set<string> listLabels;
            
            while(!in.eof()) {
                if (m->control_pressed) { break; }
                
                in >> thisListFileName; m->gobble(in);
                
                string path = m->hasPath(thisListFileName);
                if (path == "") {  thisListFileName = inputDir + thisListFileName;  }
                
                getLabels(thisListFileName, listLabels);
                listFileNames.push_back(thisListFileName);
            }
            
            ListVector* listSingle;
            map<float, int> labelBins = completeListFile(listFileNames, singleton, listLabels, listSingle);
            
            mergeLists(listFileNames, labelBins, listSingle);
        
        }else {
            
            while(!in.eof()) {
                if (m->control_pressed) { break; }
                
                in >> thiscolumn; m->gobble(in);
                in >> thisname; m->gobble(in);
                
                map<string, string> temp;
                temp[thiscolumn] = thisname;
                distName.push_back(temp);
            }
        }
        
        in.close();
        
        return singleton;
    }
    catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "readFile");
		exit(1);
	}
    
}
//**********************************************************************************************************************
int ClusterSplitCommand::getLabels(string file, set<string>& listLabels){
    try {
        ifstream in;
        m->openInputFile(file, in);

        //read headers
        m->getline(in); m->gobble(in);
        
        string label;
        while(!in.eof()) {
            if (m->control_pressed) { break; }
            
            in >> label; m->getline(in); m->gobble(in);
            
            listLabels.insert(label);
        }
        
        in.close();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterSplitCommand", "getLabels");
        exit(1);
    }
    
}
//**********************************************************************************************************************
bool ClusterSplitCommand::findVsearch(){
    try {
        
        abort = false;
        
        if (cutoffNotSet) {  m->mothurOut("\nYou did not set a cutoff, using 0.03.\n"); cutoff = 0.03; }
        
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
        
        if (m->debug) {
            m->mothurOut("[DEBUG]: vsearch location using " + vsearchLocation + "\n");
        }
        
        if (!abort) { return true; }
        
        return false;

    }
    catch(exception& e) {
        m->errorOut(e, "ClusterSplitCommand", "findVsearch");
        exit(1);
    }
    
}
//**********************************************************************************************************************
