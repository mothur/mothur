/*
 *  clustersplitcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 5/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "clustersplitcommand.h"
#include "readcluster.h"
#include "splitmatrix.h"
#include "readphylip.h"
#include "readcolumn.h"
#include "readmatrix.hpp"
#include "inputdata.h"


//**********************************************************************************************************************
vector<string> ClusterSplitCommand::setParameters(){	
	try {
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "none", "FastaTaxName",false,false); parameters.push_back(ptaxonomy);
		CommandParameter pphylip("phylip", "InputTypes", "", "", "PhylipColumnFasta", "PhylipColumnFasta", "none",false,false); parameters.push_back(pphylip);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "PhylipColumnFasta", "PhylipColumnFasta", "FastaTaxName",false,false); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "ColumnName-FastaTaxName",false,false); parameters.push_back(pname);
		CommandParameter pcolumn("column", "InputTypes", "", "", "PhylipColumnFasta", "PhylipColumnFasta", "ColumnName",false,false); parameters.push_back(pcolumn);
		CommandParameter ptaxlevel("taxlevel", "Number", "", "1", "", "", "",false,false); parameters.push_back(ptaxlevel);
		CommandParameter psplitmethod("splitmethod", "Multiple", "classify-fasta-distance", "distance", "", "", "",false,false); parameters.push_back(psplitmethod);
		CommandParameter plarge("large", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(plarge);
		CommandParameter pshowabund("showabund", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(pshowabund);
		CommandParameter ptiming("timing", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(ptiming);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "",false,false); parameters.push_back(pprocessors);
		CommandParameter pcutoff("cutoff", "Number", "", "10", "", "", "",false,false); parameters.push_back(pcutoff);
		CommandParameter pprecision("precision", "Number", "", "100", "", "", "",false,false); parameters.push_back(pprecision);
		CommandParameter pmethod("method", "Multiple", "furthest-nearest-average-weighted", "furthest", "", "", "",false,false); parameters.push_back(pmethod);
		CommandParameter phard("hard", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(phard);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
			
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
		helpString += "The cluster.split command parameter options are fasta, phylip, column, name, cutoff, precision, method, splitmethod, taxonomy, taxlevel, showabund, timing, hard, large, processors. Fasta or Phylip or column and name are required.\n";
		helpString += "The cluster.split command can split your files in 3 ways. Splitting by distance file, by classification, or by classification also using a fasta file. \n";
		helpString += "For the distance file method, you need only provide your distance file and mothur will split the file into distinct groups. \n";
		helpString += "For the classification method, you need to provide your distance file and taxonomy file, and set the splitmethod to classify.  \n";
		helpString += "You will also need to set the taxlevel you want to split by. mothur will split the sequences into distinct taxonomy groups, and split the distance file based on those groups. \n";
		helpString += "For the classification method using a fasta file, you need to provide your fasta file, names file and taxonomy file.  \n";
		helpString += "You will also need to set the taxlevel you want to split by. mothur will split the sequence into distinct taxonomy groups, and create distance files for each grouping. \n";
		helpString += "The phylip and column parameter allow you to enter your distance file. \n";
		helpString += "The fasta parameter allows you to enter your aligned fasta file. \n";
		helpString += "The name parameter allows you to enter your name file and is required if your distance file is in column format. \n";
		helpString += "The cutoff parameter allow you to set the distance you want to cluster to, default is 10.0. \n";
		helpString += "The precision parameter allows you specify the precision of the precision of the distances outputted, default=100, meaning 2 decimal places. \n";
		helpString += "The method allows you to specify what clustering algorythm you want to use, default=furthest, option furthest, nearest, or average. \n";
		helpString += "The splitmethod parameter allows you to specify how you want to split your distance file before you cluster, default=distance, options distance, classify or fasta. \n";
		helpString += "The taxonomy parameter allows you to enter the taxonomy file for your sequences, this is only valid if you are using splitmethod=classify. Be sure your taxonomy file does not include the probability scores. \n";
		helpString += "The taxlevel parameter allows you to specify the taxonomy level you want to use to split the distance file, default=1, meaning use the first taxon in each list. \n";
		helpString += "The large parameter allows you to indicate that your distance matrix is too large to fit in RAM.  The default value is false.\n";
#ifdef USE_MPI
		helpString += "When using MPI, the processors parameter is set to the number of MPI processes running. \n";
#endif
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
ClusterSplitCommand::ClusterSplitCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["list"] = tempOutNames;
		outputTypes["rabund"] = tempOutNames;
		outputTypes["sabund"] = tempOutNames;
		outputTypes["column"] = tempOutNames;
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
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
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
			}
			
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
			
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }	
			else if (fastafile == "not found") { fastafile = ""; }
			else { distfile = fastafile;  splitmethod = "fasta";  }
			
			taxFile = validParameter.validFile(parameters, "taxonomy", true);
			if (taxFile == "not open") { abort = true; }	
			else if (taxFile == "not found") { taxFile = ""; }
			
			if ((phylipfile == "") && (columnfile == "") && (fastafile == "")) { 
				//is there are current file available for either of these?
				//give priority to column, then phylip, then fasta
				columnfile = m->getColumnFile(); 
				if (columnfile != "") {  m->mothurOut("Using " + columnfile + " as input file for the column parameter."); m->mothurOutEndLine(); }
				else { 
					phylipfile = m->getPhylipFile(); 
					if (phylipfile != "") {  m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter."); m->mothurOutEndLine(); }
					else { 
						fastafile = m->getFastaFile(); 
						if (fastafile != "") {  m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
						else { 
							m->mothurOut("No valid current files. When executing a cluster.split command you must enter a phylip or a column or fastafile."); m->mothurOutEndLine(); 
							abort = true; 
						}
					}
				}
			}
			else if ((phylipfile != "") && (columnfile != "") && (fastafile != "")) { m->mothurOut("When executing a cluster.split command you must enter ONLY ONE of the following: fasta, phylip or column."); m->mothurOutEndLine(); abort = true; }
		
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
			
			if (fastafile != "") {
				if (taxFile == "") { 
					taxFile = m->getTaxonomyFile(); 
					if (taxFile != "") {  m->mothurOut("Using " + taxFile + " as input file for the taxonomy parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("You need to provide a taxonomy file if you are if you are using a fasta file to generate the split."); m->mothurOutEndLine(); 
						abort = true; 
					}	
				}
				
				if (namefile == "") { 
					namefile = m->getNameFile(); 
					if (namefile != "") {  m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("You need to provide a namefile if you are if you are using a fasta file to generate the split."); m->mothurOutEndLine(); 
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
			
			temp = validParameter.validFile(parameters, "hard", false);			if (temp == "not found") { temp = "F"; }
			hard = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "large", false);			if (temp == "not found") { temp = "F"; }
			large = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			convert(temp, processors);
			
			temp = validParameter.validFile(parameters, "splitmethod", false);	
			if (splitmethod != "fasta") {
				if (temp == "not found")  { splitmethod = "distance"; }
				else {  splitmethod = temp; }
			}
			
			temp = validParameter.validFile(parameters, "cutoff", false);		if (temp == "not found")  { temp = "10"; }
			convert(temp, cutoff); 
			cutoff += (5 / (precision * 10.0));  
			
			temp = validParameter.validFile(parameters, "taxlevel", false);		if (temp == "not found")  { temp = "1"; }
			convert(temp, taxLevelCutoff); 
			
			method = validParameter.validFile(parameters, "method", false);		if (method == "not found") { method = "furthest"; }
			
			if ((method == "furthest") || (method == "nearest") || (method == "average")) { }
			else { m->mothurOut("Not a valid clustering method.  Valid clustering algorithms are furthest, nearest or average."); m->mothurOutEndLine(); abort = true; }
			
			if ((splitmethod == "distance") || (splitmethod == "classify") || (splitmethod == "fasta")) { }
			else { m->mothurOut(splitmethod + " is not a valid splitting method.  Valid splitting algorithms are distance, classify or fasta."); m->mothurOutEndLine(); abort = true; }
			
			if ((splitmethod == "classify") && (taxFile == "")) {  m->mothurOut("You need to provide a taxonomy file if you are going to use the classify splitmethod."); m->mothurOutEndLine(); abort = true;  }

			showabund = validParameter.validFile(parameters, "showabund", false);
			if (showabund == "not found") { showabund = "T"; }

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
		set<string> labels;
		string singletonName = "";
		double saveCutoff = cutoff;

		//****************** file prep work ******************************//
		#ifdef USE_MPI
			int pid;
			int tag = 2001;
			MPI_Status status; 
			MPI_Comm_size(MPI_COMM_WORLD, &processors); //set processors to the number of mpi processes running
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
			
			if (pid == 0) { //only process 0 converts and splits
			
		#endif
		
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
			if (namefile == "") {  //you need to make a namefile for split matrix
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
		if (splitmethod == "distance")			{	split = new SplitMatrix(distfile, namefile, taxFile, cutoff, splitmethod, large);							}
		else if (splitmethod == "classify")		{	split = new SplitMatrix(distfile, namefile, taxFile, taxLevelCutoff, splitmethod, large);					}
		else if (splitmethod == "fasta")		{	split = new SplitMatrix(fastafile, namefile, taxFile, taxLevelCutoff, cutoff, splitmethod, processors, outputDir);	}
		else { m->mothurOut("Not a valid splitting method.  Valid splitting algorithms are distance, classify or fasta."); m->mothurOutEndLine(); return 0;		}
		
		split->split();
		
		if (m->control_pressed) { delete split; return 0; }
		
		singletonName = split->getSingletonNames();
		vector< map<string, string> > distName = split->getDistanceFiles();  //returns map of distance files -> namefile sorted by distance file size
		delete split;
		
		//output a merged distance file
		if (splitmethod == "fasta")		{ createMergedDistanceFile(distName); }
			
				
		if (m->control_pressed) { return 0; }
		
		m->mothurOut("It took " + toString(time(NULL) - estart) + " seconds to split the distance file."); m->mothurOutEndLine();
		estart = time(NULL);
		
		//****************** break up files between processes and cluster each file set ******************************//
	#ifdef USE_MPI
			////you are process 0 from above////
			
			vector < vector < map<string, string> > > dividedNames; //distNames[1] = vector of filenames for process 1...				
			dividedNames.resize(processors);
					
			//for each file group figure out which process will complete it
			//want to divide the load intelligently so the big files are spread between processes
			for (int i = 0; i < distName.size(); i++) { 
				int processToAssign = (i+1) % processors; 
				if (processToAssign == 0) { processToAssign = processors; }
						
				dividedNames[(processToAssign-1)].push_back(distName[i]);
			}
					
			//not lets reverse the order of ever other process, so we balance big files running with little ones
			for (int i = 0; i < processors; i++) {
				int remainder = ((i+1) % processors);
				if (remainder) {  reverse(dividedNames[i].begin(), dividedNames[i].end());  }
			}
			
			
			//send each child the list of files it needs to process
			for(int i = 1; i < processors; i++) { 
				//send number of file pairs
				int num = dividedNames[i].size();
				MPI_Send(&num, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
				
				for (int j = 0; j < num; j++) { //send filenames to process i
					char tempDistFileName[1024];
					strcpy(tempDistFileName, (dividedNames[i][j].begin()->first).c_str());
					int lengthDist = (dividedNames[i][j].begin()->first).length();
					
					MPI_Send(&lengthDist, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
					MPI_Send(tempDistFileName, 1024, MPI_CHAR, i, tag, MPI_COMM_WORLD);
					
					char tempNameFileName[1024];
					strcpy(tempNameFileName, (dividedNames[i][j].begin()->second).c_str());
					int lengthName = (dividedNames[i][j].begin()->second).length();

					MPI_Send(&lengthName, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
					MPI_Send(tempNameFileName, 1024, MPI_CHAR, i, tag, MPI_COMM_WORLD);
				}
			}
			
			//process your share
			listFileNames = cluster(dividedNames[0], labels);
			
			//receive the other processes info
			for(int i = 1; i < processors; i++) { 
				int num = dividedNames[i].size();
				
				double tempCutoff;
				MPI_Recv(&tempCutoff, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
				if (tempCutoff < cutoff) { cutoff = tempCutoff; }
				
				//send list filenames to root process
				for (int j = 0; j < num; j++) {  
					int lengthList = 0;
					char tempListFileName[1024];
				
					MPI_Recv(&lengthList, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
					MPI_Recv(tempListFileName, 1024, MPI_CHAR, i, tag, MPI_COMM_WORLD, &status); 
				
					string myListFileName = tempListFileName;
					myListFileName = myListFileName.substr(0, lengthList);
					
					listFileNames.push_back(myListFileName);
				}
				
				//send Labels to root process
				int numLabels = 0;
				MPI_Recv(&numLabels, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
				
				for (int j = 0; j < numLabels; j++) {  
					int lengthLabel = 0;
					char tempLabel[100];
				
					MPI_Recv(&lengthLabel, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
					MPI_Recv(tempLabel, 100, MPI_CHAR, i, tag, MPI_COMM_WORLD, &status); 
				
					string myLabel = tempLabel;
					myLabel = myLabel.substr(0, lengthLabel);
			
					if (labels.count(myLabel) == 0) { labels.insert(myLabel); }
				}
			}
			
		}else { //you are a child process
			vector < map<string, string> >  myNames;
			
			//recieve the files you need to process
			//receive number of file pairs
			int num = 0;
			MPI_Recv(&num, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
			
			myNames.resize(num);
	
			for (int j = 0; j < num; j++) { //receive filenames to process 
				int lengthDist = 0;
				char tempDistFileName[1024];
				
				MPI_Recv(&lengthDist, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(tempDistFileName, 1024, MPI_CHAR, 0, tag, MPI_COMM_WORLD, &status); 
				
				string myDistFileName = tempDistFileName;
				myDistFileName = myDistFileName.substr(0, lengthDist);
			
				int lengthName = 0;
				char tempNameFileName[1024];
				
				MPI_Recv(&lengthName, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(tempNameFileName, 1024, MPI_CHAR, 0, tag, MPI_COMM_WORLD, &status); 
				
				string myNameFileName = tempNameFileName;
				myNameFileName = myNameFileName.substr(0, lengthName);
				
				//save file name
				myNames[j][myDistFileName] = myNameFileName;
			}
	
			//process them
			listFileNames = cluster(myNames, labels);
			
			//send cutoff
			MPI_Send(&cutoff, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
			
			//send list filenames to root process
			for (int j = 0; j < num; j++) {  
				char tempListFileName[1024];
				strcpy(tempListFileName, listFileNames[j].c_str());
				int lengthList = listFileNames[j].length();
					
				MPI_Send(&lengthList, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
				MPI_Send(tempListFileName, 1024, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
			}
			
			//send Labels to root process
			int numLabels = labels.size();
			MPI_Send(&numLabels, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
			
			for(set<string>::iterator it = labels.begin(); it != labels.end(); ++it) {
				char tempLabel[100];
				strcpy(tempLabel, (*it).c_str());
				int lengthLabel = (*it).length();
					
				MPI_Send(&lengthLabel, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
				MPI_Send(tempLabel, 100, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
			}
		}
		
		//make everyone wait
		MPI_Barrier(MPI_COMM_WORLD);
		
	#else

		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				if(processors == 1){
					listFileNames = cluster(distName, labels); //clusters individual files and returns names of list files
				}else{
					vector < vector < map<string, string> > > dividedNames; //distNames[1] = vector of filenames for process 1...
					dividedNames.resize(processors);
					
					//for each file group figure out which process will complete it
					//want to divide the load intelligently so the big files are spread between processes
					for (int i = 0; i < distName.size(); i++) { 
						int processToAssign = (i+1) % processors; 
						if (processToAssign == 0) { processToAssign = processors; }
						
						dividedNames[(processToAssign-1)].push_back(distName[i]);
					}
					
					//not lets reverse the order of ever other process, so we balance big files running with little ones
					for (int i = 0; i < processors; i++) {
						int remainder = ((i+1) % processors);
						if (remainder) {  reverse(dividedNames[i].begin(), dividedNames[i].end());  }
					}
					
					createProcesses(dividedNames);
							
					if (m->control_pressed) { return 0; }

					//get list of list file names from each process
					for(int i=0;i<processors;i++){
						string filename = toString(processIDS[i]) + ".temp";
						ifstream in;
						m->openInputFile(filename, in);
						
						in >> tag; m->gobble(in);
						
						while(!in.eof()) {
							string tempName;
							in >> tempName; m->gobble(in);
							listFileNames.push_back(tempName);
						}
						in.close();
						remove((toString(processIDS[i]) + ".temp").c_str());
						
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
						remove((toString(processIDS[i]) + ".temp.labels").c_str());
					}
				}
		#else
				listFileNames = cluster(distName, labels); //clusters individual files and returns names of list files
		#endif
	#endif	
		if (m->control_pressed) { for (int i = 0; i < listFileNames.size(); i++) { remove(listFileNames[i].c_str()); } return 0; }
		
		if (saveCutoff != cutoff) { m->mothurOut("Cutoff was " + toString(saveCutoff) + " changed cutoff to " + toString(cutoff)); m->mothurOutEndLine();  }
		
		m->mothurOut("It took " + toString(time(NULL) - estart) + " seconds to cluster"); m->mothurOutEndLine();
		
		//****************** merge list file and create rabund and sabund files ******************************//
		estart = time(NULL);
		m->mothurOut("Merging the clustered files..."); m->mothurOutEndLine();
		
		#ifdef USE_MPI
			if (pid == 0) { //only process 0 merges
		#endif

		ListVector* listSingle;
		map<float, int> labelBins = completeListFile(listFileNames, singletonName, labels, listSingle); //returns map of label to numBins
		
		if (m->control_pressed) { if (listSingle != NULL) { delete listSingle; } for (int i = 0; i < outputNames.size(); i++) { remove(outputNames[i].c_str()); } return 0; }
		
		mergeLists(listFileNames, labelBins, listSingle);

		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) { remove(outputNames[i].c_str()); } return 0; }
		
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
		
		#ifdef USE_MPI
			} //only process 0 merges
			
			//make everyone wait
			MPI_Barrier(MPI_COMM_WORLD);
		#endif

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
			while (!in.eof()) {
				in >> firstCol >> secondCol; m->gobble(in);
				listSingle->push_back(secondCol);
			}
			in.close();
			remove(singleton.c_str());
			
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
				if (listSingle != NULL) { delete listSingle; listSingle = NULL; remove(singleton.c_str());  }
				for (int i = 0; i < listNames.size(); i++) {   remove(listNames[i].c_str());  }
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
				list->print(outFilled);
		
				//update labelBin
				labelBin[orderFloat[l]] += list->getNumBins();
									
				delete list;
									
				list = input->getListVector();
			}
			
			if (list != NULL) { delete list; }
			delete input;
			
			outFilled.close();
			remove(listNames[k].c_str());
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
		
		m->openOutputFile(fileroot+ tag + ".sabund",	outSabund);
		m->openOutputFile(fileroot+ tag + ".rabund",	outRabund);
		m->openOutputFile(fileroot+ tag + ".list",		outList);
				
		outputNames.push_back(fileroot+ tag + ".sabund");  outputTypes["list"].push_back(fileroot+ tag + ".list");
		outputNames.push_back(fileroot+ tag + ".rabund");  outputTypes["rabund"].push_back(fileroot+ tag + ".rabund");
		outputNames.push_back(fileroot+ tag + ".list");    outputTypes["sabund"].push_back(fileroot+ tag + ".sabund");
		
		map<float, int>::iterator itLabel;

		//for each label needed
		for(itLabel = userLabels.begin(); itLabel != userLabels.end(); itLabel++) {
			
			string thisLabel;
			if (itLabel->first == -1) { thisLabel = "unique"; }
			else { thisLabel = toString(itLabel->first,  length-1);  } 
			
			outList << thisLabel << '\t' << itLabel->second << '\t';

			RAbundVector* rabund = new RAbundVector();
			rabund->setLabel(thisLabel);

			//add in singletons
			if (listSingle != NULL) {
				for (int j = 0; j < listSingle->getNumBins(); j++) {
					outList << listSingle->get(j) << '\t';
					rabund->push_back(m->getNumNames(listSingle->get(j)));
				}
			}
			
			//get the list info from each file
			for (int k = 0; k < listNames.size(); k++) {
	
				if (m->control_pressed) {  if (listSingle != NULL) { delete listSingle;   } for (int i = 0; i < listNames.size(); i++) { remove(listNames[i].c_str());  } delete rabund; return 0; }
				
				InputData* input = new InputData(listNames[k], "list");
				ListVector* list = input->getListVector(thisLabel);
				
				//this file has reached the end
				if (list == NULL) { m->mothurOut("Error merging listvectors in file " + listNames[k]); m->mothurOutEndLine();  }	
				else {		
					for (int j = 0; j < list->getNumBins(); j++) {
						outList << list->get(j) << '\t';
						rabund->push_back(m->getNumNames(list->get(j)));
					}
					delete list;
				}
				delete input;
			}
			
			SAbundVector sabund = rabund->getSAbundVector();
			
			sabund.print(outSabund);
			rabund->print(outRabund);
			outList << endl;
			
			delete rabund;
		}
		
		outList.close();
		outRabund.close();
		outSabund.close();
		
		if (listSingle != NULL) { delete listSingle;  }
		
		for (int i = 0; i < listNames.size(); i++) {  remove(listNames[i].c_str());  }
		
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
	
		oldList->print(outList);
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "printData");
		exit(1);
	}
}
//**********************************************************************************************************************
int ClusterSplitCommand::createProcesses(vector < vector < map<string, string> > > dividedNames){
	try {
	
	#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 0;
		int exitCommand = 1;
		processIDS.clear();
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				set<string> labels;
				vector<string> listFileNames = cluster(dividedNames[process], labels);
				
				//write out names to file
				string filename = toString(getpid()) + ".temp";
				ofstream out;
				m->openOutputFile(filename, out);
				out << tag << endl;
				for (int j = 0; j < listFileNames.size(); j++) { out << listFileNames[j] << endl;  }
				out.close();
				
				//print out labels
				ofstream outLabels;
				filename = toString(getpid()) + ".temp.labels";
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
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		return exitCommand;
	#endif		
	
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "createProcesses");
		exit(1);
	}
}
//**********************************************************************************************************************

vector<string> ClusterSplitCommand::cluster(vector< map<string, string> > distNames, set<string>& labels){
	try {
		Cluster* cluster;
		SparseMatrix* matrix;
		ListVector* list;
		ListVector oldList;
		RAbundVector* rabund;
		
		vector<string> listFileNames;
		
		double smallestCutoff = cutoff;
		
		//cluster each distance file
		for (int i = 0; i < distNames.size(); i++) {
			if (m->control_pressed) { return listFileNames; }
			
			string thisNamefile = distNames[i].begin()->second;
			string thisDistFile = distNames[i].begin()->first;
						
			#ifdef USE_MPI
				int pid;
				MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
				
				//output your files too
				if (pid != 0) {
					cout << endl << "Reading " << thisDistFile << endl;
				}
			#endif
			
			m->mothurOutEndLine(); m->mothurOut("Reading " + thisDistFile); m->mothurOutEndLine();
			
			ReadMatrix* read = new ReadColumnMatrix(thisDistFile); 	
			read->setCutoff(cutoff);

			NameAssignment* nameMap = new NameAssignment(thisNamefile);
			nameMap->readMap();
			read->read(nameMap);
			
			if (m->control_pressed) {  delete read; delete nameMap; return listFileNames; }
			
			list = read->getListVector();
			oldList = *list;
			matrix = read->getMatrix();
			
			delete read; 
			delete nameMap; 
			
			
			#ifdef USE_MPI
				//output your files too
				if (pid != 0) {
					cout << endl << "Clustering " << thisDistFile << endl;
				}
			#endif
			
			m->mothurOutEndLine(); m->mothurOut("Clustering " + thisDistFile); m->mothurOutEndLine();
		
			rabund = new RAbundVector(list->getRAbundVector());
			
			//create cluster
			if (method == "furthest")	{	cluster = new CompleteLinkage(rabund, list, matrix, cutoff, method); }
			else if(method == "nearest"){	cluster = new SingleLinkage(rabund, list, matrix, cutoff, method); }
			else if(method == "average"){	cluster = new AverageLinkage(rabund, list, matrix, cutoff, method);	}
			tag = cluster->getTag();
		
			if (outputDir == "") { outputDir += m->hasPath(thisDistFile); }
			fileroot = outputDir + m->getRootName(m->getSimpleName(thisDistFile));
			
			ofstream listFile;
			m->openOutputFile(fileroot+ tag + ".list",	listFile);
		
			listFileNames.push_back(fileroot+ tag + ".list");
				
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
					for (int i = 0; i < listFileNames.size(); i++) {	remove(listFileNames[i].c_str()); 	}
					listFileNames.clear(); return listFileNames;
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
			listFile.close();
			
			if (m->control_pressed) { //clean up
				for (int i = 0; i < listFileNames.size(); i++) {	remove(listFileNames[i].c_str()); 	}
				listFileNames.clear(); return listFileNames;
			}
			
			remove(thisDistFile.c_str());
			remove(thisNamefile.c_str());
			
			if (saveCutoff != cutoff) { 
				if (hard)	{  saveCutoff = m->ceilDist(saveCutoff, precision);	}
				else		{	saveCutoff = m->roundDist(saveCutoff, precision);  }
			
				m->mothurOut("Cutoff was " + toString(cutoff) + " changed cutoff to " + toString(saveCutoff)); m->mothurOutEndLine();  
			}
			
			if (saveCutoff < smallestCutoff) { smallestCutoff = saveCutoff;  }
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

int ClusterSplitCommand::createMergedDistanceFile(vector< map<string, string> > distNames) {
	try{
		
#ifdef USE_MPI
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
		
		if (pid != 0) {
#endif
		
		string thisOutputDir = outputDir;
		if (outputDir == "") { thisOutputDir = m->hasPath(fastafile); }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(fastafile)) + "dist";
		remove(outputFileName.c_str());
		
		
		for (int i = 0; i < distNames.size(); i++) {
			if (m->control_pressed) {  return 0; }
			
			string thisDistFile = distNames[i].begin()->first;
			
			m->appendFiles(thisDistFile, outputFileName);
		}	
			
		outputTypes["column"].push_back(outputFileName); outputNames.push_back(outputFileName);
			
#ifdef USE_MPI
		}
#endif
				
		return 0;	
		
		
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "createMergedDistanceFile");
		exit(1);
	}
}
//**********************************************************************************************************************
