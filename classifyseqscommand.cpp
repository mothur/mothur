/*
 *  classifyseqscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 11/2/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "classifyseqscommand.h"
#include "sequence.hpp"
#include "bayesian.h"
#include "phylotree.h"
#include "phylosummary.h"
#include "knn.h"


//**********************************************************************************************************************
vector<string> ClassifySeqsCommand::getValidParameters(){	
	try {
		string AlignArray[] =  {"template","fasta","name","group","search","ksize","method","processors","taxonomy","match","mismatch","gapopen","gapextend","numwanted","cutoff","probs","iters", "outputdir","inputdir"};
		vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
ClassifySeqsCommand::ClassifySeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		vector<string> tempOutNames;
		outputTypes["taxonomy"] = tempOutNames;
		outputTypes["taxsummary"] = tempOutNames;
		outputTypes["matchdist"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "ClassifySeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> ClassifySeqsCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"fasta","template","taxonomy"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> ClassifySeqsCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************
ClassifySeqsCommand::ClassifySeqsCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			
			//valid paramters for this command
			string AlignArray[] =  {"template","fasta","name","group","search","ksize","method","processors","taxonomy","match","mismatch","gapopen","gapextend","numwanted","cutoff","probs","iters", "outputdir","inputdir"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters(); 
			
			ValidParameters validParameter("classify.seqs");
			map<string, string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["taxonomy"] = tempOutNames;
			outputTypes["taxsummary"] = tempOutNames;
			outputTypes["matchdist"] = tempOutNames;
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("template");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["template"] = inputDir + it->second;		}
				}
				
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
			}

			//check for required parameters
			templateFileName = validParameter.validFile(parameters, "template", true);
			if (templateFileName == "not found") { 
				m->mothurOut("template is a required parameter for the classify.seqs command."); 
				m->mothurOutEndLine();
				abort = true; 
			}
			else if (templateFileName == "not open") { abort = true; }	
			
						
			fastaFileName = validParameter.validFile(parameters, "fasta", false);
			if (fastaFileName == "not found") { m->mothurOut("fasta is a required parameter for the classify.seqs command."); m->mothurOutEndLine(); abort = true;  }
			else { 
				m->splitAtDash(fastaFileName, fastaFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastaFileNames.size(); i++) {
					if (inputDir != "") {
						string path = m->hasPath(fastaFileNames[i]);
						//if the user has not given a path then, add inputdir. else leave path alone.
						if (path == "") {	fastaFileNames[i] = inputDir + fastaFileNames[i];		}
					}
					
					int ableToOpen;
					
					ifstream in;
					ableToOpen = m->openInputFile(fastaFileNames[i], in, "noerror");
				
					//if you can't open it, try default location
					if (ableToOpen == 1) {
						if (m->getDefaultPath() != "") { //default path is set
							string tryPath = m->getDefaultPath() + m->getSimpleName(fastaFileNames[i]);
							m->mothurOut("Unable to open " + fastaFileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
							ifstream in2;
							ableToOpen = m->openInputFile(tryPath, in2, "noerror");
							in2.close();
							fastaFileNames[i] = tryPath;
						}
					}
					
					if (ableToOpen == 1) {
						if (m->getOutputDir() != "") { //default path is set
							string tryPath = m->getOutputDir() + m->getSimpleName(fastaFileNames[i]);
							m->mothurOut("Unable to open " + fastaFileNames[i] + ". Trying output directory " + tryPath); m->mothurOutEndLine();
							ifstream in2;
							ableToOpen = m->openInputFile(tryPath, in2, "noerror");
							in2.close();
							fastaFileNames[i] = tryPath;
						}
					}
					
					in.close();
					
					if (ableToOpen == 1) { 
						m->mothurOut("Unable to open " + fastaFileNames[i] + ". It will be disregarded."); m->mothurOutEndLine(); 
						//erase from file list
						fastaFileNames.erase(fastaFileNames.begin()+i);
						i--;
					}
					
				}
				
				//make sure there is at least one valid file left
				if (fastaFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}

			
			taxonomyFileName = validParameter.validFile(parameters, "taxonomy", true);
			if (taxonomyFileName == "not found") { 
				m->mothurOut("taxonomy is a required parameter for the classify.seqs command."); 
				m->mothurOutEndLine();
				abort = true; 
			}
			else if (taxonomyFileName == "not open") { abort = true; }	
			
			
			namefile = validParameter.validFile(parameters, "name", false);
			if (namefile == "not found") { namefile = "";  }

			else { 
				m->splitAtDash(namefile, namefileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < namefileNames.size(); i++) {
					if (inputDir != "") {
						string path = m->hasPath(namefileNames[i]);
						//if the user has not given a path then, add inputdir. else leave path alone.
						if (path == "") {	namefileNames[i] = inputDir + namefileNames[i];		}
					}
					int ableToOpen;
					
					ifstream in;
					ableToOpen = m->openInputFile(namefileNames[i], in, "noerror");
				
					//if you can't open it, try default location
					if (ableToOpen == 1) {
						if (m->getDefaultPath() != "") { //default path is set
							string tryPath = m->getDefaultPath() + m->getSimpleName(namefileNames[i]);
							m->mothurOut("Unable to open " + namefileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
							ifstream in2;
							ableToOpen = m->openInputFile(tryPath, in2, "noerror");
							in2.close();
							namefileNames[i] = tryPath;
						}
					}
					
					if (ableToOpen == 1) {
						if (m->getOutputDir() != "") { //default path is set
							string tryPath = m->getOutputDir() + m->getSimpleName(namefileNames[i]);
							m->mothurOut("Unable to open " + namefileNames[i] + ". Trying output directory " + tryPath); m->mothurOutEndLine();
							ifstream in2;
							ableToOpen = m->openInputFile(tryPath, in2, "noerror");
							in2.close();
							namefileNames[i] = tryPath;
						}
					}
					in.close();
					
					if (ableToOpen == 1) { 
						m->mothurOut("Unable to open " + namefileNames[i] + ". It will be disregarded."); m->mothurOutEndLine();  abort = true;
						//erase from file list
						namefileNames.erase(namefileNames.begin()+i);
						i--;
					}

				}
			}

			if (namefile != "") {
				if (namefileNames.size() != fastaFileNames.size()) { abort = true; m->mothurOut("If you provide a name file, you must have one for each fasta file."); m->mothurOutEndLine(); }
			}
			
			groupfile = validParameter.validFile(parameters, "group", false);
			if (groupfile == "not found") { groupfile = "";  }
			else { 
				m->splitAtDash(groupfile, groupfileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < groupfileNames.size(); i++) {
					if (inputDir != "") {
						string path = m->hasPath(groupfileNames[i]);
						//if the user has not given a path then, add inputdir. else leave path alone.
						if (path == "") {	groupfileNames[i] = inputDir + groupfileNames[i];		}
					}
					int ableToOpen;
					
					ifstream in;
					ableToOpen = m->openInputFile(groupfileNames[i], in, "noerror");
				
					//if you can't open it, try default location
					if (ableToOpen == 1) {
						if (m->getDefaultPath() != "") { //default path is set
							string tryPath = m->getDefaultPath() + m->getSimpleName(groupfileNames[i]);
							m->mothurOut("Unable to open " + groupfileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
							ifstream in2;
							ableToOpen = m->openInputFile(tryPath, in2, "noerror");
							in2.close();
							groupfileNames[i] = tryPath;
						}
					}
					
					if (ableToOpen == 1) {
						if (m->getOutputDir() != "") { //default path is set
							string tryPath = m->getOutputDir() + m->getSimpleName(groupfileNames[i]);
							m->mothurOut("Unable to open " + groupfileNames[i] + ". Trying output directory " + tryPath); m->mothurOutEndLine();
							ifstream in2;
							ableToOpen = m->openInputFile(tryPath, in2, "noerror");
							in2.close();
							groupfileNames[i] = tryPath;
						}
					}
					
					in.close();
					
					if (ableToOpen == 1) { 
						m->mothurOut("Unable to open " + groupfileNames[i] + ". It will be disregarded."); m->mothurOutEndLine(); groupfileNames[i] = "";
						//erase from file list
						groupfileNames.erase(groupfileNames.begin()+i);
						i--;
					}
				}
			}

			if (groupfile != "") {
				if (groupfileNames.size() != fastaFileNames.size()) { abort = true; m->mothurOut("If you provide a group file, you must have one for each fasta file."); m->mothurOutEndLine(); }
			}else {
				for (int i = 0; i < fastaFileNames.size(); i++) {  groupfileNames.push_back("");  }
			}
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.validFile(parameters, "ksize", false);		if (temp == "not found"){	temp = "8";				}
			convert(temp, kmerSize); 
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = "1";				}
			convert(temp, processors); 
			
			search = validParameter.validFile(parameters, "search", false);		if (search == "not found"){	search = "kmer";		}
			
			method = validParameter.validFile(parameters, "method", false);		if (method == "not found"){	method = "bayesian";	}
			
			temp = validParameter.validFile(parameters, "match", false);		if (temp == "not found"){	temp = "1.0";			}
			convert(temp, match);  
			
			temp = validParameter.validFile(parameters, "mismatch", false);		if (temp == "not found"){	temp = "-1.0";			}
			convert(temp, misMatch);  
			
			temp = validParameter.validFile(parameters, "gapopen", false);		if (temp == "not found"){	temp = "-2.0";			}
			convert(temp, gapOpen);  
			
			temp = validParameter.validFile(parameters, "gapextend", false);	if (temp == "not found"){	temp = "-1.0";			}
			convert(temp, gapExtend); 
			
			temp = validParameter.validFile(parameters, "numwanted", false);	if (temp == "not found"){	temp = "10";			}
			convert(temp, numWanted);
			
			temp = validParameter.validFile(parameters, "cutoff", false);		if (temp == "not found"){	temp = "0";				}
			convert(temp, cutoff);
			
			temp = validParameter.validFile(parameters, "probs", false);		if (temp == "not found"){	temp = "true";			}
			probs = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "iters", false);		if (temp == "not found") { temp = "100";			}
			convert(temp, iters); 


			
			if ((method == "bayesian") && (search != "kmer"))  { 
				m->mothurOut("The bayesian method requires the kmer search." + search + "will be disregarded." ); m->mothurOutEndLine();
				search = "kmer";
			}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "ClassifySeqsCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

ClassifySeqsCommand::~ClassifySeqsCommand(){	

	if (abort == false) {
		for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
	}
}

//**********************************************************************************************************************

void ClassifySeqsCommand::help(){
	try {
		m->mothurOut("The classify.seqs command reads a fasta file containing sequences and creates a .taxonomy file and a .tax.summary file.\n");
		m->mothurOut("The classify.seqs command parameters are template, fasta, name, search, ksize, method, taxonomy, processors, match, mismatch, gapopen, gapextend, numwanted and probs.\n");
		m->mothurOut("The template, fasta and taxonomy parameters are required. You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amzon.fasta \n");
		m->mothurOut("The search parameter allows you to specify the method to find most similar template.  Your options are: suffix, kmer, blast and distance. The default is kmer.\n");
		m->mothurOut("The name parameter allows you add a names file with your fasta file, if you enter multiple fasta files, you must enter matching names files for them.\n");
		m->mothurOut("The group parameter allows you add a group file so you can have the summary totals broken up by group.\n");
		m->mothurOut("The method parameter allows you to specify classification method to use.  Your options are: bayesian and knn. The default is bayesian.\n");
		m->mothurOut("The ksize parameter allows you to specify the kmer size for finding most similar template to candidate.  The default is 8.\n");
		m->mothurOut("The processors parameter allows you to specify the number of processors to use. The default is 1.\n");
		#ifdef USE_MPI
		m->mothurOut("When using MPI, the processors parameter is set to the number of MPI processes running. \n");
		#endif
		m->mothurOut("The match parameter allows you to specify the bonus for having the same base. The default is 1.0.\n");
		m->mothurOut("The mistmatch parameter allows you to specify the penalty for having different bases.  The default is -1.0.\n");
		m->mothurOut("The gapopen parameter allows you to specify the penalty for opening a gap in an alignment. The default is -2.0.\n");
		m->mothurOut("The gapextend parameter allows you to specify the penalty for extending a gap in an alignment.  The default is -1.0.\n");
		m->mothurOut("The numwanted parameter allows you to specify the number of sequence matches you want with the knn method.  The default is 10.\n");
		m->mothurOut("The cutoff parameter allows you to specify a bootstrap confidence threshold for your taxonomy.  The default is 0.\n");
		m->mothurOut("The probs parameter shuts off the bootstrapping results for the bayesian method. The default is true, meaning you want the bootstrapping to be shown.\n");
		m->mothurOut("The iters parameter allows you to specify how many iterations to do when calculating the bootstrap confidence score for your taxonomy with the bayesian method.  The default is 100.\n");
		m->mothurOut("The classify.seqs command should be in the following format: \n");
		m->mothurOut("classify.seqs(template=yourTemplateFile, fasta=yourFastaFile, method=yourClassificationMethod, search=yourSearchmethod, ksize=yourKmerSize, taxonomy=yourTaxonomyFile, processors=yourProcessors) \n");
		m->mothurOut("Example classify.seqs(fasta=amazon.fasta, template=core.filtered, method=knn, search=gotoh, ksize=8, processors=2)\n");
		m->mothurOut("The .taxonomy file consists of 2 columns: 1 = your sequence name, 2 = the taxonomy for your sequence. \n");
		m->mothurOut("The .tax.summary is a summary of the different taxonomies represented in your fasta file. \n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "help");
		exit(1);
	}
}


//**********************************************************************************************************************

int ClassifySeqsCommand::execute(){
	try {
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		if(method == "bayesian"){	classify = new Bayesian(taxonomyFileName, templateFileName, search, kmerSize, cutoff, iters);		}
		else if(method == "knn"){	classify = new Knn(taxonomyFileName, templateFileName, search, kmerSize, gapOpen, gapExtend, match, misMatch, numWanted);				}
		else {
			m->mothurOut(search + " is not a valid method option. I will run the command using bayesian.");
			m->mothurOutEndLine();
			classify = new Bayesian(taxonomyFileName, templateFileName, search, kmerSize, cutoff, iters);	
		}
		
		if (m->control_pressed) { delete classify; return 0; }
		
				
		for (int s = 0; s < fastaFileNames.size(); s++) {
		
			m->mothurOut("Classifying sequences from " + fastaFileNames[s] + " ..." ); m->mothurOutEndLine();
			
			string RippedTaxName = m->getRootName(m->getSimpleName(taxonomyFileName));
			RippedTaxName = m->getExtension(RippedTaxName.substr(0, RippedTaxName.length()-1));
			if (RippedTaxName[0] == '.') { RippedTaxName = RippedTaxName.substr(1, RippedTaxName.length()); }
			RippedTaxName +=  "."; 
		
			if (outputDir == "") { outputDir += m->hasPath(fastaFileNames[s]); }
			string newTaxonomyFile = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + RippedTaxName + "taxonomy";
			string tempTaxonomyFile = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + "taxonomy.temp";
			string taxSummary = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + RippedTaxName + "tax.summary";
			
			if ((method == "knn") && (search == "distance")) { 
				string DistName = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + "match.dist";
				classify->setDistName(DistName);  outputNames.push_back(DistName); outputTypes["matchdist"].push_back(DistName);
			}
			
			outputNames.push_back(newTaxonomyFile); outputTypes["taxonomy"].push_back(newTaxonomyFile);
			outputNames.push_back(taxSummary);	outputTypes["taxsummary"].push_back(taxSummary);
			
			int start = time(NULL);
			int numFastaSeqs = 0;
			for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
			
#ifdef USE_MPI	
				int pid, numSeqsPerProcessor; 
				int tag = 2001;
				vector<unsigned long int> MPIPos;
				
				MPI_Status status; 
				MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
				MPI_Comm_size(MPI_COMM_WORLD, &processors); 

				MPI_File inMPI;
				MPI_File outMPINewTax;
				MPI_File outMPITempTax;
							
				int outMode=MPI_MODE_CREATE|MPI_MODE_WRONLY; 
				int inMode=MPI_MODE_RDONLY; 
				
				char outNewTax[1024];
				strcpy(outNewTax, newTaxonomyFile.c_str());
				
				char outTempTax[1024];
				strcpy(outTempTax, tempTaxonomyFile.c_str());
				
				char inFileName[1024];
				strcpy(inFileName, fastaFileNames[s].c_str());

				MPI_File_open(MPI_COMM_WORLD, inFileName, inMode, MPI_INFO_NULL, &inMPI);  //comm, filename, mode, info, filepointer
				MPI_File_open(MPI_COMM_WORLD, outNewTax, outMode, MPI_INFO_NULL, &outMPINewTax);
				MPI_File_open(MPI_COMM_WORLD, outTempTax, outMode, MPI_INFO_NULL, &outMPITempTax);
				
				if (m->control_pressed) { outputTypes.clear(); MPI_File_close(&inMPI);  MPI_File_close(&outMPINewTax);   MPI_File_close(&outMPITempTax);  delete classify; return 0;  }
				
				if (pid == 0) { //you are the root process 
					
					MPIPos = m->setFilePosFasta(fastaFileNames[s], numFastaSeqs); //fills MPIPos, returns numSeqs
					
					//send file positions to all processes
					for(int i = 1; i < processors; i++) { 
						MPI_Send(&numFastaSeqs, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
						MPI_Send(&MPIPos[0], (numFastaSeqs+1), MPI_LONG, i, tag, MPI_COMM_WORLD);
					}
					
					//figure out how many sequences you have to align
					numSeqsPerProcessor = numFastaSeqs / processors;
					int startIndex =  pid * numSeqsPerProcessor;
					if(pid == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - pid * numSeqsPerProcessor; 	}
					
				
					//align your part
					driverMPI(startIndex, numSeqsPerProcessor, inMPI, outMPINewTax, outMPITempTax, MPIPos);
					
					if (m->control_pressed) {  outputTypes.clear(); MPI_File_close(&inMPI);  MPI_File_close(&outMPINewTax);   MPI_File_close(&outMPITempTax);  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());	} delete classify; return 0;  }
					
					for (int i = 1; i < processors; i++) {
						int done;
						MPI_Recv(&done, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
					}
				}else{ //you are a child process
					MPI_Recv(&numFastaSeqs, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
					MPIPos.resize(numFastaSeqs+1);
					MPI_Recv(&MPIPos[0], (numFastaSeqs+1), MPI_LONG, 0, tag, MPI_COMM_WORLD, &status);
					
					//figure out how many sequences you have to align
					numSeqsPerProcessor = numFastaSeqs / processors;
					int startIndex =  pid * numSeqsPerProcessor;
					if(pid == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - pid * numSeqsPerProcessor; 	}
					
					
					//align your part
					driverMPI(startIndex, numSeqsPerProcessor, inMPI, outMPINewTax, outMPITempTax, MPIPos);
					
					if (m->control_pressed) {  outputTypes.clear(); MPI_File_close(&inMPI);  MPI_File_close(&outMPINewTax);   MPI_File_close(&outMPITempTax);  delete classify; return 0;  }

					int done = 0;
					MPI_Send(&done, 1, MPI_INT, 0, tag, MPI_COMM_WORLD); 
				}
				
				//close files 
				MPI_File_close(&inMPI);
				MPI_File_close(&outMPINewTax);
				MPI_File_close(&outMPITempTax);
				MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
				
#else
		
			vector<unsigned long int> positions = m->divideFile(fastaFileNames[s], processors);
				
			for (int i = 0; i < (positions.size()-1); i++) {
				lines.push_back(new linePair(positions[i], positions[(i+1)]));
			}	
			
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			if(processors == 1){
				numFastaSeqs = driver(lines[0], newTaxonomyFile, tempTaxonomyFile, fastaFileNames[s]);
			}
			else{
				processIDS.resize(0);
				
				numFastaSeqs = createProcesses(newTaxonomyFile, tempTaxonomyFile, fastaFileNames[s]); 
				
			}
	#else
			numFastaSeqs = driver(lines[0], newTaxonomyFile, tempTaxonomyFile, fastaFileNames[s]);
	#endif	
#endif

		m->mothurOutEndLine();
		m->mothurOut("It took " + toString(time(NULL) - start) + " secs to classify " + toString(numFastaSeqs) + " sequences."); m->mothurOutEndLine(); m->mothurOutEndLine();
		start = time(NULL);


		#ifdef USE_MPI	
			if (pid == 0) {  //this part does not need to be paralellized
			
				if(namefile != "") { m->mothurOut("Reading " + namefileNames[s] + "..."); cout.flush();  MPIReadNamesFile(namefileNames[s]);  m->mothurOut("  Done."); m->mothurOutEndLine(); }
		#else
			//read namefile
			if(namefile != "") {
			
			    m->mothurOut("Reading " + namefileNames[s] + "..."); cout.flush();
				
				nameMap.clear(); //remove old names
				
				ifstream inNames;
				m->openInputFile(namefileNames[s], inNames);
				
				string firstCol, secondCol;
				while(!inNames.eof()) {
					inNames >> firstCol >> secondCol; m->gobble(inNames);
					
					vector<string> temp;
					m->splitAtComma(secondCol, temp);
			
					nameMap[firstCol] = temp;  
				}
				inNames.close();
				
				m->mothurOut("  Done."); m->mothurOutEndLine();
			}
		#endif

			string group = "";
			if (groupfile != "") {  group = groupfileNames[s]; }
			
			PhyloSummary taxaSum(taxonomyFileName, group);
			
			if (m->control_pressed) { outputTypes.clear();  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());	} delete classify; return 0; }
		
			if (namefile == "") {  taxaSum.summarize(tempTaxonomyFile);  }
			else {
				ifstream in;
				m->openInputFile(tempTaxonomyFile, in);
				
				//read in users taxonomy file and add sequences to tree
				string name, taxon;
				
				while(!in.eof()){
					in >> name >> taxon; m->gobble(in);
					
					itNames = nameMap.find(name);
		
					if (itNames == nameMap.end()) { 
						m->mothurOut(name + " is not in your name file please correct."); m->mothurOutEndLine(); exit(1);
					}else{
						for (int i = 0; i < itNames->second.size(); i++) { 
							taxaSum.addSeqToTree(itNames->second[i], taxon);  //add it as many times as there are identical seqs
						}
						itNames->second.clear();
						nameMap.erase(itNames->first);
					}
				}
				in.close();
			}
			remove(tempTaxonomyFile.c_str());
			
			if (m->control_pressed) {  outputTypes.clear(); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());	} delete classify; return 0; }
			
			//print summary file
			ofstream outTaxTree;
			m->openOutputFile(taxSummary, outTaxTree);
			taxaSum.print(outTaxTree);
			outTaxTree.close();
			
			//output taxonomy with the unclassified bins added
			ifstream inTax;
			m->openInputFile(newTaxonomyFile, inTax);
			
			ofstream outTax;
			string unclass = newTaxonomyFile + ".unclass.temp";
			m->openOutputFile(unclass, outTax);
			
			//get maxLevel from phylotree so you know how many 'unclassified's to add
			int maxLevel = taxaSum.getMaxLevel();
							
			//read taxfile - this reading and rewriting is done to preserve the confidence scores.
			string name, taxon;
			while (!inTax.eof()) {
				if (m->control_pressed) { outputTypes.clear();  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());	} remove(unclass.c_str()); delete classify; return 0; }

				inTax >> name >> taxon; m->gobble(inTax);
				
				string newTax = addUnclassifieds(taxon, maxLevel);
				
				outTax << name << '\t' << newTax << endl;
			}
			inTax.close();	
			outTax.close();
			
			remove(newTaxonomyFile.c_str());
			rename(unclass.c_str(), newTaxonomyFile.c_str());
			
			m->mothurOutEndLine();
			m->mothurOut("It took " + toString(time(NULL) - start) + " secs to create the summary file for " + toString(numFastaSeqs) + " sequences."); m->mothurOutEndLine(); m->mothurOutEndLine();
			
			#ifdef USE_MPI	
				}
			#endif

			m->mothurOutEndLine();
			m->mothurOut("Output File Names: "); m->mothurOutEndLine();
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
			m->mothurOutEndLine();
		}
		
		//set taxonomy file as new current taxonomyfile
		string current = "";
		itTypes = outputTypes.find("taxonomy");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setTaxonomyFile(current); }
		}
		
		delete classify;
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "execute");
		exit(1);
	}
}

/**************************************************************************************************/
string ClassifySeqsCommand::addUnclassifieds(string tax, int maxlevel) {
	try{
		string newTax, taxon;
		int level = 0;
		
		//keep what you have counting the levels
		while (tax.find_first_of(';') != -1) {
			//get taxon
			taxon = tax.substr(0,tax.find_first_of(';'))+';';
			tax = tax.substr(tax.find_first_of(';')+1, tax.length());
			newTax += taxon;
			level++;
		}
		
		//add "unclassified" until you reach maxLevel
		while (level < maxlevel) {
			newTax += "unclassified;";
			level++;
		}
		
		return newTax;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "addUnclassifieds");
		exit(1);
	}
}

/**************************************************************************************************/

int ClassifySeqsCommand::createProcesses(string taxFileName, string tempTaxFile, string filename) {
	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 1;
		int num = 0;
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driver(lines[process], taxFileName + toString(getpid()) + ".temp", tempTaxFile + toString(getpid()) + ".temp", filename);
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = filename + toString(getpid()) + ".num.temp";
				m->openOutputFile(tempFile, out);
				out << num << endl;
				out.close();

				exit(0);
			}else { 
				m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
				for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
				exit(0);
			}
		}
		
		//parent does its part
		num = driver(lines[0], taxFileName, tempTaxFile, filename);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		for (int i = 0; i < processIDS.size(); i++) {
			ifstream in;
			string tempFile =  filename + toString(processIDS[i]) + ".num.temp";
			m->openInputFile(tempFile, in);
			if (!in.eof()) { int tempNum = 0; in >> tempNum; num += tempNum; }
			in.close(); remove(tempFile.c_str());
		}
		
		for(int i=0;i<processIDS.size();i++){
			appendTaxFiles((taxFileName + toString(processIDS[i]) + ".temp"), taxFileName);
			appendTaxFiles((tempTaxFile + toString(processIDS[i]) + ".temp"), tempTaxFile);
			remove((taxFileName + toString(processIDS[i]) + ".temp").c_str());
			remove((tempTaxFile + toString(processIDS[i]) + ".temp").c_str());
		}
		
		return num;
#endif		
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "createProcesses");
		exit(1);
	}
}
/**************************************************************************************************/

void ClassifySeqsCommand::appendTaxFiles(string temp, string filename) {
	try{
		
		ofstream output;
		ifstream input;
		m->openOutputFileAppend(filename, output);
		m->openInputFile(temp, input);
		
		while(char c = input.get()){
			if(input.eof())		{	break;			}
			else				{	output << c;	}
		}
		
		input.close();
		output.close();
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "appendTaxFiles");
		exit(1);
	}
}

//**********************************************************************************************************************

int ClassifySeqsCommand::driver(linePair* filePos, string taxFName, string tempTFName, string filename){
	try {
		ofstream outTax;
		m->openOutputFile(taxFName, outTax);
		
		ofstream outTaxSimple;
		m->openOutputFile(tempTFName, outTaxSimple);
	
		ifstream inFASTA;
		m->openInputFile(filename, inFASTA);
		
		string taxonomy;

		inFASTA.seekg(filePos->start);

		bool done = false;
		int count = 0;
		
		while (!done) {
			if (m->control_pressed) { return 0; }
		
			Sequence* candidateSeq = new Sequence(inFASTA); m->gobble(inFASTA);
			
			if (candidateSeq->getName() != "") {
			
				taxonomy = classify->getTaxonomy(candidateSeq);
				
				if (m->control_pressed) { delete candidateSeq; return 0; }

				if (taxonomy != "bad seq") {
					//output confidence scores or not
					if (probs) {
						outTax << candidateSeq->getName() << '\t' << taxonomy << endl;
					}else{
						outTax << candidateSeq->getName() << '\t' << classify->getSimpleTax() << endl;
					}
					
					outTaxSimple << candidateSeq->getName() << '\t' << classify->getSimpleTax() << endl;
				}
				count++;
			}
			delete candidateSeq;
			
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				unsigned long int pos = inFASTA.tellg();
				if ((pos == -1) || (pos >= filePos->end)) { break; }
			#else
				if (inFASTA.eof()) { break; }
			#endif
			
			//report progress
			if((count) % 100 == 0){	m->mothurOut("Processing sequence: " + toString(count)); m->mothurOutEndLine();		}
			
		}
		//report progress
		if((count) % 100 != 0){	m->mothurOut("Processing sequence: " + toString(count)); m->mothurOutEndLine();		}
			
		inFASTA.close();
		outTax.close();
		outTaxSimple.close();
		
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
#ifdef USE_MPI
int ClassifySeqsCommand::driverMPI(int start, int num, MPI_File& inMPI, MPI_File& newFile, MPI_File& tempFile, vector<unsigned long int>& MPIPos){
	try {
		MPI_Status statusNew; 
		MPI_Status statusTemp; 
		MPI_Status status; 
		
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
	
		string taxonomy;
		string outputString;

		for(int i=0;i<num;i++){
		
			if (m->control_pressed) { return 0; }
		
			//read next sequence
			int length = MPIPos[start+i+1] - MPIPos[start+i];
			char* buf4 = new char[length];
			MPI_File_read_at(inMPI, MPIPos[start+i], buf4, length, MPI_CHAR, &status);
			
			string tempBuf = buf4;
			if (tempBuf.length() > length) { tempBuf = tempBuf.substr(0, length);  }
			istringstream iss (tempBuf,istringstream::in);
			delete buf4;

			Sequence* candidateSeq = new Sequence(iss);
			
			if (candidateSeq->getName() != "") {
				taxonomy = classify->getTaxonomy(candidateSeq);
				
				if (taxonomy != "bad seq") {
					//output confidence scores or not
					if (probs) {
						outputString =  candidateSeq->getName() + "\t" + taxonomy + "\n";
					}else{
						outputString =  candidateSeq->getName() + "\t" + classify->getSimpleTax() + "\n";
					}
					
					int length = outputString.length();
					char* buf2 = new char[length];
					memcpy(buf2, outputString.c_str(), length);
				
					MPI_File_write_shared(newFile, buf2, length, MPI_CHAR, &statusNew);
					delete buf2;

					outputString =  candidateSeq->getName() + "\t" + classify->getSimpleTax() + "\n";
					length = outputString.length();
					char* buf = new char[length];
					memcpy(buf, outputString.c_str(), length);
				
					MPI_File_write_shared(tempFile, buf, length, MPI_CHAR, &statusTemp);
					delete buf;
				}
			}				
			delete candidateSeq;
			
			if((i+1) % 100 == 0){	cout << "Classifying sequence " << (i+1) << endl;	}
		}
		
		if(num % 100 != 0){	cout << "Classifying sequence " << (num) << endl;	}
		
		
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "driverMPI");
		exit(1);
	}
}

//**********************************************************************************************************************
int ClassifySeqsCommand::MPIReadNamesFile(string nameFilename){
	try {
	
		nameMap.clear(); //remove old names
		
		MPI_File inMPI;
		MPI_Offset size;
		MPI_Status status;

		//char* inFileName = new char[nameFilename.length()];
		//memcpy(inFileName, nameFilename.c_str(), nameFilename.length());
		
		char inFileName[1024];
		strcpy(inFileName, nameFilename.c_str());

		MPI_File_open(MPI_COMM_WORLD, inFileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &inMPI);  
		MPI_File_get_size(inMPI, &size);
		//delete inFileName;

		char* buffer = new char[size];
		MPI_File_read(inMPI, buffer, size, MPI_CHAR, &status);

		string tempBuf = buffer;
		if (tempBuf.length() > size) { tempBuf = tempBuf.substr(0, size);  }
		istringstream iss (tempBuf,istringstream::in);
		delete buffer;
		
		string firstCol, secondCol;
		while(!iss.eof()) {
			iss >> firstCol >> secondCol; m->gobble(iss);
			
			vector<string> temp;
			m->splitAtComma(secondCol, temp);
			
			nameMap[firstCol] = temp;  
		}
	
		MPI_File_close(&inMPI);
		
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "MPIReadNamesFile");
		exit(1);
	}
}
#endif
/**************************************************************************************************/
