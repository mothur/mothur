/*
 *  classifyseqscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 11/2/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "classifyseqscommand.h"



//**********************************************************************************************************************
vector<string> ClassifySeqsCommand::setParameters(){	
	try {
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(ptaxonomy);
		CommandParameter ptemplate("reference", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(ptemplate);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pfasta);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none",false,false); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none",false,false); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none",false,false); parameters.push_back(pgroup);

		CommandParameter psearch("search", "Multiple", "kmer-blast-suffix-distance", "kmer", "", "", "",false,false); parameters.push_back(psearch);
		CommandParameter pksize("ksize", "Number", "", "8", "", "", "",false,false); parameters.push_back(pksize);
		CommandParameter pmethod("method", "Multiple", "bayesian-knn", "bayesian", "", "", "",false,false); parameters.push_back(pmethod);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "",false,false); parameters.push_back(pprocessors);
		CommandParameter pmatch("match", "Number", "", "1.0", "", "", "",false,false); parameters.push_back(pmatch);
		CommandParameter pmismatch("mismatch", "Number", "", "-1.0", "", "", "",false,false); parameters.push_back(pmismatch);
		CommandParameter pgapopen("gapopen", "Number", "", "-2.0", "", "", "",false,false); parameters.push_back(pgapopen);
		CommandParameter pgapextend("gapextend", "Number", "", "-1.0", "", "", "",false,false); parameters.push_back(pgapextend);
		//CommandParameter pflip("flip", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pflip);
		CommandParameter pcutoff("cutoff", "Number", "", "0", "", "", "",false,true); parameters.push_back(pcutoff);
		CommandParameter pprobs("probs", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(pprobs);
		CommandParameter piters("iters", "Number", "", "100", "", "", "",false,true); parameters.push_back(piters);
		CommandParameter psave("save", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(psave);
        CommandParameter pshortcuts("shortcuts", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(pshortcuts);
		CommandParameter pnumwanted("numwanted", "Number", "", "10", "", "", "",false,true); parameters.push_back(pnumwanted);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ClassifySeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The classify.seqs command reads a fasta file containing sequences and creates a .taxonomy file and a .tax.summary file.\n";
		helpString += "The classify.seqs command parameters are reference, fasta, name, group, count, search, ksize, method, taxonomy, processors, match, mismatch, gapopen, gapextend, numwanted and probs.\n";
		helpString += "The reference, fasta and taxonomy parameters are required. You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amzon.fasta \n";
		helpString += "The search parameter allows you to specify the method to find most similar template.  Your options are: suffix, kmer, blast and distance. The default is kmer.\n";
		helpString += "The name parameter allows you add a names file with your fasta file, if you enter multiple fasta files, you must enter matching names files for them.\n";
		helpString += "The group parameter allows you add a group file so you can have the summary totals broken up by group.\n";
        helpString += "The count parameter allows you add a count file so you can have the summary totals broken up by group.\n";
		helpString += "The method parameter allows you to specify classification method to use.  Your options are: bayesian and knn. The default is bayesian.\n";
		helpString += "The ksize parameter allows you to specify the kmer size for finding most similar template to candidate.  The default is 8.\n";
		helpString += "The processors parameter allows you to specify the number of processors to use. The default is 1.\n";
#ifdef USE_MPI
		helpString += "When using MPI, the processors parameter is set to the number of MPI processes running. \n";
#endif
		helpString += "If the save parameter is set to true the reference sequences will be saved in memory, to clear them later you can use the clear.memory command. Default=f.";
		helpString += "The match parameter allows you to specify the bonus for having the same base. The default is 1.0.\n";
		helpString += "The mistmatch parameter allows you to specify the penalty for having different bases.  The default is -1.0.\n";
		helpString += "The gapopen parameter allows you to specify the penalty for opening a gap in an alignment. The default is -2.0.\n";
		helpString += "The gapextend parameter allows you to specify the penalty for extending a gap in an alignment.  The default is -1.0.\n";
		helpString += "The numwanted parameter allows you to specify the number of sequence matches you want with the knn method.  The default is 10.\n";
		helpString += "The cutoff parameter allows you to specify a bootstrap confidence threshold for your taxonomy.  The default is 0.\n";
		helpString += "The probs parameter shuts off the bootstrapping results for the bayesian method. The default is true, meaning you want the bootstrapping to be shown.\n";
		helpString += "The iters parameter allows you to specify how many iterations to do when calculating the bootstrap confidence score for your taxonomy with the bayesian method.  The default is 100.\n";
		//helpString += "The flip parameter allows you shut off mothur's   The default is T.\n";
		helpString += "The classify.seqs command should be in the following format: \n";
		helpString += "classify.seqs(reference=yourTemplateFile, fasta=yourFastaFile, method=yourClassificationMethod, search=yourSearchmethod, ksize=yourKmerSize, taxonomy=yourTaxonomyFile, processors=yourProcessors) \n";
		helpString += "Example classify.seqs(fasta=amazon.fasta, reference=core.filtered, method=knn, search=gotoh, ksize=8, processors=2)\n";
		helpString += "The .taxonomy file consists of 2 columns: 1 = your sequence name, 2 = the taxonomy for your sequence. \n";
		helpString += "The .tax.summary is a summary of the different taxonomies represented in your fasta file. \n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ClassifySeqsCommand::getOutputFileNameTag(string type, string inputName=""){	
	try {
        string outputFileName = "";
		map<string, vector<string> >::iterator it;
        
        //is this a type this command creates
        it = outputTypes.find(type);
        if (it == outputTypes.end()) {  m->mothurOut("[ERROR]: this command doesn't create a " + type + " output file.\n"); }
        else {
            if (type == "taxonomy") {  outputFileName =  "taxonomy"; }
            else if (type == "accnos") {  outputFileName =  "flip.accnos"; }
            else if (type == "taxsummary") {  outputFileName =  "tax.summary"; }
            else if (type == "matchdist") {  outputFileName =  "match.dist"; }
            else { m->mothurOut("[ERROR]: No definition for type " + type + " output file tag.\n"); m->control_pressed = true;  }
        }
        return outputFileName;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "getOutputFileNameTag");
		exit(1);
	}
}
//**********************************************************************************************************************
ClassifySeqsCommand::ClassifySeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["taxonomy"] = tempOutNames;
		outputTypes["accnos"] = tempOutNames;
		outputTypes["taxsummary"] = tempOutNames;
		outputTypes["matchdist"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "ClassifySeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
ClassifySeqsCommand::ClassifySeqsCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		rdb = ReferenceDB::getInstance(); hasName = false; hasCount=false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
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
			outputTypes["accnos"] = tempOutNames;
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("reference");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["reference"] = inputDir + it->second;		}
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
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
			}

			fastaFileName = validParameter.validFile(parameters, "fasta", false);
			if (fastaFileName == "not found") { 				
				//if there is a current fasta file, use it
				string filename = m->getFastaFile(); 
				if (filename != "") { fastaFileNames.push_back(filename); m->mothurOut("Using " + filename + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else { 
				m->splitAtDash(fastaFileName, fastaFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastaFileNames.size(); i++) {
					
					bool ignore = false;
					if (fastaFileNames[i] == "current") { 
						fastaFileNames[i] = m->getFastaFile(); 
						if (fastaFileNames[i] != "") {  m->mothurOut("Using " + fastaFileNames[i] + " as input file for the fasta parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current fastafile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							fastaFileNames.erase(fastaFileNames.begin()+i);
							i--;
						}
					}
					
					if (!ignore) {
						
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
						}else {
							m->setFastaFile(fastaFileNames[i]);
						}
					}
					
				}
				
				//make sure there is at least one valid file left
				if (fastaFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}

			namefile = validParameter.validFile(parameters, "name", false);
			if (namefile == "not found") { namefile = "";  }

			else { 
				m->splitAtDash(namefile, namefileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < namefileNames.size(); i++) {
					bool ignore = false;
					if (namefileNames[i] == "current") { 
						namefileNames[i] = m->getNameFile(); 
						if (namefileNames[i] != "") {  m->mothurOut("Using " + namefileNames[i] + " as input file for the name parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current namefile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							namefileNames.erase(namefileNames.begin()+i);
							i--;
						}
					}
					
					if (!ignore) {
						
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
						}else {
							m->setNameFile(namefileNames[i]);
						}
					}
				}
			}
            
            if (namefileNames.size() != 0) { hasName = true; }
            
			if (namefile != "") {
				if (namefileNames.size() != fastaFileNames.size()) { abort = true; m->mothurOut("If you provide a name file, you must have one for each fasta file."); m->mothurOutEndLine(); }
			}
			
            //check for required parameters
			countfile = validParameter.validFile(parameters, "count", false);
			if (countfile == "not found") { 
                countfile = "";  
			}else { 
				m->splitAtDash(countfile, countfileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < countfileNames.size(); i++) {
					
					bool ignore = false;
					if (countfileNames[i] == "current") { 
						countfileNames[i] = m->getCountTableFile(); 
						if (countfileNames[i] != "") {  m->mothurOut("Using " + countfileNames[i] + " as input file for the count parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current count file, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							countfileNames.erase(countfileNames.begin()+i);
							i--;
						}
					}
					
					if (!ignore) {
						
						if (inputDir != "") {
							string path = m->hasPath(countfileNames[i]);
							//if the user has not given a path then, add inputdir. else leave path alone.
							if (path == "") {	countfileNames[i] = inputDir + countfileNames[i];		}
						}
						
						int ableToOpen;
						ifstream in;
						
						ableToOpen = m->openInputFile(countfileNames[i], in, "noerror");
						
						//if you can't open it, try default location
						if (ableToOpen == 1) {
							if (m->getDefaultPath() != "") { //default path is set
								string tryPath = m->getDefaultPath() + m->getSimpleName(countfileNames[i]);
								m->mothurOut("Unable to open " + countfileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								countfileNames[i] = tryPath;
							}
						}
						
						if (ableToOpen == 1) {
							if (m->getOutputDir() != "") { //default path is set
								string tryPath = m->getOutputDir() + m->getSimpleName(countfileNames[i]);
								m->mothurOut("Unable to open " + countfileNames[i] + ". Trying output directory " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								countfileNames[i] = tryPath;
							}
						}
						
						in.close();
						
						if (ableToOpen == 1) { 
							m->mothurOut("Unable to open " + countfileNames[i] + ". It will be disregarded."); m->mothurOutEndLine(); 
							//erase from file list
							countfileNames.erase(countfileNames.begin()+i);
							i--;
						}else {
							m->setCountTableFile(countfileNames[i]);
						}
					}
				}
			}
            
            if (countfileNames.size() != 0) { hasCount = true; if (countfileNames.size() != fastaFileNames.size()) {m->mothurOut("If you provide a count file, you must have one for each fasta file."); m->mothurOutEndLine(); } }
            
			//make sure there is at least one valid file left
            if (hasName && hasCount) { m->mothurOut("[ERROR]: You must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }

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
					}else {
						m->setGroupFile(groupfileNames[i]);
					}
				}
			}

			if (groupfile != "") {
				if (groupfileNames.size() != fastaFileNames.size()) { abort = true; m->mothurOut("If you provide a group file, you must have one for each fasta file."); m->mothurOutEndLine(); }
                if (hasCount) { m->mothurOut("[ERROR]: You must enter ONLY ONE of the following: count or group."); m->mothurOutEndLine(); abort = true; }
			}else {
				for (int i = 0; i < fastaFileNames.size(); i++) {  groupfileNames.push_back("");  }
			}
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.validFile(parameters, "ksize", false);		if (temp == "not found"){	temp = "8";				}
			m->mothurConvert(temp, kmerSize); 
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors); 
			
			temp = validParameter.validFile(parameters, "save", false);			if (temp == "not found"){	temp = "f";				}
			save = m->isTrue(temp); 
			rdb->save = save; 
			if (save) { //clear out old references
				rdb->clearMemory();	
			}
			
			//this has to go after save so that if the user sets save=t and provides no reference we abort
			templateFileName = validParameter.validFile(parameters, "reference", true);
			if (templateFileName == "not found") { 
				//check for saved reference sequences
				if (rdb->referenceSeqs.size() != 0) {
					templateFileName = "saved";
				}else {
					m->mothurOut("[ERROR]: You don't have any saved reference sequences and the reference parameter is a required for the classify.seqs command."); 
					m->mothurOutEndLine();
					abort = true; 
				}
			}else if (templateFileName == "not open") { abort = true; }	
			else {	if (save) {	rdb->setSavedReference(templateFileName);	}	}
			
			//this has to go after save so that if the user sets save=t and provides no reference we abort
			taxonomyFileName = validParameter.validFile(parameters, "taxonomy", true);
			if (taxonomyFileName == "not found") { 
				//check for saved reference sequences
				if (rdb->wordGenusProb.size() != 0) {
					taxonomyFileName = "saved";
				}else {
					m->mothurOut("[ERROR]: You don't have any saved taxonomy information and the taxonomy parameter is a required for the classify.seqs command."); 
					m->mothurOutEndLine();
					abort = true; 
				}
			}else if (taxonomyFileName == "not open") { abort = true; }	
			else {	if (save) {	rdb->setSavedTaxonomy(taxonomyFileName);	}	}
			
			search = validParameter.validFile(parameters, "search", false);		if (search == "not found"){	search = "kmer";		}
			
			method = validParameter.validFile(parameters, "method", false);		if (method == "not found"){	method = "bayesian";	}
			
			temp = validParameter.validFile(parameters, "match", false);		if (temp == "not found"){	temp = "1.0";			}
			m->mothurConvert(temp, match);  
			
			temp = validParameter.validFile(parameters, "mismatch", false);		if (temp == "not found"){	temp = "-1.0";			}
			m->mothurConvert(temp, misMatch);  
			
			temp = validParameter.validFile(parameters, "gapopen", false);		if (temp == "not found"){	temp = "-2.0";			}
			m->mothurConvert(temp, gapOpen);  
			
			temp = validParameter.validFile(parameters, "gapextend", false);	if (temp == "not found"){	temp = "-1.0";			}
			m->mothurConvert(temp, gapExtend); 
			
			temp = validParameter.validFile(parameters, "numwanted", false);	if (temp == "not found"){	temp = "10";			}
			m->mothurConvert(temp, numWanted);
			
			temp = validParameter.validFile(parameters, "cutoff", false);		if (temp == "not found"){	temp = "0";				}
			m->mothurConvert(temp, cutoff);
			
			temp = validParameter.validFile(parameters, "probs", false);		if (temp == "not found"){	temp = "true";			}
			probs = m->isTrue(temp);
            
            temp = validParameter.validFile(parameters, "shortcuts", false);	if (temp == "not found"){	temp = "true";			}
			writeShortcuts = m->isTrue(temp);
			
			//temp = validParameter.validFile(parameters, "flip", false);			if (temp == "not found"){	temp = "T";				}
			//flip = m->isTrue(temp); 
			flip = true;
			
			temp = validParameter.validFile(parameters, "iters", false);		if (temp == "not found") { temp = "100";			}
			m->mothurConvert(temp, iters); 

			
			if ((method == "bayesian") && (search != "kmer"))  { 
				m->mothurOut("The bayesian method requires the kmer search." + search + "will be disregarded." ); m->mothurOutEndLine();
				search = "kmer";
			}
			
            if (!abort) {
                if (!hasCount) {
                    if (namefileNames.size() == 0){
                        if (fastaFileNames.size() != 0) {
                            vector<string> files; files.push_back(fastaFileNames[fastaFileNames.size()-1]); 
                            parser.getNameFile(files);
                        }
                    }
                }
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

int ClassifySeqsCommand::execute(){
	try {
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
		if(method == "bayesian"){	classify = new Bayesian(taxonomyFileName, templateFileName, search, kmerSize, cutoff, iters, rand(), flip, writeShortcuts);		}
		else if(method == "knn"){	classify = new Knn(taxonomyFileName, templateFileName, search, kmerSize, gapOpen, gapExtend, match, misMatch, numWanted, rand());				}
		else {
			m->mothurOut(search + " is not a valid method option. I will run the command using bayesian.");
			m->mothurOutEndLine();
			classify = new Bayesian(taxonomyFileName, templateFileName, search, kmerSize, cutoff, iters, rand(), flip, writeShortcuts);	
		}
		
		if (m->control_pressed) { delete classify; return 0; }
				
		for (int s = 0; s < fastaFileNames.size(); s++) {
		
			m->mothurOut("Classifying sequences from " + fastaFileNames[s] + " ..." ); m->mothurOutEndLine();
			
			string baseTName = taxonomyFileName;
			if (taxonomyFileName == "saved") {baseTName = rdb->getSavedTaxonomy();	}
			
            //set rippedTaxName to 
			string RippedTaxName = "";
            bool foundDot = false;
            for (int i = baseTName.length()-1; i >= 0; i--) {
                if (foundDot && (baseTName[i] != '.')) {  RippedTaxName = baseTName[i] + RippedTaxName; }
                else if (foundDot && (baseTName[i] == '.')) {  break; }
                else if (!foundDot && (baseTName[i] == '.')) {  foundDot = true; }
            }
            if (RippedTaxName != "") { RippedTaxName +=  "."; }   
          
			if (outputDir == "") { outputDir += m->hasPath(fastaFileNames[s]); }
			string newTaxonomyFile = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + RippedTaxName + getOutputFileNameTag("taxonomy");
			string newaccnosFile = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + RippedTaxName + getOutputFileNameTag("accnos");
			string tempTaxonomyFile = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + "taxonomy.temp";
			string taxSummary = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + RippedTaxName + getOutputFileNameTag("taxsummary");
			
			if ((method == "knn") && (search == "distance")) { 
				string DistName = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + getOutputFileNameTag("matchdist");
				classify->setDistName(DistName);  outputNames.push_back(DistName); outputTypes["matchdist"].push_back(DistName);
			}
			
			outputNames.push_back(newTaxonomyFile); outputTypes["taxonomy"].push_back(newTaxonomyFile);
			outputNames.push_back(newaccnosFile); outputTypes["accnos"].push_back(newaccnosFile);
			outputNames.push_back(taxSummary);	outputTypes["taxsummary"].push_back(taxSummary);
			
			int start = time(NULL);
			int numFastaSeqs = 0;
			for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
			
#ifdef USE_MPI	
				int pid, numSeqsPerProcessor; 
				int tag = 2001;
				vector<unsigned long long> MPIPos;
				
				MPI_Status status; 
				MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
				MPI_Comm_size(MPI_COMM_WORLD, &processors); 

				MPI_File inMPI;
				MPI_File outMPINewTax;
				MPI_File outMPITempTax;
				MPI_File outMPIAcc;
							
				int outMode=MPI_MODE_CREATE|MPI_MODE_WRONLY; 
				int inMode=MPI_MODE_RDONLY; 
				
				char outNewTax[1024];
				strcpy(outNewTax, newTaxonomyFile.c_str());
				
				char outTempTax[1024];
				strcpy(outTempTax, tempTaxonomyFile.c_str());
			
				char outAcc[1024];
				strcpy(outAcc, newaccnosFile.c_str());
				
				char inFileName[1024];
				strcpy(inFileName, fastaFileNames[s].c_str());

				MPI_File_open(MPI_COMM_WORLD, inFileName, inMode, MPI_INFO_NULL, &inMPI);  //comm, filename, mode, info, filepointer
				MPI_File_open(MPI_COMM_WORLD, outNewTax, outMode, MPI_INFO_NULL, &outMPINewTax);
				MPI_File_open(MPI_COMM_WORLD, outTempTax, outMode, MPI_INFO_NULL, &outMPITempTax);
				MPI_File_open(MPI_COMM_WORLD, outAcc, outMode, MPI_INFO_NULL, &outMPIAcc);
				
				if (m->control_pressed) { outputTypes.clear(); MPI_File_close(&inMPI);  MPI_File_close(&outMPINewTax);  MPI_File_close(&outMPIAcc);   MPI_File_close(&outMPITempTax);  delete classify;  return 0;  }
				
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
					driverMPI(startIndex, numSeqsPerProcessor, inMPI, outMPINewTax, outMPITempTax, outMPIAcc, MPIPos);
					
					if (m->control_pressed) {  outputTypes.clear(); MPI_File_close(&inMPI);  MPI_File_close(&outMPINewTax); MPI_File_close(&outMPIAcc);   MPI_File_close(&outMPITempTax);  for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);	} delete classify; return 0;  }
					
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
					driverMPI(startIndex, numSeqsPerProcessor, inMPI, outMPINewTax, outMPITempTax, outMPIAcc, MPIPos);
					
					if (m->control_pressed) {  outputTypes.clear(); MPI_File_close(&inMPI);  MPI_File_close(&outMPINewTax);  MPI_File_close(&outMPIAcc);  MPI_File_close(&outMPITempTax);  delete classify; return 0;  }

					int done = 0;
					MPI_Send(&done, 1, MPI_INT, 0, tag, MPI_COMM_WORLD); 
				}
				
				//close files 
				MPI_File_close(&inMPI);
				MPI_File_close(&outMPINewTax);
				MPI_File_close(&outMPITempTax);
				MPI_File_close(&outMPIAcc); 
				MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
				
#else
		
			vector<unsigned long long> positions; 
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
			positions = m->divideFile(fastaFileNames[s], processors);
			for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(new linePair(positions[i], positions[(i+1)]));	}
#else
			if (processors == 1) {
				lines.push_back(new linePair(0, 1000));
			}else {
				positions = m->setFilePosFasta(fastaFileNames[s], numFastaSeqs); 
                if (positions.size() < processors) { processors = positions.size(); }
				
				//figure out how many sequences you have to process
				int numSeqsPerProcessor = numFastaSeqs / processors;
				for (int i = 0; i < processors; i++) {
					int startIndex =  i * numSeqsPerProcessor;
					if(i == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor; 	}
					lines.push_back(new linePair(positions[startIndex], numSeqsPerProcessor));
				}
			}
#endif
			if(processors == 1){
				numFastaSeqs = driver(lines[0], newTaxonomyFile, tempTaxonomyFile, newaccnosFile, fastaFileNames[s]);
			}else{
				numFastaSeqs = createProcesses(newTaxonomyFile, tempTaxonomyFile, newaccnosFile, fastaFileNames[s]); 
			}
#endif
			
			if (!m->isBlank(newaccnosFile)) { m->mothurOutEndLine(); m->mothurOut("[WARNING]: mothur suspects some of your sequences may be reversed, please check " + newaccnosFile + " for the list of the sequences."); m->mothurOutEndLine(); }

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
				m->readNames(namefileNames[s], nameMap);				
				m->mothurOut("  Done."); m->mothurOutEndLine();
			}
		#endif

                string group = "";
                GroupMap* groupMap = NULL;
                CountTable* ct = NULL;
                PhyloSummary* taxaSum;
                if (hasCount) { 
                    ct = new CountTable();
                    ct->readTable(countfileNames[s]);
                    taxaSum = new PhyloSummary(baseTName, ct);
                    taxaSum->summarize(tempTaxonomyFile);
                }else {
                    if (groupfile != "") {  group = groupfileNames[s]; groupMap = new GroupMap(group); groupMap->readMap(); }
                    
                    taxaSum = new PhyloSummary(baseTName, groupMap);
                    
                    if (m->control_pressed) { outputTypes.clear(); if (ct != NULL) { delete ct; }  if (groupMap != NULL) { delete groupMap; } delete taxaSum; for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);	} delete classify; return 0; }
                    
                    if (namefile == "") {  taxaSum->summarize(tempTaxonomyFile);  }
                    else {
                        ifstream in;
                        m->openInputFile(tempTaxonomyFile, in);
                        
                        //read in users taxonomy file and add sequences to tree
                        string name, taxon;
                        
                        while(!in.eof()){
                            if (m->control_pressed) { outputTypes.clear();  if (ct != NULL) { delete ct; } if (groupMap != NULL) { delete groupMap; } delete taxaSum; for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);	} delete classify; return 0; }
                            
                            in >> name >> taxon; m->gobble(in);
                            
                            itNames = nameMap.find(name);
                            
                            if (itNames == nameMap.end()) { 
                                m->mothurOut(name + " is not in your name file please correct."); m->mothurOutEndLine(); exit(1);
                            }else{
                                for (int i = 0; i < itNames->second.size(); i++) { 
                                    taxaSum->addSeqToTree(itNames->second[i], taxon);  //add it as many times as there are identical seqs
                                }
                                itNames->second.clear();
                                nameMap.erase(itNames->first);
                            }
                        }
                        in.close();
                    }
                }
			m->mothurRemove(tempTaxonomyFile);
			
			if (m->control_pressed) {  outputTypes.clear(); if (ct != NULL) { delete ct; } if (groupMap != NULL) { delete groupMap; } for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);	} delete classify; return 0; }
			
			//print summary file
			ofstream outTaxTree;
			m->openOutputFile(taxSummary, outTaxTree);
			taxaSum->print(outTaxTree);
			outTaxTree.close();
			
			//output taxonomy with the unclassified bins added
			ifstream inTax;
			m->openInputFile(newTaxonomyFile, inTax);
			
			ofstream outTax;
			string unclass = newTaxonomyFile + ".unclass.temp";
			m->openOutputFile(unclass, outTax);
			
			//get maxLevel from phylotree so you know how many 'unclassified's to add
			int maxLevel = taxaSum->getMaxLevel();
							
			//read taxfile - this reading and rewriting is done to preserve the confidence scores.
			string name, taxon;
			while (!inTax.eof()) {
				if (m->control_pressed) { outputTypes.clear(); if (ct != NULL) { delete ct; }  if (groupMap != NULL) { delete groupMap; } delete taxaSum; for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);	} m->mothurRemove(unclass); delete classify; return 0; }

				inTax >> name >> taxon; m->gobble(inTax);
				
				string newTax = addUnclassifieds(taxon, maxLevel);
				
				outTax << name << '\t' << newTax << endl;
			}
			inTax.close();	
			outTax.close();
			
            if (ct != NULL) { delete ct; }
            if (groupMap != NULL) { delete groupMap; } delete taxaSum;
			m->mothurRemove(newTaxonomyFile);
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
		
		current = "";
		itTypes = outputTypes.find("accnos");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setAccnosFile(current); }
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

int ClassifySeqsCommand::createProcesses(string taxFileName, string tempTaxFile, string accnos, string filename) {
	try {
		
		int num = 0;
		processIDS.clear();
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		int process = 1;
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driver(lines[process], taxFileName + toString(getpid()) + ".temp", tempTaxFile + toString(getpid()) + ".temp", accnos + toString(getpid()) + ".temp", filename);

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
		num = driver(lines[0], taxFileName, tempTaxFile, accnos, filename);
		
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
			in.close(); m->mothurRemove(m->getFullPathName(tempFile));
		}
#else
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the alignData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<classifyData*> pDataArray; 
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
		
		//Create processor worker threads.
		for( int i=0; i<processors-1; i++ ){
			// Allocate memory for thread data.
			string extension = "";
			if (i != 0) { extension = toString(i) + ".temp"; processIDS.push_back(i); }
			
			classifyData* tempclass = new classifyData((accnos + extension), probs, method, templateFileName, taxonomyFileName, (taxFileName + extension), (tempTaxFile + extension), filename, search, kmerSize, iters, numWanted, m, lines[i]->start, lines[i]->end, match, misMatch, gapOpen, gapExtend, cutoff, i, flip, writeShortcuts);
			pDataArray.push_back(tempclass);
			
			//MySeqSumThreadFunction is in header. It must be global or static to work with the threads.
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i] = CreateThread(NULL, 0, MyClassThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);  
			
		}
		
		//parent does its part
		num = driver(lines[processors-1], taxFileName + toString(processors-1) + ".temp", tempTaxFile + toString(processors-1) + ".temp", accnos + toString(processors-1) + ".temp", filename);
		processIDS.push_back((processors-1));
		
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
			num += pDataArray[i]->count;
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
		
	#endif	
        vector<string> nonBlankAccnosFiles;
		if (!(m->isBlank(accnos))) { nonBlankAccnosFiles.push_back(accnos); }
		else { m->mothurRemove(accnos); } //remove so other files can be renamed to it
        
		for(int i=0;i<processIDS.size();i++){
			m->appendFiles((taxFileName + toString(processIDS[i]) + ".temp"), taxFileName);
			m->appendFiles((tempTaxFile + toString(processIDS[i]) + ".temp"), tempTaxFile);
            if (!(m->isBlank(accnos + toString(processIDS[i]) + ".temp"))) {
				nonBlankAccnosFiles.push_back(accnos + toString(processIDS[i]) + ".temp");
			}else { m->mothurRemove((accnos + toString(processIDS[i]) + ".temp"));  }

			m->mothurRemove((m->getFullPathName(taxFileName) + toString(processIDS[i]) + ".temp"));
			m->mothurRemove((m->getFullPathName(tempTaxFile) + toString(processIDS[i]) + ".temp"));
		}
		
        //append accnos files
		if (nonBlankAccnosFiles.size() != 0) { 
			rename(nonBlankAccnosFiles[0].c_str(), accnos.c_str());
			
			for (int h=1; h < nonBlankAccnosFiles.size(); h++) {
				m->appendFiles(nonBlankAccnosFiles[h], accnos);
				m->mothurRemove(nonBlankAccnosFiles[h]);
			}
		}else { //recreate the accnosfile if needed
			ofstream out;
			m->openOutputFile(accnos, out);
			out.close();
		}

		return num;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "createProcesses");
		exit(1);
	}
}
//**********************************************************************************************************************

int ClassifySeqsCommand::driver(linePair* filePos, string taxFName, string tempTFName, string accnos, string filename){
	try {
		ofstream outTax;
		m->openOutputFile(taxFName, outTax);
		
		ofstream outTaxSimple;
		m->openOutputFile(tempTFName, outTaxSimple);
		
		ofstream outAcc;
		m->openOutputFile(accnos, outAcc);
	
		ifstream inFASTA;
		m->openInputFile(filename, inFASTA);
		
		string taxonomy;

		inFASTA.seekg(filePos->start);

		bool done = false;
		int count = 0;
		
		while (!done) {
			if (m->control_pressed) { 
				inFASTA.close();
				outTax.close();
				outTaxSimple.close();
				outAcc.close(); return 0; }
		
			Sequence* candidateSeq = new Sequence(inFASTA); m->gobble(inFASTA);
			
			if (candidateSeq->getName() != "") {
			
				taxonomy = classify->getTaxonomy(candidateSeq);
				
				if (m->control_pressed) { delete candidateSeq; return 0; }
				
				if (taxonomy == "unknown;") { m->mothurOut("[WARNING]: " + candidateSeq->getName() + " could not be classified. You can use the remove.lineage command with taxon=unknown; to remove such sequences."); m->mothurOutEndLine(); }
				
				//output confidence scores or not
				if (probs) {
					outTax << candidateSeq->getName() << '\t' << taxonomy << endl;
				}else{
					outTax << candidateSeq->getName() << '\t' << classify->getSimpleTax() << endl;
				}
				
				if (classify->getFlipped()) { outAcc << candidateSeq->getName() << endl; }
				
				outTaxSimple << candidateSeq->getName() << '\t' << classify->getSimpleTax() << endl;
				
				count++;
			}
			delete candidateSeq;
			
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				unsigned long long pos = inFASTA.tellg();
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
		outAcc.close();
		
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
#ifdef USE_MPI
int ClassifySeqsCommand::driverMPI(int start, int num, MPI_File& inMPI, MPI_File& newFile, MPI_File& tempFile, MPI_File& accFile, vector<unsigned long long>& MPIPos){
	try {
		MPI_Status statusNew; 
		MPI_Status statusTemp; 
		MPI_Status statusAcc; 
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
				
				if (taxonomy == "unknown;") { m->mothurOut("[WARNING]: " + candidateSeq->getName() + " could not be classified. You can use the remove.lineage command with taxon=unknown; to remove such sequences."); m->mothurOutEndLine(); }
				
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
				
				if (classify->getFlipped()) { 
					outputString =  candidateSeq->getName() + "\n";
					length = outputString.length();
					char* buf3 = new char[length];
					memcpy(buf3, outputString.c_str(), length);
					
					MPI_File_write_shared(accFile, buf3, length, MPI_CHAR, &statusAcc);
					delete buf3;
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
