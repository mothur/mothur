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
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(ptaxonomy);
		CommandParameter ptemplate("reference", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(ptemplate);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","taxonomy",false,true,true); parameters.push_back(pfasta);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
        CommandParameter poutput("output", "Multiple", "simple-detail", "detail", "", "", "","",false,false, true); parameters.push_back(poutput);
		CommandParameter psearch("search", "Multiple", "kmer-blast-suffix-distance-align", "kmer", "", "", "","",false,false); parameters.push_back(psearch);
		CommandParameter pksize("ksize", "Number", "", "8", "", "", "","",false,false); parameters.push_back(pksize);
		CommandParameter pmethod("method", "Multiple", "wang-knn-zap", "wang", "", "", "","",false,false); parameters.push_back(pmethod);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pmatch("match", "Number", "", "1.0", "", "", "","",false,false); parameters.push_back(pmatch);
        CommandParameter pprintlevel("printlevel", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pprintlevel);
		CommandParameter pmismatch("mismatch", "Number", "", "-1.0", "", "", "","",false,false); parameters.push_back(pmismatch);
		CommandParameter pgapopen("gapopen", "Number", "", "-2.0", "", "", "","",false,false); parameters.push_back(pgapopen);
		CommandParameter pgapextend("gapextend", "Number", "", "-1.0", "", "", "","",false,false); parameters.push_back(pgapextend);
		CommandParameter pcutoff("cutoff", "Number", "", "80", "", "", "","",false,true); parameters.push_back(pcutoff);
		CommandParameter pprobs("probs", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pprobs);
		CommandParameter piters("iters", "Number", "", "100", "", "", "","",false,true); parameters.push_back(piters);
        CommandParameter pshortcuts("shortcuts", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pshortcuts);
        CommandParameter prelabund("relabund", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(prelabund);
		CommandParameter pnumwanted("numwanted", "Number", "", "10", "", "", "","",false,true); parameters.push_back(pnumwanted);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
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
		helpString += "The classify.seqs command parameters are reference, fasta, name, group, count, search, ksize, method, taxonomy, processors, match, mismatch, gapopen, gapextend, numwanted, relabund and probs.\n";
		helpString += "The reference, fasta and taxonomy parameters are required. You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amzon.fasta \n";
		helpString += "The search parameter allows you to specify the method to find most similar template.  Your options are: suffix, kmer, blast, align and distance. The default is kmer.\n";
		helpString += "The name parameter allows you add a names file with your fasta file, if you enter multiple fasta files, you must enter matching names files for them.\n";
		helpString += "The group parameter allows you add a group file so you can have the summary totals broken up by group.\n";
        helpString += "The count parameter allows you add a count file so you can have the summary totals broken up by group.\n";
		helpString += "The method parameter allows you to specify classification method to use.  Your options are: wang, knn and zap. The default is wang.\n";
		helpString += "The ksize parameter allows you to specify the kmer size for finding most similar template to candidate.  The default is 8.\n";
		helpString += "The processors parameter allows you to specify the number of processors to use. The default is 1.\n";
		helpString += "The match parameter allows you to specify the bonus for having the same base. The default is 1.0.\n";
		helpString += "The mistmatch parameter allows you to specify the penalty for having different bases.  The default is -1.0.\n";
		helpString += "The gapopen parameter allows you to specify the penalty for opening a gap in an alignment. The default is -2.0.\n";
		helpString += "The gapextend parameter allows you to specify the penalty for extending a gap in an alignment.  The default is -1.0.\n";
		helpString += "The numwanted parameter allows you to specify the number of sequence matches you want with the knn method.  The default is 10.\n";
		helpString += "The cutoff parameter allows you to specify a bootstrap confidence threshold for your taxonomy.  The default is 80.\n";
		helpString += "The probs parameter shuts off the bootstrapping results for the wang and zap method. The default is true, meaning you want the bootstrapping to be shown.\n";
        helpString += "The relabund parameter allows you to indicate you want the summary file values to be relative abundances rather than raw abundances. Default=F. \n";
		helpString += "The iters parameter allows you to specify how many iterations to do when calculating the bootstrap confidence score for your taxonomy with the wang method.  The default is 100.\n";
		helpString += "The output parameter allows you to specify format of your summary file. Options are simple and detail. The default is detail.\n";
        helpString += "The printlevel parameter allows you to specify taxlevel of your summary file to print to. Options are 1 to the maz level in the file.  The default is -1, meaning max level.  If you select a level greater than the level your sequences classify to, mothur will print to the level your max level. \n";
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
string ClassifySeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "taxonomy") {  pattern = "[filename],[tag],[tag2],taxonomy"; } 
        else if (type == "taxsummary") {  pattern = "[filename],[tag],[tag2],tax.summary"; } 
        else if (type == "accnos") {  pattern =  "[filename],[tag],[tag2],flip.accnos"; }
        else if (type == "matchdist") {  pattern =  "[filename],[tag],[tag2],match.dist"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ClassifySeqsCommand", "getOutputPattern");
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
		hasName = false; hasCount=false;
		
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
					
					bool ignore = false;
					if (groupfileNames[i] == "current") { 
						groupfileNames[i] = m->getGroupFile(); 
						if (groupfileNames[i] != "") {  m->mothurOut("Using " + groupfileNames[i] + " as input file for the group parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current group file, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							groupfileNames.erase(groupfileNames.begin()+i);
							i--;
						}
					}
					
					if (!ignore) {
						
						if (inputDir != "") {
							string path = m->hasPath(groupfileNames[i]);
                            //cout << path << '\t' << inputDir << endl;
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
							m->mothurOut("Unable to open " + groupfileNames[i] + ". It will be disregarded."); m->mothurOutEndLine(); 
							//erase from file list
							groupfileNames.erase(groupfileNames.begin()+i);
							i--;
						}else {
							m->setGroupFile(groupfileNames[i]);
						}
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
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors); 

			//this has to go after save so that if the user sets save=t and provides no reference we abort
			templateFileName = validParameter.validFile(parameters, "reference", true);
			if (templateFileName == "not found") {
					m->mothurOut("[ERROR]: The reference parameter is a required for the classify.seqs command.\n"); abort = true;
			}else if (templateFileName == "not open") { abort = true; }
			
			
			//this has to go after save so that if the user sets save=t and provides no reference we abort
			taxonomyFileName = validParameter.validFile(parameters, "taxonomy", true);
			if (taxonomyFileName == "not found") {  m->mothurOut("[ERROR]: The taxonomy parameter is a required for the classify.seqs command.\n"); abort = true;
			}else if (taxonomyFileName == "not open") { abort = true; }
			
			search = validParameter.validFile(parameters, "search", false);		if (search == "not found"){	search = "kmer";		}
			
			method = validParameter.validFile(parameters, "method", false);		if (method == "not found"){	method = "wang";	}
            
            temp = validParameter.validFile(parameters, "ksize", false);		if (temp == "not found"){	
                temp = "8";	
                if (method == "zap") { temp = "7"; }
            }
			m->mothurConvert(temp, kmerSize); 
			
			temp = validParameter.validFile(parameters, "match", false);		if (temp == "not found"){	temp = "1.0";			}
			m->mothurConvert(temp, match);
            
            temp = validParameter.validFile(parameters, "printlevel", false);		if (temp == "not found"){	temp = "-1";		}
            m->mothurConvert(temp, printlevel);
			
			temp = validParameter.validFile(parameters, "mismatch", false);		if (temp == "not found"){	temp = "-1.0";			}
			m->mothurConvert(temp, misMatch);  
			
			temp = validParameter.validFile(parameters, "gapopen", false);		if (temp == "not found"){	temp = "-2.0";			}
			m->mothurConvert(temp, gapOpen);  
			
			temp = validParameter.validFile(parameters, "gapextend", false);	if (temp == "not found"){	temp = "-1.0";			}
			m->mothurConvert(temp, gapExtend); 
			
			temp = validParameter.validFile(parameters, "numwanted", false);	if (temp == "not found"){	temp = "10";			}
			m->mothurConvert(temp, numWanted);
			
			temp = validParameter.validFile(parameters, "cutoff", false);		if (temp == "not found"){	temp = "80";				}
			m->mothurConvert(temp, cutoff);
			
			temp = validParameter.validFile(parameters, "probs", false);		if (temp == "not found"){	temp = "true";			}
			probs = m->isTrue(temp);
            
            temp = validParameter.validFile(parameters, "relabund", false);		if (temp == "not found"){	temp = "false";			}
			relabund = m->isTrue(temp);
            
            temp = validParameter.validFile(parameters, "shortcuts", false);	if (temp == "not found"){	temp = "true";			}
			writeShortcuts = m->isTrue(temp);
			
			//temp = validParameter.validFile(parameters, "flip", false);			if (temp == "not found"){	temp = "T";				}
			//flip = m->isTrue(temp); 
			flip = true;
			
			temp = validParameter.validFile(parameters, "iters", false);		if (temp == "not found") { temp = "100";			}
			m->mothurConvert(temp, iters); 
            
            output = validParameter.validFile(parameters, "output", false);		if(output == "not found"){	output = "detail"; }
			if ((output != "simple") && (output != "detail")) { m->mothurOut(output + " is not a valid output form. Options are simple and detail. I will use detail."); m->mothurOutEndLine(); output = "detail"; }
            
			if ((method == "wang") && (search != "kmer"))  { 
				m->mothurOut("The wang method requires the kmer search. " + search + " will be disregarded, and kmer will be used." ); m->mothurOutEndLine();
				search = "kmer";
			}
            
            if ((method == "zap") && ((search != "kmer") && (search != "align")))  { 
				m->mothurOut("The zap method requires the kmer or align search. " + search + " will be disregarded, and kmer will be used." ); m->mothurOutEndLine();
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
        
        string outputMethodTag = method;
		if(method == "wang"){	classify = new Bayesian(taxonomyFileName, templateFileName, search, kmerSize, cutoff, iters, rand(), flip, writeShortcuts);	}
		else if(method == "knn"){	classify = new Knn(taxonomyFileName, templateFileName, search, kmerSize, gapOpen, gapExtend, match, misMatch, numWanted, rand());				}
        else if(method == "zap"){	
            outputMethodTag = search + "_" + outputMethodTag;
            if (search == "kmer") {   classify = new KmerTree(templateFileName, taxonomyFileName, kmerSize, cutoff); }
            else {  classify = new AlignTree(templateFileName, taxonomyFileName, cutoff);  }
        }
		else {
			m->mothurOut(search + " is not a valid method option. I will run the command using wang.");
			m->mothurOutEndLine();
			classify = new Bayesian(taxonomyFileName, templateFileName, search, kmerSize, cutoff, iters, rand(), flip, writeShortcuts);	
		}
		
		if (m->control_pressed) { delete classify; return 0; }
				
		for (int s = 0; s < fastaFileNames.size(); s++) {
		
			m->mothurOut("Classifying sequences from " + fastaFileNames[s] + " ..." ); m->mothurOutEndLine();
			
			string baseTName = m->getSimpleName(taxonomyFileName);
			
            //set rippedTaxName to 
			string RippedTaxName = "";
            bool foundDot = false;
            for (int i = baseTName.length()-1; i >= 0; i--) {
                if (foundDot && (baseTName[i] != '.')) {  RippedTaxName = baseTName[i] + RippedTaxName; }
                else if (foundDot && (baseTName[i] == '.')) {  break; }
                else if (!foundDot && (baseTName[i] == '.')) {  foundDot = true; }
            }
            //if (RippedTaxName != "") { RippedTaxName +=  "."; }   
          
			if (outputDir == "") { outputDir += m->hasPath(fastaFileNames[s]); }
            map<string, string> variables; 
            variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s]));
            variables["[tag]"] = RippedTaxName;
            variables["[tag2]"] = outputMethodTag;
			string newTaxonomyFile = getOutputFileName("taxonomy", variables);
			string newaccnosFile = getOutputFileName("accnos", variables);
			string tempTaxonomyFile = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + "taxonomy.temp";
			string taxSummary = getOutputFileName("taxsummary", variables);
			
			if ((method == "knn") && (search == "distance")) { 
				string DistName = getOutputFileName("matchdist", variables);
				classify->setDistName(DistName);  outputNames.push_back(DistName); outputTypes["matchdist"].push_back(DistName);
			}
			
			outputNames.push_back(newTaxonomyFile); outputTypes["taxonomy"].push_back(newTaxonomyFile);
			outputNames.push_back(taxSummary);	outputTypes["taxsummary"].push_back(taxSummary);
			
			int start = time(NULL);
			int numFastaSeqs = 0;
			for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
		
			vector<unsigned long long> positions; 
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
			positions = m->divideFile(fastaFileNames[s], processors);
			for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(new linePair(positions[i], positions[(i+1)]));	}
#else
			if (processors == 1) {
				lines.push_back(new linePair(0, 1000));
			}else {
				positions = m->setFilePosFasta(fastaFileNames[s], numFastaSeqs); 
                if (numFastaSeqs < processors) { processors = numFastaSeqs; }
				
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
			
			if (!m->isBlank(newaccnosFile)) { m->mothurOutEndLine(); m->mothurOut("[WARNING]: mothur reversed some your sequences for a better classification.  If you would like to take a closer look, please check " + newaccnosFile + " for the list of the sequences."); m->mothurOutEndLine(); 
                outputNames.push_back(newaccnosFile); outputTypes["accnos"].push_back(newaccnosFile);
            }else { m->mothurRemove(newaccnosFile); }

            m->mothurOutEndLine();
            m->mothurOut("It took " + toString(time(NULL) - start) + " secs to classify " + toString(numFastaSeqs) + " sequences."); m->mothurOutEndLine(); m->mothurOutEndLine();
            start = time(NULL);
            
            //read namefile
            if(namefile != "") {
                
                m->mothurOut("Reading " + namefileNames[s] + "..."); cout.flush();
                nameMap.clear(); //remove old names
                m->readNames(namefileNames[s], nameMap);
                m->mothurOut("  Done."); m->mothurOutEndLine();
            }
            
            //output taxonomy with the unclassified bins added
            ifstream inTax;
            m->openInputFile(newTaxonomyFile, inTax);
            
            ofstream outTax;
            string unclass = newTaxonomyFile + ".unclass.temp";
            m->openOutputFile(unclass, outTax);
            
            //get maxLevel from phylotree so you know how many 'unclassified's to add
            int maxLevel = classify->getMaxLevel();
            
            //read taxfile - this reading and rewriting is done to preserve the confidence scores.
            string name, taxon;
            string group = "";
            GroupMap* groupMap = NULL;
            CountTable* ct = NULL;
            PhyloSummary* taxaSum;
            
            if (hasCount) {
                ct = new CountTable();
                ct->readTable(countfileNames[s], true, false);
                taxaSum = new PhyloSummary(ct, relabund, printlevel);
            }else {
                if (groupfile != "") {  group = groupfileNames[s]; groupMap = new GroupMap(group); groupMap->readMap(); }
                taxaSum = new PhyloSummary(groupMap, relabund, printlevel);
            }
            
            while (!inTax.eof()) {
                if (m->control_pressed) { outputTypes.clear(); if (ct != NULL) { delete ct; }  if (groupMap != NULL) { delete groupMap; } delete taxaSum; for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);	} delete classify; return 0; }
                
                inTax >> name; m->gobble(inTax);
                taxon = m->getline(inTax); m->gobble(inTax);
                
                string newTax = m->addUnclassifieds(taxon, maxLevel, probs);
                
                outTax << name << '\t' << newTax << endl;
                
                if (namefile != "") {
                    itNames = nameMap.find(name);
                    
                    if (itNames == nameMap.end()) {
                        m->mothurOut(name + " is not in your name file please correct."); m->mothurOutEndLine(); exit(1);
                    }else{
                        for (int i = 0; i < itNames->second.size(); i++) {
                            taxaSum->addSeqToTree(itNames->second[i], newTax);  //add it as many times as there are identical seqs
                        }
                        itNames->second.clear();
                        nameMap.erase(itNames->first);
                    }
                }else {
                    taxaSum->addSeqToTree(name, newTax);
                }
            }
            inTax.close();
            outTax.close();
            
            m->mothurRemove(newTaxonomyFile);
            m->renameFile(unclass, newTaxonomyFile);
            
            if (m->control_pressed) {  outputTypes.clear(); if (ct != NULL) { delete ct; } if (groupMap != NULL) { delete groupMap; } for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);	} delete classify; return 0; }
            
            //print summary file
            ofstream outTaxTree;
            m->openOutputFile(taxSummary, outTaxTree);
            taxaSum->print(outTaxTree, output);
            outTaxTree.close();
            
            if (ct != NULL) { delete ct; }
            if (groupMap != NULL) { delete groupMap; } delete taxaSum;
            m->mothurRemove(tempTaxonomyFile);
            
            m->mothurOutEndLine();
            m->mothurOut("It took " + toString(time(NULL) - start) + " secs to create the summary file for " + toString(numFastaSeqs) + " sequences."); m->mothurOutEndLine(); m->mothurOutEndLine();
			
		}
        delete classify;
        
        m->mothurOutEndLine();
        m->mothurOut("Output File Names: "); m->mothurOutEndLine();
        for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
        m->mothurOutEndLine();
		
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
		
		
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/

int ClassifySeqsCommand::createProcesses(string taxFileName, string tempTaxFile, string accnos, string filename) {
	try {
		
		int num = 0;
		processIDS.clear();
        bool recalc = false;
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		int process = 1;
		
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driver(lines[process], taxFileName + m->mothurGetpid(process) + ".temp", tempTaxFile + m->mothurGetpid(process) + ".temp", accnos + m->mothurGetpid(process) + ".temp", filename);

				//pass numSeqs to parent
				ofstream out;
				string tempFile = filename + m->mothurGetpid(process) + ".num.temp";
				m->openOutputFile(tempFile, out);
				out << num << endl;
				out.close();

				exit(0);
            }else {
                m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(process) + "\n"); processors = process;
                for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                //wait to die
                for (int i=0;i<processIDS.size();i++) {
                    int temp = processIDS[i];
                    wait(&temp);
                }
                m->control_pressed = false;
                recalc = true;
                break;
            }
		}
        
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->control_pressed = false;  processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
            vector<unsigned long long> positions = m->divideFile(filename, processors);
            for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(new linePair(positions[i], positions[(i+1)]));	}
            
            num = 0;
            processIDS.resize(0);
            process = 1;
            
            while (process != processors) {
                pid_t pid = fork();
                
                if (pid > 0) {
                    processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                    process++;
                }else if (pid == 0){
                    num = driver(lines[process], taxFileName + m->mothurGetpid(process) + ".temp", tempTaxFile + m->mothurGetpid(process) + ".temp", accnos + m->mothurGetpid(process) + ".temp", filename);
                    
                    //pass numSeqs to parent
                    ofstream out;
                    string tempFile = filename + m->mothurGetpid(process) + ".num.temp";
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
            if (pDataArray[i]->count != pDataArray[i]->end) {
                m->mothurOut("[ERROR]: process " + toString(i) + " only processed " + toString(pDataArray[i]->count) + " of " + toString(pDataArray[i]->end) + " sequences assigned to it, quitting. \n"); m->control_pressed = true; 
            }
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
			if((count) % 100 == 0){	m->mothurOutJustToScreen("Processing sequence: " + toString(count) +"\n"); 		}
			
		}
		//report progress
		if((count) % 100 != 0){	m->mothurOutJustToScreen("Processing sequence: " + toString(count)+"\n"); 		}
			
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
/**************************************************************************************************/
