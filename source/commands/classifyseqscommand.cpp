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
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
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
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";		}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("reference");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["reference"] = inputDir + it->second;		}
				}
				
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
            }

			fastaFileName = validParameter.valid(parameters, "fasta");
			if (fastaFileName == "not found") { 				
				//if there is a current fasta file, use it
				string filename = current->getFastaFile(); 
				if (filename != "") { fastaFileNames.push_back(filename); m->mothurOut("Using " + filename + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else { 
				util.splitAtDash(fastaFileName, fastaFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastaFileNames.size(); i++) {
					
					bool ignore = false;
					if (fastaFileNames[i] == "current") { 
						fastaFileNames[i] = current->getFastaFile(); 
						if (fastaFileNames[i] != "") {  m->mothurOut("Using " + fastaFileNames[i] + " as input file for the fasta parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current fastafile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							fastaFileNames.erase(fastaFileNames.begin()+i);
							i--;
						}
					}
					
                    if (!ignore) {
                        if (util.checkLocations(fastaFileNames[i], current->getLocations())) { current->setFastaFile(fastaFileNames[i]); }
                        else { fastaFileNames.erase(fastaFileNames.begin()+i); i--; } //erase from file list
                    }
                }
				
				//make sure there is at least one valid file left
				if (fastaFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}

			namefile = validParameter.valid(parameters, "name");
			if (namefile == "not found") { namefile = "";  }
			else { 
				util.splitAtDash(namefile, namefileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < namefileNames.size(); i++) {
					bool ignore = false;
					if (namefileNames[i] == "current") { 
						namefileNames[i] = current->getNameFile(); 
						if (namefileNames[i] != "") {  m->mothurOut("Using " + namefileNames[i] + " as input file for the name parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current namefile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							namefileNames.erase(namefileNames.begin()+i);
							i--;
						}
					}
					
                    if (!ignore) {
                        if (util.checkLocations(namefileNames[i], current->getLocations())) { current->setNameFile(namefileNames[i]); }
                        else { namefileNames.erase(namefileNames.begin()+i); i--; } //erase from file list
                    }
                }
			}
            
            if (namefileNames.size() != 0) { hasName = true; }
            
			if (namefile != "") {
				if (namefileNames.size() != fastaFileNames.size()) { abort = true; m->mothurOut("If you provide a name file, you must have one for each fasta file."); m->mothurOutEndLine(); }
			}
			
            //check for required parameters
			countfile = validParameter.valid(parameters, "count");
			if (countfile == "not found") { 
                countfile = "";  
			}else { 
				util.splitAtDash(countfile, countfileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < countfileNames.size(); i++) {
					
					bool ignore = false;
					if (countfileNames[i] == "current") { 
						countfileNames[i] = current->getCountFile(); 
						if (countfileNames[i] != "") {  m->mothurOut("Using " + countfileNames[i] + " as input file for the count parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current count file, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							countfileNames.erase(countfileNames.begin()+i);
							i--;
						}
					}
					
                    if (!ignore) {
                        if (util.checkLocations(countfileNames[i], current->getLocations())) { current->setCountFile(countfileNames[i]); }
                        else { countfileNames.erase(countfileNames.begin()+i); i--; } //erase from file list
                    }
                }
			}
            
            if (countfileNames.size() != 0) { hasCount = true; if (countfileNames.size() != fastaFileNames.size()) {m->mothurOut("If you provide a count file, you must have one for each fasta file."); m->mothurOutEndLine(); } }
            
			//make sure there is at least one valid file left
            if (hasName && hasCount) { m->mothurOut("[ERROR]: You must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }

			groupfile = validParameter.valid(parameters, "group");
			if (groupfile == "not found") { groupfile = "";  }
			else { 
				util.splitAtDash(groupfile, groupfileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < groupfileNames.size(); i++) {
					
					bool ignore = false;
					if (groupfileNames[i] == "current") { 
						groupfileNames[i] = current->getGroupFile(); 
						if (groupfileNames[i] != "") {  m->mothurOut("Using " + groupfileNames[i] + " as input file for the group parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current group file, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							groupfileNames.erase(groupfileNames.begin()+i);
							i--;
						}
					}
					
                    if (!ignore) {
                        if (util.checkLocations(groupfileNames[i], current->getLocations())) { current->setGroupFile(groupfileNames[i]); }
                        else { groupfileNames.erase(groupfileNames.begin()+i); i--; } //erase from file list
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
			temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);

			//this has to go after save so that if the user sets save=t and provides no reference we abort
			templateFileName = validParameter.validFile(parameters, "reference");
			if (templateFileName == "not found") {
					m->mothurOut("[ERROR]: The reference parameter is a required for the classify.seqs command.\n"); abort = true;
			}else if (templateFileName == "not open") { abort = true; }
			
			
			//this has to go after save so that if the user sets save=t and provides no reference we abort
			taxonomyFileName = validParameter.validFile(parameters, "taxonomy");
			if (taxonomyFileName == "not found") {  m->mothurOut("[ERROR]: The taxonomy parameter is a required for the classify.seqs command.\n"); abort = true;
			}else if (taxonomyFileName == "not open") { abort = true; }
			
			search = validParameter.valid(parameters, "search");		if (search == "not found"){	search = "kmer";		}
			
			method = validParameter.valid(parameters, "method");		if (method == "not found"){	method = "wang";	}
            
            temp = validParameter.valid(parameters, "ksize");		if (temp == "not found"){	
                temp = "8";	
                if (method == "zap") { temp = "7"; }
            }
			util.mothurConvert(temp, kmerSize); 
			
			temp = validParameter.valid(parameters, "match");		if (temp == "not found"){	temp = "1.0";			}
			util.mothurConvert(temp, match);
            
            temp = validParameter.valid(parameters, "printlevel");		if (temp == "not found"){	temp = "-1";		}
            util.mothurConvert(temp, printlevel);
			
			temp = validParameter.valid(parameters, "mismatch");		if (temp == "not found"){	temp = "-1.0";			}
			util.mothurConvert(temp, misMatch);  
			
			temp = validParameter.valid(parameters, "gapopen");		if (temp == "not found"){	temp = "-2.0";			}
			util.mothurConvert(temp, gapOpen);  
			
			temp = validParameter.valid(parameters, "gapextend");	if (temp == "not found"){	temp = "-1.0";			}
			util.mothurConvert(temp, gapExtend); 
			
			temp = validParameter.valid(parameters, "numwanted");	if (temp == "not found"){	temp = "10";			}
			util.mothurConvert(temp, numWanted);
			
			temp = validParameter.valid(parameters, "cutoff");		if (temp == "not found"){	temp = "80";				}
			util.mothurConvert(temp, cutoff);
			
			temp = validParameter.valid(parameters, "probs");		if (temp == "not found"){	temp = "true";			}
			probs = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "relabund");		if (temp == "not found"){	temp = "false";			}
			relabund = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "shortcuts");	if (temp == "not found"){	temp = "true";			}
			writeShortcuts = util.isTrue(temp);
			
			flip = true;
			
			temp = validParameter.valid(parameters, "iters");		if (temp == "not found") { temp = "100";			}
			util.mothurConvert(temp, iters); 
            
            output = validParameter.valid(parameters, "output");		if(output == "not found"){	output = "detail"; }
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
                            if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
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
	if (!abort) {
		for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
	}
}
//**********************************************************************************************************************

int ClassifySeqsCommand::execute(){
	try {
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        string outputMethodTag = method;
		if(method == "wang"){	classify = new Bayesian(taxonomyFileName, templateFileName, search, kmerSize, cutoff, iters, util.getRandomNumber(), flip, writeShortcuts, current->getVersion());	}
		else if(method == "knn"){	classify = new Knn(taxonomyFileName, templateFileName, search, kmerSize, gapOpen, gapExtend, match, misMatch, numWanted, util.getRandomNumber(), current->getVersion());				}
        else if(method == "zap"){	
            outputMethodTag = search + "_" + outputMethodTag;
            if (search == "kmer") {   classify = new KmerTree(templateFileName, taxonomyFileName, kmerSize, cutoff); }
            else {  classify = new AlignTree(templateFileName, taxonomyFileName, cutoff);  }
        }
		else {
			m->mothurOut(search + " is not a valid method option. I will run the command using wang.");
			m->mothurOutEndLine();
			classify = new Bayesian(taxonomyFileName, templateFileName, search, kmerSize, cutoff, iters, util.getRandomNumber(), flip, writeShortcuts, current->getVersion());
		}
		
		if (m->getControl_pressed()) { delete classify; return 0; }
				
		for (int s = 0; s < fastaFileNames.size(); s++) {
		
			m->mothurOut("Classifying sequences from " + fastaFileNames[s] + " ..." ); m->mothurOutEndLine();
			
			string baseTName = util.getSimpleName(taxonomyFileName);
			
            //set rippedTaxName to 
			string RippedTaxName = "";
            bool foundDot = false;
            for (int i = baseTName.length()-1; i >= 0; i--) {
                if (foundDot && (baseTName[i] != '.')) {  RippedTaxName = baseTName[i] + RippedTaxName; }
                else if (foundDot && (baseTName[i] == '.')) {  break; }
                else if (!foundDot && (baseTName[i] == '.')) {  foundDot = true; }
            }
            //if (RippedTaxName != "") { RippedTaxName +=  "."; }   
          
			if (outputDir == "") { outputDir += util.hasPath(fastaFileNames[s]); }
            map<string, string> variables; 
            variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastaFileNames[s]));
            variables["[tag]"] = RippedTaxName;
            variables["[tag2]"] = outputMethodTag;
			string newTaxonomyFile = getOutputFileName("taxonomy", variables);
			string newaccnosFile = getOutputFileName("accnos", variables);
			string tempTaxonomyFile = outputDir + util.getRootName(util.getSimpleName(fastaFileNames[s])) + "taxonomy.temp";
			string taxSummary = getOutputFileName("taxsummary", variables);
            
			
			if ((method == "knn") && (search == "distance")) { 
				string DistName = getOutputFileName("matchdist", variables);
				classify->setDistName(DistName);  outputNames.push_back(DistName); outputTypes["matchdist"].push_back(DistName);
			}
			
			outputNames.push_back(newTaxonomyFile); outputTypes["taxonomy"].push_back(newTaxonomyFile);
			outputNames.push_back(taxSummary);	outputTypes["taxsummary"].push_back(taxSummary);
			
			long start = time(NULL);
			int numFastaSeqs = 0;
			for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
		
						
            numFastaSeqs = createProcesses(newTaxonomyFile, tempTaxonomyFile, newaccnosFile, fastaFileNames[s]);
			
			
			if (!util.isBlank(newaccnosFile)) { m->mothurOutEndLine(); m->mothurOut("[WARNING]: mothur reversed some your sequences for a better classification.  If you would like to take a closer look, please check " + newaccnosFile + " for the list of the sequences."); m->mothurOutEndLine(); 
                outputNames.push_back(newaccnosFile); outputTypes["accnos"].push_back(newaccnosFile);
            }else { util.mothurRemove(newaccnosFile); }

            m->mothurOutEndLine();
            m->mothurOut("It took " + toString(time(NULL) - start) + " secs to classify " + toString(numFastaSeqs) + " sequences."); m->mothurOutEndLine(); m->mothurOutEndLine();
            start = time(NULL);
            
            //read namefile
            if(namefile != "") {
                m->mothurOut("Reading " + namefileNames[s] + "..."); cout.flush();
                nameMap.clear(); //remove old names
                util.readNames(namefileNames[s], nameMap);
                m->mothurOut("  Done."); m->mothurOutEndLine();
            }
            
            //output taxonomy with the unclassified bins added
            ifstream inTax;
            util.openInputFile(newTaxonomyFile, inTax);
            
            ofstream outTax;
            string unclass = newTaxonomyFile + ".unclass.temp";
            util.openOutputFile(unclass, outTax);
            
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
                if (m->getControl_pressed()) { outputTypes.clear(); if (ct != NULL) { delete ct; }  if (groupMap != NULL) { delete groupMap; } delete taxaSum; for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);	} delete classify; return 0; }
                
                inTax >> name; util.gobble(inTax);
                taxon = util.getline(inTax); util.gobble(inTax);
                
                string newTax = util.addUnclassifieds(taxon, maxLevel, probs);
                
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
            
            util.mothurRemove(newTaxonomyFile);
            util.renameFile(unclass, newTaxonomyFile);
            
            if (m->getControl_pressed()) {  outputTypes.clear(); if (ct != NULL) { delete ct; } if (groupMap != NULL) { delete groupMap; } for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);	} delete classify; return 0; }
            
            //print summary file
            ofstream outTaxTree;
            util.openOutputFile(taxSummary, outTaxTree);
            taxaSum->print(outTaxTree, output);
            outTaxTree.close();
            
            if (ct != NULL) { delete ct; }
            if (groupMap != NULL) { delete groupMap; } delete taxaSum;
            util.mothurRemove(tempTaxonomyFile);
            
            m->mothurOutEndLine();
            m->mothurOut("It took " + toString(time(NULL) - start) + " secs to create the summary file for " + toString(numFastaSeqs) + " sequences."); m->mothurOutEndLine(); m->mothurOutEndLine();
			
		}
        delete classify;
        
        m->mothurOutEndLine();
        m->mothurOut("Output File Names: "); m->mothurOutEndLine();
        for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
        m->mothurOutEndLine();
		
		//set taxonomy file as new current taxonomyfile
		string currentName = "";
		itTypes = outputTypes.find("taxonomy");
		if (itTypes != outputTypes.end()) { if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setTaxonomyFile(currentName); } }
		
		currentName = "";
		itTypes = outputTypes.find("accnos");
		if (itTypes != outputTypes.end()) { if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setAccnosFile(currentName); } }

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
void driverClassifier(classifyData* params){
    try {
        ofstream outTax;
        params->util.openOutputFile(params->taxFName, outTax);
        
        ofstream outTaxSimple;
        params->util.openOutputFile(params->tempTFName, outTaxSimple);
        
        ofstream outAcc;
        params->util.openOutputFile(params->accnos, outAcc);
        
        ifstream inFASTA;
        params->util.openInputFile(params->filename, inFASTA);
        
        string taxonomy;
        
        inFASTA.seekg(params->start);
        
        bool done = false;
        
        while (!done) {
            if (params->m->getControl_pressed()) { break; }
            
            Sequence* candidateSeq = new Sequence(inFASTA); params->util.gobble(inFASTA);
            
            if (candidateSeq->getName() != "") {
                
                string simpleTax = ""; bool flipped = false;
                taxonomy = params->classify->getTaxonomy(candidateSeq, simpleTax, flipped);
                
                if (params->m->getControl_pressed()) { delete candidateSeq; break; }
                
                if (taxonomy == "unknown;") { params->m->mothurOut("[WARNING]: " + candidateSeq->getName() + " could not be classified. You can use the remove.lineage command with taxon=unknown; to remove such sequences."); params->m->mothurOutEndLine(); }
                
                //output confidence scores or not
                if (params->probs) {
                    outTax << candidateSeq->getName() << '\t' << taxonomy << endl;
                }else{
                    outTax << candidateSeq->getName() << '\t' << simpleTax << endl;
                }
                
                if (flipped) { outAcc << candidateSeq->getName() << endl; }
                
                outTaxSimple << candidateSeq->getName() << '\t' << simpleTax << endl;
                
                params->count++;
            }
            delete candidateSeq;
            
#if defined NON_WINDOWS
            unsigned long long pos = inFASTA.tellg();
            if ((pos == -1) || (pos >= params->end)) { break; }
#else
            if (params->count == params->end) { break; }
#endif
            
            //report progress
            if((params->count) % 100 == 0){	params->m->mothurOutJustToScreen("Processing sequence: " + toString(params->count) +"\n"); 		}
            
        }
        //report progress
        if((params->count) % 100 != 0){	params->m->mothurOutJustToScreen("Processing sequence: " + toString(params->count)+"\n"); 		}
        
        inFASTA.close();
        outTax.close();
        outTaxSimple.close();
        outAcc.close();
    }
    catch(exception& e) {
        params->m->errorOut(e, "ClassifySeqsCommand", "driver");
        exit(1);
    }
}
/**************************************************************************************************/

int ClassifySeqsCommand::createProcesses(string taxFileName, string tempTaxFile, string accnos, string filename) {
	try {
        
        //create array of worker threads
        vector<thread*> workerThreads;
        vector<classifyData*> data;
        
        long long num = 0;
        
        time_t start, end;
        time(&start);

        vector<unsigned long long> positions;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        positions = util.divideFile(filename, processors);
        for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(new linePair(positions[i], positions[(i+1)]));	}
#else
        if (processors == 1) {
            lines.push_back(new linePair(0, 1000));
        }else {
            positions = util.setFilePosFasta(filename, num);
            if (num < processors) { processors = num; }
            
            //figure out how many sequences you have to process
            int numSeqsPerProcessor = num / processors;
            for (int i = 0; i < processors; i++) {
                int startIndex =  i * numSeqsPerProcessor;
                if(i == (processors - 1)){	numSeqsPerProcessor = num - i * numSeqsPerProcessor; 	}
                lines.push_back(new linePair(positions[startIndex], numSeqsPerProcessor));
            }
        }
#endif
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            classifyData* dataBundle = new classifyData((accnos + toString(i+1) + ".temp"), probs, (taxFileName + toString(i+1) + ".temp"), (tempTaxFile + toString(i+1) + ".temp"), filename, lines[i+1]->start, lines[i+1]->end, flip, classify);
            data.push_back(dataBundle);
            
            workerThreads.push_back(new thread(driverClassifier, dataBundle));
        }
        
        classifyData* dataBundle = new classifyData(accnos, probs, taxFileName, tempTaxFile, filename, lines[0]->start, lines[0]->end, flip, classify);
        driverClassifier(dataBundle);
        num = dataBundle->count;
        delete dataBundle;
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            num += data[i]->count;
            
            delete data[i];
            delete workerThreads[i];
        }
        
        time(&end);
        m->mothurOut("It took " + toString(difftime(end, start)) + " secs to classify " + toString(num) + " sequences.\n\n");
        
        vector<string> nonBlankAccnosFiles;
        if (!(util.isBlank(accnos))) { nonBlankAccnosFiles.push_back(accnos); }
        else { util.mothurRemove(accnos); } //remove so other files can be renamed to it
        
        for (int i = 0; i < processors-1; i++) {
            util.appendFiles((taxFileName + toString(i+1) + ".temp"), taxFileName);
            util.appendFiles((tempTaxFile + toString(i+1) + ".temp"), tempTaxFile);
            if (!(util.isBlank(accnos + toString(i+1) + ".temp"))) {
                nonBlankAccnosFiles.push_back(accnos + toString(i+1) + ".temp");
            }else { util.mothurRemove((accnos + toString(i+1) + ".temp"));  }
            
            util.mothurRemove((util.getFullPathName(taxFileName) + toString(i+1) + ".temp"));
            util.mothurRemove((util.getFullPathName(tempTaxFile) + toString(i+1) + ".temp"));
        }
        
        //append accnos files
        if (nonBlankAccnosFiles.size() != 0) {
            rename(nonBlankAccnosFiles[0].c_str(), accnos.c_str());
            
            for (int h=1; h < nonBlankAccnosFiles.size(); h++) {
                util.appendFiles(nonBlankAccnosFiles[h], accnos);
                util.mothurRemove(nonBlankAccnosFiles[h]);
            }
        }else { //recreate the accnosfile if needed
            ofstream out;
            util.openOutputFile(accnos, out);
            out.close();
        }
        
        return num;
	}
	catch(exception& e) {
		m->errorOut(e, "ClassifySeqsCommand", "createProcesses");
		exit(1);
	}
}
/**************************************************************************************************/


