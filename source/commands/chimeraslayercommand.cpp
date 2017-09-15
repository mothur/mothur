/*
 *  chimeraslayercommand.cpp
 *  Mothur
 *
 *  Created by westcott on 3/31/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "chimeraslayercommand.h"
#include "deconvolutecommand.h"
#include "sequenceparser.h"
#include "counttable.h"

//**********************************************************************************************************************
vector<string> ChimeraSlayerCommand::setParameters(){	
	try {
		CommandParameter ptemplate("reference", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(ptemplate);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","chimera-accnos",false,true,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","",false,false,true); parameters.push_back(pname);
         CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
		CommandParameter pwindow("window", "Number", "", "50", "", "", "","",false,false); parameters.push_back(pwindow);
		CommandParameter pksize("ksize", "Number", "", "7", "", "", "","",false,false); parameters.push_back(pksize);
		CommandParameter pmatch("match", "Number", "", "5.0", "", "", "","",false,false); parameters.push_back(pmatch);
		CommandParameter pmismatch("mismatch", "Number", "", "-4.0", "", "", "","",false,false); parameters.push_back(pmismatch);
		CommandParameter pminsim("minsim", "Number", "", "90", "", "", "","",false,false); parameters.push_back(pminsim);
		CommandParameter pmincov("mincov", "Number", "", "70", "", "", "","",false,false); parameters.push_back(pmincov);
		CommandParameter pminsnp("minsnp", "Number", "", "10", "", "", "","",false,false); parameters.push_back(pminsnp);
		CommandParameter pminbs("minbs", "Number", "", "90", "", "", "","",false,false); parameters.push_back(pminbs);
		CommandParameter psearch("search", "Multiple", "kmer-blast", "blast", "", "", "","",false,false); parameters.push_back(psearch);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
        
		CommandParameter prealign("realign", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(prealign);
		CommandParameter ptrim("trim", "Boolean", "", "F", "", "", "","fasta",false,false); parameters.push_back(ptrim);
		CommandParameter psplit("split", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(psplit);
		CommandParameter pnumwanted("numwanted", "Number", "", "15", "", "", "","",false,false); parameters.push_back(pnumwanted);
		CommandParameter piters("iters", "Number", "", "1000", "", "", "","",false,false); parameters.push_back(piters);
		CommandParameter pdivergence("divergence", "Number", "", "1.007", "", "", "","",false,false); parameters.push_back(pdivergence);
        CommandParameter pdups("dereplicate", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pdups);
		CommandParameter pparents("parents", "Number", "", "3", "", "", "","",false,false); parameters.push_back(pparents);
		CommandParameter pincrement("increment", "Number", "", "5", "", "", "","",false,false); parameters.push_back(pincrement);
		CommandParameter pblastlocation("blastlocation", "String", "", "", "", "", "","",false,false); parameters.push_back(pblastlocation);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);

		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraSlayerCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The chimera.slayer command reads a fastafile and referencefile and outputs potentially chimeric sequences.\n";
		helpString += "This command was modeled after the chimeraSlayer written by the Broad Institute.\n";
		helpString += "The chimera.slayer command parameters are fasta, name, group, template, processors, dereplicate, trim, ksize, window, match, mismatch, divergence. minsim, mincov, minbs, minsnp, parents, search, iters, increment, numwanted, blastlocation and realign.\n";
		helpString += "The fasta parameter allows you to enter the fasta file containing your potentially chimeric sequences, and is required, unless you have a valid current fasta file. \n";
		helpString += "The name parameter allows you to provide a name file, if you are using reference=self. \n";
		helpString += "The group parameter allows you to provide a group file. The group file can be used with a namesfile and reference=self. When checking sequences, only sequences from the same group as the query sequence will be used as the reference. \n";
        helpString += "The count parameter allows you to provide a count file. The count file reference=self. If your count file contains group information, when checking sequences, only sequences from the same group as the query sequence will be used as the reference. When you use a count file with group info and dereplicate=T, mothur will create a *.pick.count_table file containing seqeunces after chimeras are removed. \n";
		helpString += "You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amazon.fasta \n";
		helpString += "The reference parameter allows you to enter a reference file containing known non-chimeric sequences, and is required. You may also set template=self, in this case the abundant sequences will be used as potential parents. \n";
		helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n";
        helpString += "If the dereplicate parameter is false, then if one group finds the seqeunce to be chimeric, then all groups find it to be chimeric, default=f.\n";
		helpString += "The trim parameter allows you to output a new fasta file containing your sequences with the chimeric ones trimmed to include only their longest piece, default=F. \n";
		helpString += "The split parameter allows you to check both pieces of non-chimeric sequence for chimeras, thus looking for trimeras and quadmeras. default=F. \n";
		helpString += "The window parameter allows you to specify the window size for searching for chimeras, default=50. \n";
		helpString += "The increment parameter allows you to specify how far you move each window while finding chimeric sequences, default=5.\n";
		helpString += "The numwanted parameter allows you to specify how many sequences you would each query sequence compared with, default=15.\n";
		helpString += "The ksize parameter allows you to input kmersize, default is 7, used if search is kmer. \n";
		helpString += "The match parameter allows you to reward matched bases in blast search, default is 5. \n";
		helpString += "The parents parameter allows you to select the number of potential parents to investigate from the numwanted best matches after rating them, default is 3. \n";
		helpString += "The mismatch parameter allows you to penalize mismatched bases in blast search, default is -4. \n";
		helpString += "The divergence parameter allows you to set a cutoff for chimera determination, default is 1.007. \n";
		helpString += "The iters parameter allows you to specify the number of bootstrap iters to do with the chimeraslayer method, default=1000.\n";
		helpString += "The minsim parameter allows you to specify a minimum similarity with the parent fragments, default=90. \n";
		helpString += "The mincov parameter allows you to specify minimum coverage by closest matches found in template. Default is 70, meaning 70%. \n";
		helpString += "The minbs parameter allows you to specify minimum bootstrap support for calling a sequence chimeric. Default is 90, meaning 90%. \n";
		helpString += "The minsnp parameter allows you to specify percent of SNPs to sample on each side of breakpoint for computing bootstrap support (default: 10) \n";
		helpString += "The search parameter allows you to specify search method for finding the closest parent. Choices are blast and kmer. Default=blast. \n";
		helpString += "The realign parameter allows you to realign the query to the potential parents. Choices are true or false, default true.  \n";
		helpString += "The blastlocation parameter allows you to specify the location of your blast executable. By default mothur will look in ./blast/bin relative to mothur's executable.  \n";
		helpString += "The chimera.slayer command should be in the following format: \n";
		helpString += "chimera.slayer(fasta=yourFastaFile, reference=yourTemplate, search=yourSearch) \n";
		helpString += "Example: chimera.slayer(fasta=AD.align, reference=core_set_aligned.imputed.fasta, search=kmer) \n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraSlayerCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "chimera") {  pattern = "[filename],slayer.chimeras"; } 
        else if (type == "accnos") {  pattern = "[filename],slayer.accnos"; } 
        else if (type == "fasta") {  pattern = "[filename],slayer.fasta"; }
        else if (type == "count") {  pattern = "[filename],slayer.pick.count_table"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraSlayerCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ChimeraSlayerCommand::ChimeraSlayerCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
		vector<string> tempOutNames;
		outputTypes["chimera"] = tempOutNames;
		outputTypes["accnos"] = tempOutNames;
		outputTypes["fasta"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "ChimeraSlayerCommand");
		exit(1);
	}
}
//***************************************************************************************************************
ChimeraSlayerCommand::ChimeraSlayerCommand(string option)  {
	try {
		abort = false; calledHelp = false;
        hasCount = false;
        hasName = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter("chimera.slayer");
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			vector<string> tempOutNames;
			outputTypes["chimera"] = tempOutNames;
			outputTypes["accnos"] = tempOutNames;
			outputTypes["fasta"] = tempOutNames;
            outputTypes["count"] = tempOutNames;
		
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
						
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", false);
			if (fastafile == "not found") { 				
				//if there is a current fasta file, use it
				string filename = m->getFastaFile(); 
				if (filename != "") { fastaFileNames.push_back(filename); m->mothurOut("Using " + filename + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { 
				m->splitAtDash(fastafile, fastaFileNames);
				
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
		
						bool ableToOpen;
						ifstream in;
						
						ableToOpen = m->openInputFile(fastaFileNames[i], in, "noerror");
					
						//if you can't open it, try default location
						if (!ableToOpen) {
							if (m->getDefaultPath() != "") { //default path is set
								string tryPath = m->getDefaultPath() + m->getSimpleName(fastaFileNames[i]);
								m->mothurOut("Unable to open " + fastaFileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								fastaFileNames[i] = tryPath;
							}
						}
						
						if (!ableToOpen) {
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
						
						if (!ableToOpen) { 
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
				if (fastaFileNames.size() == 0) { m->mothurOut("[ERROR]: no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			
			//check for required parameters
			namefile = validParameter.validFile(parameters, "name", false);
			if (namefile == "not found") { namefile = "";  	}
			else { 
				m->splitAtDash(namefile, nameFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < nameFileNames.size(); i++) {
					
					bool ignore = false;
					if (nameFileNames[i] == "current") { 
						nameFileNames[i] = m->getNameFile(); 
						if (nameFileNames[i] != "") {  m->mothurOut("Using " + nameFileNames[i] + " as input file for the name parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current namefile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							nameFileNames.erase(nameFileNames.begin()+i);
							i--;
						}
					}
					
					if (!ignore) {
						
						if (inputDir != "") {
							string path = m->hasPath(nameFileNames[i]);
							//if the user has not given a path then, add inputdir. else leave path alone.
							if (path == "") {	nameFileNames[i] = inputDir + nameFileNames[i];		}
						}
						
						bool ableToOpen;
						ifstream in;
						
						ableToOpen = m->openInputFile(nameFileNames[i], in, "noerror");
						
						//if you can't open it, try default location
						if (!ableToOpen) {
							if (m->getDefaultPath() != "") { //default path is set
								string tryPath = m->getDefaultPath() + m->getSimpleName(nameFileNames[i]);
								m->mothurOut("Unable to open " + nameFileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								nameFileNames[i] = tryPath;
							}
						}
						
						if (!ableToOpen) {
							if (m->getOutputDir() != "") { //default path is set
								string tryPath = m->getOutputDir() + m->getSimpleName(nameFileNames[i]);
								m->mothurOut("Unable to open " + nameFileNames[i] + ". Trying output directory " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								nameFileNames[i] = tryPath;
							}
						}
						
						in.close();
						
						if (!ableToOpen) { 
							m->mothurOut("Unable to open " + nameFileNames[i] + ". It will be disregarded."); m->mothurOutEndLine(); 
							//erase from file list
							nameFileNames.erase(nameFileNames.begin()+i);
							i--;
						}else {
							m->setNameFile(nameFileNames[i]);
						}
					}
				}
			}
            
            if (nameFileNames.size() != 0) { hasName = true; }
            
            //check for required parameters
            vector<string> countfileNames;
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
						if (nameFileNames[i] != "") {  m->mothurOut("Using " + countfileNames[i] + " as input file for the count parameter where you had given current."); m->mothurOutEndLine(); }
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
						
						bool ableToOpen;
						ifstream in;
						
						ableToOpen = m->openInputFile(countfileNames[i], in, "noerror");
						
						//if you can't open it, try default location
						if (!ableToOpen) {
							if (m->getDefaultPath() != "") { //default path is set
								string tryPath = m->getDefaultPath() + m->getSimpleName(countfileNames[i]);
								m->mothurOut("Unable to open " + countfileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								countfileNames[i] = tryPath;
							}
						}
						
						if (!ableToOpen) {
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
						
						if (!ableToOpen) { 
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
            
            if (countfileNames.size() != 0) { hasCount = true; }
            
			//make sure there is at least one valid file left
            if (hasName && hasCount) { m->mothurOut("[ERROR]: You must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
            
            if (!hasName && hasCount) { nameFileNames = countfileNames; }
            
			if ((hasCount || hasName) && (nameFileNames.size() != fastaFileNames.size())) { m->mothurOut("[ERROR]: The number of name or count files does not match the number of fastafiles, please correct."); m->mothurOutEndLine(); abort=true; }
			
			bool hasGroup = true;
			groupfile = validParameter.validFile(parameters, "group", false);
			if (groupfile == "not found") { groupfile = "";  hasGroup = false; }
			else { 
				m->splitAtDash(groupfile, groupFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < groupFileNames.size(); i++) {
					
					bool ignore = false;
					if (groupFileNames[i] == "current") { 
						groupFileNames[i] = m->getGroupFile(); 
						if (groupFileNames[i] != "") {  m->mothurOut("Using " + groupFileNames[i] + " as input file for the group parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current namefile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							groupFileNames.erase(groupFileNames.begin()+i);
							i--;
						}
					}
					
					if (!ignore) {
						
						if (inputDir != "") {
							string path = m->hasPath(groupFileNames[i]);
							//if the user has not given a path then, add inputdir. else leave path alone.
							if (path == "") {	groupFileNames[i] = inputDir + groupFileNames[i];		}
						}
						
						bool ableToOpen;
						ifstream in;
						
						ableToOpen = m->openInputFile(groupFileNames[i], in, "noerror");
						
						//if you can't open it, try default location
						if (!ableToOpen) {
							if (m->getDefaultPath() != "") { //default path is set
								string tryPath = m->getDefaultPath() + m->getSimpleName(groupFileNames[i]);
								m->mothurOut("Unable to open " + groupFileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								groupFileNames[i] = tryPath;
							}
						}
						
						if (!ableToOpen) {
							if (m->getOutputDir() != "") { //default path is set
								string tryPath = m->getOutputDir() + m->getSimpleName(groupFileNames[i]);
								m->mothurOut("Unable to open " + groupFileNames[i] + ". Trying output directory " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								groupFileNames[i] = tryPath;
							}
						}
						
						in.close();
						
						if (!ableToOpen) { 
							m->mothurOut("Unable to open " + groupFileNames[i] + ". It will be disregarded."); m->mothurOutEndLine(); 
							//erase from file list
							groupFileNames.erase(groupFileNames.begin()+i);
							i--;
						}else {
							m->setGroupFile(groupFileNames[i]);
						}
					}
				}
				
				//make sure there is at least one valid file left
				if (groupFileNames.size() == 0) { m->mothurOut("[ERROR]: no valid group files."); m->mothurOutEndLine(); abort = true; }
			}
			
			if (hasGroup && (groupFileNames.size() != fastaFileNames.size())) { m->mothurOut("[ERROR]: The number of groupfiles does not match the number of fastafiles, please correct."); m->mothurOutEndLine(); abort=true; }
			
            if (hasGroup && hasCount) { m->mothurOut("[ERROR]: You must enter ONLY ONE of the following: count or group."); m->mothurOutEndLine(); abort = true; }			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}
			
			string temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
			
			
			string path;
			it = parameters.find("reference");
			//user has given a template file
			if(it != parameters.end()){ 
				if (it->second == "self") { 
					templatefile = "self"; 
					if (save) {
						m->mothurOut("[WARNING]: You can't save reference=self, ignoring save."); 
						m->mothurOutEndLine();
						save = false;
					}
				}
				else {
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["reference"] = inputDir + it->second;		}
					
					templatefile = validParameter.validFile(parameters, "reference", true);
					if (templatefile == "not open") { abort = true; }
					else if (templatefile == "not found") { //check for saved reference sequences
                        m->mothurOut("[ERROR]: The reference parameter is a required, aborting.\n"); abort = true;
					}
				}
			}else if (hasName) {  templatefile = "self"; 
				if (save) {
					m->mothurOut("[WARNING]: You can't save reference=self, ignoring save."); 
					m->mothurOutEndLine();
					save = false;
				}
			}else if (hasCount) {  templatefile = "self"; 
				if (save) {
					m->mothurOut("[WARNING]: You can't save reference=self, ignoring save."); 
					m->mothurOutEndLine();
					save = false;
				}
			}
			else { 
                m->mothurOut("[ERROR]: The reference parameter is a required, aborting.\n");
					templatefile = ""; abort = true;
			}
			
			
			
			temp = validParameter.validFile(parameters, "ksize", false);			if (temp == "not found") { temp = "7"; }
			m->mothurConvert(temp, ksize);
						
			temp = validParameter.validFile(parameters, "window", false);			if (temp == "not found") { temp = "50"; }			
			m->mothurConvert(temp, window);
			
			temp = validParameter.validFile(parameters, "match", false);			if (temp == "not found") { temp = "5"; }
			m->mothurConvert(temp, match);
			
			temp = validParameter.validFile(parameters, "mismatch", false);			if (temp == "not found") { temp = "-4"; }
			m->mothurConvert(temp, mismatch);
			
			temp = validParameter.validFile(parameters, "divergence", false);		if (temp == "not found") { temp = "1.007"; }
			m->mothurConvert(temp, divR);
			
			temp = validParameter.validFile(parameters, "minsim", false);			if (temp == "not found") { temp = "90"; }
			m->mothurConvert(temp, minSimilarity);
			
			temp = validParameter.validFile(parameters, "mincov", false);			if (temp == "not found") { temp = "70"; }
			m->mothurConvert(temp, minCoverage);
			
			temp = validParameter.validFile(parameters, "minbs", false);			if (temp == "not found") { temp = "90"; }
			m->mothurConvert(temp, minBS);
			
			temp = validParameter.validFile(parameters, "minsnp", false);			if (temp == "not found") { temp = "10"; }
			m->mothurConvert(temp, minSNP);

			temp = validParameter.validFile(parameters, "parents", false);			if (temp == "not found") { temp = "3"; }
			m->mothurConvert(temp, parents); 
			
			temp = validParameter.validFile(parameters, "realign", false);			if (temp == "not found") { temp = "t"; }
			realign = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "trim", false);				if (temp == "not found") { temp = "f"; }
			trim = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "split", false);			if (temp == "not found") { temp = "f"; }
			trimera = m->isTrue(temp); 
			
			search = validParameter.validFile(parameters, "search", false);			if (search == "not found") { search = "blast"; }
			
			temp = validParameter.validFile(parameters, "iters", false);			if (temp == "not found") { temp = "1000"; }		
			m->mothurConvert(temp, iters); 
			 
			temp = validParameter.validFile(parameters, "increment", false);		if (temp == "not found") { temp = "5"; }
			m->mothurConvert(temp, increment);
			
			temp = validParameter.validFile(parameters, "numwanted", false);		if (temp == "not found") { temp = "15"; }		
			m->mothurConvert(temp, numwanted);
            
			temp = validParameter.validFile(parameters, "dereplicate", false);	
			if (temp == "not found") { temp = "false";			}
			dups = m->isTrue(temp);
			
			blastlocation = validParameter.validFile(parameters, "blastlocation", false);	
			if (blastlocation == "not found") { blastlocation = ""; }
			else {
				//add / to name if needed
				string lastChar = blastlocation.substr(blastlocation.length()-1);
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				if (lastChar != "/") { blastlocation += "/"; }
#else
				if (lastChar != "\\") { blastlocation += "\\"; }	
#endif
				blastlocation = m->getFullPathName(blastlocation);
				string formatdbCommand = "";
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				formatdbCommand = blastlocation + "formatdb";	
#else
				formatdbCommand = blastlocation + "formatdb.exe";
#endif
				
				//test to make sure formatdb exists
				ifstream in;
				formatdbCommand = m->getFullPathName(formatdbCommand);
				bool ableToOpen = m->openInputFile(formatdbCommand, in, "no error"); in.close();
				if(!ableToOpen) {	m->mothurOut("[ERROR]: " + formatdbCommand + " file does not exist. mothur requires formatdb.exe to run chimera.slayer."); m->mothurOutEndLine(); abort = true; }
				
				string blastCommand = "";
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				blastCommand = blastlocation + "megablast";	
#else
				blastCommand = blastlocation + "megablast.exe";
#endif
				//test to make sure formatdb exists
				ifstream in2;
				blastCommand = m->getFullPathName(blastCommand);
				ableToOpen = m->openInputFile(blastCommand, in2, "no error"); in2.close();
				if(!ableToOpen) {	m->mothurOut("[ERROR]: " + blastCommand + " file does not exist. mothur requires blastall.exe to run chimera.slayer."); m->mothurOutEndLine(); abort = true; }
			}

			if ((search != "blast") && (search != "kmer")) { m->mothurOut(search + " is not a valid search."); m->mothurOutEndLine(); abort = true;  }
			
			if ((hasName || hasCount) && (templatefile != "self")) { m->mothurOut("You have provided a namefile or countfile and the reference parameter is not set to self. I am not sure what reference you are trying to use, aborting."); m->mothurOutEndLine(); abort=true; }
			if (hasGroup && (templatefile != "self")) { m->mothurOut("You have provided a group file and the reference parameter is not set to self. I am not sure what reference you are trying to use, aborting."); m->mothurOutEndLine(); abort=true; }

			//until we resolve the issue 10-18-11
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
			//processors=1;
#endif
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "ChimeraSlayerCommand");
		exit(1);
	}
}
//***************************************************************************************************************

int ChimeraSlayerCommand::execute(){
	try{
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}

		for (int s = 0; s < fastaFileNames.size(); s++) {
				
			m->mothurOut("Checking sequences from " + fastaFileNames[s] + " ..." ); m->mothurOutEndLine();
		
			int start = time(NULL);	
			if (outputDir == "") { outputDir = m->hasPath(fastaFileNames[s]);  }//if user entered a file with a path then preserve it	
			map<string, string> variables; 
            variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s]));
			string outputFileName = getOutputFileName("chimera", variables);
			string accnosFileName = getOutputFileName("accnos", variables);
			string trimFastaFileName = getOutputFileName("fasta", variables);
            string newCountFile = "";
			
			//clears files
			ofstream out, out1, out2;
			m->openOutputFile(outputFileName, out); out.close(); 
			m->openOutputFile(accnosFileName, out1); out1.close();
			if (trim) { m->openOutputFile(trimFastaFileName, out2); out2.close(); }
			outputNames.push_back(outputFileName); outputTypes["chimera"].push_back(outputFileName);
			outputNames.push_back(accnosFileName); outputTypes["accnos"].push_back(accnosFileName);
			if (trim) {  outputNames.push_back(trimFastaFileName); outputTypes["fasta"].push_back(trimFastaFileName); }			
			
			//maps a filename to priority map. 
			//if no groupfile this is fastafileNames[s] -> prioirity
			//if groupfile then this is each groups seqs -> priority
			map<string, map<string, int> > fileToPriority; 
			map<string, map<string, int> >::iterator itFile;
			map<string, string> fileGroup;
			fileToPriority[fastaFileNames[s]] = priority; //default
			fileGroup[fastaFileNames[s]] = "noGroup";
            map<string, string> uniqueNames; 
			int totalChimeras = 0;
			lines.clear();
			
			if (templatefile == "self") { 
                if (hasCount) {
                    SequenceCountParser* parser = NULL;
                    setUpForSelfReference(parser, fileGroup, fileToPriority, s); 
                    if (parser != NULL) { uniqueNames = parser->getAllSeqsMap(); delete parser; }
                }else {
                    SequenceParser* parser = NULL;
                    setUpForSelfReference(parser, fileGroup, fileToPriority, s); 
                    if (parser != NULL) { uniqueNames = parser->getAllSeqsMap(); delete parser; }
                }
            }
			
			if (m->getControl_pressed()) {   for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0;	}

			if (fileToPriority.size() == 1) { //you running without a groupfile
				itFile = fileToPriority.begin();
				string thisFastaName = itFile->first;
				map<string, int> thisPriority = itFile->second;
				//break up file
				vector<unsigned long long> positions; 
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				positions = m->divideFile(thisFastaName, processors);
				for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(linePair(positions[i], positions[(i+1)]));	}
#else
				if (processors == 1) {	lines.push_back(linePair(0, 1000)); }
				else {
					positions = m->setFilePosFasta(thisFastaName, numSeqs); 
                    if (numSeqs < processors) { processors = numSeqs; }
					
					//figure out how many sequences you have to process
					int numSeqsPerProcessor = numSeqs / processors;
					for (int i = 0; i < processors; i++) {
						int startIndex =  i * numSeqsPerProcessor;
						if(i == (processors - 1)){	numSeqsPerProcessor = numSeqs - i * numSeqsPerProcessor; 	}
						lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
					}
				}
#endif
				if(processors == 1){ numSeqs = driver(lines[0], outputFileName, thisFastaName, accnosFileName, trimFastaFileName, thisPriority);  }
				else{ numSeqs = createProcesses(outputFileName, thisFastaName, accnosFileName, trimFastaFileName, thisPriority); }
				
				if (m->getControl_pressed()) {  outputTypes.clear(); if (trim) { m->mothurRemove(trimFastaFileName); } m->mothurRemove(outputFileName); m->mothurRemove(accnosFileName); for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0; }				

			}else { //you have provided a groupfile
                string countFile = "";
                if (hasCount) {
                    countFile = nameFileNames[s];
                    variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(nameFileNames[s]));
                    newCountFile = getOutputFileName("count", variables);
                }

				if (processors == 1) {
                    numSeqs = driverGroups(outputFileName, accnosFileName, trimFastaFileName, fileToPriority, fileGroup, newCountFile);
                    if (hasCount && dups) {
                        CountTable c; c.readTable(nameFileNames[s], true, false);
                        if (!m->isBlank(newCountFile)) {
                            ifstream in2;
                            m->openInputFile(newCountFile, in2);
                            
                            string name, group;
                            while (!in2.eof()) {
                                in2 >> name >> group; m->gobble(in2);
                                c.setAbund(name, group, 0);
                            }
                            in2.close();
                        }
                        m->mothurRemove(newCountFile);
                        c.printTable(newCountFile);
                    }

                }
				else {  numSeqs = createProcessesGroups(outputFileName, accnosFileName, trimFastaFileName, fileToPriority, fileGroup, newCountFile, countFile); 		} //destroys fileToPriority


                    if (!dups) {
                        totalChimeras = deconvoluteResults(uniqueNames, outputFileName, accnosFileName, trimFastaFileName);
                        m->mothurOutEndLine(); m->mothurOut(toString(totalChimeras) + " chimera found."); m->mothurOutEndLine();
                    }else {
                        if (hasCount) {
                            set<string> doNotRemove;
                            CountTable c; c.readTable(newCountFile, true, true);
                            vector<string> namesInTable = c.getNamesOfSeqs();
                            for (int i = 0; i < namesInTable.size(); i++) {
                                int temp = c.getNumSeqs(namesInTable[i]);
                                if (temp == 0) {  c.remove(namesInTable[i]);  }
                                else { doNotRemove.insert((namesInTable[i])); }
                            }
                            //remove names we want to keep from accnos file.
                            set<string> accnosNames = m->readAccnos(accnosFileName);
                            ofstream out2;
                            m->openOutputFile(accnosFileName, out2);
                            for (set<string>::iterator it = accnosNames.begin(); it != accnosNames.end(); it++) {
                                if (doNotRemove.count(*it) == 0) {  out2 << (*it) << endl; }
                            }
                            out2.close();
                            c.printTable(newCountFile);
                            outputNames.push_back(newCountFile); outputTypes["count"].push_back(newCountFile);
                        }
                    }
			}
			
            m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences.");	m->mothurOutEndLine();
		}
		
		//set accnos file as new current accnosfile
		string current = "";
		itTypes = outputTypes.find("accnos");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setAccnosFile(current); }
		}
		
		if (trim) {
			itTypes = outputTypes.find("fasta");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
			}
		}
		
        itTypes = outputTypes.find("count");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setCountTableFile(current); }
		}
        
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}	
		m->mothurOutEndLine();

		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int ChimeraSlayerCommand::deconvoluteResults(map<string, string>& uniqueNames, string outputFileName, string accnosFileName, string trimFileName){
	try {
		map<string, string>::iterator itUnique;
		int total = 0;
        
        if (trimera) { //add in more potential uniqueNames
            map<string, string> newUniqueNames = uniqueNames;
            for (map<string, string>::iterator it = uniqueNames.begin(); it != uniqueNames.end(); it++) {
                newUniqueNames[(it->first)+"_LEFT"] = (it->first)+"_LEFT";
                newUniqueNames[(it->first)+"_RIGHT"] = (it->first)+"_RIGHT";
            }
            uniqueNames = newUniqueNames;
            newUniqueNames.clear();
        }
		
		//edit accnos file
		ifstream in2; 
		m->openInputFile(accnosFileName, in2, "no error");
		
		ofstream out2;
		m->openOutputFile(accnosFileName+".temp", out2);
		
		string name; name = "";
		set<string> chimerasInFile;
		set<string>::iterator itChimeras;
		
		while (!in2.eof()) {
			if (m->getControl_pressed()) { in2.close(); out2.close(); m->mothurRemove(outputFileName); m->mothurRemove((accnosFileName+".temp")); return 0; }
			
			in2 >> name; m->gobble(in2);
			
			//find unique name
			itUnique = uniqueNames.find(name);
			
			if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing accnos results. Cannot find "+ name + "."); m->mothurOutEndLine(); m->setControl_pressed(true); }
			else {
				itChimeras = chimerasInFile.find((itUnique->second));
				
				if (itChimeras == chimerasInFile.end()) {
					out2 << itUnique->second << endl;
					chimerasInFile.insert((itUnique->second));
					total++;
				}
			}
		}
		in2.close();
		out2.close();
		
		m->mothurRemove(accnosFileName);
		rename((accnosFileName+".temp").c_str(), accnosFileName.c_str());
		
		
		//edit chimera file
		ifstream in; 
		m->openInputFile(outputFileName, in);
		
		ofstream out;
		m->openOutputFile(outputFileName+".temp", out); out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);

		string rest, parent1, parent2, line;
		set<string> namesInFile; //this is so if a sequence is found to be chimera in several samples we dont write it to the results file more than once
		set<string>::iterator itNames;
		
		//assumptions - in file each read will always look like...
		/*
		 F11Fcsw_92754	no
		 F11Fcsw_63104	F11Fcsw_33372	F11Fcsw_37007	0.89441	80.4469	0.2	1.03727	93.2961	52.2	no	0-241	243-369	
		 */
		
		//get header line
		if (!in.eof()) {
			line = m->getline(in); m->gobble(in);
			out << line << endl;
		}
		
		//for the chimera file, we want to make sure if any group finds a sequence to be chimeric then all groups do, 
		//so if this is a report that did not find it to be chimeric, but it appears in the accnos file, 
		//then ignore this report and continue until we find the report that found it to be chimeric
		
		while (!in.eof()) {
			
			if (m->getControl_pressed()) { in.close(); out.close(); m->mothurRemove((outputFileName+".temp")); return 0; }
			
			in >> name;		m->gobble(in);
			in >> parent1;	m->gobble(in);
			
			if (name == "Name") { //name = "Name" because we append the header line each time we add results from the groups
				line = m->getline(in); m->gobble(in);
			}else {
				if (parent1 == "no") {
					//find unique name
					itUnique = uniqueNames.find(name);
					
					if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find "+ name + "."); m->mothurOutEndLine(); m->setControl_pressed(true); }
					else {
						//is this sequence really not chimeric??
						itChimeras = chimerasInFile.find(itUnique->second);
						
						if (itChimeras == chimerasInFile.end()) {
							//is this sequence not already in the file
							itNames = namesInFile.find((itUnique->second));
							
							if (itNames == namesInFile.end()) { out << itUnique->second << '\t' << "no" << endl; namesInFile.insert(itUnique->second); }
						}
					}
				}else { //read the rest of the line
					double DivQLAQRB,PerIDQLAQRB,BootStrapA,DivQLBQRA,PerIDQLBQRA,BootStrapB;
					string flag, range1, range2;
					bool print = false;
					in >> parent2 >> DivQLAQRB >> PerIDQLAQRB >> BootStrapA >> DivQLBQRA >> PerIDQLBQRA >> BootStrapB >> flag >> range1 >> range2;	m->gobble(in);
					
					//find unique name
					itUnique = uniqueNames.find(name);
					
					if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find "+ name + "."); m->mothurOutEndLine(); m->setControl_pressed(true); }
					else {
						name = itUnique->second;
						//is this name already in the file
						itNames = namesInFile.find((name));
						
						if (itNames == namesInFile.end()) { //no not in file
							if (flag == "no") { //are you really a no??
								//is this sequence really not chimeric??
								itChimeras = chimerasInFile.find(name);
								
								//then you really are a no so print, otherwise skip
								if (itChimeras == chimerasInFile.end()) { print = true; }
								
							}else{ print = true; }
						}
					}
					
					if (print) {
						out << name << '\t';
						
						namesInFile.insert(name);

						//output parent1's name
						itUnique = uniqueNames.find(parent1);
						if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find parentA "+ parent1 + "."); m->mothurOutEndLine(); m->setControl_pressed(true); }
						else { out << itUnique->second << '\t'; }
						
						//output parent2's name
						itUnique = uniqueNames.find(parent2);
						if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find parentA "+ parent2 + "."); m->mothurOutEndLine(); m->setControl_pressed(true); }
						else { out << itUnique->second << '\t'; }
						
						out << DivQLAQRB << '\t' << PerIDQLAQRB << '\t' << BootStrapA << '\t' << DivQLBQRA << '\t' << PerIDQLBQRA << '\t' << BootStrapB << '\t' << flag << '\t' << range1 << '\t' << range2 << endl;
					}
				}				
			}
		}
		in.close();
		out.close();
		
		m->mothurRemove(outputFileName);
		rename((outputFileName+".temp").c_str(), outputFileName.c_str());
		
		//edit fasta file
		if (trim) {
			ifstream in3; 
			m->openInputFile(trimFileName, in3);
			
			ofstream out3;
			m->openOutputFile(trimFileName+".temp", out3);
			
			namesInFile.clear();
			
			while (!in3.eof()) {
				if (m->getControl_pressed()) { in3.close(); out3.close(); m->mothurRemove(outputFileName); m->mothurRemove(accnosFileName); m->mothurRemove((trimFileName+".temp")); return 0; }
				
				Sequence seq(in3); m->gobble(in3);
				
				if (seq.getName() != "") {
					//find unique name
					itUnique = uniqueNames.find(seq.getName());
					
					if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing accnos results. Cannot find "+ seq.getName() + "."); m->mothurOutEndLine(); m->setControl_pressed(true); }
					else {
						itNames = namesInFile.find((itUnique->second));
						
						if (itNames == namesInFile.end()) {
							seq.printSequence(out3);
						}
					}
				}
			}
			in3.close();
			out3.close();
			
			m->mothurRemove(trimFileName);
			rename((trimFileName+".temp").c_str(), trimFileName.c_str());
		}
		
		return total;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "deconvoluteResults");
		exit(1);
	}
}
//**********************************************************************************************************************
int ChimeraSlayerCommand::setUpForSelfReference(SequenceParser*& parser, map<string, string>& fileGroup, map<string, map<string, int> >& fileToPriority, int s){
	try {
		fileGroup.clear();
		fileToPriority.clear();
		
		string nameFile = "";
		if (nameFileNames.size() != 0) { //you provided a namefile and we don't need to create one
			nameFile = nameFileNames[s];
		}else {  nameFile = getNamesFile(fastaFileNames[s]); }
		
		//you provided a groupfile
		string groupFile = "";
		if (groupFileNames.size() != 0) { groupFile = groupFileNames[s]; }
		
		if (groupFile == "") { 
			if (processors != 1) { m->mothurOut("When using template=self, mothur can only use 1 processor, continuing."); m->mothurOutEndLine(); processors = 1; }
						
			//sort fastafile by abundance, returns new sorted fastafile name
			m->mothurOut("Sorting fastafile according to abundance..."); cout.flush(); 
			priority = sortFastaFile(fastaFileNames[s], nameFile);
			m->mothurOut("Done."); m->mothurOutEndLine();
			
			fileToPriority[fastaFileNames[s]] = priority;
			fileGroup[fastaFileNames[s]] = "noGroup";
		}else {
			//Parse sequences by group
            vector<string> temp;
			parser = new SequenceParser(groupFile, fastaFileNames[s], nameFile, temp);
			vector<string> groups = parser->getNamesOfGroups();
			
			for (int i = 0; i < groups.size(); i++) {
				vector<Sequence> thisGroupsSeqs = parser->getSeqs(groups[i]);
				map<string, string> thisGroupsMap = parser->getNameMap(groups[i]);
                group2NameMap[groups[i]] = thisGroupsMap;
				string newFastaFile = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + groups[i] + "-sortedTemp.fasta";
				priority = sortFastaFile(thisGroupsSeqs, thisGroupsMap, newFastaFile); 
				fileToPriority[newFastaFile] = priority;
				fileGroup[newFastaFile] = groups[i];
			}
		}
		
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "setUpForSelfReference");
		exit(1);
	}
}
//**********************************************************************************************************************
int ChimeraSlayerCommand::setUpForSelfReference(SequenceCountParser*& parser, map<string, string>& fileGroup, map<string, map<string, int> >& fileToPriority, int s){
	try {
		fileGroup.clear();
		fileToPriority.clear();
		
		string nameFile = "";
		if (nameFileNames.size() != 0) { //you provided a namefile and we don't need to create one
			nameFile = nameFileNames[s];
		}else {  m->setControl_pressed(true); return 0; }
         
		CountTable ct;
		if (!ct.testGroups(nameFile)) {  
			if (processors != 1) { m->mothurOut("When using template=self, mothur can only use 1 processor, continuing."); m->mothurOutEndLine(); processors = 1; }
            
			//sort fastafile by abundance, returns new sorted fastafile name
			m->mothurOut("Sorting fastafile according to abundance..."); cout.flush(); 
			priority = sortFastaFile(fastaFileNames[s], nameFile);
			m->mothurOut("Done."); m->mothurOutEndLine();
			
			fileToPriority[fastaFileNames[s]] = priority;
			fileGroup[fastaFileNames[s]] = "noGroup";
		}else {
			//Parse sequences by group
            vector<string> temp;
			parser = new SequenceCountParser(nameFile, fastaFileNames[s], temp);
			vector<string> groups = parser->getNamesOfGroups();
			
			for (int i = 0; i < groups.size(); i++) {
				vector<Sequence> thisGroupsSeqs = parser->getSeqs(groups[i]);
				map<string, int> thisGroupsMap = parser->getCountTable(groups[i]);
				string newFastaFile = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + groups[i] + "-sortedTemp.fasta";
				sortFastaFile(thisGroupsSeqs, thisGroupsMap, newFastaFile); 
				fileToPriority[newFastaFile] = thisGroupsMap;
				fileGroup[newFastaFile] = groups[i];
			}
		}
		
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "setUpForSelfReference");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraSlayerCommand::getNamesFile(string& inputFile){
	try {
		string nameFile = "";
		
		m->mothurOutEndLine(); m->mothurOut("No namesfile given, running unique.seqs command to generate one."); m->mothurOutEndLine(); m->mothurOutEndLine();
		
		//use unique.seqs to create new name and fastafile
		string inputString = "fasta=" + inputFile;
		m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
		m->mothurOut("Running command: unique.seqs(" + inputString + ")"); m->mothurOutEndLine(); 
		m->setMothurCalling(true);
        
		Command* uniqueCommand = new DeconvoluteCommand(inputString);
		uniqueCommand->execute();
		
		map<string, vector<string> > filenames = uniqueCommand->getOutputFiles();
		
		delete uniqueCommand;
		m->setMothurCalling(false);
		m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
		
		nameFile = filenames["name"][0];
		inputFile = filenames["fasta"][0];
		
		return nameFile;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "getNamesFile");
		exit(1);
	}
}
//**********************************************************************************************************************

int ChimeraSlayerCommand::driverGroups(string outputFName, string accnos, string fasta, map<string, map<string, int> >& fileToPriority, map<string, string>& fileGroup, string countlist){
	try {
		int totalSeqs = 0;
        ofstream outCountList;

        if (hasCount && dups) { m->openOutputFile(countlist, outCountList); }
		
		for (map<string, map<string, int> >::iterator itFile = fileToPriority.begin(); itFile != fileToPriority.end(); itFile++) {
			
			if (m->getControl_pressed()) {  return 0;  }
			
			int start = time(NULL);
			string thisFastaName = itFile->first;
			map<string, int> thisPriority = itFile->second;
			string thisoutputFileName = outputDir + m->getRootName(m->getSimpleName(thisFastaName)) + fileGroup[thisFastaName] + "slayer.chimera";
			string thisaccnosFileName = outputDir + m->getRootName(m->getSimpleName(thisFastaName)) + fileGroup[thisFastaName] + "slayer.accnos";
			string thistrimFastaFileName = outputDir + m->getRootName(m->getSimpleName(thisFastaName)) + fileGroup[thisFastaName] + "slayer.fasta";
			
			m->mothurOutEndLine(); m->mothurOut("Checking sequences from group: " + fileGroup[thisFastaName] + "."); m->mothurOutEndLine(); 
			
			lines.clear();
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
			int proc = 1;
			vector<unsigned long long> positions = m->divideFile(thisFastaName, proc);
			lines.push_back(linePair(positions[0], positions[1]));	
#else
			lines.push_back(linePair(0, 1000)); 
#endif			
			int numSeqs = driver(lines[0], thisoutputFileName, thisFastaName, thisaccnosFileName, thistrimFastaFileName, thisPriority);
			
            //if we provided a count file with group info and set dereplicate=t, then we want to create a *.pick.count_table
            //This table will zero out group counts for seqs determined to be chimeric by that group.
            if (dups) {
                if (!m->isBlank(thisaccnosFileName)) {
                    ifstream in;
                    m->openInputFile(thisaccnosFileName, in);
                    string name;
                    if (hasCount) {
                        while (!in.eof()) {
                            in >> name; m->gobble(in);
                            outCountList << name << '\t' << fileGroup[thisFastaName] << endl;
                        }
                        in.close();
                    }else {
                        map<string, map<string, string> >::iterator itGroupNameMap = group2NameMap.find(fileGroup[thisFastaName]);
                        if (itGroupNameMap != group2NameMap.end()) {
                            map<string, string> thisnamemap = itGroupNameMap->second;
                            map<string, string>::iterator itN;
                            ofstream out;
                            m->openOutputFile(thisaccnosFileName+".temp", out);
                            while (!in.eof()) {
                                in >> name; m->gobble(in);
                                itN = thisnamemap.find(name);
                                if (itN != thisnamemap.end()) {
                                    vector<string> tempNames; m->splitAtComma(itN->second, tempNames);
                                    for (int j = 0; j < tempNames.size(); j++) { out << tempNames[j] << endl; }
                                
                                }else { m->mothurOut("[ERROR]: parsing cannot find " + name + ".\n"); m->setControl_pressed(true); }
                            }
                            out.close();
                            in.close();
                            m->renameFile(thisaccnosFileName+".temp", thisaccnosFileName);
                        }else { m->mothurOut("[ERROR]: parsing cannot find " + fileGroup[thisFastaName] + ".\n"); m->setControl_pressed(true); }
                    }
                    
                }
            }

			//append files
			m->appendFiles(thisoutputFileName, outputFName); m->mothurRemove(thisoutputFileName); 
			m->appendFiles(thisaccnosFileName, accnos); m->mothurRemove(thisaccnosFileName);
			if (trim) { m->appendFiles(thistrimFastaFileName, fasta); m->mothurRemove(thistrimFastaFileName); }
			m->mothurRemove(thisFastaName);
			
			totalSeqs += numSeqs;
			
			m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences from group " + fileGroup[thisFastaName] + ".");	m->mothurOutEndLine();
		}
		
        if (hasCount && dups) { outCountList.close(); }
        
		return totalSeqs;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "driverGroups");
		exit(1);
	}
}
/**************************************************************************************************/
int ChimeraSlayerCommand::createProcessesGroups(string outputFName, string accnos, string fasta, map<string, map<string, int> >& fileToPriority, map<string, string>& fileGroup, string countlist, string countFile) {
	try {
		int process = 1;
		int num = 0;
		processIDS.clear();
        bool recalc = false;
        map<string, map<string, int> > copyFileToPriority;
        copyFileToPriority = fileToPriority;
		
		if (fileToPriority.size() < processors) { processors = fileToPriority.size(); }
        
        CountTable newCount;
        if (hasCount && dups) { newCount.readTable(countFile, true, false); }
		
		int groupsPerProcessor = fileToPriority.size() / processors;
		int remainder = fileToPriority.size() % processors;
		
		vector< map<string, map<string, int> > > breakUp;
		
		for (int i = 0; i < processors; i++) {
			map<string, map<string, int> > thisFileToPriority;
			map<string, map<string, int> >::iterator itFile;
			int count = 0;
			int enough = groupsPerProcessor;
			if (i == 0) { enough = groupsPerProcessor + remainder; }
			
			for (itFile = fileToPriority.begin(); itFile != fileToPriority.end();) {
				thisFileToPriority[itFile->first] = itFile->second;
				fileToPriority.erase(itFile++);
				count++;
				if (count == enough) { break; }
			}	
			breakUp.push_back(thisFileToPriority);
		}
				
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driverGroups(outputFName + toString(m->mothurGetpid(process)) + ".temp", accnos + m->mothurGetpid(process) + ".temp", fasta + toString(m->mothurGetpid(process)) + ".temp", breakUp[process], fileGroup, accnos + toString(m->mothurGetpid(process)) + ".byCount");
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = outputFName + toString(m->mothurGetpid(process)) + ".num.temp";
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
                m->setControl_pressed(false);
                recalc = true;
                break;
            }
		}
		
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->setControl_pressed(false);
                    processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            groupsPerProcessor = copyFileToPriority.size() / processors;
            remainder = copyFileToPriority.size() % processors;
            breakUp.clear();
            
            for (int i = 0; i < processors; i++) {
                map<string, map<string, int> > thisFileToPriority;
                map<string, map<string, int> >::iterator itFile;
                int count = 0;
                int enough = groupsPerProcessor;
                if (i == 0) { enough = groupsPerProcessor + remainder; }
                
                for (itFile = copyFileToPriority.begin(); itFile != copyFileToPriority.end();) {
                    thisFileToPriority[itFile->first] = itFile->second;
                    copyFileToPriority.erase(itFile++);
                    count++;
                    if (count == enough) { break; }
                }	
                breakUp.push_back(thisFileToPriority);
            }
            
            num = 0;
            processIDS.resize(0);
            process = 1;
            
            while (process != processors) {
                pid_t pid = fork();
                
                if (pid > 0) {
                    processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                    process++;
                }else if (pid == 0){
                    num = driverGroups(outputFName + toString(m->mothurGetpid(process)) + ".temp", accnos + m->mothurGetpid(process) + ".temp", fasta + toString(m->mothurGetpid(process)) + ".temp", breakUp[process], fileGroup, accnos + toString(m->mothurGetpid(process)) + ".byCount");
                    
                    //pass numSeqs to parent
                    ofstream out;
                    string tempFile = outputFName + toString(m->mothurGetpid(process)) + ".num.temp";
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

		num = driverGroups(outputFName, accnos, fasta, breakUp[0], fileGroup, accnos + ".byCount");

		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		for (int i = 0; i < processIDS.size(); i++) {
			ifstream in;
			string tempFile =  outputFName + toString(processIDS[i]) + ".num.temp";
			m->openInputFile(tempFile, in);
			if (!in.eof()) { int tempNum = 0;  in >> tempNum; num += tempNum; }
			in.close(); m->mothurRemove(tempFile);
		}
#else
		
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the slayerData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<slayerData*> pDataArray; 
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
		
		//Create processor worker threads.
		for(int i=1; i<processors; i++ ){
			string extension = toString(i) + ".temp";
			slayerData* tempslayer = new slayerData(group2NameMap, hasCount, dups, (accnos + toString(i) +".byCount"), (outputFName + extension), (fasta + extension), (accnos + extension), templatefile, search, blastlocation, trimera, trim, realign, m, breakUp[i], fileGroup, ksize, match, mismatch, window, minSimilarity, minCoverage, minBS, minSNP, parents, iters, increment, numwanted, divR, priority, i);
			pDataArray.push_back(tempslayer);
			processIDS.push_back(i);
			
			//MySlayerThreadFunction is in header. It must be global or static to work with the threads.
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i-1] = CreateThread(NULL, 0, MySlayerGroupThreadFunction, pDataArray[i-1], 0, &dwThreadIdArray[i-1]);   
		}
		
		num = driverGroups(outputFName, accnos, fasta, breakUp[0], fileGroup, accnos + ".byCount");
		
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
            if (pDataArray[i]->fileToPriority.size() != pDataArray[i]->end) {
                m->mothurOut("[ERROR]: process " + toString(i) + " only processed " + toString(pDataArray[i]->end) + " of " + toString(pDataArray[i]->fileToPriority.size()) + " groups assigned to it, quitting. \n"); m->setControl_pressed(true);
            }
			num += pDataArray[i]->count;
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
#endif	
        //read my own
        if (hasCount && dups) {
            if (!m->isBlank(accnos + ".byCount")) {
                ifstream in2;
                m->openInputFile(accnos + ".byCount", in2);
                
                string name, group;
                while (!in2.eof()) {
                    in2 >> name >> group; m->gobble(in2);
                    newCount.setAbund(name, group, 0);
                }
                in2.close();
            }
            m->mothurRemove(accnos + ".byCount");
        }

		
		//append output files
		for(int i=0;i<processIDS.size();i++){
			m->appendFiles((outputFName + toString(processIDS[i]) + ".temp"), outputFName);
			m->mothurRemove((outputFName + toString(processIDS[i]) + ".temp"));
			
			m->appendFiles((accnos + toString(processIDS[i]) + ".temp"), accnos);
			m->mothurRemove((accnos + toString(processIDS[i]) + ".temp"));
			
			if (trim) {
				m->appendFiles((fasta + toString(processIDS[i]) + ".temp"), fasta);
				m->mothurRemove((fasta + toString(processIDS[i]) + ".temp"));
			}
            
            if (hasCount && dups) {
                if (!m->isBlank(accnos + toString(processIDS[i]) + ".byCount")) {
                    ifstream in2;
                    m->openInputFile(accnos  + toString(processIDS[i]) + ".byCount", in2);
                    
                    string name, group;
                    while (!in2.eof()) {
                        in2 >> name >> group; m->gobble(in2);
                        newCount.setAbund(name, group, 0);
                    }
                    in2.close();
                }
                m->mothurRemove(accnos + toString(processIDS[i]) + ".byCount");
            }

		}
		
        //print new *.pick.count_table
        if (hasCount && dups) {  newCount.printTable(countlist);   }
		
		return num;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "createProcessesGroups");
		exit(1);
	}
}
//**********************************************************************************************************************

int ChimeraSlayerCommand::driver(linePair filePos, string outputFName, string filename, string accnos, string fasta, map<string, int>& priority){
	try {
		
        if (m->getDebug()) { m->mothurOut("[DEBUG]: filename = " + filename + "\n"); }
        
		MothurChimera* chimera;
		if (templatefile != "self") { //you want to run slayer with a reference template
			chimera = new ChimeraSlayer(filename, templatefile, trim, search, ksize, match, mismatch, window, divR, minSimilarity, minCoverage, minBS, minSNP, parents, iters, increment, numwanted, realign, blastlocation, m->getRandomNumber());
		}else {
			chimera = new ChimeraSlayer(filename, templatefile, trim, priority, search, ksize, match, mismatch, window, divR, minSimilarity, minCoverage, minBS, minSNP, parents, iters, increment, numwanted, realign, blastlocation, m->getRandomNumber());	
		}
		
		if (m->getControl_pressed()) { delete chimera; return 0; }
		
		if (chimera->getUnaligned()) { delete chimera; m->mothurOut("Your template sequences are different lengths, please correct."); m->mothurOutEndLine(); m->setControl_pressed(true); return 0; }
		templateSeqsLength = chimera->getLength();
		
		ofstream out;
		m->openOutputFile(outputFName, out);
		
		ofstream out2;
		m->openOutputFile(accnos, out2);
		
		ofstream out3;
		if (trim) {  m->openOutputFile(fasta, out3); }
		
		ifstream inFASTA;
		m->openInputFile(filename, inFASTA);

		inFASTA.seekg(filePos.start);
		
		if (filePos.start == 0) { chimera->printHeader(out); }

		bool done = false;
		int count = 0;
	
		while (!done) {
		
			if (m->getControl_pressed()) {	delete chimera; out.close(); out2.close(); if (trim) { out3.close(); } inFASTA.close(); return 1;	}
		
			Sequence* candidateSeq = new Sequence(inFASTA);  m->gobble(inFASTA);
			string candidateAligned = candidateSeq->getAligned();
			
			if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
				if (candidateSeq->getAligned().length() != templateSeqsLength) {  
					m->mothurOut(candidateSeq->getName() + " is not the same length as the template sequences. Skipping."); m->mothurOutEndLine();
				}else{
					//find chimeras
					chimera->getChimeras(candidateSeq);
					
					if (m->getControl_pressed()) {	delete chimera; delete candidateSeq; return 1;	}
						
					//if you are not chimeric, then check each half
					data_results wholeResults = chimera->getResults();
					
					//determine if we need to split
					bool isChimeric = false;
					
					if (wholeResults.flag == "yes") {
						string chimeraFlag = "no";
						if(  (wholeResults.results[0].bsa >= minBS && wholeResults.results[0].divr_qla_qrb >= divR)
						   ||
						   (wholeResults.results[0].bsb >= minBS && wholeResults.results[0].divr_qlb_qra >= divR) ) { chimeraFlag = "yes"; }
						
						
						if (chimeraFlag == "yes") {	
							if ((wholeResults.results[0].bsa >= minBS) || (wholeResults.results[0].bsb >= minBS)) { isChimeric = true; }
						}
					}
					
					if ((!isChimeric) && trimera) {
						
						//split sequence in half by bases
						string leftQuery, rightQuery;
						Sequence tempSeq(candidateSeq->getName(), candidateAligned);
						divideInHalf(tempSeq, leftQuery, rightQuery);
						
						//run chimeraSlayer on each piece
						Sequence* left = new Sequence(candidateSeq->getName(), leftQuery);
						Sequence* right = new Sequence(candidateSeq->getName(), rightQuery);
						
						//find chimeras
						chimera->getChimeras(left);
						data_results leftResults = chimera->getResults();
						
						chimera->getChimeras(right);
						data_results rightResults = chimera->getResults();
						
						//if either piece is chimeric then report
						Sequence trimmed = chimera->print(out, out2, leftResults, rightResults);
						if (trim) { trimmed.printSequence(out3);  }
						
						delete left; delete right;
						
					}else { //already chimeric
						//print results
						Sequence trimmed = chimera->print(out, out2);
						if (trim) { trimmed.printSequence(out3);  }
					}
					
					
				}
				count++;
			}
			
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				unsigned long long pos = inFASTA.tellg();
				if ((pos == -1) || (pos >= filePos.end)) { break; }
			#else
				if (inFASTA.eof()) { break; }
			#endif
			
			delete candidateSeq;
			//report progress
			if((count) % 100 == 0){	m->mothurOutJustToScreen("Processing sequence: " + toString(count) + "\n");		}
		}
		//report progress
		if((count) % 100 != 0){	m->mothurOutJustToScreen("Processing sequence: " + toString(count)+ "\n"); 		}
		
		int numNoParents = chimera->getNumNoParents();
		if (numNoParents == count) { m->mothurOut("[WARNING]: megablast returned 0 potential parents for all your sequences. This could be due to formatdb.exe not being setup properly, please check formatdb.log for errors."); m->mothurOutEndLine(); } 
		
		out.close();
		out2.close();
		if (trim) { out3.close(); }
		inFASTA.close();
		delete chimera;
				
		return count;
		
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "driver");
		exit(1);
	}
}
/**************************************************************************************************/

int ChimeraSlayerCommand::createProcesses(string outputFileName, string filename, string accnos, string fasta, map<string, int>& thisPriority) {
	try {
		int process = 0;
		int num = 0;
		processIDS.clear();
        bool recalc = false;
        
        if (m->getDebug()) { m->mothurOut("[DEBUG]: filename = " + filename + "\n"); }
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driver(lines[process], outputFileName + toString(m->mothurGetpid(process)) + ".temp", filename, accnos + toString(m->mothurGetpid(process)) + ".temp", fasta + toString(m->mothurGetpid(process)) + ".temp", thisPriority);
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = outputFileName + toString(m->mothurGetpid(process)) + ".num.temp";
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
                m->setControl_pressed(false);
                recalc = true;
                break;
            }
		}
        
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->setControl_pressed(false);
					  processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            lines.clear();
            vector<unsigned long long> positions = m->divideFile(filename, processors);
            for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(linePair(positions[i], positions[(i+1)]));	}
            
            num = 0;
            processIDS.resize(0);
            process = 0;
            
            while (process != processors) {
                pid_t pid = fork();
                
                if (pid > 0) {
                    processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                    process++;
                }else if (pid == 0){
                    num = driver(lines[process], outputFileName + toString(m->mothurGetpid(process)) + ".temp", filename, accnos + toString(m->mothurGetpid(process)) + ".temp", fasta + toString(m->mothurGetpid(process)) + ".temp", thisPriority);
                    
                    //pass numSeqs to parent
                    ofstream out;
                    string tempFile = outputFileName + toString(m->mothurGetpid(process)) + ".num.temp";
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

		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		for (int i = 0; i < processIDS.size(); i++) {
			ifstream in;
			string tempFile =  outputFileName + toString(processIDS[i]) + ".num.temp";
			m->openInputFile(tempFile, in);
			if (!in.eof()) { int tempNum = 0;  in >> tempNum; num += tempNum; }
			in.close(); m->mothurRemove(tempFile);
		}
#else
		
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the slayerData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<slayerData*> pDataArray; 
		DWORD   dwThreadIdArray[processors];
		HANDLE  hThreadArray[processors]; 
		
		//Create processor worker threads.
		for( int i=0; i<processors; i++ ){
			string extension = toString(i) + ".temp";
			slayerData* tempslayer = new slayerData((outputFileName + extension), (fasta + extension), (accnos + extension), filename, templatefile, search, blastlocation, trimera, trim, realign, m, lines[i].start, lines[i].end, ksize, match, mismatch, window, minSimilarity, minCoverage, minBS, minSNP, parents, iters, increment, numwanted, divR, priority, i);
			pDataArray.push_back(tempslayer);
			processIDS.push_back(i);
			
			//MySlayerThreadFunction is in header. It must be global or static to work with the threads.
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i] = CreateThread(NULL, 0, MySlayerThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);   
		}
				
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
			num += pDataArray[i]->count;
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
#endif	
		
		rename((outputFileName + toString(processIDS[0]) + ".temp").c_str(), outputFileName.c_str());
		rename((accnos + toString(processIDS[0]) + ".temp").c_str(), accnos.c_str());
		if (trim) {  rename((fasta + toString(processIDS[0]) + ".temp").c_str(), fasta.c_str()); }
		
		//append output files
		for(int i=1;i<processIDS.size();i++){
			m->appendFiles((outputFileName + toString(processIDS[i]) + ".temp"), outputFileName);
			m->mothurRemove((outputFileName + toString(processIDS[i]) + ".temp"));
			
			m->appendFiles((accnos + toString(processIDS[i]) + ".temp"), accnos);
			m->mothurRemove((accnos + toString(processIDS[i]) + ".temp"));
			
			if (trim) {
				m->appendFiles((fasta + toString(processIDS[i]) + ".temp"), fasta);
				m->mothurRemove((fasta + toString(processIDS[i]) + ".temp"));
			}
		}
		
		
		return num;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "createProcesses");
		exit(1);
	}
}

/**************************************************************************************************/

int ChimeraSlayerCommand::divideInHalf(Sequence querySeq, string& leftQuery, string& rightQuery) {
	try {
		
		string queryUnAligned = querySeq.getUnaligned();
		int numBases = int(queryUnAligned.length() * 0.5);
		
		string queryAligned = querySeq.getAligned();
		leftQuery = querySeq.getAligned();
		rightQuery = querySeq.getAligned();
		
		int baseCount = 0;
		int leftSpot = 0;
		for (int i = 0; i < queryAligned.length(); i++) {
			//if you are a base
			if (isalpha(queryAligned[i])) {		
				baseCount++; 
			}
			
			//if you have half
			if (baseCount >= numBases) {  leftSpot = i; break; } //first half
		}
		
		//blank out right side
		for (int i = leftSpot; i < leftQuery.length(); i++) { leftQuery[i] = '.'; }
		
		//blank out left side
		for (int i = 0; i < leftSpot; i++) { rightQuery[i] = '.'; }
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "divideInHalf");
		exit(1);
	}
}
/**************************************************************************************************/
map<string, int> ChimeraSlayerCommand::sortFastaFile(string fastaFile, string nameFile) {
	try {
		map<string, int> nameAbund;
		
		//read through fastafile and store info
		map<string, string> seqs;
		ifstream in;
		m->openInputFile(fastaFile, in);
		
		while (!in.eof()) {
			
			if (m->getControl_pressed()) { in.close(); return nameAbund; }
			
			Sequence seq(in); m->gobble(in);
			seqs[seq.getName()] = seq.getAligned();
		}
		
		in.close();
		
		//read namefile or countfile
		vector<seqPriorityNode> nameMapCount;
        int error = 0;
        if (hasCount) { 
            CountTable ct;
            ct.readTable(nameFile, true, false);
            
            for(map<string, string>::iterator it = seqs.begin(); it != seqs.end(); it++) {
                int num = ct.getNumSeqs(it->first);
                
                if (num == 0) { error = 1; }
                else {
                    seqPriorityNode temp(num, it->second, it->first);
                    nameMapCount.push_back(temp);
                }
            }
        }else { error = m->readNames(nameFile, nameMapCount, seqs); }
		
		if (m->getControl_pressed()) { return nameAbund; }
		
		if (error == 1) { m->setControl_pressed(true); return nameAbund; }
		if (seqs.size() != nameMapCount.size()) { m->mothurOut( "The number of sequences in your fastafile does not match the number of sequences in your namefile, aborting."); m->mothurOutEndLine(); m->setControl_pressed(true); return nameAbund; }

		sort(nameMapCount.begin(), nameMapCount.end(), compareSeqPriorityNodes);
		
		string newFasta = fastaFile + ".temp";
		ofstream out;
		m->openOutputFile(newFasta, out);
		
		//print new file in order of
		for (int i = 0; i < nameMapCount.size(); i++) {
			out << ">" << nameMapCount[i].name << endl << nameMapCount[i].seq << endl;
			nameAbund[nameMapCount[i].name] = nameMapCount[i].numIdentical;
		}
		out.close();
		
		rename(newFasta.c_str(), fastaFile.c_str());
        
		return nameAbund;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "sortFastaFile");
		exit(1);
	}
}
/**************************************************************************************************/
map<string, int> ChimeraSlayerCommand::sortFastaFile(vector<Sequence>& thisseqs, map<string, string>& nameMap, string newFile) {
	try {
		map<string, int> nameAbund;
		vector<seqPriorityNode> nameVector;
		
		//read through fastafile and store info
		map<string, string> seqs;
				
		for (int i = 0; i < thisseqs.size(); i++) {
			
			if (m->getControl_pressed()) { return nameAbund; }
			
			map<string, string>::iterator itNameMap = nameMap.find(thisseqs[i].getName());
			
			if (itNameMap == nameMap.end()){
				m->setControl_pressed(true);
				m->mothurOut("[ERROR]: " + thisseqs[i].getName() + " is in your fastafile, but is not in your namesfile, please correct."); m->mothurOutEndLine();
			}else {
				int num = m->getNumNames(itNameMap->second);
				
				seqPriorityNode temp(num, thisseqs[i].getAligned(), thisseqs[i].getName());
				nameVector.push_back(temp);
			}
		}
	
		//sort by num represented
		sort(nameVector.begin(), nameVector.end(), compareSeqPriorityNodes);
	
		if (m->getControl_pressed()) { return nameAbund; }
		
		if (thisseqs.size() != nameVector.size()) { m->mothurOut( "The number of sequences in your fastafile does not match the number of sequences in your namefile, aborting."); m->mothurOutEndLine(); m->setControl_pressed(true); return nameAbund; }
				
		ofstream out;
		m->openOutputFile(newFile, out);
		
		//print new file in order of
		for (int i = 0; i < nameVector.size(); i++) {
			out << ">" << nameVector[i].name << endl << nameVector[i].seq << endl;
			nameAbund[nameVector[i].name] = nameVector[i].numIdentical;
		}
		out.close();
		
		return nameAbund;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "sortFastaFile");
		exit(1);
	}
}
/**************************************************************************************************/
int ChimeraSlayerCommand::sortFastaFile(vector<Sequence>& thisseqs, map<string, int>& countMap, string newFile) {
	try {
		vector<seqPriorityNode> nameVector;
		
		//read through fastafile and store info
		map<string, string> seqs;
        
		for (int i = 0; i < thisseqs.size(); i++) {
			
			if (m->getControl_pressed()) { return 0; }
			
			map<string, int>::iterator itCountMap = countMap.find(thisseqs[i].getName());
			
			if (itCountMap == countMap.end()){
				m->setControl_pressed(true);
				m->mothurOut("[ERROR]: " + thisseqs[i].getName() + " is in your fastafile, but is not in your count file, please correct."); m->mothurOutEndLine();
			}else {
                seqPriorityNode temp(itCountMap->second, thisseqs[i].getAligned(), thisseqs[i].getName());
				nameVector.push_back(temp);
			}
		}
        
		//sort by num represented
		sort(nameVector.begin(), nameVector.end(), compareSeqPriorityNodes);
        
		if (m->getControl_pressed()) { return 0; }
		
		if (thisseqs.size() != nameVector.size()) { m->mothurOut( "The number of sequences in your fastafile does not match the number of sequences in your count file, aborting."); m->mothurOutEndLine(); m->setControl_pressed(true); return 0; }
        
		ofstream out;
		m->openOutputFile(newFile, out);
		
		//print new file in order of
		for (int i = 0; i < nameVector.size(); i++) {
			out << ">" << nameVector[i].name << endl << nameVector[i].seq << endl;
		}
		out.close();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "sortFastaFile");
		exit(1);
	}
}
/**************************************************************************************************/

