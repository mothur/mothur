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
#include "referencedb.h"
#include "sequenceparser.h"

//**********************************************************************************************************************
vector<string> ChimeraSlayerCommand::setParameters(){	
	try {
		CommandParameter ptemplate("reference", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(ptemplate);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pname);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pgroup);
		CommandParameter pwindow("window", "Number", "", "50", "", "", "",false,false); parameters.push_back(pwindow);
		CommandParameter pksize("ksize", "Number", "", "7", "", "", "",false,false); parameters.push_back(pksize);
		CommandParameter pmatch("match", "Number", "", "5.0", "", "", "",false,false); parameters.push_back(pmatch);
		CommandParameter pmismatch("mismatch", "Number", "", "-4.0", "", "", "",false,false); parameters.push_back(pmismatch);
		CommandParameter pminsim("minsim", "Number", "", "90", "", "", "",false,false); parameters.push_back(pminsim);
		CommandParameter pmincov("mincov", "Number", "", "70", "", "", "",false,false); parameters.push_back(pmincov);
		CommandParameter pminsnp("minsnp", "Number", "", "10", "", "", "",false,false); parameters.push_back(pminsnp);
		CommandParameter pminbs("minbs", "Number", "", "90", "", "", "",false,false); parameters.push_back(pminbs);
		CommandParameter psearch("search", "Multiple", "kmer-blast", "blast", "", "", "",false,false); parameters.push_back(psearch);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "",false,false); parameters.push_back(pprocessors);
		CommandParameter prealign("realign", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(prealign);
		CommandParameter ptrim("trim", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(ptrim);
		CommandParameter psplit("split", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(psplit);
		CommandParameter pnumwanted("numwanted", "Number", "", "15", "", "", "",false,false); parameters.push_back(pnumwanted);
		CommandParameter piters("iters", "Number", "", "1000", "", "", "",false,false); parameters.push_back(piters);
		CommandParameter pdivergence("divergence", "Number", "", "1.007", "", "", "",false,false); parameters.push_back(pdivergence);
		CommandParameter pparents("parents", "Number", "", "3", "", "", "",false,false); parameters.push_back(pparents);
		CommandParameter pincrement("increment", "Number", "", "5", "", "", "",false,false); parameters.push_back(pincrement);
		CommandParameter pblastlocation("blastlocation", "String", "", "", "", "", "",false,false); parameters.push_back(pblastlocation);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		CommandParameter psave("save", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(psave);

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
		helpString += "The chimera.slayer command parameters are fasta, name, template, processors, trim, ksize, window, match, mismatch, divergence. minsim, mincov, minbs, minsnp, parents, search, iters, increment, numwanted, blastlocation and realign.\n";
		helpString += "The fasta parameter allows you to enter the fasta file containing your potentially chimeric sequences, and is required, unless you have a valid current fasta file. \n";
		helpString += "The name parameter allows you to provide a name file, if you are using reference=self. \n";
		helpString += "The group parameter allows you to provide a group file. The group file can be used with a namesfile and reference=self. When checking sequences, only sequences from the same group as the query sequence will be used as the reference. \n";
		helpString += "You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amazon.fasta \n";
		helpString += "The reference parameter allows you to enter a reference file containing known non-chimeric sequences, and is required. You may also set template=self, in this case the abundant sequences will be used as potential parents. \n";
		helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n";
#ifdef USE_MPI
		helpString += "When using MPI, the processors parameter is set to the number of MPI processes running. \n";
#endif
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
		helpString += "The search parameter allows you to specify search method for finding the closest parent. Choices are blast, and kmer, default blast. \n";
		helpString += "The realign parameter allows you to realign the query to the potential parents. Choices are true or false, default true.  \n";
		helpString += "The blastlocation parameter allows you to specify the location of your blast executable. By default mothur will look in ./blast/bin relative to mothur's executable.  \n";
		helpString += "If the save parameter is set to true the reference sequences will be saved in memory, to clear them later you can use the clear.memory command. Default=f.";
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
ChimeraSlayerCommand::ChimeraSlayerCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
		vector<string> tempOutNames;
		outputTypes["chimera"] = tempOutNames;
		outputTypes["accnos"] = tempOutNames;
		outputTypes["fasta"] = tempOutNames;
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
		ReferenceDB* rdb = ReferenceDB::getInstance();
		
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
				if (fastaFileNames.size() == 0) { m->mothurOut("[ERROR]: no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			
			//check for required parameters
			bool hasName = true;
			namefile = validParameter.validFile(parameters, "name", false);
			if (namefile == "not found") { namefile = "";  hasName = false; }
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
						
						int ableToOpen;
						ifstream in;
						
						ableToOpen = m->openInputFile(nameFileNames[i], in, "noerror");
						
						//if you can't open it, try default location
						if (ableToOpen == 1) {
							if (m->getDefaultPath() != "") { //default path is set
								string tryPath = m->getDefaultPath() + m->getSimpleName(nameFileNames[i]);
								m->mothurOut("Unable to open " + nameFileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								nameFileNames[i] = tryPath;
							}
						}
						
						if (ableToOpen == 1) {
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
						
						if (ableToOpen == 1) { 
							m->mothurOut("Unable to open " + nameFileNames[i] + ". It will be disregarded."); m->mothurOutEndLine(); 
							//erase from file list
							nameFileNames.erase(nameFileNames.begin()+i);
							i--;
						}else {
							m->setNameFile(nameFileNames[i]);
						}
					}
				}
				
				//make sure there is at least one valid file left
				if (nameFileNames.size() == 0) { m->mothurOut("[ERROR]: no valid name files."); m->mothurOutEndLine(); abort = true; }
			}
			
			if (hasName && (nameFileNames.size() != fastaFileNames.size())) { m->mothurOut("[ERROR]: The number of namefiles does not match the number of fastafiles, please correct."); m->mothurOutEndLine(); abort=true; }
			
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
						
						int ableToOpen;
						ifstream in;
						
						ableToOpen = m->openInputFile(groupFileNames[i], in, "noerror");
						
						//if you can't open it, try default location
						if (ableToOpen == 1) {
							if (m->getDefaultPath() != "") { //default path is set
								string tryPath = m->getDefaultPath() + m->getSimpleName(groupFileNames[i]);
								m->mothurOut("Unable to open " + groupFileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								groupFileNames[i] = tryPath;
							}
						}
						
						if (ableToOpen == 1) {
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
						
						if (ableToOpen == 1) { 
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
			
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}
			
			string temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			convert(temp, processors);
			
			temp = validParameter.validFile(parameters, "save", false);			if (temp == "not found"){	temp = "f";				}
			save = m->isTrue(temp); 
			rdb->save = save; 
			if (save) { //clear out old references
				rdb->clearMemory();	
			}
			
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
						if (rdb->referenceSeqs.size() != 0) {
							templatefile = "saved";
						}else {
							m->mothurOut("[ERROR]: You don't have any saved reference sequences and the reference parameter is a required."); 
							m->mothurOutEndLine();
							abort = true; 
						}
					}else {	if (save) {	rdb->setSavedReference(templatefile);	}	}	
				}
			}else if (hasName) {  templatefile = "self"; 
				if (save) {
					m->mothurOut("[WARNING]: You can't save reference=self, ignoring save."); 
					m->mothurOutEndLine();
					save = false;
				}
			}
			else { 
				if (rdb->referenceSeqs.size() != 0) {
					templatefile = "saved";
				}else {
					m->mothurOut("[ERROR]: You don't have any saved reference sequences and the reference parameter is a required."); 
					m->mothurOutEndLine();
					templatefile = ""; abort = true; 
				} 
			}
			
			
			
			temp = validParameter.validFile(parameters, "ksize", false);			if (temp == "not found") { temp = "7"; }
			convert(temp, ksize);
						
			temp = validParameter.validFile(parameters, "window", false);			if (temp == "not found") { temp = "50"; }			
			convert(temp, window);
			
			temp = validParameter.validFile(parameters, "match", false);			if (temp == "not found") { temp = "5"; }
			convert(temp, match);
			
			temp = validParameter.validFile(parameters, "mismatch", false);			if (temp == "not found") { temp = "-4"; }
			convert(temp, mismatch);
			
			temp = validParameter.validFile(parameters, "divergence", false);		if (temp == "not found") { temp = "1.007"; }
			convert(temp, divR);
			
			temp = validParameter.validFile(parameters, "minsim", false);			if (temp == "not found") { temp = "90"; }
			convert(temp, minSimilarity);
			
			temp = validParameter.validFile(parameters, "mincov", false);			if (temp == "not found") { temp = "70"; }
			convert(temp, minCoverage);
			
			temp = validParameter.validFile(parameters, "minbs", false);			if (temp == "not found") { temp = "90"; }
			convert(temp, minBS);
			
			temp = validParameter.validFile(parameters, "minsnp", false);			if (temp == "not found") { temp = "10"; }
			convert(temp, minSNP);

			temp = validParameter.validFile(parameters, "parents", false);			if (temp == "not found") { temp = "3"; }
			convert(temp, parents); 
			
			temp = validParameter.validFile(parameters, "realign", false);			if (temp == "not found") { temp = "t"; }
			realign = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "trim", false);				if (temp == "not found") { temp = "f"; }
			trim = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "split", false);			if (temp == "not found") { temp = "f"; }
			trimera = m->isTrue(temp); 
			
			search = validParameter.validFile(parameters, "search", false);			if (search == "not found") { search = "blast"; }
			
			temp = validParameter.validFile(parameters, "iters", false);			if (temp == "not found") { temp = "1000"; }		
			convert(temp, iters); 
			 
			temp = validParameter.validFile(parameters, "increment", false);		if (temp == "not found") { temp = "5"; }
			convert(temp, increment);
			
			temp = validParameter.validFile(parameters, "numwanted", false);		if (temp == "not found") { temp = "15"; }		
			convert(temp, numwanted);
			
			blastlocation = validParameter.validFile(parameters, "blastlocation", false);	
			if (blastlocation == "not found") { blastlocation = ""; }
			else {
				//add / to name if needed
				string lastChar = blastlocation.substr(blastlocation.length()-1);
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				if (lastChar != "/") { blastlocation += "/"; }
#else
				if (lastChar != "\\") { blastlocation += "\\"; }	
#endif
				blastlocation = m->getFullPathName(blastlocation);
				string formatdbCommand = "";
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				formatdbCommand = blastlocation + "formatdb";	
#else
				formatdbCommand = blastlocation + "formatdb.exe";
#endif
				
				//test to make sure formatdb exists
				ifstream in;
				formatdbCommand = m->getFullPathName(formatdbCommand);
				int ableToOpen = m->openInputFile(formatdbCommand, in, "no error"); in.close();
				if(ableToOpen == 1) {	m->mothurOut("[ERROR]: " + formatdbCommand + " file does not exist. mothur requires formatdb.exe to run chimera.slayer."); m->mothurOutEndLine(); abort = true; }
				
				string blastCommand = "";
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				blastCommand = blastlocation + "megablast";	
#else
				blastCommand = blastlocation + "megablast.exe";
#endif
				//test to make sure formatdb exists
				ifstream in2;
				blastCommand = m->getFullPathName(blastCommand);
				ableToOpen = m->openInputFile(blastCommand, in2, "no error"); in2.close();
				if(ableToOpen == 1) {	m->mothurOut("[ERROR]: " + blastCommand + " file does not exist. mothur requires blastall.exe to run chimera.slayer."); m->mothurOutEndLine(); abort = true; }
			}

			if ((search != "blast") && (search != "kmer")) { m->mothurOut(search + " is not a valid search."); m->mothurOutEndLine(); abort = true;  }
			
			if (hasName && (templatefile != "self")) { m->mothurOut("You have provided a namefile and the reference parameter is not set to self. I am not sure what reference you are trying to use, aborting."); m->mothurOutEndLine(); abort=true; }
			if (hasGroup && (templatefile != "self")) { m->mothurOut("You have provided a group file and the reference parameter is not set to self. I am not sure what reference you are trying to use, aborting."); m->mothurOutEndLine(); abort=true; }

			//until we resolve the issue 10-18-11
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
#else
			processors=1;
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
			
			//you provided a groupfile
			string groupFile = "";
			if (groupFileNames.size() != 0) { groupFile = groupFileNames[s]; }
			
			//maps a filename to priority map. 
			//if no groupfile this is fastafileNames[s] -> prioirity
			//if groupfile then this is each groups seqs -> priority
			map<string, map<string, int> > fileToPriority; 
			map<string, map<string, int> >::iterator itFile;
			map<string, string> fileGroup;
			fileToPriority[fastaFileNames[s]] = priority; //default
			fileGroup[fastaFileNames[s]] = "noGroup";
			SequenceParser* parser = NULL;
			int totalSeqs = 0;
			int totalChimeras = 0;
			
			if ((templatefile == "self") && (groupFile == "")) { 
				fileGroup.clear();
				fileToPriority.clear();
				if (processors != 1) { m->mothurOut("When using template=self, mothur can only use 1 processor, continuing."); m->mothurOutEndLine(); processors = 1; }
				string nameFile = "";
				if (nameFileNames.size() != 0) { //you provided a namefile and we don't need to create one
					nameFile = nameFileNames[s];
				}else {  nameFile = getNamesFile(fastaFileNames[s]); }
				
				//sort fastafile by abundance, returns new sorted fastafile name
				m->mothurOut("Sorting fastafile according to abundance..."); cout.flush(); 
				priority = sortFastaFile(fastaFileNames[s], nameFile);
				m->mothurOut("Done."); m->mothurOutEndLine();
				
				fileToPriority.clear();
				fileGroup.clear();
				fileToPriority[fastaFileNames[s]] = priority;
				fileGroup[fastaFileNames[s]] = "noGroup";
				if (m->control_pressed) {  for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0;	}
			}else if ((templatefile == "self") && (groupFile != "")) {
				fileGroup.clear();
				fileToPriority.clear();
				if (processors != 1) { m->mothurOut("When using template=self, mothur can only use 1 processor, continuing."); m->mothurOutEndLine(); processors = 1; }
				string nameFile = "";
				if (nameFileNames.size() != 0) { //you provided a namefile and we don't need to create one
					nameFile = nameFileNames[s];
				}else { nameFile = getNamesFile(fastaFileNames[s]); }
				
				//Parse sequences by group
				parser = new SequenceParser(groupFile, fastaFileNames[s], nameFile);
				vector<string> groups = parser->getNamesOfGroups();
				
				for (int i = 0; i < groups.size(); i++) {
					vector<Sequence> thisGroupsSeqs = parser->getSeqs(groups[i]);
					map<string, string> thisGroupsMap = parser->getNameMap(groups[i]);
					string newFastaFile = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + groups[i] + "-sortedTemp.fasta";
					priority = sortFastaFile(thisGroupsSeqs, thisGroupsMap, newFastaFile); 
					fileToPriority[newFastaFile] = priority;
					fileGroup[newFastaFile] = groups[i];
				}
			}
			
			if (outputDir == "") { outputDir = m->hasPath(fastaFileNames[s]);  }//if user entered a file with a path then preserve it				
			string outputFileName = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + "slayer.chimera";
			string accnosFileName = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s]))  + "slayer.accnos";
			string trimFastaFileName = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s]))  + "slayer.fasta";
			
			//clears files
			ofstream out, out1, out2;
			m->openOutputFile(outputFileName, out); out.close(); 
			m->openOutputFile(accnosFileName, out1); out1.close();
			if (trim) { m->openOutputFile(trimFastaFileName, out2); out2.close(); }
			outputNames.push_back(outputFileName); outputTypes["chimera"].push_back(outputFileName);
			outputNames.push_back(accnosFileName); outputTypes["accnos"].push_back(accnosFileName);
			if (trim) {  outputNames.push_back(trimFastaFileName); outputTypes["fasta"].push_back(trimFastaFileName); }
			
			
			for (itFile = fileToPriority.begin(); itFile != fileToPriority.end(); itFile++) {
				
				string thisFastaName = itFile->first;
				map<string, int> thisPriority = itFile->second;
				
				//this is true when you have parsed by groups
				if (fileToPriority.size() > 1) { m->mothurOutEndLine(); m->mothurOut("Checking sequences from group: " + fileGroup[thisFastaName] + "."); m->mothurOutEndLine();  }
				
				string thisoutputFileName = outputDir + m->getRootName(m->getSimpleName(thisFastaName)) + fileGroup[thisFastaName] + "slayer.chimera";
				string thisaccnosFileName = outputDir + m->getRootName(m->getSimpleName(thisFastaName)) + fileGroup[thisFastaName] + "slayer.accnos";
				string thistrimFastaFileName = outputDir + m->getRootName(m->getSimpleName(thisFastaName)) + fileGroup[thisFastaName] + "slayer.fasta";
				
				//create chimera here if you are mac or linux because fork will copy for you. Create in create processes instead if you are windows.
				if (processors == 1) { templateSeqsLength = setupChimera(thisFastaName, thisPriority); }
				else {
					#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || USE_MPI
						templateSeqsLength = setupChimera(thisFastaName, thisPriority);
					#endif
				}
				
				if (m->control_pressed) {  if (parser != NULL) { delete parser; } for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0;	}
				
			#ifdef USE_MPI	
				MPIExecute(thisFastaName, thisoutputFileName, thisaccnosFileName, thistrimFastaFileName);
				if (m->control_pressed) { outputTypes.clear();  for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0;  }
			#else
					//break up file
					vector<unsigned long long> positions; 
				#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
					positions = m->divideFile(thisFastaName, processors);
					for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(new linePair(positions[i], positions[(i+1)]));	}
				#else
					if (processors == 1) {	lines.push_back(new linePair(0, 1000)); }
					else {
						positions = m->setFilePosFasta(thisFastaName, numSeqs); 
					
						//figure out how many sequences you have to process
						int numSeqsPerProcessor = numSeqs / processors;
						for (int i = 0; i < processors; i++) {
							int startIndex =  i * numSeqsPerProcessor;
							if(i == (processors - 1)){	numSeqsPerProcessor = numSeqs - i * numSeqsPerProcessor; 	}
							lines.push_back(new linePair(positions[startIndex], numSeqsPerProcessor));
						}
					}
				#endif
				
				if(processors == 1){
					numSeqs = driver(lines[0], thisoutputFileName, thisFastaName, thisaccnosFileName, thistrimFastaFileName);
					
					int numNoParents = chimera->getNumNoParents();
					if (numNoParents == numSeqs) { m->mothurOut("[WARNING]: megablast returned 0 potential parents for all your sequences. This could be due to formatdb.exe not being setup properly, please check formatdb.log for errors."); m->mothurOutEndLine(); }
					
				}else{ numSeqs = createProcesses(thisoutputFileName, thisFastaName, thisaccnosFileName, thistrimFastaFileName); }
				
				if (m->control_pressed) { if (parser != NULL) { delete parser; }  outputTypes.clear(); if (trim) { m->mothurRemove(trimFastaFileName); } m->mothurRemove(outputFileName); m->mothurRemove(accnosFileName); for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	} for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear(); delete chimera; return 0; }
				
			#endif
				
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || USE_MPI
				delete chimera;
			#endif
				
				for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
				
				//append files
				m->appendFiles(thisoutputFileName, outputFileName); m->mothurRemove(thisoutputFileName); 
				totalChimeras = m->appendFiles(thisaccnosFileName, accnosFileName); m->mothurRemove(thisaccnosFileName);
				if (trim) { m->appendFiles(thistrimFastaFileName, trimFastaFileName); m->mothurRemove(thistrimFastaFileName); }
				
				totalSeqs += numSeqs;
			}
			
			if (fileToPriority.size() > 1) { totalChimeras = deconvoluteResults(parser, outputFileName, accnosFileName, trimFastaFileName); }
			
			if (parser != NULL) { delete parser; } 
			
			m->mothurOutEndLine(); m->mothurOut(toString(totalChimeras) + " chimera found."); m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(totalSeqs) + " sequences.");	m->mothurOutEndLine();
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
int ChimeraSlayerCommand::MPIExecute(string inputFile, string outputFileName, string accnosFileName, string trimFastaFileName){
	try {
		
#ifdef USE_MPI	
		int pid, numSeqsPerProcessor; 
		int tag = 2001;
		vector<unsigned long long> MPIPos;
		
		MPI_Status status; 
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
		MPI_Comm_size(MPI_COMM_WORLD, &processors); 
		
		MPI_File inMPI;
		MPI_File outMPI;
		MPI_File outMPIAccnos;
		MPI_File outMPIFasta;
		
		int outMode=MPI_MODE_CREATE|MPI_MODE_WRONLY; 
		int inMode=MPI_MODE_RDONLY; 
		
		char outFilename[1024];
		strcpy(outFilename, outputFileName.c_str());
		
		char outAccnosFilename[1024];
		strcpy(outAccnosFilename, accnosFileName.c_str());
		
		char outFastaFilename[1024];
		strcpy(outFastaFilename, trimFastaFileName.c_str());
		
		char inFileName[1024];
		strcpy(inFileName, inputFile.c_str());
		
		MPI_File_open(MPI_COMM_WORLD, inFileName, inMode, MPI_INFO_NULL, &inMPI);  //comm, filename, mode, info, filepointer
		MPI_File_open(MPI_COMM_WORLD, outFilename, outMode, MPI_INFO_NULL, &outMPI);
		MPI_File_open(MPI_COMM_WORLD, outAccnosFilename, outMode, MPI_INFO_NULL, &outMPIAccnos);
		if (trim) { MPI_File_open(MPI_COMM_WORLD, outFastaFilename, outMode, MPI_INFO_NULL, &outMPIFasta); }
		
		if (m->control_pressed) {  MPI_File_close(&inMPI);  MPI_File_close(&outMPI); if (trim) {  MPI_File_close(&outMPIFasta);  } MPI_File_close(&outMPIAccnos);  delete chimera; return 0;  }
		
		if (pid == 0) { //you are the root process 
			m->mothurOutEndLine();
			m->mothurOut("Only reporting sequence supported by " + toString(minBS) + "% of bootstrapped results.");
			m->mothurOutEndLine();
			
			string outTemp = "Name\tLeftParent\tRightParent\tDivQLAQRB\tPerIDQLAQRB\tBootStrapA\tDivQLBQRA\tPerIDQLBQRA\tBootStrapB\tFlag\tLeftWindow\tRightWindow\n";
			
			//print header
			int length = outTemp.length();
			char* buf2 = new char[length];
			memcpy(buf2, outTemp.c_str(), length);
			
			MPI_File_write_shared(outMPI, buf2, length, MPI_CHAR, &status);
			delete buf2;
			
			MPIPos = m->setFilePosFasta(inputFile, numSeqs); //fills MPIPos, returns numSeqs
			
			if (templatefile != "self") { //if template=self we can only use 1 processor
				//send file positions to all processes
				for(int i = 1; i < processors; i++) { 
					MPI_Send(&numSeqs, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
					MPI_Send(&MPIPos[0], (numSeqs+1), MPI_LONG, i, tag, MPI_COMM_WORLD);
				}
			}
			//figure out how many sequences you have to align
			numSeqsPerProcessor = numSeqs / processors;
			int startIndex =  pid * numSeqsPerProcessor;
			if(pid == (processors - 1)){	numSeqsPerProcessor = numSeqs - pid * numSeqsPerProcessor; 	}
			
			if (templatefile == "self") { //if template=self we can only use 1 processor
				startIndex = 0;
				numSeqsPerProcessor = numSeqs;
			}
			
			//do your part
			driverMPI(startIndex, numSeqsPerProcessor, inMPI, outMPI, outMPIAccnos, outMPIFasta, MPIPos);
			
			int numNoParents = chimera->getNumNoParents();
			int temp;
			for(int i = 1; i < processors; i++) { 
				MPI_Recv(&temp, 1, MPI_INT, 1, tag, MPI_COMM_WORLD, &status);
				numNoParents += temp;
			}
			
			
			if (numSeqs == numNoParents) {  m->mothurOut("[WARNING]: megablast returned 0 potential parents for all your sequences. This could be due to formatdb.exe not being setup properly, please check formatdb.log for errors."); m->mothurOutEndLine(); }
			
			if (m->control_pressed) {  MPI_File_close(&inMPI);  MPI_File_close(&outMPI); if (trim) { MPI_File_close(&outMPIFasta); }  MPI_File_close(&outMPIAccnos);  delete chimera; return 0;  }
			
		}else{ //you are a child process
			if (templatefile != "self") { //if template=self we can only use 1 processor
				MPI_Recv(&numSeqs, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
				MPIPos.resize(numSeqs+1);
				MPI_Recv(&MPIPos[0], (numSeqs+1), MPI_LONG, 0, tag, MPI_COMM_WORLD, &status);
				
				//figure out how many sequences you have to align
				numSeqsPerProcessor = numSeqs / processors;
				int startIndex =  pid * numSeqsPerProcessor;
				if(pid == (processors - 1)){	numSeqsPerProcessor = numSeqs - pid * numSeqsPerProcessor; 	}
				
				//do your part
				driverMPI(startIndex, numSeqsPerProcessor, inMPI, outMPI, outMPIAccnos, outMPIFasta, MPIPos);
				
				int numNoParents = chimera->getNumNoParents();
				MPI_Send(&numNoParents, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
				
				if (m->control_pressed) {  MPI_File_close(&inMPI);  MPI_File_close(&outMPI); if (trim) { MPI_File_close(&outMPIFasta); }  MPI_File_close(&outMPIAccnos);  delete chimera; return 0;  }
				
			}
		}
		
		//close files 
		MPI_File_close(&inMPI);
		MPI_File_close(&outMPI);
		MPI_File_close(&outMPIAccnos); 
		if (trim) { MPI_File_close(&outMPIFasta); }
		MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
		
		
#endif		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "MPIExecute");
		exit(1);
	}
}
//**********************************************************************************************************************
int ChimeraSlayerCommand::deconvoluteResults(SequenceParser* parser, string outputFileName, string accnosFileName, string trimFileName){
	try {
		map<string, string> uniqueNames = parser->getAllSeqsMap();
		map<string, string>::iterator itUnique;
		int total = 0;
		
		//edit accnos file
		ifstream in2; 
		m->openInputFile(accnosFileName, in2, "no error");
		
		ofstream out2;
		m->openOutputFile(accnosFileName+".temp", out2);
		
		string name; name = "";
		set<string> chimerasInFile;
		set<string>::iterator itChimeras;
		
		while (!in2.eof()) {
			if (m->control_pressed) { in2.close(); out2.close(); m->mothurRemove(outputFileName); m->mothurRemove((accnosFileName+".temp")); return 0; }
			
			in2 >> name; m->gobble(in2);
			
			//find unique name
			itUnique = uniqueNames.find(name);
			
			if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing accnos results. Cannot find "+ name + "."); m->mothurOutEndLine(); m->control_pressed = true; }
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
			
			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove((outputFileName+".temp")); return 0; }
			
			in >> name;		m->gobble(in);
			in >> parent1;	m->gobble(in);
			
			if (name == "Name") { //name = "Name" because we append the header line each time we add results from the groups
				line = m->getline(in); m->gobble(in);
			}else {
				if (parent1 == "no") {
					//find unique name
					itUnique = uniqueNames.find(name);
					
					if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find "+ name + "."); m->mothurOutEndLine(); m->control_pressed = true; }
					else {
						//is this sequence really not chimeric??
						itChimeras = chimerasInFile.find(itUnique->second);
						
						if (itChimeras == chimerasInFile.end()) {
							itNames = namesInFile.find((itUnique->second));
							
							if (itNames == namesInFile.end()) {cout << itUnique->second << endl; out << itUnique->second << '\t' << "no" << endl; namesInFile.insert(itUnique->second); }
						}
					}
				}else { //read the rest of the line
					double DivQLAQRB,PerIDQLAQRB,BootStrapA,DivQLBQRA,PerIDQLBQRA,BootStrapB;
					string flag, range1, range2;
					bool print = false;
					in >> parent2 >> DivQLAQRB >> PerIDQLAQRB >> BootStrapA >> DivQLBQRA >> PerIDQLBQRA >> BootStrapB >> flag >> range1 >> range2;	m->gobble(in);
					
					//find unique name
					itUnique = uniqueNames.find(name);
					
					if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find "+ name + "."); m->mothurOutEndLine(); m->control_pressed = true; }
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
						if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find parentA "+ parent1 + "."); m->mothurOutEndLine(); m->control_pressed = true; }
						else { out << itUnique->second << '\t'; }
						
						//output parent2's name
						itUnique = uniqueNames.find(parent2);
						if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find parentA "+ parent2 + "."); m->mothurOutEndLine(); m->control_pressed = true; }
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
				if (m->control_pressed) { in3.close(); out3.close(); m->mothurRemove(outputFileName); m->mothurRemove(accnosFileName); m->mothurRemove((trimFileName+".temp")); return 0; }
				
				Sequence seq(in3); m->gobble(in3);
				
				if (seq.getName() != "") {
					//find unique name
					itUnique = uniqueNames.find(seq.getName());
					
					if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing accnos results. Cannot find "+ seq.getName() + "."); m->mothurOutEndLine(); m->control_pressed = true; }
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
int ChimeraSlayerCommand::setupChimera(string inputFile, map<string, int>& priority){
	try {
		if (templatefile != "self") { //you want to run slayer with a reference template
			chimera = new ChimeraSlayer(inputFile, templatefile, trim, search, ksize, match, mismatch, window, divR, minSimilarity, minCoverage, minBS, minSNP, parents, iters, increment, numwanted, realign, blastlocation, rand());	
		}else {
			chimera = new ChimeraSlayer(inputFile, templatefile, trim, priority, search, ksize, match, mismatch, window, divR, minSimilarity, minCoverage, minBS, minSNP, parents, iters, increment, numwanted, realign, blastlocation, rand());	
		}
		
		if (m->control_pressed) { delete chimera; return 0; }
		
		if (chimera->getUnaligned()) { delete chimera; m->mothurOut("Your template sequences are different lengths, please correct."); m->mothurOutEndLine(); m->control_pressed = true; return 0; }
		
		return (chimera->getLength());
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "setupChimera");
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
		
		Command* uniqueCommand = new DeconvoluteCommand(inputString);
		uniqueCommand->execute();
		
		map<string, vector<string> > filenames = uniqueCommand->getOutputFiles();
		
		delete uniqueCommand;
		
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

int ChimeraSlayerCommand::driver(linePair* filePos, string outputFName, string filename, string accnos, string fasta){
	try {
		ofstream out;
		m->openOutputFile(outputFName, out);
		
		ofstream out2;
		m->openOutputFile(accnos, out2);
		
		ofstream out3;
		if (trim) {  m->openOutputFile(fasta, out3); }
		
		ifstream inFASTA;
		m->openInputFile(filename, inFASTA);

		inFASTA.seekg(filePos->start);
		
		if (filePos->start == 0) { chimera->printHeader(out); }

		bool done = false;
		int count = 0;
	
		while (!done) {
		
			if (m->control_pressed) {	out.close(); out2.close(); if (trim) { out3.close(); } inFASTA.close(); return 1;	}
		
			Sequence* candidateSeq = new Sequence(inFASTA);  m->gobble(inFASTA);
			string candidateAligned = candidateSeq->getAligned();
			
			if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
				if (candidateSeq->getAligned().length() != templateSeqsLength) {  
					m->mothurOut(candidateSeq->getName() + " is not the same length as the template sequences. Skipping."); m->mothurOutEndLine();
				}else{
					//find chimeras
					chimera->getChimeras(candidateSeq);
					
					if (m->control_pressed) {	delete candidateSeq; return 1;	}
						
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
			
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				unsigned long long pos = inFASTA.tellg();
				if ((pos == -1) || (pos >= filePos->end)) { break; }
			#else
				if (inFASTA.eof()) { break; }
			#endif
			
			delete candidateSeq;
			//report progress
			if((count) % 100 == 0){	m->mothurOut("Processing sequence: " + toString(count)); m->mothurOutEndLine();		}
		}
		//report progress
		if((count) % 100 != 0){	m->mothurOut("Processing sequence: " + toString(count)); m->mothurOutEndLine();		}
		
		out.close();
		out2.close();
		if (trim) { out3.close(); }
		inFASTA.close();
				
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
#ifdef USE_MPI
int ChimeraSlayerCommand::driverMPI(int start, int num, MPI_File& inMPI, MPI_File& outMPI, MPI_File& outAccMPI, MPI_File& outFastaMPI, vector<unsigned long long>& MPIPos){
	try {				
		MPI_Status status; 
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
		
		for(int i=0;i<num;i++){
			
			if (m->control_pressed) {	return 1;	}
			
			//read next sequence
			int length = MPIPos[start+i+1] - MPIPos[start+i];

			char* buf4 = new char[length];
			MPI_File_read_at(inMPI, MPIPos[start+i], buf4, length, MPI_CHAR, &status);
	
			string tempBuf = buf4;
			if (tempBuf.length() > length) { tempBuf = tempBuf.substr(0, length);  }
			istringstream iss (tempBuf,istringstream::in);

			delete buf4;

			Sequence* candidateSeq = new Sequence(iss);  m->gobble(iss);
			string candidateAligned = candidateSeq->getAligned();
		
			if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
				
				if (candidateSeq->getAligned().length() != templateSeqsLength) {  
					m->mothurOut(candidateSeq->getName() + " is not the same length as the template sequences. Skipping."); m->mothurOutEndLine();
				}else{
		
					//find chimeras
					chimera->getChimeras(candidateSeq);
			
					if (m->control_pressed) {	delete candidateSeq; return 1;	}
					
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
						Sequence trimmed = chimera->print(outMPI, outAccMPI, leftResults, rightResults);
						if (trim) {  
							string outputString = ">" + trimmed.getName() + "\n" + trimmed.getAligned() + "\n";
							
							//write to accnos file
							int length = outputString.length();
							char* buf2 = new char[length];
							memcpy(buf2, outputString.c_str(), length);
							
							MPI_File_write_shared(outFastaMPI, buf2, length, MPI_CHAR, &status);
							delete buf2;
						}
						
						delete left; delete right;
						
					}else { 
						//print results
						Sequence trimmed = chimera->print(outMPI, outAccMPI);
						
						if (trim) {  
							string outputString = ">" + trimmed.getName() + "\n" + trimmed.getAligned() + "\n";
							
							//write to accnos file
							int length = outputString.length();
							char* buf2 = new char[length];
							memcpy(buf2, outputString.c_str(), length);
							
							MPI_File_write_shared(outFastaMPI, buf2, length, MPI_CHAR, &status);
							delete buf2;
						}
					}
					
				}
			}
			delete candidateSeq;
			
			//report progress
			if((i+1) % 100 == 0){  cout << "Processing sequence: " << (i+1) << endl;	m->mothurOutJustToLog("Processing sequence: " + toString(i+1) + "\n");		}
		}
		//report progress
		if(num % 100 != 0){		cout << "Processing sequence: " << num << endl;	m->mothurOutJustToLog("Processing sequence: " + toString(num) + "\n"); 	}
		
				
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "driverMPI");
		exit(1);
	}
}
#endif

/**************************************************************************************************/

int ChimeraSlayerCommand::createProcesses(string outputFileName, string filename, string accnos, string fasta) {
	try {
		int process = 0;
		int num = 0;
		int numNoParents = 0;
		processIDS.clear();
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driver(lines[process], outputFileName + toString(getpid()) + ".temp", filename, accnos + toString(getpid()) + ".temp", fasta + toString(getpid()) + ".temp");
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = outputFileName + toString(getpid()) + ".num.temp";
				m->openOutputFile(tempFile, out);
				out << num << '\t' << chimera->getNumNoParents() << endl;
				out.close();
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
		
		for (int i = 0; i < processIDS.size(); i++) {
			ifstream in;
			string tempFile =  outputFileName + toString(processIDS[i]) + ".num.temp";
			m->openInputFile(tempFile, in);
			if (!in.eof()) { int tempNum = 0; int tempNumParents = 0; in >> tempNum >> tempNumParents; num += tempNum; numNoParents += tempNumParents; }
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
			slayerData* tempslayer = new slayerData((outputFileName + extension), (fasta + extension), (accnos + extension), filename, templatefile, search, blastlocation, trimera, trim, realign, m, lines[i]->start, lines[i]->end, ksize, match, mismatch, window, minSimilarity, minCoverage, minBS, minSNP, parents, iters, increment, numwanted, divR, priority, i);
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
			numNoParents += pDataArray[i]->numNoParents;
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
#endif	
		if (num == numNoParents) {  m->mothurOut("[WARNING]: megablast returned 0 potential parents for all your sequences. This could be due to formatdb.exe not being setup properly, please check formatdb.log for errors."); m->mothurOutEndLine(); }
		
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
			
			if (m->control_pressed) { in.close(); return nameAbund; }
			
			Sequence seq(in); m->gobble(in);
			seqs[seq.getName()] = seq.getAligned();
		}
		
		in.close();
		
		//read namefile
		vector<seqPriorityNode> nameMapCount;
		int error = m->readNames(nameFile, nameMapCount, seqs);
		
		if (m->control_pressed) { return nameAbund; }
		
		if (error == 1) { m->control_pressed = true; return nameAbund; }
		if (seqs.size() != nameMapCount.size()) { m->mothurOut( "The number of sequences in your fastafile does not match the number of sequences in your namefile, aborting."); m->mothurOutEndLine(); m->control_pressed = true; return nameAbund; }

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
			
			if (m->control_pressed) { return nameAbund; }
			
			map<string, string>::iterator itNameMap = nameMap.find(thisseqs[i].getName());
			
			if (itNameMap == nameMap.end()){
				m->control_pressed = true;
				m->mothurOut("[ERROR]: " + thisseqs[i].getName() + " is in your fastafile, but is not in your namesfile, please correct."); m->mothurOutEndLine();
			}else {
				int num = m->getNumNames(itNameMap->second);
				
				seqPriorityNode temp(num, thisseqs[i].getAligned(), thisseqs[i].getName());
				nameVector.push_back(temp);
			}
		}
	
		//sort by num represented
		sort(nameVector.begin(), nameVector.end(), compareSeqPriorityNodes);
	
		if (m->control_pressed) { return nameAbund; }
		
		if (thisseqs.size() != nameVector.size()) { m->mothurOut( "The number of sequences in your fastafile does not match the number of sequences in your namefile, aborting."); m->mothurOutEndLine(); m->control_pressed = true; return nameAbund; }
				
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

