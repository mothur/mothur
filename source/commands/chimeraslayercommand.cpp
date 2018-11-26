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
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
						
			//check for required parameters
			fastafile = validParameter.valid(parameters, "fasta");
			if (fastafile == "not found") { 				
				//if there is a current fasta file, use it
				string filename = current->getFastaFile(); 
				if (filename != "") { fastaFileNames.push_back(filename); m->mothurOut("Using " + filename + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { 
				util.splitAtDash(fastafile, fastaFileNames);
				
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
				if (fastaFileNames.size() == 0) { m->mothurOut("[ERROR]: no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			
			//check for required parameters
			namefile = validParameter.valid(parameters, "name");
			if (namefile == "not found") { namefile = "";  	}
			else { 
				util.splitAtDash(namefile, nameFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < nameFileNames.size(); i++) {
					
					bool ignore = false;
					if (nameFileNames[i] == "current") { 
						nameFileNames[i] = current->getNameFile(); 
						if (nameFileNames[i] != "") {  m->mothurOut("Using " + nameFileNames[i] + " as input file for the name parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current namefile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							nameFileNames.erase(nameFileNames.begin()+i);
							i--;
						}
					}
					
                    if (!ignore) {
                        if (util.checkLocations(nameFileNames[i], current->getLocations())) { current->setNameFile(nameFileNames[i]); }
                        else { nameFileNames.erase(nameFileNames.begin()+i); i--; } //erase from file list
                    }
                }
			}
            
            if (nameFileNames.size() != 0) { hasName = true; }
            
            //check for required parameters
            vector<string> countfileNames;
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
						if (nameFileNames[i] != "") {  m->mothurOut("Using " + countfileNames[i] + " as input file for the count parameter where you had given current."); m->mothurOutEndLine(); }
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
            
            if (countfileNames.size() != 0) { hasCount = true; }
            
			//make sure there is at least one valid file left
            if (hasName && hasCount) { m->mothurOut("[ERROR]: You must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
            
            if (!hasName && hasCount) { nameFileNames = countfileNames; }
            
			if ((hasCount || hasName) && (nameFileNames.size() != fastaFileNames.size())) { m->mothurOut("[ERROR]: The number of name or count files does not match the number of fastafiles, please correct."); m->mothurOutEndLine(); abort=true; }
			
			bool hasGroup = true;
			groupfile = validParameter.valid(parameters, "group");
			if (groupfile == "not found") { groupfile = "";  hasGroup = false; }
			else { 
				util.splitAtDash(groupfile, groupFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < groupFileNames.size(); i++) {
					
					bool ignore = false;
					if (groupFileNames[i] == "current") { 
						groupFileNames[i] = current->getGroupFile();
						if (groupFileNames[i] != "") {  m->mothurOut("Using " + groupFileNames[i] + " as input file for the group parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current namefile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							groupFileNames.erase(groupFileNames.begin()+i);
							i--;
						}
					}
					
                    if (!ignore) {
                        if (util.checkLocations(groupFileNames[i], current->getLocations())) { current->setGroupFile(groupFileNames[i]); }
                        else { groupFileNames.erase(groupFileNames.begin()+i); i--; } //erase from file list
                    }
                }
				
				//make sure there is at least one valid file left
				if (groupFileNames.size() == 0) { m->mothurOut("[ERROR]: no valid group files."); m->mothurOutEndLine(); abort = true; }
			}
			
			if (hasGroup && (groupFileNames.size() != fastaFileNames.size())) { m->mothurOut("[ERROR]: The number of groupfiles does not match the number of fastafiles, please correct."); m->mothurOutEndLine(); abort=true; }
			
            if (hasGroup && hasCount) { m->mothurOut("[ERROR]: You must enter ONLY ONE of the following: count or group."); m->mothurOutEndLine(); abort = true; }			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";	}
			
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
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["reference"] = inputDir + it->second;		}
					
					templatefile = validParameter.validFile(parameters, "reference");
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
			
			string temp = validParameter.valid(parameters, "ksize");			if (temp == "not found") { temp = "7"; }
			util.mothurConvert(temp, ksize);
						
			temp = validParameter.valid(parameters, "window");			if (temp == "not found") { temp = "50"; }
			util.mothurConvert(temp, window);
			
			temp = validParameter.valid(parameters, "match");			if (temp == "not found") { temp = "5"; }
			util.mothurConvert(temp, match);
			
			temp = validParameter.valid(parameters, "mismatch");			if (temp == "not found") { temp = "-4"; }
			util.mothurConvert(temp, mismatch);
			
			temp = validParameter.valid(parameters, "divergence");		if (temp == "not found") { temp = "1.007"; }
			util.mothurConvert(temp, divR);
			
			temp = validParameter.valid(parameters, "minsim");			if (temp == "not found") { temp = "90"; }
			util.mothurConvert(temp, minSimilarity);
			
			temp = validParameter.valid(parameters, "mincov");			if (temp == "not found") { temp = "70"; }
			util.mothurConvert(temp, minCoverage);
			
			temp = validParameter.valid(parameters, "minbs");			if (temp == "not found") { temp = "90"; }
			util.mothurConvert(temp, minBS);
			
			temp = validParameter.valid(parameters, "minsnp");			if (temp == "not found") { temp = "10"; }
			util.mothurConvert(temp, minSNP);

			temp = validParameter.valid(parameters, "parents");			if (temp == "not found") { temp = "3"; }
			util.mothurConvert(temp, parents); 
			
			temp = validParameter.valid(parameters, "realign");			if (temp == "not found") { temp = "t"; }
			realign = util.isTrue(temp); 
			
			temp = validParameter.valid(parameters, "trim");				if (temp == "not found") { temp = "f"; }
			trim = util.isTrue(temp); 
			
			temp = validParameter.valid(parameters, "split");			if (temp == "not found") { temp = "f"; }
			trimera = util.isTrue(temp); 
			
			search = validParameter.valid(parameters, "search");			if (search == "not found") { search = "blast"; }
			
			temp = validParameter.valid(parameters, "iters");			if (temp == "not found") { temp = "1000"; }
			util.mothurConvert(temp, iters); 
			 
			temp = validParameter.valid(parameters, "increment");		if (temp == "not found") { temp = "5"; }
			util.mothurConvert(temp, increment);
			
			temp = validParameter.valid(parameters, "numwanted");		if (temp == "not found") { temp = "15"; }
			util.mothurConvert(temp, numwanted);
            
			temp = validParameter.valid(parameters, "dereplicate");
			if (temp == "not found") { temp = "false";			}
			dups = util.isTrue(temp);
			
			blastlocation = validParameter.valid(parameters, "blastlocation");	
			if (blastlocation == "not found") { blastlocation = ""; }
			else {
				//add / to name if needed
				string lastChar = blastlocation.substr(blastlocation.length()-1);
                if (lastChar != PATH_SEPARATOR) { blastlocation += PATH_SEPARATOR; }
				blastlocation = util.getFullPathName(blastlocation);
				
				//test to make sure formatdb exists
				ifstream in;
                string formatdbCommand = blastlocation + "formatdb" + EXECUTABLE_EXT;
				formatdbCommand = util.getFullPathName(formatdbCommand);
				bool ableToOpen = util.openInputFile(formatdbCommand, in, "no error"); in.close();
				if(!ableToOpen) {	m->mothurOut("[ERROR]: " + formatdbCommand + " file does not exist. mothur requires formatdb.exe to run chimera.slayer."); m->mothurOutEndLine(); abort = true; }

				//test to make sure formatdb exists
				ifstream in2;
                string blastCommand = blastlocation + "megablast" + EXECUTABLE_EXT;
				blastCommand = util.getFullPathName(blastCommand);
				ableToOpen = util.openInputFile(blastCommand, in2, "no error"); in2.close();
				if(!ableToOpen) {	m->mothurOut("[ERROR]: " + blastCommand + " file does not exist. mothur requires blastall.exe to run chimera.slayer."); m->mothurOutEndLine(); abort = true; }
			}

			if ((search != "blast") && (search != "kmer")) { m->mothurOut(search + " is not a valid search."); m->mothurOutEndLine(); abort = true;  }
			
			if ((hasName || hasCount) && (templatefile != "self")) { m->mothurOut("You have provided a namefile or countfile and the reference parameter is not set to self. I am not sure what reference you are trying to use, aborting."); m->mothurOutEndLine(); abort=true; }
			if (hasGroup && (templatefile != "self")) { m->mothurOut("You have provided a group file and the reference parameter is not set to self. I am not sure what reference you are trying to use, aborting."); m->mothurOutEndLine(); abort=true; }

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
		if (abort) { if (calledHelp) { return 0; }  return 2;	}

		for (int s = 0; s < fastaFileNames.size(); s++) {
				
			m->mothurOut("Checking sequences from " + fastaFileNames[s] + " ..." ); m->mothurOutEndLine();
		
			long start = time(NULL);	
			if (outputDir == "") { outputDir = util.hasPath(fastaFileNames[s]);  }//if user entered a file with a path then preserve it	
			map<string, string> variables; 
            variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastaFileNames[s]));
			string outputFileName = getOutputFileName("chimera", variables);
			string accnosFileName = getOutputFileName("accnos", variables);
			string trimFastaFileName = getOutputFileName("fasta", variables);
            string newCountFile = "";
			
			//clears files
			ofstream out, out1, out2;
			util.openOutputFile(outputFileName, out); out.close(); 
			util.openOutputFile(accnosFileName, out1); out1.close();
			if (trim) { util.openOutputFile(trimFastaFileName, out2); out2.close(); }
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
            int numSeqs = 0;
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
			
			if (m->getControl_pressed()) {   for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0;	}

			if (fileToPriority.size() == 1) { //you running without a groupfile
				itFile = fileToPriority.begin();
				string thisFastaName = itFile->first;
				map<string, int> thisPriority = itFile->second;
                numSeqs = driver(outputFileName, thisFastaName, accnosFileName, trimFastaFileName, thisPriority);
				
				if (m->getControl_pressed()) {  outputTypes.clear(); if (trim) { util.mothurRemove(trimFastaFileName); } util.mothurRemove(outputFileName); util.mothurRemove(accnosFileName); for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0; }				

			}else { //you have provided a groupfile
                string countFile = "";
                if (hasCount) {
                    countFile = nameFileNames[s];
                    variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(nameFileNames[s]));
                    newCountFile = getOutputFileName("count", variables);
                }

				
                numSeqs = driverGroups(outputFileName, accnosFileName, trimFastaFileName, fileToPriority, fileGroup, newCountFile);
                if (hasCount && dups) {
                    CountTable c; c.readTable(nameFileNames[s], true, false);
                    if (!util.isBlank(newCountFile)) {
                        ifstream in2;
                        util.openInputFile(newCountFile, in2);
                        
                        string name, group;
                        while (!in2.eof()) {
                            in2 >> name >> group; util.gobble(in2);
                            c.setAbund(name, group, 0);
                        }
                        in2.close();
                    }
                    util.mothurRemove(newCountFile);
                    c.printTable(newCountFile);
                }

                if (!dups) {
                    totalChimeras = deconvoluteResults(uniqueNames, outputFileName, accnosFileName, trimFastaFileName);
                    m->mothurOutEndLine(); m->mothurOut(toString(totalChimeras) + " chimera found."); m->mothurOutEndLine();
                }else {
                    if (hasCount) {
                        set<string> doNotRemove;
                        CountTable c; c.readTable(newCountFile, true, true);
                        c.eliminateZeroSeqs();
                        vector<string> namesInTable = c.getNamesOfSeqs();
                        for (int i = 0; i < namesInTable.size(); i++) { doNotRemove.insert(namesInTable[i]); }
                        //remove names we want to keep from accnos file.
                        set<string> accnosNames = util.readAccnos(accnosFileName);
                        ofstream out2;
                        util.openOutputFile(accnosFileName, out2);
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
		string currentName = "";
		itTypes = outputTypes.find("accnos");
		if (itTypes != outputTypes.end()) { if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setAccnosFile(currentName); } }
		
		if (trim) {
			itTypes = outputTypes.find("fasta");
			if (itTypes != outputTypes.end()) { if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); } }
		}
		
        itTypes = outputTypes.find("count");
		if (itTypes != outputTypes.end()) { if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCountFile(currentName); } }
        
		m->mothurOut("\nOutput File Names: \n"); 
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
		util.openInputFile(accnosFileName, in2, "no error");
		
		ofstream out2;
		util.openOutputFile(accnosFileName+".temp", out2);
		
		string name; name = "";
		set<string> chimerasInFile;
		set<string>::iterator itChimeras;
		
		while (!in2.eof()) {
			if (m->getControl_pressed()) { in2.close(); out2.close(); util.mothurRemove(outputFileName); util.mothurRemove((accnosFileName+".temp")); return 0; }
			
			in2 >> name; util.gobble(in2);
			
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
		
		util.mothurRemove(accnosFileName);
		rename((accnosFileName+".temp").c_str(), accnosFileName.c_str());
		
		
		//edit chimera file
		ifstream in; 
		util.openInputFile(outputFileName, in);
		
		ofstream out;
		util.openOutputFile(outputFileName+".temp", out); out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);

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
			line = util.getline(in); util.gobble(in);
			out << line << endl;
		}
		
		//for the chimera file, we want to make sure if any group finds a sequence to be chimeric then all groups do, 
		//so if this is a report that did not find it to be chimeric, but it appears in the accnos file, 
		//then ignore this report and continue until we find the report that found it to be chimeric
		
		while (!in.eof()) {
			
			if (m->getControl_pressed()) { in.close(); out.close(); util.mothurRemove((outputFileName+".temp")); return 0; }
			
			in >> name;		util.gobble(in);
			in >> parent1;	util.gobble(in);
			
			if (name == "Name") { //name = "Name" because we append the header line each time we add results from the groups
				line = util.getline(in); util.gobble(in);
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
					in >> parent2 >> DivQLAQRB >> PerIDQLAQRB >> BootStrapA >> DivQLBQRA >> PerIDQLBQRA >> BootStrapB >> flag >> range1 >> range2;	util.gobble(in);
					
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
		
		util.mothurRemove(outputFileName);
		rename((outputFileName+".temp").c_str(), outputFileName.c_str());
		
		//edit fasta file
		if (trim) {
			ifstream in3; 
			util.openInputFile(trimFileName, in3);
			
			ofstream out3;
			util.openOutputFile(trimFileName+".temp", out3);
			
			namesInFile.clear();
			
			while (!in3.eof()) {
				if (m->getControl_pressed()) { in3.close(); out3.close(); util.mothurRemove(outputFileName); util.mothurRemove(accnosFileName); util.mothurRemove((trimFileName+".temp")); return 0; }
				
				Sequence seq(in3); util.gobble(in3);
				
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
			
			util.mothurRemove(trimFileName);
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
				string newFastaFile = outputDir + util.getRootName(util.getSimpleName(fastaFileNames[s])) + groups[i] + "-sortedTemp.fasta";
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
				string newFastaFile = outputDir + util.getRootName(util.getSimpleName(fastaFileNames[s])) + groups[i] + "-sortedTemp.fasta";
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
		current->setMothurCalling(true);
        
		Command* uniqueCommand = new DeconvoluteCommand(inputString);
		uniqueCommand->execute();
		
		map<string, vector<string> > filenames = uniqueCommand->getOutputFiles();
		
		delete uniqueCommand;
		current->setMothurCalling(false);
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

        if (hasCount && dups) { util.openOutputFile(countlist, outCountList); }
		
		for (map<string, map<string, int> >::iterator itFile = fileToPriority.begin(); itFile != fileToPriority.end(); itFile++) {
			
			if (m->getControl_pressed()) {  return 0;  }
			
			long start = time(NULL);
			string thisFastaName = itFile->first;
			map<string, int> thisPriority = itFile->second;
			string thisoutputFileName = outputDir + util.getRootName(util.getSimpleName(thisFastaName)) + fileGroup[thisFastaName] + "slayer.chimera";
			string thisaccnosFileName = outputDir + util.getRootName(util.getSimpleName(thisFastaName)) + fileGroup[thisFastaName] + "slayer.accnos";
			string thistrimFastaFileName = outputDir + util.getRootName(util.getSimpleName(thisFastaName)) + fileGroup[thisFastaName] + "slayer.fasta";
			
			m->mothurOutEndLine(); m->mothurOut("Checking sequences from group: " + fileGroup[thisFastaName] + "."); m->mothurOutEndLine(); 
			
			lines.clear();
			int numSeqs = driver(thisoutputFileName, thisFastaName, thisaccnosFileName, thistrimFastaFileName, thisPriority);
			
            //if we provided a count file with group info and set dereplicate=t, then we want to create a *.pick.count_table
            //This table will zero out group counts for seqs determined to be chimeric by that group.
            if (dups) {
                if (!util.isBlank(thisaccnosFileName)) {
                    ifstream in;
                    util.openInputFile(thisaccnosFileName, in);
                    string name;
                    if (hasCount) {
                        while (!in.eof()) {
                            in >> name; util.gobble(in);
                            outCountList << name << '\t' << fileGroup[thisFastaName] << endl;
                        }
                        in.close();
                    }else {
                        map<string, map<string, string> >::iterator itGroupNameMap = group2NameMap.find(fileGroup[thisFastaName]);
                        if (itGroupNameMap != group2NameMap.end()) {
                            map<string, string> thisnamemap = itGroupNameMap->second;
                            map<string, string>::iterator itN;
                            ofstream out;
                            util.openOutputFile(thisaccnosFileName+".temp", out);
                            while (!in.eof()) {
                                in >> name; util.gobble(in);
                                itN = thisnamemap.find(name);
                                if (itN != thisnamemap.end()) {
                                    vector<string> tempNames; util.splitAtComma(itN->second, tempNames);
                                    for (int j = 0; j < tempNames.size(); j++) { out << tempNames[j] << endl; }
                                
                                }else { m->mothurOut("[ERROR]: parsing cannot find " + name + ".\n"); m->setControl_pressed(true); }
                            }
                            out.close();
                            in.close();
                            util.renameFile(thisaccnosFileName+".temp", thisaccnosFileName);
                        }else { m->mothurOut("[ERROR]: parsing cannot find " + fileGroup[thisFastaName] + ".\n"); m->setControl_pressed(true); }
                    }
                    
                }
            }

			//append files
			util.appendFiles(thisoutputFileName, outputFName); util.mothurRemove(thisoutputFileName); 
			util.appendFiles(thisaccnosFileName, accnos); util.mothurRemove(thisaccnosFileName);
			if (trim) { util.appendFiles(thistrimFastaFileName, fasta); util.mothurRemove(thistrimFastaFileName); }
			util.mothurRemove(thisFastaName);
			
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
//**********************************************************************************************************************

int ChimeraSlayerCommand::driver(string outputFName, string filename, string accnos, string fasta, map<string, int>& priority){
	try {
		
        if (m->getDebug()) { m->mothurOut("[DEBUG]: filename = " + filename + "\n"); }
        
		MothurChimera* chimera;
		if (templatefile != "self") { //you want to run slayer with a reference template
			chimera = new ChimeraSlayer(filename, templatefile, trim, search, ksize, match, mismatch, window, divR, minSimilarity, minCoverage, minBS, minSNP, parents, iters, increment, numwanted, realign, blastlocation, util.getRandomNumber());
		}else {
			chimera = new ChimeraSlayer(filename, templatefile, trim, priority, search, ksize, match, mismatch, window, divR, minSimilarity, minCoverage, minBS, minSNP, parents, iters, increment, numwanted, realign, blastlocation, util.getRandomNumber());	
		}
		
		if (m->getControl_pressed()) { delete chimera; return 0; }
		
		if (chimera->getUnaligned()) { delete chimera; m->mothurOut("Your template sequences are different lengths, please correct."); m->mothurOutEndLine(); m->setControl_pressed(true); return 0; }
		templateSeqsLength = chimera->getLength();
		
		ofstream out;
		util.openOutputFile(outputFName, out);
		
		ofstream out2;
		util.openOutputFile(accnos, out2);
		
		ofstream out3;
		if (trim) {  util.openOutputFile(fasta, out3); }
		
		ifstream inFASTA;
		util.openInputFile(filename, inFASTA);

		chimera->printHeader(out);

		int count = 0;
		while (!inFASTA.eof()) {
		
			if (m->getControl_pressed()) {	delete chimera; out.close(); out2.close(); if (trim) { out3.close(); } inFASTA.close(); return 1;	}
		
			Sequence* candidateSeq = new Sequence(inFASTA);  util.gobble(inFASTA);
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
		util.openInputFile(fastaFile, in);
		
		while (!in.eof()) {
			
			if (m->getControl_pressed()) { in.close(); return nameAbund; }
			
			Sequence seq(in); util.gobble(in);
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
        }else { error = util.readNames(nameFile, nameMapCount, seqs); }
		
		if (m->getControl_pressed()) { return nameAbund; }
		
		if (error == 1) { m->setControl_pressed(true); return nameAbund; }
		if (seqs.size() != nameMapCount.size()) { m->mothurOut( "The number of sequences in your fastafile does not match the number of sequences in your namefile, aborting."); m->mothurOutEndLine(); m->setControl_pressed(true); return nameAbund; }

		sort(nameMapCount.begin(), nameMapCount.end(), compareSeqPriorityNodes);
		
		string newFasta = fastaFile + ".temp";
		ofstream out;
		util.openOutputFile(newFasta, out);
		
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
				int num = util.getNumNames(itNameMap->second);
				
				seqPriorityNode temp(num, thisseqs[i].getAligned(), thisseqs[i].getName());
				nameVector.push_back(temp);
			}
		}
	
		//sort by num represented
		sort(nameVector.begin(), nameVector.end(), compareSeqPriorityNodes);
	
		if (m->getControl_pressed()) { return nameAbund; }
		
		if (thisseqs.size() != nameVector.size()) { m->mothurOut( "The number of sequences in your fastafile does not match the number of sequences in your namefile, aborting."); m->mothurOutEndLine(); m->setControl_pressed(true); return nameAbund; }
				
		ofstream out;
		util.openOutputFile(newFile, out);
		
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
		util.openOutputFile(newFile, out);
		
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

