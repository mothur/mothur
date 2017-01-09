/*
 *  chimerauchimecommand.cpp
 *  Mothur
 *
 *  Created by westcott on 5/13/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "chimerauchimecommand.h"
#include "deconvolutecommand.h"
#include "sequence.hpp"
#include "systemcommand.h"

//**********************************************************************************************************************
vector<string> ChimeraUchimeCommand::setParameters(){	
	try {
		CommandParameter ptemplate("reference", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(ptemplate);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","chimera-accnos",false,true,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
        CommandParameter pstrand("strand", "String", "", "", "", "", "","",false,false); parameters.push_back(pstrand);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		CommandParameter pabskew("abskew", "Number", "", "1.9", "", "", "","",false,false); parameters.push_back(pabskew);
		CommandParameter pchimealns("chimealns", "Boolean", "", "F", "", "", "","alns",false,false); parameters.push_back(pchimealns);
		CommandParameter pminh("minh", "Number", "", "0.3", "", "", "","",false,false); parameters.push_back(pminh);
		CommandParameter pmindiv("mindiv", "Number", "", "0.5", "", "", "","",false,false); parameters.push_back(pmindiv);
		CommandParameter pxn("xn", "Number", "", "8.0", "", "", "","",false,false); parameters.push_back(pxn);
		CommandParameter pdn("dn", "Number", "", "1.4", "", "", "","",false,false); parameters.push_back(pdn);
		CommandParameter pxa("xa", "Number", "", "1", "", "", "","",false,false); parameters.push_back(pxa);
		CommandParameter pchunks("chunks", "Number", "", "4", "", "", "","",false,false); parameters.push_back(pchunks);
		CommandParameter pminchunk("minchunk", "Number", "", "64", "", "", "","",false,false); parameters.push_back(pminchunk);
		CommandParameter pidsmoothwindow("idsmoothwindow", "Number", "", "32", "", "", "","",false,false); parameters.push_back(pidsmoothwindow);
        CommandParameter pdups("dereplicate", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pdups);

		//CommandParameter pminsmoothid("minsmoothid", "Number", "", "0.95", "", "", "",false,false); parameters.push_back(pminsmoothid);
		CommandParameter pmaxp("maxp", "Number", "", "2", "", "", "","",false,false); parameters.push_back(pmaxp);
		CommandParameter pskipgaps("skipgaps", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pskipgaps);
		CommandParameter pskipgaps2("skipgaps2", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pskipgaps2);
		CommandParameter pminlen("minlen", "Number", "", "10", "", "", "","",false,false); parameters.push_back(pminlen);
		CommandParameter pmaxlen("maxlen", "Number", "", "10000", "", "", "","",false,false); parameters.push_back(pmaxlen);
		CommandParameter pucl("ucl", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pucl);
		CommandParameter pqueryfract("queryfract", "Number", "", "0.5", "", "", "","",false,false); parameters.push_back(pqueryfract);

		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraUchimeCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraUchimeCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The chimera.uchime command reads a fastafile and referencefile and outputs potentially chimeric sequences.\n";
		helpString += "This command is a wrapper for uchime written by Robert C. Edgar.\n";
		helpString += "The chimera.uchime command parameters are fasta, name, count, reference, processors, dereplicate, abskew, chimealns, minh, mindiv, xn, dn, xa, chunks, minchunk, idsmoothwindow, minsmoothid, maxp, skipgaps, skipgaps2, minlen, maxlen, ucl, strand and queryfact.\n";
		helpString += "The fasta parameter allows you to enter the fasta file containing your potentially chimeric sequences, and is required, unless you have a valid current fasta file. \n";
		helpString += "The name parameter allows you to provide a name file, if you are using template=self. \n";
        helpString += "The count parameter allows you to provide a count file, if you are using template=self. When you use a count file with group info and dereplicate=T, mothur will create a *.pick.count_table file containing seqeunces after chimeras are removed. \n";
		helpString += "You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amazon.fasta \n";
		helpString += "The group parameter allows you to provide a group file. The group file can be used with a namesfile and reference=self. When checking sequences, only sequences from the same group as the query sequence will be used as the reference. \n";
        helpString += "If the dereplicate parameter is false, then if one group finds the seqeunce to be chimeric, then all groups find it to be chimeric, default=f.\n";
		helpString += "The reference parameter allows you to enter a reference file containing known non-chimeric sequences, and is required. You may also set template=self, in this case the abundant sequences will be used as potential parents. \n";
		helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n";
		helpString += "The abskew parameter can only be used with template=self. Minimum abundance skew. Default 1.9. Abundance skew is: min [ abund(parent1), abund(parent2) ] / abund(query).\n";
		helpString += "The chimealns parameter allows you to indicate you would like a file containing multiple alignments of query sequences to parents in human readable format. Alignments show columns with differences that support or contradict a chimeric model.\n";
		helpString += "The minh parameter - mininum score to report chimera. Default 0.3. Values from 0.1 to 5 might be reasonable. Lower values increase sensitivity but may report more false positives. If you decrease xn you may need to increase minh, and vice versa.\n";
		helpString += "The mindiv parameter - minimum divergence ratio, default 0.5. Div ratio is 100%% - %%identity between query sequence and the closest candidate for being a parent. If you don't care about very close chimeras, then you could increase mindiv to, say, 1.0 or 2.0, and also decrease minh, say to 0.1, to increase sensitivity. How well this works will depend on your data. Best is to tune parameters on a good benchmark.\n";
		helpString += "The xn parameter - weight of a no vote. Default 8.0. Decreasing this weight to around 3 or 4 may give better performance on denoised data.\n";
		helpString += "The dn parameter - pseudo-count prior on number of no votes. Default 1.4. Probably no good reason to change this unless you can retune to a good benchmark for your data. Reasonable values are probably in the range from 0.2 to 2.\n";
		helpString += "The xa parameter - weight of an abstain vote. Default 1. So far, results do not seem to be very sensitive to this parameter, but if you have a good training set might be worth trying. Reasonable values might range from 0.1 to 2.\n";
		helpString += "The chunks parameter is the number of chunks to extract from the query sequence when searching for parents. Default 4.\n";
		helpString += "The minchunk parameter is the minimum length of a chunk. Default 64.\n";
		helpString += "The idsmoothwindow parameter is the length of id smoothing window. Default 32.\n";
		//helpString += "The minsmoothid parameter - minimum factional identity over smoothed window of candidate parent. Default 0.95.\n";
		helpString += "The maxp parameter - maximum number of candidate parents to consider. Default 2. In tests so far, increasing maxp gives only a very small improvement in sensivity but tends to increase the error rate quite a bit.\n";
		helpString += "The skipgaps parameter controls how gapped columns affect counting of diffs. If skipgaps is set to T, columns containing gaps do not found as diffs. Default = T.\n";
		helpString += "The skipgaps2 parameter controls how gapped columns affect counting of diffs. If skipgaps2 is set to T, if column is immediately adjacent to a column containing a gap, it is not counted as a diff. Default = T.\n";
		helpString += "The minlen parameter is the minimum unaligned sequence length. Defaults 10. Applies to both query and reference sequences.\n";
		helpString += "The maxlen parameter is the maximum unaligned sequence length. Defaults 10000. Applies to both query and reference sequences.\n";
		helpString += "The ucl parameter - use local-X alignments. Default is global-X or false. On tests so far, global-X is always better; this option is retained because it just might work well on some future type of data.\n";
		helpString += "The queryfract parameter - minimum fraction of the query sequence that must be covered by a local-X alignment. Default 0.5. Applies only when ucl is true.\n";
		helpString += "The chimera.uchime command should be in the following format: \n";
		helpString += "chimera.uchime(fasta=yourFastaFile, reference=yourTemplate) \n";
		helpString += "Example: chimera.uchime(fasta=AD.align, reference=silva.gold.align) \n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraUchimeCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraUchimeCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "chimera") {  pattern = "[filename],[tag],uchime.chimeras"; }
        else if (type == "accnos") {  pattern = "[filename],[tag],uchime.accnos"; }
        else if (type == "alns") {  pattern = "[filename],[tag],uchime.alns"; }
        else if (type == "count") {  pattern = "[filename],[tag],uchime.pick.count_table"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraUchimeCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ChimeraUchimeCommand::ChimeraUchimeCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
		vector<string> tempOutNames;
		outputTypes["chimera"] = tempOutNames;
		outputTypes["accnos"] = tempOutNames;
		outputTypes["alns"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraUchimeCommand", "ChimeraUchimeCommand");
		exit(1);
	}
}
//***************************************************************************************************************
ChimeraUchimeCommand::ChimeraUchimeCommand(string option)  {
	try {
		abort = false; calledHelp = false; hasName=false; hasCount=false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter("chimera.uchime");
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			vector<string> tempOutNames;
			outputTypes["chimera"] = tempOutNames;
			outputTypes["accnos"] = tempOutNames;
			outputTypes["alns"] = tempOutNames;
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
			
            if (hasGroup && hasCount) { m->mothurOut("[ERROR]: You must enter ONLY ONE of the following: count or group."); m->mothurOutEndLine(); abort = true; }			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}
			
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}
			
			string path;
			it = parameters.find("reference");
			//user has given a template file
			if(it != parameters.end()){ 
				if (it->second == "self") { templatefile = "self"; }
				else {
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["reference"] = inputDir + it->second;		}
					
					templatefile = validParameter.validFile(parameters, "reference", true);
					if (templatefile == "not open") { abort = true; }
					else if (templatefile == "not found") { //check for saved reference sequences
						m->mothurOut("[ERROR]: The reference parameter is a required.\n"); abort = true;
					}
				}
			}else if (hasName) {  templatefile = "self"; }
            else if (hasCount) {  templatefile = "self"; }
			else {
                m->mothurOut("[ERROR]: The reference parameter is a required.");
					
                templatefile = ""; abort = true;
			}
				
			string temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
			
			abskew = validParameter.validFile(parameters, "abskew", false);	if (abskew == "not found"){	useAbskew = false;  abskew = "1.9";	}else{  useAbskew = true;  }
			if (useAbskew && templatefile != "self") { m->mothurOut("The abskew parameter is only valid with template=self, ignoring."); m->mothurOutEndLine(); useAbskew = false; }
			
			temp = validParameter.validFile(parameters, "chimealns", false);			if (temp == "not found") { temp = "f"; }
			chimealns = m->isTrue(temp); 
			
			minh = validParameter.validFile(parameters, "minh", false);						if (minh == "not found")			{ useMinH = false; minh = "0.3";					}	else{ useMinH = true;			}
			mindiv = validParameter.validFile(parameters, "mindiv", false);					if (mindiv == "not found")			{ useMindiv = false; mindiv = "0.5";				}	else{ useMindiv = true;			}
			xn = validParameter.validFile(parameters, "xn", false);							if (xn == "not found")				{ useXn = false; xn = "8.0";						}	else{ useXn = true;				}
			dn = validParameter.validFile(parameters, "dn", false);							if (dn == "not found")				{ useDn = false; dn = "1.4";						}	else{ useDn = true;				}
			xa = validParameter.validFile(parameters, "xa", false);							if (xa == "not found")				{ useXa = false; xa = "1";							}	else{ useXa = true;				}
			chunks = validParameter.validFile(parameters, "chunks", false);					if (chunks == "not found")			{ useChunks = false; chunks = "4";					}	else{ useChunks = true;			}
			minchunk = validParameter.validFile(parameters, "minchunk", false);				if (minchunk == "not found")		{ useMinchunk = false; minchunk = "64";				}	else{ useMinchunk = true;		}
			idsmoothwindow = validParameter.validFile(parameters, "idsmoothwindow", false);	if (idsmoothwindow == "not found")	{ useIdsmoothwindow = false; idsmoothwindow = "32";	}	else{ useIdsmoothwindow = true;	}
			//minsmoothid = validParameter.validFile(parameters, "minsmoothid", false);		if (minsmoothid == "not found")		{ useMinsmoothid = false; minsmoothid = "0.95";		}	else{ useMinsmoothid = true;	}
			maxp = validParameter.validFile(parameters, "maxp", false);						if (maxp == "not found")			{ useMaxp = false; maxp = "2";						}	else{ useMaxp = true;			}
			minlen = validParameter.validFile(parameters, "minlen", false);					if (minlen == "not found")			{ useMinlen = false; minlen = "10";					}	else{ useMinlen = true;			}
			maxlen = validParameter.validFile(parameters, "maxlen", false);					if (maxlen == "not found")			{ useMaxlen = false; maxlen = "10000";				}	else{ useMaxlen = true;			}
            
            strand = validParameter.validFile(parameters, "strand", false);	if (strand == "not found")	{  strand = "";	}
			
			temp = validParameter.validFile(parameters, "ucl", false);						if (temp == "not found") { temp = "f"; }
			ucl = m->isTrue(temp);
			
			queryfract = validParameter.validFile(parameters, "queryfract", false);			if (queryfract == "not found")		{ useQueryfract = false; queryfract = "0.5";		}	else{ useQueryfract = true;		}
			if (!ucl && useQueryfract) { m->mothurOut("queryfact may only be used when ucl=t, ignoring."); m->mothurOutEndLine(); useQueryfract = false; }
			
			temp = validParameter.validFile(parameters, "skipgaps", false);					if (temp == "not found") { temp = "t"; }
			skipgaps = m->isTrue(temp); 

			temp = validParameter.validFile(parameters, "skipgaps2", false);				if (temp == "not found") { temp = "t"; }
			skipgaps2 = m->isTrue(temp); 
            
            
			temp = validParameter.validFile(parameters, "dereplicate", false);	
			if (temp == "not found") { temp = "false";			}
			dups = m->isTrue(temp);

			
			if (hasName && (templatefile != "self")) { m->mothurOut("You have provided a namefile and the reference parameter is not set to self. I am not sure what reference you are trying to use, aborting."); m->mothurOutEndLine(); abort=true; }
            if (hasCount && (templatefile != "self")) { m->mothurOut("You have provided a countfile and the reference parameter is not set to self. I am not sure what reference you are trying to use, aborting."); m->mothurOutEndLine(); abort=true; }
			if (hasGroup && (templatefile != "self")) { m->mothurOut("You have provided a group file and the reference parameter is not set to self. I am not sure what reference you are trying to use, aborting."); m->mothurOutEndLine(); abort=true; }
			
			//look for uchime exe
			path = m->mothurProgramPath;
			
			string uchimeCommand;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
			uchimeCommand = path + "uchime";	//	format the database, -o option gives us the ability
            if (m->debug) { 
                m->mothurOut("[DEBUG]: Uchime location using \"which uchime\" = "); 
                Command* newCommand = new SystemCommand("which uchime"); m->mothurOutEndLine();
                newCommand->execute();
                delete newCommand;
                m->mothurOut("[DEBUG]: Mothur's location using \"which mothur\" = "); 
                newCommand = new SystemCommand("which mothur"); m->mothurOutEndLine();
                newCommand->execute();
                delete newCommand;
            }
#else
			uchimeCommand = path + "\\uchime.exe";
#endif
        
			//test to make sure uchime exists
			ifstream in;
			uchimeCommand = m->getFullPathName(uchimeCommand);
			int ableToOpen = m->openInputFile(uchimeCommand, in, "no error"); in.close();
			if(ableToOpen == 1) {	
                m->mothurOut(uchimeCommand + " file does not exist. Checking path... \n");
                //check to see if uchime is in the path??
                
                string uLocation = ""; 
                
                ifstream in2;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
                uLocation = m->findProgramPath("uchime");
                ableToOpen = m->openInputFile(uLocation, in2, "no error"); in2.close();
#else
                uLocation = m->findProgramPath("uchime.exe");
                ableToOpen = m->openInputFile(uLocation, in2, "no error"); in2.close();
#endif

                if(ableToOpen == 1) { m->mothurOut("[ERROR]: " + uLocation + " file does not exist. mothur requires the uchime executable."); m->mothurOutEndLine(); abort = true; } 
                else {  m->mothurOut("Found uchime in your path, using " + uLocation + "\n");uchimeLocation = uLocation; }
            }else {  uchimeLocation = uchimeCommand; }
            
            uchimeLocation = m->getFullPathName(uchimeLocation);
        }
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "ChimeraSlayerCommand");
		exit(1);
	}
}
//***************************************************************************************************************

int ChimeraUchimeCommand::execute(){
	try{
        
        if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		m->mothurOut("\nuchime by Robert C. Edgar\nhttp://drive5.com/uchime\nThis code is donated to the public domain.\n\n");
		
		for (int s = 0; s < fastaFileNames.size(); s++) {
			
			m->mothurOut("Checking sequences from " + fastaFileNames[s] + " ..." ); m->mothurOutEndLine();
			
			int start = time(NULL);	
			string nameFile = "";
			if (outputDir == "") { outputDir = m->hasPath(fastaFileNames[s]);  }//if user entered a file with a path then preserve it				
			map<string, string> variables; 
            variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s]));
            variables["[tag]"] = "denovo";
            if (templatefile != "self") { variables["[tag]"] = "ref"; }
			string outputFileName = getOutputFileName("chimera", variables);
			string accnosFileName = getOutputFileName("accnos", variables);
			string alnsFileName = getOutputFileName("alns", variables);
			string newFasta = m->getRootName(fastaFileNames[s]) + "temp";
            string newCountFile = "";
				
			//you provided a groupfile
			string groupFile = "";
            bool hasGroup = false;
			if (groupFileNames.size() != 0) { groupFile = groupFileNames[s]; hasGroup = true; }
            else if (hasCount) {
                CountTable ct;
                if (ct.testGroups(nameFileNames[s])) { hasGroup = true; }
                variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(nameFileNames[s]));
                newCountFile = getOutputFileName("count", variables);
            }
			
			if ((templatefile == "self") && (!hasGroup)) { //you want to run uchime with a template=self and no groups

				if (processors != 1) { m->mothurOut("When using template=self, mothur can only use 1 processor, continuing."); m->mothurOutEndLine(); processors = 1; }
				if (nameFileNames.size() != 0) { //you provided a namefile and we don't need to create one
					nameFile = nameFileNames[s];
				}else { nameFile = getNamesFile(fastaFileNames[s]); }
										
				map<string, string> seqs;  
				readFasta(fastaFileNames[s], seqs);  if (m->control_pressed) { for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0; }

				//read namefile
				vector<seqPriorityNode> nameMapCount;
                int error;
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
                }else {
                    error = m->readNames(nameFile, nameMapCount, seqs); if (m->control_pressed) { for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0; }
                }
				if (error == 1) { for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0; }
				if (seqs.size() != nameMapCount.size()) { m->mothurOut( "The number of sequences in your fastafile does not match the number of sequences in your namefile, aborting."); m->mothurOutEndLine(); for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0; }
				
                m->printVsearchFile(nameMapCount, newFasta, "/ab=", "/");
				fastaFileNames[s] = newFasta;
			}
			
			if (m->control_pressed) {  for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0;	}				
			
			if (hasGroup) {
				if (nameFileNames.size() != 0) { //you provided a namefile and we don't need to create one
					nameFile = nameFileNames[s];
				}else { nameFile = getNamesFile(fastaFileNames[s]); }
				
				//Parse sequences by group
                vector<string> groups;
                map<string, string> uniqueNames;
                vector<string> temp;
                if (hasCount) {
                    cparser = new SequenceCountParser(nameFile, fastaFileNames[s], temp);
                    groups = cparser->getNamesOfGroups();
                    uniqueNames = cparser->getAllSeqsMap();
                }else{
                    sparser = new SequenceParser(groupFile, fastaFileNames[s], nameFile, temp);
                    groups = sparser->getNamesOfGroups();
                    uniqueNames = sparser->getAllSeqsMap();
                }
					
				if (m->control_pressed) { for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0; }
								
				//clears files
				ofstream out, out1, out2;
				m->openOutputFile(outputFileName, out); out.close(); 
				m->openOutputFile(accnosFileName, out1); out1.close();
				if (chimealns) { m->openOutputFile(alnsFileName, out2); out2.close(); }
				int totalSeqs = 0;
				
				if(processors == 1)	{	totalSeqs = driverGroups(outputFileName, newFasta, accnosFileName, alnsFileName, newCountFile, 0, groups.size(), groups);
                    
                    if (hasCount && dups) {
                        CountTable c; c.readTable(nameFile, true, false);
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

                }else{	totalSeqs = createProcessesGroups(outputFileName, newFasta, accnosFileName, alnsFileName, newCountFile, groups, nameFile, groupFile, fastaFileNames[s]);			}

				if (m->control_pressed) {  for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0;	}				
               
                
                if (!dups) { 
                    int totalChimeras = deconvoluteResults(uniqueNames, outputFileName, accnosFileName, alnsFileName);
				
                    m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(totalSeqs) + " sequences. " + toString(totalChimeras) + " chimeras were found.");	m->mothurOutEndLine();
                    m->mothurOut("The number of sequences checked may be larger than the number of unique sequences because some sequences are found in several samples."); m->mothurOutEndLine(); 
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
                
                if (hasCount) { delete cparser; }
                else { delete sparser; }
                
				if (m->control_pressed) {  for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0;	}				
					
			}else{
				if (m->control_pressed) {  for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0;	}
			
				int numSeqs = 0;
				int numChimeras = 0;

				if(processors == 1){ numSeqs = driver(outputFileName, fastaFileNames[s], accnosFileName, alnsFileName, numChimeras); }
				else{	numSeqs = createProcesses(outputFileName, fastaFileNames[s], accnosFileName, alnsFileName, numChimeras); }
				
				//add headings
				ofstream out;
				m->openOutputFile(outputFileName+".temp", out); 
				out << "Score\tQuery\tParentA\tParentB\tIdQM\tIdQA\tIdQB\tIdAB\tIdQT\tLY\tLN\tLA\tRY\tRN\tRA\tDiv\tYN\n";
				out.close();
				
				m->appendFiles(outputFileName, outputFileName+".temp");
				m->mothurRemove(outputFileName); rename((outputFileName+".temp").c_str(), outputFileName.c_str());
				
				if (m->control_pressed) { for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	} return 0; }
			
				//remove file made for uchime
				if (templatefile == "self") {  m->mothurRemove(fastaFileNames[s]); }
			
				m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences. " + toString(numChimeras) + " chimeras were found.");	m->mothurOutEndLine();
			}
			
			outputNames.push_back(outputFileName); outputTypes["chimera"].push_back(outputFileName);
			outputNames.push_back(accnosFileName); outputTypes["accnos"].push_back(accnosFileName);
			if (chimealns) { outputNames.push_back(alnsFileName); outputTypes["alns"].push_back(alnsFileName); }
		}
	
		//set accnos file as new current accnosfile
		string current = "";
		itTypes = outputTypes.find("accnos");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setAccnosFile(current); }
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
		m->errorOut(e, "ChimeraUchimeCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int ChimeraUchimeCommand::deconvoluteResults(map<string, string>& uniqueNames, string outputFileName, string accnosFileName, string alnsFileName){
	try {
		map<string, string>::iterator itUnique;
		int total = 0;
		
		ofstream out2;
		m->openOutputFile(accnosFileName+".temp", out2);
		
		string name;
		set<string> namesInFile; //this is so if a sequence is found to be chimera in several samples we dont write it to the results file more than once
		set<string>::iterator itNames;
		set<string> chimerasInFile;
		set<string>::iterator itChimeras;

        if (!m->isBlank(accnosFileName)) {
            //edit accnos file
            ifstream in2;
            m->openInputFile(accnosFileName, in2);
            
            while (!in2.eof()) {
                if (m->control_pressed) { in2.close(); out2.close(); m->mothurRemove(outputFileName); m->mothurRemove((accnosFileName+".temp")); return 0; }
                
                in2 >> name; m->gobble(in2);
                
                //find unique name
                itUnique = uniqueNames.find(name);
                
                if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing accnos results. Cannot find " + name + "."); m->mothurOutEndLine(); m->control_pressed = true; }
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
        }
		out2.close();
		
		m->mothurRemove(accnosFileName);
		rename((accnosFileName+".temp").c_str(), accnosFileName.c_str());
		
		
		
		//edit chimera file
		ifstream in; 
		m->openInputFile(outputFileName, in);
		
		ofstream out;
		m->openOutputFile(outputFileName+".temp", out); out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		out << "Score\tQuery\tParentA\tParentB\tIdQM\tIdQA\tIdQB\tIdAB\tIdQT\tLY\tLN\tLA\tRY\tRN\tRA\tDiv\tYN\n";
		
		float temp1;
		string parent1, parent2, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, flag;
		name = "";
		namesInFile.clear();	
		//assumptions - in file each read will always look like - if uchime source is updated, revisit this code.
		/*										1	2	3	4	5	6	7	8	9	10	11	12	13	14	15
		 0.000000	F11Fcsw_33372/ab=18/		*	*	*	*	*	*	*	*	*	*	*	*	*	*	N
		 0.018300	F11Fcsw_14980/ab=16/		F11Fcsw_1915/ab=35/	F11Fcsw_6032/ab=42/	79.9	78.7	78.2	78.7	79.2	3	0	5	11	10	20	1.46	N
		*/
		
		while (!in.eof()) {
			
			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove((outputFileName+".temp")); return 0; }
			
			bool print = false;
			in >> temp1;	m->gobble(in);
			in >> name;		m->gobble(in);
			in >> parent1;	m->gobble(in);
			in >> parent2;	m->gobble(in);
			in >> temp2 >> temp3 >> temp4 >> temp5 >> temp6 >> temp7 >> temp8 >> temp9 >> temp10 >> temp11 >> temp12 >> temp13 >> flag;
			m->gobble(in);
			
			//parse name - name will look like U68590/ab=1/
			string restOfName = "";
			int pos = name.find_first_of('/');
			if (pos != string::npos) {
				restOfName = name.substr(pos);
				name = name.substr(0, pos);
			}
			
			//find unique name
			itUnique = uniqueNames.find(name);
			
			if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find "+ name + "."); m->mothurOutEndLine(); m->control_pressed = true; }
			else {
				name = itUnique->second;
				//is this name already in the file
				itNames = namesInFile.find((name));
				
				if (itNames == namesInFile.end()) { //no not in file
					if (flag == "N") { //are you really a no??
						//is this sequence really not chimeric??
						itChimeras = chimerasInFile.find(name);
						
						//then you really are a no so print, otherwise skip
						if (itChimeras == chimerasInFile.end()) { print = true; }
					}else{ print = true; }
				}
			}
			
			if (print) {
				out << temp1 << '\t' << name << restOfName << '\t';
				namesInFile.insert(name);
				
				//parse parent1 names
				if (parent1 != "*") {
					restOfName = "";
					pos = parent1.find_first_of('/');
					if (pos != string::npos) {
						restOfName = parent1.substr(pos);
						parent1 = parent1.substr(0, pos);
					}
					
					itUnique = uniqueNames.find(parent1);
					if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find parentA "+ parent1 + "."); m->mothurOutEndLine(); m->control_pressed = true; }
					else {	out << itUnique->second << restOfName << '\t';	}
				}else { out << parent1 << '\t'; }
				
				//parse parent2 names
				if (parent2 != "*") {
					restOfName = "";
					pos = parent2.find_first_of('/');
					if (pos != string::npos) {
						restOfName = parent2.substr(pos);
						parent2 = parent2.substr(0, pos);
					}
					
					itUnique = uniqueNames.find(parent2);
					if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find parentB "+ parent2 + "."); m->mothurOutEndLine(); m->control_pressed = true; }
					else {	out << itUnique->second << restOfName << '\t';	}
				}else { out << parent2 << '\t'; }
				
				out << temp2 << '\t' << temp3 << '\t' << temp4 << '\t' << temp5 << '\t' << temp6 << '\t' << temp7 << '\t' << temp8 << '\t' << temp9 << '\t' << temp10 << '\t' << temp11 << '\t' << temp12 << temp13 << '\t' << flag << endl;	
			}
		}
		in.close();
		out.close();
		
		m->mothurRemove(outputFileName);
		rename((outputFileName+".temp").c_str(), outputFileName.c_str());
		
				
		//edit anls file
		//assumptions - in file each read will always look like - if uchime source is updated, revisit this code.
		/*
		 ------------------------------------------------------------------------
		 Query   (  179 nt) F21Fcsw_11639/ab=591/
		 ParentA (  179 nt) F11Fcsw_6529/ab=1625/
		 ParentB (  181 nt) F21Fcsw_12128/ab=1827/
		 
		 A     1 AAGgAAGAtTAATACaagATGgCaTCatgAGtccgCATgTtcAcatGATTAAAG--gTaTtcCGGTagacGATGGGGATG 78
		 Q     1 AAGTAAGACTAATACCCAATGACGTCTCTAGAAGACATCTGAAAGAGATTAAAG--ATTTATCGGTGATGGATGGGGATG 78
		 B     1 AAGgAAGAtTAATcCaggATGggaTCatgAGttcACATgTccgcatGATTAAAGgtATTTtcCGGTagacGATGGGGATG 80
		 Diffs      N    N    A N?N   N N  NNN  N?NB   N ?NaNNN          B B NN    NNNN          
		 Votes      0    0    + 000   0 0  000  000+   0 00!000            + 00    0000          
		 Model   AAAAAAAAAAAAAAAAAAAAAAxBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
		 
		 A    79 CGTtccATTAGaTaGTaGGCGGGGTAACGGCCCACCtAGtCttCGATggaTAGGGGTTCTGAGAGGAAGGTCCCCCACAT 158
		 Q    79 CGTCTGATTAGCTTGTTGGCGGGGTAACGGCCCACCAAGGCAACGATCAGTAGGGGTTCTGAGAGGAAGGTCCCCCACAT 158
		 B    81 CGTtccATTAGaTaGTaGGCGGGGTAACGGCCCACCtAGtCAACGATggaTAGGGGTTCTGAGAGGAAGGTCCCCCACAT 160
		 Diffs      NNN     N N  N                   N  N BB    NNN                              
		 Votes      000     0 0  0                   0  0 ++    000                              
		 Model   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
		 
		 A   159 TGGAACTGAGACACGGTCCAA 179
		 Q   159 TGGAACTGAGACACGGTCCAA 179
		 B   161 TGGAACTGAGACACGGTCCAA 181
		 Diffs                        
		 Votes                        
		 Model   BBBBBBBBBBBBBBBBBBBBB
		 
		 Ids.  QA 76.6%, QB 77.7%, AB 93.7%, QModel 78.9%, Div. +1.5%
		 Diffs Left 7: N 0, A 6, Y 1 (14.3%); Right 35: N 1, A 30, Y 4 (11.4%), Score 0.0047
		*/
		if (chimealns) {
			ifstream in3; 
			m->openInputFile(alnsFileName, in3);
		
			ofstream out3;
			m->openOutputFile(alnsFileName+".temp", out3); out3.setf(ios::fixed, ios::floatfield); out3.setf(ios::showpoint);
		
			name = "";
			namesInFile.clear();
			string line = "";
			
			while (!in3.eof()) {
				if (m->control_pressed) { in3.close(); out3.close(); m->mothurRemove(outputFileName); m->mothurRemove((accnosFileName)); m->mothurRemove((alnsFileName+".temp")); return 0; }
				
				line = "";
				line = m->getline(in3); 
				string temp = "";
				
				if (line != "") {
					istringstream iss(line);
					iss >> temp;
					
					//are you a name line
					if ((temp == "Query") || (temp == "ParentA") || (temp == "ParentB")) {
						int spot = 0;
						for (int i = 0; i < line.length(); i++) {
							spot = i;
							if (line[i] == ')') { break; }
							else { out3 << line[i]; }
						}
						
						if (spot == (line.length() - 1)) { m->mothurOut("[ERROR]: could not line sequence name in line " + line + "."); m->mothurOutEndLine(); m->control_pressed = true; }
						else if ((spot+2) > (line.length() - 1)) { m->mothurOut("[ERROR]: could not line sequence name in line " + line + "."); m->mothurOutEndLine(); m->control_pressed = true; }
						else {
							out << line[spot] << line[spot+1];
							
							name = line.substr(spot+2);
							
							//parse name - name will either look like U68590/ab=1/ or U68590
							string restOfName = "";
							int pos = name.find_first_of('/');
							if (pos != string::npos) {
								restOfName = name.substr(pos);
								name = name.substr(0, pos);
							}
							
							//find unique name
							itUnique = uniqueNames.find(name);
							
							if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing alns results. Cannot find "+ name + "."); m->mothurOutEndLine();m->control_pressed = true;  }
							else {
								//only limit repeats on query names
								if (temp == "Query") {
									itNames = namesInFile.find((itUnique->second));
									
									if (itNames == namesInFile.end()) {
										out << itUnique->second << restOfName << endl;
										namesInFile.insert((itUnique->second));
									}
								}else { out << itUnique->second << restOfName << endl;  }
							}
							
						}
						
					}else { //not need to alter line
						out3 << line << endl;
					}
				}else { out3 << endl; }
			}
			in3.close();
			out3.close();
			
			m->mothurRemove(alnsFileName);
			rename((alnsFileName+".temp").c_str(), alnsFileName.c_str());
		}
		
		return total;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraUchimeCommand", "deconvoluteResults");
		exit(1);
	}
}		
//**********************************************************************************************************************
int ChimeraUchimeCommand::readFasta(string filename, map<string, string>& seqs){
	try {
		//create input file for uchime
		//read through fastafile and store info
		ifstream in;
		m->openInputFile(filename, in);
		
		while (!in.eof()) {
			
			if (m->control_pressed) { in.close(); return 0; }
			
			Sequence seq(in); m->gobble(in);
			seqs[seq.getName()] = seq.getAligned();
		}
		in.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraUchimeCommand", "readFasta");
		exit(1);
	}
}	
//**********************************************************************************************************************

string ChimeraUchimeCommand::getNamesFile(string& inputFile){
	try {
		string nameFile = "";
		
		m->mothurOutEndLine(); m->mothurOut("No namesfile given, running unique.seqs command to generate one."); m->mothurOutEndLine(); m->mothurOutEndLine();
		
		//use unique.seqs to create new name and fastafile
		string inputString = "fasta=" + inputFile;
		m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
		m->mothurOut("Running command: unique.seqs(" + inputString + ")"); m->mothurOutEndLine(); 
		m->mothurCalling = true;
        
		Command* uniqueCommand = new DeconvoluteCommand(inputString);
		uniqueCommand->execute();
		
		map<string, vector<string> > filenames = uniqueCommand->getOutputFiles();
		
		delete uniqueCommand;
		m->mothurCalling = false;
		m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
		
		nameFile = filenames["name"][0];
		inputFile = filenames["fasta"][0];
		
		return nameFile;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraUchimeCommand", "getNamesFile");
		exit(1);
	}
}
//**********************************************************************************************************************
int ChimeraUchimeCommand::driverGroups(string outputFName, string filename, string accnos, string alns, string countlist, int start, int end, vector<string> groups){
	try {
		
		int totalSeqs = 0;
		int numChimeras = 0;
        
        
        ofstream outCountList;
        if (hasCount && dups) { m->openOutputFile(countlist, outCountList); }
        
		for (int i = start; i < end; i++) {
			int start = time(NULL);	 if (m->control_pressed) {  outCountList.close(); m->mothurRemove(countlist); return 0; }
            
			int error;
            if (hasCount) { error = cparser->getSeqs(groups[i], filename, "/ab=", "/", true); if ((error == 1) || m->control_pressed) {  return 0; } }
            else { error = sparser->getSeqs(groups[i], filename, "/ab=", "/", true); if ((error == 1) || m->control_pressed) {  return 0; } }
			
			int numSeqs = driver((outputFName + groups[i]), filename, (accnos+groups[i]), (alns+ groups[i]), numChimeras);
			totalSeqs += numSeqs;
			
			if (m->control_pressed) { return 0; }
			
			//remove file made for uchime
			if (!m->debug) {  m->mothurRemove(filename);  }
            else { m->mothurOut("[DEBUG]: saving file: " + filename + ".\n"); }
			
            //if we provided a count file with group info and set dereplicate=t, then we want to create a *.pick.count_table
            //This table will zero out group counts for seqs determined to be chimeric by that group.
            if (dups) {
                if (!m->isBlank(accnos+groups[i])) {
                    ifstream in;
                    m->openInputFile(accnos+groups[i], in);
                    string name;
                    if (hasCount) {
                        while (!in.eof()) {
                            in >> name; m->gobble(in);
                            outCountList << name << '\t' << groups[i] << endl;
                        }
                        in.close();
                    }else {
                        map<string, string> thisnamemap = sparser->getNameMap(groups[i]);
                        map<string, string>::iterator itN;
                        ofstream out;
                        m->openOutputFile(accnos+groups[i]+".temp", out);
                        while (!in.eof()) {
                            in >> name; m->gobble(in); 
                            itN = thisnamemap.find(name);
                            if (itN != thisnamemap.end()) {
                                vector<string> tempNames; m->splitAtComma(itN->second, tempNames); 
                                for (int j = 0; j < tempNames.size(); j++) { out << tempNames[j] << endl; }
                                
                            }else { m->mothurOut("[ERROR]: parsing cannot find " + name + ".\n"); m->control_pressed = true; }
                        }
                        out.close();
                        in.close();
                        m->renameFile(accnos+groups[i]+".temp", accnos+groups[i]);
                    }
                   
                }
            }
            
			//append files
			m->appendFiles((outputFName+groups[i]), outputFName); m->mothurRemove((outputFName+groups[i]));
			m->appendFiles((accnos+groups[i]), accnos); m->mothurRemove((accnos+groups[i]));
			if (chimealns) { m->appendFiles((alns+groups[i]), alns); m->mothurRemove((alns+groups[i])); }
			
			m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences from group " + groups[i] + ".");	m->mothurOutEndLine();					
		}

        if (hasCount && dups) { outCountList.close(); }
        
        return totalSeqs;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraUchimeCommand", "driverGroups");
		exit(1);
	}
}	
//**********************************************************************************************************************

int ChimeraUchimeCommand::driver(string outputFName, string filename, string accnos, string alns, int& numChimeras){
	try {
		
		outputFName = m->getFullPathName(outputFName);
		filename = m->getFullPathName(filename);
		alns = m->getFullPathName(alns);
		
		//to allow for spaces in the path
		outputFName = "\"" + outputFName + "\"";
		filename = "\"" + filename + "\"";
		alns = "\"" + alns + "\"";
        
        if (filename.length() > 257) {
            m->mothurOut("[ERROR]: " + filename + " filename is " + toString(filename.length()) + " long. The uchime program can't handle files with a full path longer than 257 characters, please correct.\n"); m->control_pressed = true; return 0;
        }else if ((alns.length() > 257) && (chimealns)) {
            m->mothurOut("[ERROR]: " + alns + " filename is " + toString(alns.length()) + " long. The uchime program can't handle files with a full path longer than 257 characters, please correct.\n"); m->control_pressed = true; return 0;
        }else if (outputFName.length() > 257) {
            m->mothurOut("[ERROR]: " + outputFName + " filename is " + toString(outputFName.length()) + " long. The uchime program can't handle files with a full path longer than 257 characters, please correct input file name.\n"); m->control_pressed = true; return 0;
        }
				
		vector<char*> cPara;
		
		string uchimeCommand = uchimeLocation;
        uchimeCommand = "\"" + uchimeCommand + "\" ";
        
        char* tempUchime;
		tempUchime= new char[uchimeCommand.length()+1]; 
		*tempUchime = '\0';
		strncat(tempUchime, uchimeCommand.c_str(), uchimeCommand.length());
		cPara.push_back(tempUchime);
		
        //are you using a reference file
		if (templatefile != "self") {
            string outputFileName = filename.substr(1, filename.length()-2) + ".uchime_formatted";
            prepFile(filename.substr(1, filename.length()-2), outputFileName);
            filename = outputFileName;
            filename = "\"" + filename + "\"";
			//add reference file
			char* tempRef = new char[5]; 
			//strcpy(tempRef, "--db"); 
			*tempRef = '\0'; strncat(tempRef, "--db", 4);
			cPara.push_back(tempRef);  
			char* tempR = new char[templatefile.length()+1];
			//strcpy(tempR, templatefile.c_str());
			*tempR = '\0'; strncat(tempR, templatefile.c_str(), templatefile.length());
			cPara.push_back(tempR);
		}
		
		char* tempIn = new char[8]; 
		*tempIn = '\0'; strncat(tempIn, "--input", 7);
		//strcpy(tempIn, "--input"); 
		cPara.push_back(tempIn);
		char* temp = new char[filename.length()+1];
		*temp = '\0'; strncat(temp, filename.c_str(), filename.length());
		//strcpy(temp, filename.c_str());
		cPara.push_back(temp);
		
		char* tempO = new char[12]; 
		*tempO = '\0'; strncat(tempO, "--uchimeout", 11);
		//strcpy(tempO, "--uchimeout"); 
		cPara.push_back(tempO);
		char* tempout = new char[outputFName.length()+1];
		//strcpy(tempout, outputFName.c_str());
		*tempout = '\0'; strncat(tempout, outputFName.c_str(), outputFName.length());
		cPara.push_back(tempout);
		
		if (chimealns) {
			char* tempA = new char[13]; 
			*tempA = '\0'; strncat(tempA, "--uchimealns", 12);
			//strcpy(tempA, "--uchimealns"); 
			cPara.push_back(tempA);
			char* tempa = new char[alns.length()+1];
			//strcpy(tempa, alns.c_str());
			*tempa = '\0'; strncat(tempa, alns.c_str(), alns.length());
			cPara.push_back(tempa);
		}
        
        if (strand != "") {
			char* tempA = new char[9]; 
			*tempA = '\0'; strncat(tempA, "--strand", 8);
			cPara.push_back(tempA);
			char* tempa = new char[strand.length()+1];
			*tempa = '\0'; strncat(tempa, strand.c_str(), strand.length());
			cPara.push_back(tempa);
		}
		
		if (useAbskew) {
			char* tempskew = new char[9];
			*tempskew = '\0'; strncat(tempskew, "--abskew", 8);
			//strcpy(tempskew, "--abskew"); 
			cPara.push_back(tempskew);
			char* tempSkew = new char[abskew.length()+1];
			//strcpy(tempSkew, abskew.c_str());
			*tempSkew = '\0'; strncat(tempSkew, abskew.c_str(), abskew.length());
			cPara.push_back(tempSkew);
		}
		
		if (useMinH) {
			char* tempminh = new char[7]; 
			*tempminh = '\0'; strncat(tempminh, "--minh", 6);
			//strcpy(tempminh, "--minh"); 
			cPara.push_back(tempminh);
			char* tempMinH = new char[minh.length()+1];
			*tempMinH = '\0'; strncat(tempMinH, minh.c_str(), minh.length());
			//strcpy(tempMinH, minh.c_str());
			cPara.push_back(tempMinH);
		}
		
		if (useMindiv) {
			char* tempmindiv = new char[9]; 
			*tempmindiv = '\0'; strncat(tempmindiv, "--mindiv", 8);
			//strcpy(tempmindiv, "--mindiv"); 
			cPara.push_back(tempmindiv);
			char* tempMindiv = new char[mindiv.length()+1];
			*tempMindiv = '\0'; strncat(tempMindiv, mindiv.c_str(), mindiv.length());
			//strcpy(tempMindiv, mindiv.c_str());
			cPara.push_back(tempMindiv);
		}
		
		if (useXn) {
			char* tempxn = new char[5]; 
			//strcpy(tempxn, "--xn"); 
			*tempxn = '\0'; strncat(tempxn, "--xn", 4);
			cPara.push_back(tempxn);
			char* tempXn = new char[xn.length()+1];
			//strcpy(tempXn, xn.c_str());
			*tempXn = '\0'; strncat(tempXn, xn.c_str(), xn.length());
			cPara.push_back(tempXn);
		}
		
		if (useDn) {
			char* tempdn = new char[5]; 
			//strcpy(tempdn, "--dn"); 
			*tempdn = '\0'; strncat(tempdn, "--dn", 4);
			cPara.push_back(tempdn);
			char* tempDn = new char[dn.length()+1];
			*tempDn = '\0'; strncat(tempDn, dn.c_str(), dn.length());
			//strcpy(tempDn, dn.c_str());
			cPara.push_back(tempDn);
		}
		
		if (useXa) {
			char* tempxa = new char[5]; 
			//strcpy(tempxa, "--xa"); 
			*tempxa = '\0'; strncat(tempxa, "--xa", 4);
			cPara.push_back(tempxa);
			char* tempXa = new char[xa.length()+1];
			*tempXa = '\0'; strncat(tempXa, xa.c_str(), xa.length());
			//strcpy(tempXa, xa.c_str());
			cPara.push_back(tempXa);
		}
		
		if (useChunks) {
			char* tempchunks = new char[9]; 
			//strcpy(tempchunks, "--chunks"); 
			*tempchunks = '\0'; strncat(tempchunks, "--chunks", 8);
			cPara.push_back(tempchunks);
			char* tempChunks = new char[chunks.length()+1];
			*tempChunks = '\0'; strncat(tempChunks, chunks.c_str(), chunks.length());
			//strcpy(tempChunks, chunks.c_str());
			cPara.push_back(tempChunks);
		}
		
		if (useMinchunk) {
			char* tempminchunk = new char[11]; 
			//strcpy(tempminchunk, "--minchunk"); 
			*tempminchunk = '\0'; strncat(tempminchunk, "--minchunk", 10);
			cPara.push_back(tempminchunk);
			char* tempMinchunk = new char[minchunk.length()+1];
			*tempMinchunk = '\0'; strncat(tempMinchunk, minchunk.c_str(), minchunk.length());
			//strcpy(tempMinchunk, minchunk.c_str());
			cPara.push_back(tempMinchunk);
		}
		
		if (useIdsmoothwindow) {
			char* tempidsmoothwindow = new char[17]; 
			*tempidsmoothwindow = '\0'; strncat(tempidsmoothwindow, "--idsmoothwindow", 16);
			//strcpy(tempidsmoothwindow, "--idsmoothwindow"); 
			cPara.push_back(tempidsmoothwindow);
			char* tempIdsmoothwindow = new char[idsmoothwindow.length()+1];
			*tempIdsmoothwindow = '\0'; strncat(tempIdsmoothwindow, idsmoothwindow.c_str(), idsmoothwindow.length());
			//strcpy(tempIdsmoothwindow, idsmoothwindow.c_str());
			cPara.push_back(tempIdsmoothwindow);
		}
		
		/*if (useMinsmoothid) {
			char* tempminsmoothid = new char[14]; 
			//strcpy(tempminsmoothid, "--minsmoothid"); 
			*tempminsmoothid = '\0'; strncat(tempminsmoothid, "--minsmoothid", 13);
			cPara.push_back(tempminsmoothid);
			char* tempMinsmoothid = new char[minsmoothid.length()+1];
			*tempMinsmoothid = '\0'; strncat(tempMinsmoothid, minsmoothid.c_str(), minsmoothid.length());
			//strcpy(tempMinsmoothid, minsmoothid.c_str());
			cPara.push_back(tempMinsmoothid);
		}*/
		
		if (useMaxp) {
			char* tempmaxp = new char[7]; 
			//strcpy(tempmaxp, "--maxp"); 
			*tempmaxp = '\0'; strncat(tempmaxp, "--maxp", 6);
			cPara.push_back(tempmaxp);
			char* tempMaxp = new char[maxp.length()+1];
			*tempMaxp = '\0'; strncat(tempMaxp, maxp.c_str(), maxp.length());
			//strcpy(tempMaxp, maxp.c_str());
			cPara.push_back(tempMaxp);
		}
		
		if (!skipgaps) {
			char* tempskipgaps = new char[13]; 
			//strcpy(tempskipgaps, "--[no]skipgaps");
			*tempskipgaps = '\0'; strncat(tempskipgaps, "--noskipgaps", 12);
			cPara.push_back(tempskipgaps);
		}
		
		if (!skipgaps2) {
			char* tempskipgaps2 = new char[14]; 
			//strcpy(tempskipgaps2, "--[no]skipgaps2"); 
			*tempskipgaps2 = '\0'; strncat(tempskipgaps2, "--noskipgaps2", 13);
			cPara.push_back(tempskipgaps2);
		}
		
		if (useMinlen) {
			char* tempminlen = new char[9]; 
			*tempminlen = '\0'; strncat(tempminlen, "--minlen", 8);
			//strcpy(tempminlen, "--minlen"); 
			cPara.push_back(tempminlen);
			char* tempMinlen = new char[minlen.length()+1];
			//strcpy(tempMinlen, minlen.c_str());
			*tempMinlen = '\0'; strncat(tempMinlen, minlen.c_str(), minlen.length());
			cPara.push_back(tempMinlen);
		}
		
		if (useMaxlen) {
			char* tempmaxlen = new char[9]; 
			//strcpy(tempmaxlen, "--maxlen"); 
			*tempmaxlen = '\0'; strncat(tempmaxlen, "--maxlen", 8);
			cPara.push_back(tempmaxlen);
			char* tempMaxlen = new char[maxlen.length()+1];
			*tempMaxlen = '\0'; strncat(tempMaxlen, maxlen.c_str(), maxlen.length());
			//strcpy(tempMaxlen, maxlen.c_str());
			cPara.push_back(tempMaxlen);
		}
		
		if (ucl) {
			char* tempucl = new char[5]; 
			strcpy(tempucl, "--ucl"); 
			cPara.push_back(tempucl);
		}
		
		if (useQueryfract) {
			char* tempqueryfract = new char[13]; 
			*tempqueryfract = '\0'; strncat(tempqueryfract, "--queryfract", 12);
			//strcpy(tempqueryfract, "--queryfract"); 
			cPara.push_back(tempqueryfract);
			char* tempQueryfract = new char[queryfract.length()+1];
			*tempQueryfract = '\0'; strncat(tempQueryfract, queryfract.c_str(), queryfract.length());
			//strcpy(tempQueryfract, queryfract.c_str());
			cPara.push_back(tempQueryfract);
		}
		
		
		char** uchimeParameters;
		uchimeParameters = new char*[cPara.size()];
		string commandString = "";
		for (int i = 0; i < cPara.size(); i++) {  uchimeParameters[i] = cPara[i];  commandString += toString(cPara[i]) + " "; } 
		//int numArgs = cPara.size();
		
		//uchime_main(numArgs, uchimeParameters); 
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
		commandString = "\"" + commandString + "\"";
#endif
        if (m->debug) { m->mothurOut("[DEBUG]: uchime command = " + commandString + ".\n"); }
		system(commandString.c_str());
		
		//free memory
		for(int i = 0; i < cPara.size(); i++)  {  delete cPara[i];  }
		delete[] uchimeParameters; 
		
		//remove "" from filenames
		outputFName = outputFName.substr(1, outputFName.length()-2);
		filename = filename.substr(1, filename.length()-2);
		alns = alns.substr(1, alns.length()-2);
		
		if (m->control_pressed) { return 0; }
		
		//create accnos file from uchime results
		ifstream in; 
		m->openInputFile(outputFName, in);
		
		ofstream out;
		m->openOutputFile(accnos, out);
		
		int num = 0;
		numChimeras = 0;
		while(!in.eof()) {
			
			if (m->control_pressed) { break; }
			
			string name = "";
			string chimeraFlag = "";
			//in >> chimeraFlag >> name;
			
            string line = m->getline(in);
            vector<string> pieces = m->splitWhiteSpace(line);
            if (pieces.size() > 2) { 
                name = pieces[1];
                //fix name if needed
                if (templatefile == "self") { 
                    name = name.substr(0, name.length()-1); //rip off last /
                    name = name.substr(0, name.find_last_of('/'));
                }
                
                chimeraFlag = pieces[pieces.size()-1];
			}
			//for (int i = 0; i < 15; i++) {  in >> chimeraFlag; }
			m->gobble(in);
			
			if (chimeraFlag == "Y") {  out << name << endl; numChimeras++; }
			num++;
		}
		in.close();
		out.close();
		
        //if (templatefile != "self") {  m->mothurRemove(filename); }
        
		return num;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraUchimeCommand", "driver");
		exit(1);
	}
}
/**************************************************************************************************/
//uchime can't handle some of the things allowed in mothurs fasta files. This functions "cleans up" the file.
int ChimeraUchimeCommand::prepFile(string filename, string output) {
	try {
        
        ifstream in;
        m->openInputFile(filename, in);
        
        ofstream out;
        m->openOutputFile(output, out);
        
        while (!in.eof()) {
            if (m->control_pressed) { break;  }
            
            Sequence seq(in); m->gobble(in);
            
            if (seq.getName() != "") { seq.printSequence(out); }
        }
        in.close();
        out.close();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "ChimeraUchimeCommand", "prepFile");
		exit(1);
	}
}
/**************************************************************************************************/

int ChimeraUchimeCommand::createProcesses(string outputFileName, string filename, string accnos, string alns, int& numChimeras) {
	try {
		
		processIDS.clear();
		int process = 1;
		int num = 0;
		vector<string> files;
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)		
		//break up file into multiple files
		m->divideFile(filename, processors, files);
		
		if (m->control_pressed) {  return 0;  }
				
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driver(outputFileName + toString(m->mothurGetpid(process)) + ".temp", files[process], accnos + toString(m->mothurGetpid(process)) + ".temp", alns + toString(m->mothurGetpid(process)) + ".temp", numChimeras);
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = outputFileName + toString(m->mothurGetpid(process)) + ".num.temp";
				m->openOutputFile(tempFile, out);
				out << num << endl;
				out << numChimeras << endl;
				out.close();
				
				exit(0);
            }else {
                m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine();
                for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                exit(0);
            }
		}
		
		//do my part
		num = driver(outputFileName, files[0], accnos, alns, numChimeras);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		for (int i = 0; i < processIDS.size(); i++) {
			ifstream in;
			string tempFile =  outputFileName + toString(processIDS[i]) + ".num.temp";
			m->openInputFile(tempFile, in);
			if (!in.eof()) { 
				int tempNum = 0; 
				in >> tempNum; m->gobble(in);
				num += tempNum; 
				in >> tempNum;
				numChimeras += tempNum;
			}
			in.close(); m->mothurRemove(tempFile);
		}
#else
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the preClusterData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		//divide file
		int count = 0;
		int spot = 0;
        files.resize(processors, "");
        
		for (int i = 0; i < processors; i++) {
            ofstream temp;
			files[i] = filename+toString(i)+".temp";
			m->openOutputFile(filename+toString(i)+".temp", temp);
		}
		
		ifstream in;
		m->openInputFile(filename, in);
		
		while(!in.eof()) {
			Sequence tempSeq(in); m->gobble(in); 
			
			if (tempSeq.getName() != "") {
                ofstream temp;
                m->openOutputFileAppend(files[spot], temp);
                tempSeq.printSequence(temp); temp.close();
				spot++; count++;
				if (spot == processors) { spot = 0; }
			}
		}
		in.close();
		
		//sanity check for number of processors
		if (count < processors) { processors = count; }
		
		vector<uchimeData*> pDataArray; 
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
		vector<string> dummy; //used so that we can use the same struct for MyUchimeSeqsThreadFunction and MyUchimeThreadFunction
		
		//Create processor worker threads.
		for( int i=1; i<processors; i++ ){
			// Allocate memory for thread data.
			string extension = toString(i) + ".temp";
			
			uchimeData* tempUchime = new uchimeData(outputFileName+extension, uchimeLocation, templatefile, files[i], "", "", "", accnos+extension, alns+extension, "", dummy, m, 0, 0,  i);
			tempUchime->setBooleans(dups, useAbskew, chimealns, useMinH, useMindiv, useXn, useDn, useXa, useChunks, useMinchunk, useIdsmoothwindow, useMinsmoothid, useMaxp, skipgaps, skipgaps2, useMinlen, useMaxlen, ucl, useQueryfract, hasCount);
			tempUchime->setVariables(abskew, minh, mindiv, xn, dn, xa, chunks, minchunk, idsmoothwindow, minsmoothid, maxp, minlen, maxlen, queryfract, strand);
			
			pDataArray.push_back(tempUchime);
			processIDS.push_back(i);
			
			//MySeqSumThreadFunction is in header. It must be global or static to work with the threads.
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i-1] = CreateThread(NULL, 0, MyUchimeSeqsThreadFunction, pDataArray[i-1], 0, &dwThreadIdArray[i-1]);   
		}
		
		
		//using the main process as a worker saves time and memory
		num = driver(outputFileName, files[0], accnos, alns, numChimeras);
		
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
			num += pDataArray[i]->count;
			numChimeras += pDataArray[i]->numChimeras;
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
#endif		
		
		//append output files
		for(int i=0;i<processIDS.size();i++){
			m->appendFiles((outputFileName + toString(processIDS[i]) + ".temp"), outputFileName);
			m->mothurRemove((outputFileName + toString(processIDS[i]) + ".temp"));
			
			m->appendFiles((accnos + toString(processIDS[i]) + ".temp"), accnos);
			m->mothurRemove((accnos + toString(processIDS[i]) + ".temp"));
			
			if (chimealns) {
				m->appendFiles((alns + toString(processIDS[i]) + ".temp"), alns);
				m->mothurRemove((alns + toString(processIDS[i]) + ".temp"));
			}
		}
		
		//get rid of the file pieces.
		for (int i = 0; i < files.size(); i++) { m->mothurRemove(files[i]); }
        
		return num;	
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraUchimeCommand", "createProcesses");
		exit(1);
	}
}
/**************************************************************************************************/

int ChimeraUchimeCommand::createProcessesGroups(string outputFName, string filename, string accnos, string alns, string newCountFile, vector<string> groups, string nameFile, string groupFile, string fastaFile) {
	try {
		
		processIDS.clear();
		int process = 1;
		int num = 0;
        
        CountTable newCount;
        if (hasCount && dups) { newCount.readTable(nameFile, true, false); }
		
		//sanity check
		if (groups.size() < processors) { processors = groups.size(); }
		
		//divide the groups between the processors
		vector<linePair> lines;
        int remainingPairs = groups.size();
        int startIndex = 0;
        for (int remainingProcessors = processors; remainingProcessors > 0; remainingProcessors--) {
            int numPairs = remainingPairs; //case for last processor
            if (remainingProcessors != 1) { numPairs = ceil(remainingPairs / remainingProcessors); }
            lines.push_back(linePair(startIndex, (startIndex+numPairs))); //startIndex, endIndex
            startIndex = startIndex + numPairs;
            remainingPairs = remainingPairs - numPairs;
        }

#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)		
				
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driverGroups(outputFName + toString(m->mothurGetpid(process)) + ".temp", filename + toString(m->mothurGetpid(process)) + ".temp", accnos + toString(m->mothurGetpid(process)) + ".temp", alns + toString(m->mothurGetpid(process)) + ".temp", accnos + ".byCount." + toString(m->mothurGetpid(process)) + ".temp", lines[process].start, lines[process].end, groups);
				
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
            
		//do my part
		num = driverGroups(outputFName, filename, accnos, alns, accnos + ".byCount", lines[0].start, lines[0].end, groups);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
        
		for (int i = 0; i < processIDS.size(); i++) {
			ifstream in;
			string tempFile =  outputFName + toString(processIDS[i]) + ".num.temp";
			m->openInputFile(tempFile, in);
			if (!in.eof()) { int tempNum = 0; in >> tempNum; num += tempNum; }
			in.close(); m->mothurRemove(tempFile);
        }
        
#else
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the uchimeData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<uchimeData*> pDataArray; 
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
		
		//Create processor worker threads.
		for( int i=1; i<processors; i++ ){
			// Allocate memory for thread data.
			string extension = toString(i) + ".temp";
			
			uchimeData* tempUchime = new uchimeData(outputFName+extension, uchimeLocation, templatefile, filename+extension, fastaFile, nameFile, groupFile, accnos+extension, alns+extension, accnos+".byCount."+extension, groups, m, lines[i].start, lines[i].end,  i);
			tempUchime->setBooleans(dups, useAbskew, chimealns, useMinH, useMindiv, useXn, useDn, useXa, useChunks, useMinchunk, useIdsmoothwindow, useMinsmoothid, useMaxp, skipgaps, skipgaps2, useMinlen, useMaxlen, ucl, useQueryfract, hasCount);
			tempUchime->setVariables(abskew, minh, mindiv, xn, dn, xa, chunks, minchunk, idsmoothwindow, minsmoothid, maxp, minlen, maxlen, queryfract, strand);
			
			pDataArray.push_back(tempUchime);
			processIDS.push_back(i);
			
			//MyUchimeThreadFunction is in header. It must be global or static to work with the threads.
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i-1] = CreateThread(NULL, 0, MyUchimeThreadFunction, pDataArray[i-1], 0, &dwThreadIdArray[i-1]);   
		}
		
		
		//using the main process as a worker saves time and memory
		num = driverGroups(outputFName, filename, accnos, alns, accnos + ".byCount", lines[0].start, lines[0].end, groups);
		
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
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
			
			if (chimealns) {
				m->appendFiles((alns + toString(processIDS[i]) + ".temp"), alns);
				m->mothurRemove((alns + toString(processIDS[i]) + ".temp"));
			}
            
            if (hasCount && dups) {
                if (!m->isBlank(accnos + ".byCount." + toString(processIDS[i]) + ".temp")) {
                    ifstream in2;
                    m->openInputFile(accnos + ".byCount." + toString(processIDS[i]) + ".temp", in2);
                    
                    string name, group;
                    while (!in2.eof()) {
                        in2 >> name >> group; m->gobble(in2);
                        newCount.setAbund(name, group, 0);
                    }
                    in2.close();
                }
                m->mothurRemove(accnos + ".byCount." + toString(processIDS[i]) + ".temp");
            }

		}
        
        //print new *.pick.count_table
        if (hasCount && dups) {  newCount.printTable(newCountFile);   }
 		
		return num;	
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraUchimeCommand", "createProcessesGroups");
		exit(1);
	}
}
/**************************************************************************************************/

