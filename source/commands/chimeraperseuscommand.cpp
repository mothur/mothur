/*
 *  chimeraperseuscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 10/26/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "chimeraperseuscommand.h"
#include "deconvolutecommand.h"
#include "sequence.hpp"
#include "counttable.h"
#include "sequencecountparser.h"
//**********************************************************************************************************************
vector<string> ChimeraPerseusCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","chimera-accnos",false,true,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "NameCount", "NameCount", "none","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "NameCount", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
        CommandParameter pdups("dereplicate", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pdups);

		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		CommandParameter pcutoff("cutoff", "Number", "", "0.5", "", "", "","",false,false); parameters.push_back(pcutoff);
		CommandParameter palpha("alpha", "Number", "", "-5.54", "", "", "","",false,false); parameters.push_back(palpha);
		CommandParameter pbeta("beta", "Number", "", "0.33", "", "", "","",false,false); parameters.push_back(pbeta);
			
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraPerseusCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The chimera.perseus command reads a fastafile and namefile or countfile and outputs potentially chimeric sequences.\n";
		helpString += "The chimera.perseus command parameters are fasta, name, group, cutoff, processors, dereplicate, alpha and beta.\n";
		helpString += "The fasta parameter allows you to enter the fasta file containing your potentially chimeric sequences, and is required, unless you have a valid current fasta file. \n";
		helpString += "The name parameter allows you to provide a name file associated with your fasta file.\n";
        helpString += "The count parameter allows you to provide a count file associated with your fasta file. A count or name file is required. When you use a count file with group info and dereplicate=T, mothur will create a *.pick.count_table file containing seqeunces after chimeras are removed.\n";
		helpString += "You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amazon.fasta \n";
		helpString += "The group parameter allows you to provide a group file.  When checking sequences, only sequences from the same group as the query sequence will be used as the reference. \n";
		helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n";
        helpString += "If the dereplicate parameter is false, then if one group finds the seqeunce to be chimeric, then all groups find it to be chimeric, default=f.\n";
		helpString += "The alpha parameter ....  The default is -5.54. \n";
		helpString += "The beta parameter ....  The default is 0.33. \n";
		helpString += "The cutoff parameter ....  The default is 0.50. \n";
		helpString += "The chimera.perseus command should be in the following format: \n";
		helpString += "chimera.perseus(fasta=yourFastaFile, name=yourNameFile) \n";
		helpString += "Example: chimera.perseus(fasta=AD.align, name=AD.names) \n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraPerseusCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "chimera") {  pattern = "[filename],perseus.chimeras"; } 
        else if (type == "accnos") {  pattern = "[filename],perseus.accnos"; }
        else if (type == "count") {  pattern = "[filename],perseus.pick.count_table"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraPerseusCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ChimeraPerseusCommand::ChimeraPerseusCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
		vector<string> tempOutNames;
		outputTypes["chimera"] = tempOutNames;
		outputTypes["accnos"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "ChimeraPerseusCommand");
		exit(1);
	}
}
//***************************************************************************************************************
ChimeraPerseusCommand::ChimeraPerseusCommand(string option)  {
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
			
			ValidParameters validParameter("chimera.perseus");
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			vector<string> tempOutNames;
			outputTypes["chimera"] = tempOutNames;
			outputTypes["accnos"] = tempOutNames;
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
            
            if (countfileNames.size() != 0) { hasCount = true; }
            
			//make sure there is at least one valid file left
            if (hasName && hasCount) { m->mothurOut("[ERROR]: You must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
            
            if (!hasName && !hasCount) { 
                //if there is a current name file, use it, else look for current count file
				string filename = m->getNameFile(); 
				if (filename != "") { hasName = true; nameFileNames.push_back(filename); m->mothurOut("Using " + filename + " as input file for the name parameter."); m->mothurOutEndLine(); }
				else { 
                    filename = m->getCountTableFile();
                    if (filename != "") { hasCount = true; countfileNames.push_back(filename); m->mothurOut("Using " + filename + " as input file for the count parameter."); m->mothurOutEndLine(); }
                    else { m->mothurOut("[ERROR]: You must provide a count or name file."); m->mothurOutEndLine(); abort = true;  }
                }
            }
            if (!hasName && hasCount) { nameFileNames = countfileNames; }
            
			if (nameFileNames.size() != fastaFileNames.size()) { m->mothurOut("[ERROR]: The number of name or count files does not match the number of fastafiles, please correct."); m->mothurOutEndLine(); abort=true; }
			
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
			
			string temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
			
			temp = validParameter.validFile(parameters, "cutoff", false);	if (temp == "not found"){	temp = "0.50";	}
			m->mothurConvert(temp, cutoff);
			
			temp = validParameter.validFile(parameters, "alpha", false);	if (temp == "not found"){	temp = "-5.54";	}
			m->mothurConvert(temp, alpha);
			
			temp = validParameter.validFile(parameters, "cutoff", false);	if (temp == "not found"){	temp = "0.33";	}
			m->mothurConvert(temp, beta);
            
			temp = validParameter.validFile(parameters, "dereplicate", false);	
			if (temp == "not found") { temp = "false";			}
			dups = m->isTrue(temp);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "ChimeraPerseusCommand");
		exit(1);
	}
}
//***************************************************************************************************************

int ChimeraPerseusCommand::execute(){
	try{
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
				
		//process each file
		for (int s = 0; s < fastaFileNames.size(); s++) {
			
			m->mothurOut("Checking sequences from " + fastaFileNames[s] + " ..." ); m->mothurOutEndLine();
			
			int start = time(NULL);	
			if (outputDir == "") { outputDir = m->hasPath(fastaFileNames[s]);  }//if user entered a file with a path then preserve it	
			map<string, string> variables;
			variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s]));
			string outputFileName = getOutputFileName("chimera", variables);
			string accnosFileName = getOutputFileName("accnos", variables);
            string newCountFile = "";

			//string newFasta = m->getRootName(fastaFileNames[s]) + "temp";
			
			//you provided a groupfile
			string groupFile = "";
			if (groupFileNames.size() != 0) { groupFile = groupFileNames[s]; }
			
			string nameFile = "";
			if (nameFileNames.size() != 0) { //you provided a namefile and we don't need to create one
				nameFile = nameFileNames[s];
			}else { nameFile = getNamesFile(fastaFileNames[s]); }
			
			if (m->control_pressed) {  for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	} return 0;	}				
			
			int numSeqs = 0;
			int numChimeras = 0;
            
            if (hasCount) {
                CountTable* ct = new CountTable();
                ct->readTable(nameFile, true, false);
                
                if (ct->hasGroupInfo()) {
                    cparser = new SequenceCountParser(fastaFileNames[s], *ct);
                    variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(nameFile));
                    newCountFile = getOutputFileName("count", variables);
                    
                    vector<string> groups = cparser->getNamesOfGroups();
                    
                    if (m->control_pressed) { delete ct; delete cparser; for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0; }
                    
                    //clears files
                    ofstream out, out1, out2;
                    m->openOutputFile(outputFileName, out); out.close(); 
                    m->openOutputFile(accnosFileName, out1); out1.close();
                    
                    if(processors == 1)	{	numSeqs = driverGroups(outputFileName, accnosFileName, newCountFile, 0, groups.size(), groups);
                        if (dups) {
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

                    }
                    else				{	numSeqs = createProcessesGroups(outputFileName, accnosFileName, newCountFile, groups, groupFile, fastaFileNames[s], nameFile);			}
                    
                    if (m->control_pressed) {  delete ct; delete cparser; for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0;	}				
                    map<string, string> uniqueNames = cparser->getAllSeqsMap();
                    if (!dups) { 
                        numChimeras = deconvoluteResults(uniqueNames, outputFileName, accnosFileName);
                    }else {
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
                    delete cparser;

                    m->mothurOut("The number of sequences checked may be larger than the number of unique sequences because some sequences are found in several samples."); m->mothurOutEndLine(); 
                    
                    if (m->control_pressed) {  delete ct; for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0;  }	
                    
                }else {
                    if (processors != 1) { m->mothurOut("Your count file does not contain group information, mothur can only use 1 processor, continuing."); m->mothurOutEndLine(); processors = 1; }
                    
                    //read sequences and store sorted by frequency
                    vector<seqData> sequences = readFiles(fastaFileNames[s], ct);
                    
                    if (m->control_pressed) { delete ct; for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	} return 0; }
                    
                    numSeqs = driver(outputFileName, sequences, accnosFileName, numChimeras);   
                }
                delete ct;
            }else {
                if (groupFile != "") {
                    //Parse sequences by group
                    parser = new SequenceParser(groupFile, fastaFileNames[s], nameFile);
                    vector<string> groups = parser->getNamesOfGroups();
                    
                    if (m->control_pressed) { delete parser; for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0; }
                    
                    //clears files
                    ofstream out, out1, out2;
                    m->openOutputFile(outputFileName, out); out.close(); 
                    m->openOutputFile(accnosFileName, out1); out1.close();
                    
                    if(processors == 1)	{	numSeqs = driverGroups(outputFileName, accnosFileName, "", 0, groups.size(), groups);	}
                    else				{	numSeqs = createProcessesGroups(outputFileName, accnosFileName, "", groups, groupFile, fastaFileNames[s], nameFile);			}
                    
                    if (m->control_pressed) {  delete parser; for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0;	}				
                    map<string, string> uniqueNames = parser->getAllSeqsMap();
                    if (!dups) { 
                        numChimeras = deconvoluteResults(uniqueNames, outputFileName, accnosFileName);
                    }
                    delete parser;
                    
                    m->mothurOut("The number of sequences checked may be larger than the number of unique sequences because some sequences are found in several samples."); m->mothurOutEndLine(); 
                    
                    if (m->control_pressed) {  for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0;  }		
                }else{
                    if (processors != 1) { m->mothurOut("Without a groupfile, mothur can only use 1 processor, continuing."); m->mothurOutEndLine(); processors = 1; }
                    
                    //read sequences and store sorted by frequency
                    vector<seqData> sequences = readFiles(fastaFileNames[s], nameFile);
                    
                    if (m->control_pressed) { for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	} return 0; }
                    
                    numSeqs = driver(outputFileName, sequences, accnosFileName, numChimeras); 
                }
			}
            
			if (m->control_pressed) { for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	} return 0; }
			
			m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences. " + toString(numChimeras) + " chimeras were found.");	m->mothurOutEndLine();
			outputNames.push_back(outputFileName); outputTypes["chimera"].push_back(outputFileName);
			outputNames.push_back(accnosFileName); outputTypes["accnos"].push_back(accnosFileName);
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
		m->errorOut(e, "ChimeraPerseusCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraPerseusCommand::getNamesFile(string& inputFile){
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
		m->errorOut(e, "ChimeraPerseusCommand", "getNamesFile");
		exit(1);
	}
}
//**********************************************************************************************************************
int ChimeraPerseusCommand::driverGroups(string outputFName, string accnos, string countlist, int start, int end, vector<string> groups){
	try {
		
		int totalSeqs = 0;
		int numChimeras = 0;
        
        ofstream outCountList;
        if (hasCount && dups) { m->openOutputFile(countlist, outCountList); }
		
		for (int i = start; i < end; i++) {
			
			m->mothurOutEndLine(); m->mothurOut("Checking sequences from group " + groups[i] + "...");	m->mothurOutEndLine();					
			
			int start = time(NULL);	 if (m->control_pressed) {  return 0; }
			
			vector<seqData> sequences = loadSequences(groups[i]);
			
			if (m->control_pressed) { return 0; }
			
			int numSeqs = driver((outputFName + groups[i]), sequences, (accnos+groups[i]), numChimeras);
			totalSeqs += numSeqs;
			
			if (m->control_pressed) { return 0; }
            
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
                        map<string, string> thisnamemap = parser->getNameMap(groups[i]);
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
			
			m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences from group " + groups[i] + ".");	m->mothurOutEndLine();					
		}	
		
        if (hasCount && dups) { outCountList.close(); }
        
		return totalSeqs;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "driverGroups");
		exit(1);
	}
}	
//**********************************************************************************************************************
vector<seqData> ChimeraPerseusCommand::loadSequences(string group){
	try {
        bool error = false;
		alignLength = 0;
        vector<seqData> sequences;
        if (hasCount) {
            vector<Sequence> thisGroupsSeqs = cparser->getSeqs(group);
            map<string, int> counts = cparser->getCountTable(group);
            map<string, int>::iterator it;
            
            for (int i = 0; i < thisGroupsSeqs.size(); i++) {
                
                if (m->control_pressed) {  return sequences; }
                
                it = counts.find(thisGroupsSeqs[i].getName());
                if (it == counts.end()) { error = true; m->mothurOut("[ERROR]: " + thisGroupsSeqs[i].getName() + " is in your fasta file and not in your count file, please correct."); m->mothurOutEndLine(); }
                else {
                    thisGroupsSeqs[i].setAligned(removeNs(thisGroupsSeqs[i].getUnaligned()));
                    sequences.push_back(seqData(thisGroupsSeqs[i].getName(), thisGroupsSeqs[i].getUnaligned(), it->second));
                    if (thisGroupsSeqs[i].getUnaligned().length() > alignLength) { alignLength = thisGroupsSeqs[i].getUnaligned().length(); }
                }
            }
        }else{
            vector<Sequence> thisGroupsSeqs = parser->getSeqs(group);
            map<string, string> nameMap = parser->getNameMap(group);
            map<string, string>::iterator it;
           
            for (int i = 0; i < thisGroupsSeqs.size(); i++) {
                
                if (m->control_pressed) {  return sequences; }
                
                it = nameMap.find(thisGroupsSeqs[i].getName());
                if (it == nameMap.end()) { error = true; m->mothurOut("[ERROR]: " + thisGroupsSeqs[i].getName() + " is in your fasta file and not in your namefile, please correct."); m->mothurOutEndLine(); }
                else {
                    int num = m->getNumNames(it->second);
                    thisGroupsSeqs[i].setAligned(removeNs(thisGroupsSeqs[i].getUnaligned()));
                    sequences.push_back(seqData(thisGroupsSeqs[i].getName(), thisGroupsSeqs[i].getUnaligned(), num));
                    if (thisGroupsSeqs[i].getUnaligned().length() > alignLength) { alignLength = thisGroupsSeqs[i].getUnaligned().length(); }
                }
            }
            
		}
		
        if (error) { m->control_pressed = true; }
		//sort by frequency
		sort(sequences.rbegin(), sequences.rend());
		
		return sequences;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "loadSequences");
		exit(1);
	}
}

//**********************************************************************************************************************
vector<seqData> ChimeraPerseusCommand::readFiles(string inputFile, string name){
	try {
		map<string, int>::iterator it;
		map<string, int> nameMap = m->readNames(name);
		
		//read fasta file and create sequenceData structure - checking for file mismatches
		vector<seqData> sequences;
		bool error = false;
		ifstream in;
		m->openInputFile(inputFile, in);
		alignLength = 0;
        
		while (!in.eof()) {
			
			if (m->control_pressed) { in.close(); return sequences; }
			
			Sequence temp(in); m->gobble(in);
			
			it = nameMap.find(temp.getName());
			if (it == nameMap.end()) { error = true; m->mothurOut("[ERROR]: " + temp.getName() + " is in your fasta file and not in your namefile, please correct."); m->mothurOutEndLine(); }
			else {
                temp.setAligned(removeNs(temp.getUnaligned()));
				sequences.push_back(seqData(temp.getName(), temp.getUnaligned(), it->second));
                if (temp.getUnaligned().length() > alignLength) { alignLength = temp.getUnaligned().length(); }
			}
		}
		in.close();
		
		if (error) { m->control_pressed = true; }
		
		//sort by frequency
		sort(sequences.rbegin(), sequences.rend());
		
		return sequences;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "readFiles");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraPerseusCommand::removeNs(string seq){
	try {
        string newSeq = "";
        for (int i = 0; i < seq.length(); i++) {
            if (seq[i] != 'N') {  newSeq += seq[i]; }
        }
        return newSeq;
    }
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "removeNs");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<seqData> ChimeraPerseusCommand::readFiles(string inputFile, CountTable* ct){
	try {		
		//read fasta file and create sequenceData structure - checking for file mismatches
		vector<seqData> sequences;
		ifstream in;
		m->openInputFile(inputFile, in);
		alignLength = 0;
        
		while (!in.eof()) {
            Sequence temp(in); m->gobble(in);
			
			int count = ct->getNumSeqs(temp.getName());
			if (m->control_pressed) { break; }
			else {
                temp.setAligned(removeNs(temp.getUnaligned()));
				sequences.push_back(seqData(temp.getName(), temp.getUnaligned(), count));
                if (temp.getUnaligned().length() > alignLength) { alignLength = temp.getUnaligned().length(); }
			}
		}
		in.close();
		
		//sort by frequency
		sort(sequences.rbegin(), sequences.rend());
		
		return sequences;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "getNamesFile");
		exit(1);
	}
}
//**********************************************************************************************************************
int ChimeraPerseusCommand::driver(string chimeraFileName, vector<seqData>& sequences, string accnosFileName, int& numChimeras){
	try {
		
		vector<vector<double> > correctModel(4);	//could be an option in the future to input own model matrix
		for(int i=0;i<4;i++){	correctModel[i].resize(4);	}
		
		correctModel[0][0] = 0.000000;	//AA
		correctModel[1][0] = 11.619259;	//CA
		correctModel[2][0] = 11.694004;	//TA
		correctModel[3][0] = 7.748623;	//GA
		
		correctModel[1][1] = 0.000000;	//CC
		correctModel[2][1] = 7.619657;	//TC
		correctModel[3][1] = 12.852562;	//GC
		
		correctModel[2][2] = 0.000000;	//TT
		correctModel[3][2] = 10.964048;	//TG
		
		correctModel[3][3] = 0.000000;	//GG
		
		for(int i=0;i<4;i++){
			for(int j=0;j<i;j++){
				correctModel[j][i] = correctModel[i][j];
			}
		}
		
		int numSeqs = sequences.size();
		//int alignLength = sequences[0].sequence.size();
		
		ofstream chimeraFile;
		ofstream accnosFile;
		m->openOutputFile(chimeraFileName, chimeraFile); 
		m->openOutputFile(accnosFileName, accnosFile); 
		
		Perseus myPerseus;
		vector<vector<double> > binMatrix = myPerseus.binomial(alignLength);
		
		chimeraFile << "SequenceIndex\tName\tDiffsToBestMatch\tBestMatchIndex\tBestMatchName\tDiffstToChimera\tIndexofLeftParent\tIndexOfRightParent\tNameOfLeftParent\tNameOfRightParent\tDistanceToBestMatch\tcIndex\t(cIndex - singleDist)\tloonIndex\tMismatchesToChimera\tMismatchToTrimera\tChimeraBreakPoint\tLogisticProbability\tTypeOfSequence\n";
		
		vector<bool> chimeras(numSeqs, 0);
		
		for(int i=0;i<numSeqs;i++){	
			if (m->control_pressed) { chimeraFile.close(); m->mothurRemove(chimeraFileName); accnosFile.close(); m->mothurRemove(accnosFileName); return 0; }
    
			vector<bool> restricted = chimeras;
			
			vector<vector<int> > leftDiffs(numSeqs);
			vector<vector<int> > leftMaps(numSeqs);
			vector<vector<int> > rightDiffs(numSeqs);
			vector<vector<int> > rightMaps(numSeqs);
			
			vector<int> singleLeft, bestLeft;
			vector<int> singleRight, bestRight;
			
			int bestSingleIndex, bestSingleDiff;
			vector<pwAlign> alignments(numSeqs);
			
			int comparisons = myPerseus.getAlignments(i, sequences, alignments, leftDiffs, leftMaps, rightDiffs, rightMaps, bestSingleIndex, bestSingleDiff, restricted);
			if (m->control_pressed) { chimeraFile.close(); m->mothurRemove(chimeraFileName); accnosFile.close(); m->mothurRemove(accnosFileName); return 0; }

			int minMismatchToChimera, leftParentBi, rightParentBi, breakPointBi;
			
			string dummyA, dummyB;
			
            if (sequences[i].sequence.size() < 3) { 
                chimeraFile << i << '\t' << sequences[i].seqName << "\t0\t0\tNull\t0\t0\t0\tNull\tNull\t0.0\t0.0\t0.0\t0\t0\t0\t0.0\t0.0\tgood" << endl;
            }else if(comparisons >= 2){	
				minMismatchToChimera = myPerseus.getChimera(sequences, leftDiffs, rightDiffs, leftParentBi, rightParentBi, breakPointBi, singleLeft, bestLeft, singleRight, bestRight, restricted);
				if (m->control_pressed) { chimeraFile.close(); m->mothurRemove(chimeraFileName); accnosFile.close(); m->mothurRemove(accnosFileName); return 0; }

				int minMismatchToTrimera = numeric_limits<int>::max();
				int leftParentTri, middleParentTri, rightParentTri, breakPointTriA, breakPointTriB;
				
				if(minMismatchToChimera >= 3 && comparisons >= 3){
					minMismatchToTrimera = myPerseus.getTrimera(sequences, leftDiffs, leftParentTri, middleParentTri, rightParentTri, breakPointTriA, breakPointTriB, singleLeft, bestLeft, singleRight, bestRight, restricted);
					if (m->control_pressed) { chimeraFile.close(); m->mothurRemove(chimeraFileName); accnosFile.close(); m->mothurRemove(accnosFileName); return 0; }
				}
				
				double singleDist = myPerseus.modeledPairwiseAlignSeqs(sequences[i].sequence, sequences[bestSingleIndex].sequence, dummyA, dummyB, correctModel);
				
				if (m->control_pressed) { chimeraFile.close(); m->mothurRemove(chimeraFileName); accnosFile.close(); m->mothurRemove(accnosFileName); return 0; }

				string type;
				string chimeraRefSeq;
				
				if(minMismatchToChimera - minMismatchToTrimera >= 3){
					type = "trimera";
					chimeraRefSeq = myPerseus.stitchTrimera(alignments, leftParentTri, middleParentTri, rightParentTri, breakPointTriA, breakPointTriB, leftMaps, rightMaps);
				}
				else{
					type = "chimera";
					chimeraRefSeq = myPerseus.stitchBimera(alignments, leftParentBi, rightParentBi, breakPointBi, leftMaps, rightMaps);
				}

                if (m->control_pressed) { chimeraFile.close(); m->mothurRemove(chimeraFileName); accnosFile.close(); m->mothurRemove(accnosFileName); return 0; }
				
				double chimeraDist = myPerseus.modeledPairwiseAlignSeqs(sequences[i].sequence, chimeraRefSeq, dummyA, dummyB, correctModel);
				
				if (m->control_pressed) { chimeraFile.close(); m->mothurRemove(chimeraFileName); accnosFile.close(); m->mothurRemove(accnosFileName); return 0; }

				double cIndex = chimeraDist;//modeledPairwiseAlignSeqs(sequences[i].sequence, chimeraRefSeq);
				double loonIndex = myPerseus.calcLoonIndex(sequences[i].sequence, sequences[leftParentBi].sequence, sequences[rightParentBi].sequence, breakPointBi, binMatrix);		
				
				if (m->control_pressed) { chimeraFile.close(); m->mothurRemove(chimeraFileName); accnosFile.close(); m->mothurRemove(accnosFileName); return 0; }

				chimeraFile << i << '\t' << sequences[i].seqName << '\t' << bestSingleDiff << '\t' << bestSingleIndex << '\t' << sequences[bestSingleIndex].seqName << '\t';
				chimeraFile << minMismatchToChimera << '\t' << leftParentBi << '\t' << rightParentBi << '\t' << sequences[leftParentBi].seqName << '\t' << sequences[rightParentBi].seqName << '\t';
				chimeraFile << singleDist << '\t' << cIndex << '\t' << (cIndex - singleDist) << '\t' << loonIndex << '\t';
				chimeraFile << minMismatchToChimera << '\t' << minMismatchToTrimera << '\t' << breakPointBi << '\t';
				
				double probability = myPerseus.classifyChimera(singleDist, cIndex, loonIndex, alpha, beta);
				
				chimeraFile << probability << '\t';
				
				if(probability > cutoff){ 
					chimeraFile << type << endl;
					accnosFile << sequences[i].seqName << endl;
					chimeras[i] = 1;
					numChimeras++;
				}
				else{
					chimeraFile << "good" << endl;
				}
				
			}
			else{
				chimeraFile << i << '\t' << sequences[i].seqName << "\t0\t0\tNull\t0\t0\t0\tNull\tNull\t0.0\t0.0\t0.0\t0\t0\t0\t0.0\t0.0\tgood" << endl;
			}
	
			//report progress
			if((i+1) % 100 == 0){ 	m->mothurOutJustToScreen("Processing sequence: " + toString(i+1) + "\n");		}
		}
		
		if((numSeqs) % 100 != 0){ 	m->mothurOutJustToScreen("Processing sequence: " + toString(numSeqs) + "\n");		}
		
		chimeraFile.close();
		accnosFile.close();
		
		return numSeqs;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "driver");
		exit(1);
	}
}
/**************************************************************************************************/
int ChimeraPerseusCommand::createProcessesGroups(string outputFName, string accnos, string newCountFile, vector<string> groups, string group, string fasta, string name) {
	try {
		
		vector<int> processIDS;
		int process = 1;
		int num = 0;
		
        CountTable newCount;
        if (hasCount && dups) { newCount.readTable(name, true, false); }
        
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
				num = driverGroups(outputFName + toString(m->mothurGetpid(process)) + ".temp", accnos + toString(m->mothurGetpid(process)) + ".temp", accnos + ".byCount." + toString(m->mothurGetpid(process)) + ".temp", lines[process].start, lines[process].end, groups);
				
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
		num = driverGroups(outputFName, accnos, accnos + ".byCount", lines[0].start, lines[0].end, groups);
		
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
		//Windows version shared memory, so be careful when passing variables through the preClusterData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<perseusData*> pDataArray; 
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
		
		//Create processor worker threads.
		for( int i=1; i<processors; i++ ){
			// Allocate memory for thread data.
			string extension = toString(i) + ".temp";
			
			perseusData* tempPerseus = new perseusData(dups, hasName, hasCount, alpha, beta, cutoff, outputFName+extension, fasta, name, group, accnos+extension,  accnos+".byCount."+extension, groups, m, lines[i].start, lines[i].end, i);
			
			pDataArray.push_back(tempPerseus);
			processIDS.push_back(i);
			
			//MyPerseusThreadFunction is in header. It must be global or static to work with the threads.
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i-1] = CreateThread(NULL, 0, MyPerseusThreadFunction, pDataArray[i-1], 0, &dwThreadIdArray[i-1]);   
		}
		
		
		//using the main process as a worker saves time and memory
		num = driverGroups(outputFName, accnos, accnos + ".byCount", lines[0].start, lines[0].end, groups);
		
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
		m->errorOut(e, "ChimeraPerseusCommand", "createProcessesGroups");
		exit(1);
	}
}
//**********************************************************************************************************************
int ChimeraPerseusCommand::deconvoluteResults(map<string, string>& uniqueNames, string outputFileName, string accnosFileName){
	try {
		map<string, string>::iterator itUnique;
		int total = 0;
		
		//edit accnos file
		ifstream in2; 
		m->openInputFile(accnosFileName, in2);
		
		ofstream out2;
		m->openOutputFile(accnosFileName+".temp", out2);
		
		string name;
		set<string> namesInFile; //this is so if a sequence is found to be chimera in several samples we dont write it to the results file more than once
		set<string>::iterator itNames;
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
		
		int DiffsToBestMatch, BestMatchIndex, DiffstToChimera, IndexofLeftParent, IndexOfRightParent;
		float temp1,temp2, temp3, temp4, temp5, temp6, temp7, temp8;
		string index, BestMatchName, parent1, parent2, flag;
		name = "";
		namesInFile.clear();	
		//assumptions - in file each read will always look like 
		/*										
		 SequenceIndex	Name	DiffsToBestMatch	BestMatchIndex	BestMatchName	DiffstToChimera	IndexofLeftParent	IndexOfRightParent	NameOfLeftParent	NameOfRightParent	DistanceToBestMatch	cIndex	(cIndex - singleDist)	loonIndex	MismatchesToChimera	MismatchToTrimera	ChimeraBreakPoint	LogisticProbability	TypeOfSequence
		 0	F01QG4L02JVBQY	0	0	Null	0	0	0	Null	Null	0.0	0.0	0.0	0.0	0	0	0	0.0	0.0	good
		 1	F01QG4L02ICTC6	0	0	Null	0	0	0	Null	Null	0.0	0.0	0.0	0.0	0	0	0	0.0	0.0	good
		 2	F01QG4L02JZOEC	48	0	F01QG4L02JVBQY	47	0	0	F01QG4L02JVBQY	F01QG4L02JVBQY	2.0449	2.03545	-0.00944493	0	47	2147483647	138	0	good
		 3	F01QG4L02G7JEC	42	0	F01QG4L02JVBQY	40	1	0	F01QG4L02ICTC6	F01QG4L02JVBQY	1.87477	1.81113	-0.0636404	5.80145	40	2147483647	25	0	good
		 */
		
		//get and print headers
		BestMatchName = m->getline(in); m->gobble(in);
		out << BestMatchName << endl;
		
		while (!in.eof()) {
			
			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove((outputFileName+".temp")); return 0; }
			
			bool print = false;
			in >> index;	m->gobble(in);
			
			if (index != "SequenceIndex") { //if you are not a header line, there will be a header line for each group if group file is given
				in >> name;		m->gobble(in);
				in >> DiffsToBestMatch; m->gobble(in);
				in >> BestMatchIndex; m->gobble(in);
				in >> BestMatchName; m->gobble(in);
				in >> DiffstToChimera; m->gobble(in);
				in >> IndexofLeftParent; m->gobble(in);
				in >> IndexOfRightParent; m->gobble(in);
				in >> parent1;	m->gobble(in);
				in >> parent2;	m->gobble(in);
				in >> temp1 >> temp2 >> temp3 >> temp4 >> temp5 >> temp6 >> temp7 >> temp8 >> flag; m->gobble(in);
				
				//find unique name
				itUnique = uniqueNames.find(name);
				
				if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find "+ name + "."); m->mothurOutEndLine(); m->control_pressed = true; }
				else {
					name = itUnique->second;
					//is this name already in the file
					itNames = namesInFile.find((name));
					
					if (itNames == namesInFile.end()) { //no not in file
						if (flag == "good") { //are you really a no??
							//is this sequence really not chimeric??
							itChimeras = chimerasInFile.find(name);
							
							//then you really are a no so print, otherwise skip
							if (itChimeras == chimerasInFile.end()) { print = true; }
						}else{ print = true; }
					}
				}
				
				if (print) {
					out << index << '\t' << name  << '\t' << DiffsToBestMatch << '\t' << BestMatchIndex << '\t';
					namesInFile.insert(name);
					
					if (BestMatchName != "Null") {
						itUnique = uniqueNames.find(BestMatchName);
						if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find BestMatchName "+ BestMatchName + "."); m->mothurOutEndLine(); m->control_pressed = true; }
						else {	out << itUnique->second << '\t';	}					
					}else { out << "Null" << '\t'; }
					
					out << DiffstToChimera << '\t' << IndexofLeftParent << '\t' << IndexOfRightParent << '\t';
					
					if (parent1 != "Null") {
						itUnique = uniqueNames.find(parent1);
						if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find parent1 "+ parent1 + "."); m->mothurOutEndLine(); m->control_pressed = true; }
						else {	out << itUnique->second << '\t';	}
					}else { out << "Null" << '\t'; }
					
					if (parent1 != "Null") {
						itUnique = uniqueNames.find(parent2);
						if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find parent2 "+ parent2 + "."); m->mothurOutEndLine(); m->control_pressed = true; }
						else {	out << itUnique->second << '\t';	}
					}else { out << "Null" << '\t'; }
					
					out << temp1 << '\t' << temp2 << '\t' << temp3 << '\t' << temp4 << '\t' << temp5 << '\t' << temp6 << '\t' << temp7 << '\t' << temp8 << '\t' << flag << endl;	
				}
			}else { index = m->getline(in); m->gobble(in); }
		}
		in.close();
		out.close();
		
		m->mothurRemove(outputFileName);
		rename((outputFileName+".temp").c_str(), outputFileName.c_str());
		
		return total;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "deconvoluteResults");
		exit(1);
	}
}	
//**********************************************************************************************************************


