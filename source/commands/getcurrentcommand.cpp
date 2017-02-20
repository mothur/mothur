/*
 *  getcurrentcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 3/16/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "getcurrentcommand.h"

//**********************************************************************************************************************
vector<string> GetCurrentCommand::setParameters(){	
	try {
		CommandParameter pclear("clear", "String", "", "", "", "", "","",false,false); parameters.push_back(pclear);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetCurrentCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetCurrentCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The get.current command outputs the current files saved by mothur.\n";
		helpString += "The get.current command has one parameter: clear.\n";
		helpString += "The clear parameter is used to indicate which file types you would like to clear values for, multiple types can be separated by dashes.\n";
		helpString += "The get.current command should be in the following format: \n";
		helpString += "get.current() or get.current(clear=fasta-name-accnos)\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "GetCurrentCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
GetCurrentCommand::GetCurrentCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
        vector<string> tempOutNames;
        outputTypes["summary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "GetCurrentCommand", "GetCurrentCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetCurrentCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "summary") {  pattern = "[filename],current_files.summary"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "GetCurrentCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
GetCurrentCommand::GetCurrentCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
            vector<string> tempOutNames;
            outputTypes["summary"] = tempOutNames;
            
			clearTypes = validParameter.validFile(parameters, "clear", false);			
			if (clearTypes == "not found") { clearTypes = ""; }
			else { m->splitAtDash(clearTypes, types);	}
            
            //if the user changes the output directory command factory will send this info to us in the output parameter
            outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){  outputDir = ""; }
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "GetCurrentCommand", "GetCurrentCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetCurrentCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
		cFactory = CommandFactory::getInstance();
        
		//user wants to clear a type
		if (types.size() != 0) {
			for (int i = 0; i < types.size(); i++) {
				
				if (m->control_pressed) { break; }
				
				//look for file types
				if (types[i] == "fasta") {
					m->setFastaFile("");
				}else if (types[i] == "qfile") {
					m->setQualFile("");
				}else if (types[i] == "phylip") {
					m->setPhylipFile("");
				}else if (types[i] == "column") {
					m->setColumnFile("");
				}else if (types[i] == "list") {
					m->setListFile("");
				}else if (types[i] == "rabund") {
					m->setRabundFile("");
				}else if (types[i] == "sabund") {
					m->setSabundFile("");
				}else if (types[i] == "name") {
					m->setNameFile("");
				}else if (types[i] == "group") {
					m->setGroupFile("");
				}else if (types[i] == "order") {
					m->setOrderFile("");
				}else if (types[i] == "ordergroup") {
					m->setOrderGroupFile("");
				}else if (types[i] == "tree") {
					m->setTreeFile("");
				}else if (types[i] == "shared") {
					m->setSharedFile("");
				}else if (types[i] == "relabund") {
					m->setRelAbundFile("");
				}else if (types[i] == "design") {
					m->setDesignFile("");
				}else if (types[i] == "sff") {
					m->setSFFFile("");
				}else if (types[i] == "oligos") {
					m->setOligosFile("");
				}else if (types[i] == "accnos") {
					m->setAccnosFile("");
				}else if (types[i] == "taxonomy") {
					m->setTaxonomyFile("");
                }else if (types[i] == "constaxonomy") {
                    m->setConsTaxonomyFile("");
                }else if (types[i] == "contigsreport") {
                    m->setContigsReportFile("");
				}else if (types[i] == "flow") {
					m->setFlowFile("");
                }else if (types[i] == "biom") {
					m->setBiomFile("");
                }else if (types[i] == "count") {
					m->setCountTableFile("");
                }else if (types[i] == "summary") {
					m->setSummaryFile("");
                }else if (types[i] == "file") {
                    m->setFileFile("");
				}else if (types[i] == "processors") {
					m->setProcessors("1");
				}else if (types[i] == "all") {
					m->clearCurrentFiles();
				}else {
					m->mothurOut("[ERROR]: mothur does not save a current file for " + types[i]); m->mothurOutEndLine();
				}
			}
		}
		
        unsigned long long ramUsed, total;
        ramUsed = m->getRAMUsed(); total = m->getTotalRAM();
        m->mothurOut("\nCurrent RAM usage: " + toString(ramUsed/(double)GIG) + " Gigabytes. Total Ram: " + toString(total/(double)GIG) + " Gigabytes.\n");
        
		if (m->hasCurrentFiles()) {
            map<string, string> variables;
            variables["[filename]"] = m->getFullPathName(outputDir);
            string filename = getOutputFileName("summary", variables);
            
			m->mothurOutEndLine(); m->mothurOut("Current files saved by mothur:"); m->mothurOutEndLine();
			m->printCurrentFiles(filename);
            outputNames.push_back(filename); outputTypes["summary"].push_back(filename);
		}
        
        string inputDir = cFactory->getInputDir();
        if (inputDir != "") {
            m->mothurOutEndLine(); m->mothurOut("Current input directory saved by mothur: " + inputDir); m->mothurOutEndLine();
        }
        
        string outputDir = cFactory->getOutputDir();
        if (outputDir != "") {
            m->mothurOutEndLine(); m->mothurOut("Current output directory saved by mothur: " + outputDir); m->mothurOutEndLine();
        }
        string defaultPath = m->getDefaultPath();
        if (defaultPath != "") {
            m->mothurOutEndLine(); m->mothurOut("Current default directory saved by mothur: " + defaultPath); m->mothurOutEndLine();
        }
        
        
        string temp = "./";
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
            temp = ".\\";
#endif
        temp = m->getFullPathName(temp);
        m->mothurOutEndLine(); m->mothurOut("Current working directory: " + temp); m->mothurOutEndLine();
        
        if (m->hasCurrentFiles()) {
            m->mothurOutEndLine();
            m->mothurOut("Output File Names: "); m->mothurOutEndLine();
            for (int i = 0; i < outputNames.size(); i++) { m->mothurOut(outputNames[i]); m->mothurOutEndLine(); }
            m->mothurOutEndLine();
        }
        
        return 0;
        
	}
	
	catch(exception& e) {
		m->errorOut(e, "GetCurrentCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************



