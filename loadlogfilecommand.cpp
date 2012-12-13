//
//  loadlogfilecommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 6/13/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "loadlogfilecommand.h"
#include "commandoptionparser.hpp"
#include "commandfactory.hpp"
#include "setcurrentcommand.h"

//**********************************************************************************************************************
vector<string> LoadLogfileCommand::setParameters(){	
	try {
        CommandParameter plogfile("logfile", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(plogfile);		
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "LoadLogfileCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string LoadLogfileCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The load.logfile command extracts the current file names from a logfile.\n";
		helpString += "The load.logfile parameter is logfile, and it is required.\n";
		helpString += "The load.logfile command should be in the following format: \n";
		helpString += "load.logfile(logfile=yourLogFile)\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "LoadLogfileCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
LoadLogfileCommand::LoadLogfileCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
	}
	catch(exception& e) {
		m->errorOut(e, "LoadLogfileCommand", "LoadLogfileCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
LoadLogfileCommand::LoadLogfileCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			//valid paramters for this command
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
            
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
                
                string path;
				it = parameters.find("logfile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["logfile"] = inputDir + it->second;		}
				}
				
            }
                        
            //get shared file, it is required
			logfile = validParameter.validFile(parameters, "logfile", true);
			if (logfile == "not open") { logfile = ""; abort = true; }	
			else if (logfile == "not found") { m->mothurOut("The logfile parameter is required."); m->mothurOutEndLine();abort = true;  }
            
            //if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = m->hasPath(logfile); //if user entered a file with a path then preserve it	
			}
        }
		
	}
	catch(exception& e) {
		m->errorOut(e, "NewCommand", "NewCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int LoadLogfileCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        m->mothurOutEndLine(); 
        m->mothurOut("Extracting current files names..."); m->mothurOutEndLine(); 
        m->mothurOutEndLine(); 
        
        CommandFactory* cFactory;
        cFactory = CommandFactory::getInstance();
        
        ifstream in;
        m->openInputFile(logfile, in);
        
        set<string> currentTypes = m->getCurrentTypes();
        map<string, string> currentFiles;
        string commandName = "";
        bool skip = false;
        string line = "";
        
        while (!in.eof()) {
            if (m->control_pressed) { break; }

            if (!skip) {  line = m->getline(in); m->gobble(in); }
            m->gobble(in);
             
            //look for "mothur >"
            int pos = line.find("mothur > "); //command line
            int pos2 = line.find("Output File "); //indicates command completed and we can update the current file
            int pos3 = line.find("/*****************"); 
            
            //skipping over parts where a command runs another command
            if (pos3 != string::npos) {
                while (!in.eof()) {
                    if (m->control_pressed) { break; }
                    line = m->getline(in); m->gobble(in); 
                    int posTemp = line.find("/*****************");
                    if (posTemp != string::npos) { break; }
                }
            }
            
            if (pos != string::npos) { 
                skip=false;
                //extract command name and option string
                string input = line.substr(pos+9);
                CommandOptionParser parser(input);
                commandName = parser.getCommandString();
                string options = parser.getOptionString();
                
                //parse out parameters in option string
                map<string,string> parameters;
                OptionParser optionParser(options, parameters);
                
                for (map<string, string>::iterator it = parameters.begin(); it != parameters.end(); it++) {
                    if (currentTypes.count((it->first)) != 0) { //if this is a type we save
                        if (it->second != "current") { currentFiles[it->first] = it->second; }//save the input file name as current
                    }
                }
            }else if (pos2 != string::npos) { 
                //read file output file names
                vector<string> theseOutputNames;
                while (!in.eof()) {
                    if (m->control_pressed) { break; }
                    line = m->getline(in); m->gobble(in); 
                    int pos = line.find("mothur > ");
                    if (pos != string::npos) { skip = true; break;  }
                    else {  theseOutputNames.push_back(line);  }
                }
                //ask command for the output names for each type based on inputs
                Command* command = cFactory->getCommand(commandName);
                map<string, vector<string> > thisOutputTypes = command->getOutputFiles();
                
                
                for (map<string, vector<string> >::iterator it = thisOutputTypes.begin(); it != thisOutputTypes.end(); it++) {
                    if (currentTypes.count((it->first)) != 0) {  //do we save this type
                        //if yes whats its tag
                        map<string, string>::iterator itCurrentFiles = currentFiles.find(it->first);
                        string thisTypesCurrentFile = "";
                        if (itCurrentFiles != currentFiles.end()) { thisTypesCurrentFile = itCurrentFiles->second;  }
                        
                        //outputfilename pattern for this input type
                        string pattern = command->getOutputPattern(it->first);
                        updateCurrent(pattern, it->first, thisTypesCurrentFile, theseOutputNames, currentFiles);
                        
                        //cout << "current=\n\n";
                        //for (map<string, string>::iterator itcc = currentFiles.begin(); itcc != currentFiles.end(); itcc++) {
                          //  cout << itcc->first << '\t' << itcc->second << endl;
                       // }
                    }
                }
            }
        }
        in.close();
        
        if (m->control_pressed) { return 0; }
                
        //output results
        string inputString = "";
        for (map<string, string>::iterator it = currentFiles.begin(); it != currentFiles.end(); it++) { inputString += it->first + "=" + it->second + ","; }
        
        if (inputString != "") { 
            inputString = inputString.substr(0, inputString.length()-1); 
            m->mothurOutEndLine(); 
			m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
			m->mothurOut("Running command: set.current(" + inputString + ")"); m->mothurOutEndLine(); 
			m->mothurCalling = true;
            
			Command* currentCommand = new SetCurrentCommand(inputString);
			currentCommand->execute();
			
            delete currentCommand;
			m->mothurCalling = false;
			m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
			

        } 
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "LoadLogfileCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int LoadLogfileCommand::updateCurrent(string pattern, string type, string thisTypesCurrentFile, vector<string> filenames, map<string, string>& currentFiles){
	try {
        
        vector<string> patterns; m->splitAtChar(pattern, patterns, '-');
        
        for (int i = 0; i < patterns.size(); i++) {
            
            vector<string> peices; m->splitAtChar(patterns[i], peices, ',');
            //cout << "patterns i = " << patterns[i] << endl;
            if (peices.size() != 0) {
                string tag = peices[peices.size()-1];
                //cout << "tag = " << tag << endl;
                if (peices[peices.size()-1] == "[extension]") {  tag = m->getExtension(thisTypesCurrentFile);  }
               
                //search for the tag in the list of output files
                for (int h = 0; h < filenames.size(); h++) {
                    
                    if (m->control_pressed) { return 0; }
                    //cout << "filename h = " << filenames[h]<< endl;
                    string ending = filenames[h].substr(filenames[h].length()-tag.length(), tag.length());
                    if (ending == tag) { //if it's there and this is a type we save a current version of, save it
                        if ((type == "column") || (type == "phylip")) { //check for format
                            string RippedName = "";
                            bool foundDot = false;
                            for (int i = filenames[h].length()-1; i >= 0; i--) {
                                if (foundDot && (filenames[h][i] != '.')) {  RippedName = filenames[h][i] + RippedName; }
                                else if (foundDot && (filenames[h][i] == '.')) {  break; }
                                else if (!foundDot && (filenames[h][i] == '.')) {  foundDot = true; }
                            }
                            if ((RippedName == "phylip") || (RippedName == "lt") || (RippedName == "square"))  {  currentFiles["phylip"] = filenames[h];  }
                            else {  currentFiles["column"] = filenames[h]; }
                        }else {  currentFiles[type] = filenames[h]; }
                        break;
                    }
                }
            }
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "LoadLogfileCommand", "updateCurrent");
		exit(1);
	}
}
//**********************************************************************************************************************


