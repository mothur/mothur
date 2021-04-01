/*
 *  interface.cpp
 *  
 *
 *  Created by Pat Schloss on 8/14/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

#include "batchengine.hpp"
#include "scriptengine.hpp"
#include "interactengine.hpp"
#include "mothurout.h"

/**************************************************************************************************/

CommandFactory* CommandFactory::_uniqueInstance = 0;
MothurOut* MothurOut::_uniqueInstance = 0;
CurrentFile* CurrentFile::instance = 0;
/***********************************************************************/
volatile int ctrlc_pressed = 0;
void ctrlc_handler ( int sig ) {
	MothurOut* m = MothurOut::getInstance();
    ctrlc_pressed = 1;
	m->setControl_pressed(ctrlc_pressed);
	
	if (m->getExecuting()) { //if mid command quit execution, else quit mothur
        m->mothurOut("\nquitting command...\n");
	}else{
		m->mothurOut("quitting mothur\n");
		exit(1);
	}
}
/***********************************************************************/
int main(int argc, char *argv[], char *envp[]){
	MothurOut* m = MothurOut::getInstance();
	try {
        CurrentFile* current = CurrentFile::getInstance();
        Utils util;
        bool createLogFile = true;
        
		signal(SIGINT, ctrlc_handler );
        
        string defaultPath, mothurVersion, releaseDate, OS, toolsPath;
        util.mothurInitialPrep(defaultPath, toolsPath, mothurVersion, releaseDate, OS);
        
        current->setReleaseDate(releaseDate);
        current->setVersion(mothurVersion);
		
		#ifdef MOTHUR_FILES
			current->setDefaultPath(defaultPath);
		#endif
        
        #ifdef MOTHUR_TOOLS
            current->setToolsPath(toolsPath);
        #endif
        
		if (argc>1) {
            if (argc > 2) { //is one of these -q for quiet mode?
                if (argc > 3) { m->appendLogBuffer("[ERROR]: mothur only allows command inputs and the -q command line options.\n  i.e. ./mothur \"#summary.seqs(fasta=final.fasta);\" -q\n or ./mothur -q \"#summary.seqs(fasta=final.fasta);\"\n"); return 0; }
                else {
                    string argv1 = argv[1];
                    string argv2 = argv[2];
                    if ((argv1 == "--quiet") || (argv1 == "-q")) {
                        m->setQuietMode(true);
                        argv[1] = argv[2];
                    }else if ((argv2 == "--quiet") || (argv2 == "-q")) {
                         m->setQuietMode(true);
                    }else {
                        m->appendLogBuffer("[ERROR]: mothur only allows command inputs and the -q command line options.\n");
                        m->appendLogBuffer("[ERROR]: Unrecognized options: " + argv1 + " " + argv2 + "\n");
                        return 0;
                    }
                }
            }
		}
        
        map<string, string> environmentalVariables;
        for (char **env = envp; *env != 0; env++){
            string thisEvn = *env;
            string key, value; value = thisEvn;
            util.splitAtEquals(key, value);
            
            map<string, string>::iterator it = environmentalVariables.find(key);
            if (it == environmentalVariables.end())     { environmentalVariables[key] = value;  }
            else                                        { it->second = value;                   }
            //m->mothurOut("[DEBUG]: Setting environment variable " + key + " to " + value + "\n"); 
            if (m->getDebug()) { m->mothurOut("[DEBUG]: Setting environment variable " + key + " to " + value + "\n"); }
        }
        
		Engine* mothur = NULL;
		bool bail = false;
		string input;
 
		if(argc>1){
			input = argv[1];
			if (input[0] == '#') {
				m->appendLogBuffer("Script Mode\n\n");
				mothur = new ScriptEngine(argv[0], argv[1], environmentalVariables);
			}else if ((input == "--version") || (input == "-v")) {
				cout << (OS + "\nMothur version=" + mothurVersion + "\nRelease Date=" + releaseDate + "\n\n"); return 0;
            }else if ((input == "--help") || (input == "-h")) {
                createLogFile = false;
                m->appendLogBuffer("Script Mode\n\n");

                string helpQuit = "#help();quit();";
               
                argv[1] = util.mothurConvert(helpQuit);
                mothur = new ScriptEngine(argv[0], argv[1], environmentalVariables);
			}else{
				m->appendLogBuffer("Batch Mode\n\n");
                mothur = new BatchEngine(argv[0], argv[1], environmentalVariables);
			}
		}else{
			m->appendLogBuffer("Interactive Mode\n\n");
            mothur = new InteractEngine(argv[0], environmentalVariables);
		}
		
		while(!bail)	{   bail = mothur->getInput();	}
		
		string newlogFileName = mothur->getLogFileName();
		        
        if (!createLogFile) { util.mothurRemove(newlogFileName); }
				
		if (mothur != NULL) { delete mothur; }
        
        int returnCode = 0;
        if (m->getNumErrors() != 0) { returnCode = 1; }
        m->closeLog();
        
		return returnCode;
	}
	catch(exception& e) {
		m->errorOut(e, "mothur", "main");
		exit(1);
	}
}

/**************************************************************************************************/

