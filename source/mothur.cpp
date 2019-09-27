/*
 *  interface.cpp
 *  
 *
 *  Created by Pat Schloss on 8/14/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */
 
#include "mothur.h"
#include "engine.hpp"
#include "mothurout.h"
#include "currentfile.h"

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
int main(int argc, char *argv[]){
	MothurOut* m = MothurOut::getInstance();
	try {
        CurrentFile* current = CurrentFile::getInstance();
        Utils util;
        bool createLogFile = true;
        
		signal(SIGINT, ctrlc_handler );
				
		#if defined NON_WINDOWS
			system("clear");
		#else
			system("CLS");
		#endif
		
		#ifdef MOTHUR_FILES
			string temp = MOTHUR_FILES; 
		
			//add / to name if needed
			string lastChar = temp.substr(temp.length()-1);
			if (lastChar != PATH_SEPARATOR) { temp += PATH_SEPARATOR; }
		
			temp = util.getFullPathName(temp);
			current->setDefaultPath(temp);
		#endif
        
        #ifdef LOGFILE_NAME
            string logfilename = LOGFILE_NAME;
            logfilename = util.getFullPathName(logfilename);
        
            m->appendLogBuffer("Using Static Logfile " + logfilename +  "\n");
        
            m->setLogFileName(logfilename, false);
            m->mothurOut("\n");
        #endif
        
        string releaseDate = "";
        #ifdef RELEASE_DATE
            releaseDate = RELEASE_DATE;
        #else
            string year, month, day;
            util.getCurrentDate(year, month, day);
            releaseDate = month + "/" + day + "/" + year;
        #endif
        
		//get releaseDate from Make
		 
		string mothurVersion = VERSION; 
		current->setReleaseDate(releaseDate);
		current->setVersion(mothurVersion);
		
		//will make the gui output "pretty"
		bool outputHeader = true;
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
		
		if (outputHeader)  {
			//version
			#if defined NON_WINDOWS
				#if defined (__APPLE__) || (__MACH__)
					m->appendLogBuffer("Mac version\n\n");
				#else
					m->appendLogBuffer("Linux version\n\n");
				#endif
			#else
				m->appendLogBuffer("Windows version\n\n");
			#endif		
			
            string packagesUsed = "";
			#ifdef USE_READLINE
                packagesUsed += "ReadLine,";
			#endif
            
            #ifdef USE_BOOST
                packagesUsed += "Boost,";
            #endif
            
            #ifdef USE_HDF5
                packagesUsed += "HDF5,";
            #endif
            
            if (packagesUsed != "") {
                //remove last comma
                packagesUsed = packagesUsed.substr(0,packagesUsed.length()-1);
                m->appendLogBuffer("Using " + packagesUsed + "\n");
            }
			
			#ifdef MOTHUR_FILES
				m->appendLogBuffer("\nUsing default file location " + temp + "\n\n");
			#endif
			
			//header
			m->appendLogBuffer("mothur v." + mothurVersion + "\n");
			m->appendLogBuffer("Last updated: " + releaseDate + "\n");
			m->appendLogBuffer("by\n");
			m->appendLogBuffer("Patrick D. Schloss\n\n");
			m->appendLogBuffer("Department of Microbiology & Immunology\n\n");
			m->appendLogBuffer("University of Michigan\n");
			m->appendLogBuffer("http://www.mothur.org\n\n");
			m->appendLogBuffer("When using, please cite:\n");
			m->appendLogBuffer("Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.\n\n");
			m->appendLogBuffer("Distributed under the GNU General Public License\n\n");
			m->appendLogBuffer("Type 'help()' for information on the commands that are available\n\n");
			m->appendLogBuffer("For questions and analysis support, please visit our forum at https://forum.mothur.org\n\n");
			m->appendLogBuffer("Type 'quit()' to exit program\n\n");
		}
		
        m->setRandomSeed(19760620);
        m->appendLogBuffer("[NOTE]: Setting random seed to 19760620.\n\n");
        
		Engine* mothur = NULL;
		bool bail = 0;
		string input;
 
		if(argc>1){
			input = argv[1];
			if (input[0] == '#') {
				m->appendLogBuffer("Script Mode\n\n");
				mothur = new ScriptEngine(argv[0], argv[1]);
			}else if ((input == "--version") || (input == "-v")) {
                createLogFile = false;
                string OS = "";
                //version
                #if defined NON_WINDOWS
                #if defined (__APPLE__) || (__MACH__)
                OS = "Mac ";
                #else
                OS = "Linux ";
                #endif
                
                #else
                OS = "Windows ";
                #endif
                
				cout << (OS + "\nMothur version=" + mothurVersion + "\nRelease Date=" + releaseDate + "\n\n");
				return 0;
                
            }else if ((input == "--help") || (input == "-h")) {
                createLogFile = false;
                m->appendLogBuffer("Script Mode\n\n");

                char* temp = new char[16];
                *temp = '\0'; strncat(temp, "#help();quit();", 15);
                
                argv[1] = temp;
                mothur = new ScriptEngine(argv[0], argv[1]);
			}else{
				m->appendLogBuffer("Batch Mode\n\n");
                mothur = new BatchEngine(argv[0], argv[1]);
			}
		}else{
			m->appendLogBuffer("Interactive Mode\n\n");
            mothur = new InteractEngine(argv[0]);
		}
		
		while(bail == 0)	{	bail = mothur->getInput();	}
		
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

