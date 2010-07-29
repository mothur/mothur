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
#include "globaldata.hpp"
#include "mothurout.h"

/**************************************************************************************************/

GlobalData* GlobalData::_uniqueInstance = 0;
CommandFactory* CommandFactory::_uniqueInstance = 0;
MothurOut* MothurOut::_uniqueInstance = 0;

/***********************************************************************/
volatile int ctrlc_pressed = 0;
void ctrlc_handler ( int sig ) {
	MothurOut* m = MothurOut::getInstance();
    ctrlc_pressed = 1;
	m->control_pressed = ctrlc_pressed;
	
	if (m->executing) { //if mid command quit execution, else quit mothur
		m->mothurOutEndLine(); m->mothurOut("quitting command...");  m->mothurOutEndLine();
	}else{
		m->mothurOut("quitting mothur");  m->mothurOutEndLine();
		exit(1);
	}
}
/***********************************************************************/
int main(int argc, char *argv[]){
	MothurOut* m = MothurOut::getInstance();
	try {
		
		signal(SIGINT, ctrlc_handler );
				
		time_t ltime = time(NULL); /* calendar time */  
		string logFileName = "mothur." + toString(ltime) + ".logfile";
		
		#ifdef USE_MPI
			MPI_Init(&argc, &argv); 
		#endif

		m->setFileName(logFileName);
		
				
		//version
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			system("clear");
			#if defined (__APPLE__) || (__MACH__)
				m->mothurOutJustToLog("Mac version");
				m->mothurOutEndLine(); m->mothurOutEndLine();
			#else
				m->mothurOutJustToLog("Linux version");
				m->mothurOutEndLine(); m->mothurOutEndLine();
			#endif

		#else
			system("CLS");
			m->mothurOutJustToLog("Windows version");
			m->mothurOutEndLine(); m->mothurOutEndLine();
		#endif		
		
		#ifdef USE_READLINE
			m->mothurOutJustToLog("Using ReadLine");
			m->mothurOutEndLine(); m->mothurOutEndLine();
		#endif
		
		#ifdef MOTHUR_FILES
			string temp = MOTHUR_FILES; 
			
			//add / to name if needed
			string lastChar = temp.substr(temp.length()-1);
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				if (lastChar != "/") { temp += "/"; }
			#else
				if (lastChar != "\\") { temp += "\\"; }	
			#endif
			
			temp = getFullPathName(temp);
			m->setDefaultPath(temp);
			
			m->mothurOutJustToLog("Using default file location " + temp);
			m->mothurOutEndLine(); m->mothurOutEndLine();
		#endif
		
		#ifdef BIT_VERSION
			m->mothurOutJustToLog("Running 64Bit Version");
			m->mothurOutEndLine(); m->mothurOutEndLine();
		#else
			m->mothurOutJustToLog("Running 32Bit Version");
			m->mothurOutEndLine(); m->mothurOutEndLine();
		#endif
		
		//header
		m->mothurOut("mothur v.1.12.1");
		m->mothurOutEndLine();		
		m->mothurOut("Last updated: 7/29/2010");
		m->mothurOutEndLine();	
		m->mothurOutEndLine();		
		m->mothurOut("by");
		m->mothurOutEndLine();		
		m->mothurOut("Patrick D. Schloss");
		m->mothurOutEndLine();
		m->mothurOutEndLine();			
		m->mothurOut("Department of Microbiology & Immunology");
		m->mothurOutEndLine();	
		m->mothurOut("University of Michigan");
		m->mothurOutEndLine();			
		m->mothurOut("pschloss@umich.edu");
		m->mothurOutEndLine();		
		m->mothurOut("http://www.mothur.org");
		m->mothurOutEndLine();
		m->mothurOutEndLine();
		m->mothurOut("When using, please cite:");
		m->mothurOutEndLine();
		m->mothurOut("Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.");
		m->mothurOutEndLine();	
		m->mothurOutEndLine();		
		m->mothurOut("Distributed under the GNU General Public License");
		m->mothurOutEndLine();
		m->mothurOutEndLine();			
		m->mothurOut("Type 'help()' for information on the commands that are available");
		m->mothurOutEndLine();
		m->mothurOutEndLine();			
		m->mothurOut("Type 'quit()' to exit program");
		m->mothurOutEndLine();	
		
		#ifdef USE_MPI
			m->mothurOutJustToLog("Using MPI\tversion ");
			int version, subversion;
			MPI_Get_version(&version, &subversion);
			m->mothurOutJustToLog(toString(version) + "." + toString(subversion) + "\n");
		#endif
		
		//srand(54321);
		srand( (unsigned)time( NULL ) );
		
		Engine* mothur;
		bool bail = 0;
		string input;

		if(argc>1){
			input = argv[1];

			if (input[0] == '#') {
				m->mothurOutJustToLog("Script Mode");
				m->mothurOutEndLine(); m->mothurOutEndLine();

				mothur = new ScriptEngine(argv[0], argv[1]);
			}else{
				m->mothurOutJustToLog("Batch Mode");
				m->mothurOutEndLine(); m->mothurOutEndLine();
				
				mothur = new BatchEngine(argv[0], argv[1]);
			}
		}
		else{
			m->mothurOutJustToLog("Interactive Mode");
			m->mothurOutEndLine(); m->mothurOutEndLine();
			
			mothur = new InteractEngine(argv[0]);	
		}
		
		while(bail == 0)	{	bail = mothur->getInput();	}
		
		//closes logfile so we can rename
		m->closeLog();
		
		string outputDir = mothur->getOutputDir();
		string tempLog = mothur->getLogFileName();
		bool append = mothur->getAppend();
		
		string newlogFileName;
		if (tempLog != "") {
			newlogFileName = outputDir + tempLog;
			
			if (!append) {	
				//need this because m->mothurOut makes the logfile, but doesn't know where to put it
				rename(logFileName.c_str(), newlogFileName.c_str()); //logfile with timestamp

			}else {
				ofstream outNewLog;
				openOutputFileAppend(newlogFileName, outNewLog);
				outNewLog << endl << endl << "*********************************************************************************" << endl << endl;
				outNewLog.close();
				
				appendFiles(logFileName, newlogFileName);
				remove(logFileName.c_str());
			}
		}else{  
			newlogFileName = outputDir + logFileName;
			//need this because m->mothurOut makes the logfile, but doesn't know where to put it
			rename(logFileName.c_str(), newlogFileName.c_str()); //logfile with timestamp
		}
		
				
		delete mothur;
		
		#ifdef USE_MPI
			MPI_Finalize();
		#endif
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "mothur", "main");
		exit(1);
	}
}

/**************************************************************************************************/

