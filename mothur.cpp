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

/**************************************************************************************************/

GlobalData* GlobalData::_uniqueInstance = 0;
CommandFactory* CommandFactory::_uniqueInstance = 0;

int main(int argc, char *argv[]){
	try {
		
		string log = "mothur.logFile";
		remove(log.c_str());
		
		time_t ltime = time(NULL); /* calendar time */  
		string logFileName = "mothur." + toString(ltime) + ".logfile";
		
		//version
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			system("clear");
			#if defined (__APPLE__) || (__MACH__)
				mothurOutJustToLog("Mac version");
				mothurOutEndLine(); mothurOutEndLine();
			#else
				mothurOutJustToLog("Linux version");
				mothurOutEndLine(); mothurOutEndLine();
			#endif

		#else
			system("CLS");
			mothurOutJustToLog("Windows version");
			mothurOutEndLine(); mothurOutEndLine();
		#endif		
		
		#ifdef USE_READLINE
			mothurOutJustToLog("Using ReadLine");
			mothurOutEndLine(); mothurOutEndLine();
		#endif
		
		//header
		mothurOut("mothur v.1.8");
		mothurOutEndLine();		
		mothurOut("Last updated: 2/02/2010");
		mothurOutEndLine();	
		mothurOutEndLine();		
		mothurOut("by");
		mothurOutEndLine();		
		mothurOut("Patrick D. Schloss");
		mothurOutEndLine();
		mothurOutEndLine();			
		mothurOut("Department of Microbiology & Immunology");
		mothurOutEndLine();	
		mothurOut("University of Michigan");
		mothurOutEndLine();			
		mothurOut("pschloss@umich.edu");
		mothurOutEndLine();		
		mothurOut("http://www.mothur.org");
		mothurOutEndLine();
		mothurOutEndLine();
		mothurOut("When using, please cite:");
		mothurOutEndLine();
		mothurOut("Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41.");
		mothurOutEndLine();	
		mothurOutEndLine();		
		mothurOut("Distributed under the GNU General Public License");
		mothurOutEndLine();
		mothurOutEndLine();			
		mothurOut("Type 'help()' for information on the commands that are available");
		mothurOutEndLine();
		mothurOutEndLine();			
		mothurOut("Type 'quit()' to exit program");
		mothurOutEndLine();	

				
		//srand(54321);
		srand( (unsigned)time( NULL ) );
		
		Engine* mothur;
		bool bail = 0;
		string input;

		if(argc>1){
			input = argv[1];

			if (input[0] == '#') {
				mothur = new ScriptEngine(argv[0], argv[1]);
			}else{
				mothur = new BatchEngine(argv[0], argv[1]);
			}
		}
		else{
			mothur = new InteractEngine(argv[0]);	
		}
		
		while(bail == 0)		{	bail = mothur->getInput();			}
		
		string outputDir = mothur->getOutputDir();
		logFileName = outputDir + logFileName;
	
		//need this because mothur.h makes the logfile, but doesn't know where to put it
		rename(log.c_str(), logFileName.c_str()); //logfile with timestamp
		
		delete mothur;

		return 0;
	}
	catch(exception& e) {
		errorOut(e, "mothur", "main");
		exit(1);
	}
}

/**************************************************************************************************/

