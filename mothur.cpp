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

int main(int argc, char *argv[]){
	try {
		
		//remove old logfile
		string logFileName = "mothur.logFile";
		remove(logFileName.c_str());
		
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

		
		//header
		mothurOut("mothur v.1.4.1");
		mothurOutEndLine();		
		mothurOut("Last updated: 6/23/2009");
		mothurOutEndLine();	
		mothurOutEndLine();		
		mothurOut("by");
		mothurOutEndLine();		
		mothurOut("Patrick D. Schloss");
		mothurOutEndLine();
		mothurOutEndLine();			
		mothurOut("Department of Microbiology");
		mothurOutEndLine();		
		mothurOut("pschloss@micro.umass.edu");
		mothurOutEndLine();		
		mothurOut("http://schloss.micro.umass.edu/mothur");
		mothurOutEndLine();	
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
	
		delete mothur;
	
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "mothur", "main");
		exit(1);
	}
}

/**************************************************************************************************/

