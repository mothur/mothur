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
		cout << "Standard Error: " << e.what() << " has occurred in the dotur class Function main. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the dotur class function main. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/**************************************************************************************************/

