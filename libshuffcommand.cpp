/*
 *  libshuffcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "libshuffcommand.h"

//**********************************************************************************************************************


LibShuffCommand::LibShuffCommand(){
	try {
		//globaldata = GlobalData::getInstance();
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the LibShuffCommand class Function LibShuffCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the LibShuffCommand class function LibShuffCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
			
}

//**********************************************************************************************************************

LibShuffCommand::~LibShuffCommand(){
	
}

//**********************************************************************************************************************

int LibShuffCommand::execute(){
	try {
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the LibShuffCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the LibShuffCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************
