/*
 *  progress.cpp
 *  
 *
 *  Created by Pat Schloss on 8/14/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */


#include "progress.hpp"

const int totalTicks = 50;
const char marker = '|';


/***********************************************************************/

Progress::Progress(){
	try {
		m = MothurOut::getInstance();
		m->mothurOut("********************#****#****#****#****#****#****#****#****#****#****#");
		
		nTicks = 0;
		finalPos = 0;
	}
	catch(exception& e) {
		m->errorOut(e, "Progress", "Progress");
		exit(1);
	}
}

/***********************************************************************/

Progress::Progress(string job, int end){
	try {
		m = MothurOut::getInstance();
		
		m->mothurOut("********************#****#****#****#****#****#****#****#****#****#****#\n");
		cout << setw(20) << left << job << setw(1) << marker;
		m->mothurOutJustToLog(job);
		m->mothurOut(toString(marker));
		cout.flush();

		nTicks = 0;
		finalPos = end;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Progress", "Progress");
		exit(1);
	}
}

/***********************************************************************/

void Progress::newLine(string job, int end){
	try {
		m->mothurOutEndLine();
		cout << setw(20) << left << job << setw(1) << marker;
		m->mothurOutJustToLog(job);
		m->mothurOut(toString(marker));
		cout.flush();
		
		nTicks = 0;
		finalPos = end;
	}
	catch(exception& e) {
		m->errorOut(e, "Progress", "newLine");
		exit(1);
	}
}
	
/***********************************************************************/

void Progress::update(const int currentPos){
	try {
		int ratio = int(totalTicks * (float)currentPos / finalPos);
	
		if(ratio > nTicks){
			for(int i=nTicks;i<ratio;i++){
				m->mothurOut(toString(marker));
				cout.flush();
			}
			nTicks = ratio;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "Progress", "update");
		exit(1);
	}
}

/***********************************************************************/

void Progress::finish(){
	try {
		for(int i=nTicks;i<totalTicks;i++){
			m->mothurOut(toString(marker));
			cout.flush();
		}
	
	
		m->mothurOutEndLine();
		m->mothurOut("***********************************************************************\n");
		cout.flush();
	}
	catch(exception& e) {
		m->errorOut(e, "Progress", "finish");
		exit(1);
	}
}

/***********************************************************************/
