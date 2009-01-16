/*
 *  progress.cpp
 *  
 *
 *  Created by Pat Schloss on 8/14/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

#include <iostream>
#include <string>
#include "progress.hpp"
#include <exception>

using namespace std;

const int totalTicks = 50;
const char marker = '|';


/***********************************************************************/

Progress::Progress(string job, int end){
	try {
		cout << "*******************#****#****#****#****#****#****#****#****#****#****#\n";
		cout << job << marker;
		cout.flush();
	
		nTicks = 0;
		finalPos = end;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Progress class Function Progress. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Progress class function Progress. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

void Progress::update(const int currentPos){
	try {
		int ratio = int(totalTicks * (float)currentPos / finalPos);
	
		if(ratio > nTicks){
			for(int i=nTicks;i<ratio;i++){
				cout << marker;
				cout.flush();
			}
			nTicks = ratio;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Progress class Function update. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Progress class function update. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

void Progress::finish(){
	try {
		for(int i=nTicks;i<totalTicks;i++){
			cout << marker;
			cout.flush();
		}
	
	
		cout << endl;
		cout << "**********************************************************************\n";
		cout.flush();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Progress class Function finish. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Progress class function finish. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/
