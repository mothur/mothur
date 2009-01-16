/*
 *  raredisplay.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 11/18/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "raredisplay.h"

/***********************************************************************/

void RareDisplay::init(string label){
	try {
		this->label = label;
		if(nIters != 1){
			openOutputFile(tempOutName, tempOutFile);
			openInputFile(tempInName, tempInFile);
		}
		else{
			openOutputFile(tempOutName, tempOutFile);
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the RareDisplay class Function init. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the RareDisplay class function init. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

void RareDisplay::update(SAbundVector* rank){
	try {
		int newNSeqs = rank->getNumSeqs();
		vector<double> data = estimate->getValues(rank);

		if(nIters != 1){

			double oldNSeqs, oldMean, oldVar;
		
			tempInFile >> oldNSeqs >> oldMean >> oldVar;
		
			double oldS = oldVar * ( nIters - 2 );
			double delta = data[0] - oldMean;
			double newMean = oldMean + delta / nIters;
			double newS = oldS + delta * ( data[0] - newMean );
			double newVar = newS / ( nIters - 1 );

			tempOutFile << newNSeqs << '\t' << newMean << '\t' << newVar << endl;
		}
		else{
			tempOutFile << newNSeqs << '\t' << data[0] << '\t' << 0 << endl;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the RareDisplay class Function update. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the RareDisplay class function update. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
};

/***********************************************************************/
void RareDisplay::update(SharedRAbundVector* shared1, SharedRAbundVector* shared2, int numSeqs, int numGroupComb) {
	try {
		vector<double> data = estimate->getValues(shared1, shared2); 
		double newNSeqs = data[0];
		
		if(nIters != 1){
			double oldNSeqs, oldMean, oldVar;
		
			tempInFile >> oldNSeqs >> oldMean >> oldVar;
		
			double oldS = oldVar * ( nIters - 2 );
			double delta = data[0] - oldMean;
			double newMean = oldMean + delta / nIters;
			double newS = oldS + delta * ( data[0] - newMean );
			double newVar = newS / ( nIters - 1 );

			tempOutFile << newNSeqs << '\t' << newMean << '\t' << newVar << endl;
		}
		else{
			tempOutFile << newNSeqs << '\t' << data[0] << '\t' << 0 << endl;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the RareDisplay class Function update. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the RareDisplay class function update. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
};

/***********************************************************************/

void RareDisplay::reset(){
	try {
		if(nIters != 1){
			tempOutFile.close();
			tempInFile.close();
		}
		else{
			tempOutFile.close();
		}
	
		nIters++;
	
		remove(tempInName.c_str());
		rename(tempOutName.c_str(), tempInName.c_str());	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the RareDisplay class Function reset. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the RareDisplay class function reset. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

void RareDisplay::close(){
	try {
		output->initFile(label);
	
		openInputFile(tempInName, tempInFile);
	
		while(!tempInFile.eof()){
			int nSeqs;
			tempInFile >> nSeqs;
		
			vector<double> data(3,0);
			double variance = 0;
		
			tempInFile >> data[0];
			tempInFile >> variance;
		
			double ci = 1.96 * pow(variance, 0.5);
			data[1] = data[0] - ci;
			data[2] = data[0] + ci;
		
			output->output(nSeqs, data);
		
			gobble(tempInFile);
		}
		tempInFile.close();
	
		remove(tempInName.c_str());
		remove(tempOutName.c_str());
	
		nIters = 1;
		output->resetFile();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the RareDisplay class Function close. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the RareDisplay class function close. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

