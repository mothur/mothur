/*
 *  raredisplay.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 11/18/08.
 *  Copyright 2008 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "raredisplay.h"

/***********************************************************************/

void RareDisplay::init(string label){
	try {
		this->label = label;
	}
	catch(exception& e) {
		m->errorOut(e, "RareDisplay", "init");
		exit(1);
	}
}

/***********************************************************************/

void RareDisplay::update(SAbundVector* rank){
	try {
		int newNSeqs = rank->getNumSeqs();
		vector<double> data = estimate->getValues(rank);

		if(nIters != 1){

			double oldS = var[index] * ( nIters - 2 );
			double delta = data[0] - results[index];
			double newMean = results[index] + delta / nIters;
			double newS = oldS + delta * ( data[0] - newMean );
			double newVar = newS / ( nIters - 1 );

			seqs[index] = newNSeqs;
			results[index] = newMean; 
			var[index] = newVar;
			
			index++;  
		}
		else{
			seqs.push_back(newNSeqs); 
			results.push_back(data[0]);
			var.push_back(0.0);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "RareDisplay", "update");
		exit(1);
	}
}

/***********************************************************************/
void RareDisplay::update(vector<SharedRAbundVector*> shared, int numSeqs, int numGroupComb) {
	try {
		vector<double> data = estimate->getValues(shared); 
		double newNSeqs = data[0];
		
		if(nIters != 1){
		
			double oldS = var[index] * ( nIters - 2 );
			double delta = data[0] - results[index];
			double newMean = results[index] + delta / nIters;
			double newS = oldS + delta * ( data[0] - newMean );
			double newVar = newS / ( nIters - 1 );
			seqs[index] = newNSeqs;
			results[index] = newMean; 
			var[index] = newVar;
			
			index++;  
		}
		else{
			
			seqs.push_back(newNSeqs); 
			results.push_back(data[0]);
			var.push_back(0.0);

		}
	}
	catch(exception& e) {
		m->errorOut(e, "RareDisplay", "update");
		exit(1);
	}
}

/***********************************************************************/

void RareDisplay::reset(){
	try {
		nIters++;
		index = 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RareDisplay", "reset");
		exit(1);
	}
}

/***********************************************************************/

void RareDisplay::close(){
	try {
		
		output->initFile(label);
	
		for (int i = 0; i < seqs.size(); i++) {
		
			vector<double> data(3,0);
			double variance = var[i];
			
			data[0] = results[i];
			
			double ci = 1.96 * pow(variance, 0.5);
			data[1] = data[0] - ci;
			data[2] = data[0] + ci;
		
			output->output(seqs[i], data);
		}
		
		nIters = 1;
		index = 0;
		
		seqs.clear();
		results.clear();
		var.clear();
		
		output->resetFile();
	}
	catch(exception& e) {
		m->errorOut(e, "RareDisplay", "close");
		exit(1);
	}
}

/***********************************************************************/

