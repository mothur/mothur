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

		map<int, vector<double> >::iterator it = results.find(newNSeqs);
        if (it == results.end()) { //first iter for this count
            vector<double> temp;
            temp.push_back(data[0]);
            results[newNSeqs] = temp;
        }else {
            it->second.push_back(data[0]);
        }
	}
	catch(exception& e) {
		m->errorOut(e, "RareDisplay", "update");
		exit(1);
	}
}

/***********************************************************************/
void RareDisplay::update(vector<SharedRAbundVector*> shared, int numSeqs, int numGroupComb, vector<string> g) {
	try {
		vector<double> data = estimate->getValues(shared);
        Groups = g;
		
		map<int, vector<double> >::iterator it = results.find(numSeqs);
        if (it == results.end()) { //first iter for this count
            vector<double> temp;
            temp.push_back(data[0]);
            results[numSeqs] = temp;
        }else {
            it->second.push_back(data[0]);
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
	}
	catch(exception& e) {
		m->errorOut(e, "RareDisplay", "reset");
		exit(1);
	}
}

/***********************************************************************/

void RareDisplay::close(){
	try {
		output->setLabelName(label);
	
		for (map<int, vector<double> >::iterator it = results.begin(); it != results.end(); it++) {
		
			vector<double> data(3,0);
            
            sort((it->second).begin(), (it->second).end());
            
            vector<double> thisResults = it->second;
            double meanResults = util.getAverage(thisResults);
			data[0] = meanResults;
			data[1] = (it->second)[(int)(0.025*(nIters-1))];
			data[2] = (it->second)[(int)(0.975*(nIters-1))];
		
			output->updateOutput(it->first, data);
		}
		
		nIters = 1;
        results.clear();
		
		output->resetFile();
	}
	catch(exception& e) {
		m->errorOut(e, "RareDisplay", "close");
		exit(1);
	}
}
/***********************************************************************/
