/*
 *  distancecommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "distancecommand.h"
#include "ignoregaps.h"
#include "eachgapdist.h"
#include "eachgapignore.h"
#include "onegapdist.h"
#include "onegapignore.h"

//**********************************************************************************************************************

DistanceCommand::DistanceCommand(){
	try {
		globaldata = GlobalData::getInstance();
		validCalculator = new ValidCalculators();
		countends = globaldata->getCountEnds();
		convert(globaldata->getProcessors(), processors);
		convert(globaldata->getCutOff(), cutoff);
		
		int i;
		if (countends == "T") {
			for (i=0; i<globaldata->Estimators.size(); i++) {
				if (validCalculator->isValidCalculator("distance", globaldata->Estimators[i]) == true) { 
					if (globaldata->Estimators[i] == "nogaps") { 
						distCalculator = new ignoreGaps();
					}else if (globaldata->Estimators[i] == "eachgap") { 
						distCalculator = new eachGapDist();	
					}else if (globaldata->Estimators[i] == "onegap") {
					distCalculator = new oneGapDist();					}
				}
			}
		}else {
			for (i=0; i<globaldata->Estimators.size(); i++) {
				if (validCalculator->isValidCalculator("distance", globaldata->Estimators[i]) == true) { 
					if (globaldata->Estimators[i] == "nogaps") { 
						distCalculator = new ignoreGaps();	
					}else if (globaldata->Estimators[i] == "eachgap") { 
						distCalculator = new eachGapIgnoreTermGapDist();
					}else if (globaldata->Estimators[i] == "onegap") { 
						distCalculator = new oneGapIgnoreTermGapDist();	
					}
				}
			}
		}
		
		//reset calc for next command
		globaldata->setCalc("");
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the DistanceCommand class Function DistanceCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the DistanceCommand class function DistanceCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************

int DistanceCommand::execute(){
	try {
		
		//read file
		string filename = globaldata->inputFileName;
		
		if(globaldata->getFastaFile() != "") {
		readSeqs =  new ReadFasta(filename); }
		else if(globaldata->getNexusFile() != "") {
		readSeqs = new ReadNexus(filename); }
		else if(globaldata->getClustalFile() != "") {
		readSeqs = new ReadClustal(filename); }
		else if(globaldata->getPhylipFile() != "") {
		readSeqs = new ReadPhylip(filename); }
		
		readSeqs->read();
		seqDB = readSeqs->getDB();
		
		int numSeqs = seqDB->getNumSeqs();
		cutoff += 0.005;
		
		string distFile = getRootName(globaldata->getFastaFile()) + "dist";
		
		remove(distFile.c_str());
		
		//#	if defined (_WIN32)
		//figure out how to implement the fork and wait commands in windows
		//	driver(distCalculator, seqDB, 0, numSeqs, distFile, cutoff);
		//#	endif
		
#if defined (__APPLE__) || (__MACH__)
		if(processors == 1){
			driver(distCalculator, seqDB, 0, numSeqs, distFile, cutoff);	
		}	
		else if(processors == 2){
			
			int pid = fork();
			if(pid > 0){
				driver(distCalculator, seqDB, 0, (numSeqs/sqrt(2)), distFile + "tempa", cutoff);	
				appendFiles((distFile+"tempa"), distFile);
				remove((distFile + "tempa").c_str());	
			}
			else{
				driver(distCalculator, seqDB, (numSeqs/sqrt(2)), numSeqs, distFile + "tempb", cutoff);	
				appendFiles((distFile+"tempb"), distFile);
				remove((distFile + "tempb").c_str());	
			}
			wait(NULL);
			
		}
		else if(processors == 3){
			int pid1 = fork();
			if(pid1 > 0){
				int pid2 = fork();
				if(pid2 > 0){
					driver(distCalculator, seqDB, 0, sqrt(3) * numSeqs / 3, distFile + "tempa", cutoff);
					appendFiles(distFile+"tempa", distFile);
					remove((distFile + "tempa").c_str());	
				}
				else{
					driver(distCalculator, seqDB, sqrt(3) * numSeqs / 3, sqrt(6) * numSeqs / 3, distFile + "tempb", cutoff);	
					appendFiles(distFile+"tempb", distFile);
					remove((distFile + "tempb").c_str());				
				}
				wait(NULL);
			}
			else{
				driver(distCalculator, seqDB, sqrt(6) * numSeqs / 3, numSeqs, distFile + "tempc", cutoff);	
				appendFiles(distFile+"tempc", distFile);
				remove((distFile + "tempc").c_str());			
			}
			wait(NULL);
		}
		else if(processors == 4){
			int pid1 = fork();
			if(pid1 > 0){
				int pid2 = fork();
				if(pid2 > 0){
					driver(distCalculator, seqDB, 0, numSeqs / 2, distFile + "tempa", cutoff);	
					appendFiles(distFile+"tempa", distFile);
					remove((distFile + "tempa").c_str());			
				}
				else{
					driver(distCalculator, seqDB, numSeqs / 2, (numSeqs/sqrt(2)), distFile + "tempb", cutoff);	
					appendFiles(distFile+"tempb", distFile);
					remove((distFile + "tempb").c_str());				
				}
				wait(NULL);
			}
			else{
				int pid3 = fork();
				if(pid3 > 0){
					driver(distCalculator, seqDB, (numSeqs/sqrt(2)), (sqrt(3) * numSeqs / 2), distFile + "tempc", cutoff);	
					appendFiles(distFile+"tempc", distFile);
					remove((distFile + "tempc").c_str());				
				}
				else{
					driver(distCalculator, seqDB, (sqrt(3) * numSeqs / 2), numSeqs, distFile + "tempd", cutoff);	
					appendFiles(distFile+"tempd", distFile);
					remove((distFile + "tempd").c_str());				
				}
				wait(NULL);
			}
			wait(NULL);
		}
		wait(NULL);
#elif (linux) || (__linux)
		if(processors == 1){
			driver(distCalculator, seqDB, 0, numSeqs, distFile, cutoff);	
		}	
		else if(processors == 2){
			
			int pid = fork();
			if(pid > 0){
				driver(distCalculator, seqDB, 0, (numSeqs/sqrt(2)), distFile + "tempa", cutoff);	
				appendFiles((distFile+"tempa"), distFile);
				remove((distFile + "tempa").c_str());	
			}
			else{
				driver(distCalculator, seqDB, (numSeqs/sqrt(2)), numSeqs, distFile + "tempb", cutoff);	
				appendFiles((distFile+"tempb"), distFile);
				remove((distFile + "tempb").c_str());	
			}
			wait();
			
		}
		else if(processors == 3){
			int pid1 = fork();
			if(pid1 > 0){
				int pid2 = fork();
				if(pid2 > 0){
					driver(distCalculator, seqDB, 0, sqrt(3) * numSeqs / 3, distFile + "tempa", cutoff);
					appendFiles(distFile+"tempa", distFile);
					remove((distFile + "tempa").c_str());	
				}
				else{
					driver(distCalculator, seqDB, sqrt(3) * numSeqs / 3, sqrt(6) * numSeqs / 3, distFile + "tempb", cutoff);	
					appendFiles(distFile+"tempb", distFile);
					remove((distFile + "tempb").c_str());				
				}
				wait();
			}
			else{
				driver(distCalculator, seqDB, sqrt(6) * numSeqs / 3, numSeqs, distFile + "tempc", cutoff);	
				appendFiles(distFile+"tempc", distFile);
				remove((distFile + "tempc").c_str());			
			}
			wait();
		}
		else if(processors == 4){
			int pid1 = fork();
			if(pid1 > 0){
				int pid2 = fork();
				if(pid2 > 0){
					driver(distCalculator, seqDB, 0, numSeqs / 2, distFile + "tempa", cutoff);	
					appendFiles(distFile+"tempa", distFile);
					remove((distFile + "tempa").c_str());			
				}
				else{
					driver(distCalculator, seqDB, numSeqs / 2, (numSeqs/sqrt(2)), distFile + "tempb", cutoff);	
					appendFiles(distFile+"tempb", distFile);
					remove((distFile + "tempb").c_str());				
				}
				wait();
			}
			else{
				int pid3 = fork();
				if(pid3 > 0){
					driver(distCalculator, seqDB, (numSeqs/sqrt(2)), (sqrt(3) * numSeqs / 2), distFile + "tempc", cutoff);	
					appendFiles(distFile+"tempc", distFile);
					remove((distFile + "tempc").c_str());				
				}
				else{
					driver(distCalculator, seqDB, (sqrt(3) * numSeqs / 2), numSeqs, distFile + "tempd", cutoff);	
					appendFiles(distFile+"tempd", distFile);
					remove((distFile + "tempd").c_str());				
				}
				wait();
			}
			wait();
		}
		wait();
		
#else
		driver(distCalculator, seqDB, 0, numSeqs, distFile, cutoff);
#endif
		
		delete distCalculator;
		
		return 0;
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the DistanceCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the DistanceCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/**************************************************************************************************/
/////// need to fix to work with calcs and sequencedb
int DistanceCommand::driver(Dist* distCalculator, SequenceDB* align, int startLine, int endLine, string dFileName, float cutoff){
	try {
		int startTime = time(NULL);
		
		ofstream distFile(dFileName.c_str(), ios::trunc);
		distFile.setf(ios::fixed, ios::showpoint);
		distFile << setprecision(4);
		
		for(int i=startLine;i<endLine;i++){
			
			for(int j=0;j<i;j++){
				distCalculator->calcDist(align->get(i), align->get(j));
				double dist = distCalculator->getDist();
				
				if(dist <= cutoff){
					distFile << align->get(i).getName() << ' ' << align->get(j).getName() << ' ' << dist << endl;
				}
				
			}
			if(i % 100 == 0){
				cout << i << '\t' << time(NULL) - startTime << endl;
			}
			
		}
		cout << endLine-1 << '\t' << time(NULL) - startTime << endl;
		
		return 1;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the DistanceCommand class Function driver. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the DistanceCommand class function driver. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
	
}

/**************************************************************************************************/
void DistanceCommand::appendFiles(string temp, string filename) {
	try{
		ofstream output;
		ifstream input;
		
		//open output file in append mode
		openOutputFileAppend(filename, output);
		
		//open temp file for reading
		openInputFile(temp, input);
		
		string line;
		//read input file and write to output file
		while(input.eof() != true) {
			getline(input, line); //getline removes the newline char
			if (line != "") {
				output << line << endl;   // Appending back newline char 
			}
		}	
		
		input.close();
		output.close();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the DistanceCommand class Function appendFiles. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the DistanceCommand class function appendFiles. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
/**************************************************************************************************/