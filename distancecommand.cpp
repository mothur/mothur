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
		
		//open file
		string filename = globaldata->getFastaFile();
		openInputFile(filename, in);

		
		int i;
		if (isTrue(countends) == true) {
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
		
		//reads fasta file and fills sequenceDB
		if(globaldata->getFastaFile() != "") {  seqDB = new SequenceDB(in);  }
		else { cout << "Error no fasta file." << endl; return 0; }
				
		int numSeqs = seqDB->getNumSeqs();
		cutoff += 0.005;
		
		string distFile = getRootName(globaldata->getFastaFile()) + "dist";
		string phylipFile = getRootName(globaldata->getFastaFile()) + "phylip.dist";
		
		remove(phylipFile.c_str());
		remove(distFile.c_str());
		
		//output numSeqs to phylip formatted dist file
		openOutputFile(phylipFile, phylipOut);
		phylipOut << numSeqs << endl;
		phylipOut.close();
		
		//#	if defined (_WIN32)
		//figure out how to implement the fork and wait commands in windows
		//	driver(distCalculator, seqDB, 0, numSeqs, distFile, phylipFile, cutoff);
		//#	endif
		
#if defined (__APPLE__) || (__MACH__)
		if(processors == 1){
			driver(distCalculator, seqDB, 0, numSeqs, distFile, phylipFile + "tempPhylipA", cutoff);
			appendFiles((phylipFile + "tempPhylipA"), phylipFile);
			remove((phylipFile + "tempPhylipA").c_str());	
		}	
		else if(processors == 2){
			
			int pid = fork();
			if(pid > 0){
				driver(distCalculator, seqDB, 0, (numSeqs/sqrt(2)), distFile + "tempa", phylipFile + "tempPhylipA", cutoff);	
				appendFiles((distFile+"tempa"), distFile);
				remove((distFile + "tempa").c_str());
				appendFiles((phylipFile + "tempPhylipA"), phylipFile);
				remove((phylipFile + "tempPhylipA").c_str());
			}
			else{
				driver(distCalculator, seqDB, (numSeqs/sqrt(2)), numSeqs, distFile + "tempb", phylipFile + "tempPhylipB", cutoff);	
				appendFiles((distFile+"tempb"), distFile);
				remove((distFile + "tempb").c_str());
				appendFiles((phylipFile + "tempPhylipB"), phylipFile);
				remove((phylipFile + "tempPhylipB").c_str());	
			}
			wait(NULL);
			
		}
		else if(processors == 3){
			int pid1 = fork();
			if(pid1 > 0){
				int pid2 = fork();
				if(pid2 > 0){
					driver(distCalculator, seqDB, 0, sqrt(3) * numSeqs / 3, distFile + "tempa", phylipFile + "tempPhylipA", cutoff);
					appendFiles(distFile+"tempa", distFile);
					appendFiles((phylipFile + "tempPhylipA"), phylipFile);
					remove((distFile + "tempa").c_str());
					remove((phylipFile + "tempPhylipA").c_str());	
				}
				else{
					driver(distCalculator, seqDB, sqrt(3) * numSeqs / 3, sqrt(6) * numSeqs / 3, distFile + "tempb", phylipFile + "tempPhylipB", cutoff);	
					appendFiles(distFile+"tempb", distFile);
					appendFiles((phylipFile + "tempPhylipB"), phylipFile);
					remove((distFile + "tempb").c_str());	
					remove((phylipFile + "tempPhylipB").c_str());			
				}
				wait(NULL);
			}
			else{
				driver(distCalculator, seqDB, sqrt(6) * numSeqs / 3, numSeqs, distFile + "tempc", phylipFile + "tempPhylipC", cutoff);	
				appendFiles(distFile+"tempc", distFile);
				appendFiles((phylipFile + "tempPhylipC"), phylipFile);
				remove((distFile + "tempc").c_str());	
				remove((phylipFile + "tempPhylipC").c_str());		
			}
			wait(NULL);
		}
		else if(processors == 4){
			int pid1 = fork();
			if(pid1 > 0){
				int pid2 = fork();
				if(pid2 > 0){
					driver(distCalculator, seqDB, 0, numSeqs / 2, distFile + "tempa", phylipFile + "tempPhylipA", cutoff);	
					appendFiles(distFile+"tempa", distFile);
					appendFiles((phylipFile + "tempPhylipA"), phylipFile);
					remove((distFile + "tempa").c_str());
					remove((phylipFile + "tempPhylipA").c_str());			
				}
				else{
					driver(distCalculator, seqDB, numSeqs / 2, (numSeqs/sqrt(2)), distFile + "tempb", phylipFile + "tempPhylipB", cutoff);	
					appendFiles(distFile+"tempb", distFile);
					appendFiles((phylipFile + "tempPhylipB"), phylipFile);
					remove((distFile + "tempb").c_str());
					remove((phylipFile + "tempPhylipB").c_str());				
				}
				wait(NULL);
			}
			else{
				int pid3 = fork();
				if(pid3 > 0){
					driver(distCalculator, seqDB, (numSeqs/sqrt(2)), (sqrt(3) * numSeqs / 2), distFile + "tempc", phylipFile + "tempPhylipC", cutoff);	
					appendFiles(distFile+"tempc", distFile);
					appendFiles((phylipFile + "tempPhylipC"), phylipFile);
					remove((distFile + "tempc").c_str());
					remove((phylipFile + "tempPhylipC").c_str());	
								
				}
				else{
					driver(distCalculator, seqDB, (sqrt(3) * numSeqs / 2), numSeqs, distFile + "tempd", phylipFile + "tempPhylipD", cutoff);	
					appendFiles(distFile+"tempd", distFile);
					appendFiles((phylipFile + "tempPhylipD"), phylipFile);
					remove((distFile + "tempd").c_str());	
					remove((phylipFile + "tempPhylipD").c_str());			
				}
				wait(NULL);
			}
			wait(NULL);
		}
		wait(NULL);
#elif (linux) || (__linux)
		if(processors == 1){
			driver(distCalculator, seqDB, 0, numSeqs, distFile, phylipFile + "tempPhylipA", cutoff);
			appendFiles((phylipFile + "tempPhylipA"), phylipFile);	
			remove((phylipFile + "tempPhylipA").c_str());	
		}	
		else if(processors == 2){
			
			int pid = fork();
			if(pid > 0){
				driver(distCalculator, seqDB, 0, (numSeqs/sqrt(2)), distFile + "tempa", phylipFile + "tempPhylipA", cutoff);	
				appendFiles((distFile+"tempa"), distFile);
				appendFiles((phylipFile + "tempPhylipA"), phylipFile);
				remove((distFile + "tempa").c_str());	
				remove((phylipFile + "tempPhylipA").c_str());

			}
			else{
				driver(distCalculator, seqDB, (numSeqs/sqrt(2)), numSeqs, distFile + "tempb", phylipFile + "tempPhylipB", cutoff);	
				appendFiles((distFile+"tempb"), distFile);
				appendFiles((phylipFile + "tempPhylipB"), phylipFile);
				remove((distFile + "tempb").c_str());
				remove((phylipFile + "tempPhylipB").c_str());	
			}
			wait();
			
		}
		else if(processors == 3){
			int pid1 = fork();
			if(pid1 > 0){
				int pid2 = fork();
				if(pid2 > 0){
					driver(distCalculator, seqDB, 0, (numSeqs/sqrt(2)), distFile + "tempa", phylipFile + "tempPhylipA", cutoff);	
					appendFiles((distFile+"tempa"), distFile);
					appendFiles((phylipFile + "tempPhylipA"), phylipFile);
					remove((distFile + "tempa").c_str());	
					remove((phylipFile + "tempPhylipA").c_str());

				}
				else{
					driver(distCalculator, seqDB, (numSeqs/sqrt(2)), numSeqs, distFile + "tempb", phylipFile + "tempPhylipB", cutoff);	
					appendFiles((distFile+"tempb"), distFile);
					appendFiles((phylipFile + "tempPhylipB"), phylipFile);
					remove((distFile + "tempb").c_str());
					remove((phylipFile + "tempPhylipB").c_str());	
				}
				wait();
			}
			else{
				driver(distCalculator, seqDB, sqrt(6) * numSeqs / 3, numSeqs, distFile + "tempc", phylipFile + "tempPhylipC", cutoff);	
				appendFiles(distFile+"tempc", distFile);
				appendFiles((phylipFile + "tempPhylipC"), phylipFile);
				remove((distFile + "tempc").c_str());	
				remove((phylipFile + "tempPhylipC").c_str());				
			}
			wait();
		}
		else if(processors == 4){
			int pid1 = fork();
			if(pid1 > 0){
				int pid2 = fork();
				if(pid2 > 0){
					driver(distCalculator, seqDB, 0, (numSeqs/sqrt(2)), distFile + "tempa", phylipFile + "tempPhylipA", cutoff);	
					appendFiles((distFile+"tempa"), distFile);
					appendFiles((phylipFile + "tempPhylipA"), phylipFile);
					remove((distFile + "tempa").c_str());	
					remove((phylipFile + "tempPhylipA").c_str());
				}
				else{
					driver(distCalculator, seqDB, (numSeqs/sqrt(2)), numSeqs, distFile + "tempb", phylipFile + "tempPhylipB", cutoff);	
					appendFiles((distFile+"tempb"), distFile);
					appendFiles((phylipFile + "tempPhylipB"), phylipFile);
					remove((distFile + "tempb").c_str());
					remove((phylipFile + "tempPhylipB").c_str());	
				}
				wait();
			}
			else{
				int pid3 = fork();
				if(pid3 > 0){
					driver(distCalculator, seqDB, sqrt(6) * numSeqs / 3, numSeqs, distFile + "tempc", phylipFile + "tempPhylipC", cutoff);	
					appendFiles(distFile+"tempc", distFile);
					appendFiles((phylipFile + "tempPhylipC"), phylipFile);
					remove((distFile + "tempc").c_str());	
					remove((phylipFile + "tempPhylipC").c_str());					
				}
				else{
					driver(distCalculator, seqDB, (sqrt(3) * numSeqs / 2), numSeqs, distFile + "tempd", phylipFile + "tempPhylipD", cutoff);	
					appendFiles(distFile+"tempd", distFile);
					appendFiles((phylipFile + "tempPhylipD"), phylipFile);
					remove((distFile + "tempd").c_str());	
					remove((phylipFile + "tempPhylipD").c_str());				
				}
				wait();
			}
			wait();
		}
		wait();
		
#else
		driver(distCalculator, seqDB, 0, numSeqs, distFile, phylipFile + "tempPhylipA", cutoff);
		appendFiles((phylipFile + "tempPhylipA"), phylipFile);	
		remove((phylipFile + "tempPhylipA").c_str());
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
int DistanceCommand::driver(Dist* distCalculator, SequenceDB* align, int startLine, int endLine, string dFileName, string pFilename, float cutoff){
	try {
		int startTime = time(NULL);
		
		//column file
		ofstream distFile(dFileName.c_str(), ios::trunc);
		distFile.setf(ios::fixed, ios::showpoint);
		distFile << setprecision(4);
		
		//column file
		ofstream philFile(pFilename.c_str(), ios::trunc);
		philFile.setf(ios::fixed, ios::showpoint);
		philFile << setprecision(4);
		
		for(int i=startLine;i<endLine;i++){
			
			for(int j=0;j<i;j++){
				distCalculator->calcDist(*(align->get(i)), *(align->get(j)));
				double dist = distCalculator->getDist();
				
				if(dist <= cutoff){
					distFile << align->get(i)->getName() << ' ' << align->get(j)->getName() << ' ' << dist << endl;
				}
				philFile << dist << '\t';
			}
			
			philFile << endl;
			
			if(i % 100 == 0){
				cout << i << '\t' << time(NULL) - startTime << endl;
			}
			
		}
		cout << endLine-1 << '\t' << time(NULL) - startTime << endl;
		
		//philFile.close();
		//distFile.close();
		
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