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
		phylip = globaldata->getPhylipFile();
		
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
		
		string outputFile;
		
		//doses the user want the phylip formatted file as well
		if (isTrue(phylip) == true) {
			outputFile = getRootName(globaldata->getFastaFile()) + "phylip.dist";
			remove(outputFile.c_str());
			
			//output numSeqs to phylip formatted dist file
			openOutputFile(outputFile, outFile);
			outFile << numSeqs << endl;
			outFile.close();
		}else { //user wants column format
			outputFile = getRootName(globaldata->getFastaFile()) + "dist";
			remove(outputFile.c_str());
		}
				
		//#	if defined (_WIN32)
		//figure out how to implement the fork and wait commands in windows
		//	driver(distCalculator, seqDB, 0, numSeqs, distFile, phylipFile, cutoff);
		//#	endif
		
				
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		//if you don't need to fork anything
		if(processors == 1){
			driver(distCalculator, seqDB, 0, numSeqs, outputFile + ".temp", cutoff);
			appendFiles((outputFile + ".temp"), outputFile);
			remove((outputFile + ".temp").c_str());
		}else{ //you have multiple processors
			
			for (int i = 0; i < processors; i++) {
				lines.push_back(new linePair());
				lines[i]->start = int (sqrt(float(i)/float(processors)) * numSeqs);
				lines[i]->end = int (sqrt(float(i+1)/float(processors)) * numSeqs);
			}

			cout << lines[0]->start << '\t' << lines[0]->end << endl;
			cout << lines[1]->start << '\t' << lines[1]->end << endl;

			createProcesses(outputFile); 
		
			//append and remove temp files
			for (it = processIDS.begin(); it != processIDS.end(); it++) {
				appendFiles((outputFile + toString(it->second) + ".temp"), outputFile);
				remove((outputFile + toString(it->second) + ".temp").c_str());
			}
		}
#else
		driver(distCalculator, seqDB, 0, numSeqs, outputFile + ".temp", cutoff);
		appendFiles((outputFile + ".temp"), outputFile);
		remove((outputFile + ".temp").c_str());
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
void DistanceCommand::createProcesses(string filename) {
	try {
		int process = 0;
		processIDS.clear();
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS[lines[process]->end] = pid;  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				driver(distCalculator, seqDB, lines[process]->start, lines[process]->end, filename + toString(getpid()) + ".temp", cutoff);
				exit(0);
			}else { cout << "unable to spawn the necessary processes." << endl; exit(0); }
		}
	
		//force parent to wait until all the processes are done
		for (it = processIDS.begin(); it != processIDS.end(); it++) { 
			int temp = it->second;
			wait(&temp);
		}
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the DistanceCommand class Function createProcesses. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the DistanceCommand class function createProcesses. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/**************************************************************************************************/
/////// need to fix to work with calcs and sequencedb
int DistanceCommand::driver(Dist* distCalculator, SequenceDB* align, int startLine, int endLine, string dFileName, float cutoff){
	try {

		int startTime = time(NULL);
		
		//column file
		ofstream outFile(dFileName.c_str(), ios::trunc);
		outFile.setf(ios::fixed, ios::showpoint);
		outFile << setprecision(4);
		
		for(int i=startLine;i<endLine;i++){
			
			for(int j=0;j<i;j++){
				distCalculator->calcDist(*(align->get(i)), *(align->get(j)));
				double dist = distCalculator->getDist();
				
				if(dist <= cutoff){
					if (isTrue(phylip) != true) { outFile << align->get(i)->getName() << ' ' << align->get(j)->getName() << ' ' << dist << endl; }
				}
				if (isTrue(phylip) == true) {  outFile << dist << '\t'; }
				
			}
			
			if (isTrue(phylip) == true) { outFile << endl; }
			
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