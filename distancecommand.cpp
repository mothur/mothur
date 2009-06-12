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

DistanceCommand::DistanceCommand(string option){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		Estimators.clear();
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta", "phylip", "calc", "countends", "cutoff", "processors"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			parser = new OptionParser();
			parser->parse(option, parameters);  delete parser;
			
			ValidParameters* validParameter = new ValidParameters();
		
			//check to make sure all parameters are valid for command
			for (it2 = parameters.begin(); it2 != parameters.end(); it2++) { 
				if (validParameter->isValidParameter(it2->first, myArray, it2->second) != true) {  abort = true;  }
			}
			
			//check for required parameters
			fastafile = validParameter->validFile(parameters, "fasta", true);
			if (fastafile == "not found") { cout << "fasta is a required parameter for the dist.seqs command." << endl; abort = true; }
			else if (fastafile == "not open") { abort = true; }	
			else { 
				globaldata->setFastaFile(fastafile);
				openInputFile(fastafile, in);
			}
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			calc = validParameter->validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "onegap";  }
			else { 
				 if (calc == "default")  {  calc = "onegap";  }
			}
			splitAtDash(calc, Estimators);

			string temp;
			temp = validParameter->validFile(parameters, "countends", false);			if (temp == "not found") { temp = "T"; }
			convert(temp, countends); 
			
			temp = validParameter->validFile(parameters, "cutoff", false);				if (temp == "not found") { temp = "1.0"; }
			convert(temp, cutoff); 
			
			temp = validParameter->validFile(parameters, "processors", false);			if (temp == "not found") { temp = "1"; }
			convert(temp, processors); 
			
			phylip = validParameter->validFile(parameters, "phylip", false);			if (phylip == "not found") { phylip = "F"; }
	
			delete validParameter;
			
			validCalculator = new ValidCalculators();
			
			int i;
			if (isTrue(countends) == true) {
				for (i=0; i<Estimators.size(); i++) {
					if (validCalculator->isValidCalculator("distance", Estimators[i]) == true) { 
						if (Estimators[i] == "nogaps") { 
							distCalculator = new ignoreGaps();
						}else if (Estimators[i] == "eachgap") { 
							distCalculator = new eachGapDist();	
						}else if (Estimators[i] == "onegap") {
						distCalculator = new oneGapDist();					}
					}
				}
			}else {
				for (i=0; i<Estimators.size(); i++) {
					if (validCalculator->isValidCalculator("distance", Estimators[i]) == true) { 
						if (Estimators[i] == "nogaps") { 
							distCalculator = new ignoreGaps();	
						}else if (Estimators[i] == "eachgap") { 
							distCalculator = new eachGapIgnoreTermGapDist();
						}else if (Estimators[i] == "onegap") { 
							distCalculator = new oneGapIgnoreTermGapDist();	
						}
					}
				}
			}

			delete validCalculator;
		}
				
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

void DistanceCommand::help(){
	try {
		cout << "The dist.seqs command reads a file containing sequences and creates a distance file." << "\n";
		cout << "The dist.seqs command parameters are fasta, calc, countends, cutoff and processors.  " << "\n";
		cout << "The fasta parameter is required." << "\n";
		cout << "The calc parameter allows you to specify the method of calculating the distances.  Your options are: nogaps, onegap or eachgap. The default is onegap." << "\n";
		cout << "The countends parameter allows you to specify whether to include terminal gaps in distance.  Your options are: T or F. The default is T." << "\n";
		cout << "The cutoff parameter allows you to specify maximum distance to keep. The default is 1.0." << "\n";
		cout << "The processors parameter allows you to specify number of processors to use.  The default is 1." << "\n";
		cout << "The dist.seqs command should be in the following format: " << "\n";
		cout << "dist.seqs(fasta=yourFastaFile, calc=yourCalc, countends=yourEnds, cutoff= yourCutOff, processors=yourProcessors) " << "\n";
		cout << "Example dist.seqs(fasta=amazon.fasta, calc=eachgap, countends=F, cutoff= 2.0, processors=3)." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. calc), '=' and parameters (i.e.yourCalc)." << "\n" << "\n";
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the DistanceCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the DistanceCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}


//**********************************************************************************************************************

int DistanceCommand::execute(){
	try {
		
		if (abort == true) { return 0; }
		
		//reads fasta file and fills sequenceDB
		seqDB = new SequenceDB(in); 
				
		int numSeqs = seqDB->getNumSeqs();
		cutoff += 0.005;
		
		string outputFile;
		
		//doses the user want the phylip formatted file as well
		if (isTrue(phylip) == true) {
			outputFile = getRootName(fastafile) + "phylip.dist";
			remove(outputFile.c_str());
			
			//output numSeqs to phylip formatted dist file
			openOutputFile(outputFile, outFile);
			outFile << numSeqs << endl;
			outFile.close();
		}else { //user wants column format
			outputFile = getRootName(fastafile) + "dist";
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