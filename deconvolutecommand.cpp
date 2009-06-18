/*
 *  deconvolute.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/21/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "deconvolutecommand.h"

/**************************************************************************************/
DeconvoluteCommand::DeconvoluteCommand(string option) {	
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string, string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//check for required parameters
			filename = validParameter.validFile(parameters, "fasta", true);
			if (filename == "not open") { abort = true; }
			else if (filename == "not found") { filename = ""; cout << "fasta is a required parameter for the unique.seqs command." << endl; abort = true;  }	
			
		}

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the DeconvoluteCommand class Function DeconvoluteCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the DeconvoluteCommand class function DeconvoluteCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
//**********************************************************************************************************************

void DeconvoluteCommand::help(){
	try {
		cout << "The unique.seqs command reads a fastafile and creates a namesfile." << "\n";
		cout << "It creates a file where the first column is the groupname and the second column is a list of sequence names who have the same sequence. " << "\n";
		cout << "If the sequence is unique the second column will just contain its name. " << "\n";
		cout << "The unique.seqs command parameter is fasta and it is required." << "\n";
		cout << "The unique.seqs command should be in the following format: " << "\n";
		cout << "unique.seqs(fasta=yourFastaFile) " << "\n";	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the DeconvoluteCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the DeconvoluteCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/**************************************************************************************/
int DeconvoluteCommand::execute() {	
	try {
		
		if (abort == true) { return 0; }
	
		//prepare filenames and open files
		outputFileName = (getRootName(filename) + "names");
		outFastafile = (getRootName(filename) + "unique.fasta");
		
		openInputFile(filename, in);
		openOutputFile(outputFileName, out);
		openOutputFile(outFastafile, outFasta);

		//constructor reads in file and store internally
		fastamap = new FastaMap();
	
		//two columns separated by tabs sequence name and then sequence
		fastamap->readFastaFile(in);
		
		//print out new names file 
		//file contains 2 columns separated by tabs.  the first column is the groupname(name of first sequence found.
		//the second column is the list of names of identical sequences separated by ','.
		fastamap->printNamesFile(out);
		fastamap->printCondensedFasta(outFasta);
		
		in.close();
		out.close();
		outFasta.close();
	
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the DeconvoluteCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the DeconvoluteCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/**************************************************************************************/
