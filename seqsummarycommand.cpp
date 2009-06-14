/*
 *  seqcoordcommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 5/30/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "seqsummarycommand.h"
#include "sequence.hpp"

//***************************************************************************************************************

SeqSummaryCommand::SeqSummaryCommand(string option){
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
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") { fastafile = ""; cout << "fasta is a required parameter for the summary.seqs command." << endl; abort = true;  }	
			
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SeqSummaryCommand class Function SeqSummaryCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SeqSummaryCommand class function SeqSummaryCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************

void SeqSummaryCommand::help(){
	try {
		cout << "The summary.seqs command reads a fastafile and ...." << "\n";
		cout << "The summary.seqs command parameter is fasta and it is required." << "\n";
		cout << "The summary.seqs command should be in the following format: " << "\n";
		cout << "summary.seqs(fasta=yourFastaFile) " << "\n";
		cout << "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile)." << "\n" << "\n";	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SeqSummaryCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SeqSummaryCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//***************************************************************************************************************

SeqSummaryCommand::~SeqSummaryCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************

int SeqSummaryCommand::execute(){
	try{
		
		if (abort == true) { return 0; }
		
		ifstream inFASTA;
		openInputFile(fastafile, inFASTA);
		int numSeqs = 0;

		ofstream outSummary;
		string summaryFile = fastafile + ".summary";
		openOutputFile(summaryFile, outSummary);
		
		vector<int> startPosition;
		vector<int> endPosition;
		vector<int> seqLength;
		vector<int> ambigBases;
		vector<int> longHomoPolymer;
		
		outSummary << "seqname\tstart\tend\tnbases\tambigs\tpolymer" << endl;			

		while(!inFASTA.eof()){
			Sequence current(inFASTA);
			startPosition.push_back(current.getStartPos());
			endPosition.push_back(current.getEndPos());
			seqLength.push_back(current.getNumBases());
			ambigBases.push_back(current.getAmbigBases());
			longHomoPolymer.push_back(current.getLongHomoPolymer());

			outSummary << current.getName() << '\t';
			outSummary << current.getStartPos() << '\t' << current.getEndPos() << '\t';
			outSummary << current.getNumBases() << '\t' << current.getAmbigBases() << '\t';
			outSummary << current.getLongHomoPolymer() << endl;
			
			numSeqs++;
			gobble(inFASTA);
		}
		inFASTA.close();
		
		sort(startPosition.begin(), startPosition.end());
		sort(endPosition.begin(), endPosition.end());
		sort(seqLength.begin(), seqLength.end());
		sort(ambigBases.begin(), ambigBases.end());
		sort(longHomoPolymer.begin(), longHomoPolymer.end());
		
		int ptile0_25	= int(numSeqs * 0.025);
		int ptile25		= int(numSeqs * 0.250);
		int ptile50		= int(numSeqs * 0.500);
		int ptile75		= int(numSeqs * 0.750);
		int ptile97_5	= int(numSeqs * 0.975);
		int ptile100	= numSeqs - 1;
		
		cout << endl;
		cout << "\t\tStart\tEnd\tNBases\tAmbigs\tPolymer" << endl;
		cout << "Minimum:\t" << startPosition[0] << '\t' << endPosition[0] << '\t' << seqLength[0] << '\t' << ambigBases[0] << '\t' << longHomoPolymer[0] << endl;
		cout << "2.5%-tile:\t" << startPosition[ptile0_25] << '\t' << endPosition[ptile0_25] << '\t' << seqLength[ptile0_25] << '\t' << ambigBases[ptile0_25] << '\t' << longHomoPolymer[ptile0_25] << endl;
		cout << "25%-tile:\t" << startPosition[ptile25] << '\t' << endPosition[ptile25] << '\t' << seqLength[ptile25] << '\t' << ambigBases[ptile25] << '\t' << longHomoPolymer[ptile25] << endl;
		cout << "Median: \t" << startPosition[ptile50] << '\t' << endPosition[ptile50] << '\t' << seqLength[ptile50] << '\t' << ambigBases[ptile50] << '\t' << longHomoPolymer[ptile50] << endl;
		cout << "75%-tile:\t" << startPosition[ptile75] << '\t' << endPosition[ptile75] << '\t' << seqLength[ptile75] << '\t' << ambigBases[ptile75] << '\t' << longHomoPolymer[ptile75] << endl;
		cout << "97.5%-tile:\t" << startPosition[ptile97_5] << '\t' << endPosition[ptile97_5] << '\t' << seqLength[ptile97_5] << '\t' << ambigBases[ptile97_5] << '\t' << longHomoPolymer[ptile97_5] << endl;
		cout << "Maximum:\t" << startPosition[ptile100] << '\t' << endPosition[ptile100] << '\t' << seqLength[ptile100] << '\t' << ambigBases[ptile100] << '\t' << longHomoPolymer[ptile100] << endl;
		cout << "# of Seqs:\t" << numSeqs << endl;
		
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SeqSummaryCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SeqSummaryCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	
}

//***************************************************************************************************************


