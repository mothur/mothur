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
			else if (fastafile == "not found") { fastafile = ""; mothurOut("fasta is a required parameter for the summary.seqs command."); mothurOutEndLine(); abort = true;  }	
			
		}
	}
	catch(exception& e) {
		errorOut(e, "SeqSummaryCommand", "SeqSummaryCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void SeqSummaryCommand::help(){
	try {
		mothurOut("The summary.seqs command reads a fastafile and ....\n");
		mothurOut("The summary.seqs command parameter is fasta and it is required.\n");
		mothurOut("The summary.seqs command should be in the following format: \n");
		mothurOut("summary.seqs(fasta=yourFastaFile) \n");
		mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n\n");	
	}
	catch(exception& e) {
		errorOut(e, "SeqSummaryCommand", "help");
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
			if (current.getName() != "") {
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
			}
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
		
		//to compensate for blank sequences that would result in startPosition and endPostion equalling -1
		if (startPosition[0] == -1) {  startPosition[0] = 0;	}
		if (endPosition[0] == -1)	{  endPosition[0] = 0;		}
		
		mothurOutEndLine();
		mothurOut("\t\tStart\tEnd\tNBases\tAmbigs\tPolymer"); mothurOutEndLine();
		mothurOut("Minimum:\t" + toString(startPosition[0]) + "\t" + toString(endPosition[0]) + "\t" + toString(seqLength[0]) + "\t" + toString(ambigBases[0]) + "\t" + toString(longHomoPolymer[0])); mothurOutEndLine();
		mothurOut("2.5%-tile:\t" + toString(startPosition[ptile0_25]) + "\t" + toString(endPosition[ptile0_25]) + "\t" + toString(seqLength[ptile0_25]) + "\t" + toString(ambigBases[ptile0_25]) + "\t"+ toString(longHomoPolymer[ptile0_25])); mothurOutEndLine();
		mothurOut("25%-tile:\t" + toString(startPosition[ptile25]) + "\t" + toString(endPosition[ptile25]) + "\t" + toString(seqLength[ptile25]) + "\t" + toString(ambigBases[ptile25]) + "\t" + toString(longHomoPolymer[ptile25])); mothurOutEndLine();
		mothurOut("Median: \t" + toString(startPosition[ptile50]) + "\t" + toString(endPosition[ptile50]) + "\t" + toString(seqLength[ptile50]) + "\t" + toString(ambigBases[ptile50]) + "\t" + toString(longHomoPolymer[ptile50])); mothurOutEndLine();
		mothurOut("75%-tile:\t" + toString(startPosition[ptile75]) + "\t" + toString(endPosition[ptile75]) + "\t" + toString(seqLength[ptile75]) + "\t" + toString(ambigBases[ptile75]) + "\t" + toString(longHomoPolymer[ptile75])); mothurOutEndLine();
		mothurOut("97.5%-tile:\t" + toString(startPosition[ptile97_5]) + "\t" + toString(endPosition[ptile97_5]) + "\t" + toString(seqLength[ptile97_5]) + "\t" + toString(ambigBases[ptile97_5]) + "\t" + toString(longHomoPolymer[ptile97_5])); mothurOutEndLine();
		mothurOut("Maximum:\t" + toString(startPosition[ptile100]) + "\t" + toString(endPosition[ptile100]) + "\t" + toString(seqLength[ptile100]) + "\t" + toString(ambigBases[ptile100]) + "\t" + toString(longHomoPolymer[ptile100])); mothurOutEndLine();
		mothurOut("# of Seqs:\t" + toString(numSeqs)); mothurOutEndLine();
		
		outSummary.close();
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "SeqSummaryCommand", "execute");
		exit(1);
	}
}

//***************************************************************************************************************


