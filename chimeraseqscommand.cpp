/*
 *  chimeraseqscommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/29/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "chimeraseqscommand.h"

//***************************************************************************************************************

ChimeraSeqsCommand::ChimeraSeqsCommand(string option){
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta", "filter", "correction", "processors", "method" };
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
			else if (fastafile == "not found") { fastafile = ""; mothurOut("fasta is a required parameter for the chimera.seqs command."); mothurOutEndLine(); abort = true;  }	
			
			string temp;
			temp = validParameter.validFile(parameters, "filter", false);			if (temp == "not found") { temp = "F"; }
			filter = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "correction", false);		if (temp == "not found") { temp = "T"; }
			correction = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "processors", true);		if (temp == "not found") { temp = "1"; }
			convert(temp, processors);
			
			method = validParameter.validFile(parameters, "method", false);		if (method == "not found") { method = "bellerophon"; }

		}
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSeqsCommand", "ChimeraSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void ChimeraSeqsCommand::help(){
	try {
		mothurOut("The chimera.seqs command reads a fastafile and creates a sorted priority score list of potentially chimeric sequences (ideally, the sequences should already be aligned).\n");
		mothurOut("The chimera.seqs command parameters are fasta, filter, correction, processors and method.  fasta is required.\n");
		mothurOut("The filter parameter allows you to specify if you would like to apply a 50% soft filter.  The default is false. \n");
		mothurOut("The correction parameter allows you to .....  The default is true. \n");
		mothurOut("The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n");
		mothurOut("The method parameter allows you to specify the method for finding chimeric sequences.  The default is bellerophon. \n");
		mothurOut("The chimera.seqs command should be in the following format: \n");
		mothurOut("chimera.seqs(fasta=yourFastaFile, filter=yourFilter, correction=yourCorrection, processors=yourProcessors, method=bellerophon) \n");
		mothurOut("Example: chimera.seqs(fasta=AD.align, filter=True, correction=true, processors=2, method=yourMethod) \n");
		mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n\n");	
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSeqsCommand", "help");
		exit(1);
	}
}

//***************************************************************************************************************

ChimeraSeqsCommand::~ChimeraSeqsCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************

int ChimeraSeqsCommand::execute(){
	try{
		
		if (abort == true) { return 0; }
		
		//do soft filter
		if (filter)  {
			string optionString = "fasta=" + fastafile + ", soft=50.0, vertical=F";
			filterSeqs = new FilterSeqsCommand(optionString);
			filterSeqs->execute();
			delete filterSeqs;
			
			//reset fastafile to filtered file
			fastafile = getRootName(fastafile) + "filter.fasta";
		}
		
		//read in sequences
		readSeqs();
		
		//int numSeqs = seqs.size();
		
		//find average midpoint of seqs
		midpoint = findAverageMidPoint();
				
		//this should be parallelized
		//generatePreferences();
		
				
		//output results to screen						
		mothurOutEndLine();
		mothurOut("\t\t"); mothurOutEndLine();
		//mothurOut("Minimum:\t" + toString(startPosition[0]) + "\t" + toString(endPosition[0]) + "\t" + toString(seqLength[0]) + "\t" + toString(ambigBases[0]) + "\t" + toString(longHomoPolymer[0])); mothurOutEndLine();
		//mothurOut("2.5%-tile:\t" + toString(startPosition[ptile0_25]) + "\t" + toString(endPosition[ptile0_25]) + "\t" + toString(seqLength[ptile0_25]) + "\t" + toString(ambigBases[ptile0_25]) + "\t"+ toString(longHomoPolymer[ptile0_25])); mothurOutEndLine();
		//mothurOut("25%-tile:\t" + toString(startPosition[ptile25]) + "\t" + toString(endPosition[ptile25]) + "\t" + toString(seqLength[ptile25]) + "\t" + toString(ambigBases[ptile25]) + "\t" + toString(longHomoPolymer[ptile25])); mothurOutEndLine();
		//mothurOut("Median: \t" + toString(startPosition[ptile50]) + "\t" + toString(endPosition[ptile50]) + "\t" + toString(seqLength[ptile50]) + "\t" + toString(ambigBases[ptile50]) + "\t" + toString(longHomoPolymer[ptile50])); mothurOutEndLine();
		//mothurOut("75%-tile:\t" + toString(startPosition[ptile75]) + "\t" + toString(endPosition[ptile75]) + "\t" + toString(seqLength[ptile75]) + "\t" + toString(ambigBases[ptile75]) + "\t" + toString(longHomoPolymer[ptile75])); mothurOutEndLine();
		//mothurOut("97.5%-tile:\t" + toString(startPosition[ptile97_5]) + "\t" + toString(endPosition[ptile97_5]) + "\t" + toString(seqLength[ptile97_5]) + "\t" + toString(ambigBases[ptile97_5]) + "\t" + toString(longHomoPolymer[ptile97_5])); mothurOutEndLine();
		//mothurOut("Maximum:\t" + toString(startPosition[ptile100]) + "\t" + toString(endPosition[ptile100]) + "\t" + toString(seqLength[ptile100]) + "\t" + toString(ambigBases[ptile100]) + "\t" + toString(longHomoPolymer[ptile100])); mothurOutEndLine();
		//mothurOut("# of Seqs:\t" + toString(numSeqs)); mothurOutEndLine();
		
		//outSummary.close();
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSeqsCommand", "execute");
		exit(1);
	}
}

//***************************************************************************************************************
void ChimeraSeqsCommand::readSeqs(){
	try {
		ifstream inFASTA;
		openInputFile(fastafile, inFASTA);
		
		//read in seqs and store in vector
		while(!inFASTA.eof()){
			Sequence current(inFASTA);
			
			seqs.push_back(current);
			
			gobble(inFASTA);
		}
		inFASTA.close();

	}
	catch(exception& e) {
		errorOut(e, "ChimeraSeqsCommand", "readSeqs");
		exit(1);
	}
}


//***************************************************************************************************************
int ChimeraSeqsCommand::findAverageMidPoint(){
	try {
		int totalMids = 0;
		int averageMid = 0;
		
		//loop through the seqs and find midpoint
		for (int i = 0; i < seqs.size(); i++) {
			
			//get unaligned sequence
			seqs[i].setUnaligned(seqs[i].getUnaligned());  //if you read an aligned file the unaligned is really aligned, so we need to make sure its unaligned
			
			string unaligned = seqs[i].getUnaligned();
			string aligned = seqs[i].getAligned();
			
			//find midpoint of this seq
			int count = 0;
			int thismid = 0;
			for (int j = 0; j < aligned.length(); j++) {
				
				thismid++;
				
				//if you are part of the unaligned sequence increment
				if (isalpha(aligned[j])) {  count++;  }
				
				//if you have reached the halfway point stop
				if (count >= (unaligned.length() / 2)) { break; }
			}
			
			//add this mid to total
			totalMids += thismid;
		
		}
		
		averageMid = (totalMids / seqs.size());
		
		return averageMid; 
	
	
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSeqsCommand", "findAverageMidPoint");
		exit(1);
	}
}

/***************************************************************************************************************
int ChimeraSeqsCommand::createSparseMatrix(int startLine, int endLine, SparseMatrix* sparse){
	try {

		for(int i=startLine; i<endLine; i++){
			
			for(int j=0;j<i;j++){
			
				distCalculator->calcDist(seqs.get(i), seqs.get(j));
				double dist = distCalculator->getDist();
				
				
				
			}
			
	
		return 1;
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSeqsCommand", "createSparseMatrix");
		exit(1);
	}
}
/**************************************************************************************************/

