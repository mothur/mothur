/*
 *  seqcoordcommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 5/30/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "seqsummarycommand.h"

//***************************************************************************************************************

SeqSummaryCommand::SeqSummaryCommand(){
	try {
		globaldata = GlobalData::getInstance();
		
		if(globaldata->getFastaFile() != "")		{	readSeqs = new ReadFasta(globaldata->inputFileName);	}
		else if(globaldata->getNexusFile() != "")	{	readSeqs = new ReadNexus(globaldata->inputFileName);	}
		else if(globaldata->getClustalFile() != "") {	readSeqs = new ReadClustal(globaldata->inputFileName);	}
		else if(globaldata->getPhylipFile() != "")	{	readSeqs = new ReadPhylip(globaldata->inputFileName);	}
		
		readSeqs->read();
		db = readSeqs->getDB();
		numSeqs = db->size();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SeqCoordCommand class Function SeqCoordCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SeqCoordCommand class function SeqCoordCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//***************************************************************************************************************

SeqSummaryCommand::~SeqSummaryCommand(){
}

//***************************************************************************************************************

int SeqSummaryCommand::execute(){
	try{
		
		ofstream outfile;
		string summaryFile = getRootName(globaldata->inputFileName) + "summary";
		openOutputFile(summaryFile, outfile);

		vector<int> startPosition(numSeqs, 0);
		vector<int> endPosition(numSeqs, 0);
		vector<int> seqLength(numSeqs, 0);
		vector<int> ambigBases(numSeqs, 0);
		vector<int> longHomoPolymer(numSeqs, 0);
		
		if(db->get(0).getIsAligned() == 1){
			outfile << "seqname\tstart\tend\tlength\tambiguities\tlonghomopolymer" << endl;			
			for(int i = 0; i < numSeqs; i++) {
				Sequence current = db->get(i);
				startPosition[i] = current.getStartPos();
				endPosition[i] = current.getEndPos();
				seqLength[i] = current.getNumBases();
				ambigBases[i] = current.getAmbigBases();
				longHomoPolymer[i] = current.getLongHomoPolymer();
				outfile << current.getName() << '\t' << startPosition[i] << '\t' << endPosition[i] << '\t' << seqLength[i] << '\t' << ambigBases[i] << '\t' << longHomoPolymer[i] << endl;
			}
		}
		else{
			outfile << "seqname\tlength\tambiguities\tlonghomopolymer" << endl;
			for(int i=0;i<numSeqs;i++){
				Sequence current = db->get(i);
				seqLength[i] = current.getNumBases();
				ambigBases[i] = current.getAmbigBases();
				longHomoPolymer[i] = current.getLongHomoPolymer();
				outfile << current.getName() << '\t' << seqLength[i] << '\t' << ambigBases[i] << '\t' << longHomoPolymer[i] << endl;
			}
		}
		
		sort(seqLength.begin(), seqLength.end());
		sort(ambigBases.begin(), ambigBases.end());
		sort(longHomoPolymer.begin(), longHomoPolymer.end());
		
		int median			= int(numSeqs * 0.500);
		int lowestPtile		= int(numSeqs * 0.025);
		int lowPtile		= int(numSeqs * 0.250);
		int highPtile		= int(numSeqs * 0.750);
		int highestPtile	= int(numSeqs * 0.975);
		int max				= numSeqs - 1;
		
		cout << endl;
		if(db->get(0).getIsAligned() == 1){
			sort(startPosition.begin(), startPosition.end());
			sort(endPosition.begin(), endPosition.end());
					
			cout << "\t\tStart\tEnd\tLength\tN's\tPolymer" << endl;
			cout << "Minimum:\t" << startPosition[0] << '\t' << endPosition[0] << '\t' << seqLength[0] << '\t' << ambigBases[0] << '\t' << longHomoPolymer[0] << endl;
			cout << "2.5%-tile:\t" << startPosition[lowestPtile] << '\t' << endPosition[lowestPtile] << '\t' << seqLength[lowestPtile] << '\t' << ambigBases[lowestPtile] << '\t' << longHomoPolymer[lowestPtile] << endl;
			cout << "25%-tile:\t" << startPosition[lowPtile] << '\t' << endPosition[lowPtile] << '\t' << seqLength[lowPtile] << '\t' << ambigBases[lowPtile] << '\t' << longHomoPolymer[lowPtile] << endl;
			cout << "Median: \t" << startPosition[median] << '\t' << endPosition[median] << '\t' << seqLength[median] << '\t' << ambigBases[median] << '\t' << longHomoPolymer[median] << endl;
			cout << "75%-tile:\t" << startPosition[highPtile] << '\t' << endPosition[highPtile] << '\t' << seqLength[highPtile] << '\t' << ambigBases[highPtile] << '\t' << longHomoPolymer[highPtile] << endl;
			cout << "97.5%-tile:\t" << startPosition[highestPtile] << '\t' << endPosition[highestPtile] << '\t' << seqLength[highestPtile] << '\t' << ambigBases[highestPtile] << '\t' << longHomoPolymer[highestPtile] << endl;
			cout << "Maximum:\t" << startPosition[max] << '\t' << endPosition[max] << '\t' << seqLength[max] << '\t' << ambigBases[max] << '\t' << longHomoPolymer[max] << endl;
		}
		else{
			cout << "\t\tLength\tN's\tPolymer" << endl;
			cout << "Minimum:\t" << seqLength[0] << '\t' << ambigBases[0] << '\t' << longHomoPolymer[0] << endl;
			cout << "2.5%-tile:\t" << seqLength[lowestPtile] << '\t' << ambigBases[lowestPtile] << '\t' << longHomoPolymer[lowestPtile] << endl;
			cout << "25%-tile:\t" << seqLength[lowPtile] << '\t' << ambigBases[lowPtile] << '\t' << longHomoPolymer[lowPtile] << endl;
			cout << "Median: \t" << seqLength[median] << '\t' << ambigBases[median] << '\t' << longHomoPolymer[median] << endl;
			cout << "75%-tile:\t"<< seqLength[highPtile] << '\t' << ambigBases[highPtile] << '\t' << longHomoPolymer[highPtile] << endl;
			cout << "97.5%-tile:\t"<< seqLength[highestPtile] << '\t' << ambigBases[highestPtile] << '\t' << longHomoPolymer[highestPtile] << endl;
			cout << "Maximum:\t" << seqLength[max] << '\t' << ambigBases[max] << '\t' << longHomoPolymer[max] << endl;
		}
		cout << "# of Seqs:\t" << numSeqs << endl;
		
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the FilterSeqsCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the FilterSeqsCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	
}

//***************************************************************************************************************


