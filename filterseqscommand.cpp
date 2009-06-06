/*
 *  filterseqscommand.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/4/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "filterseqscommand.h"

/**************************************************************************************/

FilterSeqsCommand::FilterSeqsCommand(){

	globaldata = GlobalData::getInstance();
	
	if(globaldata->getFastaFile() == "")		{	cout << "You must enter a fasta formatted file" << endl;	}
	trump = globaldata->getTrump()[0];
	numSeqs = 0;

}

/**************************************************************************************/

void FilterSeqsCommand::doHard() {
	
	string hardName = globaldata->getHard();
	string hardFilter = "";
		
	ifstream fileHandle;
	openInputFile(hardName, fileHandle);
	
	fileHandle >> filter;

}

/**************************************************************************************/

void FilterSeqsCommand::doTrump(Sequence seq) {
	
	string curAligned = seq.getAligned();
	
	for(int j = 0; j < alignmentLength; j++) {
		if(curAligned[j] == trump){
			filter[j] = '0';
		}
	}

}

/**************************************************************************************/

void FilterSeqsCommand::doVertical() {

	for(int i=0;i<alignmentLength;i++){
		if(gap[i] == numSeqs)	{	filter[i] = '0';	}
	}
	
}

/**************************************************************************************/

void FilterSeqsCommand::doSoft() {
	
	int threshold = int (soft * numSeqs);
	bool keep = 0;
	
	for(int i=0;i<alignmentLength;i++){
		if(a[i] >= threshold)		{	keep = 1;	}
		else if(t[i] >= threshold)	{	keep = 1;	}
		else if(g[i] >= threshold)	{	keep = 1;	}
		else if(c[i] >= threshold)	{	keep = 1;	}
		
		if(keep == 0)	{	filter[i] = 0;		}
	}
}

/**************************************************************************************/

void FilterSeqsCommand::getFreqs(Sequence seq) {

	string curAligned = seq.getAligned();;
	
	for(int j=0;j<alignmentLength;j++){
		if(toupper(curAligned[j]) == 'A')										{	a[j]++;		}
		else if(toupper(curAligned[j]) == 'T' || toupper(curAligned[j]) == 'U')	{	t[j]++;		}
		else if(toupper(curAligned[j]) == 'G')									{	g[j]++;		}
		else if(toupper(curAligned[j]) == 'C')									{	c[j]++;		}
		else if(curAligned[j] == '-' || curAligned[j] == '.')					{	gap[j]++;	}
	}
	
}

/**************************************************************************************/

int FilterSeqsCommand::execute() {	
	try {
		ifstream inFASTA;
		openInputFile(globaldata->getFastaFile(), inFASTA);
		
		Sequence testSeq(inFASTA);
		alignmentLength = testSeq.getAlignLength();
		inFASTA.seekg(0);
		
		if(globaldata->getSoft() != "" || isTrue(globaldata->getVertical())){
			a.assign(alignmentLength, 0);
			t.assign(alignmentLength, 0);
			g.assign(alignmentLength, 0);
			c.assign(alignmentLength, 0);
			gap.assign(alignmentLength, 0);
		}
		if(globaldata->getSoft() != ""){
			soft = (float)atoi(globaldata->getSoft().c_str()) / 100.0;
		}
		
		if(globaldata->getHard().compare("") != 0)	{	doHard();								}
		else										{	filter = string(alignmentLength, '1');	}

		if(globaldata->getTrump().compare("") != 0 || isTrue(globaldata->getVertical()) || globaldata->getSoft().compare("") != 0){
		
			while(!inFASTA.eof()){
				Sequence seq(inFASTA);
				if(globaldata->getTrump().compare("") != 0)	{	doTrump(seq);		}
				if(isTrue(globaldata->getVertical()) || globaldata->getSoft().compare("") != 0){	getFreqs(seq);	}
				numSeqs++;
				cout.flush();
			}
		
		}
		inFASTA.close();
		
		if(isTrue(globaldata->getVertical()) == 1)	{	doVertical();	}
		if(globaldata->getSoft().compare("") != 0)	{	doSoft();		}			

		ofstream outFilter;
		string filterFile = getRootName(globaldata->inputFileName) + "filter";
		openOutputFile(filterFile, outFilter);
		outFilter << filter << endl;
		outFilter.close();
		

		openInputFile(globaldata->getFastaFile(), inFASTA);
		string filteredFasta = getRootName(globaldata->inputFileName) + "filter.fasta";
		ofstream outFASTA;
		openOutputFile(filteredFasta, outFASTA);

		numSeqs = 0;
		while(!inFASTA.eof()){
			Sequence seq(inFASTA);
			string align = seq.getAligned();
			string filterSeq = "";
	
			for(int j=0;j<alignmentLength;j++){
				if(filter[j] == '1'){
					filterSeq += align[j];
				}
			}

			outFASTA << '>' << seq.getName() << endl << filterSeq << endl;
			numSeqs++;
			gobble(inFASTA);
		}
		outFASTA.close();
		inFASTA.close();
		
		
		int filteredLength = 0;
		for(int i=0;i<alignmentLength;i++){
			if(filter[i] == '1'){	filteredLength++;	}
		}
		
		cout << endl;
		cout << "Length of filtered alignment: " << filteredLength << endl;
		cout << "Number of columns removed: " << alignmentLength-filteredLength << endl;
		cout << "Length of the original alignment: " << alignmentLength << endl;
		cout << "Number of sequences used to construct filter: " << numSeqs << endl;
		
		globaldata->clear();
		
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

/**************************************************************************************/
