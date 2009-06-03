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
	
	if(globaldata->getFastaFile() != "")		{	readSeqs =  new ReadFasta(globaldata->inputFileName);	}
	else if(globaldata->getNexusFile() != "")	{	readSeqs = new ReadNexus(globaldata->inputFileName);	}
	else if(globaldata->getClustalFile() != "") {	readSeqs = new ReadClustal(globaldata->inputFileName);	}
	else if(globaldata->getPhylipFile() != "")	{	readSeqs = new ReadPhylip(globaldata->inputFileName);	}
	
	readSeqs->read();
	db = readSeqs->getDB();
	numSeqs = db->size();
	
	alignmentLength = db->get(0).getAlignLength();

	filter = string(alignmentLength, '1');
}

/**************************************************************************************/

void FilterSeqsCommand::doHard() {
	
	string hardName = globaldata->getHard();
	string hardFilter = "";
		
	ifstream fileHandle;
	openInputFile(hardName, fileHandle);
	
	fileHandle >> hardFilter;
	
	if(hardFilter.length() != filter.length()){
		cout << "The hard filter is not the same length as the alignment: Hard filter will not be applied." << endl;
	}
	else{
		filter = hardFilter;
	}
	
}

/**************************************************************************************/

void FilterSeqsCommand::doTrump() {

	char trump = globaldata->getTrump()[0];
	
	for(int i = 0; i < numSeqs; i++) {
		string curAligned = db->get(i).getAligned();;

		for(int j = 0; j < alignmentLength; j++) {
			if(curAligned[j] == trump){
				filter[j] = '0';
			}
		}
	}

}

/**************************************************************************************/

void FilterSeqsCommand::doVertical() {

	vector<int> counts(alignmentLength, 0);
	
	for(int i = 0; i < numSeqs; i++) {
		string curAligned = db->get(i).getAligned();;
		
		for(int j = 0; j < alignmentLength; j++) {
			if(curAligned[j] == '-' || curAligned[j] == '.'){
				counts[j]++;
			}
		}
	}
	for(int i=0;i<alignmentLength;i++){
		if(counts[i] == numSeqs)	{	filter[i] = '0';		}
	}
}

/**************************************************************************************/

void FilterSeqsCommand::doSoft() {

	int softThreshold = numSeqs * (float)atoi(globaldata->getSoft().c_str()) / 100.0;

	vector<int> a(alignmentLength, 0);
	vector<int> t(alignmentLength, 0);
	vector<int> g(alignmentLength, 0);
	vector<int> c(alignmentLength, 0);
	vector<int> x(alignmentLength, 0);
	
	for(int i=0;i<numSeqs;i++){
		string curAligned = db->get(i).getAligned();;

		for(int j=0;j<alignmentLength;j++){
			if(toupper(curAligned[j]) == 'A')										{	a[j]++;	}
			else if(toupper(curAligned[j]) == 'T' || toupper(curAligned[i]) == 'U')	{	t[j]++;	}
			else if(toupper(curAligned[j]) == 'G')									{	g[j]++;	}
			else if(toupper(curAligned[j]) == 'C')									{	c[j]++;	}
		}
	}

	for(int i=0;i<alignmentLength;i++){
		if(a[i] < softThreshold && t[i] < softThreshold && g[i] < softThreshold && c[i] < softThreshold){
			filter[i] = '0';			
		}
	}
}

/**************************************************************************************/

int FilterSeqsCommand::execute() {	
	try {
						
		if(globaldata->getHard().compare("") != 0)		{	doHard();		}	//	has to be applied first!
		if(globaldata->getTrump().compare("") != 0)		{	doTrump();		}
		if(globaldata->getVertical() == "T")			{	doVertical();	}
		if(globaldata->getSoft().compare("") != 0)		{	doSoft();		}

		ofstream outfile;
		string filterFile = getRootName(globaldata->inputFileName) + "filter";
		openOutputFile(filterFile, outfile);

		outfile << filter << endl;
		outfile.close();
		
		string filteredFasta = getRootName(globaldata->inputFileName) + "filter.fasta";
		openOutputFile(filteredFasta, outfile);

		for(int i=0;i<numSeqs;i++){
			string curAligned = db->get(i).getAligned();
			outfile << '>' << db->get(i).getName() << endl;
			for(int j=0;j<alignmentLength;j++){
				if(filter[j] == '1'){
					outfile << curAligned[j];
				}
			}
			outfile << endl;
		}
		outfile.close();
		
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
