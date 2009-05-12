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
void FilterSeqsCommand::doTrump() {
	trump = globaldata->getTrump();
	for(int i = 0; i < db->size(); i++) {
		Sequence cur = db->get(i);
		string curAligned = cur.getUnaligned();
		for(int j = 0; j < curAligned.length(); j++) {
			string curChar = curAligned.substr(j, 1);
			if(curChar.compare(trump) == 0) 
				columnsToRemove[j] = true;
		}
	}
}

/**************************************************************************************/
void FilterSeqsCommand::doSoft() {
	soft = atoi(globaldata->getSoft().c_str());
	vector<vector<int> > columnSymbolSums;
	vector<vector<string> > columnSymbols;
	for(int i = 0; i < db->get(0).getLength(); i++) {
		vector<string> symbols;
		vector<int> sums;
		columnSymbols.push_back(symbols);
		columnSymbolSums.push_back(sums);
	}
	
	for(int i = 0; i < db->size(); i++) {
		Sequence cur = db->get(i);
		string curAligned = cur.getUnaligned();
		
		for(int j = 0; j < curAligned.length(); j++) {
			string curChar = curAligned.substr(j, 1);
			vector<string> curColumnSymbols = columnSymbols[j];
			bool newSymbol = true;
			
			for(int k = 0; k < curColumnSymbols.size(); k++) 
				if(curChar.compare(curColumnSymbols[k]) == 0) {
					newSymbol = false;
					columnSymbolSums[j][k]++;
				}
			
			if(newSymbol) {
				columnSymbols[j].push_back(curChar);
				columnSymbolSums[j].push_back(1);
			}
		}
	}
	
	
	for(int i = 0; i < columnSymbolSums.size(); i++) {
		int totalSum = 0;
		int max = 0;
		vector<int> curColumnSymbols = columnSymbolSums[i];
		
		for(int j = 0; j < curColumnSymbols.size(); j++) {
			int curSum = curColumnSymbols[j];
			//cout << columnSymbols[i][j] << ": " << curSum << "\n";
			if(curSum > max)
				max = curSum;
			totalSum += curSum;
		}
		//cout << "\n";
		
		if((double)max/(double)totalSum * 100 < soft)
			columnsToRemove[i] = true;
	}
}

/**************************************************************************************/
void FilterSeqsCommand::doFilter() {
	filter = globaldata->getFilter();
	ifstream filehandle;
	openInputFile(filter, filehandle);
	
	char c;
	int count = 0;
	while(!filehandle.eof()) {
		c = filehandle.get();
		if(c == '0') 
			columnsToRemove[count] = true;
		count++;
	}
}

/**************************************************************************************/
int FilterSeqsCommand::execute() {	
	try {
		globaldata = GlobalData::getInstance();
		filename = globaldata->inputFileName;
		
		if(globaldata->getFastaFile() != "") {
			readSeqs =  new ReadFasta(filename); }
		else if(globaldata->getNexusFile() != "") {
			readSeqs = new ReadNexus(filename); }
		else if(globaldata->getClustalFile() != "") {
			readSeqs = new ReadClustal(filename); }
		else if(globaldata->getPhylipFile() != "") {
			readSeqs = new ReadPhylip(filename); }
			
		readSeqs->read();
		db = readSeqs->getDB();
		
		//for(int i = 0; i < db->size(); i++) {
//			cout << db->get(i).getLength() << "\n" << db->get(i).getName() << ": " << db->get(i).getUnaligned() << "\n\n";
//		}

		for(int i = 0; i < db->get(0).getLength(); i++) 
			columnsToRemove.push_back(false);
		
				
		if(globaldata->getTrump().compare("") != 0) 
			doTrump();
		else if(globaldata->getSoft().compare("") != 0)
			doSoft();
		else if(globaldata->getFilter().compare("") != 0) 
			doFilter();
		
		//for(int i = 0; i < columnsToRemove.size(); i++)
//		{
//			cout << "Remove Column " << i << " = ";
//			if(columnsToRemove[i])
//				cout << "true\n";
//			else
//				cout << "false\n";
//		}


		//Creating the new SequenceDB 
		SequenceDB newDB;
		for(int i = 0; i < db->size(); i++) {
			Sequence curSeq = db->get(i);
			string curAligned = curSeq.getUnaligned();
			string curName = curSeq.getName();
			string newAligned = "";
			for(int j = 0; j < curAligned.length(); j++) 
				if(!columnsToRemove[j]) 
					newAligned += curAligned.substr(j, 1);
			
			Sequence newSeq(curName, newAligned);
			newDB.add(newSeq);
		}
		
		string newFileName = getRootName(filename) + "filter.fa";
		ofstream outfile;
		outfile.open(newFileName.c_str());
		newDB.print(outfile);
		outfile.close();
			
		globaldata->clear();
		//delete db;
		//delete newDB;
		
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