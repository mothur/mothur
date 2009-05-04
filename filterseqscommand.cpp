/*
 *  filterseqscommand.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/4/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "filterseqscommand.h"
#include <iostream>
#include <fstream>

/**************************************************************************************/
void FilterSeqsCommand::doTrump() {
	//trump = globaldata->getTrump();
//	
//	for(int i = 0; i < db->size(); i++) {
//		Sequence cur = db->get(i);
//		string curAligned = cur.getAligned();
//		
//		for(int j = 0; j < curAligned.length-1; j++) {
//			string curChar = curAligned.substr(j, j+1);
//			
//			if(curChar.compare(trump) == 0) 
//				columnsToRemove[j] = true;
//		}
//	}
}

/**************************************************************************************/
void FilterSeqsCommand::doSoft() {
	//soft = atoi(globaldata->getSoft().c_str());
//	vector<vector<int> > columnSymbolSums;
//	vector<vector<string> > columnSymbols;
//	for(int i = 0; i < db->get(0).getLength(); i++) {
//		vector<string> symbols;
//		vector<int> sums;
//		columnSymbols[i] = symbols;
//		columnSymbolSums[i] = sums;
//	}
//	
//	for(int i = 0; i < db->size(); i++) {
//		Sequence cur = db->get(i);
//		string curAligned = cur.getAligned();
//		
//		for(int j = 0; j < curAligned.length-1; j++) {
//			string curChar = curAligned.substr(j, j+1);
//			vector<string> curColumnSymbols = columnSymbols[j];
//			
//			bool newSymbol = true;
//			
//			for(int k = 0; j < curColumnSymbols.size(); j++) 
//				if(curChar.compare(curColumnSymbols[k]) == 0) {
//					newSymbol = false;
//					columnSymbolSums[j][k]++;
//				}
//			
//			if(newSymbol) {
//				columnSymbols.push_back(curChar);
//				columnSymbolSums[j].push_back(1);
//			}
//		}
//	}
//	
//	for(int i = 0; i < columnSymbolSums.size(); i++) {
//		int totalSum = 0;
//		int max = 0;
//		vector<int> curColumn = columnSymbolSums[i];
//		
//		for(int j = 0; j < curColumn.size(); j++) {
//			int curSum = curColumn[j];
//			if(curSum > max)
//				max = curSum;
//			totalSum += curSum;
//		}
//		
//		if((double)max/(double)totalSum * 100 < soft)
//			columnsToRemove[i] = true;
//	}
}
void FilterSeqsCommand::doFilter() {}
/**************************************************************************************/
int FilterSeqsCommand::execute() {	
	try {
		globaldata = GlobalData::getInstance();
		filename = globaldata->inputFileName;
		
		if(globaldata->getFastaFile().compare("") != 0) {
			readFasta = new ReadFasta(filename);
			readFasta->read();
			db = readFasta->getDB();
		}
		
		else if(globaldata->getNexusFile().compare("") != 0) {
			readNexus = new ReadNexus(filename);
			readNexus->read();
			db = readNexus->getDB();
		}
		
		else if(globaldata->getClustalFile().compare("") != 0) {
			readClustal = new ReadClustal(filename);
			readClustal->read();
			db = readClustal->getDB();
		}

		else if(globaldata->getPhylipFile().compare("") != 0) {
			readPhylip = new ReadPhylip(filename);
			readPhylip->read();
			db = readPhylip->getDB();
		}
	
		for(int i = 0; i < db->get(0).getLength(); i++) 
			columnsToRemove[i] = false;
			
		// Trump
		if(globaldata->getTrump().compare("") != 0) {
		
			
		}
		
		// Soft
		if(globaldata->getSoft().compare("") != 0) {}

		
			
		
		// Filter
		//if(globaldata->getFilter().compare("") != 0) {
//
//			filter = globaldata->getFilter();
//			ifstream filehandle;
//			openInputFile(filter, filehandle);
//			
//			char c;
//			int count = 0;
//			while(!filehandle.eof()) {
//				c = filehandle.get();
//				if(c == '0') 
//					columnsToRemove[count] = true;
//				count++;
//			}
//		}
		
		
		
			
			
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
