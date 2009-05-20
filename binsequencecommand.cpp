/*
 *  binsequencecommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/3/09.
 *  Copyright 2009 Schloss Lab UMASS Amhers. All rights reserved.
 *
 */

#include "binsequencecommand.h"

//**********************************************************************************************************************
BinSeqCommand::BinSeqCommand(){
	try {
		globaldata = GlobalData::getInstance();
		fastafile = globaldata->getFastaFile();
		namesfile = globaldata->getNameFile();
		openInputFile(fastafile, in);
		
		fasta = new FastaMap();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the BinSeqCommand class Function BinSeqCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the BinSeqCommand class function BinSeqCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************

BinSeqCommand::~BinSeqCommand(){
	delete input;
	delete read;
	delete fasta;
	delete list;
}

//**********************************************************************************************************************

int BinSeqCommand::execute(){
	try {
		int count = 1;
		string binnames, name, sequence;
		
		//read fastafile
		fasta->readFastaFile(in);
		
		//set format to list so input can get listvector
		globaldata->setFormat("list");
		
		//if user gave a namesfile then use it
		if (namesfile != "") {
			readNamesFile();
		}
		
		//read list file
		read = new ReadOTUFile(globaldata->getListFile());	
		read->read(&*globaldata); 
		
		input = globaldata->ginput;
		list = globaldata->gListVector;
				
		while(list != NULL){
			
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(list->getLabel()) == 1){
				
				string outputFileName = getRootName(globaldata->getListFile()) + list->getLabel() + ".fasta";
				openOutputFile(outputFileName, out);

				cout << list->getLabel() << '\t' << count << endl;
				
				//for each bin in the list vector
				for (int i = 0; i < list->size(); i++) {

					binnames = list->get(i);
					while (binnames.find_first_of(',') != -1) { 
						name = binnames.substr(0,binnames.find_first_of(','));
						binnames = binnames.substr(binnames.find_first_of(',')+1, binnames.length());
						
						//do work for that name
						sequence = fasta->getSequence(name);
						if (sequence != "not found") {
							name = name + "|" + toString(i+1);
							out << ">" << name << endl;
							out << sequence << endl;
						}else { 
							cout << name << " is missing from your fasta or name file. Please correct. " << endl; 
							remove(outputFileName.c_str());
							return 0;
						}
						
					}
					
					//get last name
					sequence = fasta->getSequence(binnames);
					if (sequence != "not found") {
						name = binnames + '|' + toString(i+1);
						out << ">" << name << endl;
						out << sequence << endl;
					}else { 
						cout << binnames << " is missing from your fasta or name file. Please correct. " << endl; 
						remove(outputFileName.c_str());
						return 0;
					}
					
				}
				out.close();
			}
			
			delete list;
			list = input->getListVector();
			count++;
		}
		
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the BinSeqCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the BinSeqCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************
void BinSeqCommand::readNamesFile() {
	try {
		vector<string> dupNames;
		openInputFile(namesfile, inNames);
		
		string name, names, sequence;
	
		while(inNames){
			inNames >> name;			//read from first column  A
			inNames >> names;		//read from second column  A,B,C,D
			
			dupNames.clear();
			
			//parse names into vector
			splitAtComma(names, dupNames);
			
			//store names in fasta map
			sequence = fasta->getSequence(name);
			for (int i = 0; i < dupNames.size(); i++) {
				fasta->push_back(dupNames[i], sequence);
			}
		
			gobble(inNames);
		}
		inNames.close();

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the BinSeqCommand class Function readNamesFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the BinSeqCommand class function readNamesFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************



