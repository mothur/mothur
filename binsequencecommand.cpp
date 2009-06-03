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
		groupfile = globaldata->getGroupFile();
		openInputFile(fastafile, in);
		
		if (groupfile != "") {
			//read in group map info.
			groupMap = new GroupMap(groupfile);
			groupMap->readMap();
		}
		
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
	if (groupfile != "") {
		delete groupMap;
	}
}

//**********************************************************************************************************************

int BinSeqCommand::execute(){
	try {
		int count = 1;
		int error = 0;
		
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
		ListVector* lastList = list;
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = globaldata->labels;
		set<int> userLines = globaldata->lines;

				
		while((list != NULL) && ((globaldata->allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
			
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(list->getLabel()) == 1){
				
				error = process(list, count);	
				if (error == 1) { return 0; }	
							
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				userLines.erase(count);
			}
			
			if ((anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastList->getLabel()) != 1)) {
				
				error = process(lastList, count);	
				if (error == 1) { return 0; }
													
				processedLabels.insert(lastList->getLabel());
				userLabels.erase(lastList->getLabel());
				
			}
			
			if (count != 1) { delete lastList; }
			lastList = list;			

			list = input->getListVector();
			count++;
		}
		
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			cout << "Your file does not include the label "<< *it; 
			if (processedLabels.count(lastList->getLabel()) != 1) {
				cout << ". I will use " << lastList->getLabel() << "." << endl;
				needToRun = true;
			}else {
				cout << ". Please refer to " << lastList->getLabel() << "." << endl;
			}
		}
		
		//run last line if you need to
		if (needToRun == true)  {
			error = process(lastList, count);	
			if (error == 1) { return 0; }			
		}
		
		delete lastList;
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
//return 1 if error, 0 otherwise
int BinSeqCommand::process(ListVector* list, int count) {
	try {
				string binnames, name, sequence;
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
							//if you don't have groups
							if (groupfile == "") {
								name = name + "|" + toString(i+1);
								out << ">" << name << endl;
								out << sequence << endl;
							}else {//if you do have groups
								string group = groupMap->getGroup(name);
								if (group == "not found") {  
									cout << name << " is missing from your group file. Please correct. " << endl;
									remove(outputFileName.c_str());
									return 1;
								}else{
									name = name + "|" + group + "|" + toString(i+1);
									out << ">" << name << endl;
									out << sequence << endl;
								}
							}
						}else { 
							cout << name << " is missing from your fasta or name file. Please correct. " << endl; 
							remove(outputFileName.c_str());
							return 1;
						}
						
					}
					
					//get last name
					sequence = fasta->getSequence(binnames);
					if (sequence != "not found") {
						//if you don't have groups
						if (groupfile == "") {
							binnames = binnames + "|" + toString(i+1);
							out << ">" << binnames << endl;
							out << sequence << endl;
						}else {//if you do have groups
							string group = groupMap->getGroup(binnames);
							if (group == "not found") {  
								cout << binnames << " is missing from your group file. Please correct. " << endl;
								remove(outputFileName.c_str());
								return 1;
							}else{
								binnames = binnames + "|" + group + "|" + toString(i+1);
								out << ">" << binnames << endl;
								out << sequence << endl;
							}
						}
					}else { 
						cout << binnames << " is missing from your fasta or name file. Please correct. " << endl; 
						remove(outputFileName.c_str());
						return 1;
					}
				}
					
				out.close();
				return 0;

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


