/*
 *  getrepseqscommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/19/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "getrepseqscommand.h"

//**********************************************************************************************************************
GetRepSeqsCommand::GetRepSeqsCommand(){
	try {
		globaldata = GlobalData::getInstance();
		fastafile = globaldata->getFastaFile();
		namesfile = globaldata->getNameFile();
		openInputFile(fastafile, in);
		
		fasta = new FastaMap();
		
		//read in group map info.
		groupMap = new GroupMap(globaldata->getGroupFile());
		groupMap->readMap();
			
		//fill filehandles with neccessary ofstreams
		int i;
		ofstream* temp;
		//one for each group
		for (i=0; i<groupMap->getNumGroups(); i++) {
			temp = new ofstream;
			filehandles[groupMap->namesOfGroups[i]] = temp;
		}
		
		//one for shared
		temp = new ofstream;
		string s = "shared";
		filehandles[s] = temp;
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetRepSeqsCommand class Function GetRepSeqsCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetRepSeqsCommand class function GetRepSeqsCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************

GetRepSeqsCommand::~GetRepSeqsCommand(){
	delete input;
	delete read;
	delete fasta;
	delete list;
}

//**********************************************************************************************************************

int GetRepSeqsCommand::execute(){
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
				
				cout << list->getLabel() << '\t' << count << endl;
				
				//open output list files
				for (int i=0; i<groupMap->getNumGroups(); i++) {//opens an output file for each group
					openOutputFile(fastafile + groupMap->namesOfGroups[i] + list->getLabel() + ".fasta", *(filehandles[groupMap->namesOfGroups[i]]));
					used[groupMap->namesOfGroups[i]] = false;
				}
				string s = "shared";
				openOutputFile(fastafile + s + list->getLabel() + ".fasta", *(filehandles[s]));
				used[s] = false;
				
				
				//for each bin in the list vector
				for (int i = 0; i < list->size(); i++) {
					seq.clear();
					//uses this to determine if the bin is unique to one group or if it is shared
					map<string, string> groups;

					//determine if this otu is unique to one group or not
					binnames = list->get(i);
					while (binnames.find_first_of(',') != -1) { 
						//parse out each name in bin
						name = binnames.substr(0,binnames.find_first_of(','));
						binnames = binnames.substr(binnames.find_first_of(',')+1, binnames.length());
						
						//do work for that name
						sequence = fasta->getSequence(name);
						if (sequence != "not found") {
							string group = groupMap->getGroup(name);
							if (group != "not found") {  groups[group] = group;	}  //add group to list of groups in this bin
							else {	
								cout << "error sequence " << name << " is not assigned a group in your groupfile. Please correct." << endl;
								removeFiles(list->getLabel());
								return 0;
							}
							name = ">" + name + "|" + toString(i+1);
							seq[name] = sequence;
						}else { 
							cout << name << " is missing from your fasta or name file. Please correct. " << endl; 
							removeFiles(list->getLabel());
							return 0;
						}
						
					}
					
					//get last name
					sequence = fasta->getSequence(binnames);
					if (sequence != "not found") {
						string group = groupMap->getGroup(binnames);
						if (group != "not found") {  groups[group] = group;	}  //add group to list of groups in this bin
						else {	
							cout << "error sequence " << binnames << " is not assigned a group in your groupfile. Please correct." << endl;
							removeFiles(list->getLabel());
							return 0;
						}
						binnames = ">" + binnames + "|" + toString(i+1);  //attach bin number to name
						seq[binnames] = sequence;
					}else { 
						cout << binnames << " is missing from your fasta or name file. Please correct. " << endl; 
						removeFiles(list->getLabel());
						return 0;
					}
					
					//output each bin to files
					//what file does this bin need to be outputted to 
					if (groups.size() == 1) { //this bin is unique to one group
						it3 = groups.begin();
						string uniqueGroup = it3->first;
						used[uniqueGroup] = true;
						//print out sequences from that bin to shared file
						for (it3 = seq.begin(); it3 != seq.end(); it3++){
							*(filehandles[uniqueGroup]) << it3->first << endl;
							*(filehandles[uniqueGroup]) << it3->second << endl;
						}
					}else {//this bin has sequences from multiple groups in it
						used[s] = true;
						//print out sequences from that bin to shared file
						for (it3 = seq.begin(); it3 != seq.end(); it3++){
							*(filehandles[s]) << it3->first << endl;
							*(filehandles[s]) << it3->second << endl;
						}
					}
				}
				
				//close ostreams and remove unused files
				for (it = filehandles.begin(); it != filehandles.end(); it++) {
					it->second->close();
					if (used[it->first] == false) { string filename = fastafile + it->first + list->getLabel() + ".fasta";  remove(filename.c_str());  }
				}

			}
			
			delete list;
			list = input->getListVector();
			count++;
		}
		
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetRepSeqsCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetRepSeqsCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************
void GetRepSeqsCommand::readNamesFile() {
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
		cout << "Standard Error: " << e.what() << " has occurred in the GetRepSeqsCommand class Function readNamesFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetRepSeqsCommand class function readNamesFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************
void GetRepSeqsCommand::removeFiles(string label) {
	try {
			//close ostreams
			for (it = filehandles.begin(); it != filehandles.end(); it++) {
				it->second->close();
			}

			//remove output files because there was an error
			for (int i=0; i<groupMap->getNumGroups(); i++) {
				string outputFileName = fastafile + groupMap->namesOfGroups[i] + label + ".fasta";
				remove(outputFileName.c_str());
			}
			string outputFileName = fastafile + "shared"+ label + ".fasta";
			remove(outputFileName.c_str());

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetRepSeqsCommand class Function removeFiles. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetRepSeqsCommand class function removeFiles. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************

