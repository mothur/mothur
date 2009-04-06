/*
 *  getoturepcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/6/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "getoturepcommand.h"

//**********************************************************************************************************************
GetOTURepCommand::GetOTURepCommand(){
	try{
		globaldata = GlobalData::getInstance();
	
		if(globaldata->gSparseMatrix != NULL)	{	matrix = new SparseMatrix(*globaldata->gSparseMatrix);		}
		
		//listOfNames bin 0 = first name read in distance matrix, listOfNames bin 1 = second name read in distance matrix
		if(globaldata->gListVector != NULL)		{	
			listOfNames = new ListVector(*globaldata->gListVector);	
			
			//map names to rows in sparsematrix
			for (int i = 0; i < listOfNames->size(); i++) {
				nameToIndex[listOfNames->get(i)] = i;
			}
		}else { cout << "error" << endl; }

		
		fastafile = globaldata->getFastaFile();
		namesfile = globaldata->getNameFile();
		openInputFile(fastafile, in);
		
		fasta = new FastaMap();

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetOTURepCommand class Function GetOTURepCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetOTURepCommand class function GetOTURepCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
//**********************************************************************************************************************

GetOTURepCommand::~GetOTURepCommand(){
	delete matrix;
	delete list;
	delete input;
	delete read;
	delete fasta;
}

//**********************************************************************************************************************

int GetOTURepCommand::execute(){
	try {
		int count = 1;
		string nameRep, name, sequence;
		
		//read fastafile
		fasta->readFastaFile(in);
		
		//set format to list so input can get listvector
		globaldata->setFormat("list");
		
		//if user gave a namesfile then use it
		if (namesfile != "") {
			readNamesFile();
		}
		
		//read list file
		read = new ReadPhilFile(globaldata->getListFile());	
		read->read(&*globaldata); 
		
		input = globaldata->ginput;
		list = globaldata->gListVector;
		
		while(list != NULL){
			
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(list->getLabel()) == 1){
				
				//create output file
				string outputFileName = getRootName(globaldata->getListFile()) + list->getLabel() + ".fastarep";
				openOutputFile(outputFileName, out);

				cout << list->getLabel() << '\t' << count << endl;
				
				//for each bin in the list vector
				for (int i = 0; i < list->size(); i++) {
					nameRep = FindRep(i);
					
					//print out name and sequence for that bin
					sequence = fasta->getSequence(nameRep);

					if (sequence != "not found") {
						nameRep = nameRep + "bin" + toString(i+1);
						out << ">" << nameRep << endl;
						out << sequence << endl;
					}else { 
						cout << nameRep << " is missing from your fasta or name file. Please correct. " << endl; 
						remove(outputFileName.c_str());
						return 0;
					}
				}
			}
			
			list = input->getListVector();
			count++;
		}

		
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetOTURepCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetOTURepCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

//**********************************************************************************************************************
void GetOTURepCommand::readNamesFile() {
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
		cout << "Standard Error: " << e.what() << " has occurred in the GetOTURepCommand class Function readNamesFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetOTURepCommand class function readNamesFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************
string GetOTURepCommand::FindRep(int bin) {
	try{
		vector<string> names;
		map<string, float> sums;
		map<string, float>::iterator it4;
		map<int, string> binMap; //subset of namesToIndex - just member of this bin
		string binnames;
		float min = 10000;
		string minName;
		
		binnames = list->get(bin);
		
		//parse names into vector
		splitAtComma(binnames, names);
		
		//if only 1 sequence in bin then that's the rep
		if (names.size() == 1) { return names[0]; }
		else {
			//fill binMap
			for (int i = 0; i < names.size(); i++) {
				for (it3 = nameToIndex.begin(); it3 != nameToIndex.end(); it3++) {
					if (it3->first == names[i]) {  
						binMap[it3->second] = it3->first;

						//initialize sums map
						sums[it3->first] = 0.0;
						break;
					}
				}
			}
			
			//go through each cell in the sparsematrix
			for(MatData currentCell = matrix->begin(); currentCell != matrix->end(); currentCell++){
				//is this a distance between 2 members of this bin
				it = binMap.find(currentCell->row);
				it2 = binMap.find(currentCell->column);
				
				//sum the distance of the sequences in the bin to eachother
				if ((it != binMap.end()) && (it2 != binMap.end())) {
					//this is a cell that repesents the distance between to of this bins members
					sums[it->second] += currentCell->dist;
					sums[it2->second] += currentCell->dist;
				}
			}
			
			//smallest sum is the representative
			for (it4 = sums.begin(); it4 != sums.end(); it4++) {
				if (it4->second < min) {
					min = it4->second;
					minName = it4->first;
				}

			}
			
			return minName;
		}
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetOTURepCommand class Function FindRep. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetOTURepCommand class function FindRep. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}





