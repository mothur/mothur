/*
 *  fastamap.cpp
 *  mothur
 *
 *  Created by Sarah Westcott on 1/16/09.
 *  Copyright 2009 Schloss Lab UMASS AMherst. All rights reserved.
 *
 */

#include "fastamap.h"

/*******************************************************************************/
void FastaMap::readFastaFile(ifstream& in) {
	try {
		string name, sequence, line;
		sequence = "";
		int c;
		string temp;
		
		
		//read through file
		while ((c = in.get()) != EOF) {
			name = ""; sequence = ""; 
			//is this a name
			if (c == '>') { 
				name = readName(in); 
				sequence = readSequence(in); 
			}else {  cout << "Error fasta in your file. Please correct." << endl; }

			//store info in map
			//input sequence info into map
			seqmap[name] = sequence;  
			it = data.find(sequence);
			if (it == data.end()) { 	//it's unique.
				data[sequence].groupname = name;  //group name will be the name of the first duplicate sequence found.
				data[sequence].groupnumber = 1;
				data[sequence].names = name;
			}else { // its a duplicate.
				data[sequence].names += "," + name;
				data[sequence].groupnumber++;
			}	
			
			gobble(in);
		}
					
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the FastaMap class Function readFastaFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the FastaMap class function readFastaFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/*******************************************************************************/
string FastaMap::readName(ifstream& in) {
	try{
		string name = "";
		int c;
		string temp;
		
		while ((c = in.get()) != EOF) {
			//if c is not a line return
			if (c != 10) {
				name += c;
			}else { break;  }
		}
			
		return name;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the FastaMap class Function readName. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the FastaMap class function readName. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/*******************************************************************************/
string FastaMap::readSequence(ifstream& in) {
	try{
		string sequence = "";
		string line;
		int pos, c;
		
		while (!in.eof()) {
			//save position in file in case next line is a new name.
			pos = in.tellg();
			line = "";
			in >> line;			
			//if you are at a new name
			if (line[0] == '>') {
				//put file pointer back since you are now at a new name
				in.seekg(pos, ios::beg);
				c = in.get();  //because you put it back to a newline char
				break;
			}else {  sequence += line;	}
		}
			
		return sequence;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the FastaMap class Function readSequence. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the FastaMap class function readSequence. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/*******************************************************************************/
string FastaMap::getGroupName(string seq) {  //pass a sequence name get its group
	return data[seq].groupname;
}
/*******************************************************************************/
string FastaMap::getNames(string seq) {	//pass a sequence get the string of names in the group separated by ','s.
	return data[seq].names;
}
/*******************************************************************************/
int FastaMap::getGroupNumber(string seq) {	//pass a sequence get the number of identical sequences.
	return data[seq].groupnumber;
}
/*******************************************************************************/
string FastaMap::getSequence(string name) {
	it2 = seqmap.find(name);
	if (it2 == seqmap.end()) { 	//it's not found
		return "not found";
	}else { // found it
		return it2->second;
	}
}	
/*******************************************************************************/
void FastaMap::push_back(string name, string seq) {
	it = data.find(seq);
	if (it == data.end()) { 	//it's unique.
		data[seq].groupname = name;  //group name will be the name of the first duplicate sequence found.
		data[seq].groupnumber = 1;
		data[seq].names = name;
	}else { // its a duplicate.
		data[seq].names += "," + name;
		data[seq].groupnumber++;
	}
	
	seqmap[name] = seq;
}
/*******************************************************************************/
int FastaMap::sizeUnique(){ //returns datas size which is the number of unique sequences
	return data.size();
}
/*******************************************************************************/
void FastaMap::printNamesFile(ostream& out){ //prints data
	try {
		// two column file created with groupname and them list of identical sequence names
		for (it = data.begin(); it != data.end(); it++) {
			out << it->second.groupname << '\t' << it->second.names << endl;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the FastaMap class Function print. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the FastaMap class function print. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/*******************************************************************************/
void FastaMap::printCondensedFasta(ostream& out){ //prints data
	try {
		//creates a fasta file
		for (it = data.begin(); it != data.end(); it++) {
			out << ">" << it->second.groupname << endl;
			out << it->first << endl;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the FastaMap class Function print. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the FastaMap class function print. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/*******************************************************************************/

