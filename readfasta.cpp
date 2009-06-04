/*
 *  readfasta.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/21/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "readfasta.h"
#include <iostream>
#include <fstream>

/*******************************************************************************/
ReadFasta::ReadFasta(string file) : ReadSeqs(file) {}
/*******************************************************************************/
ReadFasta::~ReadFasta(){
	//for(int i = 0; i < sequencedb.getNumSeqs(); i++)
		//delete sequencedb.get(i);
}
/*******************************************************************************/
void ReadFasta::read() {
	try {
	/*string name = "";
	string sequence = "";
	string temp;
	int count = 0;
	
	while(!filehandle.eof()){
		if(count == 0)
			filehandle >> temp;
		if(temp.substr(0,1).compare(">") == 0) {
			if(count != 0) {
				Sequence newSequence(name, sequence);
				sequencedb.add(newSequence);
				sequence = "";
			}
			else
				count++;
			name = temp.substr(1,temp.length()-1);
		}
		else {
			sequence += temp;
		}
		
		filehandle >> temp;
		gobble(filehandle);
		
		if(filehandle.eof())
			sequence += temp;
			
	}
	Sequence newSequence(name, sequence);
	sequencedb.add(newSequence); */

		string name, sequence, line;
		sequence = "";
		int c;
		string temp;
		
		
		//read through file
		while ((c = filehandle.get()) != EOF) {
			name = ""; sequence = ""; 
			//is this a name
			if (c == '>') { 
				name = readName(filehandle); 
				sequence = readSequence(filehandle); 
			}else {  cout << "Error fasta in your file. Please correct." << endl; }

			//input sequence info into sequencedb
			Sequence newSequence(name, sequence);
			sequencedb.add(newSequence);
			
			//takes care of white space
			gobble(filehandle);
		}

		filehandle.close();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadFasta class Function read. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadFasta class function read. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/*********************************************************************************/
SequenceDB* ReadFasta::getDB() {
	return &sequencedb;
}
/*******************************************************************************/
string ReadFasta::readName(ifstream& in) {
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
		cout << "Standard Error: " << e.what() << " has occurred in the ReadFasta class Function readName. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadFasta class function readName. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/*******************************************************************************/
string ReadFasta::readSequence(ifstream& in) {
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
		cout << "Standard Error: " << e.what() << " has occurred in the ReadFasta class Function readSequence. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadFasta class function readSequence. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/*******************************************************************************/

