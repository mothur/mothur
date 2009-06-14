/*
 *  sequencedb.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/13/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sequencedb.h"
#include "sequence.hpp"
#include "mothur.h"
#include "calculator.h"


/***********************************************************************/

SequenceDB::SequenceDB() {}
/***********************************************************************/
//the clear function free's the memory
SequenceDB::~SequenceDB() { clear(); }

/***********************************************************************/

SequenceDB::SequenceDB(int newSize) {
	data.resize(newSize, Sequence());
}

/***********************************************************************/

SequenceDB::SequenceDB(ifstream& filehandle) {
	try{
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
			data.push_back(newSequence);
			
			//takes care of white space
			gobble(filehandle);
		}

		filehandle.close();
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SequenceDB class Function SequenceDB. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SequenceDB class function SequenceDB. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/*******************************************************************************/
string SequenceDB::readName(ifstream& in) {
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
		cout << "Standard Error: " << e.what() << " has occurred in the SequenceDB class Function readName. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SequenceDB class function readName. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/*******************************************************************************/
string SequenceDB::readSequence(ifstream& in) {
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
		cout << "Standard Error: " << e.what() << " has occurred in the SequenceDB class Function readSequence. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SequenceDB class function readSequence. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
	
/***********************************************************************/

int SequenceDB::getNumSeqs() {
	return data.size();
}

/***********************************************************************/

void SequenceDB::set(int index, string newUnaligned) {
	try {
		data[index] = Sequence(data[index].getName(), newUnaligned);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SequenceDB class Function set. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SequenceDB class function set. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

void SequenceDB::set(int index, Sequence newSeq) {
	try {
		data[index] = newSeq;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SequenceDB class Function set. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SequenceDB class function set. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

Sequence SequenceDB::get(int index) {
	return data[index];
}

/***********************************************************************/

void SequenceDB::resize(int newSize) {
	try {
		data.resize(newSize);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SequenceDB class Function resize. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SequenceDB class function resize. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

void SequenceDB::clear() {
	try {
		data.clear();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SequenceDB class Function clear. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SequenceDB class function clear. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

int SequenceDB::size() {
	return data.size();
}

/***********************************************************************/

void SequenceDB::print(ostream& out) {
	try {
		for(int i = 0; i < data.size(); i++) {
			data[i].printSequence(out);
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SequenceDB class Function print. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SequenceDB class function print. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
	
/***********************************************************************/

void SequenceDB::push_back(Sequence newSequence) {
	try {
		data.push_back(newSequence);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SequenceDB class Function add. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SequenceDB class function add. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

