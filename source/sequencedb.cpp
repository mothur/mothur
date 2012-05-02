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

SequenceDB::SequenceDB() {  m = MothurOut::getInstance();  length = 0; samelength = true; }
/***********************************************************************/
//the clear function free's the memory
SequenceDB::~SequenceDB() { clear(); }

/***********************************************************************/

SequenceDB::SequenceDB(int newSize) {
	data.resize(newSize, Sequence());
	length = 0; samelength = true;
}

/***********************************************************************/

SequenceDB::SequenceDB(ifstream& filehandle) {
	try{
		length = 0; samelength = true;
				
		//read through file
		while (!filehandle.eof()) {
			//input sequence info into sequencedb
			Sequence newSequence(filehandle);
			
			if (newSequence.getName() != "") {   
				if (length == 0) { length = newSequence.getAligned().length(); }
				if (length != newSequence.getAligned().length()) { samelength = false; }
				data.push_back(newSequence);  
			}
			
			//takes care of white space
			m->gobble(filehandle);
		}

		filehandle.close();
		
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceDB", "SequenceDB");
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
		m->errorOut(e, "SequenceDB", "readName");
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
		
		if (length == 0) { length = sequence.length(); }
		if (length != sequence.length()) { samelength = false; }
		
		return sequence;
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceDB", "readSequence");
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
		if (length == 0) { length = newUnaligned.length(); }
		if (length != newUnaligned.length()) { samelength = false; }
		
		data[index] = Sequence(data[index].getName(), newUnaligned);
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceDB", "set");
		exit(1);
	}
}

/***********************************************************************/

void SequenceDB::set(int index, Sequence newSeq) {
	try {
		if (length == 0) { length = newSeq.getAligned().length(); }
		if (length != newSeq.getAligned().length()) { samelength = false; }

		data[index] = newSeq;
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceDB", "set");
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
		m->errorOut(e, "SequenceDB", "resize");
		exit(1);
	}
}

/***********************************************************************/

void SequenceDB::clear() {
	try {
		data.clear();
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceDB", "clear");
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
		m->errorOut(e, "SequenceDB", "print");
		exit(1);
	}
}
	
/***********************************************************************/

void SequenceDB::push_back(Sequence newSequence) {
	try {
		if (length == 0) { length = newSequence.getAligned().length(); }
		if (length != newSequence.getAligned().length()) { samelength = false; }

		data.push_back(newSequence);
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceDB", "push_back");
		exit(1);
	}
}

/***********************************************************************/

