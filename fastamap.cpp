/*
 *  fastamap.cpp
 *  mothur
 *
 *  Created by Sarah Westcott on 1/16/09.
 *  Copyright 2009 Schloss Lab UMASS AMherst. All rights reserved.
 *
 */

#include "fastamap.h"
#include "sequence.hpp"

/*******************************************************************************/

void FastaMap::readFastaFile(string inFileName) {
	try {
		ifstream in;
		openInputFile(inFileName, in);
		string name, sequence, line;
		sequence = "";
		string temp;

		while(!in.eof()){
			Sequence currSeq(in);
			name = currSeq.getName();
			
			if (name != "") {
				if(currSeq.getIsAligned())	{	sequence = currSeq.getAligned();	}
				else						{	sequence = currSeq.getUnaligned();	}
				
				seqmap[name] = sequence;  
				map<string,group>::iterator it = data.find(sequence);
				if (it == data.end()) { 	//it's unique.
					data[sequence].groupname = name;  //group name will be the name of the first duplicate sequence found.
					//				data[sequence].groupnumber = 1;
					data[sequence].names = name;
				}else { // its a duplicate.
					data[sequence].names += "," + name;
					//				data[sequence].groupnumber++;
				}	
			}
			gobble(in);
		}
		in.close();		
	}
	catch(exception& e) {
		m->errorOut(e, "FastaMap", "readFastaFile");
		exit(1);
	}
}

/*******************************************************************************/

void FastaMap::readFastaFile(string inFastaFile, string oldNameFileName){ //prints data
	
	ifstream oldNameFile;
	openInputFile(oldNameFileName, oldNameFile);
	
	map<string,string> oldNameMap;
	string name, list;
	while(!oldNameFile.eof()){
		oldNameFile >> name >> list;
		oldNameMap[name] = list;
		gobble(oldNameFile);
	}
	oldNameFile.close();
	
	ifstream inFASTA;
	openInputFile(inFastaFile, inFASTA);
	string sequence;
	while(!inFASTA.eof()){
		Sequence currSeq(inFASTA);
		name = currSeq.getName();
		
		if (name != "") {
			if(currSeq.getIsAligned())	{	sequence = currSeq.getAligned();	}
			else						{	sequence = currSeq.getUnaligned();	}
			
			seqmap[name] = sequence;  
			map<string,group>::iterator it = data.find(sequence);
			if (it == data.end()) { 	//it's unique.
				data[sequence].groupname = name;  //group name will be the name of the first duplicate sequence found.
				//			data[sequence].groupnumber = 1;
				data[sequence].names = oldNameMap[name];
			}else { // its a duplicate.
				data[sequence].names += "," + oldNameMap[name];
				//			data[sequence].groupnumber++;
			}	
		}
		gobble(inFASTA);
	}
	
	
	inFASTA.close();
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

string FastaMap::getSequence(string name) {
	
	map<string,string>::iterator it = seqmap.find(name);
	if (it == seqmap.end()) { 	return "not found";		}
	else					{	return it->second;		}
	
}	

/*******************************************************************************/

void FastaMap::push_back(string name, string seq) {
	
	map<string,group>::iterator it = data.find(seq);
	if (it == data.end()) { 	//it's unique.
		data[seq].groupname = name;  //group name will be the name of the first duplicate sequence found.
		data[seq].names = name;
	}else { // its a duplicate.
		data[seq].names += "," + name;
	}
	seqmap[name] = seq;
}

/*******************************************************************************/

int FastaMap::sizeUnique(){ //returns datas size which is the number of unique sequences
	return data.size();
}

/*******************************************************************************/

void FastaMap::printNamesFile(string outFileName){ //prints data
	try {
		ofstream outFile;
		openOutputFile(outFileName, outFile);
		
		// two column file created with groupname and them list of identical sequence names
		for (map<string,group>::iterator it = data.begin(); it != data.end(); it++) {
			outFile << it->second.groupname << '\t' << it->second.names << endl;
		}
		outFile.close();
	}
	catch(exception& e) {
		m->errorOut(e, "FastaMap", "printNamesFile");
		exit(1);
	}
}

/*******************************************************************************/

void FastaMap::printCondensedFasta(string outFileName){ //prints data
	try {
		ofstream out;
		openOutputFile(outFileName, out);
		//creates a fasta file
		for (map<string,group>::iterator it = data.begin(); it != data.end(); it++) {
			out << ">" << it->second.groupname << endl;
			out << it->first << endl;
		}
		out.close();
	}
	catch(exception& e) {
		m->errorOut(e, "FastaMap", "printCondensedFasta");
		exit(1);
	}
}

/*******************************************************************************/

