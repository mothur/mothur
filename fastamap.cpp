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
 FastaMap::FastaMap(ifstream& in) {
	//int numberOfSequences = 0;
	
	string name, sequence, line;
	sequence = "";
	
	getline(in, line);
	name = line.substr(1, line.length());  //rips off '>'
	
	//read through file
	while (getline(in, line)) {
		if (isalnum(line.at(0))){  //if it's a sequence line
			sequence += line;
		}
		else{
			//input sequence info into map
			it = data.find(sequence);
			if (it == data.end()) { 	//it's unique.
				data[sequence].groupname = name;  //group name will be the name of the first duplicate sequence found.
				data[sequence].groupnumber = 1;
				data[sequence].names = name;
			}else { // its a duplicate.
				data[sequence].names += "," + name;
				data[sequence].groupnumber++;	
			}
			name = (line.substr(1, (line.npos))); //The line you just read is a new name so rip off '>'
		}
	}
	
	//store last sequence and name info.
	it = data.find(sequence);
	if (it == data.end()) { 	//it's unique.
		data[sequence].groupname = name;  //group name will be the name of the first duplicate sequence found.
		data[sequence].groupnumber = 1;
		data[sequence].names = name;
	}else { // its a duplicate.
		data[sequence].names += "," + name;
		data[sequence].groupnumber++;	
	}	
}
/*******************************************************************************/
string FastaMap::getGroupName(string seq) {  //pass a sequence name get its group
	return data[seq].groupname;
}
/*******************************************************************************/
int FastaMap::getGroupNumber(string seq) {  //pass a sequence name get number of sequence in its group
	return data[seq].groupnumber;
}
/*******************************************************************************/
string FastaMap::getNames(string seq) {	//pass a sequence get the string of names in the group separated by ','s.
	return data[seq].names;
}
/*******************************************************************************/
void FastaMap::push_back(string seq, string Name) {//sequencename, name
	data[seq].groupname = Name;
	data[seq].groupnumber = 1;
	data[seq].names = Name;
}
/*******************************************************************************/
void FastaMap::clear() { //clears out data
	data.clear();
}
/*******************************************************************************/
int FastaMap::size(){ //returns datas size which is the number of unique sequences
	return data.size();
}
/*******************************************************************************/
void FastaMap::print(ostream&){ //prints data

}
/*******************************************************************************/
