#ifndef GROUPMAP_H
#define GROUPMAP_H
/*
 *  groupmap.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 12/1/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <Carbon/Carbon.h>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include "utilities.hpp"


class GroupMap {
public:
	GroupMap(string);
	~GroupMap();
	void readMap();
	int getNumGroups();
	string getGroup(string);
	vector<string> namesOfGroups;
	
	
private:
	ifstream fileHandle;
	string groupFileName;
	int numGroups;
	map<string, string>::iterator it;
	void setNamesOfGroups(string); 
	map<string, string> groupmap; //sequence name and groupname
};

#endif