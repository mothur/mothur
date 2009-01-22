#ifndef SHARED_H
#define SHARED_H

/*
 *  shared.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 12/5/08.
 *  Copyright 2008 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
 
using namespace std;

#include <Carbon/Carbon.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include "sharedrabundvector.h"
#include "listvector.hpp"
#include "globaldata.hpp"

class Shared {
	public:
		Shared();
		~Shared();
		void getSharedVectors(int, ListVector*);
		map<string, SharedRAbundVector*> sharedGroups; //string is groupname, SharedVector* is out info for that group
		
	private:
		GlobalData* globaldata;
		void parse(int, ListVector*);
		vector< map<string, SharedRAbundVector*> > sharedRAbund;  //contains all the info needed to create the .shared file not sure if we will need 
};

#endif