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

#include "mothur.h"
#include "sharedrabundvector.h"
#include "sharedlistvector.h"
#include "globaldata.hpp"

class Shared {
	public:
		Shared();
		~Shared() {};
		void getSharedVectors(SharedListVector*);
		map<string, SharedRAbundVector*> sharedGroups; //string is groupname, SharedVector* is out info for that group
		
	private:
		GlobalData* globaldata;
		map<string, SharedRAbundVector*>::iterator it;
		void parse(int, SharedListVector*);
		//vector< map<string, SharedRAbundVector*> > sharedRAbund;  //contains all the info needed to create the .shared file not sure if we will need 
};

#endif
