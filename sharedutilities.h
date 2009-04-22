#ifndef SHAREDUTIL_H
#define SHAREDUTIL_H
/*
 *  sharedutilities.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "mothur.h"

class SharedRAbundVector;
class SharedOrderVector;

/**************************************************************************************************/

class SharedUtil {
	public:
		SharedUtil() {};
		~SharedUtil() {};
		
		void getSharedVectors(vector<string>, vector<SharedRAbundVector*>&, SharedOrderVector*);
		void getSharedVectorswithReplacement(vector<string>, vector<SharedRAbundVector*>&, SharedOrderVector*);
		void setGroups(vector<string>&, vector<string>&);  //globaldata->Groups, your tree or group map
		void setGroups(vector<string>&, vector<string>&, string);  //globaldata->Groups, your tree or group map, mode
		void setGroups(vector<string>&, vector<string>&, string&, int&, string);  //globaldata->Groups, your tree or group map, allgroups, numGroups, mode
		void getCombos(vector<string>&, vector<string>, int&); //groupcomb, globaldata->Groups, numcomb
		void updateGroupIndex(vector<string>&, map<string, int>&); //globaldata->Groups, groupmap->groupIndex
		
	private:
	
		bool isValidGroup(string, vector<string>);
};

/**************************************************************************************************/



#endif