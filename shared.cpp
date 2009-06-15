/*
 *  shared.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 12/5/08.
 *  Copyright 2008 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "shared.h"

/**************************************************************************************************/

Shared::Shared(){
   globaldata = GlobalData::getInstance();
}

/**************************************************************************************************/
void Shared::getSharedVectors(SharedListVector* list) {
		string label, group;
		int i,j;
		label = list->getLabel();
		
		for (it = sharedGroups.begin(); it != sharedGroups.end(); it++) {  delete it->second;  }
		sharedGroups.clear();  //removes old info.
		
		//initalize sharedGroups 
		for (j=0; j<globaldata->gGroupmap->getNumGroups(); j++) {//for each group
			group = globaldata->gGroupmap->namesOfGroups[j];
			sharedGroups[group] = new SharedRAbundVector();
			sharedGroups[group]->setLabel(label);
			for (i = 0; i<list->size(); i++) { //for each otu
				sharedGroups[group]->push_back(0, i, group); //initialize to 0.
			}
		}
			
		//fills sharedGroups
		for (i = 0; i<list->size(); i++) {
			parse(i, list);
		}
		
		//updates sharedVector
		sharedRAbund.push_back(sharedGroups);
}



/***********************************************************************/
void Shared::parse(int index, SharedListVector* list) {

		string prefix, suffix, groupsName;
		suffix = list->get(index);
	
		while (suffix.find_first_of(',') != -1) {//while you still have sequences
			prefix = suffix.substr(0,suffix.find_first_of(','));
			if ((suffix.find_first_of(',')+1) <= suffix.length()) {  //checks to make sure you don't have comma at end of string
				suffix = suffix.substr(suffix.find_first_of(',')+1, suffix.length());
			}
			groupsName = globaldata->gGroupmap->getGroup(prefix);
			if (groupsName != "not found") {
				sharedGroups[groupsName]->set(index, (sharedGroups[groupsName]->getAbundance(index) + 1), groupsName); //increment shared vector for that group
			}else {
				cerr << "Error: Sequence '" << prefix << "' was not found in the group file, please correct\n";
			}
		}
		
		//save last name after comma
		groupsName = globaldata->gGroupmap->getGroup(suffix);
		if (groupsName != "not found") {
			sharedGroups[groupsName]->set(index, ((sharedGroups[groupsName]->getAbundance(index)) + 1), groupsName); //increment shared vector for that group
		}else {
			cerr << "Error: Sequence '" << suffix << "' was not found in the group file, please correct\n";
		}
	}
