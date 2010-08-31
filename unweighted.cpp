/*
 *  unweighted.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "unweighted.h"

/**************************************************************************************************/

EstOutput Unweighted::getValues(Tree* t, int p, string o) {
	try {
		globaldata = GlobalData::getInstance();
		processors = p;
		outputDir = o;
			
		//if the users enters no groups then give them the score of all groups
		int numGroups = globaldata->Groups.size();
		
		//calculate number of comparsions
		int numComp = 0;
		vector< vector<string> > namesOfGroupCombos;
		for (int r=0; r<numGroups; r++) { 
			for (int l = r+1; l < numGroups; l++) {
				numComp++;
				vector<string> groups; groups.push_back(globaldata->Groups[r]); groups.push_back(globaldata->Groups[l]);
				namesOfGroupCombos.push_back(groups);
			}
		}
		
		if (numComp != 1) {
			vector<string> groups;
			if (numGroups == 0) {
				//get score for all users groups
				for (int i = 0; i < tmap->namesOfGroups.size(); i++) {
					if (tmap->namesOfGroups[i] != "xxx") {
						groups.push_back(tmap->namesOfGroups[i]);
					}
				}
				namesOfGroupCombos.push_back(groups);
			}else {
				for (int i = 0; i < globaldata->Groups.size(); i++) {
					groups.push_back(globaldata->Groups[i]);
				}
				namesOfGroupCombos.push_back(groups);
			}
		}

		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			if(processors == 1){
				data = driver(t, namesOfGroupCombos, 0, namesOfGroupCombos.size());
			}else{
				int numPairs = namesOfGroupCombos.size();
				
				int numPairsPerProcessor = numPairs / processors;
				
				for (int i = 0; i < processors; i++) {
					int startPos = i * numPairsPerProcessor;
					if(i == processors - 1){
						numPairsPerProcessor = numPairs - i * numPairsPerProcessor;
					}
					lines.push_back(new linePair(startPos, numPairsPerProcessor));
				}

				data = createProcesses(t, namesOfGroupCombos);
				
				for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
			}
		#else
			data = driver(t, namesOfGroupCombos, 0, namesOfGroupCombos.size());
		#endif
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Unweighted", "getValues");
		exit(1);
	}
}
/**************************************************************************************************/

EstOutput Unweighted::createProcesses(Tree* t, vector< vector<string> > namesOfGroupCombos) {
	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 1;
		int num = 0;
		vector<int> processIDS;
		
		EstOutput results;
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				EstOutput myresults;
				myresults = driver(t, namesOfGroupCombos, lines[process]->start, lines[process]->num);
				
				if (m->control_pressed) { exit(0); }
				
				m->mothurOut("Merging results."); m->mothurOutEndLine();
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = outputDir + toString(getpid()) + ".unweighted.results.temp";
				m->openOutputFile(tempFile, out);
				out << myresults.size() << endl;
				for (int i = 0; i < myresults.size(); i++) {  out << myresults[i] << '\t';  } out << endl;
				out.close();
				
				exit(0);
			}else { m->mothurOut("unable to spawn the necessary processes."); m->mothurOutEndLine(); exit(0); }
		}
		
		results = driver(t, namesOfGroupCombos, lines[0]->start, lines[0]->num);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<(processors-1);i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		if (m->control_pressed) { return results; }
		
		//get data created by processes
		for (int i=0;i<(processors-1);i++) { 
			ifstream in;
			string s = outputDir + toString(processIDS[i]) + ".unweighted.results.temp";
			m->openInputFile(s, in);
			
			//get quantiles
			if (!in.eof()) {
				int num;
				in >> num; m->gobble(in);
				
				if (m->control_pressed) { break; }
				
				double w; 
				for (int j = 0; j < num; j++) {
					in >> w;
					results.push_back(w);
				}
				m->gobble(in);
			}
			in.close();
			remove(s.c_str());
		}
		
		m->mothurOut("DONE."); m->mothurOutEndLine(); m->mothurOutEndLine();
		
		return results;
#endif		
	}
	catch(exception& e) {
		m->errorOut(e, "Unweighted", "createProcesses");
		exit(1);
	}
}
/**************************************************************************************************/
EstOutput Unweighted::driver(Tree* t, vector< vector<string> > namesOfGroupCombos, int start, int num) { 
 try {
		
		EstOutput results; results.resize(num);
		
		int count = 0;
		int total = num;
		int twentyPercent = (total * 0.20);

		for (int h = start; h < (start+num); h++) {
		
			if (m->control_pressed) { return results; }
		
			double UniqueBL=0.0000;  //a branch length is unique if it's chidren are from the same group
			double totalBL = 0.00;	//all branch lengths
			double UW = 0.00;		//Unweighted Value = UniqueBL / totalBL;
				
			for(int i=0;i<t->getNumNodes();i++){
			
				if (m->control_pressed) {  return data; }
				
				//pcountSize = 0, they are from a branch that is entirely from a group the user doesn't want
				//pcountSize = 2, not unique to one group
				//pcountSize = 1, unique to one group
				
				int pcountSize = 0;
				for (int j = 0; j < namesOfGroupCombos[h].size(); j++) {
					map<string, int>::iterator itGroup = t->tree[i].pcount.find(namesOfGroupCombos[h][j]);
					if (itGroup != t->tree[i].pcount.end()) { pcountSize++; if (pcountSize > 1) { break; } } 
				}
				
				if (pcountSize == 0) { }
				else if ((t->tree[i].getBranchLength() != -1) && (pcountSize == 1)) {  UniqueBL += abs(t->tree[i].getBranchLength());	}
					
				if ((t->tree[i].getBranchLength() != -1) && (pcountSize != 0)) {  
					totalBL += abs(t->tree[i].getBranchLength()); 
				}
			}
		
			UW = (UniqueBL / totalBL);  
	
			if (isnan(UW) || isinf(UW)) { UW = 0; }
	
			results[count] = UW;
			count++;

			//report progress
			if((count) % twentyPercent == 0){	m->mothurOut("Percentage complete: " + toString(int((count / (float)total) * 100.0))); m->mothurOutEndLine();		}
		}
		
		m->mothurOut("Percentage complete: 100"); m->mothurOutEndLine();
		
		return results; 
	}
	catch(exception& e) {
		m->errorOut(e, "Unweighted", "driver");
		exit(1);
	}
}
/**************************************************************************************************/

EstOutput Unweighted::getValues(Tree* t, string groupA, string groupB) { 
 try {
	globaldata = GlobalData::getInstance();
		
		vector<string> groups;
		double UniqueBL;  //a branch length is unique if it's chidren are from the same group
		double totalBL;	//all branch lengths
		double UW;		//Unweighted Value = UniqueBL / totalBL;
		copyTree = new Tree;

		//if the users enters no groups then give them the score of all groups
		int numGroups = globaldata->Groups.size();

		//calculate number of comparsions
		int numComp = 0;
		for (int r=0; r<numGroups; r++) { 
			for (int l = r+1; l < numGroups; l++) {
				numComp++;
			}
		}

		//numComp+1 for AB, AC, BC, ABC
		data.resize(numComp+1,0);
		
		int count = 0;
		for (int a=0; a<numGroups; a++) { 
			for (int l = a+1; l < numGroups; l++) {
				UniqueBL=0.0000;  //a branch length is unique if it's chidren are from the same group
				totalBL = 0.00;	//all branch lengths
				UW = 0.00;		//Unweighted Value = UniqueBL / totalBL;
				
				//copy random tree passed in
				copyTree->getCopy(t);
								
				//groups in this combo
				groups.push_back(globaldata->Groups[a]); groups.push_back(globaldata->Groups[l]);
				
				//swap labels in the groups you want to compare
				copyTree->assembleRandomUnifracTree(groups[0], groups[1]);
				
				if (m->control_pressed) { delete copyTree; return data; }
				
				for(int i=0;i<copyTree->getNumNodes();i++){
			
					if (m->control_pressed) {  return data; }
					
					//pcountSize = 0, they are from a branch that is entirely from a group the user doesn't want
					//pcountSize = 2, not unique to one group
					//pcountSize = 1, unique to one group
					
					int pcountSize = 0;
					for (int j = 0; j < groups.size(); j++) {
						map<string, int>::iterator itGroup = copyTree->tree[i].pcount.find(groups[j]);
						if (itGroup != copyTree->tree[i].pcount.end()) { pcountSize++; } 
					}
					
					if (pcountSize == 0) { }
					else if ((copyTree->tree[i].getBranchLength() != -1) && (pcountSize == 1)) {  UniqueBL += abs(copyTree->tree[i].getBranchLength());	}
						
					if ((copyTree->tree[i].getBranchLength() != -1) && (pcountSize != 0)) {  
						totalBL += abs(copyTree->tree[i].getBranchLength()); 
					}
				}

				
				UW = (UniqueBL / totalBL);  
	
				if (isnan(UW) || isinf(UW)) { UW = 0; }
	
				data[count] = UW;
				count++;
				groups.clear();
			}
		}
		
		
		if (numComp != 1) {
			if (numGroups == 0) {
				//get score for all users groups
				for (int i = 0; i < tmap->namesOfGroups.size(); i++) {
					if (tmap->namesOfGroups[i] != "xxx") {
						groups.push_back(tmap->namesOfGroups[i]);
					}
				}
			}else {
				for (int i = 0; i < globaldata->Groups.size(); i++) {
					groups.push_back(globaldata->Groups[i]);
				}
			}
		
			UniqueBL=0.0000;  //a branch length is unique if it's chidren are from the same group
			totalBL = 0.00;	//all branch lengths
			UW = 0.00;		//Unweighted Value = UniqueBL / totalBL;
		
			//copy random tree passed in
			copyTree->getCopy(t);
				
			//swap labels in all the groups you want to compare
			copyTree->assembleRandomUnifracTree(groups);
			
			if (m->control_pressed) { delete copyTree; return data; }

			for(int i=0;i<copyTree->getNumNodes();i++){
			
				if (m->control_pressed) {  return data; }
				
				//pcountSize = 0, they are from a branch that is entirely from a group the user doesn't want
				//pcountSize = 2, not unique to one group
				//pcountSize = 1, unique to one group
				
				int pcountSize = 0;
				for (int j = 0; j < groups.size(); j++) {
					map<string, int>::iterator itGroup = copyTree->tree[i].pcount.find(groups[j]);
					if (itGroup != copyTree->tree[i].pcount.end()) { pcountSize++; if (pcountSize > 1) { break; } } 
				}
				
				if (pcountSize == 0) { }
				else if ((copyTree->tree[i].getBranchLength() != -1) && (pcountSize == 1)) {  UniqueBL += abs(copyTree->tree[i].getBranchLength());	}
					
				if ((copyTree->tree[i].getBranchLength() != -1) && (pcountSize != 0)) {  
					totalBL += abs(copyTree->tree[i].getBranchLength()); 
				}
			}
		
			UW = (UniqueBL / totalBL);  
	
			if (isnan(UW) || isinf(UW)) { UW = 0; }
	
			data[count] = UW;
		}
		
		delete copyTree;
		
		return data;
	
	}
	catch(exception& e) {
		m->errorOut(e, "Unweighted", "getValues");
		exit(1);
	}
}

/**************************************************************************************************/


