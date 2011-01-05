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
			for (int l = 0; l < r; l++) {
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
		
					lines.push_back(linePair(startPos, numPairsPerProcessor));
				}

				data = createProcesses(t, namesOfGroupCombos);
				lines.clear();
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
				myresults = driver(t, namesOfGroupCombos, lines[process].start, lines[process].num);
				
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
			}else { 
				m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
				for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
				exit(0); 
			}
		}
		
		results = driver(t, namesOfGroupCombos, lines[0].start, lines[0].num);
		
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
		int numLeaves = t->getNumLeaves();
		int total = num;
		int twentyPercent = (total * 0.20);
		if (twentyPercent == 0) { twentyPercent = 1; }
		
		
		for (int h = start; h < (start+num); h++) {
	//cout << namesOfGroupCombos[h][0] << '\t' << namesOfGroupCombos[h][1] << endl;		
			if (m->control_pressed) { return results; }
		
			double UniqueBL=0.0000;  //a branch length is unique if it's chidren are from the same group
			double totalBL = 0.00;	//all branch lengths
			double UW = 0.00;		//Unweighted Value = UniqueBL / totalBL;
			map<int, double> tempTotals; //maps node to total Branch Length
			map<int, int> nodePcountSize; //maps node to pcountSize
			map<int, int>::iterator itCount;
				
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

				nodePcountSize[i] = pcountSize;
				
				//cout << i << '\t' << t->tree[i].getName() << " br = " << abs(t->tree[i].getBranchLength()) << '\t';		
				if (pcountSize == 0) { }
				else if ((t->tree[i].getBranchLength() != -1) && (pcountSize == 1)) {  UniqueBL += abs(t->tree[i].getBranchLength()); }
				
				//if you are a leaf from a users group add to total
				if (i < numLeaves) {
					if ((t->tree[i].getBranchLength() != -1) && pcountSize != 0) { 
					//cout << "added to total" << endl; 
						totalBL += abs(t->tree[i].getBranchLength()); 
					}
					tempTotals[i] = 0.0;  //we don't care about you, or we have already added you
				}else{ //if you are not a leaf 
					//do both your chidren have have descendants from the users groups? 
					int lc = t->tree[i].getLChild();
					int rc = t->tree[i].getRChild();
					
					//if yes, add your childrens tempTotals
					if ((nodePcountSize[lc] != 0) && (nodePcountSize[rc] != 0)) {
						totalBL += tempTotals[lc] + tempTotals[rc]; 
						//cout << "added to total " << tempTotals[lc] << '\t' << tempTotals[rc] << endl;
						if (t->tree[i].getBranchLength() != -1) {
							tempTotals[i] = abs(t->tree[i].getBranchLength());
						}else {
							tempTotals[i] = 0.0;
						}
					}else if ((nodePcountSize[lc] == 0) && (nodePcountSize[rc] == 0)) { tempTotals[i] = 0.0;  //we don't care about you
					}else { //if no, your tempTotal is your childrens temp totals + your branch length
						tempTotals[i] = tempTotals[lc] + tempTotals[rc] + abs(t->tree[i].getBranchLength()); 
					}
					//cout << "temptotal = "<< tempTotals[i] << endl;
				}

			}
	//cout << UniqueBL << '\t' << totalBL << endl;		
			UW = (UniqueBL / totalBL);  
	
			if (isnan(UW) || isinf(UW)) { UW = 0; }
	
			results[count] = UW;
			count++;

			//report progress
			if((count % twentyPercent) == 0) {	float tempOut = (count / (float)total); if (isnan(tempOut) || isinf(tempOut)) { tempOut = 0.0; } m->mothurOut("Percentage complete: " + toString((int(tempOut) * 100.0))); m->mothurOutEndLine();	}
		}
		
		//report progress
		if((count % twentyPercent) != 0) {	float tempOut = (count / (float)total); if (isnan(tempOut) || isinf(tempOut)) { tempOut = 0.0; } m->mothurOut("Percentage complete: " + toString((int(tempOut) * 100.0))); m->mothurOutEndLine();	}
		
		return results; 
	}
	catch(exception& e) {
		m->errorOut(e, "Unweighted", "driver");
		exit(1);
	}
}
/**************************************************************************************************/

EstOutput Unweighted::getValues(Tree* t, string groupA, string groupB, int p, string o) { 
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
			for (int l = 0; l < r; l++) {
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
				data = driver(t, namesOfGroupCombos, 0, namesOfGroupCombos.size(), true);
			}else{
				int numPairs = namesOfGroupCombos.size();
				
				int numPairsPerProcessor = numPairs / processors;
				
				for (int i = 0; i < processors; i++) {
					int startPos = i * numPairsPerProcessor;
					if(i == processors - 1){
						numPairsPerProcessor = numPairs - i * numPairsPerProcessor;
					}
					lines.push_back(linePair(startPos, numPairsPerProcessor));
				}

				data = createProcesses(t, namesOfGroupCombos, true);
				
				lines.clear();
			}
		#else
			data = driver(t, namesOfGroupCombos, 0, namesOfGroupCombos.size(), true);
		#endif
	
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Unweighted", "getValues");
		exit(1);
	}
}
/**************************************************************************************************/

EstOutput Unweighted::createProcesses(Tree* t, vector< vector<string> > namesOfGroupCombos, bool usingGroups) {
	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 1;
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
				myresults = driver(t, namesOfGroupCombos, lines[process].start, lines[process].num, usingGroups);
				
				if (m->control_pressed) { exit(0); }
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = outputDir + toString(getpid()) + ".unweighted.results.temp";
				m->openOutputFile(tempFile, out);
				out << myresults.size() << endl;
				for (int i = 0; i < myresults.size(); i++) {  out << myresults[i] << '\t';  } out << endl;
				out.close();
				
				exit(0);
			}else { 
				m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
				for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
				exit(0); 
			}
		}
		
		results = driver(t, namesOfGroupCombos, lines[0].start, lines[0].num, usingGroups);
		
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
		
		return results;
#endif		
	}
	catch(exception& e) {
		m->errorOut(e, "Unweighted", "createProcesses");
		exit(1);
	}
}
/**************************************************************************************************/
EstOutput Unweighted::driver(Tree* t, vector< vector<string> > namesOfGroupCombos, int start, int num, bool usingGroups) { 
 try {
		
		EstOutput results; results.resize(num);
		
		int count = 0;
		int numLeaves = t->getNumLeaves();
		
		Tree* copyTree = new Tree;
		
		for (int h = start; h < (start+num); h++) {
		
			if (m->control_pressed) { return results; }
		
			//copy random tree passed in
			copyTree->getCopy(t);
				
			//swap labels in the groups you want to compare
			copyTree->assembleRandomUnifracTree(namesOfGroupCombos[h]);
			
			double UniqueBL=0.0000;  //a branch length is unique if it's chidren are from the same group
			double totalBL = 0.00;	//all branch lengths
			double UW = 0.00;		//Unweighted Value = UniqueBL / totalBL;
			map<int, double> tempTotals; //maps node to total Branch Length
			map<int, int> nodePcountSize; //maps node to pcountSize
				
			for(int i=0;i<copyTree->getNumNodes();i++){
			
				if (m->control_pressed) {  return data; }
				
				//pcountSize = 0, they are from a branch that is entirely from a group the user doesn't want
				//pcountSize = 2, not unique to one group
				//pcountSize = 1, unique to one group
				
				int pcountSize = 0;
				for (int j = 0; j < namesOfGroupCombos[h].size(); j++) {
					map<string, int>::iterator itGroup = copyTree->tree[i].pcount.find(namesOfGroupCombos[h][j]);
					if (itGroup != copyTree->tree[i].pcount.end()) { pcountSize++; if (pcountSize > 1) { break; } } 
				}
				
				nodePcountSize[i] = pcountSize;
			
				if (pcountSize == 0) { }
				else if ((copyTree->tree[i].getBranchLength() != -1) && (pcountSize == 1)) {  UniqueBL += abs(copyTree->tree[i].getBranchLength());	 }
				
				//if you are a leaf from a users group add to total
				if (i < numLeaves) {
					if ((copyTree->tree[i].getBranchLength() != -1) && pcountSize != 0) { 
						totalBL += abs(copyTree->tree[i].getBranchLength()); 
					}
					tempTotals[i] = 0.0;  //we don't care about you, or we have already added you
				}else{ //if you are not a leaf 
					//do both your chidren have have descendants from the users groups? 
					int lc = copyTree->tree[i].getLChild();
					int rc = copyTree->tree[i].getRChild();
					
					//if yes, add your childrens tempTotals
					if ((nodePcountSize[lc] != 0) && (nodePcountSize[rc] != 0)) {
						totalBL += tempTotals[lc] + tempTotals[rc]; 
						
						if (copyTree->tree[i].getBranchLength() != -1) {
							tempTotals[i] = abs(copyTree->tree[i].getBranchLength());
						}else {
							tempTotals[i] = 0.0;
						}
					}else if ((nodePcountSize[lc] == 0) && (nodePcountSize[rc] == 0)) { tempTotals[i] = 0.0;  //we don't care about you
					}else { //if no, your tempTotal is your childrens temp totals + your branch length
						tempTotals[i] = tempTotals[lc] + tempTotals[rc] + abs(copyTree->tree[i].getBranchLength()); 
					}
					
				}

			}
		
			UW = (UniqueBL / totalBL);  
	
			if (isnan(UW) || isinf(UW)) { UW = 0; }
	
			results[count] = UW;
			count++;

		}
		
		delete copyTree;
		
		return results; 
	}
	catch(exception& e) {
		m->errorOut(e, "Unweighted", "driver");
		exit(1);
	}
}
/**************************************************************************************************/


