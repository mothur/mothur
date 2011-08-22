/*
 *  parsimony.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/26/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "parsimony.h"

/**************************************************************************************************/

EstOutput Parsimony::getValues(Tree* t, int p, string o) {
	try {
		processors = p;
		outputDir = o;
		
		//if the users enters no groups then give them the score of all groups
		vector<string> mGroups = m->getGroups();
		int numGroups = mGroups.size();
		
		//calculate number of comparsions
		int numComp = 0;
		vector< vector<string> > namesOfGroupCombos;
		for (int r=0; r<numGroups; r++) { 
			for (int l = 0; l < r; l++) {
				numComp++;
				vector<string> groups; groups.push_back(mGroups[r]); groups.push_back(mGroups[l]);
				//cout << globaldata->Groups[r] << '\t' << globaldata->Groups[l] << endl;
				namesOfGroupCombos.push_back(groups);
			}
		}

		//numComp+1 for AB, AC, BC, ABC
		if (numComp != 1) {
			vector<string> groups;
			if (numGroups == 0) {
				//get score for all users groups
				vector<string> tGroups = tmap->getNamesOfGroups();
				for (int i = 0; i < tGroups.size(); i++) {
					if (tGroups[i] != "xxx") {
						groups.push_back(tGroups[i]);
						//cout << tmap->namesOfGroups[i] << endl;
					}
				}
				namesOfGroupCombos.push_back(groups);
			}else {
				for (int i = 0; i < mGroups.size(); i++) {
					groups.push_back(mGroups[i]);
					//cout << globaldata->Groups[i] << endl;
				}
				namesOfGroupCombos.push_back(groups);
			}
		}
		
	#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		if(processors == 1){
			data = driver(t, namesOfGroupCombos, 0, namesOfGroupCombos.size());
		}else{
			lines.clear();
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
		}
	#else
		data = driver(t, namesOfGroupCombos, 0, namesOfGroupCombos.size());
	#endif
		
		return data;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Parsimony", "getValues");
		exit(1);
	}
}
/**************************************************************************************************/

EstOutput Parsimony::createProcesses(Tree* t, vector< vector<string> > namesOfGroupCombos) {
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
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = outputDir + toString(getpid()) + ".pars.results.temp";
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
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		if (m->control_pressed) { return results; }
			
		//get data created by processes
		for (int i=0;i<processIDS.size();i++) { 
			ifstream in;
			string s = outputDir + toString(processIDS[i]) + ".pars.results.temp";
			m->openInputFile(s, in);
			
			//get scores
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
			m->mothurRemove(s);
		}
		
		return results;
#endif		
	}
	catch(exception& e) {
		m->errorOut(e, "Parsimony", "createProcesses");
		exit(1);
	}
}
/**************************************************************************************************/
EstOutput Parsimony::driver(Tree* t, vector< vector<string> > namesOfGroupCombos, int start, int num) { 
	try {
		
		EstOutput results; results.resize(num);
		
		Tree* copyTree = new Tree(tmap);
		int count = 0;
		
		for (int h = start; h < (start+num); h++) {
					
			if (m->control_pressed) { delete copyTree; return results; }
	
			int score = 0;
			
			//groups in this combo
			vector<string> groups = namesOfGroupCombos[h];
			
			//copy users tree so that you can redo pgroups 
			copyTree->getCopy(t);
			
			//create pgroups that reflect the groups the user want to use
			for(int i=copyTree->getNumLeaves();i<copyTree->getNumNodes();i++){
				copyTree->tree[i].pGroups = (copyTree->mergeUserGroups(i, groups));
			}
			
			for(int i=copyTree->getNumLeaves();i<copyTree->getNumNodes();i++){
				
				if (m->control_pressed) { return data; }
				
				int lc = copyTree->tree[i].getLChild();
				int rc = copyTree->tree[i].getRChild();
				
				int iSize = copyTree->tree[i].pGroups.size();
				int rcSize = copyTree->tree[rc].pGroups.size();
				int lcSize = copyTree->tree[lc].pGroups.size();
				
				//if isize are 0 then that branch is to be ignored
				if (iSize == 0) { }
				else if ((rcSize == 0) || (lcSize == 0)) { }
				//if you have more groups than either of your kids then theres been a change.
				else if(iSize > rcSize || iSize > lcSize){
					score++;
				}
			} 
			
			results[count] = score;
			count++;
		}
					
		delete copyTree;
			
		return results; 
	}
	catch(exception& e) {
		m->errorOut(e, "Parsimony", "driver");
		exit(1);
	}
}

/**************************************************************************************************/

