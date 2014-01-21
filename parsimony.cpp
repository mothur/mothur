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
        CountTable* ct = t->getCountTable();
		
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
				vector<string> tGroups = ct->getNamesOfGroups();
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
        
        lines.clear();
        int remainingPairs = namesOfGroupCombos.size();
        int startIndex = 0;
        for (int remainingProcessors = processors; remainingProcessors > 0; remainingProcessors--) {
            int numPairs = remainingPairs; //case for last processor
            if (remainingProcessors != 1) { numPairs = ceil(remainingPairs / remainingProcessors); }
            lines.push_back(linePair(startIndex, numPairs)); //startIndex, numPairs
            startIndex = startIndex + numPairs;
            remainingPairs = remainingPairs - numPairs;
        }
        
        data = createProcesses(t, namesOfGroupCombos, ct);
		
		return data;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Parsimony", "getValues");
		exit(1);
	}
}
/**************************************************************************************************/

EstOutput Parsimony::createProcesses(Tree* t, vector< vector<string> > namesOfGroupCombos, CountTable* ct) {
	try {
        int process = 1;
		vector<int> processIDS;
		
		EstOutput results;

#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				EstOutput myresults;
				myresults = driver(t, namesOfGroupCombos, lines[process].start, lines[process].num, ct);
				
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
		
		results = driver(t, namesOfGroupCombos, lines[0].start, lines[0].num, ct);
		
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
#else
        //fill in functions
        vector<parsData*> pDataArray;
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1];
        vector<CountTable*> cts;
        vector<Tree*> trees;
		
		//Create processor worker threads.
		for( int i=1; i<processors; i++ ){
            CountTable* copyCount = new CountTable();
            copyCount->copy(ct);
            Tree* copyTree = new Tree(copyCount);
            copyTree->getCopy(t);
            
            cts.push_back(copyCount);
            trees.push_back(copyTree);
            
            parsData* temppars = new parsData(m, lines[i].start, lines[i].num, namesOfGroupCombos, copyTree, copyCount);
			pDataArray.push_back(temppars);
			processIDS.push_back(i);
            
			hThreadArray[i-1] = CreateThread(NULL, 0, MyParsimonyThreadFunction, pDataArray[i-1], 0, &dwThreadIdArray[i-1]);
		}
		
		results = driver(t, namesOfGroupCombos, lines[0].start, lines[0].num, ct);
		
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
            for (int j = 0; j < pDataArray[i]->results.size(); j++) {  results.push_back(pDataArray[i]->results[j]);  }
			delete cts[i];
            delete trees[i];
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
		
#endif		
        return results;
	}
	catch(exception& e) {
		m->errorOut(e, "Parsimony", "createProcesses");
		exit(1);
	}
}
/**************************************************************************************************/
EstOutput Parsimony::driver(Tree* t, vector< vector<string> > namesOfGroupCombos, int start, int num, CountTable* ct) { 
	try {
		
		EstOutput results; results.resize(num);
		
		Tree* copyTree = new Tree(ct);
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

