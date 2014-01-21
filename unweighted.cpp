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
		processors = p;
		outputDir = o;
        
        CountTable* ct = t->getCountTable();
        
		//if the users enters no groups then give them the score of all groups
		int numGroups = m->getNumGroups();
		
		//calculate number of comparsions
		int numComp = 0;
		vector< vector<string> > namesOfGroupCombos;
		for (int r=0; r<numGroups; r++) { 
			for (int l = 0; l < r; l++) {
				numComp++;
				vector<string> groups; groups.push_back((m->getGroups())[r]); groups.push_back((m->getGroups())[l]);
				namesOfGroupCombos.push_back(groups);
			}
		}
		
		if (numComp != 1) {
			vector<string> groups;
			if (numGroups == 0) {
				//get score for all users groups
				for (int i = 0; i < (ct->getNamesOfGroups()).size(); i++) {
					if ((ct->getNamesOfGroups())[i] != "xxx") {
						groups.push_back((ct->getNamesOfGroups())[i]);
					}
				}
				namesOfGroupCombos.push_back(groups);
			}else {
				for (int i = 0; i < m->getNumGroups(); i++) {
					groups.push_back((m->getGroups())[i]);
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

        lines.clear();
        
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Unweighted", "getValues");
		exit(1);
	}
}
/**************************************************************************************************/

EstOutput Unweighted::createProcesses(Tree* t, vector< vector<string> > namesOfGroupCombos, CountTable* ct) {
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
				
				//m->mothurOut("Merging results."); m->mothurOutEndLine();
				
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
		
		results = driver(t, namesOfGroupCombos, lines[0].start, lines[0].num, ct);
		
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
			m->mothurRemove(s);
		}
#else
		//fill in functions
        vector<unweightedData*> pDataArray;
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
            
            unweightedData* tempweighted = new unweightedData(m, lines[i].start, lines[i].num, namesOfGroupCombos, copyTree, copyCount, includeRoot);
			pDataArray.push_back(tempweighted);
			processIDS.push_back(i);
            
			hThreadArray[i-1] = CreateThread(NULL, 0, MyUnWeightedThreadFunction, pDataArray[i-1], 0, &dwThreadIdArray[i-1]);
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
		m->errorOut(e, "Unweighted", "createProcesses");
		exit(1);
	}
}
/**************************************************************************************************/
EstOutput Unweighted::driver(Tree* t, vector< vector<string> > namesOfGroupCombos, int start, int num, CountTable* ct) { 
 try {
	
	 
		EstOutput results; results.resize(num);
		
		int count = 0;
		int total = num;
					
		for (int h = start; h < (start+num); h++) {
				
			if (m->control_pressed) { return results; }
		
			double UniqueBL=0.0000;  //a branch length is unique if it's chidren are from the same group
			double totalBL = 0.00;	//all branch lengths
			double UW = 0.00;		//Unweighted Value = UniqueBL / totalBL;
						
			//find a node that belongs to one of the groups in this combo
			int nodeBelonging = -1;
			for (int g = 0; g < namesOfGroupCombos[h].size(); g++) {
				if (t->groupNodeInfo[namesOfGroupCombos[h][g]].size() != 0) { nodeBelonging = t->groupNodeInfo[namesOfGroupCombos[h][g]][0]; break; }
			}
			
			//sanity check
			if (nodeBelonging == -1) {
				m->mothurOut("[WARNING]: cannot find a nodes in the tree from grouping "); 
				for (int g = 0; g < namesOfGroupCombos[h].size()-1; g++) { m->mothurOut(namesOfGroupCombos[h][g] + "-"); }
				m->mothurOut(namesOfGroupCombos[h][namesOfGroupCombos[h].size()-1]);
				m->mothurOut(", skipping."); m->mothurOutEndLine(); results[count] = UW;
			}else{
				//cout << "trying to get root" << endl;	
				//if including the root this clears rootForGrouping[namesOfGroupCombos[h]]
				getRoot(t, nodeBelonging, namesOfGroupCombos[h]);
				//cout << "here" << endl;	
				for(int i=0;i<t->getNumNodes();i++){
					
					if (m->control_pressed) {  return data; }
					//cout << i << endl;	
					//pcountSize = 0, they are from a branch that is entirely from a group the user doesn't want
					//pcountSize = 2, not unique to one group
					//pcountSize = 1, unique to one group
					
					int pcountSize = 0;
					for (int j = 0; j < namesOfGroupCombos[h].size(); j++) {
						map<string, int>::iterator itGroup = t->tree[i].pcount.find(namesOfGroupCombos[h][j]);
						if (itGroup != t->tree[i].pcount.end()) { pcountSize++; if (pcountSize > 1) { break; } } 
					}
					
					
					//unique calc
					if (pcountSize == 0) { }
					else if ((t->tree[i].getBranchLength() != -1) && (pcountSize == 1) && (rootForGrouping[namesOfGroupCombos[h]].count(i) == 0)) { //you have a unique branch length and you are not the root 
						UniqueBL += abs(t->tree[i].getBranchLength()); 
					}
						
					//total calc
					if (pcountSize == 0) { }
					else if ((t->tree[i].getBranchLength() != -1) && (pcountSize != 0) && (rootForGrouping[namesOfGroupCombos[h]].count(i) == 0)) { //you have a branch length and you are not the root 
						totalBL += abs(t->tree[i].getBranchLength()); 
					}
				}
	//cout << UniqueBL << '\t' << totalBL << endl;		
				UW = (UniqueBL / totalBL);  
	
				if (isnan(UW) || isinf(UW)) { UW = 0; }
	
				results[count] = UW;
			}
			count++;

        }
		
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
		processors = p;
		outputDir = o;
		
        CountTable* ct = t->getCountTable();
     
		//if the users enters no groups then give them the score of all groups
		int numGroups = m->getNumGroups();
		
		//calculate number of comparsions
		int numComp = 0;
		vector< vector<string> > namesOfGroupCombos;
		for (int r=0; r<numGroups; r++) { 
			for (int l = 0; l < r; l++) {
				numComp++;
				vector<string> groups; groups.push_back((m->getGroups())[r]); groups.push_back((m->getGroups())[l]);
				namesOfGroupCombos.push_back(groups);
			}
		}
		
		if (numComp != 1) {
			vector<string> groups;
			if (numGroups == 0) {
				//get score for all users groups
				for (int i = 0; i < (ct->getNamesOfGroups()).size(); i++) {
					if ((ct->getNamesOfGroups())[i] != "xxx") {
						groups.push_back((ct->getNamesOfGroups())[i]);
					}
				}
				namesOfGroupCombos.push_back(groups);
			}else {
				for (int i = 0; i < m->getNumGroups(); i++) {
					groups.push_back((m->getGroups())[i]);
				}
				namesOfGroupCombos.push_back(groups);
			}
		}
     
        lines.clear();
        int numPairs = namesOfGroupCombos.size();
        int numPairsPerProcessor = ceil(numPairs / processors);
     
        for (int i = 0; i < processors; i++) {
            int startPos = i * numPairsPerProcessor;
            if(i == processors - 1){ numPairsPerProcessor = numPairs - i * numPairsPerProcessor; }
            lines.push_back(linePair(startPos, numPairsPerProcessor));
        }
     
        data = createProcesses(t, namesOfGroupCombos, true, ct);
        lines.clear();

		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Unweighted", "getValues");
		exit(1);
	}
}
/**************************************************************************************************/

EstOutput Unweighted::createProcesses(Tree* t, vector< vector<string> > namesOfGroupCombos, bool usingGroups, CountTable* ct) {
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
				myresults = driver(t, namesOfGroupCombos, lines[process].start, lines[process].num, usingGroups, ct);
				
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
		
		results = driver(t, namesOfGroupCombos, lines[0].start, lines[0].num, usingGroups, ct);
		
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
			m->mothurRemove(s);
		}
#else
       //for some reason it doesn't seem to be calculating hte random trees scores.  all scores are the same even though copytree appears to be randomized.
        
        /*
        //fill in functions
        vector<unweightedData*> pDataArray;
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
            
            unweightedData* tempweighted = new unweightedData(m, lines[i].start, lines[i].num, namesOfGroupCombos, copyTree, copyCount, includeRoot);
			pDataArray.push_back(tempweighted);
			processIDS.push_back(i);
            
			hThreadArray[i-1] = CreateThread(NULL, 0, MyUnWeightedRandomThreadFunction, pDataArray[i-1], 0, &dwThreadIdArray[i-1]);
		}
		
		results = driver(t, namesOfGroupCombos, lines[0].start, lines[0].num, usingGroups, ct);
		
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
            for (int j = 0; j < pDataArray[i]->results.size(); j++) {  results.push_back(pDataArray[i]->results[j]);  }
			delete cts[i];
            delete trees[i];
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}	*/
        
        results = driver(t, namesOfGroupCombos, 0, namesOfGroupCombos.size(), usingGroups, ct);
#endif
        return results;
	}
	catch(exception& e) {
		m->errorOut(e, "Unweighted", "createProcesses");
		exit(1);
	}
}
/**************************************************************************************************/
EstOutput Unweighted::driver(Tree* t, vector< vector<string> > namesOfGroupCombos, int start, int num, bool usingGroups, CountTable* ct) { 
 try {
		
		EstOutput results; results.resize(num);
		
		int count = 0;
		
		Tree* copyTree = new Tree(ct);
		
		for (int h = start; h < (start+num); h++) {
		
			if (m->control_pressed) { return results; }
		
			//copy random tree passed in
			copyTree->getCopy(t);
				
			//swap labels in the groups you want to compare
			copyTree->assembleRandomUnifracTree(namesOfGroupCombos[h]);
			
			double UniqueBL=0.0000;  //a branch length is unique if it's chidren are from the same group
			double totalBL = 0.00;	//all branch lengths
			double UW = 0.00;		//Unweighted Value = UniqueBL / totalBL;
			//find a node that belongs to one of the groups in this combo
			int nodeBelonging = -1;
			for (int g = 0; g < namesOfGroupCombos[h].size(); g++) {
				if (copyTree->groupNodeInfo[namesOfGroupCombos[h][g]].size() != 0) { nodeBelonging = copyTree->groupNodeInfo[namesOfGroupCombos[h][g]][0]; break; }
			}
			
			//sanity check
			if (nodeBelonging == -1) {
				m->mothurOut("[WARNING]: cannot find a nodes in the tree from grouping "); 
				for (int g = 0; g < namesOfGroupCombos[h].size()-1; g++) { m->mothurOut(namesOfGroupCombos[h][g] + "-"); }
				m->mothurOut(namesOfGroupCombos[h][namesOfGroupCombos[h].size()-1]);
				m->mothurOut(", skipping."); m->mothurOutEndLine(); results[count] = UW;
			}else{
				
				//if including the root this clears rootForGrouping[namesOfGroupCombos[h]]
				getRoot(copyTree, nodeBelonging, namesOfGroupCombos[h]);
				
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
					
					//unique calc
					if (pcountSize == 0) { }
					else if ((copyTree->tree[i].getBranchLength() != -1) && (pcountSize == 1) && (rootForGrouping[namesOfGroupCombos[h]].count(i) == 0)) { //you have a unique branch length and you are not the root 
						UniqueBL += abs(copyTree->tree[i].getBranchLength()); 
					}
					
					//total calc
					if (pcountSize == 0) { }
					else if ((copyTree->tree[i].getBranchLength() != -1) && (pcountSize != 0) && (rootForGrouping[namesOfGroupCombos[h]].count(i) == 0)) { //you have a branch length and you are not the root 
						totalBL += abs(copyTree->tree[i].getBranchLength()); 
					}
					
				}
				//cout << UniqueBL << '\t' << totalBL << endl;		
				UW = (UniqueBL / totalBL);  
				
				if (isnan(UW) || isinf(UW)) { UW = 0; }
				
				results[count] = UW;
            }
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
int Unweighted::getRoot(Tree* t, int v, vector<string> grouping) { 
	try {
		//you are a leaf so get your parent
		int index = t->tree[v].getParent();
		
		if (includeRoot) { 
			rootForGrouping[grouping].clear();
		}else {
			
			//my parent is a potential root
			rootForGrouping[grouping].insert(index);
			
			//while you aren't at root
			while(t->tree[index].getParent() != -1){
				//cout << index << endl;	
				if (m->control_pressed) {  return 0; }
				
				//am I the root for this grouping? if so I want to stop "early"
				//does my sibling have descendants from the users groups? 
				//if so I am not the root
				int parent = t->tree[index].getParent();
				int lc = t->tree[parent].getLChild();
				int rc = t->tree[parent].getRChild();
				
				int sib = lc;
				if (lc == index) { sib = rc; }
				
				map<string, int>::iterator itGroup;
				int pcountSize = 0;
				for (int j = 0; j < grouping.size(); j++) {
					map<string, int>::iterator itGroup = t->tree[sib].pcount.find(grouping[j]);
					if (itGroup != t->tree[sib].pcount.end()) { pcountSize++; if (pcountSize > 1) { break; } } 
				}
				
				//if yes, I am not the root
				if (pcountSize != 0) {
					rootForGrouping[grouping].clear();
					rootForGrouping[grouping].insert(parent);
				}
				
				index = parent;	
			}
			
			//get all nodes above the root to add so we don't add their u values above
			index = *(rootForGrouping[grouping].begin());
			while(t->tree[index].getParent() != -1){
				int parent = t->tree[index].getParent();
				rootForGrouping[grouping].insert(parent);
				//cout << parent << " in root" << endl;
				index = parent;
			}
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "Unweighted", "getRoot");
		exit(1);
	}
}
/**************************************************************************************************/


