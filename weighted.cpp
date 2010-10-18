/*
 *  weighted.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "weighted.h"

/**************************************************************************************************/

EstOutput Weighted::getValues(Tree* t, int p, string o) {
    try {
		globaldata = GlobalData::getInstance();
		
		data.clear(); //clear out old values
		int numGroups;
		vector<double> D;
		processors = p;
		outputDir = o;
		
		numGroups = globaldata->Groups.size();
		
		if (m->control_pressed) { return data; }
		
		//calculate number of comparisons i.e. with groups A,B,C = AB, AC, BC = 3;
		vector< vector<string> > namesOfGroupCombos;
		for (int i=0; i<numGroups; i++) { 
			for (int l = 0; l < i; l++) {	
				//initialize weighted scores
				//WScore[globaldata->Groups[i]+globaldata->Groups[l]] = 0.0;
				vector<string> groups; groups.push_back(globaldata->Groups[i]); groups.push_back(globaldata->Groups[l]);
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
		m->errorOut(e, "Weighted", "getValues");
		exit(1);
	}
}
/**************************************************************************************************/

EstOutput Weighted::createProcesses(Tree* t, vector< vector<string> > namesOfGroupCombos) {
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
	
				EstOutput Myresults;
				Myresults = driver(t, namesOfGroupCombos, lines[process].start, lines[process].num);
			
				m->mothurOut("Merging results."); m->mothurOutEndLine();
				
				//pass numSeqs to parent
				ofstream out;

				string tempFile = outputDir + toString(getpid()) + ".weighted.results.temp";
	
				m->openOutputFile(tempFile, out);
	
				out << Myresults.size() << endl;
				for (int i = 0; i < Myresults.size(); i++) {  out << Myresults[i] << '\t';  } out << endl;
				out.close();
				
				exit(0);
			}else { m->mothurOut("unable to spawn the necessary processes."); m->mothurOutEndLine(); exit(0); }
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
			string s = outputDir + toString(processIDS[i]) + ".weighted.results.temp";
			m->openInputFile(s, in);
			
			//get quantiles
			while (!in.eof()) {
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
		m->errorOut(e, "Weighted", "createProcesses");
		exit(1);
	}
}
/**************************************************************************************************/
EstOutput Weighted::driver(Tree* t, vector< vector<string> > namesOfGroupCombos, int start, int num) { 
 try {
		EstOutput results;
		vector<double> D;
		
		int count = 0;
		for (int h = start; h < (start+num); h++) {
		
			if (m->control_pressed) { return results; }
		
			//initialize weighted score
			string groupA = namesOfGroupCombos[h][0]; 
			string groupB = namesOfGroupCombos[h][1];
			
			set<int> validBranches;
			WScore[groupA+groupB] = 0.0;
			D.push_back(0.0000); //initialize a spot in D for each combination
			
			//adding the wieghted sums from group i
			for (int j = 0; j < t->groupNodeInfo[groupA].size(); j++) { //the leaf nodes that have seqs from group i
				map<string, int>::iterator it = t->tree[t->groupNodeInfo[groupA][j]].pcount.find(groupA);
				int numSeqsInGroupI = it->second;
				
				double sum = getLengthToRoot(t, t->groupNodeInfo[groupA][j], groupA, groupB);
				double weightedSum = ((numSeqsInGroupI * sum) / (double)tmap->seqsPerGroup[groupA]);
			
				D[count] += weightedSum;
			}
			
			//adding the wieghted sums from group l
			for (int j = 0; j < t->groupNodeInfo[groupB].size(); j++) { //the leaf nodes that have seqs from group l
				map<string, int>::iterator it = t->tree[t->groupNodeInfo[groupB][j]].pcount.find(groupB);
				int numSeqsInGroupL = it->second;
				
				double sum = getLengthToRoot(t, t->groupNodeInfo[groupB][j], groupA, groupB);
				double weightedSum = ((numSeqsInGroupL * sum) / (double)tmap->seqsPerGroup[groupB]);
			
				D[count] += weightedSum;
			}
			count++;
		}
		
		//calculate u for the group comb 
		
		for (int h = start; h < (start+num); h++) {	
			//report progress
			m->mothurOut("Processing combo: " + toString(h)); m->mothurOutEndLine();
			
			int numLeaves = t->getNumLeaves();
			map<int, double> tempTotals; //maps node to total Branch Length
			map<int, int> nodePcountSize; //maps node to pcountSize
			
			string groupA = namesOfGroupCombos[h][0]; 
			string groupB = namesOfGroupCombos[h][1];
			
			//calculate u for the group comb 
			for(int i=0;i<t->getNumNodes();i++){
				
				if (m->control_pressed) { return data; }
				
				double u;
				int pcountSize = 0;
				//does this node have descendants from groupA
				it = t->tree[i].pcount.find(groupA);
				//if it does u = # of its descendants with a certain group / total number in tree with a certain group
				if (it != t->tree[i].pcount.end()) {
					u = (double) t->tree[i].pcount[groupA] / (double) tmap->seqsPerGroup[groupA];
					pcountSize++;
				}else { u = 0.00; }
				
				
				//does this node have descendants from group l
				it = t->tree[i].pcount.find(groupB);
				//if it does subtract their percentage from u
				if (it != t->tree[i].pcount.end()) {
					u -= (double) t->tree[i].pcount[groupB] / (double) tmap->seqsPerGroup[groupB];
					pcountSize++;
				}
				
				u = abs(u * t->tree[i].getBranchLength());
				
				nodePcountSize[i] = pcountSize;
				
				//if you are a leaf from a users group add to total
				if (i < numLeaves) {
					if ((t->tree[i].getBranchLength() != -1) && pcountSize != 0) { 
						//cout << "added to total" << endl; 
						WScore[(groupA+groupB)] += u; 
					}
					tempTotals[i] = 0.0;  //we don't care about you, or we have already added you
				}else{ //if you are not a leaf 
					//do both your chidren have have descendants from the users groups? 
					int lc = t->tree[i].getLChild();
					int rc = t->tree[i].getRChild();
					
					//if yes, add your childrens tempTotals
					if ((nodePcountSize[lc] != 0) && (nodePcountSize[rc] != 0)) {
						WScore[(groupA+groupB)] += tempTotals[lc] + tempTotals[rc]; 
						//cout << "added to total " << tempTotals[lc] << '\t' << tempTotals[rc] << endl;
						if (t->tree[i].getBranchLength() != -1) {
							tempTotals[i] = u;
						}else {
							tempTotals[i] = 0.0;
						}
					}else if ((nodePcountSize[lc] == 0) && (nodePcountSize[rc] == 0)) { tempTotals[i] = 0.0;  //we don't care about you
					}else { //if no, your tempTotal is your childrens temp totals + your branch length
						tempTotals[i] = tempTotals[lc] + tempTotals[rc] + u; 
					}
					//cout << "temptotal = "<< tempTotals[i] << endl;
				}
			}
		}
		
		/********************************************************/
		//calculate weighted score for the group combination
		double UN;	
		count = 0;
		for (int h = start; h < (start+num); h++) {
			UN = (WScore[namesOfGroupCombos[h][0]+namesOfGroupCombos[h][1]] / D[count]);
		
			if (isnan(UN) || isinf(UN)) { UN = 0; } 
			results.push_back(UN);
			count++;
		}
				
		return results; 
	}
	catch(exception& e) {
		m->errorOut(e, "Weighted", "driver");
		exit(1);
	}
}
/**************************************************************************************************/
EstOutput Weighted::getValues(Tree* t, string groupA, string groupB) { 
 try {
		globaldata = GlobalData::getInstance();
		
		data.clear(); //clear out old values
		
		if (m->control_pressed) { return data; }
		
		//initialize weighted score
		WScore[(groupA+groupB)] = 0.0;
		double D = 0.0;
		set<int> validBranches;
		
		vector<string> groups; groups.push_back(groupA); groups.push_back(groupB);
		
		//adding the wieghted sums from group i
		for (int j = 0; j < t->groupNodeInfo[groups[0]].size(); j++) { //the leaf nodes that have seqs from group i
			map<string, int>::iterator it = t->tree[t->groupNodeInfo[groups[0]][j]].pcount.find(groups[0]);
			int numSeqsInGroupI = it->second;
			
			double sum = getLengthToRoot(t, t->groupNodeInfo[groups[0]][j], groups[0], groups[1]);
			double weightedSum = ((numSeqsInGroupI * sum) / (double)tmap->seqsPerGroup[groups[0]]);
		
			D += weightedSum;
		}
		
		//adding the wieghted sums from group l
		for (int j = 0; j < t->groupNodeInfo[groups[1]].size(); j++) { //the leaf nodes that have seqs from group l
			map<string, int>::iterator it = t->tree[t->groupNodeInfo[groups[1]][j]].pcount.find(groups[1]);
			int numSeqsInGroupL = it->second;
			
			double sum = getLengthToRoot(t, t->groupNodeInfo[groups[1]][j], groups[0], groups[1]);
			double weightedSum = ((numSeqsInGroupL * sum) / (double)tmap->seqsPerGroup[groups[1]]);
		
			D += weightedSum;
		}
		
		int numLeaves = t->getNumLeaves();
		map<int, double> tempTotals; //maps node to total Branch Length
		map<int, int> nodePcountSize; //maps node to pcountSize
		
		//calculate u for the group comb 
		for(int i=0;i<t->getNumNodes();i++){
		
			if (m->control_pressed) { return data; }
			
			double u;
			int pcountSize = 0;
			//does this node have descendants from groupA
			it = t->tree[i].pcount.find(groupA);
			//if it does u = # of its descendants with a certain group / total number in tree with a certain group
			if (it != t->tree[i].pcount.end()) {
				u = (double) t->tree[i].pcount[groupA] / (double) tmap->seqsPerGroup[groupA];
				pcountSize++;
			}else { u = 0.00; }
		
						
			//does this node have descendants from group l
			it = t->tree[i].pcount.find(groupB);
			//if it does subtract their percentage from u
			if (it != t->tree[i].pcount.end()) {
				u -= (double) t->tree[i].pcount[groupB] / (double) tmap->seqsPerGroup[groupB];
				pcountSize++;
			}
						
			u = abs(u * t->tree[i].getBranchLength());
			
			nodePcountSize[i] = pcountSize;
				
			//if you are a leaf from a users group add to total
			if (i < numLeaves) {
				if ((t->tree[i].getBranchLength() != -1) && pcountSize != 0) { 
				//cout << "added to total" << endl; 
					WScore[(groupA+groupB)] += u; 
				}
				tempTotals[i] = 0.0;  //we don't care about you, or we have already added you
			}else{ //if you are not a leaf 
				//do both your chidren have have descendants from the users groups? 
				int lc = t->tree[i].getLChild();
				int rc = t->tree[i].getRChild();
				
				//if yes, add your childrens tempTotals
				if ((nodePcountSize[lc] != 0) && (nodePcountSize[rc] != 0)) {
					WScore[(groupA+groupB)] += tempTotals[lc] + tempTotals[rc]; 
					//cout << "added to total " << tempTotals[lc] << '\t' << tempTotals[rc] << endl;
					if (t->tree[i].getBranchLength() != -1) {
						tempTotals[i] = u;
					}else {
						tempTotals[i] = 0.0;
					}
				}else if ((nodePcountSize[lc] == 0) && (nodePcountSize[rc] == 0)) { tempTotals[i] = 0.0;  //we don't care about you
				}else { //if no, your tempTotal is your childrens temp totals + your branch length
					tempTotals[i] = tempTotals[lc] + tempTotals[rc] + u; 
				}
				//cout << "temptotal = "<< tempTotals[i] << endl;
			}
		}
		
		/********************************************************/
		
		//calculate weighted score for the group combination
		double UN;	
		UN = (WScore[(groupA+groupB)] / D);
		
		if (isnan(UN) || isinf(UN)) { UN = 0; } 
		data.push_back(UN);
				
		return data; 
	}
	catch(exception& e) {
		m->errorOut(e, "Weighted", "getValues");
		exit(1);
	}
}
/**************************************************************************************************/
double Weighted::getLengthToRoot(Tree* t, int v, string groupA, string groupB) { 
	try {
		
		double sum = 0.0;
		map<int, double> tempTotals; //maps node to total Branch Length
		map<int, int> nodePcountSize; //maps node to pcountSize
		map<int, int>::iterator itCount;

		int index = v;
	
		//you are a leaf
		if(t->tree[index].getBranchLength() != -1){	sum += abs(t->tree[index].getBranchLength());	}
		tempTotals[index] = 0.0;
		index = t->tree[index].getParent();	
			
		//while you aren't at root
		while(t->tree[index].getParent() != -1){

			if (m->control_pressed) {  return sum; }
				
			int pcountSize = 0;
			map<string, int>::iterator itGroup = t->tree[index].pcount.find(groupA);
			if (itGroup != t->tree[index].pcount.end()) { pcountSize++;  } 
			itGroup = t->tree[index].pcount.find(groupB);
			if (itGroup != t->tree[index].pcount.end()) { pcountSize++;  } 
		
			nodePcountSize[index] = pcountSize;
			
			//do both your chidren have have descendants from the users groups? 
			int lc = t->tree[index].getLChild();
			int rc = t->tree[index].getRChild();
			
			itCount = nodePcountSize.find(lc); 
			if (itCount == nodePcountSize.end()) { 
				int LpcountSize = 0;
				itGroup = t->tree[lc].pcount.find(groupA);
				if (itGroup != t->tree[lc].pcount.end()) { LpcountSize++;  } 
				itGroup = t->tree[lc].pcount.find(groupB);
				if (itGroup != t->tree[lc].pcount.end()) { LpcountSize++;  } 
				nodePcountSize[lc] = LpcountSize;
			}
			
			itCount = nodePcountSize.find(rc); 
			if (itCount == nodePcountSize.end()) { 
				int RpcountSize = 0;
				itGroup = t->tree[rc].pcount.find(groupA);
				if (itGroup != t->tree[rc].pcount.end()) { RpcountSize++;  } 
				itGroup = t->tree[rc].pcount.find(groupB);
				if (itGroup != t->tree[rc].pcount.end()) { RpcountSize++;  } 
				nodePcountSize[rc] = RpcountSize;
			}
				
			//if yes, add your childrens tempTotals
			if ((nodePcountSize[lc] != 0) && (nodePcountSize[rc] != 0)) {
				sum += tempTotals[lc] + tempTotals[rc]; 

				//cout << "added to total " << tempTotals[lc] << '\t' << tempTotals[rc] << endl;
				if (t->tree[index].getBranchLength() != -1) {
					tempTotals[index] = abs(t->tree[index].getBranchLength());
				}else {
					tempTotals[index] = 0.0;
				}
			}else { //if no, your tempTotal is your childrens temp totals + your branch length
				tempTotals[index] = tempTotals[lc] + tempTotals[rc] + abs(t->tree[index].getBranchLength()); 
			}
			//cout << "temptotal = "<< tempTotals[i] << endl;
			
			index = t->tree[index].getParent();	
		}

		return sum;
	}
	catch(exception& e) {
		m->errorOut(e, "Weighted", "getBranchLengthSums");
		exit(1);
	}
}
/**************************************************************************************************/



