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
	
		vector<double> sums = getBranchLengthSums(t);
		
		if (m->control_pressed) { return data; }
		
		//calculate number of comparisons i.e. with groups A,B,C = AB, AC, BC = 3;
		vector< vector<string> > namesOfGroupCombos;
		for (int i=0; i<numGroups; i++) { 
			for (int l = i+1; l < numGroups; l++) {	
				//initialize weighted scores
				//WScore[globaldata->Groups[i]+globaldata->Groups[l]] = 0.0;
				vector<string> groups; groups.push_back(globaldata->Groups[i]); groups.push_back(globaldata->Groups[l]);
				namesOfGroupCombos.push_back(groups);
			}
		}
		
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			if(processors == 1){
				data = driver(t, namesOfGroupCombos, 0, namesOfGroupCombos.size(), sums);
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

				data = createProcesses(t, namesOfGroupCombos, sums);
				
				for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
			}
		#else
			data = driver(t, namesOfGroupCombos, 0, namesOfGroupCombos.size(), sums);
		#endif
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Weighted", "getValues");
		exit(1);
	}
}
/**************************************************************************************************/

EstOutput Weighted::createProcesses(Tree* t, vector< vector<string> > namesOfGroupCombos, vector<double>& sums) {
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
				Myresults = driver(t, namesOfGroupCombos, lines[process]->start, lines[process]->num, sums);
				
				if (m->control_pressed) { exit(0); }
				
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
	
		results = driver(t, namesOfGroupCombos, lines[0]->start, lines[0]->num, sums);
		
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
EstOutput Weighted::driver(Tree* t, vector< vector<string> > namesOfGroupCombos, int start, int num, vector<double>& sums) { 
 try {
		globaldata = GlobalData::getInstance();
		
		EstOutput results;
		vector<double> D;
		
		int count = 0;
		int total = num;
		int twentyPercent = (total * 0.20);

		for (int h = start; h < (start+num); h++) {
		
			if (m->control_pressed) { return results; }
		
			//initialize weighted score
			string groupA = namesOfGroupCombos[h][0]; 
			string groupB = namesOfGroupCombos[h][1];
			
			WScore[groupA+groupB] = 0.0;
			D.push_back(0.0000); //initialize a spot in D for each combination
			
			//adding the wieghted sums from group i
			for (int j = 0; j < t->groupNodeInfo[groupA].size(); j++) { //the leaf nodes that have seqs from group i
				map<string, int>::iterator it = t->tree[t->groupNodeInfo[groupA][j]].pcount.find(groupA);
				int numSeqsInGroupI = it->second;
				
				double weightedSum = ((numSeqsInGroupI * sums[t->groupNodeInfo[groupA][j]]) / (double)tmap->seqsPerGroup[groupA]);
			
				D[count] += weightedSum;
			}
			
			//adding the wieghted sums from group l
			for (int j = 0; j < t->groupNodeInfo[groupB].size(); j++) { //the leaf nodes that have seqs from group l
				map<string, int>::iterator it = t->tree[t->groupNodeInfo[groupB][j]].pcount.find(groupB);
				int numSeqsInGroupL = it->second;
				
				double weightedSum = ((numSeqsInGroupL * sums[t->groupNodeInfo[groupB][j]]) / (double)tmap->seqsPerGroup[groupB]);
			
				D[count] += weightedSum;
			}
			count++;
			
			//report progress
			if((count) % twentyPercent == 0){	m->mothurOut("Percentage complete: " + toString(int((count / (float)total) * 100.0))); m->mothurOutEndLine();		}
		}
		
		m->mothurOut("Percentage complete: 100"); m->mothurOutEndLine();
		
		//calculate u for the group comb 
		for(int i=0;i<t->getNumNodes();i++){
			for (int h = start; h < (start+num); h++) {
				
				string groupA = namesOfGroupCombos[h][0]; 
				string groupB = namesOfGroupCombos[h][1];

				if (m->control_pressed) { return results; }
				
				double u;
				//does this node have descendants from groupA
				it = t->tree[i].pcount.find(groupA);
				//if it does u = # of its descendants with a certain group / total number in tree with a certain group
				if (it != t->tree[i].pcount.end()) {
					u = (double) t->tree[i].pcount[groupA] / (double) tmap->seqsPerGroup[groupA];
				}else { u = 0.00; }
			
							
				//does this node have descendants from group l
				it = t->tree[i].pcount.find(groupB);
				//if it does subtract their percentage from u
				if (it != t->tree[i].pcount.end()) {
					u -= (double) t->tree[i].pcount[groupB] / (double) tmap->seqsPerGroup[groupB];
				}
							
				u = abs(u * t->tree[i].getBranchLength());
						
				//save groupcombs u value
				WScore[(groupA+groupB)] += u;
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
EstOutput Weighted::getValues(Tree* t, string groupA, string groupB, vector<double>& sums) { 
 try {
		globaldata = GlobalData::getInstance();
		
		data.clear(); //clear out old values
		
		if (m->control_pressed) { return data; }
		
		//initialize weighted score
		WScore[(groupA+groupB)] = 0.0;
		double D = 0.0;
		
		vector<string> groups; groups.push_back(groupA); groups.push_back(groupB);
		
		//adding the wieghted sums from group i
		for (int j = 0; j < t->groupNodeInfo[groups[0]].size(); j++) { //the leaf nodes that have seqs from group i
			map<string, int>::iterator it = t->tree[t->groupNodeInfo[groups[0]][j]].pcount.find(groups[0]);
			int numSeqsInGroupI = it->second;
			
			double weightedSum = ((numSeqsInGroupI * sums[t->groupNodeInfo[groups[0]][j]]) / (double)tmap->seqsPerGroup[groups[0]]);
		
			D += weightedSum;
		}
		
		//adding the wieghted sums from group l
		for (int j = 0; j < t->groupNodeInfo[groups[1]].size(); j++) { //the leaf nodes that have seqs from group l
			map<string, int>::iterator it = t->tree[t->groupNodeInfo[groups[1]][j]].pcount.find(groups[1]);
			int numSeqsInGroupL = it->second;
			
			double weightedSum = ((numSeqsInGroupL * sums[t->groupNodeInfo[groups[1]][j]]) / (double)tmap->seqsPerGroup[groups[1]]);
		
			D += weightedSum;
		}
		
		//calculate u for the group comb 
		for(int i=0;i<t->getNumNodes();i++){
		
			if (m->control_pressed) { return data; }
			
			double u;
			//does this node have descendants from groupA
			it = t->tree[i].pcount.find(groupA);
			//if it does u = # of its descendants with a certain group / total number in tree with a certain group
			if (it != t->tree[i].pcount.end()) {
				u = (double) t->tree[i].pcount[groupA] / (double) tmap->seqsPerGroup[groupA];
			}else { u = 0.00; }
		
						
			//does this node have descendants from group l
			it = t->tree[i].pcount.find(groupB);
			//if it does subtract their percentage from u
			if (it != t->tree[i].pcount.end()) {
				u -= (double) t->tree[i].pcount[groupB] / (double) tmap->seqsPerGroup[groupB];
			}
						
			u = abs(u * t->tree[i].getBranchLength());
					
			//save groupcombs u value
			WScore[(groupA+groupB)] += u;
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
vector<double> Weighted::getBranchLengthSums(Tree* t) { 
	try {
		
		vector<double> sums;
		
		for(int v=0;v<t->getNumLeaves();v++){
			
			if (m->control_pressed) { return sums; }
				
			int index = v;
			double sum = 0.0000;
	
			//while you aren't at root
			while(t->tree[index].getParent() != -1){
						
				//if you have a BL
				if(t->tree[index].getBranchLength() != -1){
					sum += abs(t->tree[index].getBranchLength());
				}
				index = t->tree[index].getParent();
			}
					
			//get last breanch length added
			if(t->tree[index].getBranchLength() != -1){
				sum += abs(t->tree[index].getBranchLength());
			}
			
			sums.push_back(sum);
		}
		
		return sums;
	}
	catch(exception& e) {
		m->errorOut(e, "Weighted", "getBranchLengthSums");
		exit(1);
	}
}
/**************************************************************************************************/



