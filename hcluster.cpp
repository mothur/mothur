/*
 *  hcluster.cpp
 *  Mothur
 *
 *  Created by westcott on 10/13/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "hcluster.h"
#include "rabundvector.hpp"
#include "listvector.hpp"
#include "sparsematrix.hpp"

/***********************************************************************/
HCluster::HCluster(RAbundVector* rav, ListVector* lv, string ms, string d, NameAssignment* n, float c) :  rabund(rav), list(lv), method(ms), distfile(d), nameMap(n), cutoff(c) {
	try {
		m = MothurOut::getInstance();
		mapWanted = false;
		exitedBreak = false; 
		numSeqs = list->getNumSeqs();
		
		//initialize cluster array
		for (int i = 0; i < numSeqs; i++) {
			clusterNode temp(1, -1, i);
			clusterArray.push_back(temp);
		}
		
		if ((method == "furthest") || (method == "nearest")) {
			m->openInputFile(distfile, filehandle);
		}else{  
			processFile();  
		}
	}
	catch(exception& e) {
		m->errorOut(e, "HCluster", "HCluster");
		exit(1);
	}
}
/***********************************************************************/

void HCluster::clusterBins(){
	try {
		//cout << smallCol << '\t' << smallRow << '\t' << smallDist << '\t' << rabund->get(clusterArray[smallRow].smallChild) << '\t' << rabund->get(clusterArray[smallCol].smallChild);

		rabund->set(clusterArray[smallCol].smallChild, rabund->get(clusterArray[smallRow].smallChild)+rabund->get(clusterArray[smallCol].smallChild));	
		rabund->set(clusterArray[smallRow].smallChild, 0);	
		rabund->setLabel(toString(smallDist));

		//cout << '\t' << rabund->get(clusterArray[smallRow].smallChild) << '\t' << rabund->get(clusterArray[smallCol].smallChild) << endl;
	}
	catch(exception& e) {
		m->errorOut(e, "HCluster", "clusterBins");
		exit(1);
	}


}

/***********************************************************************/

void HCluster::clusterNames(){
	try {
		///cout << smallCol << '\t' << smallRow << '\t' << smallDist << '\t' << list->get(clusterArray[smallRow].smallChild) << '\t' << list->get(clusterArray[smallCol].smallChild);
		if (mapWanted) {  updateMap();  }
		
		list->set(clusterArray[smallCol].smallChild, list->get(clusterArray[smallRow].smallChild)+','+list->get(clusterArray[smallCol].smallChild));
		list->set(clusterArray[smallRow].smallChild, "");	
		list->setLabel(toString(smallDist));
	
		//cout << '\t' << list->get(clusterArray[smallRow].smallChild) << '\t' << list->get(clusterArray[smallCol].smallChild) << endl;

    }
	catch(exception& e) {
		m->errorOut(e, "HCluster", "clusterNames");
		exit(1);
	}

}
/***********************************************************************/
int HCluster::getUpmostParent(int node){
	try {
		while (clusterArray[node].parent != -1) {
			node = clusterArray[node].parent;
		}
		
		return node;
	}
	catch(exception& e) {
		m->errorOut(e, "HCluster", "getUpmostParent");
		exit(1);
	}
}
/***********************************************************************/
void HCluster::printInfo(){
	try {
		
		cout << "link table" << endl;
		for (itActive = activeLinks.begin(); itActive!= activeLinks.end(); itActive++) {
			cout << itActive->first << " = " << itActive->second << endl;
		}
		cout << endl;
		for (int i = 0; i < linkTable.size(); i++) {
			cout << i << '\t';
			for (it = linkTable[i].begin(); it != linkTable[i].end(); it++) {
				cout << it->first << '-' << it->second << '\t' ;
			}
			cout << endl;
		}
		cout << endl << "clusterArray" << endl;
		
		for (int i = 0; i < clusterArray.size(); i++) {
			cout << i << '\t' << clusterArray[i].numSeq << '\t' << clusterArray[i].parent << '\t' << clusterArray[i].smallChild << endl;
		}
		cout << endl;
		
		
	}
	catch(exception& e) {
		m->errorOut(e, "HCluster", "getUpmostParent");
		exit(1);
	}
}
/***********************************************************************/
int HCluster::makeActive() {
	try {
		int linkValue = 1; 

		itActive = activeLinks.find(smallRow);
		it2Active = activeLinks.find(smallCol);
		
		if ((itActive == activeLinks.end()) && (it2Active == activeLinks.end())) { //both are not active so add them
			int size = linkTable.size();
			map<int, int> temp; map<int, int> temp2;
			
			//add link to eachother
			temp[smallRow] = 1;							//	   1	2  								
			temp2[smallCol] = 1;						// 1   0	1 
														// 2   1	0				
			linkTable.push_back(temp);
			linkTable.push_back(temp2);
			
			//add to activeLinks
			activeLinks[smallRow] = size;
			activeLinks[smallCol] = size+1;

		}else if ((itActive != activeLinks.end()) && (it2Active == activeLinks.end())) {  //smallRow is active, smallCol is not
			 int size = linkTable.size();
			 int alreadyActiveRow = itActive->second;
			 map<int, int> temp; 
			
			//add link to eachother
			temp[smallRow] = 1;									//	   6	2	3	5
			linkTable.push_back(temp);							// 6   0	1	2	0
			linkTable[alreadyActiveRow][smallCol] = 1;			// 2   1	0	1	1
																// 3   2	1	0	0
																// 5   0    1   0   0	
			//add to activeLinks
			activeLinks[smallCol] = size;
	
		}else if ((itActive == activeLinks.end()) && (it2Active != activeLinks.end())) {  //smallCol is active, smallRow is not
			 int size = linkTable.size();
			 int alreadyActiveCol = it2Active->second;
			 map<int, int> temp; 
			
			//add link to eachother
			temp[smallCol] = 1;									//	   6	2	3	5
			linkTable.push_back(temp);							// 6   0	1	2	0
			linkTable[alreadyActiveCol][smallRow] = 1;			// 2   1	0	1	1
																// 3   2	1	0	0
																// 5   0    1   0   0	
			//add to activeLinks
			activeLinks[smallRow] = size;

		}else { //both are active so add one
			int row = itActive->second;
			int col = it2Active->second;
			
			
			linkTable[row][smallCol]++;
			linkTable[col][smallRow]++;
			linkValue = linkTable[row][smallCol];
		}
		
		return linkValue;
	}
	catch(exception& e) {
		m->errorOut(e, "HCluster", "makeActive");
		exit(1);
	}
}
/***********************************************************************/
void HCluster::updateArrayandLinkTable() {
	try {
		//if cluster was made update clusterArray and linkTable
		int size = clusterArray.size();
		
		//add new node
		clusterNode temp(clusterArray[smallRow].numSeq + clusterArray[smallCol].numSeq, -1, clusterArray[smallCol].smallChild);
		clusterArray.push_back(temp);
		
		//update child nodes
		clusterArray[smallRow].parent = size;
		clusterArray[smallCol].parent = size;
		
		if (method == "furthest") {
			
			//update linkTable by merging clustered rows and columns
			int rowSpot = activeLinks[smallRow];
			int colSpot = activeLinks[smallCol];
			
			//fix old rows
			for (int i = 0; i < linkTable.size(); i++) {
				//check if they are in map
				it = linkTable[i].find(smallRow);
				it2 = linkTable[i].find(smallCol);
				
				if ((it!=linkTable[i].end()) && (it2!=linkTable[i].end())) { //they are both there
					linkTable[i][size] = linkTable[i][smallRow]+linkTable[i][smallCol];
					linkTable[i].erase(smallCol); //delete col row
					linkTable[i].erase(smallRow); //delete col row
				}else if ((it==linkTable[i].end()) && (it2!=linkTable[i].end())) { //only col
					linkTable[i][size] = linkTable[i][smallCol];
					linkTable[i].erase(smallCol); //delete col 
				}else if ((it!=linkTable[i].end()) && (it2==linkTable[i].end())) { //only row
					linkTable[i][size] = linkTable[i][smallRow];
					linkTable[i].erase(smallRow); //delete col 
				}
			}
			
			//merge their values
			for (it = linkTable[rowSpot].begin(); it != linkTable[rowSpot].end(); it++) {
				it2 = linkTable[colSpot].find(it->first);  //does the col also have this
				
				if (it2 == linkTable[colSpot].end()) { //not there so add it
					linkTable[colSpot][it->first] = it->second;
				}else { //merge them
					linkTable[colSpot][it->first] = it->second + it2->second;
				}
			}
			
			linkTable[colSpot].erase(size);
			linkTable.erase(linkTable.begin()+rowSpot);  //delete row
			
			//update activerows
			activeLinks.erase(smallRow);
			activeLinks.erase(smallCol);
			activeLinks[size] = colSpot;
			
			//adjust everybody elses spot since you deleted - time vs. space
			for (itActive = activeLinks.begin(); itActive != activeLinks.end(); itActive++) {
				if (itActive->second > rowSpot) {  activeLinks[itActive->first]--;	}
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "HCluster", "updateArrayandLinkTable");
		exit(1);
	}
}
/***********************************************************************/
double HCluster::update(int row, int col, float distance){
	try {
		bool cluster = false;
		smallRow = row;
		smallCol = col;
		smallDist = distance;
		
		//find upmost parent of row and col
		smallRow = getUpmostParent(smallRow);
		smallCol = getUpmostParent(smallCol);

		//you don't want to cluster with yourself
		if (smallRow != smallCol) {
			
			if ((method == "furthest") || (method == "nearest")) {
				//can we cluster???
				if (method == "nearest") { cluster = true;  }
				else{ //assume furthest
					//are they active in the link table
					int linkValue = makeActive(); //after this point this nodes info is active in linkTable
					if (linkValue == (clusterArray[smallRow].numSeq * clusterArray[smallCol].numSeq)) {		cluster = true;		}
				}
				
				if (cluster) { 
					updateArrayandLinkTable();
					clusterBins();
					clusterNames();
				}
			}else {
				cluster = true;
				updateArrayandLinkTable();
				clusterBins();
				clusterNames();
				combineFile();
			}
		}
		
		return cutoff;
		//printInfo();
	}
	catch(exception& e) {
		m->errorOut(e, "HCluster", "update");
		exit(1);
	}
}
/***********************************************************************/
void HCluster::setMapWanted(bool ms)  {  
	try {
		mapWanted = ms;
		
		//initialize map
		for (int i = 0; i < list->getNumBins(); i++) {
			
			//parse bin 
			string names = list->get(i);
			while (names.find_first_of(',') != -1) { 
				//get name from bin
				string name = names.substr(0,names.find_first_of(','));
				//save name and bin number
				seq2Bin[name] = i;
				names = names.substr(names.find_first_of(',')+1, names.length());
			}
			
			//get last name
			seq2Bin[names] = i;
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "HCluster", "setMapWanted");
		exit(1);
	}
}
/***********************************************************************/
void HCluster::updateMap() {
try {
		//update location of seqs in smallRow since they move to smallCol now
		string names = list->get(clusterArray[smallRow].smallChild);
		while (names.find_first_of(',') != -1) { 
			//get name from bin
			string name = names.substr(0,names.find_first_of(','));
			//save name and bin number
			seq2Bin[name] = clusterArray[smallCol].smallChild;
			names = names.substr(names.find_first_of(',')+1, names.length());
		}
			
		//get last name
		seq2Bin[names] = clusterArray[smallCol].smallChild;
	}
	catch(exception& e) {
		m->errorOut(e, "HCluster", "updateMap");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<seqDist> HCluster::getSeqs(){
	try {
		vector<seqDist> sameSeqs;
		
		if ((method == "furthest") || (method == "nearest")) {
			sameSeqs = getSeqsFNNN();
		}else{
			sameSeqs = getSeqsAN();	
		}
				
		return sameSeqs;
	}
	catch(exception& e) {
		m->errorOut(e, "HCluster", "getSeqs");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<seqDist> HCluster::getSeqsFNNN(){
	try {
		string firstName, secondName;
		float distance, prevDistance;
		vector<seqDist> sameSeqs;
		prevDistance = -1;
	
		//if you are not at the beginning of the file
		if (exitedBreak) { 
			sameSeqs.push_back(next);
			prevDistance = next.dist;
			exitedBreak = false;
		}
	
		//get entry
		while (!filehandle.eof()) {
			
			filehandle >> firstName >> secondName >> distance;    m->gobble(filehandle); 
	
			//save first one
			if (prevDistance == -1) { prevDistance = distance; }
			
			map<string,int>::iterator itA = nameMap->find(firstName);
			map<string,int>::iterator itB = nameMap->find(secondName);
			if(itA == nameMap->end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the names file, please correct\n"); exit(1);  }
			if(itB == nameMap->end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the names file, please correct\n"); exit(1);  }
		
			//using cutoff
			if (distance > cutoff) { break; }
		
			if (distance != -1) { //-1 means skip me
			
				//are the distances the same
				if (distance == prevDistance) { //save in vector
					seqDist temp(itA->second, itB->second, distance);
					sameSeqs.push_back(temp);
					exitedBreak = false;
				}else{ 
					next.seq1 = itA->second;
					next.seq2 = itB->second;
					next.dist = distance;
					exitedBreak = true;
					break;
				}
			}
		}
		
		//rndomize matching dists
		random_shuffle(sameSeqs.begin(), sameSeqs.end());
		
		return sameSeqs;
	}
	catch(exception& e) {
		m->errorOut(e, "HCluster", "getSeqsFNNN");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<seqDist> HCluster::getSeqsAN(){
	try {
		int firstName, secondName;
		float prevDistance;
		vector<seqDist> sameSeqs;
		prevDistance = -1;
		
		m->openInputFile(distfile, filehandle, "no error"); 
		
		//is the smallest value in mergedMin or the distfile?
		float mergedMinDist = 10000;
		float distance = 10000;
		if (mergedMin.size() > 0) { mergedMinDist = mergedMin[0].dist;  }
			
		if (!filehandle.eof()) {  
			filehandle >> firstName >> secondName >> distance;    m->gobble(filehandle);
			//save first one
			if (prevDistance == -1) { prevDistance = distance; } 
			if (distance != -1) { //-1 means skip me
				seqDist temp(firstName, secondName, distance);
				sameSeqs.push_back(temp);
			}else{ distance = 10000; }
		}
		
		if (mergedMinDist < distance) { //get minimum distance from mergedMin
			//remove distance we saved from file
			sameSeqs.clear();
			prevDistance = mergedMinDist;
			
			for (int i = 0; i < mergedMin.size(); i++) {
				if (mergedMin[i].dist == prevDistance) {
					sameSeqs.push_back(mergedMin[i]);
				}else { break; }
			}
		}else{ //get minimum from file
			//get entry
			while (!filehandle.eof()) {
				
				filehandle >> firstName >> secondName >> distance;    m->gobble(filehandle); 
				
				if (prevDistance == -1) { prevDistance = distance; }
				
				if (distance != -1) { //-1 means skip me
					//are the distances the same
					if (distance == prevDistance) { //save in vector
						seqDist temp(firstName, secondName, distance);
						sameSeqs.push_back(temp);
					}else{	
						break;
					}
				}
			}
		}
		filehandle.close();
		
		//randomize matching dists
		random_shuffle(sameSeqs.begin(), sameSeqs.end());
		
		//can only return one value since once these are merged the other distances in sameSeqs may have changed
		vector<seqDist> temp;
		if (sameSeqs.size() > 0) {  temp.push_back(sameSeqs[0]);  }
		
		return temp;
	}
	catch(exception& e) {
		m->errorOut(e, "HCluster", "getSeqsAN");
		exit(1);
	}
}

/***********************************************************************/
int HCluster::combineFile() {
	try {
		//int bufferSize = 64000;  //512k - this should be a variable that the user can set to optimize code to their hardware
		//char* inputBuffer;
		//inputBuffer = new char[bufferSize];
		//size_t numRead;
		
		string tempDistFile = distfile + ".temp";
		ofstream out;
		m->openOutputFile(tempDistFile, out);
		
		//FILE* in;
		//in = fopen(distfile.c_str(), "rb");
	
		ifstream in;
		m->openInputFile(distfile, in, "no error");
		
		int first, second;
		float dist;
		
		vector< map<int, float> > smallRowColValues;
		smallRowColValues.resize(2);  //0 = row, 1 = col
		int count = 0;
				
		//go through file pulling out distances related to rows merging
		//if mergedMin contains distances add those back into file
		//bool done = false;
		//partialDist = "";
		//while ((numRead = fread(inputBuffer, 1, bufferSize, in)) != 0) {
//cout << "number of char read = " << numRead << endl;
//cout << inputBuffer << endl;
			//if (numRead < bufferSize) { done = true; }
			
			//parse input into individual distances
			//int spot = 0;
			//string outputString = "";
			//while(spot < numRead) {
	//cout << "spot = " << spot << endl;
			 //  seqDist nextDist = getNextDist(inputBuffer, spot, bufferSize);
			   
			   //you read a partial distance
			  // if (nextDist.seq1 == -1) { break;  }
			while (!in.eof()) {
				//first = nextDist.seq1; second = nextDist.seq2; dist = nextDist.dist;
	//cout << "next distance = " << first << '\t' << second << '\t' << dist << endl;		   
			   //since file is sorted and mergedMin is sorted 
			   //you can put the smallest distance from each through the code below and keep the file sorted
			   
			   in >> first >> second >> dist; m->gobble(in);
			   
			   if (m->control_pressed) { in.close(); out.close(); m->mothurRemove(tempDistFile); return 0; }
			   
			   //while there are still values in mergedMin that are smaller than the distance read from file
			   while (count < mergedMin.size())  {
			   
					//is the distance in mergedMin smaller than from the file
					if (mergedMin[count].dist < dist) {
					//is this a distance related to the columns merging?
					//if yes, save in memory
						if ((mergedMin[count].seq1 == smallRow) && (mergedMin[count].seq2 == smallCol)) { //do nothing this is the smallest distance from last time
						}else if (mergedMin[count].seq1 == smallCol) {
							smallRowColValues[1][mergedMin[count].seq2] = mergedMin[count].dist;
						}else if (mergedMin[count].seq2 == smallCol) {
							smallRowColValues[1][mergedMin[count].seq1] = mergedMin[count].dist;
						}else if (mergedMin[count].seq1 == smallRow) {
							smallRowColValues[0][mergedMin[count].seq2] = mergedMin[count].dist;
						}else if (mergedMin[count].seq2 == smallRow) {
							smallRowColValues[0][mergedMin[count].seq1] = mergedMin[count].dist;
						}else { //if no, write to temp file
							//outputString += toString(mergedMin[count].seq1) + '\t' + toString(mergedMin[count].seq2) + '\t' + toString(mergedMin[count].dist) + '\n';
							//if (mergedMin[count].dist < cutoff) { 
								out << mergedMin[count].seq1 << '\t' << mergedMin[count].seq2 << '\t' << mergedMin[count].dist << endl;
							//}
						}
						count++;
					}else{   break;	}
			   }
			   
			   //is this a distance related to the columns merging?
			   //if yes, save in memory
			   if ((first == smallRow) && (second == smallCol)) { //do nothing this is the smallest distance from last time
			   }else if (first == smallCol) {
					smallRowColValues[1][second] = dist;
			   }else if (second == smallCol) {
					smallRowColValues[1][first] = dist;
			   }else if (first == smallRow) {
					smallRowColValues[0][second] = dist;
			   }else if (second == smallRow) {
					smallRowColValues[0][first] = dist;
			   
			   }else { //if no, write to temp file
					//outputString += toString(first) + '\t' + toString(second) + '\t' + toString(dist) + '\n';
				   //if (dist < cutoff) {
					   out << first << '\t' << second << '\t' << dist << endl;
				   //}
			   }
			}
			
			//out << outputString;
			//if(done) { break; }
		//}
		//fclose(in);
		in.close();
		
		//if values in mergedMin are larger than the the largest in file then
		while (count < mergedMin.size())  {  
			//is this a distance related to the columns merging?
			//if yes, save in memory
			if ((mergedMin[count].seq1 == smallRow) && (mergedMin[count].seq2 == smallCol)) { //do nothing this is the smallest distance from last time
			}else if (mergedMin[count].seq1 == smallCol) {
				smallRowColValues[1][mergedMin[count].seq2] = mergedMin[count].dist;
			}else if (mergedMin[count].seq2 == smallCol) {
				smallRowColValues[1][mergedMin[count].seq1] = mergedMin[count].dist;
			}else if (mergedMin[count].seq1 == smallRow) {
				smallRowColValues[0][mergedMin[count].seq2] = mergedMin[count].dist;
			}else if (mergedMin[count].seq2 == smallRow) {
				smallRowColValues[0][mergedMin[count].seq1] = mergedMin[count].dist;
				
			}else { //if no, write to temp file
				//if (mergedMin[count].dist < cutoff) {
					out << mergedMin[count].seq1 << '\t' << mergedMin[count].seq2 << '\t' << mergedMin[count].dist << endl;
				//}
			}
			count++;
		}
		out.close();
		mergedMin.clear();
			
		//rename tempfile to distfile
		m->mothurRemove(distfile);
		rename(tempDistFile.c_str(), distfile.c_str());
//cout << "remove = "<< renameOK << " rename = " << ok << endl;	

		//merge clustered rows averaging the distances
		map<int, float>::iterator itMerge;
		map<int, float>::iterator it2Merge;
		for(itMerge = smallRowColValues[0].begin(); itMerge != smallRowColValues[0].end(); itMerge++) {			
			//does smallRowColValues[1] have a distance to this seq too?
			it2Merge = smallRowColValues[1].find(itMerge->first);
			
			float average;
			if (it2Merge != smallRowColValues[1].end()) { //if yes, then average
				//average
				if (method == "average") {
					int total = clusterArray[smallRow].numSeq + clusterArray[smallCol].numSeq;
					average = ((clusterArray[smallRow].numSeq * itMerge->second) + (clusterArray[smallCol].numSeq * it2Merge->second)) / (float) total;
				}else { //weighted
					average = ((itMerge->second * 1.0) + (it2Merge->second * 1.0)) / (float) 2.0;				
				}
				
				smallRowColValues[1].erase(it2Merge);
				
				seqDist temp(clusterArray[smallRow].parent, itMerge->first, average);
				mergedMin.push_back(temp);
			}else {  
				//can't find value so update cutoff
				if (cutoff > itMerge->second) { cutoff = itMerge->second; }
			}
		}
		
		//update cutoff
		for(itMerge = smallRowColValues[1].begin(); itMerge != smallRowColValues[1].end(); itMerge++) {	
			if (cutoff > itMerge->second) { cutoff = itMerge->second; }
		}
		
		//sort merged values
		sort(mergedMin.begin(), mergedMin.end(), compareSequenceDistance);	
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "HCluster", "combineFile");
		exit(1);
	}
}
/***********************************************************************
seqDist HCluster::getNextDist(char* buffer, int& index, int size){
	try {
		seqDist next;
		int indexBefore = index;
		string first, second, distance;
		first = ""; second = ""; distance = "";
		int tabCount = 0;
//cout << "partial = " << partialDist << endl;		
		if (partialDist != "") { //read what you can, you know it is less than a whole distance.
			for (int i = 0; i < partialDist.size(); i++) {
				if (tabCount == 0) {
					if (partialDist[i] == '\t') { tabCount++; }
					else {  first += partialDist[i];	}
				}else if (tabCount == 1) {
					if (partialDist[i] == '\t') { tabCount++; }
					else {  second += partialDist[i];	}
				}else if (tabCount == 2) {
					distance +=  partialDist[i];
				}
			}
			partialDist = "";
		}
	
		//try to get another distance
		bool gotDist = false;
		while (index < size) {
			if ((buffer[index] == 10) || (buffer[index] == 13)) { //newline in unix or windows
				gotDist = true;
				
				//m->gobble space
				while (index < size) {		
					if (isspace(buffer[index])) { index++; }
					else { break; }		
				}
				break;
			}else{
				if (tabCount == 0) {
					if (buffer[index] == '\t') { tabCount++; }
					else {  first += buffer[index];	}
				}else if (tabCount == 1) {
					if (buffer[index] == '\t') { tabCount++; }
					else {  second += buffer[index];	}
				}else if (tabCount == 2) {
					distance +=  buffer[index];
				}
				index++;
			}
		}
		
		//there was not a whole distance in the buffer, ie. buffer = "1	2	0.01	2	3	0."
		//then you want to save the partial distance.
		if (!gotDist) {
			for (int i = indexBefore; i < size; i++) {
				partialDist += buffer[i];
			}
			index = size + 1;
			next.seq1 = -1; next.seq2 = -1; next.dist = 0.0;
		}else{
			int firstname, secondname;
			float dist;
			
			convert(first, firstname);
			convert(second, secondname);
			convert(distance, dist);
			
			next.seq1 = firstname; next.seq2 = secondname; next.dist = dist;
		}
						
		return next;
	}
	catch(exception& e) {
		m->errorOut(e, "HCluster", "getNextDist");
		exit(1);
	}
}
/***********************************************************************/
int HCluster::processFile() {
	try {
		string firstName, secondName;
		float distance;
		
		ifstream in;
		m->openInputFile(distfile, in, "no error");
		
		ofstream out;
		string outTemp = distfile + ".temp";
		m->openOutputFile(outTemp, out);
	
		//get entry
		while (!in.eof()) {
			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove(outTemp); return 0; }
			
			in >> firstName >> secondName >> distance;    m->gobble(in);		
			
			map<string,int>::iterator itA = nameMap->find(firstName);
			map<string,int>::iterator itB = nameMap->find(secondName);
			if(itA == nameMap->end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the names file, please correct\n"); exit(1);  }
			if(itB == nameMap->end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the names file, please correct\n"); exit(1);  }
		
			//using cutoff
			if (distance > cutoff) { break; }
		
			if (distance != -1) { //-1 means skip me
				out << itA->second << '\t' << itB->second << '\t' << distance << endl;
			}
		}
		
		in.close();
		out.close();
		
		m->mothurRemove(distfile);
		rename(outTemp.c_str(), distfile.c_str());
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "HCluster", "processFile");
		exit(1);
	}
}
/***********************************************************************/








