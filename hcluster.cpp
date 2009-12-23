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

HCluster::HCluster(RAbundVector* rav, ListVector* lv, string m) :  rabund(rav), list(lv), method(m){
	try {
		mapWanted = false;
		exitedBreak = false; 
		numSeqs = list->getNumSeqs();
		
		//initialize cluster array
		for (int i = 0; i < numSeqs; i++) {
			clusterNode temp(1, -1, i);
			clusterArray.push_back(temp);
		}
		
	}
	catch(exception& e) {
		errorOut(e, "HCluster", "HCluster");
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
		errorOut(e, "HCluster", "clusterBins");
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
		errorOut(e, "HCluster", "clusterNames");
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
		errorOut(e, "HCluster", "getUpmostParent");
		exit(1);
	}
}
/***********************************************************************/
void HCluster::printInfo(){
	try {
		
		cout << "link table" << endl;
		for (it = activeLinks.begin(); it!= activeLinks.end(); it++) {
			cout << it->first << " = " << it->second << endl;
		}
		cout << endl;
		for (int i = 0; i < linkTable.size(); i++) {
			cout << i << '\t';
			for (it = linkTable[i].begin(); it != linkTable[i].end(); it++) {
				cout << it->first << '-' << it->second << '\t';
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
		errorOut(e, "HCluster", "getUpmostParent");
		exit(1);
	}
}
/***********************************************************************/
int HCluster::makeActive() {
	try {
	
		int linkValue = 1; 
//cout << "active - here" << endl;		
		it = activeLinks.find(smallRow);
		it2 = activeLinks.find(smallCol);
		
		if ((it == activeLinks.end()) && (it2 == activeLinks.end())) { //both are not active so add them
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
//cout << "active - here1" << endl;
		}else if ((it != activeLinks.end()) && (it2 == activeLinks.end())) {  //smallRow is active, smallCol is not
			 int size = linkTable.size();
			 int alreadyActiveRow = it->second;
			 map<int, int> temp; 
			
			//add link to eachother
			temp[smallRow] = 1;							//	   6	2	3	5
			linkTable.push_back(temp);					// 6   0	1	2	0
			linkTable[alreadyActiveRow][smallCol] = 1;	// 2   1	0	1	1
														// 3   2	1	0	0
														// 5   0    1   0   0	
			//add to activeLinks
			activeLinks[smallCol] = size;
//cout << "active - here2" << endl;			
		}else if ((it == activeLinks.end()) && (it2 != activeLinks.end())) {  //smallCol is active, smallRow is not
			 int size = linkTable.size();
			 int alreadyActiveCol = it2->second;
			 map<int, int> temp; 
			
			//add link to eachother
			temp[smallCol] = 1;							//	   6	2	3	5
			linkTable.push_back(temp);					// 6   0	1	2	0
			linkTable[alreadyActiveCol][smallRow] = 1;	// 2   1	0	1	1
														// 3   2	1	0	0
														// 5   0    1   0   0	
			//add to activeLinks
			activeLinks[smallRow] = size;
//cout << "active - here3" << endl;
		}else { //both are active so add one
			int row = it->second;
			int col = it2->second;
//cout << row << '\t' << col << endl;			
			
			linkTable[row][smallCol]++;
			linkTable[col][smallRow]++;
			linkValue = linkTable[row][smallCol];
//cout << "active - here4" << endl;
		}
		
		return linkValue;
	}
	catch(exception& e) {
		errorOut(e, "HCluster", "makeActive");
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
			
			//update linkTable by merging clustered rows and columns
			int rowSpot = activeLinks[smallRow];
			int colSpot = activeLinks[smallCol];
	//cout << "here" << endl;		
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
	//printInfo();
//cout << "here2" << endl;
			//merge their values
			for (it = linkTable[rowSpot].begin(); it != linkTable[rowSpot].end(); it++) {
				it2 = linkTable[colSpot].find(it->first);  //does the col also have this
				
				if (it2 == linkTable[colSpot].end()) { //not there so add it
					linkTable[colSpot][it->first] = it->second;
				}else { //merge them
					linkTable[colSpot][it->first] = it->second+it2->second;
				}
			}
//cout << "here3" << endl;			
			linkTable[colSpot].erase(size);
			linkTable.erase(linkTable.begin()+rowSpot);  //delete row
	//printInfo();		
			//update activerows
			activeLinks.erase(smallRow);
			activeLinks.erase(smallCol);
			activeLinks[size] = colSpot;
			
			//adjust everybody elses spot since you deleted - time vs. space
			for (it = activeLinks.begin(); it != activeLinks.end(); it++) {
				if (it->second > rowSpot) {  activeLinks[it->first]--;	}
			}
			
//cout << "here4" << endl;
	
	}
	catch(exception& e) {
		errorOut(e, "HCluster", "updateArrayandLinkTable");
		exit(1);
	}
}
/***********************************************************************/
void HCluster::update(int row, int col, float distance){
	try {
		
		smallRow = row;
		smallCol = col;
		smallDist = distance;
		
		//find upmost parent of row and col
		smallRow = getUpmostParent(smallRow);
		smallCol = getUpmostParent(smallCol);
	//cout << "row = " << row << " smallRow = " << smallRow <<  " col = " << col << " smallCol = " << smallCol << " dist = " << distance << endl;
		//you don't want to cluster with yourself
		if (smallRow != smallCol) {
			//are they active in the link table
			int linkValue = makeActive(); //after this point this nodes info is active in linkTable
			//printInfo();			
			//cout << "linkValue = " << linkValue << " times = " << (clusterArray[smallRow].numSeq * clusterArray[smallCol].numSeq) << endl;
			//can we cluster???
			bool cluster = false;
			
			if (method == "nearest") { cluster = true;  }
			else if (method == "average") { cout << "still working on this... " << endl; //got to figure this out 
			}else{ //assume furthest
				if (linkValue == (clusterArray[smallRow].numSeq * clusterArray[smallCol].numSeq)) { cluster = true; }
			}
			
			if (cluster) { 
				updateArrayandLinkTable();
				clusterBins();
				clusterNames();
			}
		}
		//printInfo();
	}
	catch(exception& e) {
		errorOut(e, "HCluster", "update");
		exit(1);
	}
}
/***********************************************************************/
void HCluster::setMapWanted(bool m)  {  
	try {
		mapWanted = m;
		
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
		errorOut(e, "HCluster", "setMapWanted");
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
		errorOut(e, "HCluster", "updateMap");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<seqDist> HCluster::getSeqs(ifstream& filehandle, NameAssignment* nameMap, float cutoff){
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
			
			filehandle >> firstName >> secondName >> distance;    gobble(filehandle); 
	
			//save first one
			if (prevDistance == -1) { prevDistance = distance; }
			
			map<string,int>::iterator itA = nameMap->find(firstName);
			map<string,int>::iterator itB = nameMap->find(secondName);
			if(itA == nameMap->end()){  cerr << "AAError: Sequence '" << firstName << "' was not found in the names file, please correct\n"; exit(1);  }
			if(itB == nameMap->end()){  cerr << "ABError: Sequence '" << secondName << "' was not found in the names file, please correct\n"; exit(1);  }
		
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
		errorOut(e, "HCluster", "getSeqs");
		exit(1);
	}
}
/***********************************************************************/



