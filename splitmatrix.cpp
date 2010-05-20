/*
 *  splitmatrix.cpp
 *  Mothur
 *
 *  Created by westcott on 5/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "splitmatrix.h"

/***********************************************************************/

SplitMatrix::SplitMatrix(string distfile, string name, float c){
	m = MothurOut::getInstance();
	distFile = distfile;
	cutoff = c;
	namefile = name;
}

/***********************************************************************/

int SplitMatrix::split(){
	try {
        
		vector<set<string> > groups;
		int numGroups = 0;

		ofstream outFile;
		ifstream dFile;
		openInputFile(distFile, dFile);
	
		while(dFile){
			string seqA, seqB;
			float dist;

			dFile >> seqA >> seqB >> dist;
			
			if(dist < cutoff){
				//cout << "in cutoff: " << dist << endl;
				int groupIDA = -1;
				int groupIDB = -1;
				int groupID = -1;
				int prevGroupID = -1;
				
				for(int i=0;i<numGroups;i++){
					set<string>::iterator aIt = groups[i].find(seqA);
					set<string>::iterator bIt = groups[i].find(seqB);
					
					if(groupIDA == -1 && aIt != groups[i].end()){//seqA is not already assigned to a group and is in group[i], so assign seqB to group[i]
						groups[i].insert(seqB);
						groupIDA = i;
						groupID = groupIDA;

						//cout << "in aIt: " << groupID << endl;
	//					break;
					}
					else if(groupIDB == -1 && bIt != groups[i].end()){//seqB is not already assigned to a group and is in group[i], so assign seqA to group[i]
						groups[i].insert(seqA);
						groupIDB = i;
						groupID = groupIDB;

					//	cout << "in bIt: " << groupID << endl;
	//					break;
					}
				
					if(groupIDA != -1 && groupIDB != -1){//both ifs above have been executed, so we need to decide who to assign them to
						if(groupIDA < groupIDB){
						//	cout << "A: " << groupIDA << "\t" << groupIDB << endl;
							groups[groupIDA].insert(groups[groupIDB].begin(), groups[groupIDB].end()); //merge two groups into groupIDA
							groups[groupIDB].clear(); 
							groupID = groupIDA;
						}
						else{
						//	cout << "B: " << groupIDA << "\t" << groupIDB << endl;
							groups[groupIDB].insert(groups[groupIDA].begin(), groups[groupIDA].end()); //merge two groups into groupIDB
							groups[groupIDA].clear();  
							groupID = groupIDB;
						}
						break;
					}
				}
				
	//windows is gonna gag on the reuse of outFile, will need to make it local...
				
				if(groupIDA == -1 && groupIDB == -1){ //we need a new group
					set<string> newGroup;
					newGroup.insert(seqA);
					newGroup.insert(seqB);
					groups.push_back(newGroup);
									
					outFile.close();
					string fileName = distFile + "." + toString(numGroups) + ".temp";
					outFile.open(fileName.c_str(), ios::ate);

					outFile << seqA << '\t' << seqB << '\t' << dist << endl;
					numGroups++;
				}
				else{
					string fileName = distFile + "." + toString(groupID) + ".temp";
					if(groupID != prevGroupID){
						outFile.close();
						outFile.open(fileName.c_str(), ios::app);
						prevGroupID	= groupID;
					}
					outFile << seqA << '\t' << seqB << '\t' << dist << endl;
					
					if(groupIDA != -1 && groupIDB != -1){ //merge distance files of two groups you merged above
						string row, column, distance;
						if(groupIDA<groupIDB){
							string fileName = distFile + "." + toString(groupIDB) + ".temp";
							ifstream fileB(fileName.c_str());
							while(fileB){
								fileB >> row >> column >> distance;
								outFile << row << '\t' << column << '\t' << distance << endl;
								gobble(fileB);
							}
							fileB.close();
							remove(fileName.c_str());
						}
						else{
							string fileName = distFile + "." + toString(groupIDA) + ".temp";
							ifstream fileA(fileName.c_str());
							while(fileA){
								fileA >> row >> column >> distance;
								outFile << row << '\t' << column << '\t' << distance << endl;
								gobble(fileA);
							}
							fileA.close();
							remove(fileName.c_str());
						}					
					}
				}
			}
			gobble(dFile);
		}
		outFile.close();
		dFile.close();
	
		ifstream bigNameFile(namefile.c_str());
		if(!bigNameFile){
			cerr << "Error: We can't open the name file\n";
			exit(1);
		}
		
		map<string, string> nameMap;
		string name, nameList;
		while(bigNameFile){
			bigNameFile >> name >> nameList;
			nameMap[name] = nameList;
			gobble(bigNameFile);
		}
		bigNameFile.close();
			
		for(int i=0;i<numGroups;i++){  //parse names file to match distance files
			int numSeqsInGroup = groups[i].size();
			
			if(numSeqsInGroup > 0){
				string fileName = namefile + "." + toString(i) + ".temp";
				ofstream smallNameFile(fileName.c_str(), ios::ate);
				
				for(set<string>::iterator gIt=groups[i].begin();gIt!=groups[i].end();gIt++){
					map<string,string>::iterator nIt = nameMap.find(*gIt);
					
					if (nIt != nameMap.end()) {
						smallNameFile << nIt->first << '\t' << nIt->second << endl;
						nameMap.erase(nIt);
					}else{
						m->mothurOut((*gIt) + " is in your distance file and not in your namefile.  Please correct."); m->mothurOutEndLine(); exit(1);
					}
				}
				smallNameFile.close();
			}
		}
			
		//names of singletons
		if (nameMap.size() != 0) {
			singleton = namefile + ".extra.temp";
			ofstream remainingNames(singleton.c_str(), ios::ate);
			for(map<string,string>::iterator nIt=nameMap.begin();nIt!=nameMap.end();nIt++){
				remainingNames << nIt->first << '\t' << nIt->second << endl;
			}
			remainingNames.close();
		}else { singleton = "none"; }
			
		for(int i=0;i<numGroups;i++){
			if(groups[i].size() > 0){
				string tempNameFile = namefile + "." + toString(i) + ".temp";
				string tempDistFile = distFile + "." + toString(i) + ".temp";
				
				map<string, string> temp;
				temp[tempDistFile] = tempNameFile;
				dists.push_back(temp);
			}
		}
				
		return 0;
			
	}
	catch(exception& e) {
		m->errorOut(e, "SplitMatrix", "split");
		exit(1);
	}
}
//********************************************************************************************************************
//sorts biggest to smallest
inline bool compareFileSizes(map<string, string> left, map<string, string> right){
	
	FILE * pFile;
	long leftsize = 0;
		
	//get num bytes in file
	string filename = left.begin()->first;
	pFile = fopen (filename.c_str(),"rb");
	string error = "Error opening " + filename;
	if (pFile==NULL) perror (error.c_str());
	else{
		fseek (pFile, 0, SEEK_END);
		leftsize=ftell (pFile);
		fclose (pFile);
	}

	FILE * pFile2;
	long rightsize = 0;
		
	//get num bytes in file
	filename = right.begin()->first;
	pFile2 = fopen (filename.c_str(),"rb");
	error = "Error opening " + filename;
	if (pFile2==NULL) perror (error.c_str());
	else{
		fseek (pFile2, 0, SEEK_END);
		rightsize=ftell (pFile2);
		fclose (pFile2);
	}

	return (leftsize > rightsize);	
} 
/***********************************************************************/
//returns map of distance files -> namefile sorted by distance file size
vector< map< string, string> > SplitMatrix::getDistanceFiles(){
	try {	
		
		sort(dists.begin(), dists.end(), compareFileSizes);
		
		return dists;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitMatrix", "getDistanceFiles");
		exit(1);
	}
}
/***********************************************************************/
SplitMatrix::~SplitMatrix(){}
/***********************************************************************/

