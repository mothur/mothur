/*
 *  clusterclassic.cpp
 *  Mothur
 *
 *  Created by westcott on 10/29/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "clusterclassic.h"
#include "progress.hpp"

/***********************************************************************/
ClusterClassic::ClusterClassic(float c, string f, bool s) : method(f), smallDist(1e6), nseqs(0), sim(s) {
	try {
		mapWanted = false;  //set to true by mgcluster to speed up overlap merge
	
		//save so you can modify as it changes in average neighbor
		cutoff = c;
		aboveCutoff = cutoff + 10000.0;
		m = MothurOut::getInstance();
        if(method == "furthest")        {   tag = "fn";   }
        else if (method == "average")   {   tag = "an";   }
        else if (method == "weighted")  {   tag = "wn";   }        
        else if (method == "nearest")   {   tag = "nn";   }
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterClassic", "ClusterClassic");
		exit(1);
	}
}
/***********************************************************************/
int ClusterClassic::readPhylipFile(string filename, NameAssignment* nameMap) {
	try {
		double distance;
		int square;
		string name;
		vector<string> matrixNames;
		
		ifstream fileHandle;
		m->openInputFile(filename, fileHandle);
		
		fileHandle >> nseqs >> name;

		matrixNames.push_back(name);

		if(nameMap == NULL){
				list = new ListVector(nseqs);
				list->set(0, name);
		}
		else{
				list = new ListVector(nameMap->getListVector());
				if(nameMap->count(name)==0){        m->mothurOut("Error: Sequence '" + name + "' was not found in the names file, please correct"); m->mothurOutEndLine(); }
		}
		
		//initialize distance matrix to cutoff
		dMatrix.resize(nseqs);
		//colDist temp(0, 0, aboveCutoff);
		//rowSmallDists.resize(nseqs, temp);
		for (int i = 1; i < nseqs; i++) {			
			dMatrix[i].resize(i, aboveCutoff);		
		}												
														
		
		char d;
		while((d=fileHandle.get()) != EOF){

				if(isalnum(d)){
						square = 1;
						fileHandle.putback(d);
						for(int i=0;i<nseqs;i++){
								fileHandle >> distance;
						}
						break;
				}
				if(d == '\n'){
						square = 0;
						break;
				}
		}

		Progress* reading;

		if(square == 0){

				reading = new Progress("Reading matrix:     ", nseqs * (nseqs - 1) / 2);

				int        index = 0;

				for(int i=1;i<nseqs;i++){
						if (m->control_pressed) {  fileHandle.close();  delete reading; return 0; }
						
						fileHandle >> name;
						matrixNames.push_back(name);
		

						//there's A LOT of repeated code throughout this method...
						if(nameMap == NULL){
								list->set(i, name);
						
								for(int j=0;j<i;j++){
								
										if (m->control_pressed) { delete reading; fileHandle.close(); return 0;  }
										
										fileHandle >> distance;
							
										if (distance == -1) { distance = 1000000; }
										else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
								
										//if(distance < cutoff){
											dMatrix[i][j] = distance;
											if (distance < smallDist) { smallDist = distance; }
											//if (rowSmallDists[i].dist > distance) {  rowSmallDists[i].dist = distance; rowSmallDists[i].col = j; rowSmallDists[i].row = i; }
											//if (rowSmallDists[j].dist > distance) {  rowSmallDists[j].dist = distance; rowSmallDists[j].col = i; rowSmallDists[j].row = j; }
										//}
										index++;
										reading->update(index);
								}
				
						}
						else{
								if(nameMap->count(name)==0){        m->mothurOut("Error: Sequence '" + name + "' was not found in the names file, please correct"); m->mothurOutEndLine(); }
				
								for(int j=0;j<i;j++){
										fileHandle >> distance;
										
										if (m->control_pressed) { delete reading; fileHandle.close(); return 0;  }
				
										if (distance == -1) { distance = 1000000; }
										else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
										
										//if(distance < cutoff){
											if (distance < smallDist) { smallDist = distance; }
											
											int row = nameMap->get(matrixNames[i]);
											int col = nameMap->get(matrixNames[j]);
											
											if (row < col) {  dMatrix[col][row] = distance; }
											else { dMatrix[row][col] = distance; }
											
											//if (rowSmallDists[row].dist > distance) {  rowSmallDists[row].dist = distance; rowSmallDists[row].col = col; rowSmallDists[row].row = row; }
											//if (rowSmallDists[col].dist > distance) {  rowSmallDists[col].dist = distance; rowSmallDists[col].col = row; rowSmallDists[col].row = col; }
										//}
										index++;
										reading->update(index);
								}
						}
				}
		}
		else{

				reading = new Progress("Reading matrix:     ", nseqs * nseqs);
		
				int index = nseqs;

				for(int i=1;i<nseqs;i++){
						fileHandle >> name;                
						matrixNames.push_back(name);
						
						if(nameMap == NULL){
								list->set(i, name);
								for(int j=0;j<nseqs;j++){
										fileHandle >> distance;
										
										if (m->control_pressed) {  fileHandle.close();  delete reading; return 0; }
										
										if (distance == -1) { distance = 1000000; }
										else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
										
										if(j < i){
											if (distance < smallDist) { smallDist = distance; }
											
											dMatrix[i][j] = distance;
											//if (rowSmallDists[i].dist > distance) {  rowSmallDists[i].dist = distance; rowSmallDists[i].col = j; rowSmallDists[i].row = i; }
											//if (rowSmallDists[j].dist > distance) {  rowSmallDists[j].dist = distance; rowSmallDists[j].col = i; rowSmallDists[j].row = j; }
										}
										index++;
										reading->update(index);
								}
						
						}
						else{
								if(nameMap->count(name)==0){        m->mothurOut("Error: Sequence '" + name + "' was not found in the names file, please correct"); m->mothurOutEndLine(); }
				
								for(int j=0;j<nseqs;j++){
										fileHandle >> distance;
										
										if (m->control_pressed) {  fileHandle.close();  delete reading; return 0; }
										
									   if (distance == -1) { distance = 1000000; }
										else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.                                                        
										
										if(j < i){
											if (distance < smallDist) { smallDist = distance; }
											
											int row = nameMap->get(matrixNames[i]);
											int col = nameMap->get(matrixNames[j]);
											
											if (row < col) {  dMatrix[col][row] = distance; }
											else { dMatrix[row][col] = distance; }
											
											//if (rowSmallDists[row].dist > distance) {  rowSmallDists[row].dist = distance; rowSmallDists[row].col = col; rowSmallDists[row].row = row; }
											//if (rowSmallDists[col].dist > distance) {  rowSmallDists[col].dist = distance; rowSmallDists[col].col = row; rowSmallDists[col].row = col; }
										}
										index++;
										reading->update(index);
								}
						}
				}
		}
		
		if (m->control_pressed) {  fileHandle.close();  delete reading; return 0; }
		
		reading->finish();
		delete reading;

		list->setLabel("0");
		rabund = new RAbundVector(list->getRAbundVector());
		
		fileHandle.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterClassic", "readPhylipFile");
		exit(1);
	}

}
/***********************************************************************/
//sets smallCol and smallRow, returns distance
double ClusterClassic::getSmallCell() {
	try {
			
		smallDist = aboveCutoff;
		smallRow = 1;
		smallCol = 0;
		
		vector<colDist> mins;
		
		for(int i=1;i<nseqs;i++){
			for(int j=0;j<i;j++){ 
				if (dMatrix[i][j] < smallDist) {
					mins.clear();
					colDist temp(i, j, dMatrix[i][j]);
					mins.push_back(temp); 
					smallDist = dMatrix[i][j];
				}else if (dMatrix[i][j] == smallDist) {
					colDist temp(i, j, dMatrix[i][j]);
					mins.push_back(temp); 

				}
			}
		}
		
		if (mins.size() > 0) {
			int zrand = 0;
			if (mins.size() > 1) {
				//pick random number between 0 and mins.size()
				zrand = (int)((float)(rand()) / (RAND_MAX / (mins.size()-1) + 1));
			}
			
			smallRow = mins[zrand].row;
			smallCol = mins[zrand].col;
		
		}
	//cout << smallRow << '\t' << smallCol << '\t' << smallDist << endl;
		//eliminate smallCell
		if (smallRow < smallCol) { dMatrix[smallCol][smallRow] = aboveCutoff; }
		else { dMatrix[smallRow][smallCol] = aboveCutoff; }
		
		return smallDist;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterClassic", "getSmallCell");
		exit(1);
	}
}
/***********************************************************************/
void ClusterClassic::clusterBins(){
	try {
	//	cout << smallCol << '\t' << smallRow << '\t' << smallDist << '\t' << rabund->get(smallRow) << '\t' << rabund->get(smallCol);

		rabund->set(smallRow, rabund->get(smallRow)+rabund->get(smallCol));	
		rabund->set(smallCol, 0);	
		/*for (int i = smallCol+1; i < rabund->size(); i++) {
			rabund->set((i-1), rabund->get(i));
		}
		rabund->resize((rabund->size()-1));*/
		rabund->setLabel(toString(smallDist));

	//	cout << '\t' << rabund->get(smallRow) << '\t' << rabund->get(smallCol) << endl;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterClassic", "clusterBins");
		exit(1);
	}
}
/***********************************************************************/
void ClusterClassic::clusterNames(){
	try {
	//	cout << smallCol << '\t' << smallRow << '\t' << smallDist << '\t' << list->get(smallRow) << '\t' << list->get(smallCol);
		if (mapWanted) {  updateMap();  }
		
		list->set(smallRow, list->get(smallRow)+','+list->get(smallCol));
		list->set(smallCol, "");	
		/*for (int i = smallCol+1; i < list->size(); i++) {
			list->set((i-1), list->get(i));
		}
		list->resize((list->size()-1));*/
		list->setLabel(toString(smallDist));
	
	//	cout << '\t' << list->get(smallRow) << '\t' << list->get(smallCol) << endl;
    }
	catch(exception& e) {
		m->errorOut(e, "ClusterClassic", "clusterNames");
		exit(1);
	}
}
/***********************************************************************/
void ClusterClassic::update(double& cutOFF){
	try {
//print();		
		getSmallCell();
		
		int r, c;
		r = smallRow; c = smallCol;
				
		for(int i=0;i<nseqs;i++){
			if(i != r && i != c){
				double distRow, distCol, newDist;
				if (i > r) { distRow = dMatrix[i][r]; }
				else { distRow =  dMatrix[r][i]; }

				if (i > c) { distCol = dMatrix[i][c]; dMatrix[i][c] = aboveCutoff; } //like removeCell
				else { distCol =  dMatrix[c][i]; dMatrix[c][i] = aboveCutoff; }
				
				if(method == "furthest"){
					newDist = max(distRow, distCol);
				}
				else if (method == "average"){
					int rowBin = rabund->get(r);
					int colBin = rabund->get(c);
					newDist = (colBin * distCol + rowBin * distRow) / (rowBin + colBin);
				}
				else if (method == "weighted"){
					newDist = (distCol + distRow) / 2.0;
				}
				else if (method == "nearest"){
					newDist = min(distRow, distCol);
				}
				//cout << "newDist = " << newDist << endl;	
				if (i > r) { dMatrix[i][r] = newDist; }
				else { dMatrix[r][i] = newDist; }
				
			}
		}
			
		clusterBins();
		clusterNames();
		
		//resize each row
		/*for(int i=0;i<nseqs;i++){
			for(int j=c+1;j<dMatrix[i].size();j++){
				dMatrix[i][j-1]=dMatrix[i][j];
			}
		}			
		
		//resize each col
		for(int i=c+1;i<nseqs;i++){
			for(int j=0;j < dMatrix[i-1].size();j++){
				dMatrix[i-1][j]=dMatrix[i][j];
			}
		}	
		
		nseqs--;
		dMatrix.pop_back();*/

	}
	catch(exception& e) {
		m->errorOut(e, "ClusterClassic", "update");
		exit(1);
	}
}
/***********************************************************************/
void ClusterClassic::setMapWanted(bool f)  {  
	try {
		mapWanted = f;
		
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
		m->errorOut(e, "ClusterClassic", "setMapWanted");
		exit(1);
	}
}
/***********************************************************************/
void ClusterClassic::updateMap() {
try {
		//update location of seqs in smallRow since they move to smallCol now
		string names = list->get(smallRow);
		while (names.find_first_of(',') != -1) { 
			//get name from bin
			string name = names.substr(0,names.find_first_of(','));
			//save name and bin number
			seq2Bin[name] = smallCol;
			names = names.substr(names.find_first_of(',')+1, names.length());
		}
			
		//get last name
		seq2Bin[names] = smallCol;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterClassic", "updateMap");
		exit(1);
	}
}
/***********************************************************************/
void ClusterClassic::print() {
try {
		//update location of seqs in smallRow since they move to smallCol now
		for (int i = 0; i < dMatrix.size(); i++) {
			m->mothurOut("row = " + toString(i) + "\t");
			for (int j = 0; j < dMatrix[i].size(); j++) {
				m->mothurOut(toString(dMatrix[i][j]) + "\t");
			}
			m->mothurOutEndLine();
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterClassic", "updateMap");
		exit(1);
	}
}
/***********************************************************************/

