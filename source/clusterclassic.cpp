/*
 *  clusterclassic.cpp
 *  Mothur
 *
 *  Created by westcott on 10/29/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "clusterclassic.h"
#include "utils.hpp"

/***********************************************************************/
ClusterClassic::ClusterClassic(float c, string f, bool s) : method(f),  nseqs(0), sim(s) {
	try {
        smallDist = MOTHURMAX;
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
        Utils util; util.openInputFile(filename, fileHandle);
		
        string numTest;
		fileHandle >> numTest >> name;
        
        if (!util.isContainingOnlyDigits(numTest)) { m->mothurOut("[ERROR]: expected a number and got " + numTest + ", quitting."); m->mothurOutEndLine(); exit(1); }
        else { convert(numTest, nseqs); }


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
		for (int i = 1; i < nseqs; i++) {	 dMatrix[i].resize(i, aboveCutoff);	 }
        
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

		

		if(square == 0){


				int        index = 0;

				for(int i=1;i<nseqs;i++){
						if (m->getControl_pressed()) {  fileHandle.close();   return 0; }
						
						fileHandle >> name;
						matrixNames.push_back(name);
		

						//there's A LOT of repeated code throughout this method...
						if(nameMap == NULL){
								list->set(i, name);
						
								for(int j=0;j<i;j++){
								
										if (m->getControl_pressed()) {  fileHandle.close(); return 0;  }
										
										fileHandle >> distance;
							
										if (util.isEqual(distance,-1)) { distance = 1000000; }
										else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
								
										
                                        dMatrix[i][j] = distance;
                                        if (distance < smallDist) { smallDist = distance; }
											
										index++;
										
								}
				
						}
						else{
								if(nameMap->count(name)==0){        m->mothurOut("Error: Sequence '" + name + "' was not found in the names file, please correct"); m->mothurOutEndLine(); }
				
								for(int j=0;j<i;j++){
										fileHandle >> distance;
										
										if (m->getControl_pressed()) {  fileHandle.close(); return 0;  }
				
										if (util.isEqual(distance,-1)) { distance = 1000000; }
										else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
										
                                    
                                        if (distance < smallDist) { smallDist = distance; }
                                    
                                        int row = nameMap->get(matrixNames[i]);
                                        int col = nameMap->get(matrixNames[j]);
                                    
                                        if (row < col) {  dMatrix[col][row] = distance; }
                                        else { dMatrix[row][col] = distance; }
											
										index++;
										
								}
						}
				}
		}
		else{

				
		
				int index = nseqs;

				for(int i=1;i<nseqs;i++){
						fileHandle >> name;                
						matrixNames.push_back(name);
						
						if(nameMap == NULL){
								list->set(i, name);
								for(int j=0;j<nseqs;j++){
										fileHandle >> distance;
										
										if (m->getControl_pressed()) {  fileHandle.close();  return 0; }
										
										if (util.isEqual(distance,-1)) { distance = 1000000; }
										else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
										
										if(j < i){
											if (distance < smallDist) { smallDist = distance; }
											
											dMatrix[i][j] = distance;
										}
										index++;
										
								}
						
						}
						else{
								if(nameMap->count(name)==0){        m->mothurOut("Error: Sequence '" + name + "' was not found in the names file, please correct"); m->mothurOutEndLine(); }
				
								for(int j=0;j<nseqs;j++){
										fileHandle >> distance;
										
										if (m->getControl_pressed()) {  fileHandle.close();   return 0; }
										
									   if (distance == -1) { distance = 1000000; }
										else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.                                                        
										
										if(j < i){
											if (distance < smallDist) { smallDist = distance; }
											
											int row = nameMap->get(matrixNames[i]);
											int col = nameMap->get(matrixNames[j]);
											
											if (row < col) {  dMatrix[col][row] = distance; }
											else { dMatrix[row][col] = distance; }
										}
										index++;
										
								}
						}
				}
		}
		
		if (m->getControl_pressed()) {  fileHandle.close();  return 0; }
		
		

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
int ClusterClassic::readPhylipFile(string filename, CountTable* countTable) {
	try {
		double distance;
		int square;
		string name;
		vector<string> matrixNames;
		
		ifstream fileHandle;
		Utils util; util.openInputFile(filename, fileHandle);
		
        string numTest;
		fileHandle >> numTest >> name;
        
        if (!util.isContainingOnlyDigits(numTest)) { m->mothurOut("[ERROR]: expected a number and got " + numTest + ", quitting."); m->mothurOutEndLine(); exit(1); }
        else { convert(numTest, nseqs); }
        
        
		matrixNames.push_back(name);
        
		if(countTable == NULL){
            list = new ListVector(nseqs);
            list->set(0, name);
        }
        else{  list = new ListVector(countTable->getListVector()); }

		
		//initialize distance matrix to cutoff
		dMatrix.resize(nseqs);
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
        
		
        
		if(square == 0){
            
            
            int        index = 0;
            
            for(int i=1;i<nseqs;i++){
                if (m->getControl_pressed()) {  fileHandle.close(); return 0; }
                
                fileHandle >> name;
                matrixNames.push_back(name);
                
                
                //there's A LOT of repeated code throughout this method...
                 if(countTable == NULL){
                    list->set(i, name);
                    
                    for(int j=0;j<i;j++){
                        
                        if (m->getControl_pressed()) { fileHandle.close(); return 0;  }
                        
                        fileHandle >> distance;
                        
                        if (util.isEqual(distance,-1)) { distance = 1000000; }
                        else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                        
                        
                        dMatrix[i][j] = distance;
                        if (distance < smallDist) { smallDist = distance; }
                                                index++;
                        
                    }
                    
                }
                else{
                    for(int j=0;j<i;j++){
                        fileHandle >> distance;
                        
                        if (m->getControl_pressed()) {  fileHandle.close(); return 0;  }
                        
                        if (util.isEqual(distance,-1)) { distance = 1000000; }
                        else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                        
                        if (distance < smallDist) { smallDist = distance; }
                        
                        int row = countTable->get(matrixNames[i]);
                        int col = countTable->get(matrixNames[j]);
                       
                        if (row < col) {  dMatrix[col][row] = distance; }
                        else { dMatrix[row][col] = distance; }
                        
                        index++;
                        
                    }
                }
            }
		}
		else{
            
            
            int index = nseqs;
            
            for(int i=1;i<nseqs;i++){
                fileHandle >> name;                
                matrixNames.push_back(name);
                
                if(countTable == NULL){
                    list->set(i, name);
                    for(int j=0;j<nseqs;j++){
                        fileHandle >> distance;
                        
                        if (m->getControl_pressed()) {  fileHandle.close();   return 0; }
                        
                        if (util.isEqual(distance,-1)) { distance = 1000000; }
                        else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                        
                        if(j < i){
                            if (distance < smallDist) { smallDist = distance; }
                            
                            dMatrix[i][j] = distance;
                        }
                        index++;
                    }
                    
                }
                else{
                    
                    for(int j=0;j<nseqs;j++){
                        fileHandle >> distance;
                        
                        if (m->getControl_pressed()) {  fileHandle.close();   return 0; }
                        
                        if (util.isEqual(distance,-1)) { distance = 1000000; }
                        else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.                                                        
                        
                        if(j < i){
                            if (distance < smallDist) { smallDist = distance; }
                            
                            int row = countTable->get(matrixNames[i]);
                            int col = countTable->get(matrixNames[j]);
                            
                            if (row < col) {  dMatrix[col][row] = distance; }
                            else { dMatrix[row][col] = distance; }
                        }
                        index++;
                    }
                }
            }
		}
		
		if (m->getControl_pressed()) {  fileHandle.close();   return 0; }
		
		
        
		list->setLabel("0");
		rabund = new RAbundVector();
        rabund->setLabel(list->getLabel());  
        
        for(int i = 0; i < list->getNumBins(); i++) { 
            if (m->getControl_pressed()) { break; }
            vector<string> binNames;
            string bin = list->get(i);
            util.splitAtComma(bin, binNames);
            int total = 0;
            for (int j = 0; j < binNames.size(); j++) { total += countTable->getNumSeqs(binNames[j]);  }
            rabund->push_back(total);   
        }
		
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
                    //cout << "adding " << i << '\t' << j << " to mins\n";
					smallDist = dMatrix[i][j];
				}else if (dMatrix[i][j] == smallDist) {
					colDist temp(i, j, dMatrix[i][j]);
					mins.push_back(temp);
                    //cout << "adding " << i << '\t' << j << " to mins\n";

				}
			}
		}
        
        if (mins.size() != 0) {
            Utils util; util.mothurRandomShuffle(mins);  //randomize the order of the iterators in the mins vector
        
            smallRow = mins[0].row;
            smallCol = mins[0].col;
        }
        
        //cout << "num mins = " << mins.size() << endl;
        
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
		rabund->set(smallRow, rabund->get(smallRow)+rabund->get(smallCol));	
		rabund->set(smallCol, 0);
		rabund->setLabel(toString(smallDist));
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterClassic", "clusterBins");
		exit(1);
	}
}
/***********************************************************************/
void ClusterClassic::clusterNames(){
	try {
		if (mapWanted) {  updateMap();  }
		
		list->set(smallRow, list->get(smallRow)+','+list->get(smallCol));
		list->set(smallCol, "");
		list->setLabel(toString(smallDist));
    }
	catch(exception& e) {
		m->errorOut(e, "ClusterClassic", "clusterNames");
		exit(1);
	}
}
/***********************************************************************/
void ClusterClassic::update(double& cutOFF){
	try {
		
		getSmallCell();
		
		int r, c;
		r = smallRow; c = smallCol;
        
        //cout << "smallest row / col = " << smallRow << '\t' << smallCol << endl;
				
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
					
				if (i > r) { dMatrix[i][r] = newDist; }
				else { dMatrix[r][i] = newDist; }
				
			}
		}
			
		clusterBins();
		clusterNames();
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
        Utils util;
		
		//initialize map
		for (int i = 0; i < list->getNumBins(); i++) {
			
			//parse bin 
			string names = list->get(i);
			vector<string> binnames;
            util.splitAtComma(names, binnames);
            for (int j = 0; j < binnames.size(); j++) {
				//save name and bin number
				seq2Bin[binnames[j]] = i;
			}
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
        vector<string> binnames;
        Utils util; util.splitAtComma(names, binnames);
        for (int j = 0; j < binnames.size(); j++) {
            //save name and bin number
            seq2Bin[binnames[j]] = smallCol;
        }
		
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

