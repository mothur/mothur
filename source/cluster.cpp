/*
 *  cluster.cpp
 *  
 *
 *  Created by Pat Schloss on 8/14/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

#include "cluster.hpp"
#include "rabundvector.hpp"
#include "listvector.hpp"

/***********************************************************************/

Cluster::Cluster(RAbundVector* rav, ListVector* lv, SparseDistanceMatrix* dm, float c, string f, float cs) :
rabund(rav), list(lv), dMatrix(dm), method(f), adjust(cs)
{
	try {
        
        mapWanted = false;  //set to true by mgcluster to speed up overlap merge
        
        //save so you can modify as it changes in average neighbor
        cutoff = c;
        m = MothurOut::getInstance();
	}
	catch(exception& e) {
		m->errorOut(e, "Cluster", "Cluster");
		exit(1);
	}
}
/***********************************************************************/
void Cluster::clusterBins(){
	try {
 		rabund->set(smallCol, rabund->get(smallRow)+rabund->get(smallCol));	
		rabund->set(smallRow, 0);	
		rabund->setLabel(toString(smallDist));
	}
	catch(exception& e) {
		m->errorOut(e, "Cluster", "clusterBins");
		exit(1);
	}
}
/***********************************************************************/

void Cluster::clusterNames(){
	try {
		if (mapWanted) {  updateMap();  }
		
		list->set(smallCol, list->get(smallRow)+','+list->get(smallCol));
		list->set(smallRow, "");	
		list->setLabel(toString(smallDist));
    }
	catch(exception& e) {
		m->errorOut(e, "Cluster", "clusterNames");
		exit(1);
	}
}
/***********************************************************************/
bool Cluster::update(double& cutOFF){
	try {
        smallCol = dMatrix->getSmallestCell(smallRow);
        nColCells = dMatrix->seqVec[smallCol].size();
        nRowCells = dMatrix->seqVec[smallRow].size();
        
		vector<int> foundCol(nColCells, 0);
        
		int search;
		bool changed = false;
        
		for (int i=nRowCells-1;i>=0;i--) {  //matrix indexes sorted from largest to smallest, so start at smallest index
            if (m->getControl_pressed()) { break; }
             
			//if you are not the smallCell
			if (dMatrix->seqVec[smallRow][i].index != smallCol) { 
                search = dMatrix->seqVec[smallRow][i].index;
                
				bool merged = false;
				for (int j=0;j<nColCells;j++) {  //go through each distance the smallCol has looking for matching distance to find
                    
					if (dMatrix->seqVec[smallCol][j].index != smallRow) {  //if you are not the smallest distance
						if (dMatrix->seqVec[smallCol][j].index == search) {  //we found a distance for the merge
							foundCol[j] = 1;
							merged = true;
							changed = updateDistance(dMatrix->seqVec[smallCol][j], dMatrix->seqVec[smallRow][i]);
                            dMatrix->updateCellCompliment(smallCol, j);
							break;
						}else if (dMatrix->seqVec[smallCol][j].index < search) { //we don't have a distance for this cell
                            if (!util.isEqual(adjust, -1)) { //adjust
                                merged = true;
                                PDistCell value(search, adjust); //create a distance for the missing value
                                int location = dMatrix->addCellSorted(smallCol, value);
                                changed = updateDistance(dMatrix->seqVec[smallCol][location], dMatrix->seqVec[smallRow][i]);
                                dMatrix->updateCellCompliment(smallCol, location);
                                nColCells++;
                                foundCol.push_back(0); //add a new found column
                                //adjust value
                                for (int k = foundCol.size()-1; k > location; k--) { foundCol[k] = foundCol[k-1]; }
                                foundCol[location] = 1;
                            }
                            j+=nColCells;  //jump out of loop and remove cell below
                        } 
                    }
				}
				//if not merged it you need it for warning 
				if ((!merged) && (method == "average" || method == "weighted")) {   if (cutOFF > dMatrix->seqVec[smallRow][i].dist) {   cutOFF = dMatrix->seqVec[smallRow][i].dist; } }
                if ((method == "nearest") && (!merged)) { //you are a row dist without a column dist, add you as a column dist
                    PDistCell value(search, dMatrix->seqVec[smallRow][i].dist); //create a distance for the missing value
                    int location = dMatrix->addCellSorted(smallCol, value); nColCells++;
                    foundCol.push_back(0); //add a new found column
                    //adjust value
                    for (int k = foundCol.size()-1; k > location; k--) { foundCol[k] = foundCol[k-1]; }
                    foundCol[location] = 1;
                }
                dMatrix->rmCell(smallRow, i);  
			}
		}
		clusterBins();
		clusterNames();
        
        if (method == "nearest") {
            for (int i=nColCells-1;i>=0;i--) { //remove any unfound dists from merged column, need special case for nn, since unfound dists mean above the cutoff -> keep smaller dist in col
                if (foundCol[i] == 0) {  //not found
                    if (dMatrix->seqVec[smallCol][i].index == smallRow) { //you are smallest distance
                        dMatrix->rmCell(smallCol, i);
                        break;
                    }
                }
            }
        }else {
            for (int i=nColCells-1;i>=0;i--) { //remove any unfound dists from merged column, need special case for nn, since unfound dists mean above the cutoff -> keep smaller dist in col
                if (foundCol[i] == 0) {  //not found
                    if (!util.isEqual(adjust, -1)) { //adjust
                        PDistCell value(smallCol, adjust); //create a distance for the missing value
                        changed = updateDistance(dMatrix->seqVec[smallCol][i], value);
                        dMatrix->updateCellCompliment(smallCol, i);
                    }else {
                        if (method == "average" || method == "weighted") {
                            if (dMatrix->seqVec[smallCol][i].index != smallRow) { //if you are not hte smallest distance
                                if (cutOFF > dMatrix->seqVec[smallCol][i].dist) {   cutOFF = dMatrix->seqVec[smallCol][i].dist;   }
                            }
                        }
                    }
                    dMatrix->rmCell(smallCol, i);
                }
            }
        }
        //dMatrix->print();
        return changed;
	}
	catch(exception& e) {
		m->errorOut(e, "Cluster", "update");
		exit(1);
	}
}
/***********************************************************************/
void Cluster::setMapWanted(bool f)  {  
	try {
		mapWanted = f;
		
        //initialize map
		for (int k = 0; k < list->getNumBins(); k++) {
            
            string names = list->get(k);
            
            //parse bin
            string individual = "";
            int binNameslength = names.size();
            for(int j=0;j<binNameslength;j++){
                if(names[j] == ','){
                    seq2Bin[individual] = k;
                    individual = "";				
                }
                else{  individual += names[j];  }
            }
            //get last name
            seq2Bin[individual] = k;
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "Cluster", "setMapWanted");
		exit(1);
	}
}
/***********************************************************************/
void Cluster::updateMap() {
    try {
		//update location of seqs in smallRow since they move to smallCol now
		string names = list->get(smallRow);
		
        string individual = "";
        int binNameslength = names.size();
        for(int j=0;j<binNameslength;j++){
            if(names[j] == ','){
                seq2Bin[individual] = smallCol;
                individual = "";				
            }
            else{  individual += names[j];  }
        }
        //get last name
        seq2Bin[individual] = smallCol;		
	
	}
	catch(exception& e) {
		m->errorOut(e, "Cluster", "updateMap");
		exit(1);
	}
}
/***********************************************************************/



