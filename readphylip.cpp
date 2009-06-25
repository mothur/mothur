/*
 *  readphylip.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/21/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "readphylip.h"
#include "progress.hpp"

/***********************************************************************/

ReadPhylipMatrix::ReadPhylipMatrix(string distFile){
        
        successOpen = openInputFile(distFile, fileHandle);
        
}

/***********************************************************************/

void ReadPhylipMatrix::read(NameAssignment* nameMap){
        try {
        
                        float distance;
                        int square, nseqs;
                        string name;
                        vector<string> matrixNames;
        
                        fileHandle >> nseqs >> name;

                        matrixNames.push_back(name);

                        if(nameMap == NULL){
                                list = new ListVector(nseqs);
                                list->set(0, name);
                        }
                        else{
                                list = new ListVector(nameMap->getListVector());
                                if(nameMap->count(name)==0){        mothurOut("Error: Sequence '" + name + "' was not found in the names file, please correct"); mothurOutEndLine(); }
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
                                        fileHandle >> name;
                                        matrixNames.push_back(name);
        
                                        //there's A LOT of repeated code throughout this method...
                                        if(nameMap == NULL){
                                                list->set(i, name);
                                        
                                                for(int j=0;j<i;j++){
                                                        fileHandle >> distance;
                                                
                                                        if (distance == -1) { distance = 1000000; }
                                                
                                                        if(distance < cutoff){
                                                                PCell value(i, j, distance);
                                                                D->addCell(value);
                                                        }
                                                        index++;
                                                        reading->update(index);
                                                }
                                
                                        }
                                        else{
                                                if(nameMap->count(name)==0){        mothurOut("Error: Sequence '" + name + "' was not found in the names file, please correct"); mothurOutEndLine(); }
                                
                                                for(int j=0;j<i;j++){
                                                        fileHandle >> distance;
                                
                                                        if (distance == -1) { distance = 1000000; }
                                                        
                                                        if(distance < cutoff){
                                                                PCell value(nameMap->get(matrixNames[i]), nameMap->get(matrixNames[j]), distance);
                                                                D->addCell(value);
                                                        }
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
                                        
                                                        if (distance == -1) { distance = 1000000; }
                                                        
                                                        if(distance < cutoff && j < i){
                                                                PCell value(i, j, distance);
                                                                D->addCell(value);
                                                        }
                                                        index++;
                                                        reading->update(index);
                                                }
                                        
                                        }
                                        else{
                                                if(nameMap->count(name)==0){        mothurOut("Error: Sequence '" + name + "' was not found in the names file, please correct"); mothurOutEndLine(); }
                                
                                                for(int j=0;j<nseqs;j++){
                                                        fileHandle >> distance;
                        
                                                        if (distance == -1) { distance = 1000000; }
                                                        
                                                        if(distance < cutoff && j < i){
                                                                PCell value(nameMap->get(matrixNames[i]), nameMap->get(matrixNames[j]), distance);
                                                                D->addCell(value);
                                                        }
                                                        index++;
                                                        reading->update(index);
                                                }
                                        }
                                }
                        }
                        reading->finish();
                        delete reading;

                        list->setLabel("0");
                        fileHandle.close();

                        if(nameMap != NULL){
                                for(int i=0;i<matrixNames.size();i++){
                                        nameMap->erase(matrixNames[i]);
                                }
                                if(nameMap->size() > 0){
                                        //should probably tell them what is missing if we missed something
                                        mothurOut("missed something\t" + toString(nameMap->size())); mothurOutEndLine();
                                }
                        }

                }
        catch(exception& e) {
               errorOut(e, "ReadPhylipMatrix", "read");
                exit(1);
        }
	}

/***********************************************************************/

ReadPhylipMatrix::~ReadPhylipMatrix(){
       // delete D;
       // delete list;
}
