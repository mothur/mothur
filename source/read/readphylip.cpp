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
        
        successOpen = m->openInputFile(distFile, fileHandle);
		sim=false;
        
}
/***********************************************************************/

ReadPhylipMatrix::ReadPhylipMatrix(string distFile, bool s){
	
	successOpen = m->openInputFile(distFile, fileHandle);
	sim=s;
}


/***********************************************************************/

int ReadPhylipMatrix::read(NameAssignment* nameMap){
        try {
        
                        float distance;
                        int square, nseqs; 
                        string name;
                        vector<string> matrixNames;
						
						string numTest;
						fileHandle >> numTest >> name;
			
						if (!m->isContainingOnlyDigits(numTest)) { m->mothurOut("[ERROR]: expected a number and got " + numTest + ", quitting."); m->mothurOutEndLine(); exit(1); }
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
                        DMatrix->resize(nseqs);
        
                        if(square == 0){

                                reading = new Progress("Reading matrix:     ", nseqs * (nseqs - 1) / 2);
                
                                int        index = 0;
               
                                for(int i=1;i<nseqs;i++){
										if (m->getControl_pressed()) {  fileHandle.close();  delete reading; return 0; }
										
                                        fileHandle >> name;
                                        matrixNames.push_back(name);
						
        
                                        //there's A LOT of repeated code throughout this method...
                                        if(nameMap == NULL){
                                                list->set(i, name);
                                        
                                                for(int j=0;j<i;j++){
												
														if (m->getControl_pressed()) { delete reading; fileHandle.close(); return 0;  }
														
                                                        fileHandle >> distance;
                                                
                                                        if (distance == -1) { distance = 1000000; }
														else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                                                
                                                        if(distance < cutoff){
                                                            PDistCell value(i, distance);
                                                            DMatrix->addCell(j, value);
                                                        }
                                                        index++;
                                                        reading->update(index);
                                                }
                                
                                        }
                                        else{
                                                if(nameMap->count(name)==0){        m->mothurOut("Error: Sequence '" + name + "' was not found in the names file, please correct"); m->mothurOutEndLine(); }
                                
                                                for(int j=0;j<i;j++){
                                                        fileHandle >> distance;
														
														if (m->getControl_pressed()) { delete reading; fileHandle.close(); return 0;  }
                                
                                                        if (distance == -1) { distance = 1000000; }
														else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                                                        
                                                        if(distance < cutoff){
                                                            PDistCell value(nameMap->get(matrixNames[i]), distance);
                                                            DMatrix->addCell(nameMap->get(matrixNames[j]), value);
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
														
														if (m->getControl_pressed()) {  fileHandle.close();  delete reading; return 0; }
														
                                                        if (distance == -1) { distance = 1000000; }
														else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                                                        
                                                        if(distance < cutoff && j < i){
                                                            PDistCell value(i, distance);
                                                            DMatrix->addCell(j, value);
                                                        }
                                                        index++;
                                                        reading->update(index);
                                                }
                                        
                                        }
                                        else{
                                                if(nameMap->count(name)==0){        m->mothurOut("Error: Sequence '" + name + "' was not found in the names file, please correct"); m->mothurOutEndLine(); }
                                
                                                for(int j=0;j<nseqs;j++){
                                                        fileHandle >> distance;
														
														if (m->getControl_pressed()) {  fileHandle.close();  delete reading; return 0; }
														
                                                       if (distance == -1) { distance = 1000000; }
														else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.                                                        
                                                        
														if(distance < cutoff && j < i){
                                                            PDistCell value(nameMap->get(matrixNames[i]), distance);
                                                            DMatrix->addCell(nameMap->get(matrixNames[j]), value);
                                                        }
                                                        index++;
                                                        reading->update(index);
                                                }
                                        }
                                }
                        }
						
						if (m->getControl_pressed()) {  fileHandle.close();  delete reading; return 0; }
						
                        reading->finish();
                        delete reading;
            
                        list->setLabel("0");
                        fileHandle.close();

                    						
						return 1;

                }
        catch(exception& e) {
               m->errorOut(e, "ReadPhylipMatrix", "read");
                exit(1);
        }
	}
/***********************************************************************/

int ReadPhylipMatrix::read(CountTable* countTable){
    try {
        
        float distance;
        int square, nseqs; 
        string name;
        vector<string> matrixNames;
        
        string numTest;
        fileHandle >> numTest >> name;
        
        if (!m->isContainingOnlyDigits(numTest)) { m->mothurOut("[ERROR]: expected a number and got " + numTest + ", quitting."); m->mothurOutEndLine(); exit(1); }
        else { convert(numTest, nseqs); }
        
        matrixNames.push_back(name);
        
        if(countTable == NULL){
            list = new ListVector(nseqs);
            list->set(0, name);
        }
        else{  list = new ListVector(countTable->getListVector()); }
        
        if (m->getControl_pressed()) { return 0; }
        
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
        DMatrix->resize(nseqs);
        
        if(square == 0){
            
            reading = new Progress("Reading matrix:     ", nseqs * (nseqs - 1) / 2);
            
            int        index = 0;
            
            for(int i=1;i<nseqs;i++){
                if (m->getControl_pressed()) {  fileHandle.close();  delete reading; return 0; }
                
                fileHandle >> name;
                matrixNames.push_back(name);
                
                
                //there's A LOT of repeated code throughout this method...
                if(countTable == NULL){
                    list->set(i, name);
                    
                    for(int j=0;j<i;j++){
                        
                        if (m->getControl_pressed()) { delete reading; fileHandle.close(); return 0;  }
                        
                        fileHandle >> distance;
                        
                        if (distance == -1) { distance = 1000000; }
                        else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                        
                        if(distance < cutoff){
                            PDistCell value(i, distance);
                            DMatrix->addCell(j, value);
                        }
                        index++;
                        reading->update(index);
                    }
                    
                }
                else{
                    for(int j=0;j<i;j++){
                        fileHandle >> distance;
                        
                        if (m->getControl_pressed()) { delete reading; fileHandle.close(); return 0;  }
                        
                        if (distance == -1) { distance = 1000000; }
                        else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                        
                        if(distance < cutoff){
                            int iIndex = countTable->get(matrixNames[i]);
                            int jIndex = countTable->get(matrixNames[j]);
                            
                            if (m->getControl_pressed()) { delete reading; fileHandle.close(); return 0;  }
                            if (iIndex < jIndex) { 
                                PDistCell value(jIndex, distance);
                                DMatrix->addCell(iIndex, value);
                            }else {
                                PDistCell value(iIndex, distance);
                                DMatrix->addCell(jIndex, value);

                            }
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
                
                
                
                if(countTable == NULL){
                    list->set(i, name);
                    for(int j=0;j<nseqs;j++){
                        fileHandle >> distance;
                        
                        if (m->getControl_pressed()) {  fileHandle.close();  delete reading; return 0; }
                        
                        if (distance == -1) { distance = 1000000; }
                        else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.
                        
                        if(distance < cutoff && j < i){
                            PDistCell value(i, distance);
                            DMatrix->addCell(j, value);
                        }
                        index++;
                        reading->update(index);
                    }
                    
                }
                else{
                    for(int j=0;j<nseqs;j++){
                        fileHandle >> distance;
                        
                        if (m->getControl_pressed()) {  fileHandle.close();  delete reading; return 0; }
                        
                        if (distance == -1) { distance = 1000000; }
                        else if (sim) { distance = 1.0 - distance;  }  //user has entered a sim matrix that we need to convert.                                                        
                        
                        if(distance < cutoff && j < i){
                            int iIndex = countTable->get(matrixNames[i]);
                            int jIndex = countTable->get(matrixNames[j]);
                            
                            if (m->getControl_pressed()) { delete reading; fileHandle.close(); return 0;  }
                            if (iIndex < jIndex) { 
                                PDistCell value(jIndex, distance);
                                DMatrix->addCell(iIndex, value);
                            }else {
                                PDistCell value(iIndex, distance);
                                DMatrix->addCell(jIndex, value);
                                
                            }
                        }
                        index++;
                        reading->update(index);
                    }
                }
            }
        }
        
        if (m->getControl_pressed()) {  fileHandle.close();  delete reading; return 0; }
        
        reading->finish();
        delete reading;
        
        list->setLabel("0");
        fileHandle.close();
        
        
        return 1;
        
    }
    catch(exception& e) {
        m->errorOut(e, "ReadPhylipMatrix", "read");
        exit(1);
    }
}
/***********************************************************************/
ReadPhylipMatrix::~ReadPhylipMatrix(){}
/***********************************************************************/

