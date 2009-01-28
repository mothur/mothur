/*
 *  readmatrix.cpp
 *  
 *
 *  Created by Pat Schloss on 8/13/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

using namespace std;

#include <string>
#include <map>
#include "utilities.hpp"
#include "sparsematrix.hpp"
#include "progress.hpp"
#include "listvector.hpp"
#include "rabundvector.hpp"
#include <exception>

#include "readmatrix.hpp"


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
				if(nameMap->count(name)==0){	cout << "Error: Sequence '" << name << "' was not found in the names file, please correct" << endl; }
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

				reading = new Progress("Reading matrix:    ", nseqs * (nseqs - 1) / 2);
		
				int	index = 0;
		
				for(int i=1;i<nseqs;i++){
					fileHandle >> name;
					matrixNames.push_back(name);
	
					//there's A LOT of repeated code throughout this method...
					if(nameMap == NULL){
						list->set(i, name);
					
						for(int j=0;j<i;j++){
							fileHandle >> distance;
						
							if(distance < cutoff){
								PCell value(i, j, distance);
								D->addCell(value);
							}
							index++;
							reading->update(index);
						}
				
					}
					else{
						if(nameMap->count(name)==0){	cout << "Error: Sequence '" << name << "' was not found in the names file, please correct" << endl; }
				
						for(int j=0;j<i;j++){
							fileHandle >> distance;
						
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

				reading = new Progress("Reading matrix:    ", nseqs * nseqs);
			
				int index = nseqs;
		
				for(int i=1;i<nseqs;i++){
					fileHandle >> name;		
					matrixNames.push_back(name);
	
					if(nameMap == NULL){
						list->set(i, name);
						for(int j=0;j<nseqs;j++){
							fileHandle >> distance;
					
							if(distance < cutoff && j < i){
								PCell value(i, j, distance);
								D->addCell(value);
							}
							index++;
							reading->update(index);
						}
					
					}
					else{
						if(nameMap->count(name)==0){	cout << "Error: Sequence '" << name << "' was not found in the names file, please correct" << endl; }
				
						for(int j=0;j<nseqs;j++){
							fileHandle >> distance;
					
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
					cout << "missed something" << '\t' << nameMap->size() << endl;
				}
			}

		}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadPhylipMatrix class Function read. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadPhylipMatrix class function read. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

ReadPhylipMatrix::~ReadPhylipMatrix(){
	delete D;
	delete list;
}

/***********************************************************************/

ReadColumnMatrix::ReadColumnMatrix(string df) : distFile(df){
	
	successOpen = openInputFile(distFile, fileHandle);
	
}

/***********************************************************************/

void ReadColumnMatrix::read(NameAssignment* nameMap){
	try {		
	
			string firstName, secondName;
			float distance;
			int nseqs = nameMap->size();

			list = new ListVector(nameMap->getListVector());
		
			Progress* reading = new Progress("Reading matrix:    ", nseqs * nseqs);
	
			int lt = 1;
			int refRow = 0;	//we'll keep track of one cell - Cell(refRow,refCol) - and see if it's transpose
			int refCol = 0; //shows up later - Cell(refCol,refRow).  If it does, then its a square matrix
	
			//need to see if this is a square or a triangular matrix...
			while(fileHandle && lt == 1){  //let's assume it's a triangular matrix...
			
				fileHandle >> firstName >> secondName >> distance;	// get the row and column names and distance
		
				if(nameMap->count(firstName)==0){
					cerr << "AError: Sequence '" << firstName << "' was not found in the names file, please correct\n";
				}
				if(nameMap->count(secondName)==0){
					cerr << "AError: Sequence '" << secondName << "' was not found in the names file, please correct\n";
				}
		
				if(distance < cutoff && nameMap->get(firstName) != nameMap->get(secondName)){
					if(nameMap->get(firstName) > nameMap->get(secondName)){
						PCell value(nameMap->get(firstName), nameMap->get(secondName), distance);
				
						if(refRow == refCol){		// in other words, if we haven't loaded refRow and refCol...
							refRow = nameMap->get(firstName);
							refCol = nameMap->get(secondName);
							D->addCell(value);
						}
						else if(refRow == nameMap->get(firstName) && refCol == nameMap->get(secondName)){
							lt = 0;
						}
						else{
							D->addCell(value);
						}
					}
					else if(nameMap->get(firstName) < nameMap->get(secondName)){
						PCell value(nameMap->get(secondName), nameMap->get(firstName), distance);
				
						if(refRow == refCol){		// in other words, if we haven't loaded refRow and refCol...
							refRow = nameMap->get(firstName);
							refCol = nameMap->get(secondName);
							D->addCell(value);
						}
						else if(refRow == nameMap->get(secondName) && refCol == nameMap->get(firstName)){
							lt = 0;
						}
						else{
							D->addCell(value);
						}
					}
					reading->update(nameMap->get(firstName) * nseqs);
				}
				gobble(fileHandle);
			}

			if(lt == 0){  // oops, it was square
				fileHandle.close();  //let's start over
				D->clear();  //let's start over
			   
				openInputFile(distFile, fileHandle);  //let's start over

				while(fileHandle){
					fileHandle >> firstName >> secondName >> distance;
			
					if(nameMap->count(firstName)==0){
						cerr << "BError: Sequence '" << firstName << "' was not found in the names file, please correct\n";
					}
					if(nameMap->count(secondName)==0){
						cerr << "BError: Sequence '" << secondName << "' was not found in the names file, please correct\n";
					}
			
					if(distance < cutoff && nameMap->get(firstName) > nameMap->get(secondName)){
						PCell value(nameMap->get(firstName), nameMap->get(secondName), distance);
						D->addCell(value);
						reading->update(nameMap->get(firstName) * nseqs);
					}
			
					gobble(fileHandle);
				}
			}
		//	else if(lt == 0){
		//		while(fileHandle){
		//			fileHandle >> firstName >> secondName >> distance;
		//			
		//			if(nameMap->count(firstName)==0){
		//				cerr << "CError: Sequence '" << firstName << "' was not found in the names file, please correct\n";
		//			}
		//			if(nameMap->count(secondName)==0){
		//				cerr << "CError: Sequence '" << secondName << "' was not found in the names file, please correct\n";
		//			}
		//			
		//			if(distance < cutoff && (*nameMap)[firstName].second < (*nameMap)[secondName].second){
		////				cout << (*nameMap)[secondName] << ' ' << (*nameMap)[firstName] << ' ' << distance << endl;
		//				D->addCell(Cell((*nameMap)[secondName].second, (*nameMap)[firstName].second, distance));
		//				reading->update((*nameMap)[secondName].second * nseqs);
		//			}
		//
		//			gobble(fileHandle);
		//		}
		//	}	
			reading->finish();
			fileHandle.close();
	
			list->setLabel("0");
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadColumnMatrix class Function read. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadColumnMatrix class function read. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

/***********************************************************************/

ReadColumnMatrix::~ReadColumnMatrix(){
	delete D;
	delete list;
}


/***********************************************************************/

ReadPhilFile::ReadPhilFile(string pf): philFile(pf){
	
	successOpen = openInputFile(philFile, fileHandle);
	
}

/***********************************************************************/
//This function reads the list, rabund or sabund files to be used by collect and rarefact command.
void ReadPhilFile::read(GlobalData* globaldata){
	try {
		if (globaldata->getOrderFile() == "") {
			//you have two inputs because in the next if statement if you only have one then it moves ahead in the same file.  
			//So when you run the collect or summary commands you miss a line.
			input = new InputData(philFile, globaldata->getFormat()); //format tells you whether philFile is list, rabund, sabund.
			inputSabund = new InputData(philFile, globaldata->getFormat()); //format tells you whether philFile is list, rabund, sabund.
		}else {//there is an orderfile
			input = new InputData(philFile, globaldata->getOrderFile(), globaldata->getFormat());
		}
		globaldata->ginput = input;	//saving to be used by collector and rarefact commands.
		
		if ((globaldata->getFormat() == "list") || (globaldata->getFormat() == "rabund") || (globaldata->getFormat() == "sabund")) {//you are reading a list, rabund or sabund file for collect, rarefaction or summary.
			order = input->getOrderVector();
			globaldata->gorder = order;	//saving to be used by collect and rarefact commands.
			sabund = inputSabund->getSAbundVector(); 
			globaldata->sabund = sabund; //saving to be used by summary command.
		}else if (globaldata->getFormat() == "shared") {
			SharedList = input->getSharedListVector(); //you are reading for parselist command, or shared commands.
			globaldata->gSharedList = SharedList;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadPhilFile class Function read. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadPhilFile class function read. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

ReadPhilFile::~ReadPhilFile(){
//	delete input;
//	delete order;
}

/***********************************************************************/

