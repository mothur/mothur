/*
 *  formatcolumn.cpp
 *  Mothur
 *
 *  Created by westcott on 1/13/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "formatcolumn.h"
#include "progress.hpp"

/***********************************************************************/
FormatColumnMatrix::FormatColumnMatrix(string df) : filename(df){
	openInputFile(filename, fileHandle);
}
/***********************************************************************/

void FormatColumnMatrix::read(NameAssignment* nameMap){
	try {		

		string firstName, secondName;
		float distance;
		int nseqs = nameMap->size();

		list = new ListVector(nameMap->getListVector());
	
		Progress* reading = new Progress("Formatting matrix:     ", nseqs * nseqs);

		int lt = 1;
		int refRow = 0;	//we'll keep track of one cell - Cell(refRow,refCol) - and see if it's transpose
		int refCol = 0; //shows up later - Cell(refCol,refRow).  If it does, then its a square matrix

		//need to see if this is a square or a triangular matrix...
		
		ofstream out;
		string tempOutFile = filename + ".temp";
		openOutputFile(tempOutFile, out);
	
		while(fileHandle && lt == 1){  //let's assume it's a triangular matrix...
		
			fileHandle >> firstName >> secondName >> distance;	// get the row and column names and distance
	
			map<string,int>::iterator itA = nameMap->find(firstName);
			map<string,int>::iterator itB = nameMap->find(secondName);
			if(itA == nameMap->end()){	cerr << "AAError: Sequence '" << firstName << "' was not found in the names file, please correct\n"; exit(1);	}
			if(itB == nameMap->end()){	cerr << "ABError: Sequence '" << secondName << "' was not found in the names file, please correct\n"; exit(1);	}

			if (distance == -1) { distance = 1000000; }
		
			if((distance < cutoff) && (itA != itB)){
				if(refRow == refCol){		// in other words, if we haven't loaded refRow and refCol...
					refRow = itA->second;
					refCol = itB->second;
					
					//making it square
					out << itA->second << '\t' << itB->second << '\t' << distance << endl;
					out << itB->second << '\t' << itA->second << '\t' << distance << endl;
				}
				else if(refRow == itA->second && refCol == itB->second){	lt = 0;		} //you are square
				else if(refRow == itB->second && refCol == itA->second){	lt = 0;		}  //you are square
				else{	//making it square
					out << itA->second << '\t' << itB->second << '\t' << distance << endl;
					out << itB->second << '\t' << itA->second << '\t' << distance << endl;
				}
				
				reading->update(itA->second * nseqs / 2);
			}
			gobble(fileHandle);
		}
		out.close();
		fileHandle.close();
	
		string squareFile;
		if(lt == 0){  // oops, it was square
			squareFile = filename;
		}else{ squareFile = tempOutFile; }
		
		//sort file by first column so the distances for each row are together
		string outfile = getRootName(squareFile) + "sorted.dist.temp";
		
		//use the unix sort 
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			string command = "sort -n " + squareFile + " -o " + outfile;
			system(command.c_str());
		#else //sort using windows sort
			string command = "sort " + squareFile + " /O " + outfile;
			system(command.c_str());
		#endif
		

		//output to new file distance for each row and save positions in file where new row begins
		ifstream in;
		openInputFile(outfile, in);
		
		distFile = outfile + ".rowFormatted";
		openOutputFile(distFile, out);
		
		rowPos.resize(nseqs, -1);
		int currentRow;
		int first, second;
		float dist;
		map<int, float> rowMap;
		map<int, float>::iterator itRow;
		
		//get first currentRow
		in >> first;
		currentRow = first;
		
		string firstString = toString(first);
		for(int k = 0; k < firstString.length(); k++)  {   in.putback(firstString[k]);  }
		
		while(!in.eof()) {
			in >> first >> second >> dist; gobble(in);
			
			if (first != currentRow) {
				//save position in file of each new row
				rowPos[currentRow] = out.tellp();
				
				out << currentRow << '\t' << rowMap.size() << '\t';
				
				for (itRow = rowMap.begin(); itRow != rowMap.end(); itRow++) {
					out << itRow->first << '\t' << itRow->second << '\t';
				}
				out << endl;
				
				currentRow = first;
				rowMap.clear();
				
				//save row you just read
				if (dist < cutoff) {
					rowMap[second] = dist;
				}
			}else{
				if (dist < cutoff) {
					rowMap[second] = dist;
				}
			}
		}
		
		//print last Row
		//save position in file of each new row
		rowPos[currentRow] = out.tellp();
		
		out << currentRow << '\t' << rowMap.size() << '\t';
		
		for (itRow = rowMap.begin(); itRow != rowMap.end(); itRow++) {
			out << itRow->first << '\t' << itRow->second << '\t';
		}
		out << endl;
		
		
		in.close();
		out.close();
		
		
		remove(tempOutFile.c_str());
		remove(outfile.c_str());
		
		reading->finish();
		list->setLabel("0");

	}
	catch(exception& e) {
		m->errorOut(e, "FormatColumnMatrix", "read");
		exit(1);
	}
}

/***********************************************************************/
FormatColumnMatrix::~FormatColumnMatrix(){}
/***********************************************************************/



