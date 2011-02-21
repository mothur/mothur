/*
 *  readphylipvector.cpp
 *  mothur
 *
 *  Created by westcott on 1/11/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "readphylipvector.h"

/***********************************************************************/
ReadPhylipVector::ReadPhylipVector(string d) {
	try {
		m = MothurOut::getInstance();
		distFile = d;
	}
	catch(exception& e) {
		m->errorOut(e, "ReadPhylipVector", "ReadPhylipVector");
		exit(1);
	}
}
/***********************************************************************/
vector<string> ReadPhylipVector::read(vector< vector<double> >& matrix) {
	try {
		vector<string> names;
		
		ifstream in;
		m->openInputFile(distFile, in);
		
		//check whether matrix is square
		char d;
		int square = 1;
		int numSeqs;
		string name;
		
		in >> numSeqs >> name; 
		
		while((d=in.get()) != EOF){
			
			//is d a number meaning its square
			if(isalnum(d)){ 
				square = 1; 
				break; 
			}
			
			//is d a line return meaning its lower triangle
			if(d == '\n'){
				square = 2;
				break;
			}
		}
		in.close();
		
		
		//reopen and read now that you know whether you are square
		ifstream f;
		m->openInputFile(distFile, f);
		
		int rank;
		f >> rank;
		
		names.resize(rank);
		matrix.resize(rank);
		if(square == 1){
			for(int i=0;i<rank;i++)
				matrix[i].resize(rank);
			for(int i=0;i<rank;i++) {
				f >> names[i];
				for(int j=0;j<rank;j++) {
					if (m->control_pressed) { return names; }
					
					f >> matrix[i][j];
					if (matrix[i][j] == -0.0000)
						matrix[i][j] = 0.0000;
				}
			}
		}
		else if(square == 2){
			for(int i=0;i<rank;i++){
				matrix[i].resize(rank);
			}
			matrix[0][0] = 0.0000;
			f >> names[0];
			for(int i=1;i<rank;i++){
				f >> names[i];
				matrix[i][i]=0.0000;
				for(int j=0;j<i;j++){
					if (m->control_pressed) { return names; }
					f >> matrix[i][j];
					if (matrix[i][j] == -0.0000)
						matrix[i][j] = 0.0000;
					matrix[j][i]=matrix[i][j];
				}
			}
		}
		f.close();
		
		return names;
	}
	catch(exception& e) {
		m->errorOut(e, "ReadPhylipVector", "read");
		exit(1);
	}
}
/***********************************************************************/
vector<string> ReadPhylipVector::read(vector<seqDist>& matrix) {
	try {
		vector<string> names;
		
		ifstream in;
		m->openInputFile(distFile, in);
		
		//check whether matrix is square
		char d;
		int square = 1;
		int numSeqs;
		string name;
		
		in >> numSeqs >> name; 
		
		while((d=in.get()) != EOF){
			
			//is d a number meaning its square
			if(isalnum(d)){ 
				square = 1; 
				break; 
			}
			
			//is d a line return meaning its lower triangle
			if(d == '\n'){
				square = 2;
				break;
			}
		}
		in.close();
		
		
		//reopen and read now that you know whether you are square
		ifstream f;
		m->openInputFile(distFile, f);
		
		int rank;
		float temp;
		f >> rank;
		
		names.resize(rank);
		if(square == 1){
			for(int i=0;i<rank;i++) {
				f >> names[i];
				for(int j=0;j<rank;j++) {
					if (m->control_pressed) { return names; }
					
					f >> temp;
					
					if (j < i) { //only save lt
						seqDist dist(i, j, temp);
						matrix.push_back(dist);
					}
				}
			}
		}
		else if(square == 2){
			f >> names[0];
			for(int i=1;i<rank;i++){
				f >> names[i];
				for(int j=0;j<i;j++){
					if (m->control_pressed) { return names; }
					f >> temp;
					seqDist dist(i, j, temp);
					matrix.push_back(dist);
				}
			}
		}
		f.close();
		
		return names;
	}
	catch(exception& e) {
		m->errorOut(e, "ReadPhylipVector", "read");
		exit(1);
	}
}
/***********************************************************************/


