/*
 *  readcluster.cpp
 *  Mothur
 *
 *  Created by westcott on 10/28/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "readcluster.h"

/***********************************************************************/

ReadCluster::ReadCluster(string distfile, float c){
		globaldata = GlobalData::getInstance();
        distFile = distfile;
		cutoff = c;
}

/***********************************************************************/

void ReadCluster::read(NameAssignment* nameMap){
	try {
        
		if (format == "phylip") { convertPhylip2Column(nameMap); }
		else { list = new ListVector(nameMap->getListVector());  }
		
		OutPutFile = sortFile(distFile);
			
	}
	catch(exception& e) {
		errorOut(e, "ReadCluster", "read");
		exit(1);
	}
}
/***********************************************************************/

void ReadCluster::convertPhylip2Column(NameAssignment* nameMap){
	try {	
		//convert phylip file to column file
		map<int, string> rowToName;
		map<int, string>::iterator it;
		
		ifstream in;
		ofstream out;
		string tempFile = distFile + ".column.temp";
		
		openInputFile(distFile, in);
		openOutputFile(tempFile, out);
		
		float distance;
		int square, nseqs;
		string name;
		vector<string> matrixNames;
        
		in >> nseqs >> name;
		rowToName[0] = name;
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
		while((d=in.get()) != EOF){
			
			if(isalnum(d)){
				square = 1;
				in.putback(d);
				for(int i=0;i<nseqs;i++){
					in >> distance;
				}
				break;
			}
			if(d == '\n'){
				square = 0;
				break;
			}
		}
        
		if(square == 0){
					
			for(int i=1;i<nseqs;i++){
				in >> name;
				rowToName[i] = name;
				matrixNames.push_back(name);
				
				//there's A LOT of repeated code throughout this method...
				if(nameMap == NULL){
					list->set(i, name);
					
					for(int j=0;j<i;j++){
						in >> distance;
						
						if (distance == -1) { distance = 1000000; }
						
						if(distance < cutoff){
							out << i << '\t' << j << '\t' << distance << endl;
						}
					}
					
				}
				else{
					if(nameMap->count(name)==0){        mothurOut("Error: Sequence '" + name + "' was not found in the names file, please correct"); mothurOutEndLine(); }
					
					for(int j=0;j<i;j++){
						in >> distance;
						
						if (distance == -1) { distance = 1000000; }
						
						if(distance < cutoff){
							out << i << '\t' << j << '\t' << distance << endl;
						}
					}
				}
			}
		}
		else{
			for(int i=1;i<nseqs;i++){
				in >> name;                
				rowToName[i] = name;
				matrixNames.push_back(name);
		
				if(nameMap == NULL){
					list->set(i, name);
					for(int j=0;j<nseqs;j++){
						in >> distance;
					
						if (distance == -1) { distance = 1000000; }
						
						if(distance < cutoff && j < i){
							out << i << '\t' << j << '\t' << distance << endl;
						}
					}
				}
				else{
					if(nameMap->count(name)==0){        mothurOut("Error: Sequence '" + name + "' was not found in the names file, please correct"); mothurOutEndLine(); }
					
					for(int j=0;j<nseqs;j++){
						in >> distance;
                        
						if (distance == -1) { distance = 1000000; }
						
						if(distance < cutoff && j < i){
							out << i << '\t' << j << '\t' << distance << endl;
						}
						
					}
				}
			}
		}
		
		list->setLabel("0");
		in.close();
		out.close();
		
		if(nameMap == NULL){
			nameMap = new NameAssignment();
			for(int i=0;i<matrixNames.size();i++){
				nameMap->push_back(matrixNames[i]);
			}
			globaldata->nameMap = nameMap;
		}
		
	
		ifstream in2;
		ofstream out2;
		
		string outputFile = getRootName(distFile) + "column.dist";
		openInputFile(tempFile, in2);
		openOutputFile(outputFile, out2);
		
		int first, second;
		float dist;
		
		while (in2) {
			in2 >> first >> second >> dist;
			out2 << rowToName[first] << '\t' << rowToName[second] << '\t' << dist << endl;
			gobble(in2);
		}
		in2.close();
		out2.close();
		
		remove(tempFile.c_str());
		distFile = outputFile;
	}
	catch(exception& e) {
		errorOut(e, "ReadCluster", "convertPhylip2Column");
		exit(1);
	}
}
/***********************************************************************/

ReadCluster::~ReadCluster(){}
/***********************************************************************/

