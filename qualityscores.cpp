/*
 *  qualityscores.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 7/12/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "qualityscores.h"

/**************************************************************************************************/

QualityScores::QualityScores(){
	try {
		m = MothurOut::getInstance();
		seqName = "";
		seqLength = -1;
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "QualityScores");
		exit(1);
	}							
}

/**************************************************************************************************/

QualityScores::QualityScores(ifstream& qFile, int l){
	try {
		
		m = MothurOut::getInstance();

		seqName = "";
		seqLength = l;
		int score;
		
		string line;
		getline(qFile, line); gobble(qFile);
		istringstream nameStream(line);
	
		nameStream >> seqName;
		seqName = seqName.substr(1); 

		//getline(qFile, line);
		//istringstream qualStream(line);
	
		//while(qualStream){
		//	qualStream >> score;
		//	qScores.push_back(score);
		//}
		//qScores.pop_back();
		
		//seqLength = qScores.size();	
		
		for(int i=0;i<seqLength;i++){
			qFile >> score;
			qScores.push_back(score);
		}
		gobble(qFile);

	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "QualityScores");
		exit(1);
	}							

}

/**************************************************************************************************/

string QualityScores::getName(){
	
	try {
		return seqName;
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "getName");
		exit(1);
	}							
}

/**************************************************************************************************/

void QualityScores::printQScores(ofstream& qFile){
	try {
		
		double aveQScore = calculateAverage();
		
		qFile << '>' << seqName << '\t' << aveQScore << endl;
		
		for(int i=0;i<seqLength;i++){
			qFile << qScores[i] << ' ';
		}
		qFile << endl;
				
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "printQScores");
		exit(1);
	}							
}

/**************************************************************************************************/

void QualityScores::trimQScores(int start, int end){
	try {
		vector<int> hold;
		
		if(end == -1){		
			hold = vector<int>(qScores.begin()+start, qScores.end());
			qScores = hold;		
		}
		if(start == -1){
			hold = vector<int>(qScores.begin(), qScores.begin()+end);	//not sure if indexing is correct
			qScores = hold;		
		}

		seqLength = qScores.size();
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "trimQScores");
		exit(1);
	}							
}

/**************************************************************************************************/

void QualityScores::flipQScores(){
	try {
		
		vector<int> temp = qScores;
		for(int i=0;i<seqLength;i++){
			qScores[seqLength - i - 1] = temp[i];
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "flipQScores");
		exit(1);
	}							
}

/**************************************************************************************************/

bool QualityScores::stripQualThreshold(Sequence& sequence, double qThreshold){
	try {
		string rawSequence = sequence.getUnaligned();
		int seqLength = sequence.getNumBases();
		
		if(seqName != sequence.getName()){
			m->mothurOut("sequence name mismatch btwn fasta: " + sequence.getName() + " and qual file: " + seqName);
			m->mothurOutEndLine();	
		}
		
		int end;
		for(int i=0;i<seqLength;i++){
			end = i;
			if(qScores[i] < qThreshold){
				break;
			}
		}
		
		sequence.setUnaligned(rawSequence.substr(0,end));
		trimQScores(-1, end);
	
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "flipQScores");
		exit(1);
	}							
	
}

/**************************************************************************************************/

bool QualityScores::stripQualRollingAverage(Sequence& sequence, double qThreshold){
	try {
		string rawSequence = sequence.getUnaligned();
		int seqLength = sequence.getNumBases();
		
		if(seqName != sequence.getName()){
			m->mothurOut("sequence name mismatch btwn fasta: " + sequence.getName() + " and qual file: " + seqName);
			m->mothurOutEndLine();	
		}
		
		int end = -1;
		double rollingSum = 0.0000;
		
		for(int i=0;i<seqLength;i++){

			rollingSum += (double)qScores[i];
			
			if(rollingSum / (double)(i+1) < qThreshold){
				end = i;
//				cout << i+1 << '\t' << seqName << '\t' << rollingSum / (double)(i+1) << endl;
				break;
			}
		}
		
		if(end == -1){	end = seqLength;	}
		
		sequence.setUnaligned(rawSequence.substr(0,end));
		trimQScores(-1, end);
		
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "flipQScores");
		exit(1);
	}							
	
}

/**************************************************************************************************/

bool QualityScores::stripQualWindowAverage(Sequence& sequence, int stepSize, int windowSize, double qThreshold){
	try {
		string rawSequence = sequence.getUnaligned();
		int seqLength = sequence.getNumBases();
		
		if(seqName != sequence.getName()){
			m->mothurOut("sequence name mismatch btwn fasta: " + sequence.getName() + " and qual file: " + seqName);
			m->mothurOutEndLine();	
		}
		
		int end = windowSize;
		int start = 0;
		

		while(start < seqLength){
			double windowSum = 0.0000;

			for(int i=start;i<end;i++){
				windowSum += qScores[i];				
			}
			double windowAverage = windowSum / (double)(end-start);
			
			if(windowAverage < qThreshold){
				end = start;
				break;
			}
			start += stepSize;
			end = start + windowSize;
			if(end >= seqLength){	end = seqLength - 1;	}
		}
		
		
		if(end == -1){	end = seqLength;	}
		
		sequence.setUnaligned(rawSequence.substr(0,end));
		trimQScores(-1, end);
		
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "flipQScores");
		exit(1);
	}							
	
}

/**************************************************************************************************/

double QualityScores::calculateAverage(){
	
	double aveQScore = 0.0000;
	
	for(int i=0;i<seqLength;i++){
		aveQScore += (double) qScores[i];
	}
	aveQScore /= (double) seqLength;
	
	return aveQScore;
}

/**************************************************************************************************/

bool QualityScores::cullQualAverage(Sequence& sequence, double qAverage){
	try {
		string rawSequence = sequence.getUnaligned();
		bool success = 0;	//guilty until proven innocent
		
		if(seqName != sequence.getName())	{
			m->mothurOut("sequence name mismatch btwn fasta: " + sequence.getName() + " and qual file: " + seqName);
			m->mothurOutEndLine();	
		} 
			
		double aveQScore = calculateAverage();
		
		if(aveQScore >= qAverage)	{	success = 1;	}
		else						{	success = 0;	}
		
		return success;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "cullQualAverage");
		exit(1);
	}
}

/**************************************************************************************************/
