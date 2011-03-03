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

QualityScores::QualityScores(ifstream& qFile){
	try {
		
		m = MothurOut::getInstance();
		
		seqName = "";
		int score;
		
		qFile >> seqName; 
		m->getline(qFile);
		
		if (seqName == "")	{
			m->mothurOut("Error reading quality file, name blank at position, " + toString(qFile.tellg()));
			m->mothurOutEndLine(); 
		}
		else{
			seqName = seqName.substr(1);
		}
		
		string qScoreString = m->getline(qFile);
		
		istringstream qScoreStringStream(qScoreString);
		while(!qScoreStringStream.eof()){
			qScoreStringStream >> score;
			qScores.push_back(score);
		}
		qScores.pop_back();
		seqLength = qScores.size();
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
			if(qScores.size() > end){
				hold = vector<int>(qScores.begin(), qScores.begin()+end);
				qScores = hold;		
			}
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
		
		//every score passed
		if (end == (seqLength-1)) { end = seqLength; }
		
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
			m->mothurOut("sequence name mismatch between fasta: " + sequence.getName() + " and qual file: " + seqName);
			m->mothurOutEndLine();
		}
		
		int end = windowSize;
		int start = 0;
		
		if(seqLength < windowSize) {	return 0;	}
			
		while(start < seqLength){
			double windowSum = 0.0000;

			for(int i=start;i<end;i++){
				windowSum += qScores[i];
			}
			double windowAverage = windowSum / (double)(end-start);
			
			if(windowAverage < qThreshold){
				end = end - stepSize;
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
		m->errorOut(e, "QualityScores", "cullQualAverage");
		exit(1);
	}
}

/**************************************************************************************************/

void QualityScores::updateQScoreErrorMap(map<char, vector<int> >& qualErrorMap, string errorSeq, int start, int stop, int weight){
	try {

		int seqLength = errorSeq.size();
		
		int qIndex = start - 1;
		for(int i=0;i<seqLength;i++){
			
			if(errorSeq[i] == 'm')		{	qualErrorMap['m'][qScores[qIndex]] += weight;	}
			else if(errorSeq[i] == 's')	{	qualErrorMap['s'][qScores[qIndex]] += weight;	}
			else if(errorSeq[i] == 'i')	{	qualErrorMap['i'][qScores[qIndex]] += weight;	}
			else if(errorSeq[i] == 'a')	{	qualErrorMap['a'][qScores[qIndex]] += weight;	}
			else if(errorSeq[i] == 'd')	{	/*	there are no qScores for deletions	*/		}

			if(errorSeq[i] != 'd')		{	qIndex++;	}
			if(qIndex > stop){	break;	}
		}	
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "updateQScoreErrorMap");
		exit(1);
	}
}

/**************************************************************************************************/

void QualityScores::updateForwardMap(vector<vector<int> >& forwardMap, int start, int stop, int weight){
	try {
		
		int index = 0;
		for(int i=start-1;i<stop;i++){
			forwardMap[index++][qScores[i]] += weight;
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "updateForwardMap");
		exit(1);
	}
}

/**************************************************************************************************/

void QualityScores::updateReverseMap(vector<vector<int> >& reverseMap, int start, int stop, int weight){
	try {
		
		int index = 0;
		for(int i=stop-1;i>=start;i--){
			reverseMap[index++][qScores[i]] += weight;
		}
		
	}	
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "updateForwardMap");
		exit(1);
	}
}

/**************************************************************************************************/
