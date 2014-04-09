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

QualityScores::QualityScores(string n, vector<int> s){
	try {
		m = MothurOut::getInstance();
		setName(n);
        setScores(s);
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

		int score;
		seqName = getSequenceName(qFile); m->gobble(qFile);
		
        if (m->debug) { m->mothurOut("[DEBUG]: name = '" + seqName + "'\n.");  }
        
		if (!m->control_pressed) {
            string qScoreString = m->getline(qFile); m->gobble(qFile);
            
            if (m->debug) { m->mothurOut("[DEBUG]: scores = '" + qScoreString + "'\n.");  }
            
            while(qFile.peek() != '>' && qFile.peek() != EOF){
                if (m->control_pressed) { break; }
                string temp = m->getline(qFile); m->gobble(qFile);
                //if (m->debug) { m->mothurOut("[DEBUG]: scores = '" + temp + "'\n.");  }
                qScoreString +=  ' ' + temp;
            }
            //cout << "done reading " << endl;
            istringstream qScoreStringStream(qScoreString);
            int count = 0;
            while(!qScoreStringStream.eof()){
                if (m->control_pressed) { break; }
                string temp;
                qScoreStringStream >> temp;  m->gobble(qScoreStringStream);
                
                //if (m->debug) { m->mothurOut("[DEBUG]: score " + toString(qScores.size()) + " = '" + temp + "'\n.");  }
                
                //check temp to make sure its a number
                if (!m->isContainingOnlyDigits(temp)) { m->mothurOut("[ERROR]: In sequence " + seqName + "'s quality scores, expected a number and got " + temp + ", setting score to 0."); m->mothurOutEndLine(); temp = "0"; }
                convert(temp, score);
                
                //cout << count << '\t' << score << endl;
                qScores.push_back(score);
                count++;
            }
        }
		
		seqLength = qScores.size();
		//cout << "seqlength = " << seqLength  << endl;
		
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "QualityScores");
		exit(1);
	}							
	
}
//********************************************************************************************************************
string QualityScores::getSequenceName(ifstream& qFile) {
	try {
		string name = "";
		
        qFile >> name;
        m->getline(qFile);
		
		if (name.length() != 0) { 
            
			name = name.substr(1); 
            
            m->checkName(name);
            
        }else{ m->mothurOut("Error in reading your qfile, at position " + toString(qFile.tellg()) + ". Blank name."); m->mothurOutEndLine(); m->control_pressed = true;  }
        
		return name;
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "getSequenceName");
		exit(1);
	}
}
//********************************************************************************************************************
void QualityScores::setName(string name) {
	try {
      
        m->checkName(name);   
        seqName = name;
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "setName");
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
		
		double aveQScore = calculateAverage(false);
		
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
		

		//cout << seqName << '\t' << start << '\t' << end << '\t' << qScores.size() << endl;
		//for (int i = 0; i < qScores.size(); i++) { cout << qScores[i] << end; }
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
			m->mothurOutEndLine();	m->control_pressed = true;
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

bool QualityScores::stripQualRollingAverage(Sequence& sequence, double qThreshold, bool logTransform){
	try {
		string rawSequence = sequence.getUnaligned();
		int seqLength = sequence.getNumBases();
		
		if(seqName != sequence.getName()){
			m->mothurOut("sequence name mismatch btwn fasta: " + sequence.getName() + " and qual file: " + seqName);
			m->mothurOutEndLine();	
		}
		
		int end = -1;
		double rollingSum = 0.0000;
        double value = 0.0;
		
		for(int i=0;i<seqLength;i++){
            
            if (logTransform)   {
                rollingSum += (double)pow(10.0, qScores[i]);
                value = log10(rollingSum / (double)(i+1));
                
            } //Sum 10^Q
            else                {
                rollingSum += (double)qScores[i];
                value = rollingSum / (double)(i+1);
            }
			
			
			if(value < qThreshold){
				end = i;
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

bool QualityScores::stripQualWindowAverage(Sequence& sequence, int stepSize, int windowSize, double qThreshold, bool logTransform){
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
			
		while((start+windowSize) < seqLength){
			double windowSum = 0.0000;

			for(int i=start;i<end;i++){
                if (logTransform)   {  windowSum += pow(10.0, qScores[i]);  }
                else                {  windowSum += qScores[i];             }
			}
			double windowAverage = 0.0;
            if (logTransform)   { windowAverage = log10(windowSum / (double)(end-start)); }
            else                { windowAverage = windowSum / (double)(end-start);      }
				
			if(windowAverage < qThreshold){
				end = end - stepSize;
				break;
			}
			
			start += stepSize;
			end = start + windowSize;
				
			if(end >= seqLength){	end = seqLength;	}
				
		}
	
		if(end == -1){	end = seqLength;	}
		
		//failed first window
		if (end < windowSize) { return 0; }
			
		sequence.setUnaligned(rawSequence.substr(0,end));
		trimQScores(-1, end);
		
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "stripQualWindowAverage");
		exit(1);
	}							
	
}

/**************************************************************************************************/

double QualityScores::calculateAverage(bool logTransform){
	
	double aveQScore = 0.0000;
	
	for(int i=0;i<seqLength;i++){
        if (logTransform)   {  aveQScore += pow(10.0, qScores[i]);  }
        else                {  aveQScore += qScores[i];             }
	}
    
    if (logTransform)   {  aveQScore = log10(aveQScore /(double) seqLength);    }
    else                {  aveQScore /= (double) seqLength;                     }
	
	return aveQScore;
}

/**************************************************************************************************/

bool QualityScores::cullQualAverage(Sequence& sequence, double qAverage, bool logTransform){
	try {
		string rawSequence = sequence.getUnaligned();
		bool success = 0;	//guilty until proven innocent
		
		if(seqName != sequence.getName())	{
			m->mothurOut("sequence name mismatch btwn fasta: " + sequence.getName() + " and qual file: " + seqName);
			m->mothurOutEndLine();	
		} 
			
		double aveQScore = calculateAverage(logTransform);
		
        if (m->debug) { m->mothurOut("[DEBUG]: " + sequence.getName() + " average = " + toString(aveQScore) + "\n"); }
        
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
			else if(errorSeq[i] == 'a')	{	qualErrorMap['a'][qScores[qIndex]] += weight;	/*if(qScores[qIndex] != 0){	cout << qIndex << '\t';		}*/	}
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
		for(int i=stop-1;i>=start-1;i--){
			reverseMap[index++][qScores[i]] += weight;
		}
		
	}	
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "updateReverseMap");
		exit(1);
	}
}

/**************************************************************************************************/
