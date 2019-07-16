/*
 *  pds.seqdist.cpp
 *  
 *
 *  Created by Pat Schloss on 8/12/11.
 *  Copyright 2011 Patrick D. Schloss. All rights reserved.
 *
 */

#include "myseqdist.h"
#include "sequence.hpp"

/**************************************************************************************************/
correctDist::correctDist(int p) : processors(p) {
	try {
		m = MothurOut::getInstance();
	}
	catch(exception& e) {
		m->errorOut(e, "correctDist", "correctDist");
		exit(1);
	}
}
/**************************************************************************************************/
correctDist::correctDist(string sequenceFileName, int p) : processors(p) {
	try {
		m = MothurOut::getInstance();
		getSequences(sequenceFileName);
	}
	catch(exception& e) {
		m->errorOut(e, "correctDist", "correctDist");
		exit(1);
	}
}
/**************************************************************************************************/
int correctDist::addSeq(string seqName, string seqSeq){
	try {
		names.push_back(seqName);
		sequences.push_back(fixSequence(seqSeq));
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "correctDist", "addSeq");
		exit(1);
	}
}
/**************************************************************************************************/
void correctDist::execute(string distanceFileName){
	try {
        createProcess(distanceFileName);
	}
	catch(exception& e) {
		m->errorOut(e, "correctDist", "execute");
		exit(1);
	}
}
/**************************************************************************************************/
int correctDist::getSequences(string sequenceFileName){
	try {
		ifstream sequenceFile;
        Utils util; util.openInputFile(sequenceFileName, sequenceFile);
		string seqName, seqSeq;
		
		while(!sequenceFile.eof()){
			if (m->getControl_pressed()) { break; }
			
			Sequence temp(sequenceFile); util.gobble(sequenceFile);
			
			if (temp.getName() != "") {
				names.push_back(temp.getName());
				sequences.push_back(fixSequence(temp.getAligned()));
			}
		}
		sequenceFile.close();
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "correctDist", "getSequences");
		exit(1);
	}	
}

/**************************************************************************************************/
vector<int> correctDist::fixSequence(string sequence){
	try {
		int alignLength = sequence.length();
		
		vector<int> seqVector;
		
		for(int i=0;i<alignLength;i++){
			if(sequence[i] == 'A')		{	seqVector.push_back(0);		}
			else if(sequence[i] == 'C')	{	seqVector.push_back(1);		}
			else if(sequence[i] == 'T')	{	seqVector.push_back(2);		}
			else if(sequence[i] == 'G')	{	seqVector.push_back(3);		}
			else if(sequence[i] == 'N')	{	seqVector.push_back(0);		}//hmmmm....
		}
		
		return seqVector;
	}
	catch(exception& e) {
		m->errorOut(e, "correctDist", "fixSequence");
		exit(1);
	}	
}
/**************************************************************************************************/
struct correctData {
    string outputFileName;
    long long startLine, endLine;
    vector<vector<double> > correctMatrix;
    vector<vector<int> > sequences;
    MothurOut* m;
    Utils util;
    
    correctData(string ofn, vector<vector<int> > seqs, long long sLine, long long eLine) {
        outputFileName = ofn;
        startLine = sLine;
        endLine = eLine;
        sequences = seqs;
        m = MothurOut::getInstance();
        
        correctMatrix.resize(4);
        for(int i=0;i<4;i++){	correctMatrix[i].resize(4);	}
        
        correctMatrix[0][0] = 0.000000;		//AA
        correctMatrix[1][0] = 11.619259;	//CA
        correctMatrix[2][0] = 11.694004;	//TA
        correctMatrix[3][0] = 7.748623;		//GA
        
        correctMatrix[1][1] = 0.000000;		//CC
        correctMatrix[2][1] = 7.619657;		//TC
        correctMatrix[3][1] = 12.852562;	//GC
        
        correctMatrix[2][2] = 0.000000;		//TT
        correctMatrix[3][2] = 10.964048;	//TG
        
        correctMatrix[3][3] = 0.000000;		//GG
        
        for(int i=0;i<4;i++){
            for(int j=0;j<i;j++){
                correctMatrix[j][i] = correctMatrix[i][j];
            }
        }
    }
};
/**************************************************************************************************/
int getLastMatch(char direction, vector<vector<char> >& alignMoves, int i, int j, vector<int>& seqA, vector<int>& seqB, MothurOut* m){
    try {
        
        char nullReturn = -1;
        
        while(i>=1 && j>=1){
            
            if (m->getControl_pressed()) { return nullReturn; }
            
            if(direction == 'd'){
                if(seqA[i-1] == seqB[j-1])	{	return seqA[i-1];	}
                else						{	return nullReturn;	}
            }
            
            else if(direction == 'l')		{	j--;				}
            else							{	i--;				}
            
            direction = alignMoves[i][j];
        }
        
        return nullReturn;
    }
    catch(exception& e) {
        m->errorOut(e, "correctDist", "getLastMatch");
        exit(1);
    }	
}
/**************************************************************************************************/
double getDist(vector<int>& seqA, vector<int>& seqB, vector<vector<double> >& correctMatrix, MothurOut* m){
    try {
        
        int lengthA = seqA.size();
        int lengthB = seqB.size();
        
        vector<vector<double> > alignMatrix(lengthA+1);
        vector<vector<char> > alignMoves(lengthA+1);
        
        for(int i=0;i<=lengthA;i++){
            alignMatrix[i].resize(lengthB+1, 0);
            alignMoves[i].resize(lengthB+1, 'x');
        }
        
        for(int i=0;i<=lengthA;i++){
            alignMatrix[i][0] = 15.0 * i;
            alignMoves[i][0] = 'u';
        }
        for(int i=0;i<=lengthB;i++){
            alignMatrix[0][i] = 15.0 * i;
            alignMoves[0][i] = 'l';
        }
        
        for(int i=1;i<=lengthA;i++){
            for(int j=1;j<=lengthB;j++){
                
                if (m->getControl_pressed()) {  return 0;  }
                
                double nogap;
                nogap = alignMatrix[i-1][j-1] + correctMatrix[seqA[i-1]][seqB[j-1]];

                double gap;
                double left;
                if(i == lengthA){ left = alignMatrix[i][j-1]; } //terminal gap
                else{
                    if(seqB[j-1] == getLastMatch('l', alignMoves, i, j, seqA, seqB, m))     { gap = 4.0;    }
                    else                                                                    {  gap = 15.0;  }
                    
                    left = alignMatrix[i][j-1] + gap;
                }
                
                
                double up;
                if(j == lengthB){ up = alignMatrix[i-1][j]; }  //terminal gap
                else{
                    if(seqA[i-1] == getLastMatch('u', alignMoves, i, j, seqA, seqB, m))     {  gap = 4.0;   }
                    else                                                                    {  gap = 15.0;  }
                    
                    up = alignMatrix[i-1][j] + gap;
                }

                if(nogap < left){
                    if(nogap < up){
                        alignMoves[i][j] = 'd';
                        alignMatrix[i][j] = nogap;
                    }
                    else{
                        alignMoves[i][j] = 'u';
                        alignMatrix[i][j] = up;
                    }
                }
                else{
                    if(left < up){
                        alignMoves[i][j] = 'l';
                        alignMatrix[i][j] = left;
                    }
                    else{
                        alignMoves[i][j] = 'u';
                        alignMatrix[i][j] = up;
                    }
                }
            }
        }
        
        int i = lengthA;
        int j = lengthB;
        int count = 0;
        
        while(i > 0 && j > 0){
            
            if (m->getControl_pressed()) {  return 0;  }
            
            if(alignMoves[i][j] == 'd')         { count++; i--; j--;                    }
            else if(alignMoves[i][j] == 'u')    { if(j != lengthB){ count++; } i--;     }
            else if(alignMoves[i][j] == 'l')    {  if(i != lengthA){ count++; } j--;    }
        }
        
        return alignMatrix[lengthA][lengthB] / (double)count;
    }
    catch(exception& e) {
        m->errorOut(e, "correctDist", "getDist");
        exit(1);
    }	
}
/**************************************************************************************************/
int driverCorrect(correctData* params){
    try {
        ofstream distFile;
        params->util.openOutputFile(params->outputFileName, distFile);
        distFile << setprecision(9);
        
        if(params->startLine == 0){ distFile << params->sequences.size() << endl; }
        
        int startTime = time(NULL);
        params->m->mothurOut("\nCalculating distances for (" + toString(params->startLine+1) + " to " + toString(params->endLine+1) + ")... \n");
        
        for(int i = params->startLine;i < params->endLine; i++){
            if (params->m->getControl_pressed()) { distFile.close(); return 0; }
            
            distFile << i;
            for(int j=0;j<i;j++){ distFile << ' ' << getDist(params->sequences[i], params->sequences[j], params->correctMatrix, params->m); }
            distFile << endl;
            
            if(i % 100 == 0){ params->m->mothurOutJustToScreen(toString(i) + "\t" + toString(time(NULL) - startTime)+"\n"); }
        }
        distFile.close();
        
        if((params->endLine-1) % 100 != 0){ params->m->mothurOutJustToScreen(toString(params->endLine-1) + "\t" + toString(time(NULL) - startTime)+"\n"); }
        params->m->mothurOut("Done.\n");
        
        return 0;
    }
    catch(exception& e) {
        params->m->errorOut(e, "correctDist", "driverCorrect");
        exit(1);
    }	
}
/**************************************************************************************************/

int correctDist::createProcess(string distanceFileName){
	try {
        vector<linePair> lines;
        for(int i=0;i<processors;i++){
            linePair thisLine(int (sqrt(float(i)/float(processors)) * sequences.size()), int (sqrt(float(i+1)/float(processors)) * sequences.size()));
            lines.push_back(thisLine);
        }

        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<correctData*> data;
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            string extension = toString(i+1) + ".temp";
            correctData* dataBundle = new correctData(distanceFileName+extension, sequences, lines[i+1].start, lines[i+1].end);
            data.push_back(dataBundle);
            
            std::thread* thisThread = new std::thread(driverCorrect, dataBundle);
            workerThreads.push_back(thisThread);
        }
        
        correctData* dataBundle = new correctData(distanceFileName, sequences, lines[0].start, lines[0].end);
        driverCorrect(dataBundle);
        delete dataBundle;
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            
            string extension = toString(i+1) + ".temp";
            util.appendFiles((distanceFileName+extension), distanceFileName);
            util.mothurRemove(distanceFileName+extension);
            
            delete data[i];
            delete workerThreads[i];
        }

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "correctDist", "createProcess");
		exit(1);
	}	
}
/**************************************************************************************************/



