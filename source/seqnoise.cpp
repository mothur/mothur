/*
 *  mySeqNoise.cpp
 *  
 *
 *  Created by Pat Schloss on 8/31/11.
 *  Copyright 2011 Patrick D. Schloss. All rights reserved.
 *
 */

#include "seqnoise.h"
#include "sequence.hpp"
#include "listvector.hpp"
#include "inputdata.h"

#define MIN_DELTA 1.0e-6
#define MIN_ITER 20
#define MAX_ITER 1000
#define MIN_COUNT 0.1
#define MIN_TAU   1.0e-4
#define MIN_WEIGHT 0.1


/**************************************************************************************************/
int seqNoise::getSequenceData(string sequenceFileName, vector<string>& sequences){
	try {
		
		ifstream sequenceFile;
		util.openInputFile(sequenceFileName, sequenceFile);
		
		while(!sequenceFile.eof()){
			
			if (m->getControl_pressed()) { break; }
			
			Sequence temp(sequenceFile); util.gobble(sequenceFile);
			
			if (temp.getName() != "") {
				sequences.push_back(temp.getAligned());
			}
		}
		sequenceFile.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "seqNoise", "getSequenceData");
		exit(1);
	}
}
/**************************************************************************************************/
int seqNoise::addSeq(string seq, vector<string>& sequences){ 
	try {
		sequences.push_back(seq);
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "seqNoise", "addSeq");
		exit(1);
	}
}	
/**************************************************************************************************/
//no checks for file mismatches
int seqNoise::getRedundantNames(string namesFileName, vector<string>& uniqueNames, vector<string>& redundantNames, vector<int>& seqFreq){
	try {
		string unique, redundant;
		ifstream namesFile;
		util.openInputFile(namesFileName, namesFile);
		
		for(int i=0;i<redundantNames.size();i++){
			
			if (m->getControl_pressed()) { break; }
			
			namesFile >> uniqueNames[i]; util.gobble(namesFile);
			namesFile >> redundantNames[i]; util.gobble(namesFile);
			
			seqFreq[i] = util.getNumNames(redundantNames[i]);
		}
		namesFile.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "seqNoise", "getRedundantNames");
		exit(1);
	}
}
/**************************************************************************************************/
int seqNoise::addRedundantName(string uniqueName, string redundantName, vector<string>& uniqueNames, vector<string>& redundantNames, vector<int>& seqFreq){ 
	try {
		
		uniqueNames.push_back(uniqueName);
		redundantNames.push_back(redundantName);
		seqFreq.push_back(util.getNumNames(redundantName));
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "seqNoise", "addRedundantName");
		exit(1);
	}
}	
/**************************************************************************************************/
int seqNoise::getDistanceData(string distFileName, vector<double>& distances){
	try {
		
		ifstream distFile;
		util.openInputFile(distFileName, distFile);
		
		int numSeqs = 0;
		string name = "";
		
		distFile >> numSeqs;
		
		for(int i=0;i<numSeqs;i++){
			
			if (m->getControl_pressed()) {  break; }
			
			distances[i * numSeqs + i] = 0.0000;
			
			distFile >> name;
			
			for(int j=0;j<i;j++){
				distFile >> distances[i * numSeqs + j];
				distances[j * numSeqs + i] = distances[i * numSeqs + j];
			}
		}
		
		distFile.close();	
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "seqNoise", "getDistanceData");
		exit(1);
	}
}

/**************************************************************************************************/
int seqNoise::getListData(string listFileName, double cutOff, vector<int>& otuData, vector<int>& otuFreq, vector<vector<int> >& otuBySeqLookUp){
	try {
		
		ifstream listFile;
		util.openInputFile(listFileName, listFile);
		
		bool adjustCutoff = true;
        string lastLabel = "";
        string readLabels = "";
        
		while(!listFile.eof()){
            
            ListVector list(listFile, readLabels, lastLabel); util.gobble(listFile); //10/18/13 - change to reading with listvector to accomodate changes to the listfiel format. ie. adding header labels.
            
            string thisLabel = list.getLabel();
            lastLabel = thisLabel;
            
            if (thisLabel == "unique") {} //skip to next label in listfile
            else {
                double threshold;
                util.mothurConvert(thisLabel, threshold);
                
                if(threshold < cutOff){} //skip to next label in listfile
                else{
                    adjustCutoff = false;
                    int numOTUs = list.getNumBins();
                    otuFreq.resize(numOTUs, 0);
                    
                    for(int i=0;i<numOTUs;i++){
                        
                        if (m->getControl_pressed()) { return 0; }
                        
                        string otu = list.get(i);
                        int count = 0;
                        string number = "";
                        
                        for(int j=0;j<otu.size();j++){
                            if(otu[j] != ','){
                                number += otu[j];
                            }
                            else{
                                int index = atoi(number.c_str());
                                otuData[index] = i;
                                count++;
                                number = "";
                            }
                        }
                        
                        int index = atoi(number.c_str());
                        otuData[index] = i;
                        count++;
                        
                        otuFreq[i] = count;
                    }
                    
                    otuBySeqLookUp.resize(numOTUs);
                    
                    int numSeqs = otuData.size();
                    
                    for(int i=0;i<numSeqs;i++){
                        if (m->getControl_pressed()) { return 0; }
                        otuBySeqLookUp[otuData[i]].push_back(i);
                    }
                    for(int i=0;i<numOTUs;i++){
                        if (m->getControl_pressed()) { return 0; }
                        for(int j=otuBySeqLookUp[i].size();j<numSeqs;j++){
                            otuBySeqLookUp[i].push_back(0);
                        }
                    }
                    
                    break;
                }
            }
		}
		
		listFile.close();
		
		//the listfile does not contain a threshold greater than the cutoff so use highest value
		if (adjustCutoff) {
            
            InputData input(listFileName, "list", nullVector);
            ListVector* list = input.getListVector(lastLabel);
            
            int numOTUs = list->getNumBins();
            otuFreq.resize(numOTUs, 0);
			
			for(int i=0;i<numOTUs;i++){
				
				if (m->getControl_pressed()) { return 0; }
				
				string otu = list->get(i);
				
				int count = 0;
				string number = "";
				
				for(int j=0;j<otu.size();j++){
					if(otu[j] != ','){
						number += otu[j];
					}
					else{
						int index = atoi(number.c_str());
						otuData[index] = i;
						count++;
						number = "";
					}
				}
				
				int index = atoi(number.c_str());
				otuData[index] = i;
				count++;
				
				otuFreq[i] = count;
			}
			
			otuBySeqLookUp.resize(numOTUs);
			
			int numSeqs = otuData.size();
			
			for(int i=0;i<numSeqs;i++){
				if (m->getControl_pressed()) { return 0; }
				otuBySeqLookUp[otuData[i]].push_back(i);
			}
			for(int i=0;i<numOTUs;i++){
				if (m->getControl_pressed()) { return 0; }
				for(int j=otuBySeqLookUp[i].size();j<numSeqs;j++){
					otuBySeqLookUp[i].push_back(0);
				}
			}
			
            delete list;
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "seqNoise", "getListData");
		exit(1);
	}
}

/**************************************************************************************************/
int seqNoise::updateOTUCountData(vector<int> otuFreq,
								 vector<vector<int> > otuBySeqLookUp,
								 vector<vector<int> > aanI,
								 vector<int>& anP,
								 vector<int>& anI,
								 vector<int>& cumCount
								 ){
	try {
		int numOTUs = otuFreq.size();
		
		int count = 0;
		
		for(int i=0;i<numOTUs;i++){
			cumCount[i] = count;
			
			if (m->getControl_pressed()) { return 0; }
			
			for(int j=0;j<otuFreq[i];j++){
				anP[count] = otuBySeqLookUp[i][j];
				anI[count] = aanI[i][j];
				
				count++;
			}
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "seqNoise", "updateOTUCountData");
		exit(1);
	}
}
/**************************************************************************************************/
double seqNoise::calcNewWeights(
					  vector<double>& weights,	//
					  vector<int> seqFreq,		//
					  vector<int> anI,			//
					  vector<int> cumCount,		//
					  vector<int> anP,			//
					  vector<int> otuFreq,		//
					  vector<double> tau		//
					  ){
	try {
		
		int numOTUs = weights.size();
		double maxChange = -1;
		
		cout.flush();
		
		for(int i=0;i<numOTUs;i++){
			
			if (m->getControl_pressed()) { return 0; }
			
			double change = weights[i];
			
			weights[i] = 0.0000;
			
			for(int j=0;j<otuFreq[i];j++){
				
				int index1 = cumCount[i] + j;
				int index2 = anI[index1];
				
				double currentTau = tau[anP[index1]];
				double freq = double(seqFreq[index2]);
				
				weights[i] += currentTau * freq;
			}
			change = fabs(weights[i] - change);
			
			if(change > maxChange){	maxChange = change;	}
			cout.flush();
		}
		return maxChange;
		
	}
	catch(exception& e) {
		m->errorOut(e, "seqNoise", "calcNewWeights");
		exit(1);
	}
}

/**************************************************************************************************/

int seqNoise::calcCentroids(
				   vector<int> anI,
				   vector<int> anP,
				   vector<int>& change, 
				   vector<int>& centroids, 
				   vector<int> cumCount,
				   vector<double> distances,///
				   vector<int> seqFreq, 
				   vector<int> otuFreq, 
				   vector<double> tau 
				   ){
	try {
		int numOTUs = change.size();
		int numSeqs = seqFreq.size();
		
		for(int i=0;i<numOTUs;i++){
			
			if (m->getControl_pressed()) { return 0; }
			
			int minFIndex = -1;
			double minFValue = 1e10;
			
			change[i] = 0;
			double count = 0.00000;
			
			int freqOfOTU = otuFreq[i];
			
			for(int j=0;j<freqOfOTU;j++){
				int index = cumCount[i] + j;
				count += seqFreq[anI[index]]*tau[anP[index]];
			}
			
			if(freqOfOTU > 0 && count > MIN_COUNT){
				
				vector<double> adF(freqOfOTU);
				vector<int> anL(freqOfOTU);
				
				for(int j=0;j<freqOfOTU;j++){
					anL[j] = anI[cumCount[i] + j];
					adF[j] = 0.0000;
				}
				
				for(int j=0;j<freqOfOTU;j++){		
					int index = cumCount[i] + j;
					double curTau = tau[anP[index]];
					
					for(int k=0;k<freqOfOTU;k++){
						double dist = distances[anL[j]*numSeqs + anL[k]];
						
						adF[k] += dist * curTau * seqFreq[anL[j]];
					}
				}
				
				for(int j=0;j<freqOfOTU;j++){
					if(adF[j] < minFValue){
						minFIndex = j;
						minFValue = adF[j];
					}
				}
				
				if(centroids[i] != anL[minFIndex]){
					change[i] = 1;
					centroids[i] = anL[minFIndex];
				}
			}
			else if(centroids[i] != -1){
				change[i] = 1;
				centroids[i] = -1;			
			}
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "seqNoise", "calcCentroids");
		exit(1);
	}
}

/**************************************************************************************************/

int seqNoise::checkCentroids(vector<double>& weights, vector<int> centroids){
	try {
		int numOTUs = centroids.size();
		vector<int> unique(numOTUs, 1);
		
		double minWeight = MIN_WEIGHT;
		for(int i=0;i<numOTUs;i++){
			if (m->getControl_pressed()) { return 0; }
			if(weights[i] < minWeight){	unique[i] = -1;	}
		}
		
		for(int i=0;i<numOTUs;i++){
			if (m->getControl_pressed()) { return 0; }
			if(unique[i] == 1){
				for(int j=i+1; j<numOTUs;j++){
					if(unique[j] == 1){
						if(centroids[i] == centroids[j]){
							unique[j] = 0;
							weights[i] += weights[j];
							weights[j] = 0;
						}
					}
				}
			}		
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "seqNoise", "checkCentroids");
		exit(1);
	}
}

/**************************************************************************************************/

int seqNoise::setUpOTUData(vector<int>& otuData, vector<double>& percentage, vector<int> cumCount, vector<double> tau, vector<int> otuFreq, vector<int> anP, vector<int> anI){
	try {

		int numOTUs = cumCount.size();
		int numSeqs = otuData.size();
		
		vector<double> bestTau(numSeqs, 0);
		vector<double> bestIndex(numSeqs, -1);
		
		for(int i=0;i<numOTUs;i++){
			if (m->getControl_pressed()) { return 0; }
			for(int j=0;j<otuFreq[i];j++){
				
				int index1 = cumCount[i] + j;
				double thisTau = tau[anP[index1]];
				int index2 = anI[index1];
				
				if(thisTau > bestTau[index2]){
					bestTau[index2] = thisTau;
					bestIndex[index2] = i;
				}
			}		
		}
		
		for(int i=0;i<numSeqs;i++){
			if (m->getControl_pressed()) { return 0; }
			otuData[i] = bestIndex[i];
			percentage[i] = 1 - bestTau[i];
		}
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "seqNoise", "setUpOTUData");
		exit(1);
	}
}

/**************************************************************************************************/

int seqNoise::finishOTUData(vector<int> otuData, vector<int>& otuFreq, vector<int>& anP, vector<int>& anI, vector<int>& cumCount, vector<vector<int> >& otuBySeqLookUp, vector<vector<int> >& aanI, vector<double>& tau){
	try {
		int numSeqs = otuData.size();
		int numOTUs = otuFreq.size();
		int total = numSeqs;
		
		otuFreq.assign(numOTUs, 0);
		tau.assign(numSeqs, 1);
		anP.assign(numSeqs, 0);
		anI.assign(numSeqs, 0);
		
		for(int i=0;i<numSeqs;i++){
			if (m->getControl_pressed()) { return 0; }
			int otu = otuData[i];
			total++;
			
			otuBySeqLookUp[otu][otuFreq[otu]] = i;
			aanI[otu][otuFreq[otu]] = i;
			otuFreq[otu]++;
		}
		updateOTUCountData(otuFreq, otuBySeqLookUp, aanI, anP, anI, cumCount);
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "seqNoise", "finishOTUData");
		exit(1);
	}
}

/**************************************************************************************************/

int seqNoise::getLastMatch(char direction, vector<vector<char> >& alignMoves, int i, int j, vector<int>& seqA, vector<int>& seqB){
	try{
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
		m->errorOut(e, "seqNoise", "getLastMatch");
		exit(1);
	}
}

/**************************************************************************************************/

int seqNoise::countDiffs(vector<int> query, vector<int> ref){
	try {
		//double MATCH = 5.0;
		//double MISMATCH = -2.0;
		//double GAP = -2.0;
		
		vector<vector<double> > correctMatrix(4);
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
		
		int queryLength = query.size();
		int refLength = ref.size();
		
		vector<vector<double> > alignMatrix(queryLength + 1);
		vector<vector<char> > alignMoves(queryLength + 1);
		
		for(int i=0;i<=queryLength;i++){
			if (m->getControl_pressed()) { return 0; }
			alignMatrix[i].resize(refLength + 1, 0);
			alignMoves[i].resize(refLength + 1, 'x');
		}
		
		for(int i=0;i<=queryLength;i++){
			if (m->getControl_pressed()) { return 0; }
			alignMatrix[i][0] = 15.0 * i;
			alignMoves[i][0] = 'u';
		}
		
		for(int i=0;i<=refLength;i++){
			if (m->getControl_pressed()) { return 0; }
			alignMatrix[0][i] = 15.0 * i;
			alignMoves[0][i] = 'l';
		}
		
		for(int i=1;i<=queryLength;i++){
			if (m->getControl_pressed()) { return 0; }
			for(int j=1;j<=refLength;j++){
				
				double nogap;		
				nogap = alignMatrix[i-1][j-1] + correctMatrix[query[i-1]][ref[j-1]];
				
				
				double gap;
				double left;
				if(i == queryLength){ //terminal gap
					left = alignMatrix[i][j-1];
				}
				else{
					if(ref[j-1] == getLastMatch('l', alignMoves, i, j, query, ref)){
						gap = 4.0;
					}
					else{
						gap = 15.0;
					}
					
					left = alignMatrix[i][j-1] + gap;
				}
				
				
				double up;
				if(j == refLength){ //terminal gap
					up = alignMatrix[i-1][j];
				}
				else{
					
					if(query[i-1] == getLastMatch('u', alignMoves, i, j, query, ref)){
						gap = 4.0;
					}
					else{
						gap = 15.0;
					}
					
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
		
		int i = queryLength;
		int j = refLength;
		int diffs = 0;
		
		//	string alignA = "";
		//	string alignB = "";
		//	string bases = "ACTG";
		
		while(i > 0 && j > 0){
			if (m->getControl_pressed()) { return 0; }
			if(alignMoves[i][j] == 'd'){
				//			alignA = bases[query[i-1]] + alignA;
				//			alignB = bases[ref[j-1]] + alignB;
				
				if(query[i-1] != ref[j-1])	{	diffs++;	}
				
				i--;
				j--;
			}
			else if(alignMoves[i][j] == 'u'){
				if(j != refLength){
					//				alignA = bases[query[i-1]] + alignA;
					//				alignB = '-' + alignB;
					
					diffs++;
				}
				
				i--;
			}
			else if(alignMoves[i][j] == 'l'){
				if(i != queryLength){
					//				alignA = '-' + alignA;
					//				alignB = bases[ref[j-1]] + alignB;
					
					diffs++;
				}
				
				j--;
			}
		}
				
		return diffs;
	}
	catch(exception& e) {
		m->errorOut(e, "seqNoise", "countDiffs");
		exit(1);
	}
	
}

/**************************************************************************************************/

vector<int> seqNoise::convertSeq(string bases){
	try {
		vector<int> numbers(bases.length(), -1);
		
		for(int i=0;i<bases.length();i++){
			if (m->getControl_pressed()) { return numbers; }
			
			char b = bases[i];
			
			if(b == 'A')	{	numbers[i] = 0;	}
			else if(b=='C')	{	numbers[i] = 1;	}
			else if(b=='T')	{	numbers[i] = 2;	}
			else if(b=='G')	{	numbers[i] = 3;	}
			else			{	numbers[i] = 0;	}
		}
		
		return numbers;
	}
	catch(exception& e) {
		m->errorOut(e, "seqNoise", "convertSeq");
		exit(1);
	}
}

/**************************************************************************************************/

string seqNoise::degapSeq(string aligned){
	try {
		string unaligned = "";
		
		for(int i=0;i<aligned.length();i++){
			
			if (m->getControl_pressed()) { return ""; }
			
			if(aligned[i] != '-' && aligned[i] != '.'){
				unaligned += aligned[i];
			}
		}
		
		return unaligned;
	}
	catch(exception& e) {
		m->errorOut(e, "seqNoise", "degapSeq");
		exit(1);
	}
}

/**************************************************************************************************/

int seqNoise::writeOutput(string fastaFileName, string namesFileName, string uMapFileName, vector<int> finalTau, vector<int> centroids, vector<int> otuData, vector<string> sequences, vector<string> uniqueNames, vector<string> redundantNames, vector<int> seqFreq, vector<double>& distances){
	try {
		int numOTUs = finalTau.size();
		int numSeqs = uniqueNames.size();
		
		ofstream fastaFile(fastaFileName.c_str());
		ofstream namesFile(namesFileName.c_str());
		ofstream uMapFile(uMapFileName.c_str());
		
		vector<int> maxSequenceAbund(numOTUs, 0);
		vector<int> maxSequenceIndex(numOTUs, 0);
		
		for(int i=0;i<numSeqs;i++){
			if (m->getControl_pressed()) { return 0; }
			if(maxSequenceAbund[otuData[i]] < seqFreq[i]){
				maxSequenceAbund[otuData[i]] = seqFreq[i];
				maxSequenceIndex[otuData[i]] = i;
			}
		}
		
		int count = 1;
		
		for(int i=0;i<numOTUs;i++){
			if (m->getControl_pressed()) { return 0; }
			
			if(finalTau[i] > 0){
				
				if((maxSequenceIndex[i] != centroids[i]) && util.isEqual(distances[maxSequenceIndex[i]*numSeqs + centroids[i]], 0)){
					centroids[i] = maxSequenceIndex[i];
				}
				
				int index = centroids[i];
				
				fastaFile << '>' << uniqueNames[index] << endl << sequences[index] << endl;
				namesFile << uniqueNames[index] << '\t';
				
				string refSeq = sequences[index];
				string redundantSeqs = redundantNames[index];;
				
				
				vector<freqData> frequencyData;
				
				for(int j=0;j<numSeqs;j++){
					if(otuData[j] == i && j != index){
						frequencyData.push_back(freqData(j, seqFreq[j]));
					}
				}
				sort(frequencyData.rbegin(), frequencyData.rend());
				
				string refDegap = degapSeq(refSeq);
				vector<int> rUnalign = convertSeq(refDegap);
				
				uMapFile << "ideal_seq_" << count << '\t' << finalTau[i] << endl;
				uMapFile << uniqueNames[index] << '\t' << seqFreq[index] << "\t0\t" << refDegap << endl;
				
				
				for(int j=0;j<frequencyData.size();j++){
					if (m->getControl_pressed()) { return 0; }
					redundantSeqs += ',' + redundantNames[frequencyData[j].index];
					
					uMapFile << uniqueNames[frequencyData[j].index] << '\t' << seqFreq[frequencyData[j].index] << '\t';
					
					string querySeq = sequences[frequencyData[j].index];
					
					string queryDegap = degapSeq(querySeq);
					vector<int> qUnalign = convertSeq(queryDegap);
					
					int udiffs = countDiffs(qUnalign, rUnalign);
					uMapFile << udiffs << '\t' << queryDegap << endl;
					
				}					
				
				uMapFile << endl;
				namesFile << redundantSeqs << endl;
				count++;
				
			}
		}
		fastaFile.close();
		namesFile.close();
		uMapFile.close();
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "seqNoise", "writeOutput");
		exit(1);
	}
}

/**************************************************************************************************

int main(int argc, char *argv[]){
	
	double sigma = 100;
	sigma = atof(argv[5]);
	
	double cutOff = 0.08;
	int minIter = 10;
	int maxIter = 1000;
	double minDelta = 1e-6;
	
	string sequenceFileName = argv[1];
	string fileNameStub = sequenceFileName.substr(0,sequenceFileName.find_last_of('.')) + ".shhh";
	
	vector<string> sequences;
	getSequenceData(sequenceFileName, sequences);
	
	int numSeqs = sequences.size();
	
	vector<string> uniqueNames(numSeqs);
	vector<string> redundantNames(numSeqs);
	vector<int> seqFreq(numSeqs);
	
	string namesFileName = argv[4];
	getRedundantNames(namesFileName, uniqueNames, redundantNames, seqFreq);
	
	string distFileName = argv[2];
	vector<double> distances(numSeqs * numSeqs);
	getDistanceData(distFileName, distances);
	
	string listFileName = argv[3];
	vector<int> otuData(numSeqs);
	vector<int> otuFreq;
	vector<vector<int> > otuBySeqLookUp;
	
	getListData(listFileName, cutOff, otuData, otuFreq, otuBySeqLookUp);
	
	int numOTUs = otuFreq.size();
	
	vector<double> weights(numOTUs, 0);
	vector<int> change(numOTUs, 1);
	vector<int> centroids(numOTUs, -1);
	vector<int> cumCount(numOTUs, 0);
	
	vector<double> tau(numSeqs, 1);
	vector<int> anP(numSeqs, 0);
	vector<int> anI(numSeqs, 0);
	vector<int> anN(numSeqs, 0);
	vector<vector<int> > aanI = otuBySeqLookUp;
	
	int numIters = 0;
	double maxDelta = 1e6;
	
	while(numIters < minIter || ((maxDelta > minDelta) && (numIters < maxIter))){
		
		updateOTUCountData(otuFreq, otuBySeqLookUp, aanI, anP, anI, cumCount);
		maxDelta = calcNewWeights(weights, seqFreq, anI, cumCount, anP, otuFreq, tau);
		
		calcCentroids(anI, anP, change, centroids, cumCount, distances, seqFreq, otuFreq, tau);
		checkCentroids(weights, centroids);
		
		otuFreq.assign(numOTUs, 0);
		
		int total = 0;
		
		for(int i=0;i<numSeqs;i++){
			double offset = 1e6;
			double norm = 0.0000;
			double minWeight = MIN_WEIGHT;
			vector<double> currentTau(numOTUs);
			
			for(int j=0;j<numOTUs;j++){
				if(weights[j] > minWeight && distances[i * numSeqs+centroids[j]] < offset){
					offset = distances[i * numSeqs+centroids[j]];
				}
			}
			
			for(int j=0;j<numOTUs;j++){
				if(weights[j] > minWeight){
					currentTau[j] = exp(sigma * (-distances[(i * numSeqs + centroids[j])] + offset)) * weights[j];
					norm += currentTau[j];
				}
				else{
					currentTau[j] = 0.0000;
				}
			}			
			
			for(int j=0;j<numOTUs;j++){
				currentTau[j] /= norm;
			}
			
			for(int j=0;j<numOTUs;j++){
				
				if(currentTau[j] > MIN_TAU){
					int oldTotal = total;
					total++;
					
					tau.resize(oldTotal+1);
					tau[oldTotal] = currentTau[j];
					otuBySeqLookUp[j][otuFreq[j]] = oldTotal;
					aanI[j][otuFreq[j]] = i;
					otuFreq[j]++;
					
				}
			}
			
			anP.resize(total);
			anI.resize(total);
		}
		
		numIters++;
	}
	
	updateOTUCountData(otuFreq, otuBySeqLookUp, aanI, anP, anI, cumCount);
	
	vector<double> percentage(numSeqs);
	setUpOTUData(otuData, percentage, cumCount, tau, otuFreq, anP, anI);
	finishOTUData(otuData, otuFreq, anP, anI, cumCount, otuBySeqLookUp, aanI, tau);
	
	change.assign(numOTUs, 1);
	calcCentroids(anI, anP, change, centroids, cumCount, distances, seqFreq, otuFreq, tau);
	
	
	vector<int> finalTau(numOTUs, 0);
	for(int i=0;i<numSeqs;i++){
		finalTau[otuData[i]] += int(seqFreq[i]);
	}
	
	writeOutput(fileNameStub, finalTau, centroids, otuData, sequences, uniqueNames, redundantNames, seqFreq, distances);
	
	return 0;
}

**************************************************************************************************/
