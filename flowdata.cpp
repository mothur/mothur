/*
 *  flowdata.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 12/22/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "flowdata.h"

//**********************************************************************************************************************

FlowData::FlowData(){}

//**********************************************************************************************************************

FlowData::FlowData(ifstream& flowFile, float signal, float noise, int maxHomoP){

	try {
		m = MothurOut::getInstance();

		baseFlow = "TACG";
		seqName = "";
		numFlows = 0;
		locationString = "";
		seqLength = 0;
		
		string lengthString;
		string flowString;
		
		flowFile >> seqName >> locationString >> lengthString >> flowString;

		convert(lengthString.substr(7), seqLength);
		convert(flowString.substr(9), numFlows);

		flowData.resize(numFlows);
		
		if (seqName == "") {
			m->mothurOut("Error reading quality file, name blank at position, " + toString(flowFile.tellg()));
			m->mothurOutEndLine(); 
		}
		else{
			seqName = seqName.substr(1);
			for(int i=0;i<numFlows;i++)	{	flowFile >> flowData[i];	}
		}
		
		findDeadSpot(signal, noise, maxHomoP);
		translateFlow();
		
		m->gobble(flowFile);
	}
	catch(exception& e) {
		m->errorOut(e, "FlowData", "FlowData");
		exit(1);
	}
	
}

//**********************************************************************************************************************

void FlowData::findDeadSpot(float signalIntensity, float noiseIntensity, int maxHomoP){
	try{
		
		int currLength = 0;
		float maxIntensity = (float) maxHomoP + 0.49;
		
		deadSpot = 0;
		while(currLength < seqLength + 4){
			int signal = 0;
			int noise = 0;
			
			for(int i=0;i<4;i++){
				float intensity = flowData[i + 4 * deadSpot];
				if(intensity > signalIntensity){
					signal++;

					if(intensity  < noiseIntensity || intensity > maxIntensity){
						noise++;
					}
				}
				currLength += (int)(intensity+0.5);
			}

			if(noise > 0 || signal == 0){
				break;
			}
		
			deadSpot++;
		}
		deadSpot *= 4;
		seqLength = currLength;
		
	}
	catch(exception& e) {
		m->errorOut(e, "FlowData", "findDeadSpot");
		exit(1);
	}
}

//**********************************************************************************************************************

void FlowData::translateFlow(){
	
	try{
		sequence = "";
		for(int i=0;i<deadSpot;i++){
			int intensity = (int)(flowData[i] + 0.5);
			char base = baseFlow[i % 4];
			
			for(int j=0;j<intensity;j++){
				sequence += base;
			}
		}

		if(sequence.size() > 4){
			sequence = sequence.substr(4);
		}
		else{
			sequence = "NNNN";
		}
	}
	catch(exception& e) {
		m->errorOut(e, "FlowData", "translateFlow");
		exit(1);
	}
}

//**********************************************************************************************************************

void FlowData::capFlows(int maxFlows){
	
	try{
		
		numFlows = maxFlows;
		if(deadSpot > maxFlows){	deadSpot = maxFlows;	}
		
	}
	catch(exception& e) {
		m->errorOut(e, "FlowData", "capFlows");
		exit(1);
	}
}

//**********************************************************************************************************************

bool FlowData::hasMinFlows(int minFlows){
	
	try{
		bool pastMin = 0;
		
		if(deadSpot >= minFlows){	pastMin = 1;	}
		return pastMin;
	}
	catch(exception& e) {
		m->errorOut(e, "FlowData", "hasMinFlows");
		exit(1);
	}
}

//**********************************************************************************************************************

Sequence FlowData::getSequence(){
	
	try{
		return Sequence(seqName, sequence);
	}
	catch(exception& e) {
		m->errorOut(e, "FlowData", "getSequence");
		exit(1);
	}
}

//**********************************************************************************************************************

int FlowData::getSeqLength(){
	
	try{
		return seqLength;		
	}
	catch(exception& e) {
		m->errorOut(e, "FlowData", "getSeqLength");
		exit(1);
	}
}

//**********************************************************************************************************************

void FlowData::printFlows(ofstream& outFlowFile){
	try{
	//	outFlowFile << '>' << seqName << locationString << " length=" << seqLength << " numflows=" << maxFlows << endl;
		outFlowFile << seqName << ' ' << deadSpot << ' ' << setprecision(2);

		for(int i=0;i<numFlows;i++){
			outFlowFile << flowData[i] << ' ';
		}
		outFlowFile << endl;
	}
	catch(exception& e) {
		m->errorOut(e, "FlowData", "printFlows");
		exit(1);
	}
}

//**********************************************************************************************************************

void FlowData::printFlows(ofstream& outFlowFile, string scrapCode){
	try{
	
	//	outFlowFile << '>' << seqName << locationString << " length=" << seqLength << " numflows=" << maxFlows << endl;

		outFlowFile << seqName << '|' << scrapCode << ' ' << deadSpot << ' ' << setprecision(2);
		
		for(int i=0;i<numFlows;i++){
			outFlowFile << flowData[i] << ' ';
		}
		outFlowFile << endl;
	}
	catch(exception& e) {
		m->errorOut(e, "FlowData", "printFlows");
		exit(1);
	}
}

//**********************************************************************************************************************

void FlowData::printFASTA(ofstream& outFASTA){
	try{
		
		outFASTA << '>' << seqName << endl;
		outFASTA << sequence << endl;

	}
	catch(exception& e) {
		m->errorOut(e, "FlowData", "printFlows");
		exit(1);
	}
}


//**********************************************************************************************************************
