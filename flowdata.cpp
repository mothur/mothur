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

FlowData::~FlowData(){	/*	do nothing	*/	}

//**********************************************************************************************************************

FlowData::FlowData(int numFlows, float signal, float noise, int maxHomoP, string baseFlow) : 
			numFlows(numFlows), signalIntensity(signal), noiseIntensity(noise), maxHomoP(maxHomoP), baseFlow(baseFlow){

	try {
		m = MothurOut::getInstance();

		flowData.assign(numFlows, 0);
//		baseFlow = "TACG";
		seqName = "";
		locationString = "";
	}
	catch(exception& e) {
		m->errorOut(e, "FlowData", "FlowData");
		exit(1);
	}
	
}

//**********************************************************************************************************************

bool FlowData::getNext(ifstream& flowFile){
	
	try {
        seqName = getSequenceName(flowFile);
		flowFile >> endFlow;	
        if (!m->control_pressed) {
            for(int i=0;i<numFlows;i++)	{	flowFile >> flowData[i]; 	}
            updateEndFlow(); 
            translateFlow();
            m->gobble(flowFile);
		}
            
		if(flowFile){	return 1;	}
		else		{	return 0;	}
	}
	catch(exception& e) {
		m->errorOut(e, "FlowData", "getNext");
		exit(1);
	}
	
}
//********************************************************************************************************************
string FlowData::getSequenceName(ifstream& flowFile) {
	try {
		string name = "";
		
        flowFile >> name;
		
		if (name.length() != 0) { 
            for (int i = 0; i < name.length(); i++) {
                if (name[i] == ':') { name[i] = '_'; m->changedSeqNames = true; }
            }
        }else{ m->mothurOut("Error in reading your flowfile, at position " + toString(flowFile.tellg()) + ". Blank name."); m->mothurOutEndLine(); m->control_pressed = true;  }
        
		return name;
	}
	catch(exception& e) {
		m->errorOut(e, "FlowData", "getSequenceName");
		exit(1);
	}
}

//**********************************************************************************************************************

void FlowData::updateEndFlow(){
	try{
		
		//int currLength = 0;
		float maxIntensity = (float) maxHomoP + 0.49;
		
		int deadSpot = 0;
				
		while(deadSpot < endFlow){
			int signal = 0;
			int noise = 0;
			
			for(int i=0;i<4;i++){
				float intensity = flowData[i + deadSpot];
				if(intensity > signalIntensity){
					signal++;

					if(intensity  < noiseIntensity || intensity > maxIntensity){
						noise++;
					}
				}
			}

			if(noise > 0 || signal == 0){
				break;
			}
		
			deadSpot += 4;
		}
		endFlow = deadSpot;

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
		for(int i=0;i<endFlow;i++){
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

void FlowData::capFlows(int mF){
	
	try{
		
		maxFlows = mF;
		if(endFlow > maxFlows){	endFlow = maxFlows;	}	
        translateFlow();
		
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
		if(endFlow >= minFlows){	pastMin = 1;	}

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

void FlowData::printFlows(ofstream& outFlowFile){
	try{
//	outFlowFile << '>' << seqName << locationString << " length=" << seqLength << " numflows=" << maxFlows << endl;
		outFlowFile << seqName << ' ' << endFlow << ' ' << setprecision(2);

		for(int i=0;i<maxFlows;i++){
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
		outFlowFile << seqName << '|' << scrapCode << ' ' << endFlow << ' ' << setprecision(2);
		
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

string FlowData::getName(){
	
	try{
		return seqName;
	}
	catch(exception& e) {
		m->errorOut(e, "FlowData", "getName");
		exit(1);
	}
}

//**********************************************************************************************************************
