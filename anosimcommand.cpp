/*
 *  anosimcommand.cpp
 *  mothur
 *
 *  Created by westcott on 2/14/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "anosimcommand.h"
#include "inputdata.h"
#include "readphylipvector.h"

//**********************************************************************************************************************
vector<string> AnosimCommand::getValidParameters(){	
	try {
		string Array[] =  {"outputdir","iters","phylip","design", "alpha","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
AnosimCommand::AnosimCommand(){	
	try {
		abort = true; calledHelp = true; 
		vector<string> tempOutNames;
		outputTypes["anosim"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "AnosimCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> AnosimCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"design"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> AnosimCommand::getRequiredFiles(){	
	try {
		string Array[] =  {};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************

AnosimCommand::AnosimCommand(string option) {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {"outputdir","iters","phylip","design", "alpha","inputdir"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			map<string,string>::iterator it;
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["anosim"] = tempOutNames;
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("design");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["design"] = inputDir + it->second;		}
				}
				
				it = parameters.find("phylip");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["phylip"] = inputDir + it->second;		}
				}
			}
			
			phylipFileName = validParameter.validFile(parameters, "phylip", true);
			if (phylipFileName == "not open") { phylipFileName = ""; abort = true; }
			else if (phylipFileName == "not found") { 
				phylipFileName = ""; 	
			
				//check currentFiles for a phylip file
				//if (currentFiles->getPhylipFile() != "") {  phylipFileName = currentFiles->getPhylipFile(); m->mothurOut("Using " + phylipFileName + " as phylip file."); m->mothurOutEndLine();
				//}else { m->mothurOut("You must provide an phylip file."); m->mothurOutEndLine(); abort = true;  }
			}	
			
			//check for required parameters
			designFileName = validParameter.validFile(parameters, "design", true);
			if (designFileName == "not open") { abort = true; }
			else if (designFileName == "not found") {
				designFileName = "";
				
				//check currentFiles for a design file
				//if (currentFiles->getDesignFile() != "") {  designFileName = currentFiles->getDesignFile(); m->mothurOut("Using " + designFileName + " as design file."); m->mothurOutEndLine();
				//}else { m->mothurOut("You must provide an design file."); m->mothurOutEndLine(); abort = true;  }
				
			}	
			
			string temp = validParameter.validFile(parameters, "iters", false);
			if (temp == "not found") { temp = "1000"; }
			convert(temp, iters); 
			
			temp = validParameter.validFile(parameters, "alpha", false);
			if (temp == "not found") { temp = "0.05"; }
			convert(temp, experimentwiseAlpha); 
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "AnosimCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void AnosimCommand::help(){
	try {
		m->mothurOut("Referenced: Clarke, K. R. (1993). Non-parametric multivariate analysis of changes in community structure.   _Australian Journal of Ecology_ 18, 117-143.\n");
		m->mothurOut("The anosim command outputs a .anosim file. \n");
		m->mothurOut("The anosim command parameters are phylip, iters, and alpha.  The phylip and design parameters are required.\n");
		m->mothurOut("The design parameter allows you to assign your samples to groups when you are running anosim. It is required. \n");
		m->mothurOut("The design file looks like the group file.  It is a 2 column tab delimited file, where the first column is the sample name and the second column is the group the sample belongs to.\n");
		m->mothurOut("The iters parameter allows you to set number of randomization for the P value.  The default is 1000. \n");
		m->mothurOut("The anosim command should be in the following format: anosim(phylip=file.dist, design=file.design).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. iters), '=' and parameters (i.e. 1000).\n\n");
		
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

AnosimCommand::~AnosimCommand(){}

//**********************************************************************************************************************

int AnosimCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//read design file
		designMap = new GroupMap(designFileName);
		designMap->readDesignMap();
		
		if (outputDir == "") { outputDir = m->hasPath(phylipFileName); }
		
		//read in distance matrix and square it
		ReadPhylipVector readMatrix(phylipFileName);
		vector<string> sampleNames = readMatrix.read(distanceMatrix);
		
		for(int i=0;i<distanceMatrix.size();i++){
			for(int j=0;j<i;j++){
				distanceMatrix[i][j] *= distanceMatrix[i][j];	
			}
		}
		
		//link designMap to rows/columns in distance matrix
		map<string, vector<int> > origGroupSampleMap;
		for(int i=0;i<sampleNames.size();i++){
			origGroupSampleMap[designMap->getGroup(sampleNames[i])].push_back(i);
		}
		int numGroups = origGroupSampleMap.size();
		
		//create a new filename
		ofstream ANOSIMFile;
		string ANOSIMFileName = outputDir + m->getRootName(m->getSimpleName(phylipFileName))  + "anosim";				
		m->openOutputFile(ANOSIMFileName, ANOSIMFile);
		outputNames.push_back(ANOSIMFileName); outputTypes["anosim"].push_back(ANOSIMFileName);
		m->mothurOut("\ncomparison\tR-value\tP-value\n");
		ANOSIMFile << "comparison\tR-value\tP-value\n";
		
		
		double fullANOSIMPValue = runANOSIM(ANOSIMFile, distanceMatrix, origGroupSampleMap, experimentwiseAlpha);
		
		
		if(fullANOSIMPValue <= experimentwiseAlpha && numGroups > 2){

			int numCombos = numGroups * (numGroups-1) / 2;
			double pairwiseAlpha = experimentwiseAlpha / (double) numCombos;

			for(map<string, vector<int> >::iterator itA=origGroupSampleMap.begin();itA!=origGroupSampleMap.end();itA++){
				map<string, vector<int> >::iterator itB = itA;
				itB++;
				for(itB;itB!=origGroupSampleMap.end();itB++){
					
					map<string, vector<int> > subGroupSampleMap;
					
					subGroupSampleMap[itA->first] = itA->second;	string groupA = itA->first;
					subGroupSampleMap[itB->first] = itB->second;	string groupB = itB->first;
			
					vector<int> subIndices;
					for(map<string, vector<int> >::iterator it=subGroupSampleMap.begin();it!=subGroupSampleMap.end();it++){
						subIndices.insert(subIndices.end(), it->second.begin(), it->second.end());
					}
					int subNumSamples = subIndices.size();

					sort(subIndices.begin(), subIndices.end());		
					
					vector<vector<double> > subDistMatrix(distanceMatrix.size());
					for(int i=0;i<distanceMatrix.size();i++){
						subDistMatrix[i].assign(distanceMatrix.size(), -1);
					}

					for(int i=0;i<subNumSamples;i++){
						for(int j=0;j<i;j++){
							subDistMatrix[subIndices[i]][subIndices[j]] = distanceMatrix[subIndices[i]][subIndices[j]];
						}
					}

					runANOSIM(ANOSIMFile, subDistMatrix, subGroupSampleMap, pairwiseAlpha);

				}
			}
			
			m->mothurOut("\nExperiment-wise error rate: " + toString(experimentwiseAlpha) + '\n');
			m->mothurOut("Pair-wise error rate (Bonferroni): " + toString(pairwiseAlpha) + '\n');
		}
		else{
			m->mothurOut("\nExperiment-wise error rate: " + toString(experimentwiseAlpha) + '\n');
		}
		m->mothurOut("If you have borderline P-values, you should try increasing the number of iterations\n");
		ANOSIMFile.close();
		
			
		delete designMap;
				
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

double AnosimCommand::runANOSIM(ofstream& ANOSIMFile, vector<vector<double> > dMatrix, map<string, vector<int> > groupSampleMap, double alpha) {
	try {

		
		vector<vector<double> > rankMatrix = convertToRanks(dMatrix);
		double RValue = calcR(rankMatrix, groupSampleMap);
		
		int pCount = 0;
		for(int i=0;i<iters;i++){
			map<string, vector<int> > randGroupSampleMap = getRandomizedGroups(groupSampleMap);
			double RValueRand = calcR(rankMatrix, randGroupSampleMap);
			if(RValue <= RValueRand){	pCount++;	}
		}

		double pValue = (double)pCount / (double) iters;
		string pString = "";
		if(pValue < 1/(double)iters){	pString = '<' + toString(1/(double)iters);	}
		else						{	pString = toString(pValue);					}
		
		
		map<string, vector<int> >::iterator it=groupSampleMap.begin();
		m->mothurOut(it->first);
		ANOSIMFile << it->first;
		it++;
		for(it;it!=groupSampleMap.end();it++){
			m->mothurOut('-' + it->first);
			ANOSIMFile << '-' << it->first;
		
		}
		m->mothurOut('\t' + toString(RValue) + '\t' + pString);
		ANOSIMFile << '\t' << RValue << '\t' << pString;

		if(pValue < alpha){
			ANOSIMFile << "*";
			m->mothurOut("*");
		}
		ANOSIMFile << endl;
		m->mothurOutEndLine();
		
		return pValue;
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "calcAnisom");
		exit(1);
	}
}

//**********************************************************************************************************************

double AnosimCommand::calcR(vector<vector<double> > rankMatrix, map<string, vector<int> > groupSampleMap){
	try {

		int numSamples = 0;
		for(map<string, vector<int> >::iterator it=groupSampleMap.begin();it!=groupSampleMap.end();it++){
			numSamples += it->second.size();
		}
		
		
		double within = 0.0;
		int numWithinComps = 0;		
		
		for(map<string, vector<int> >::iterator it=groupSampleMap.begin();it!=groupSampleMap.end();it++){
			vector<int> indices = it->second;
			for(int i=0;i<indices.size();i++){
				for(int j=0;j<i;j++){
					if(indices[i] > indices[j])	{	within += rankMatrix[indices[i]][indices[j]];	}
					else						{	within += rankMatrix[indices[j]][indices[i]];	}
					numWithinComps++;
				}
			}
		}
		
		within /= (float) numWithinComps;
		
		double between = 0.0;
		int numBetweenComps = 0;

		map<string, vector<int> >::iterator itB;
		
		for(map<string, vector<int> >::iterator itA=groupSampleMap.begin();itA!=groupSampleMap.end();itA++){

			for(int i=0;i<itA->second.size();i++){
				int A = itA->second[i];
				map<string, vector<int> >::iterator itB = itA;
				itB++;
				for(itB;itB!=groupSampleMap.end();itB++){
					for(int j=0;j<itB->second.size();j++){
						int B = itB->second[j];
						if(A>B)	{	between += rankMatrix[A][B];	}
						else	{	between += rankMatrix[B][A];	}
						numBetweenComps++;
					}					
				}
				
			}
		}
		
		
		between /= (float) numBetweenComps;
		
		double Rvalue = (between - within)/(numSamples * (numSamples-1) / 4.0);
				
		return Rvalue;
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "calcWithinBetween");
		exit(1);
	}
}

//**********************************************************************************************************************

vector<vector<double> > AnosimCommand::convertToRanks(vector<vector<double> > dist) {
	try {
		vector<seqDist> cells;
		vector<vector<double> > ranks = dist;
		
		for (int i = 0; i < dist.size(); i++) {
			for (int j = 0; j < i; j++) {
				if(dist[i][j] != -1){
					seqDist member(i, j, dist[i][j]);
					cells.push_back(member);
				}
			}
		}
		
		
		//sort distances
		sort(cells.begin(), cells.end(), compareSequenceDistance); 	

		//find ranks of distances
		int index = 0;
		int indexSum = 0;
		for(int i=0;i<cells.size()-1;i++){

			index = i;
			indexSum = i + 1;
			while(dist[cells[index].seq1][cells[index].seq2] == dist[cells[index+1].seq1][cells[index+1].seq2]){
				index++;				
				indexSum += index + 1;
			}
			
			if(index == i){
				ranks[cells[i].seq1][cells[i].seq2] = i+1;
			}
			else{
				double aveIndex = (double)indexSum / (double)(index - i + 1);
				for(int j=i;j<=index;j++){
					ranks[cells[j].seq1][cells[j].seq2] = aveIndex;
				}					
				i = index;
			}
		}
		
		if(indexSum == cells.size() - 1){
			ranks[cells[cells.size()-1].seq1][cells[cells.size()-1].seq2] = indexSum + 1;
		}

		return ranks;
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "convertToRanks");
		exit(1);
	}
}

//**********************************************************************************************************************

map<string, vector<int> > AnosimCommand::getRandomizedGroups(map<string, vector<int> > origMapping){
	try{
		vector<int> sampleIndices;
		vector<int> samplesPerGroup;
		
		map<string, vector<int> >::iterator it;
		for(it=origMapping.begin();it!=origMapping.end();it++){
			vector<int> indices = it->second;
			samplesPerGroup.push_back(indices.size());
			sampleIndices.insert(sampleIndices.end(), indices.begin(), indices.end());
		}
		
		random_shuffle(sampleIndices.begin(), sampleIndices.end());
		
		int index = 0;
		map<string, vector<int> > randomizedGroups = origMapping;
		for(it=randomizedGroups.begin();it!=randomizedGroups.end();it++){
			for(int i=0;i<it->second.size();i++){
				it->second[i] = sampleIndices[index++];				
			}
		}
		
		return randomizedGroups;		
	}
	catch (exception& e) {
		m->errorOut(e, "AnosimCommand", "randomizeGroups");
		exit(1);
	}
}

//**********************************************************************************************************************



