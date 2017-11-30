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
#include "designmap.h"

//**********************************************************************************************************************
vector<string> AnosimCommand::setParameters(){	
	try {
		CommandParameter pdesign("design", "InputTypes", "", "", "none", "none", "none","anosim",false,true,true); parameters.push_back(pdesign);
		CommandParameter pphylip("phylip", "InputTypes", "", "", "none", "none", "none","anosim",false,true,true); parameters.push_back(pphylip);
		CommandParameter piters("iters", "Number", "", "1000", "", "", "","",false,false); parameters.push_back(piters);
		CommandParameter palpha("alpha", "Number", "", "0.05", "", "", "","",false,false); parameters.push_back(palpha);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);

		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string AnosimCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "Referenced: Clarke, K. R. (1993). Non-parametric multivariate analysis of changes in community structure.   _Australian Journal of Ecology_ 18, 117-143.\n";
		helpString += "The anosim command outputs a .anosim file. \n";
		helpString += "The anosim command parameters are phylip, iters, and alpha.  The phylip and design parameters are required, unless you have valid current files.\n";
		helpString += "The design parameter allows you to assign your samples to groups when you are running anosim. It is required. \n";
		helpString += "The design file looks like the group file.  It is a 2 column tab delimited file, where the first column is the sample name and the second column is the group the sample belongs to.\n";
		helpString += "The iters parameter allows you to set number of randomization for the P value.  The default is 1000. \n";
		helpString += "The anosim command should be in the following format: anosim(phylip=file.dist, design=file.design).\n";
		helpString += "Note: No spaces between parameter labels (i.e. iters), '=' and parameters (i.e. 1000).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string AnosimCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "anosim") {  pattern = "[filename],anosim"; } //makes file like: amazon.align
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "AnosimCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
AnosimCommand::AnosimCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
		vector<string> tempOutNames;
		outputTypes["anosim"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "AnosimCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

AnosimCommand::AnosimCommand(string option) {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
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
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";	}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("design");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["design"] = inputDir + it->second;		}
				}
				
				it = parameters.find("phylip");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["phylip"] = inputDir + it->second;		}
				}
			}
			
			phylipFileName = validParameter.validFile(parameters, "phylip");
			if (phylipFileName == "not open") { phylipFileName = ""; abort = true; }
			else if (phylipFileName == "not found") { 
				//if there is a current phylip file, use it
				phylipFileName = current->getPhylipFile(); 
				if (phylipFileName != "") { m->mothurOut("Using " + phylipFileName + " as input file for the phylip parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current phylip file and the phylip parameter is required."); m->mothurOutEndLine(); abort = true; }
				
			}else { current->setPhylipFile(phylipFileName); }	
			
			//check for required parameters
			designFileName = validParameter.validFile(parameters, "design");
			if (designFileName == "not open") { designFileName = ""; abort = true; }
			else if (designFileName == "not found") {
				//if there is a current design file, use it
				designFileName = current->getDesignFile(); 
				if (designFileName != "") { m->mothurOut("Using " + designFileName + " as input file for the design parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current design file and the design parameter is required."); m->mothurOutEndLine(); abort = true; }								
			}else { current->setDesignFile(designFileName); }	
			
			string temp = validParameter.valid(parameters, "iters");
			if (temp == "not found") { temp = "1000"; }
			util.mothurConvert(temp, iters); 
			
			temp = validParameter.valid(parameters, "alpha");
			if (temp == "not found") { temp = "0.05"; }
			util.mothurConvert(temp, experimentwiseAlpha); 
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "AnosimCommand", "AnosimCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int AnosimCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		//read design file
		designMap = new DesignMap(designFileName);
		
		if (outputDir == "") { outputDir = util.hasPath(phylipFileName); }
		
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
			string group = designMap->get(sampleNames[i]);
			
			if (group == "not found") {
				m->mothurOut("[ERROR]: " + sampleNames[i] + " is not in your design file, please correct."); m->mothurOutEndLine(); m->setControl_pressed(true);
			}else { origGroupSampleMap[group].push_back(i); }
		}
		int numGroups = origGroupSampleMap.size();
		
		if (m->getControl_pressed()) { delete designMap; return 0; }
		
		//create a new filename
		ofstream ANOSIMFile;
        map<string, string> variables; variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(phylipFileName));
		string ANOSIMFileName = getOutputFileName("anosim", variables);	
        
		util.openOutputFile(ANOSIMFileName, ANOSIMFile);
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
		
		
		//map<string, vector<int> >::iterator it=groupSampleMap.begin();
        vector<string> sampleNames;
        for(map<string, vector<int> >::iterator it = groupSampleMap.begin();it!=groupSampleMap.end();it++){ sampleNames.push_back(it->first); }
        string output = util.getStringFromVector(sampleNames, "-");
        
		m->mothurOut(output + '\t' + toString(RValue) + '\t' + pString);
		ANOSIMFile << output << '\t' << RValue << '\t' << pString;

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
		
		m->mothurRandomShuffle(sampleIndices);
		
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



