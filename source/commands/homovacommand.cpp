/*
 *  homovacommand.cpp
 *  mothur
 *
 *  Created by westcott on 2/8/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "homovacommand.h"
#include "groupmap.h"
#include "readphylipvector.h"
#include "sharedutilities.h"
#include "designmap.h"

//**********************************************************************************************************************
vector<string> HomovaCommand::setParameters(){	
	try {
		CommandParameter pdesign("design", "InputTypes", "", "", "none", "none", "none","homova",false,true,true); parameters.push_back(pdesign);
		CommandParameter pphylip("phylip", "InputTypes", "", "", "none", "none", "none","homova",false,true,true); parameters.push_back(pphylip);
        CommandParameter psets("sets", "String", "", "", "", "", "","",false,false); parameters.push_back(psets);
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
		m->errorOut(e, "HomovaCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string HomovaCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "Referenced: Stewart CN, Excoffier L (1996). Assessing population genetic structure and variability with RAPD data: Application to Vaccinium macrocarpon (American Cranberry). J Evol Biol 9: 153-71.\n";
		helpString += "The homova command outputs a .homova file. \n";
		helpString += "The homova command parameters are phylip, iters, sets and alpha.  The phylip and design parameters are required, unless valid current files exist.\n";
		helpString += "The design parameter allows you to assign your samples to groups when you are running homova. It is required. \n";
		helpString += "The design file looks like the group file.  It is a 2 column tab delimited file, where the first column is the sample name and the second column is the group the sample belongs to.\n";
        helpString += "The sets parameter allows you to specify which of the sets in your designfile you would like to analyze. The set names are separated by dashes. THe default is all sets in the designfile.\n";
		helpString += "The iters parameter allows you to set number of randomization for the P value.  The default is 1000. \n";
		helpString += "The homova command should be in the following format: homova(phylip=file.dist, design=file.design).\n";
		helpString += "Note: No spaces between parameter labels (i.e. iters), '=' and parameters (i.e. 1000).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "HomovaCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string HomovaCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "homova") {  pattern = "[filename],homova"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "HomovaCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
HomovaCommand::HomovaCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["homova"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "HomovaCommand", "HomovaCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

HomovaCommand::HomovaCommand(string option) {
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
			outputTypes["homova"] = tempOutNames;
			
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
				//if there is a current phylip file, use it
				phylipFileName = m->getPhylipFile(); 
				if (phylipFileName != "") { m->mothurOut("Using " + phylipFileName + " as input file for the phylip parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current phylip file and the phylip parameter is required."); m->mothurOutEndLine(); abort = true; }
				
			}else { m->setPhylipFile(phylipFileName); }	
			
			//check for required parameters
			designFileName = validParameter.validFile(parameters, "design", true);
			if (designFileName == "not open") { abort = true; }
			else if (designFileName == "not found") {
				//if there is a current design file, use it
				designFileName = m->getDesignFile(); 
				if (designFileName != "") { m->mothurOut("Using " + designFileName + " as input file for the design parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current design file and the design parameter is required."); m->mothurOutEndLine(); abort = true; }								
			}else { m->setDesignFile(designFileName); }	
			
			string temp = validParameter.validFile(parameters, "iters", false);
			if (temp == "not found") { temp = "1000"; }
			m->mothurConvert(temp, iters); 
			
			temp = validParameter.validFile(parameters, "alpha", false);
			if (temp == "not found") { temp = "0.05"; }
			m->mothurConvert(temp, experimentwiseAlpha); 
            
            string sets = validParameter.validFile(parameters, "sets", false);			
			if (sets == "not found") { sets = ""; }
			else { 
				m->splitAtDash(sets, Sets);
			}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "HomovaCommand", "HomovaCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int HomovaCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//read design file
		designMap = new DesignMap(designFileName);
		
		if (outputDir == "") { outputDir = m->hasPath(phylipFileName); }
		
		//read in distance matrix and square it
		ReadPhylipVector readMatrix(phylipFileName);
		vector<string> sampleNames = readMatrix.read(distanceMatrix);
		
        if (Sets.size() != 0) { //user selected sets, so we want to remove the samples not in those sets
            SharedUtil util; 
            vector<string> dGroups = designMap->getCategory();
            util.setGroups(Sets, dGroups);  
            
            for(int i=0;i<distanceMatrix.size();i++){
                
                if (m->control_pressed) { delete designMap; return 0; }
                
                string group = designMap->get(sampleNames[i]);
                
                if (group == "not found") {
                    m->mothurOut("[ERROR]: " + sampleNames[i] + " is not in your design file, please correct."); m->mothurOutEndLine(); m->control_pressed = true;
                }else if (!m->inUsersGroups(group, Sets)){  //not in set we want remove it
                    //remove from all other rows
                    for(int j=0;j<distanceMatrix.size();j++){
                        distanceMatrix[j].erase(distanceMatrix[j].begin()+i);
                    }
                    distanceMatrix.erase(distanceMatrix.begin()+i);
                    sampleNames.erase(sampleNames.begin()+i);
                    i--;
                }
            }
        }

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
				m->mothurOut("[ERROR]: " + sampleNames[i] + " is not in your design file, please correct."); m->mothurOutEndLine(); m->control_pressed = true;
			}else { origGroupSampleMap[group].push_back(i); }
		}
		int numGroups = origGroupSampleMap.size();
		
		if (m->control_pressed) { delete designMap; return 0; }
		
		//create a new filename
		ofstream HOMOVAFile;
        map<string, string> variables; 
		variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(phylipFileName));
		string HOMOVAFileName = getOutputFileName("homova", variables);				
		m->openOutputFile(HOMOVAFileName, HOMOVAFile);
		outputNames.push_back(HOMOVAFileName); outputTypes["homova"].push_back(HOMOVAFileName);
		
		HOMOVAFile << "HOMOVA\tBValue\tP-value\tSSwithin/(Ni-1)_values" << endl;
		m->mothurOut("HOMOVA\tBValue\tP-value\tSSwithin/(Ni-1)_values\n");
		
		double fullHOMOVAPValue = runHOMOVA(HOMOVAFile, origGroupSampleMap, experimentwiseAlpha);

		if(fullHOMOVAPValue <= experimentwiseAlpha && numGroups > 2){
			
			int numCombos = numGroups * (numGroups-1) / 2;
			double pairwiseAlpha = experimentwiseAlpha / (double) numCombos;
			
			map<string, vector<int> >::iterator itA;
			map<string, vector<int> >::iterator itB;
			
			for(itA=origGroupSampleMap.begin();itA!=origGroupSampleMap.end();itA++){
				itB = itA;itB++;
				for(;itB!=origGroupSampleMap.end();itB++){
					map<string, vector<int> > pairwiseGroupSampleMap;
					pairwiseGroupSampleMap[itA->first] = itA->second;
					pairwiseGroupSampleMap[itB->first] = itB->second;
					
					runHOMOVA(HOMOVAFile, pairwiseGroupSampleMap, pairwiseAlpha);
				}			
			}
			HOMOVAFile << endl;
			m->mothurOutEndLine();
			
			m->mothurOut("Experiment-wise error rate: " + toString(experimentwiseAlpha) + '\n');
			m->mothurOut("Pair-wise error rate (Bonferroni): " + toString(pairwiseAlpha) + '\n');
		}
		else{
			m->mothurOut("Experiment-wise error rate: " + toString(experimentwiseAlpha) + '\n');
		}
		
		m->mothurOut("If you have borderline P-values, you should try increasing the number of iterations\n");
		
		delete designMap;
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "HomovaCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************

double HomovaCommand::runHOMOVA(ofstream& HOMOVAFile, map<string, vector<int> > groupSampleMap, double alpha){
	try {
		map<string, vector<int> >::iterator it;
		int numGroups = groupSampleMap.size();
		
		vector<double> ssWithinOrigVector;
		double bValueOrig = calcBValue(groupSampleMap, ssWithinOrigVector);
		
		double counter = 0;
		for(int i=0;i<iters;i++){
			vector<double> ssWithinRandVector;
			map<string, vector<int> > randomizedGroup = getRandomizedGroups(groupSampleMap);
			double bValueRand = calcBValue(randomizedGroup, ssWithinRandVector);
			if(bValueRand >= bValueOrig){	counter++;	}
		}
		
		double pValue = (double) counter / (double) iters;
		string pString = "";
		if(pValue < 1/(double)iters){	pString = '<' + toString(1/(double)iters);	}
		else						{	pString = toString(pValue);					}
		
		
		//print homova table
		it = groupSampleMap.begin();
		HOMOVAFile << it->first;
		m->mothurOut(it->first);
		it++;
		for(;it!=groupSampleMap.end();it++){
			HOMOVAFile << '-' << it->first;
			m->mothurOut('-' + it->first);
		}

		HOMOVAFile << '\t' << bValueOrig << '\t' << pString;
		m->mothurOut('\t' + toString(bValueOrig) + '\t' + pString);
		
		if(pValue < alpha){
			HOMOVAFile << "*";
			m->mothurOut("*");
		}

		for(int i=0;i<numGroups;i++){
			HOMOVAFile << '\t' << ssWithinOrigVector[i];
			m->mothurOut('\t' + toString(ssWithinOrigVector[i]));
		}
		HOMOVAFile << endl;
		m->mothurOutEndLine();
		
		return pValue;	
	}
	catch(exception& e) {
		m->errorOut(e, "HomovaCommand", "runHOMOVA");
		exit(1);
	}
}

//**********************************************************************************************************************

double HomovaCommand::calcSigleSSWithin(vector<int> sampleIndices) {
	try {
		double ssWithin = 0.0;
		int numSamplesInGroup = sampleIndices.size();
		
		for(int i=0;i<numSamplesInGroup;i++){
			int row = sampleIndices[i];
			
			for(int j=0;j<numSamplesInGroup;j++){
				int col = sampleIndices[j];
				
				if(col < row){
					ssWithin += distanceMatrix[row][col];
				}
				
			}
		}
		
		ssWithin /= numSamplesInGroup;
		return ssWithin;
	}
	catch(exception& e) {
		m->errorOut(e, "HomovaCommand", "calcSigleSSWithin");
		exit(1);
	}
}

//**********************************************************************************************************************

double HomovaCommand::calcBValue(map<string, vector<int> > groupSampleMap, vector<double>& ssWithinVector) {
	try {

		map<string, vector<int> >::iterator it;
		
		double numGroups = (double)groupSampleMap.size();
		ssWithinVector.resize(numGroups, 0);
		
		double totalNumSamples = 0;
		double ssWithinFull;
		double secondTermSum = 0;
		double inverseOneMinusSum = 0;
		int index = 0;
		
		ssWithinVector.resize(numGroups, 0);
		for(it = groupSampleMap.begin();it!=groupSampleMap.end();it++){
			int numSamplesInGroup = it->second.size();
			totalNumSamples += numSamplesInGroup;
			
			ssWithinVector[index] = calcSigleSSWithin(it->second);
			ssWithinFull += ssWithinVector[index];
			
			secondTermSum += (numSamplesInGroup - 1) * log(ssWithinVector[index] / (double)(numSamplesInGroup - 1));
			inverseOneMinusSum += 1.0 / (double)(numSamplesInGroup - 1);
			
			ssWithinVector[index] /= (double)(numSamplesInGroup - 1); //this line is only for output purposes to scale SSw by the number of samples in the group
			index++;
		}
		
		double B = (totalNumSamples - numGroups) * log(ssWithinFull/(totalNumSamples-numGroups)) - secondTermSum;
		double denomintor = 1 + 1.0/(3.0 * (numGroups - 1.0)) * (inverseOneMinusSum - 1.0 / (double) (totalNumSamples - numGroups));
		B /= denomintor;
		
		return B;
		
	}
	catch(exception& e) {
		m->errorOut(e, "HomovaCommand", "calcBValue");
		exit(1);
	}
}

//**********************************************************************************************************************

map<string, vector<int> > HomovaCommand::getRandomizedGroups(map<string, vector<int> > origMapping){
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
		m->errorOut(e, "AmovaCommand", "randomizeGroups");
		exit(1);
	}
}

//**********************************************************************************************************************


