/*
 *  amovacommand.cpp
 *  mothur
 *
 *  Created by westcott on 2/7/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "amovacommand.h"
#include "readphylipvector.h"
#include "designmap.h"
#include "sharedutilities.h"


//**********************************************************************************************************************
vector<string> AmovaCommand::setParameters(){	
	try {
		CommandParameter pdesign("design", "InputTypes", "", "", "none", "none", "none","amova",false,true,true); parameters.push_back(pdesign);
        CommandParameter psets("sets", "String", "", "", "", "", "","",false,false); parameters.push_back(psets);
		CommandParameter pphylip("phylip", "InputTypes", "", "", "none", "none", "none","amova",false,true,true); parameters.push_back(pphylip);
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
		m->errorOut(e, "AmovaCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string AmovaCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "Referenced: Anderson MJ (2001). A new method for non-parametric multivariate analysis of variance. Austral Ecol 26: 32-46.";
		helpString += "The amova command outputs a .amova file.";
		helpString += "The amova command parameters are phylip, iters, sets and alpha.  The phylip and design parameters are required, unless you have valid current files.";
		helpString += "The design parameter allows you to assign your samples to groups when you are running amova. It is required.";
		helpString += "The design file looks like the group file.  It is a 2 column tab delimited file, where the first column is the sample name and the second column is the group the sample belongs to.";
        helpString += "The sets parameter allows you to specify which of the sets in your designfile you would like to analyze. The set names are separated by dashes. THe default is all sets in the designfile.\n";
		helpString += "The iters parameter allows you to set number of randomization for the P value.  The default is 1000.";
		helpString += "The amova command should be in the following format: amova(phylip=file.dist, design=file.design).";
		helpString += "Note: No spaces between parameter labels (i.e. iters), '=' and parameters (i.e. 1000).";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "AmovaCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string AmovaCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "amova") {  pattern = "[filename],amova"; } //makes file like: amazon.align
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "AmovaCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
AmovaCommand::AmovaCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["amova"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "AmovaCommand", "AmovaCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
AmovaCommand::AmovaCommand(string option) {
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
			outputTypes["amova"] = tempOutNames;
			
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
			if (designFileName == "not open") { designFileName = ""; abort = true; }
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
		m->errorOut(e, "AmovaCommand", "AmovaCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int AmovaCommand::execute(){
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
		ofstream AMOVAFile;
        map<string, string> variables; variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(phylipFileName));
		string AMOVAFileName = getOutputFileName("amova", variables);	
        
		m->openOutputFile(AMOVAFileName, AMOVAFile);
		outputNames.push_back(AMOVAFileName); outputTypes["amova"].push_back(AMOVAFileName);
		
		double fullANOVAPValue = runAMOVA(AMOVAFile, origGroupSampleMap, experimentwiseAlpha);
		if(fullANOVAPValue <= experimentwiseAlpha && numGroups > 2){
			
			int numCombos = numGroups * (numGroups-1) / 2;
			double pairwiseAlpha = experimentwiseAlpha / (double) numCombos;
			
			map<string, vector<int> >::iterator itA;
			map<string, vector<int> >::iterator itB;
			
			for(itA=origGroupSampleMap.begin();itA!=origGroupSampleMap.end();itA++){
				itB = itA;itB++;
				for(itB;itB!=origGroupSampleMap.end();itB++){
					
					map<string, vector<int> > pairwiseGroupSampleMap;
					pairwiseGroupSampleMap[itA->first] = itA->second;
					pairwiseGroupSampleMap[itB->first] = itB->second;
					
					runAMOVA(AMOVAFile, pairwiseGroupSampleMap, pairwiseAlpha);
				}			
			}
			m->mothurOut("Experiment-wise error rate: " + toString(experimentwiseAlpha) + '\n');
			m->mothurOut("Pair-wise error rate (Bonferroni): " + toString(pairwiseAlpha) + '\n');
		}
		else{
			m->mothurOut("Experiment-wise error rate: " + toString(experimentwiseAlpha) + '\n');
		}
		m->mothurOut("If you have borderline P-values, you should try increasing the number of iterations\n");
		AMOVAFile.close();
		
		delete designMap;
	 
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "AmovaCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************

double AmovaCommand::runAMOVA(ofstream& AMOVAFile, map<string, vector<int> > groupSampleMap, double alpha) {
	try {
		map<string, vector<int> >::iterator it;

		int numGroups = groupSampleMap.size();
		int totalNumSamples = 0;

		for(it = groupSampleMap.begin();it!=groupSampleMap.end();it++){
			totalNumSamples += it->second.size();			
		}

		double ssTotalOrig = calcSSTotal(groupSampleMap);
		double ssWithinOrig = calcSSWithin(groupSampleMap);
		double ssAmongOrig = ssTotalOrig - ssWithinOrig;
		
		double counter = 0;
		for(int i=0;i<iters;i++){
			map<string, vector<int> > randomizedGroup = getRandomizedGroups(groupSampleMap);
			double ssWithinRand = calcSSWithin(randomizedGroup);
			if(ssWithinRand <= ssWithinOrig){	counter++;	}
		}
		
		double pValue = (double)counter / (double) iters;
		string pString = "";
		if(pValue < 1/(double)iters){	pString = '<' + toString(1/(double)iters);	}
		else						{	pString = toString(pValue);					}
		
		
		//print anova table
		it = groupSampleMap.begin();
		AMOVAFile << it->first;
		m->mothurOut(it->first);
		it++;
		for(it;it!=groupSampleMap.end();it++){
			AMOVAFile << '-' << it->first;
			m->mothurOut('-' + it->first);
		}
		
		AMOVAFile << "\tAmong\tWithin\tTotal" << endl;
		m->mothurOut("\tAmong\tWithin\tTotal\n");
		
		AMOVAFile << "SS\t" << ssAmongOrig << '\t' << ssWithinOrig << '\t' << ssTotalOrig << endl;
		m->mothurOut("SS\t" + toString(ssAmongOrig) + '\t' + toString(ssWithinOrig) + '\t' + toString(ssTotalOrig) + '\n');
		
		int dfAmong = numGroups - 1;				double MSAmong = ssAmongOrig / (double) dfAmong;
		int dfWithin = totalNumSamples - numGroups;	double MSWithin = ssWithinOrig / (double) dfWithin;
		int dfTotal = totalNumSamples - 1;			double Fs = MSAmong / MSWithin;
		
		AMOVAFile << "df\t" << dfAmong << '\t' << dfWithin << '\t' << dfTotal << endl;
		m->mothurOut("df\t" + toString(dfAmong) + '\t' + toString(dfWithin) + '\t' + toString(dfTotal) + '\n');

		AMOVAFile << "MS\t" << MSAmong << '\t' << MSWithin << endl << endl;
		m->mothurOut("MS\t" + toString(MSAmong) + '\t' + toString(MSWithin) + "\n\n");

		AMOVAFile << "Fs:\t" << Fs << endl;
		m->mothurOut("Fs:\t" + toString(Fs) + '\n');
		
		AMOVAFile << "p-value: " << pString;
		m->mothurOut("p-value: " + pString);

		if(pValue < alpha){
			AMOVAFile << "*";
			m->mothurOut("*");
		}
		AMOVAFile << endl << endl;
		m->mothurOutEndLine();m->mothurOutEndLine();

		return pValue;
	}
	catch(exception& e) {
		m->errorOut(e, "AmovaCommand", "runAMOVA");
		exit(1);
	}
}

//**********************************************************************************************************************

map<string, vector<int> > AmovaCommand::getRandomizedGroups(map<string, vector<int> > origMapping){
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
		m->errorOut(e, "AmovaCommand", "getRandomizedGroups");
		exit(1);
	}
}

//**********************************************************************************************************************

double AmovaCommand::calcSSTotal(map<string, vector<int> >& groupSampleMap) {
	try {
		
		vector<int> indices;
		map<string, vector<int> >::iterator it;
		for(it=groupSampleMap.begin();it!=groupSampleMap.end();it++){
			indices.insert(indices.end(), it->second.begin(), it->second.end());			
		}
		sort(indices.begin(), indices.end());
			
 		int numIndices =indices.size();
		double ssTotal = 0.0;
		
		for(int i=1;i<numIndices;i++){
			int row = indices[i];
			
			for(int j=0;j<i;j++){
				ssTotal += distanceMatrix[row][indices[j]];
			}
		}
		ssTotal /= numIndices;
			
		return ssTotal;
	}
	catch(exception& e) {
		m->errorOut(e, "AmovaCommand", "calcSSTotal");
		exit(1);
	}
}

//**********************************************************************************************************************

double AmovaCommand::calcSSWithin(map<string, vector<int> >& groupSampleMap) {
	try {

		double ssWithin = 0.0;
		
		map<string, vector<int> >::iterator it;
		for(it=groupSampleMap.begin();it!=groupSampleMap.end();it++){
			
			double withinGroup = 0;
			
			vector<int> samples = it->second;
			
			for(int i=0;i<samples.size();i++){
				int row = samples[i];

				for(int j=0;j<samples.size();j++){
					int col = samples[j];

					if(col < row){
						withinGroup += distanceMatrix[row][col];
					}
					
				}
			}

			ssWithin += withinGroup / samples.size();
		}

		return ssWithin;
	}
	catch(exception& e) {
		m->errorOut(e, "AmovaCommand", "calcSSWithin");
		exit(1);
	}
}

//**********************************************************************************************************************
