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
        
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["amova"] = tempOutNames;
	
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
		helpString += "Referenced: Anderson MJ (2001). A new method for non-parametric multivariate analysis of variance. Austral Ecol 26: 32-46.\n";
		helpString += "The amova command outputs a .amova file.\n";
		helpString += "The amova command parameters are " + getCommandParameters() + ". The phylip and design parameters are required, unless you have valid current files.\n";
		helpString += "The design parameter allows you to assign your samples to groups when you are running amova. It is required.\n";
        helpString += "The sets parameter allows you to specify which of the sets in your designfile you would like to analyze. The set names are separated by dashes. The default is all sets in the design file.\n";
		helpString += "The iters parameter allows you to set number of randomization for the P value.  The default is 1000.\n";
		helpString += "The amova command should be in the following format: amova(phylip=file.dist, design=file.design).\n";
		
        getCommonQuestions();
        
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "AmovaCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string AmovaCommand::getCommonQuestions(){
    try {
        vector<string> questions, issues, qanswers, ianswers, howtos, hanswers;
        
        string issue = "...XXX is not in your design file, please correct."; issues.push_back(issue);
        string ianswer = "\tMothur expects the design file to be 2 column with a header line. The first column should contain the names of the samples in the distance matrix. The second column should contain the treatment each sample is assigned to. \n"; ianswers.push_back(ianswer);
        
        string commonQuestions = util.getFormattedHelp(questions, qanswers, issues, ianswers, howtos, hanswers);
        
        return commonQuestions;
    }
    catch(exception& e) {
        m->errorOut(e, "AmovaCommand", "getCommonQuestions");
        exit(1);
    }
}
//**********************************************************************************************************************
string AmovaCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "amova") {  pattern = "[filename],amova"; } //makes file like: amazon.align
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "AmovaCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
AmovaCommand::AmovaCommand(string option) {
	try {
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			//if the user changes the output directory command factory will send this info to us in the output parameter
            ValidParameters validParameter;
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";	}
			
			phylipFileName = validParameter.validFile(parameters, "phylip");
			if (phylipFileName == "not open") { phylipFileName = ""; abort = true; }
			else if (phylipFileName == "not found") { 
				//if there is a current phylip file, use it
				phylipFileName = current->getPhylipFile(); 
				if (phylipFileName != "") { m->mothurOut("Using " + phylipFileName + " as input file for the phylip parameter.\n");  }
				else { 	m->mothurOut("You have no current phylip file and the phylip parameter is required.\n");  abort = true; }
			}else { current->setPhylipFile(phylipFileName); }
			
			//check for required parameters
			designFileName = validParameter.validFile(parameters, "design");
			if (designFileName == "not open") { designFileName = ""; abort = true; }
			else if (designFileName == "not found") {
				//if there is a current design file, use it
				designFileName = current->getDesignFile(); 
				if (designFileName != "") { m->mothurOut("Using " + designFileName + " as input file for the design parameter.\n");  }
				else { 	m->mothurOut("You have no current design file and the design parameter is required.\n");  abort = true; }
			}else { current->setDesignFile(designFileName); }	

			string temp = validParameter.valid(parameters, "iters");
			if (temp == "not found") { temp = "1000"; }
			util.mothurConvert(temp, iters); 
			
			temp = validParameter.valid(parameters, "alpha");
			if (temp == "not found") { temp = "0.05"; }
			util.mothurConvert(temp, experimentwiseAlpha); 
            
            string sets = validParameter.valid(parameters, "sets");			
			if (sets == "not found") { sets = ""; }
			else { 
				util.splitAtDash(sets, Sets);
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
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		//read design file
		designMap = new DesignMap(designFileName);

		if (outputDir == "") { outputDir = util.hasPath(phylipFileName); }
						
		//read in distance matrix and square it
		ReadPhylipVector readMatrix(phylipFileName);
		vector<string> sampleNames = readMatrix.read(distanceMatrix);
        
        if (Sets.size() != 0) { //user selected sets, so we want to remove the samples not in those sets
            for(int i=0;i<distanceMatrix.size();i++){
                
                if (m->getControl_pressed()) { delete designMap; return 0; }
                
                string group = designMap->get(sampleNames[i]);
                
                if (group == "not found") {
                    m->mothurOut("[ERROR]: " + sampleNames[i] + " is not in your design file, please correct.\n");  m->setControl_pressed(true);
                }else if (!util.inUsersGroups(group, Sets)){  //not in set we want remove it
                    //remove from all other rows
                    for(int j=0;j<distanceMatrix.size();j++){ distanceMatrix[j].erase(distanceMatrix[j].begin()+i); }
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
				m->mothurOut("[ERROR]: " + sampleNames[i] + " is not in your design file, please correct.\n");  m->setControl_pressed(true);
			}else { origGroupSampleMap[group].push_back(i); }
			
		}
		int numGroups = origGroupSampleMap.size();
		
		if (m->getControl_pressed()) { delete designMap; return 0; }
		
		//create a new filename
		ofstream AMOVAFile;
        map<string, string> variables; variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(phylipFileName));
		string AMOVAFileName = getOutputFileName("amova", variables);	
        
		util.openOutputFile(AMOVAFileName, AMOVAFile);
		outputNames.push_back(AMOVAFileName); outputTypes["amova"].push_back(AMOVAFileName);
		
		double fullANOVAPValue = runAMOVA(AMOVAFile, origGroupSampleMap, experimentwiseAlpha);
		if(fullANOVAPValue <= experimentwiseAlpha && numGroups > 2){
			
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
	 
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
		
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
		
        vector<string> sampleNames;
		for(it = groupSampleMap.begin();it!=groupSampleMap.end();it++){ sampleNames.push_back(it->first); }
        string output = util.getStringFromVector(sampleNames, "-");
		
		AMOVAFile << output << "\tAmong\tWithin\tTotal" << endl;
		m->mothurOut(output + "\tAmong\tWithin\tTotal\n");
		
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
		
		util.mothurRandomShuffle(sampleIndices);
		
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
