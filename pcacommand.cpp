/*
 *  pcacommand.cpp
 *  mothur
 *
 *  Created by westcott on 1/7/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "pcacommand.h"
#include "inputdata.h"

//**********************************************************************************************************************
vector<string> PCACommand::getValidParameters(){	
	try {
		string Array[] =  {"label", "groups","metric","outputdir","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PCACommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
PCACommand::PCACommand(){	
	try {
		abort = true; calledHelp = true; 
		vector<string> tempOutNames;
		outputTypes["pca"] = tempOutNames;
		outputTypes["loadings"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "PCACommand", "PCACommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> PCACommand::getRequiredParameters(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PCACommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> PCACommand::getRequiredFiles(){	
	try {
		string Array[] =  {"shared","relabund","or"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PCACommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************

PCACommand::PCACommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		globaldata = GlobalData::getInstance();
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"label","groups","metric","outputdir", "inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser. getParameters();
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		if (inputDir == "not found"){	inputDir = "";		}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["pca"] = tempOutNames;
			outputTypes["loadings"] = tempOutNames;
			
			//make sure the user has already run the read.otu command
			if ((globaldata->getSharedFile() == "") && (globaldata->getRelAbundFile() == "")) {
				m->mothurOut("You must read a list and a group, shared or relabund file before you can use the pca command."); m->mothurOutEndLine(); abort = true; 
			}
			
			if (globaldata->getSharedFile() != "")		{ mode = "shared"; inputFile = globaldata->getSharedFile();		}
			if (globaldata->getRelAbundFile() != "")	{ mode = "relabund"; inputFile = globaldata->getRelAbundFile(); }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(inputFile); //if user entered a file with a path then preserve it	
			}
						
			string temp = validParameter.validFile(parameters, "metric", false);	if (temp == "not found"){	temp = "T";				}
			metric = m->isTrue(temp); 
			
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; labels = globaldata->labels; if(labels.size() == 0) {  m->mothurOut("You did not provide a label, I will use the first label in your inputfile."); m->mothurOutEndLine(); } }
			else { m->splitAtDash(label, labels); }
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = "";  }
			else { m->splitAtDash(groups, Groups);	}			
			globaldata->Groups = Groups;			
			
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "PCACommand", "PCACommand");
		exit(1);
	}
}
//**********************************************************************************************************************
void PCACommand::help(){
	try {
		m->mothurOut("The pca command can only be run after a successful read.otu command of a shared or relabund file."); m->mothurOutEndLine();
		m->mothurOut("The pca command parameters are label, groups and metric. No parameters are required."); m->mothurOutEndLine();
		m->mothurOut("The label parameter is used to analyze specific labels in your input. Default is the first label in your shared or relabund file. Multiple labels may be separated by dashes.\n");
		m->mothurOut("The groups parameter allows you to specify which groups you would like analyzed. Groupnames are separated by dashes.\n");
		m->mothurOut("The metric parameter allows indicate you if would like the pearson correlation coefficient calculated. Default=True"); m->mothurOutEndLine();
		m->mothurOut("Example pca(groups=yourGroups).\n");
		m->mothurOut("Example pca(groups=A-B-C).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "PCACommand", "help");
		exit(1);
	}
}
//**********************************************************************************************************************
PCACommand::~PCACommand(){}
//**********************************************************************************************************************
int PCACommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::showpoint);
		cerr.setf(ios::fixed, ios::floatfield);
		cerr.setf(ios::showpoint);
		
		//get first line of shared file
		vector< vector<double> > matrix;
		InputData* input;
		if (mode == "shared")			{  
			input = new InputData(inputFile, "sharedfile");
		}else if (mode == "relabund")	{ 
			input = new InputData(inputFile, "relabund");
		}else {  m->mothurOut("[ERROR]: filetype not recognized."); m->mothurOutEndLine();  return 0; }
		
		vector<SharedRAbundFloatVector*> lookupFloat = input->getSharedRAbundFloatVectors();
		string lastLabel = lookupFloat[0]->getLabel();
			
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//if the user gave no labels, then use the first one read
		if (labels.size() == 0) { 
			label = lastLabel;  
			
			process(lookupFloat);
		}
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookupFloat[0] != NULL) && (userLabels.size() != 0)) {
			
			if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } delete input; for (int i = 0; i < lookupFloat.size(); i++) {  delete lookupFloat[i];  }  lookupFloat.clear(); return 0;  }
			
			if(labels.count(lookupFloat[0]->getLabel()) == 1){
				processedLabels.insert(lookupFloat[0]->getLabel());
				userLabels.erase(lookupFloat[0]->getLabel());
				
				process(lookupFloat);
			}
			
			if ((m->anyLabelsToProcess(lookupFloat[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookupFloat[0]->getLabel();
				
				for (int i = 0; i < lookupFloat.size(); i++) {  delete lookupFloat[i];  }  lookupFloat.clear();
				lookupFloat = input->getSharedRAbundFloatVectors(lastLabel);
				
				process(lookupFloat);
				
				processedLabels.insert(lookupFloat[0]->getLabel());
				userLabels.erase(lookupFloat[0]->getLabel());
				
				//restore real lastlabel to save below
				lookupFloat[0]->setLabel(saveLabel);
			}
			
			lastLabel = lookupFloat[0]->getLabel();			
			
			//get next line to process
			//prevent memory leak
			for (int i = 0; i < lookupFloat.size(); i++) {  delete lookupFloat[i];  } lookupFloat.clear();
			lookupFloat = input->getSharedRAbundFloatVectors();
		}
		
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } delete input; for (int i = 0; i < lookupFloat.size(); i++) {  delete lookupFloat[i];  } lookupFloat.clear(); return 0;  }
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			m->mothurOut("Your file does not include the label " + *it); 
			if (processedLabels.count(lastLabel) != 1) {
				m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
				needToRun = true;
			}else {
				m->mothurOut(". Please refer to " + lastLabel + "."); m->mothurOutEndLine();
			}
		}
		
		//run last label if you need to
		if (needToRun == true)  {
			for (int i = 0; i < lookupFloat.size(); i++) {  if (lookupFloat[i] != NULL) {	delete lookupFloat[i];	} }  lookupFloat.clear();
			lookupFloat = input->getSharedRAbundFloatVectors(lastLabel);
			
			process(lookupFloat);
			
			for (int i = 0; i < lookupFloat.size(); i++) {  if (lookupFloat[i] != NULL) {	delete lookupFloat[i];	} } lookupFloat.clear();
		}	
		
		for (int i = 0; i < lookupFloat.size(); i++) {  if (lookupFloat[i] != NULL) {	delete lookupFloat[i];	} } lookupFloat.clear();
		delete input;
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "PCACommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
vector< vector<double> > PCACommand::createMatrix(vector<SharedRAbundFloatVector*> lookupFloat){
	try {
		vector< vector<double> > matrix; matrix.resize(lookupFloat.size());
		
		//fill matrix with shared files relative abundances
		for (int i = 0; i < lookupFloat.size(); i++) {
			for (int j = 0; j < lookupFloat[i]->getNumBins(); j++) {
				matrix[i].push_back(lookupFloat[i]->getAbundance(j));
			}
		}
		
		vector< vector<double> > transposeMatrix; transposeMatrix.resize(matrix[0].size());
		for (int i = 0; i < transposeMatrix.size(); i++) {
			for (int j = 0; j < matrix.size(); j++) {
				transposeMatrix[i].push_back(matrix[j][i]);
			}
		}
		
		matrix = linearCalc.matrix_mult(matrix, transposeMatrix);
		
		return matrix;
	}
	catch(exception& e) {
		m->errorOut(e, "PCACommand", "createMatrix");	
		exit(1);
	}
}
//**********************************************************************************************************************
int PCACommand::process(vector<SharedRAbundFloatVector*>& lookupFloat){
	try {
		m->mothurOut("\nProcessing " + lookupFloat[0]->getLabel()); m->mothurOutEndLine();
		
		vector< vector<double> > matrix; matrix.resize(lookupFloat.size());
		
		//fill matrix with shared files relative abundances
		for (int i = 0; i < lookupFloat.size(); i++) {
			for (int j = 0; j < lookupFloat[i]->getNumBins(); j++) {
				matrix[i].push_back(lookupFloat[i]->getAbundance(j));
			}
		}
		
		vector< vector<double> > transposeMatrix; transposeMatrix.resize(matrix[0].size());
		for (int i = 0; i < transposeMatrix.size(); i++) {
			for (int j = 0; j < matrix.size(); j++) {
				transposeMatrix[i].push_back(matrix[j][i]);
			}
		}
		
		matrix = linearCalc.matrix_mult(matrix, transposeMatrix);		
		
		double offset = 0.0000;
		vector<double> d;
		vector<double> e;
		vector<vector<double> > G = matrix;
		vector<vector<double> > copy_G;
			
		for(int count=0;count<2;count++){
			linearCalc.tred2(G, d, e);				if (m->control_pressed) { return 0; }
			linearCalc.qtli(d, e, G);				if (m->control_pressed) { return 0; }
			offset = d[d.size()-1];
			if(offset > 0.0) break;
		} 
		
		if (m->control_pressed) { return 0; }
		
		string fbase = outputDir + m->getRootName(m->getSimpleName(inputFile));
		string outputFileName = fbase + lookupFloat[0]->getLabel();
		output(outputFileName, globaldata->Groups, G, d);
		
		if (metric) {   
			
			for (int i = 1; i < 4; i++) {
				
				vector< vector<double> > EuclidDists = linearCalc.calculateEuclidianDistance(G, i); //G is the pcoa file
				
				if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } return 0; }
				
				double corr = linearCalc.calcPearson(EuclidDists, matrix); //G is the pcoa file, D is the users distance matrix
				
				m->mothurOut("Pearson's coefficient using " + toString(i) + " axis: " + toString(corr)); m->mothurOutEndLine();
				
				if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } return 0; }
			}
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "PCACommand", "process");	
		exit(1);
	}
}
/*********************************************************************************************************************************/

void PCACommand::output(string fnameRoot, vector<string> name_list, vector<vector<double> >& G, vector<double> d) {
	try {
		int rank = name_list.size();
		double dsum = 0.0000;
		for(int i=0;i<rank;i++){
			dsum += d[i];
			for(int j=0;j<rank;j++){
				if(d[j] >= 0)	{	G[i][j] *= pow(d[j],0.5);	}
				else			{	G[i][j] = 0.00000;			}
			}
		}
		
		ofstream pcaData((fnameRoot+".pca.axes").c_str(), ios::trunc);
		pcaData.setf(ios::fixed, ios::floatfield);
		pcaData.setf(ios::showpoint);	
		outputNames.push_back(fnameRoot+".pca.axes");
		outputTypes["pca"].push_back(fnameRoot+".pca.axes");
		
		ofstream pcaLoadings((fnameRoot+".pca.loadings").c_str(), ios::trunc);
		pcaLoadings.setf(ios::fixed, ios::floatfield);
		pcaLoadings.setf(ios::showpoint);
		outputNames.push_back(fnameRoot+".pca.loadings");
		outputTypes["loadings"].push_back(fnameRoot+".pca.loadings");	
		
		pcaLoadings << "axis\tloading\n";
		for(int i=0;i<rank;i++){
			pcaLoadings << i+1 << '\t' << d[i] * 100.0 / dsum << endl;
		}
		
		pcaData << "group";
		for(int i=0;i<rank;i++){
			pcaData << '\t' << "axis" << i+1;
		}
		pcaData << endl;
		
		for(int i=0;i<rank;i++){
			pcaData << name_list[i] << '\t';
			for(int j=0;j<rank;j++){
				pcaData << G[i][j] << '\t';
			}
			pcaData << endl;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "PCACommand", "output");
		exit(1);
	}
}
/*********************************************************************************************************************************/


