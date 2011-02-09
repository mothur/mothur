/*
 *  mantelcommand.cpp
 *  mothur
 *
 *  Created by westcott on 2/9/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "mantelcommand.h"
#include "readphylipvector.h"

//**********************************************************************************************************************
vector<string> MantelCommand::getValidParameters(){	
	try {
		string Array[] =  {"phylip1","phylip2","method","iters","outputdir","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MantelCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> MantelCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"phylip1", "phylip2"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MantelCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
MantelCommand::MantelCommand(){	
	try {
		abort = true; calledHelp = true; 
		vector<string> tempOutNames;
		outputTypes["mantel"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "MantelCommand", "MantelCommand");
		exit(1);
	}
}

//**********************************************************************************************************************
vector<string> MantelCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MantelCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************
MantelCommand::MantelCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"phylip1","phylip2","method","iters","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			vector<string> tempOutNames;
			outputTypes["mantel"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("phylip1");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["phylip1"] = inputDir + it->second;		}
				}
				
				it = parameters.find("phylip2");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["phylip2"] = inputDir + it->second;		}
				}
			}
			
			
			//check for required parameters
			phylipfile1 = validParameter.validFile(parameters, "phylip1", true);
			if (phylipfile1 == "not open") { phylipfile1 = ""; abort = true; }
			else if (phylipfile1 == "not found") { phylipfile1 = ""; m->mothurOut("phylip1 is a required parameter for the mantel command."); m->mothurOutEndLine(); abort = true;  }	
			
			phylipfile2 = validParameter.validFile(parameters, "phylip2", true);
			if (phylipfile2 == "not open") { phylipfile2 = ""; abort = true; }
			else if (phylipfile2 == "not found") { phylipfile2 = ""; m->mothurOut("phylip2 is a required parameter for the mantel command."); m->mothurOutEndLine(); abort = true;  }	
			
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(phylipfile1);	}
			
			method = validParameter.validFile(parameters, "method", false);		if (method == "not found"){	method = "pearson";		}
			
			string temp = validParameter.validFile(parameters, "iters", false);			if (temp == "not found") { temp = "1000"; }
			convert(temp, iters);
			
			if ((method != "pearson") && (method != "spearman") && (method != "kendall")) { m->mothurOut(method + " is not a valid method. Valid methods are pearson, spearman, and kendall."); m->mothurOutEndLine(); abort = true; }
		}
	}
	catch(exception& e) {
		m->errorOut(e, "MantelCommand", "MantelCommand");		
		exit(1);
	}
}
//**********************************************************************************************************************

void MantelCommand::help(){
	try {
		m->mothurOut("Sokal, R. R., & Rohlf, F. J. (1995). Biometry, 3rd edn. New York: Freeman.\n");
		m->mothurOut("The mantel command reads two distance matrices and calculates the mantel correlation coefficient.\n");
		m->mothurOut("The mantel command parameters are phylip1, phylip2, iters and method.  The phylip1 and phylip2 parameters are required.  Matrices must be the same size and contain the same names.\n");
		m->mothurOut("The method parameter allows you to select what method you would like to use. Options are pearson, spearman and kendall. Default=pearson.\n");
		m->mothurOut("The iters parameter allows you to set number of randomization for the P value.  The default is 1000. \n");
		m->mothurOut("The mantel command should be in the following format: mantel(phylip1=veg.dist, phylip2=env.dist).\n");
		m->mothurOut("The mantel command outputs a .mantel file.\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. phylip1), '=' and parameters (i.e. veg.dist).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "MantelCommand", "help");	
		exit(1);
	}
}

//**********************************************************************************************************************

MantelCommand::~MantelCommand(){}

//**********************************************************************************************************************

int MantelCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		/***************************************************/
		//    reading distance files					   //
		/***************************************************/
		
		//read phylip1
		ReadPhylipVector readMatrix(phylipfile1);
		vector< vector<double> > matrix1;
		vector<string> names1 = readMatrix.read(matrix1);
		
		if (m->control_pressed) { return 0; }
		
		//read phylip2
		ReadPhylipVector readMatrix2(phylipfile2);
		vector< vector<double> > matrix2;
		vector<string> names2 = readMatrix2.read(matrix2);
		
		if (m->control_pressed) { return 0; }
		
		//make sure matrix2 and matrix1 are in the same order
		if (names1 == names2) { //then everything is in same order and same size
		}else if (names1.size() != names2.size()) { //wrong size no need to order, abort
			m->mothurOut("[ERROR]: distance matrices are not the same size, aborting."); m->mothurOutEndLine();
			m->control_pressed = true;
		}else { //sizes are the same, but either the names are different or they are in different order
			m->mothurOut("[WARNING]: Names do not match between distance files. Comparing based on order in files."); m->mothurOutEndLine();
		}	
		
		if (m->control_pressed) { return 0; }
		
		/***************************************************/
		//    calculating mantel and signifigance		   //
		/***************************************************/
		
		//calc mantel coefficient
		LinearAlgebra linear;
		double mantel = 0.0;
		if (method == "pearson")		{  mantel = linear.calcPearson(matrix1, matrix2);	}
		else if (method == "spearman")	{  mantel = linear.calcSpearman(matrix1, matrix2);	}
		else if (method == "kendall")	{  mantel = linear.calcKendall(matrix1, matrix2);	}
		
		
		//calc signifigance
		int count = 0;
		for (int i = 0; i < iters; i++) {
			
			if (m->control_pressed) { return 0; }
			
			//randomize matrix2
			vector< vector<double> > matrix2Copy = matrix2;
			random_shuffle(matrix2Copy.begin(), matrix2Copy.end());
		
			//calc random mantel
			double randomMantel = 0.0;
			if (method == "pearson")		{  randomMantel = linear.calcPearson(matrix1, matrix2Copy);	}
			else if (method == "spearman")	{  randomMantel = linear.calcSpearman(matrix1, matrix2Copy);	}
			else if (method == "kendall")	{  randomMantel = linear.calcKendall(matrix1, matrix2Copy);	}
			
			if (randomMantel >= mantel) { count++; }
		}
		
		double pValue = count / (float) iters;
		
		if (m->control_pressed) { return 0; }
		
		string outputFile = outputDir + m->getRootName(m->getSimpleName(phylipfile1)) + "mantel";
		outputNames.push_back(outputFile); outputTypes["mantel"].push_back(outputFile);
		ofstream out;
		
		m->openOutputFile(outputFile, out);
		
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		cout.setf(ios::fixed, ios::floatfield); cout.setf(ios::showpoint);
		
		out << "Mantel\tpValue" << endl;
		out << mantel << '\t' << pValue << endl;
		
		out.close();
	
		cout << "\nmantel = " << mantel << "\tpValue = " << pValue << endl;
		m->mothurOutJustToLog("\nmantel = " + toString(mantel) + "\tpValue = " + toString(pValue) + "\n"); 
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();		
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "MantelCommand", "execute");	
		exit(1);
	}
}

//**********************************************************************************************************************


