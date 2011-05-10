
/*
 *  pcacommand.cpp
 *  Mothur
 *
 *  Created by westcott on 1/4/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "pcoacommand.h"
#include "readphylipvector.h"


//**********************************************************************************************************************
vector<string> PCOACommand::setParameters(){	
	try {
		CommandParameter pphylip("phylip", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pphylip);
		CommandParameter pmetric("metric", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(pmetric);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string PCOACommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The pcoa command parameters are phylip and metric"; 
		helpString += "The phylip parameter allows you to enter your distance file.";
		helpString += "The metric parameter allows indicate you if would like the pearson correlation coefficient calculated. Default=True"; 
		helpString += "Example pcoa(phylip=yourDistanceFile).\n";
		helpString += "Note: No spaces between parameter labels (i.e. phylip), '=' and parameters (i.e.yourDistanceFile).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "getHelpString");
		exit(1);
	}
}


//**********************************************************************************************************************
PCOACommand::PCOACommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["pcoa"] = tempOutNames;
		outputTypes["loadings"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "PCOACommand");
		exit(1);
	}
}
//**********************************************************************************************************************

PCOACommand::PCOACommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser. getParameters();
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("phylip");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["phylip"] = inputDir + it->second;		}
				}
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["pcoa"] = tempOutNames;
			outputTypes["loadings"] = tempOutNames;
			
			//required parameters
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not open") { abort = true; }
			else if (phylipfile == "not found") { 			
				//if there is a current phylip file, use it
				phylipfile = m->getPhylipFile(); 
				if (phylipfile != "") { m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current phylip file and the phylip parameter is required."); m->mothurOutEndLine(); abort = true; }
			}	
			
			filename = phylipfile;  
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(phylipfile); //if user entered a file with a path then preserve it	
			}
			
			string temp = validParameter.validFile(parameters, "metric", false);	if (temp == "not found"){	temp = "T";				}
			metric = m->isTrue(temp); 
		}

	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "PCOACommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int PCOACommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		cout.setf(ios::fixed, ios::floatfield);
		cout.setf(ios::showpoint);
		cerr.setf(ios::fixed, ios::floatfield);
		cerr.setf(ios::showpoint);
		
		vector<string> names;
		vector<vector<double> > D;
	
		fbase = outputDir + m->getRootName(m->getSimpleName(filename));
		
		ReadPhylipVector readFile(filename);
		names = readFile.read(D);
		
		if (m->control_pressed) { return 0; }
   	
		double offset = 0.0000;
		vector<double> d;
		vector<double> e;
		vector<vector<double> > G = D;
		//vector<vector<double> > copy_G;
				
		m->mothurOut("\nProcessing...\n");
		
		for(int count=0;count<2;count++){
			linearCalc.recenter(offset, D, G);		if (m->control_pressed) { return 0; }
			linearCalc.tred2(G, d, e);				if (m->control_pressed) { return 0; }
			linearCalc.qtli(d, e, G);				if (m->control_pressed) { return 0; }
			offset = d[d.size()-1];
			if(offset > 0.0) break;
		} 
		
		if (m->control_pressed) { return 0; }
		
		output(fbase, names, G, d);
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } return 0; }
		
		if (metric) {   
			
			for (int i = 1; i < 4; i++) {
							
				vector< vector<double> > EuclidDists = linearCalc.calculateEuclidianDistance(G, i); //G is the pcoa file
				
				if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } return 0; }
				
				double corr = linearCalc.calcPearson(EuclidDists, D); //G is the pcoa file, D is the users distance matrix
				
				m->mothurOut("Rsq " + toString(i) + " axis: " + toString(corr * corr)); m->mothurOutEndLine();
				
				if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } return 0; }
			}
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "execute");
		exit(1);
	}
}
/*********************************************************************************************************************************/

void PCOACommand::get_comment(istream& f, char begin, char end){
	try {
		char d=f.get();
		while(d != end){	d = f.get();	}
		d = f.peek();
	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "get_comment");
		exit(1);
	}
}	
/*********************************************************************************************************************************/

void PCOACommand::output(string fnameRoot, vector<string> name_list, vector<vector<double> >& G, vector<double> d) {
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
		
		ofstream pcaData((fnameRoot+"pcoa.axes").c_str(), ios::trunc);
		pcaData.setf(ios::fixed, ios::floatfield);
		pcaData.setf(ios::showpoint);	
		outputNames.push_back(fnameRoot+"pcoa.axes");
		outputTypes["pcoa"].push_back(fnameRoot+"pcoa.axes");
		
		ofstream pcaLoadings((fnameRoot+"pcoa.loadings").c_str(), ios::trunc);
		pcaLoadings.setf(ios::fixed, ios::floatfield);
		pcaLoadings.setf(ios::showpoint);
		outputNames.push_back(fnameRoot+"pcoa.loadings");
		outputTypes["loadings"].push_back(fnameRoot+"pcoa.loadings");	
		
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
		m->errorOut(e, "PCOACommand", "output");
		exit(1);
	}
}

/*********************************************************************************************************************************/

