
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
vector<string> PCOACommand::getValidParameters(){	
	try {
		string Array[] =  {"phylip", "metric","outputdir","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
PCOACommand::PCOACommand(){	
	try {
		abort = true; calledHelp = true; 
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
vector<string> PCOACommand::getRequiredParameters(){	
	try {
		string Array[] =  {"phylip"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> PCOACommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************

PCOACommand::PCOACommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"phylip","metric","outputdir", "inputdir"};
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
			else if (phylipfile == "not found") { phylipfile = ""; abort = true; }	
			else {	filename = phylipfile;  }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(phylipfile); //if user entered a file with a path then preserve it	
			}
			
			//error checking on files	
			if (phylipfile == "")	{ m->mothurOut("You must provide a distance file before running the pcoa command."); m->mothurOutEndLine(); abort = true; }		
		
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
void PCOACommand::help(){
	try {
	
		m->mothurOut("The pcoa command parameters are phylip and metric"); m->mothurOutEndLine();
		m->mothurOut("The phylip parameter allows you to enter your distance file."); m->mothurOutEndLine();
		m->mothurOut("The metric parameter allows indicate you if would like the pearson correlation coefficient calculated. Default=True"); m->mothurOutEndLine();
		m->mothurOut("Example pcoa(phylip=yourDistanceFile).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. phylip), '=' and parameters (i.e.yourDistanceFile).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "help");
		exit(1);
	}
}
//**********************************************************************************************************************
PCOACommand::~PCOACommand(){}
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
		vector<vector<double> > copy_G;
				
		m->mothurOut("\nProcessing...\n\n");
		
		for(int count=0;count<2;count++){
			recenter(offset, D, G);					if (m->control_pressed) { return 0; }
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
				
				m->mothurOut("Pearson's coefficient using " + toString(i) + " axis: " + toString(corr)); m->mothurOutEndLine();
				
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

void PCOACommand::recenter(double offset, vector<vector<double> > D, vector<vector<double> >& G){
	try {
		int rank = D.size();
		
		vector<vector<double> > A(rank);
		vector<vector<double> > C(rank);
		for(int i=0;i<rank;i++){
			A[i].resize(rank);
			C[i].resize(rank);
		}
		
		double scale = -1.0000 / (double) rank;
		
		for(int i=0;i<rank;i++){
			A[i][i] = 0.0000;
			C[i][i] = 1.0000 + scale;
			for(int j=i+1;j<rank;j++){
				A[i][j] = A[j][i] = -0.5 * D[i][j] * D[i][j] + offset;
				C[i][j] = C[j][i] = scale;
			}
		}
		
		A = linearCalc.matrix_mult(C,A);
		G = linearCalc.matrix_mult(A,C);
	}
	catch(exception& e) {
		m->errorOut(e, "PCOACommand", "recenter");
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

