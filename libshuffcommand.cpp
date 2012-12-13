/*
 *  libshuffcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class is designed to implement an integral form of the Cramer-von Mises statistic.
	you may refer to the "Integration of Microbial Ecology and Statistics: A Test To Compare Gene Libraries" 
	paper in Applied and Environmental Microbiology, Sept. 2004, p. 5485-5492 0099-2240/04/$8.00+0  
	DOI: 10.1128/AEM.70.9.5485-5492.2004 Copyright 2004 American Society for Microbiology for more information. */


#include "libshuffcommand.h"
#include "libshuff.h"
#include "slibshuff.h"
#include "dlibshuff.h"



//**********************************************************************************************************************
vector<string> LibShuffCommand::setParameters(){	
	try {
		CommandParameter pphylip("phylip", "InputTypes", "", "", "none", "none", "none","coverage-libshuffsummary",false,true,true); parameters.push_back(pphylip);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(pgroup);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter piters("iters", "Number", "", "10000", "", "", "","",false,false); parameters.push_back(piters);
		CommandParameter pstep("step", "Number", "", "0.01", "", "", "","",false,false); parameters.push_back(pstep);
		CommandParameter pcutoff("cutoff", "Number", "", "1.0", "", "", "","",false,false); parameters.push_back(pcutoff);
		CommandParameter pform("form", "Multiple", "discrete-integral", "integral", "", "", "","",false,false); parameters.push_back(pform);
		CommandParameter psim("sim", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(psim);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "LibShuffCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string LibShuffCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The libshuff command parameters are phylip, group, sim, groups, iters, step, form and cutoff.  phylip and group parameters are required, unless you have valid current files.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 2 valid groups.\n";
		helpString += "The group names are separated by dashes.  The iters parameter allows you to specify how many random matrices you would like compared to your matrix.\n";
		helpString += "The step parameter allows you to specify change in distance you would like between each output if you are using the discrete form.\n";
		helpString += "The form parameter allows you to specify if you would like to analyze your matrix using the discrete or integral form. Your options are integral or discrete.\n";
		helpString += "The libshuff command should be in the following format: libshuff(groups=yourGroups, iters=yourIters, cutOff=yourCutOff, form=yourForm, step=yourStep).\n";
		helpString += "Example libshuff(groups=A-B-C, iters=500, form=discrete, step=0.01, cutOff=2.0).\n";
		helpString += "The default value for groups is all the groups in your groupfile, iters is 10000, cutoff is 1.0, form is integral and step is 0.01.\n";
		helpString += "The libshuff command output two files: .coverage and .slsummary their descriptions are in the manual.\n";
		helpString += "Note: No spaces between parameter labels (i.e. iters), '=' and parameters (i.e.yourIters).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "LibShuffCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string LibShuffCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "coverage") {  pattern = "[filename],libshuff.coverage"; } 
        else if (type == "libshuffsummary") {  pattern = "[filename],libshuff.summary"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "LibShuffCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
LibShuffCommand::LibShuffCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["coverage"] = tempOutNames;
		outputTypes["libshuffsummary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "LibShuffCommand", "LibShuffCommand");
		exit(1);
	}
}

//**********************************************************************************************************************
LibShuffCommand::LibShuffCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			map<string,string>::iterator it;
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["coverage"] = tempOutNames;
			outputTypes["libshuffsummary"] = tempOutNames;
			
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
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not open") { phylipfile = ""; abort = true; }
			else if (phylipfile == "not found") { 
				phylipfile = m->getPhylipFile(); 
				if (phylipfile != "") {  m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter."); m->mothurOutEndLine(); }
				else { 
					m->mothurOut("You must provide a phylip file."); m->mothurOutEndLine(); 
					abort = true;
				} 
			}else { m->setPhylipFile(phylipfile); }	
			
			//check for required parameters
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { groupfile = ""; abort = true; }
			else if (groupfile == "not found") { 
				groupfile = m->getGroupFile(); 
				if (groupfile != "") {  m->mothurOut("Using " + groupfile + " as input file for the group parameter."); m->mothurOutEndLine(); }
				else { 
					m->mothurOut("You must provide a group file."); m->mothurOutEndLine(); 
					abort = true;
				} 
			}else { m->setGroupFile(groupfile); }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(phylipfile); //if user entered a file with a path then preserve it	
			}
						
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; savegroups = groups; }
			else { 
				savegroups = groups;
				m->splitAtDash(groups, Groups);
				m->setGroups(Groups);
			}
				
			string temp;
			temp = validParameter.validFile(parameters, "iters", false);				if (temp == "not found") { temp = "10000"; }
			m->mothurConvert(temp, iters); 
			
			temp = validParameter.validFile(parameters, "cutoff", false);				if (temp == "not found") { temp = "1.0"; }
			m->mothurConvert(temp, cutOff); 
			
			temp = validParameter.validFile(parameters, "step", false);				if (temp == "not found") { temp = "0.01"; }
			m->mothurConvert(temp, step); 
			
			temp = validParameter.validFile(parameters, "sim", false);				if (temp == "not found") { temp = "F"; }
			sim = m->isTrue(temp); 
			
			userform = validParameter.validFile(parameters, "form", false);			if (userform == "not found") { userform = "integral"; }
			
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "LibShuffCommand", "LibShuffCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int LibShuffCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//read files
		groupMap = new GroupMap(groupfile);
		int error = groupMap->readMap();
		if (error == 1) { delete groupMap; return 0; }
		
		ifstream in;
		m->openInputFile(phylipfile, in);
		matrix = new FullMatrix(in, groupMap, sim); //reads the matrix file
		in.close();
		
		if (m->control_pressed) { delete groupMap; delete matrix; return 0; }
		
		//if files don't match...
		if (matrix->getNumSeqs() < groupMap->getNumSeqs()) {  
			m->mothurOut("Your distance file contains " + toString(matrix->getNumSeqs()) + " sequences, and your group file contains " + toString(groupMap->getNumSeqs()) + " sequences.");  m->mothurOutEndLine();				
			//create new group file
			if(outputDir == "") { outputDir += m->hasPath(groupfile); }
			
			string newGroupFile = outputDir + m->getRootName(m->getSimpleName(groupfile)) + "editted.groups";
			outputNames.push_back(newGroupFile);
			ofstream outGroups;
			m->openOutputFile(newGroupFile, outGroups);
			
			for (int i = 0; i < matrix->getNumSeqs(); i++) {
				if (m->control_pressed) { delete groupMap; delete matrix; outGroups.close(); m->mothurRemove(newGroupFile); return 0; }
				
				Names temp = matrix->getRowInfo(i);
				outGroups << temp.seqName << '\t' << temp.groupName << endl;
			}
			outGroups.close();
			
			m->mothurOut(newGroupFile + " is a new group file containing only the sequence that are in your distance file. I will read this file instead."); m->mothurOutEndLine();
			
			//read new groupfile
			delete groupMap; 
			groupfile = newGroupFile;
			
			groupMap = new GroupMap(groupfile);
			groupMap->readMap();
			
			if (m->control_pressed) { delete groupMap; delete matrix; m->mothurRemove(newGroupFile); return 0; }
		}
		
			
		setGroups();								//set the groups to be analyzed and sorts them
		
		if (numGroups < 2) { m->mothurOut("[ERROR]: libshuff requires at least 2 groups, you only have " + toString(numGroups) + ", aborting."); m->mothurOutEndLine(); m->control_pressed = true; }
		
		if (m->control_pressed) { delete groupMap; delete matrix; return 0; }
		
		/********************************************************************************************/
		//this is needed because when we read the matrix we sort it into groups in alphabetical order
		//the rest of the command and the classes used in this command assume specific order
		/********************************************************************************************/
		matrix->setGroups(groupMap->getNamesOfGroups());
		vector<int> sizes;
		for (int i = 0; i < (groupMap->getNamesOfGroups()).size(); i++) {   sizes.push_back(groupMap->getNumSeqs((groupMap->getNamesOfGroups())[i]));  }
		matrix->setSizes(sizes);
			
			
		if(userform == "discrete"){
			form = new DLibshuff(matrix, iters, step, cutOff);
		}
		else{
			form = new SLibshuff(matrix, iters, cutOff);
		}
	
		savedDXYValues = form->evaluateAll();
		savedMinValues = form->getSavedMins();
		
		if (m->control_pressed) {  delete form; m->clearGroups(); delete matrix; delete groupMap; return 0; }
	
		pValueCounts.resize(numGroups);
		for(int i=0;i<numGroups;i++){
			pValueCounts[i].assign(numGroups, 0);
		}
	
		if (m->control_pressed) {  outputTypes.clear(); delete form; m->clearGroups(); delete matrix; delete groupMap; return 0; }
				
		Progress* reading = new Progress();
		
		for(int i=0;i<numGroups-1;i++) {
			for(int j=i+1;j<numGroups;j++) {
				
				if (m->control_pressed) {  outputTypes.clear();  delete form; m->clearGroups(); delete matrix; delete groupMap; delete reading; return 0; }

				reading->newLine(groupNames[i]+'-'+groupNames[j], iters);
				int spoti = groupMap->groupIndex[groupNames[i]]; //neccessary in case user selects groups so you know where they are in the matrix
				int spotj = groupMap->groupIndex[groupNames[j]];
	
				for(int p=0;p<iters;p++) {	
					
					if (m->control_pressed) {  outputTypes.clear(); delete form; m->clearGroups(); delete matrix; delete groupMap; delete reading; return 0; }
					
					form->randomizeGroups(spoti,spotj); 
					if(form->evaluatePair(spoti,spotj) >= savedDXYValues[spoti][spotj])	{	pValueCounts[i][j]++;	}
					if(form->evaluatePair(spotj,spoti) >= savedDXYValues[spotj][spoti])	{	pValueCounts[j][i]++;	}
					
					if (m->control_pressed) {  outputTypes.clear(); delete form; m->clearGroups(); delete matrix; delete groupMap; delete reading; return 0; }
					
					reading->update(p);			
				}
				form->resetGroup(spoti);
				form->resetGroup(spotj);
			}
		}
		
		if (m->control_pressed) { outputTypes.clear();  delete form; m->clearGroups(); delete matrix; delete groupMap; delete reading; return 0; }
	
		reading->finish();
		delete reading;

		m->mothurOutEndLine();
		printSummaryFile();
		printCoverageFile();
				
		//clear out users groups
		m->clearGroups();
		delete form;
		
		delete matrix; delete groupMap;
		
		if (m->control_pressed) {  outputTypes.clear(); for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }

		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "LibShuffCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************

int LibShuffCommand::printCoverageFile() {
	try {

		ofstream outCov;
        map<string, string> variables; 
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(phylipfile));
		summaryFile = getOutputFileName("coverage", variables);
		m->openOutputFile(summaryFile, outCov);
		outputNames.push_back(summaryFile); outputTypes["coverage"].push_back(summaryFile);
		outCov.setf(ios::fixed, ios::floatfield); outCov.setf(ios::showpoint);
		//cout.setf(ios::fixed, ios::floatfield); cout.setf(ios::showpoint);
		
		map<double,vector<int> > allDistances;
		map<double,vector<int> >::iterator it;

		vector<vector<int> > indices(numGroups);
		int numIndices = numGroups * numGroups;
		
		int index = 0;
		for(int i=0;i<numGroups;i++){
			indices[i].assign(numGroups,0);
			for(int j=0;j<numGroups;j++){
				indices[i][j] = index++;
				
				int spoti = groupMap->groupIndex[groupNames[i]]; //neccessary in case user selects groups so you know where they are in the matrix
				int spotj = groupMap->groupIndex[groupNames[j]];
				
				for(int k=0;k<savedMinValues[spoti][spotj].size();k++){
					
					if(m->control_pressed)  { outCov.close(); return 0; }
					
					if(allDistances[savedMinValues[spoti][spotj][k]].size() != 0){
						allDistances[savedMinValues[spoti][spotj][k]][indices[i][j]]++;
					}
					else{
						allDistances[savedMinValues[spoti][spotj][k]].assign(numIndices, 0);
						allDistances[savedMinValues[spoti][spotj][k]][indices[i][j]] = 1;
					}
				}
			}
		}
		it=allDistances.begin();
		
		//cout << setprecision(8);

		vector<int> prevRow = it->second;
		it++;
		
		for(;it!=allDistances.end();it++){
			for(int i=0;i<it->second.size();i++){
				it->second[i] += prevRow[i];
			}
			prevRow = it->second;
		}
		
		vector<int> lastRow = allDistances.rbegin()->second;
		outCov << setprecision(8);
		
		outCov << "dist";
		for (int i = 0; i < numGroups; i++){
			outCov << '\t' << groupNames[i];
		}
		for (int i=0;i<numGroups;i++){
			for(int j=i+1;j<numGroups;j++){
				if(m->control_pressed)  { outCov.close(); return 0; }
				outCov << '\t' << groupNames[i] << '-' << groupNames[j] << '\t';
				outCov << groupNames[j] << '-' << groupNames[i];
			}
		}
		outCov << endl;
		
		for(it=allDistances.begin();it!=allDistances.end();it++){
			outCov << it->first << '\t';
			for(int i=0;i<numGroups;i++){
				outCov << it->second[indices[i][i]]/(float)lastRow[indices[i][i]] << '\t';
			}
			for(int i=0;i<numGroups;i++){
				for(int j=i+1;j<numGroups;j++){
					if(m->control_pressed)  { outCov.close(); return 0; }
					
					outCov << it->second[indices[i][j]]/(float)lastRow[indices[i][j]] << '\t';
					outCov << it->second[indices[j][i]]/(float)lastRow[indices[j][i]] << '\t';
				}
			}
			outCov << endl;
		}
		outCov.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "LibShuffCommand", "printCoverageFile");
		exit(1);
	}
} 

//**********************************************************************************************************************

int LibShuffCommand::printSummaryFile() {
	try {

		ofstream outSum;
        map<string, string> variables; 
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(phylipfile));
		summaryFile = getOutputFileName("libshuffsummary",variables);
		m->openOutputFile(summaryFile, outSum);
		outputNames.push_back(summaryFile); outputTypes["libshuffsummary"].push_back(summaryFile);

		outSum.setf(ios::fixed, ios::floatfield); outSum.setf(ios::showpoint);
		cout.setf(ios::fixed, ios::floatfield); cout.setf(ios::showpoint);
		
		cout << setw(20) << left << "Comparison" << '\t' << setprecision(8) << "dCXYScore" << '\t' << "Significance" << endl;
		m->mothurOutJustToLog("Comparison\tdCXYScore\tSignificance"); m->mothurOutEndLine();
		outSum << setw(20) << left << "Comparison" << '\t' << setprecision(8) << "dCXYScore" << '\t' << "Significance" << endl;
	
		int precision = (int)log10(iters);
		for(int i=0;i<numGroups;i++){
			for(int j=i+1;j<numGroups;j++){
				if(m->control_pressed)  { outSum.close(); return 0; }
				
				int spoti = groupMap->groupIndex[groupNames[i]]; //neccessary in case user selects groups so you know where they are in the matrix
				int spotj = groupMap->groupIndex[groupNames[j]];
				
				if(pValueCounts[i][j]){
					cout << setw(20) << left << groupNames[i]+'-'+groupNames[j] << '\t' << setprecision(8) << savedDXYValues[spoti][spotj] << '\t' << setprecision(precision) << pValueCounts[i][j]/(float)iters << endl;
					m->mothurOutJustToLog(groupNames[i]+"-"+groupNames[j] + "\t" + toString(savedDXYValues[spoti][spotj]) + "\t" + toString((pValueCounts[i][j]/(float)iters))); m->mothurOutEndLine();
					outSum << setw(20) << left << groupNames[i]+'-'+groupNames[j] << '\t' << setprecision(8) << savedDXYValues[spoti][spotj] << '\t' << setprecision(precision) << pValueCounts[i][j]/(float)iters << endl;
				}
				else{
					cout << setw(20) << left << groupNames[i]+'-'+groupNames[j] << '\t' << setprecision(8) << savedDXYValues[spoti][spotj] << '\t' << '<' <<setprecision(precision) << 1/(float)iters << endl;
					m->mothurOutJustToLog(groupNames[i]+"-"+groupNames[j] + "\t" + toString(savedDXYValues[spoti][spotj]) + "\t" + toString((1/(float)iters))); m->mothurOutEndLine();
					outSum << setw(20) << left << groupNames[i]+'-'+groupNames[j] << '\t' << setprecision(8) << savedDXYValues[spoti][spotj] << '\t' << '<' <<setprecision(precision) << 1/(float)iters << endl;
				}
				if(pValueCounts[j][i]){
					cout << setw(20) << left << groupNames[j]+'-'+groupNames[i] << '\t' << setprecision(8) << savedDXYValues[spotj][spoti] << '\t' << setprecision (precision) << pValueCounts[j][i]/(float)iters << endl;
					m->mothurOutJustToLog(groupNames[j]+"-"+groupNames[i] + "\t" + toString(savedDXYValues[spotj][spoti]) + "\t" + toString((pValueCounts[j][i]/(float)iters))); m->mothurOutEndLine();
					outSum << setw(20) << left << groupNames[j]+'-'+groupNames[i] << '\t' << setprecision(8) << savedDXYValues[spotj][spoti] << '\t' << setprecision (precision) << pValueCounts[j][i]/(float)iters << endl;
				}
				else{
					cout << setw(20) << left << groupNames[j]+'-'+groupNames[i] << '\t' << setprecision(8) << savedDXYValues[spotj][spoti] << '\t' << '<' <<setprecision (precision) << 1/(float)iters << endl;
					m->mothurOutJustToLog(groupNames[j]+"-"+groupNames[i] + "\t" + toString(savedDXYValues[spotj][spoti]) + "\t" + toString((1/(float)iters))); m->mothurOutEndLine();
					outSum << setw(20) << left << groupNames[j]+'-'+groupNames[i] << '\t' << setprecision(8) << savedDXYValues[spotj][spoti] << '\t' << '<' <<setprecision (precision) << 1/(float)iters << endl;
				}
			}
		}
		
		outSum.close();
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "LibShuffCommand", "printSummaryFile");
		exit(1);
	}
} 

//**********************************************************************************************************************

void LibShuffCommand::setGroups() {
	try {
		vector<string> myGroups = m->getGroups();
		//if the user has not entered specific groups to analyze then do them all
		if (m->getNumGroups() == 0) {
			numGroups = groupMap->getNumGroups();
			for (int i=0; i < numGroups; i++) { 
				myGroups.push_back((groupMap->getNamesOfGroups())[i]);
			}
		} else {
			if (savegroups != "all") {
				//check that groups are valid
				for (int i = 0; i < myGroups.size(); i++) {
					if (groupMap->isValidGroup(myGroups[i]) != true) {
						m->mothurOut(myGroups[i] + " is not a valid group, and will be disregarded."); m->mothurOutEndLine();
						// erase the invalid group from globaldata->Groups
						myGroups.erase(myGroups.begin()+i);
					}
				}
			
				//if the user only entered invalid groups
				if ((myGroups.size() == 0) || (myGroups.size() == 1)) { 
					numGroups = groupMap->getNumGroups();
					for (int i=0; i < numGroups; i++) { 
						myGroups.push_back((groupMap->getNamesOfGroups())[i]);
					}
					m->mothurOut("When using the groups parameter you must have at least 2 valid groups. I will run the command using all the groups in your groupfile."); m->mothurOutEndLine();
				} else { numGroups = myGroups.size(); }
			} else { //users wants all groups
				numGroups = groupMap->getNumGroups();
				myGroups.clear();
				for (int i=0; i < numGroups; i++) { 
					myGroups.push_back((groupMap->getNamesOfGroups())[i]);
				}
			}
		}

		//sort so labels match
		sort(myGroups.begin(), myGroups.end());
		
		//sort
		//sort(groupMap->namesOfGroups.begin(), groupMap->namesOfGroups.end());
		
		for (int i = 0; i < (groupMap->getNamesOfGroups()).size(); i++) {  groupMap->groupIndex[(groupMap->getNamesOfGroups())[i]] = i;  }

		groupNames = myGroups;
		m->setGroups(myGroups);

	}
	catch(exception& e) {
		m->errorOut(e, "LibShuffCommand", "setGroups");
		exit(1);
	}
}

/***********************************************************/
