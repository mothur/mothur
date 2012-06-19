/*
 *  heatmapsimcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "heatmapsimcommand.h"
#include "sharedjabund.h"
#include "sharedsorabund.h"
#include "sharedjclass.h"
#include "sharedsorclass.h"
#include "sharedjest.h"
#include "sharedsorest.h"
#include "sharedthetayc.h"
#include "sharedthetan.h"
#include "sharedmorisitahorn.h"
#include "sharedbraycurtis.h"


//**********************************************************************************************************************
vector<string> HeatMapSimCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "PhylipColumnShared", "PhylipColumnShared", "none",false,false); parameters.push_back(pshared);	
		CommandParameter pphylip("phylip", "InputTypes", "", "", "PhylipColumnShared", "PhylipColumnShared", "none",false,false); parameters.push_back(pphylip);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "ColumnName",false,false); parameters.push_back(pname);
		CommandParameter pcolumn("column", "InputTypes", "", "", "PhylipColumnShared", "PhylipColumnShared", "ColumnName",false,false); parameters.push_back(pcolumn);		
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pcalc("calc", "Multiple", "jabund-sorabund-jclass-sorclass-jest-sorest-thetayc-thetan-morisitahorn-braycurtis", "jest-thetayc", "", "", "",true,false); parameters.push_back(pcalc);
		CommandParameter pfontsize("fontsize", "Number", "", "24", "", "", "",false,false); parameters.push_back(pfontsize);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMapSimCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string HeatMapSimCommand::getHelpString(){	
	try {
		string helpString = "";
		ValidCalculators validCalculator;
		helpString += "The heatmap.sim command parameters are shared, phylip, column, name, groups, calc, fontsize and label.  shared or phylip or column and name are required unless valid current files exist.\n";
		helpString += "There are two ways to use the heatmap.sim command. The first is with the read.otu command. \n";
		helpString += "With the read.otu command you may use the groups, label and calc parameters. \n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like included in your heatmap.\n";
		helpString += "The group names are separated by dashes. The label parameter allows you to select what distance levels you would like a heatmap created for, and is also separated by dashes.\n";
		helpString += "The fontsize parameter allows you to adjust the font size of the picture created, default=24.\n";
		helpString += "The heatmap.sim command should be in the following format: heatmap.sim(groups=yourGroups, calc=yourCalc, label=yourLabels).\n";
		helpString += "Example heatmap.sim(groups=A-B-C, calc=jabund).\n";
		helpString += "The default value for groups is all the groups in your groupfile, and all labels in your inputfile will be used.\n";
		helpString +=  validCalculator.printCalc("heat");
		helpString += "The default value for calc is jclass-thetayc.\n";
		helpString += "The heatmap.sim command outputs a .svg file for each calculator you choose at each label you specify.\n";
		helpString += "The second way to use the heatmap.sim command is with a distance file representing the distance bewteen your groups. \n";
		helpString += "Using the command this way, the phylip or column parameter are required, and only one may be used.  If you use a column file the name filename is required. \n";
		helpString += "The heatmap.sim command should be in the following format: heatmap.sim(phylip=yourDistanceFile).\n";
		helpString += "Example heatmap.sim(phylip=amazonGroups.dist).\n";
		helpString += "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMapSimCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************

string HeatMapSimCommand::getOutputFileNameTag(string type, string inputName=""){	
	try {
        string outputFileName = "";
		map<string, vector<string> >::iterator it;
        
        //is this a type this command creates
        it = outputTypes.find(type);
        if (it == outputTypes.end()) {  m->mothurOut("[ERROR]: this command doesn't create a " + type + " output file.\n"); }
        else {
            if (type == "svg")            {   outputFileName =  "svg";   }
            else { m->mothurOut("[ERROR]: No definition for type " + type + " output file tag.\n"); m->control_pressed = true;  }
        }
        return outputFileName;
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMapSimCommand", "getOutputFileNameTag");
		exit(1);
	}
}

//**********************************************************************************************************************
HeatMapSimCommand::HeatMapSimCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["svg"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMapSimCommand", "HeatMapSimCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

HeatMapSimCommand::HeatMapSimCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
			
		//allow user to run help
		if(option == "help") {  help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["svg"] = tempOutNames;
			
			format = "";
				
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
				
				it = parameters.find("column");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["column"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
				
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
			}

			//required parameters
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not open") { abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; }	
			else {  format = "phylip"; 	inputfile = phylipfile; m-> setPhylipFile(phylipfile); if (outputDir == "") { outputDir += m->hasPath(phylipfile); }  }
			
			columnfile = validParameter.validFile(parameters, "column", true);
			if (columnfile == "not open") { abort = true; }	
			else if (columnfile == "not found") { columnfile = ""; }
			else {  format = "column";	inputfile = columnfile; m->setColumnFile(columnfile); if (outputDir == "") { outputDir += m->hasPath(columnfile); } }
			
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { abort = true; }	
			else if (sharedfile == "not found") { sharedfile = ""; }
			else {  format = "shared";	inputfile = sharedfile; m->setSharedFile(sharedfile); if (outputDir == "") { outputDir += m->hasPath(sharedfile); } }
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			else { m->setNameFile(namefile); }
			
			
			//error checking on files			
			if ((sharedfile == "") && ((phylipfile == "") && (columnfile == "")))	{ 
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") { format = "shared"; inputfile = sharedfile; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 
					//is there are current file available for either of these?
					//give priority to column, then phylip
					columnfile = m->getColumnFile(); 
					if (columnfile != "") {  format = "column"; inputfile = columnfile; m->mothurOut("Using " + columnfile + " as input file for the column parameter."); m->mothurOutEndLine(); }
					else { 
						phylipfile = m->getPhylipFile(); 
						if (phylipfile != "") { format = "phylip";  inputfile = phylipfile; m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter."); m->mothurOutEndLine(); }
						else { 
							m->mothurOut("No valid current files. You must provide a shared or phylip or column file."); m->mothurOutEndLine(); 
							abort = true;
						}
					}
				}
			}
			else if ((phylipfile != "") && (columnfile != "")) { m->mothurOut("When running the heatmap.sim command with a distance file you may not use both the column and the phylip parameters."); m->mothurOutEndLine(); abort = true; }
			
			if (columnfile != "") {
				if (namefile == "") { 
					namefile = m->getNameFile(); 
					if (namefile != "") {  m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("You need to provide a namefile if you are going to use the column format."); m->mothurOutEndLine(); 
						abort = true; 
					}	
				}
			}
			
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(inputfile);		}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
							
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
				
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "jest-thetayc";  }
			else { 
				if (calc == "default")  {  calc = "jest-thetayc";  }
			}
			m->splitAtDash(calc, Estimators);
			if (m->inUsersGroups("citation", Estimators)) { 
				ValidCalculators validCalc; validCalc.printCitations(Estimators); 
				//remove citation from list of calcs
				for (int i = 0; i < Estimators.size(); i++) { if (Estimators[i] == "citation") {  Estimators.erase(Estimators.begin()+i); break; } }
			}
				
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				m->splitAtDash(groups, Groups);
				m->setGroups(Groups);
			}
			
			string temp = validParameter.validFile(parameters, "fontsize", false);				if (temp == "not found") { temp = "24"; }
			m->mothurConvert(temp, fontsize);
			
			if (abort == false) {
				ValidCalculators validCalculator;
			
				int i;
				for (i=0; i<Estimators.size(); i++) {
					if (validCalculator.isValidCalculator("heat", Estimators[i]) == true) { 
						if (Estimators[i] == "jabund") { 	
							heatCalculators.push_back(new JAbund());
						}else if (Estimators[i] == "sorabund") { 
							heatCalculators.push_back(new SorAbund());
						}else if (Estimators[i] == "jclass") { 
							heatCalculators.push_back(new Jclass());
						}else if (Estimators[i] == "sorclass") { 
							heatCalculators.push_back(new SorClass());
						}else if (Estimators[i] == "jest") { 
							heatCalculators.push_back(new Jest());
						}else if (Estimators[i] == "sorest") { 
							heatCalculators.push_back(new SorEst());
						}else if (Estimators[i] == "thetayc") { 
							heatCalculators.push_back(new ThetaYC());
						}else if (Estimators[i] == "thetan") { 
							heatCalculators.push_back(new ThetaN());
						}else if (Estimators[i] == "morisitahorn") { 
							heatCalculators.push_back(new MorHorn());
						}else if (Estimators[i] == "braycurtis") { 
							heatCalculators.push_back(new BrayCurtis());
						}
					}
				}
				
			}
		}

				

	}
	catch(exception& e) {
		m->errorOut(e, "HeatMapSimCommand", "HeatMapSimCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int HeatMapSimCommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		heatmap = new HeatMapSim(outputDir, inputfile, fontsize);
		
		if (format == "shared") {
			runCommandShared();
		}else{	runCommandDist();	}
		
		delete heatmap;
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } outputTypes.clear(); return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMapSimCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
int HeatMapSimCommand::runCommandShared() {
	try {
		//if the users entered no valid calculators don't execute command
		if (heatCalculators.size() == 0) { m->mothurOut("No valid calculators."); m->mothurOutEndLine(); return 0; }
		
		input = new InputData(sharedfile, "sharedfile");
		lookup = input->getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
			
		if (lookup.size() < 2) { m->mothurOut("You have not provided enough valid groups.  I cannot run the command."); m->mothurOutEndLine(); return 0;}
				
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		if (m->control_pressed) {  delete input;  for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  m->clearGroups(); return 0; }
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (m->control_pressed) { delete input;  for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } m->clearGroups(); return 0; }

			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			
	
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				vector<string> outfilenames = heatmap->getPic(lookup, heatCalculators);
				for(int i = 0; i < outfilenames.size(); i++) { outputNames.push_back(outfilenames[i]);  outputTypes["svg"].push_back(outfilenames[i]); }
					
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}
				
			if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup[0]->getLabel();
			
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
				lookup = input->getSharedRAbundVectors(lastLabel);				

				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				vector<string> outfilenames = heatmap->getPic(lookup, heatCalculators);
				for(int i = 0; i < outfilenames.size(); i++) { outputNames.push_back(outfilenames[i]); outputTypes["svg"].push_back(outfilenames[i]);  }
					
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				
				//restore real lastlabel to save below
				lookup[0]->setLabel(saveLabel);
			}
				
			//prevent memory leak
			 
			lastLabel = lookup[0]->getLabel();			

			//get next line to process
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
			lookup = input->getSharedRAbundVectors();

		}
		
			
		if (m->control_pressed) {  delete input;  m->clearGroups();  return 0; }

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
		
		if (m->control_pressed) {  delete input;  m->clearGroups(); return 0; }
		
		//run last label if you need to
		if (needToRun == true)  {
			for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) { delete lookup[i]; } } 
			lookup = input->getSharedRAbundVectors(lastLabel);				

			m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
			vector<string> outfilenames = heatmap->getPic(lookup, heatCalculators);
			for(int i = 0; i < outfilenames.size(); i++) { outputNames.push_back(outfilenames[i]); outputTypes["svg"].push_back(outfilenames[i]);  }
			
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
		}
		
		if (m->control_pressed) {  delete input;  m->clearGroups(); return 0; }
			
		//reset groups parameter
		m->clearGroups();  
			
		delete input;  
	
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMapSimCommand", "runCommandShared");
		exit(1);
	}
}
//**********************************************************************************************************************
int HeatMapSimCommand::runCommandDist() {
	try {
	
		vector< vector<double> > matrix;
		vector<string> names;
		ifstream in;
		
		//read distance file and create distance vector and names vector
		if (format == "phylip") {
			//read phylip file
			m->openInputFile(phylipfile, in);
			
			string name;
			int numSeqs;
			in >> numSeqs >> name; 
			
			//save name
			names.push_back(name);
		
			//resize the matrix and fill with zeros
			matrix.resize(numSeqs); 
			for(int i = 0; i < numSeqs; i++) {
				matrix[i].resize(numSeqs, 0.0);
			}
					
			//determine if matrix is square or lower triangle
			//if it is square read the distances for the first sequence
			char d;
			bool square;
			while((d=in.get()) != EOF){
				
				//is d a number meaning its square
				if(isalnum(d)){ 
					square = true;
					in.putback(d);
					
					for(int i=0;i<numSeqs;i++){
						in >> matrix[0][i];
					}
					break;
				}
				
				//is d a line return meaning its lower triangle
				if(d == '\n'){
					square = false;
					break;
				}
			}
			
			//read rest of matrix
			if (square == true) { 
				for(int i=1;i<numSeqs;i++){
					in >> name;		
					names.push_back(name);
					
					if (m->control_pressed) { return 0; }
					
					for(int j=0;j<numSeqs;j++) { in >> matrix[i][j];  }
					m->gobble(in);
				}
			}else { 
				double dist;
				for(int i=1;i<numSeqs;i++){
					in >> name;	
					names.push_back(name);	
					
					if (m->control_pressed) { return 0; }
					
					for(int j=0;j<i;j++){
						in >> dist;
						matrix[i][j] = dist;  matrix[j][i] = dist;
					}
					m->gobble(in);
				}
			}
			in.close();
		}else {
			//read names file
			NameAssignment* nameMap = new NameAssignment(namefile);
			nameMap->readMap();
			
			//put names in order in vector
			for (int i = 0; i < nameMap->size(); i++) {
				names.push_back(nameMap->get(i));
			}
			
			//resize matrix
			matrix.resize(nameMap->size());
			for (int i = 0; i < nameMap->size(); i++) {
				matrix[i].resize(nameMap->size(), 0.0);
			}
			
			//read column file
			string first, second;
			double dist;
			m->openInputFile(columnfile, in);
			
			while (!in.eof()) {
				in >> first >> second >> dist; m->gobble(in);
				
				if (m->control_pressed) { return 0; }
				
				map<string, int>::iterator itA = nameMap->find(first);
				map<string, int>::iterator itB = nameMap->find(second);
				
				if(itA == nameMap->end()){  m->mothurOut("AAError: Sequence '" + first + "' was not found in the names file, please correct\n"); exit(1);  }
				if(itB == nameMap->end()){  m->mothurOut("ABError: Sequence '" + second + "' was not found in the names file, please correct\n"); exit(1);  }
				
				//save distance
				matrix[itA->second][itB->second] = dist;
				matrix[itB->second][itA->second] = dist;
			}
			in.close();
			
			delete nameMap;
		}
		
		
		string outputFileName = heatmap->getPic(matrix, names);
		outputNames.push_back(outputFileName); //vector<vector<double>>, vector<string>
		outputTypes["svg"].push_back(outputFileName);
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMapSimCommand", "runCommandDist");
		exit(1);
	}
}
//**********************************************************************************************************************






