/*
 *  rarefactcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "rarefactcommand.h"
#include "ace.h"
#include "sobs.h"
#include "nseqs.h"
#include "chao1.h"
#include "bootstrap.h"
#include "simpson.h"
#include "simpsoneven.h"
#include "heip.h"
#include "smithwilson.h"
#include "invsimpson.h"
#include "npshannon.h"
#include "shannoneven.h"
#include "shannon.h"
#include "jackknife.h"
#include "coverage.h"


//**********************************************************************************************************************
vector<string> RareFactCommand::setParameters(){	
	try {
		CommandParameter plist("list", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(plist);
		CommandParameter prabund("rabund", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(prabund);
		CommandParameter psabund("sabund", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(psabund);
		CommandParameter pshared("shared", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(pshared);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pfreq("freq", "Number", "", "100", "", "", "",false,false); parameters.push_back(pfreq);
		CommandParameter piters("iters", "Number", "", "1000", "", "", "",false,false); parameters.push_back(piters);
		CommandParameter pcalc("calc", "Multiple", "sobs-chao-nseqs-coverage-ace-jack-shannon-shannoneven-npshannon-heip-smithwilson-simpson-simpsoneven-invsimpson-bootstrap", "sobs", "", "", "",true,false); parameters.push_back(pcalc);
		CommandParameter pabund("abund", "Number", "", "10", "", "", "",false,false); parameters.push_back(pabund);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "",false,false); parameters.push_back(pprocessors);
		CommandParameter pgroupmode("groupmode", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(pgroupmode);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string RareFactCommand::getHelpString(){	
	try {
		ValidCalculators validCalculator;
		string helpString = "";
		helpString += "The rarefaction.single command parameters are list, sabund, rabund, shared, label, iters, freq, calc, processors and abund.  list, sabund, rabund or shared is required unless you have a valid current file. \n";
		helpString += "The freq parameter is used indicate when to output your data, by default it is set to 100. But you can set it to a percentage of the number of sequence. For example freq=0.10, means 10%. \n";
		helpString += "The processors parameter allows you to specify the number of processors to use. The default is 1.\n";
		helpString += "The rarefaction.single command should be in the following format: \n";
		helpString += "rarefaction.single(label=yourLabel, iters=yourIters, freq=yourFreq, calc=yourEstimators).\n";
		helpString += "Example rarefaction.single(label=unique-.01-.03, iters=10000, freq=10, calc=sobs-rchao-race-rjack-rbootstrap-rshannon-rnpshannon-rsimpson).\n";
		helpString += "The default values for iters is 1000, freq is 100, and calc is rarefaction which calculates the rarefaction curve for the observed richness.\n";
		validCalculator.printCalc("rarefaction");
		helpString += "If you are running rarefaction.single with a shared file and would like your results collated in one file, set groupmode=t. (Default=true).\n";
		helpString += "The label parameter is used to analyze specific labels in your input.\n";
		helpString += "Note: No spaces between parameter labels (i.e. freq), '=' and parameters (i.e.yourFreq).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
RareFactCommand::RareFactCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["rarefaction"] = tempOutNames;
		outputTypes["r_chao"] = tempOutNames;
		outputTypes["r_ace"] = tempOutNames;
		outputTypes["r_jack"] = tempOutNames;
		outputTypes["r_shannon"] = tempOutNames;
		outputTypes["r_shannoneven"] = tempOutNames;
		outputTypes["r_heip"] = tempOutNames;
		outputTypes["r_smithwilson"] = tempOutNames;
		outputTypes["r_npshannon"] = tempOutNames;
		outputTypes["r_simpson"] = tempOutNames;
		outputTypes["r_simpsoneven"] = tempOutNames;
		outputTypes["r_invsimpson"] = tempOutNames;
		outputTypes["r_bootstrap"] = tempOutNames;
		outputTypes["r_coverage"] = tempOutNames;
		outputTypes["r_nseqs"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactCommand", "RareFactCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
RareFactCommand::RareFactCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
						
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			map<string,string>::iterator it;
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["rarefaction"] = tempOutNames;
			outputTypes["r_chao"] = tempOutNames;
			outputTypes["r_ace"] = tempOutNames;
			outputTypes["r_jack"] = tempOutNames;
			outputTypes["r_shannon"] = tempOutNames;
			outputTypes["r_shannoneven"] = tempOutNames;
			outputTypes["r_heip"] = tempOutNames;
			outputTypes["r_smithwilson"] = tempOutNames;
			outputTypes["r_npshannon"] = tempOutNames;
			outputTypes["r_simpson"] = tempOutNames;
			outputTypes["r_simpsoneven"] = tempOutNames;
			outputTypes["r_invsimpson"] = tempOutNames;
			outputTypes["r_bootstrap"] = tempOutNames;
			outputTypes["r_coverage"] = tempOutNames;
			outputTypes["r_nseqs"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
				
				it = parameters.find("rabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["rabund"] = inputDir + it->second;		}
				}
				
				it = parameters.find("sabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["sabund"] = inputDir + it->second;		}
				}
				
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { listfile = ""; abort = true; }
			else if (listfile == "not found") { listfile = ""; }
			else {  format = "list"; inputfile = listfile; m->setListFile(listfile); }
			
			sabundfile = validParameter.validFile(parameters, "sabund", true);
			if (sabundfile == "not open") { sabundfile = ""; abort = true; }	
			else if (sabundfile == "not found") { sabundfile = ""; }
			else {  format = "sabund"; inputfile = sabundfile; m->setSabundFile(sabundfile); }
			
			rabundfile = validParameter.validFile(parameters, "rabund", true);
			if (rabundfile == "not open") { rabundfile = ""; abort = true; }	
			else if (rabundfile == "not found") { rabundfile = ""; }
			else {  format = "rabund"; inputfile = rabundfile; m->setRabundFile(rabundfile); }
			
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }	
			else if (sharedfile == "not found") { sharedfile = ""; }
			else {  format = "sharedfile"; inputfile = sharedfile; m->setSharedFile(sharedfile); }
				
			if ((sharedfile == "") && (listfile == "") && (rabundfile == "") && (sabundfile == "")) { 
				//is there are current file available for any of these?
				//give priority to shared, then list, then rabund, then sabund
				//if there is a current shared file, use it
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") { inputfile = sharedfile; format = "sharedfile"; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 
					listfile = m->getListFile(); 
					if (listfile != "") { inputfile = listfile; format = "list"; m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
					else { 
						rabundfile = m->getRabundFile(); 
						if (rabundfile != "") { inputfile = rabundfile; format = "rabund"; m->mothurOut("Using " + rabundfile + " as input file for the rabund parameter."); m->mothurOutEndLine(); }
						else { 
							sabundfile = m->getSabundFile(); 
							if (sabundfile != "") { inputfile = sabundfile; format = "sabund"; m->mothurOut("Using " + sabundfile + " as input file for the sabund parameter."); m->mothurOutEndLine(); }
							else { 
								m->mothurOut("No valid current files. You must provide a list, sabund, rabund or shared file before you can use the collect.single command."); m->mothurOutEndLine(); 
								abort = true;
							}
						}
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
			if (calc == "not found") { calc = "sobs";  }
			else { 
				 if (calc == "default")  {  calc = "sobs";  }
			}
			m->splitAtDash(calc, Estimators);
			if (m->inUsersGroups("citation", Estimators)) { 
				ValidCalculators validCalc; validCalc.printCitations(Estimators); 
				//remove citation from list of calcs
				for (int i = 0; i < Estimators.size(); i++) { if (Estimators[i] == "citation") {  Estimators.erase(Estimators.begin()+i); break; } }
			}

			string temp;
			temp = validParameter.validFile(parameters, "freq", false);			if (temp == "not found") { temp = "100"; }
			m->mothurConvert(temp, freq); 
			
			temp = validParameter.validFile(parameters, "abund", false);			if (temp == "not found") { temp = "10"; }
			m->mothurConvert(temp, abund); 
			
			temp = validParameter.validFile(parameters, "iters", false);			if (temp == "not found") { temp = "1000"; }
			m->mothurConvert(temp, nIters); 
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
			
			temp = validParameter.validFile(parameters, "groupmode", false);		if (temp == "not found") { temp = "T"; }
			groupMode = m->isTrue(temp);
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactCommand", "RareFactCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int RareFactCommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		if ((format != "sharedfile")) { inputFileNames.push_back(inputfile);  }
		else {  inputFileNames = parseSharedFile(sharedfile);  format = "rabund"; }
				
		if (m->control_pressed) { return 0; }
		
		map<int, string> file2Group; //index in outputNames[i] -> group
		for (int p = 0; p < inputFileNames.size(); p++) {
			
			string fileNameRoot = outputDir + m->getRootName(m->getSimpleName(inputFileNames[p]));
						
			if (m->control_pressed) {  outputTypes.clear(); for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	}  m->clearGroups();  return 0; }
			
			if (inputFileNames.size() > 1) {
				m->mothurOutEndLine(); m->mothurOut("Processing group " + groups[p]); m->mothurOutEndLine(); m->mothurOutEndLine();
			}
			int i;
			ValidCalculators validCalculator;
			
			  
			for (i=0; i<Estimators.size(); i++) {
				if (validCalculator.isValidCalculator("rarefaction", Estimators[i]) == true) { 
					if (Estimators[i] == "sobs") { 
						rDisplays.push_back(new RareDisplay(new Sobs(), new ThreeColumnFile(fileNameRoot+"rarefaction")));
						outputNames.push_back(fileNameRoot+"rarefaction"); outputTypes["rarefaction"].push_back(fileNameRoot+"rarefaction");
					}else if (Estimators[i] == "chao") { 
						rDisplays.push_back(new RareDisplay(new Chao1(), new ThreeColumnFile(fileNameRoot+"r_chao")));
						outputNames.push_back(fileNameRoot+"r_chao"); outputTypes["r_chao"].push_back(fileNameRoot+"r_chao");
					}else if (Estimators[i] == "ace") { 
						if(abund < 5)
							abund = 10;
						rDisplays.push_back(new RareDisplay(new Ace(abund), new ThreeColumnFile(fileNameRoot+"r_ace")));
						outputNames.push_back(fileNameRoot+"r_ace"); outputTypes["r_ace"].push_back(fileNameRoot+"r_ace");
					}else if (Estimators[i] == "jack") { 
						rDisplays.push_back(new RareDisplay(new Jackknife(), new ThreeColumnFile(fileNameRoot+"r_jack")));
						outputNames.push_back(fileNameRoot+"r_jack"); outputTypes["r_jack"].push_back(fileNameRoot+"r_jack");
					}else if (Estimators[i] == "shannon") { 
						rDisplays.push_back(new RareDisplay(new Shannon(), new ThreeColumnFile(fileNameRoot+"r_shannon")));
						outputNames.push_back(fileNameRoot+"r_shannon"); outputTypes["r_shannon"].push_back(fileNameRoot+"r_shannon");
					}else if (Estimators[i] == "shannoneven") { 
						rDisplays.push_back(new RareDisplay(new ShannonEven(), new ThreeColumnFile(fileNameRoot+"r_shannoneven")));
						outputNames.push_back(fileNameRoot+"r_shannoneven"); outputTypes["r_shannoneven"].push_back(fileNameRoot+"r_shannoneven");
					}else if (Estimators[i] == "heip") { 
						rDisplays.push_back(new RareDisplay(new Heip(), new ThreeColumnFile(fileNameRoot+"r_heip")));
						outputNames.push_back(fileNameRoot+"r_heip"); outputTypes["r_heip"].push_back(fileNameRoot+"r_heip");
					}else if (Estimators[i] == "smithwilson") { 
						rDisplays.push_back(new RareDisplay(new SmithWilson(), new ThreeColumnFile(fileNameRoot+"r_smithwilson")));
						outputNames.push_back(fileNameRoot+"r_smithwilson"); outputTypes["r_smithwilson"].push_back(fileNameRoot+"r_smithwilson");
					}else if (Estimators[i] == "npshannon") { 
						rDisplays.push_back(new RareDisplay(new NPShannon(), new ThreeColumnFile(fileNameRoot+"r_npshannon")));
						outputNames.push_back(fileNameRoot+"r_npshannon"); outputTypes["r_npshannon"].push_back(fileNameRoot+"r_npshannon");
					}else if (Estimators[i] == "simpson") { 
						rDisplays.push_back(new RareDisplay(new Simpson(), new ThreeColumnFile(fileNameRoot+"r_simpson")));
						outputNames.push_back(fileNameRoot+"r_simpson"); outputTypes["r_simpson"].push_back(fileNameRoot+"r_simpson");
					}else if (Estimators[i] == "simpsoneven") { 
						rDisplays.push_back(new RareDisplay(new SimpsonEven(), new ThreeColumnFile(fileNameRoot+"r_simpsoneven")));
						outputNames.push_back(fileNameRoot+"r_simpsoneven"); outputTypes["r_simpsoneven"].push_back(fileNameRoot+"r_simpsoneven");
					}else if (Estimators[i] == "invsimpson") { 
						rDisplays.push_back(new RareDisplay(new InvSimpson(), new ThreeColumnFile(fileNameRoot+"r_invsimpson")));
						outputNames.push_back(fileNameRoot+"r_invsimpson"); outputTypes["r_invsimpson"].push_back(fileNameRoot+"r_invsimpson");
					}else if (Estimators[i] == "bootstrap") { 
						rDisplays.push_back(new RareDisplay(new Bootstrap(), new ThreeColumnFile(fileNameRoot+"r_bootstrap")));
						outputNames.push_back(fileNameRoot+"r_bootstrap"); outputTypes["r_bootstrap"].push_back(fileNameRoot+"r_bootstrap");
					}else if (Estimators[i] == "coverage") { 
						rDisplays.push_back(new RareDisplay(new Coverage(), new ThreeColumnFile(fileNameRoot+"r_coverage")));
						outputNames.push_back(fileNameRoot+"r_coverage"); outputTypes["r_coverage"].push_back(fileNameRoot+"r_coverage");
					}else if (Estimators[i] == "nseqs") { 
						rDisplays.push_back(new RareDisplay(new NSeqs(), new ThreeColumnFile(fileNameRoot+"r_nseqs")));
						outputNames.push_back(fileNameRoot+"r_nseqs"); outputTypes["r_nseqs"].push_back(fileNameRoot+"r_nseqs");
					}
                    file2Group[outputNames.size()-1] = groups[p];
				}
			}
			
			
			//if the users entered no valid calculators don't execute command
			if (rDisplays.size() == 0) { for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	}  return 0; }
			
			input = new InputData(inputFileNames[p], format);			
			order = input->getOrderVector();
			string lastLabel = order->getLabel();
			
			//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
			set<string> processedLabels;
			set<string> userLabels = labels;
			
			if (m->control_pressed) { for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	}  delete input;  delete order;  for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
			
			//as long as you are not at the end of the file or done wih the lines you want
			while((order != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
				
				if (m->control_pressed) { for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	}  delete input;  delete order;  for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }

				
				if(allLines == 1 || labels.count(order->getLabel()) == 1){
					
					m->mothurOut(order->getLabel()); m->mothurOutEndLine();
					rCurve = new Rarefact(order, rDisplays, processors);
					rCurve->getCurve(freq, nIters);
					delete rCurve;
					
					processedLabels.insert(order->getLabel());
					userLabels.erase(order->getLabel());
				}
				
				if ((m->anyLabelsToProcess(order->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = order->getLabel();
					
					delete order;
					order = (input->getOrderVector(lastLabel));
					
					m->mothurOut(order->getLabel()); m->mothurOutEndLine();
					rCurve = new Rarefact(order, rDisplays, processors);
					rCurve->getCurve(freq, nIters);
					delete rCurve;
					
					processedLabels.insert(order->getLabel());
					userLabels.erase(order->getLabel());
					
					//restore real lastlabel to save below
					order->setLabel(saveLabel);
				}
				
				lastLabel = order->getLabel();		
				
				delete order;
				order = (input->getOrderVector());
			}
			
			if (m->control_pressed) { for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	}  delete input;   for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }

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
			
			if (m->control_pressed) { for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	}  delete input;   for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }

			//run last label if you need to
			if (needToRun == true)  {
				if (order != NULL) {	delete order;	}
				order = (input->getOrderVector(lastLabel));
				
				m->mothurOut(order->getLabel()); m->mothurOutEndLine();
				rCurve = new Rarefact(order, rDisplays, processors);
				rCurve->getCurve(freq, nIters);
				delete rCurve;
				
				delete order;
			}
			
			
			for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	}	
			rDisplays.clear();
			delete input;  
		}
		
		
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }

		//create summary file containing all the groups data for each label - this function just combines the info from the files already created.
		if ((sharedfile != "") && (groupMode)) {   outputNames = createGroupFile(outputNames, file2Group);  }

		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }

		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> RareFactCommand::createGroupFile(vector<string>& outputNames, map<int, string> file2Group) {
	try {
		
		vector<string> newFileNames;
		
		//find different types of files
		map<string, vector<string> > typesFiles;
		for (int i = 0; i < outputNames.size(); i++) {
			string extension = m->getExtension(outputNames[i]);
			
			ifstream in;
			m->openInputFile(outputNames[i], in);
			
			string labels = m->getline(in);
			string newLine = labels.substr(0, labels.find_first_of('\t'));
			
			newLine += "\tGroup" + labels.substr(labels.find_first_of('\t'));
			
			typesFiles[extension].push_back(outputNames[i]+"_"+file2Group[i]);
			
			string combineFileName = outputDir + m->getRootName(m->getSimpleName(sharedfile)) + "groups" + extension;
			
			//print headers
			ofstream out;
			m->openOutputFile(combineFileName, out);
			out << newLine << endl;
			out.close();
			
		}
		
		//for each type create a combo file
		map<int, int> lineToNumber; 
		for (map<string, vector<string> >::iterator it = typesFiles.begin(); it != typesFiles.end(); it++) {
			
			ofstream out;
			string combineFileName = outputDir + m->getRootName(m->getSimpleName(sharedfile)) + "groups" + it->first;
			m->openOutputFileAppend(combineFileName, out);
			newFileNames.push_back(combineFileName);
			
			vector<string> thisTypesFiles = it->second;
		
			//open each type summary file
			map<string, vector<string> > files; //maps file name to lines in file
			int maxLines = 0;
			int numColumns = 0;
			for (int i=0; i<thisTypesFiles.size(); i++) {
                
                string thisfilename = thisTypesFiles[i].substr(0, thisTypesFiles[i].find_last_of('_'));
                string group = thisTypesFiles[i].substr(thisTypesFiles[i].find_last_of('_')+1);
                
				ifstream temp;
				m->openInputFile(thisfilename, temp);
				
				//read through first line - labels
				m->getline(temp);	m->gobble(temp);
				
				vector<string> thisFilesLines;
				//string fileNameRoot = m->getRootName(thisTypesFiles[i]);
				//map<string, string>::iterator itName = nameMap.find(fileNameRoot);
				//string group = "";
				//if (itName != nameMap.end()) {
				//	group = itName->second;
				//}else {
				//	group = "not found" + i;
				//	m->mothurOut("[ERROR]: can't parse filename."); m->mothurOutEndLine();
				//}
				
				thisFilesLines.push_back(group);
				int count = 1;
				while (!temp.eof()){
				
					string thisLine = m->getline(temp);
					
					string numSampled = thisLine.substr(0, thisLine.find_first_of('\t'));
					int num = 0;
					convert(numSampled, num);
					numColumns = m->getNumChar(thisLine, '\t');
					lineToNumber[count] = num;
					count++;
									
					thisFilesLines.push_back(thisLine);
					
					m->gobble(temp);
				}
				
				files[thisTypesFiles[i]] = thisFilesLines;
				
				//save longest file for below
				if (maxLines < thisFilesLines.size()) { maxLines = thisFilesLines.size(); }
				
				temp.close();
				m->mothurRemove(thisfilename);
			}
			
			
			//for each label
			for (int k = 1; k < maxLines; k++) {
				
				//grab data for each group
				for (int i=0; i<thisTypesFiles.size(); i++) {
					
					map<int, int>::iterator itLine = lineToNumber.find(k);
					if (itLine != lineToNumber.end()) {
						string output = toString(itLine->second);
						if (k < files[thisTypesFiles[i]].size()) {
							string line = files[thisTypesFiles[i]][k];
							output = line.substr(0, line.find_first_of('\t'));
							output += '\t' + files[thisTypesFiles[i]][0] + '\t' + line.substr(line.find_first_of('\t'));
						}else{
							output += '\t' + files[thisTypesFiles[i]][0] + '\t';
							for (int h = 0; h < numColumns; h++) {
								output += "NA\t";
							}
						}
						out << output << endl;
					}else { m->mothurOut("[ERROR]: parsing results, cant find " + toString(k)); m->mothurOutEndLine(); }
				}
			}	
			
			out.close();
			
		}
		
		//return combine file name
		return newFileNames;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactCommand", "createGroupFile");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> RareFactCommand::parseSharedFile(string filename) {
	try {
		vector<string> filenames;
		
		map<string, ofstream*> filehandles;
		map<string, ofstream*>::iterator it3;
		
		input = new InputData(filename, "sharedfile");
		vector<SharedRAbundVector*> lookup = input->getSharedRAbundVectors();
		
		string sharedFileRoot = m->getRootName(filename);
		
		//clears file before we start to write to it below
		for (int i=0; i<lookup.size(); i++) {
			m->mothurRemove((sharedFileRoot + lookup[i]->getGroup() + ".rabund"));
			filenames.push_back((sharedFileRoot + lookup[i]->getGroup() + ".rabund"));
		}
		
		ofstream* temp;
		for (int i=0; i<lookup.size(); i++) {
			temp = new ofstream;
			filehandles[lookup[i]->getGroup()] = temp;
			groups.push_back(lookup[i]->getGroup());
		}

		while(lookup[0] != NULL) {
		
			for (int i = 0; i < lookup.size(); i++) {
				RAbundVector rav = lookup[i]->getRAbundVector();
				m->openOutputFileAppend(sharedFileRoot + lookup[i]->getGroup() + ".rabund", *(filehandles[lookup[i]->getGroup()]));
				rav.print(*(filehandles[lookup[i]->getGroup()]));
				(*(filehandles[lookup[i]->getGroup()])).close();
			}
		
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
			lookup = input->getSharedRAbundVectors();
		}
		
		//free memory
		for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {
			delete it3->second;
		}
		
		delete input;
		m->clearGroups();

		return filenames;
	}
	catch(exception& e) {
		m->errorOut(e, "RareFactCommand", "parseSharedFile");
		exit(1);
	}
}
//**********************************************************************************************************************



