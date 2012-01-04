/*
 *  collectcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "collectcommand.h"
#include "ace.h"
#include "sobs.h"
#include "nseqs.h"
#include "chao1.h"
#include "bootstrap.h"
#include "simpson.h"
#include "simpsoneven.h"
#include "invsimpson.h"
#include "npshannon.h"
#include "shannon.h"
#include "smithwilson.h"
#include "heip.h"
#include "shannoneven.h"
#include "jackknife.h"
#include "geom.h"
#include "qstat.h"
#include "logsd.h"
#include "bergerparker.h"
#include "bstick.h"
#include "goodscoverage.h"
#include "efron.h"
#include "boneh.h"
#include "solow.h"
#include "shen.h"
#include "coverage.h"


//**********************************************************************************************************************
vector<string> CollectCommand::setParameters(){	
	try {
		CommandParameter plist("list", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(plist);
		CommandParameter prabund("rabund", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(prabund);
		CommandParameter psabund("sabund", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(psabund);
		CommandParameter pshared("shared", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(pshared);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pfreq("freq", "Number", "", "100", "", "", "",false,false); parameters.push_back(pfreq);
		CommandParameter pcalc("calc", "Multiple", "sobs-chao-nseqs-coverage-ace-jack-shannon-shannoneven-np_shannon-heip-smithwilson-simpson-simpsoneven-invsimpson-bootstrap-geometric-qstat-logseries-bergerparker-bstick-goodscoverage-efron-boneh-solow-shen", "sobs-chao-ace-jack-shannon-npshannon-simpson", "", "", "",true,false); parameters.push_back(pcalc);
		CommandParameter pabund("abund", "Number", "", "10", "", "", "",false,false); parameters.push_back(pabund);
		CommandParameter psize("size", "Number", "", "0", "", "", "",false,false); parameters.push_back(psize);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "CollectCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string CollectCommand::getHelpString(){	
	try {
		string helpString = "";
		ValidCalculators validCalculator;
		helpString += "The collect.single command parameters are list, sabund, rabund, shared, label, freq, calc and abund.  list, sabund, rabund or shared is required unless you have a valid current file. \n";
		helpString += "The collect.single command should be in the following format: \n";
		helpString += "The freq parameter is used indicate when to output your data, by default it is set to 100. But you can set it to a percentage of the number of sequence. For example freq=0.10, means 10%. \n";
		helpString += "collect.single(label=yourLabel, iters=yourIters, freq=yourFreq, calc=yourEstimators).\n";
		helpString += "Example collect(label=unique-.01-.03, iters=10000, freq=10, calc=sobs-chao-ace-jack).\n";
		helpString += "The default values for freq is 100, and calc are sobs-chao-ace-jack-shannon-npshannon-simpson.\n";
		helpString += validCalculator.printCalc("single");
		helpString += "The label parameter is used to analyze specific labels in your input.\n";
		helpString += "Note: No spaces between parameter labels (i.e. freq), '=' and parameters (i.e.yourFreq).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "CollectCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
CollectCommand::CollectCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["sobs"] = tempOutNames;
		outputTypes["chao"] = tempOutNames;
		outputTypes["nseqs"] = tempOutNames;
		outputTypes["coverage"] = tempOutNames;
		outputTypes["ace"] = tempOutNames;
		outputTypes["jack"] = tempOutNames;
		outputTypes["shannon"] = tempOutNames;
		outputTypes["shannoneven"] = tempOutNames;
		outputTypes["np_shannon"] = tempOutNames;
		outputTypes["heip"] = tempOutNames;
		outputTypes["smithwilson"] = tempOutNames;
		outputTypes["simpson"] = tempOutNames;
		outputTypes["simpsoneven"] = tempOutNames;
		outputTypes["invsimpson"] = tempOutNames;
		outputTypes["bootstrap"] = tempOutNames;
		outputTypes["geometric"] = tempOutNames;
		outputTypes["qstat"] = tempOutNames;
		outputTypes["logseries"] = tempOutNames;
		outputTypes["bergerparker"] = tempOutNames;
		outputTypes["bstick"] = tempOutNames;
		outputTypes["goodscoverage"] = tempOutNames;
		outputTypes["efron"] = tempOutNames;
		outputTypes["boneh"] = tempOutNames;
		outputTypes["solow"] = tempOutNames;
		outputTypes["shen"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "CollectCommand", "CollectCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
CollectCommand::CollectCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
		
		//allow user to run help
		if(option == "help") { help(); calledHelp = true; abort = true; }
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
			outputTypes["sobs"] = tempOutNames;
			outputTypes["chao"] = tempOutNames;
			outputTypes["nseqs"] = tempOutNames;
			outputTypes["coverage"] = tempOutNames;
			outputTypes["ace"] = tempOutNames;
			outputTypes["jack"] = tempOutNames;
			outputTypes["shannon"] = tempOutNames;
			outputTypes["shannoneven"] = tempOutNames;
			outputTypes["np_shannon"] = tempOutNames;
			outputTypes["heip"] = tempOutNames;
			outputTypes["smithwilson"] = tempOutNames;
			outputTypes["simpson"] = tempOutNames;
			outputTypes["simpsoneven"] = tempOutNames;
			outputTypes["invsimpson"] = tempOutNames;
			outputTypes["bootstrap"] = tempOutNames;
			outputTypes["geometric"] = tempOutNames;
			outputTypes["qstat"] = tempOutNames;
			outputTypes["logseries"] = tempOutNames;
			outputTypes["bergerparker"] = tempOutNames;
			outputTypes["bstick"] = tempOutNames;
			outputTypes["goodscoverage"] = tempOutNames;
			outputTypes["efron"] = tempOutNames;
			outputTypes["boneh"] = tempOutNames;
			outputTypes["solow"] = tempOutNames;
			outputTypes["shen"] = tempOutNames;
			
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
			
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
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
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			//NOTE: if you add new calc options, don't forget to add them to the parameter initialize in setParameters or the gui won't be able to use them
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "sobs-chao-ace-jack-shannon-npshannon-simpson";  }
			else { 
				 if (calc == "default")  {  calc = "sobs-chao-ace-jack-shannon-npshannon-simpson";  }
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
			
			temp = validParameter.validFile(parameters, "abund", false);		if (temp == "not found") { temp = "10"; }
			m->mothurConvert(temp, abund); 
			
			temp = validParameter.validFile(parameters, "size", false);			if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, size); 
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "CollectCommand", "CollectCommand");
		exit(1);
	}			
}
//**********************************************************************************************************************

int CollectCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
	
		if ((format != "sharedfile")) { inputFileNames.push_back(inputfile);  }
		else {  inputFileNames = parseSharedFile(sharedfile);  format = "rabund"; }
	
		for (int p = 0; p < inputFileNames.size(); p++) {
			
			if (m->control_pressed) {  outputTypes.clear(); for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	}  m->clearGroups();  return 0; }
			
			if (outputDir == "") { outputDir += m->hasPath(inputFileNames[p]); }
			string fileNameRoot = outputDir + m->getRootName(m->getSimpleName(inputFileNames[p]));
			//globaldata->inputFileName = inputFileNames[p];
		
			if (inputFileNames.size() > 1) {
				m->mothurOutEndLine(); m->mothurOut("Processing group " + groups[p]); m->mothurOutEndLine(); m->mothurOutEndLine();
			}
		
			ValidCalculators validCalculator;
			
			for (int i=0; i<Estimators.size(); i++) {
				if (validCalculator.isValidCalculator("single", Estimators[i]) == true) { 
					if (Estimators[i] == "sobs") { 
						cDisplays.push_back(new CollectDisplay(new Sobs(), new OneColumnFile(fileNameRoot+"sobs")));
						outputNames.push_back(fileNameRoot+"sobs"); outputTypes["sobs"].push_back(fileNameRoot+"sobs");
					}else if (Estimators[i] == "chao") { 
						cDisplays.push_back(new CollectDisplay(new Chao1(), new ThreeColumnFile(fileNameRoot+"chao")));
						outputNames.push_back(fileNameRoot+"chao"); outputTypes["chao"].push_back(fileNameRoot+"chao");
					}else if (Estimators[i] == "nseqs") { 
						cDisplays.push_back(new CollectDisplay(new NSeqs(), new OneColumnFile(fileNameRoot+"nseqs")));
						outputNames.push_back(fileNameRoot+"nseqs"); outputTypes["nseqs"].push_back(fileNameRoot+"nseqs");
					}else if (Estimators[i] == "coverage") { 
						cDisplays.push_back(new CollectDisplay(new Coverage(), new OneColumnFile(fileNameRoot+"coverage")));
						outputNames.push_back(fileNameRoot+"coverage"); outputTypes["coverage"].push_back(fileNameRoot+"coverage");
					}else if (Estimators[i] == "ace") { 
						cDisplays.push_back(new CollectDisplay(new Ace(abund), new ThreeColumnFile(fileNameRoot+"ace")));
						outputNames.push_back(fileNameRoot+"ace"); outputTypes["ace"].push_back(fileNameRoot+"ace");
					}else if (Estimators[i] == "jack") { 
						cDisplays.push_back(new CollectDisplay(new Jackknife(), new ThreeColumnFile(fileNameRoot+"jack")));
						outputNames.push_back(fileNameRoot+"jack"); outputTypes["jack"].push_back(fileNameRoot+"jack");
					}else if (Estimators[i] == "shannon") { 
						cDisplays.push_back(new CollectDisplay(new Shannon(), new ThreeColumnFile(fileNameRoot+"shannon")));
						outputNames.push_back(fileNameRoot+"shannon"); outputTypes["shannon"].push_back(fileNameRoot+"shannon");
					}else if (Estimators[i] == "shannoneven") { 
						cDisplays.push_back(new CollectDisplay(new ShannonEven(), new OneColumnFile(fileNameRoot+"shannoneven")));
						outputNames.push_back(fileNameRoot+"shannoneven"); outputTypes["shannoneven"].push_back(fileNameRoot+"shannoneven");
					}else if (Estimators[i] == "npshannon") { 
						cDisplays.push_back(new CollectDisplay(new NPShannon(), new OneColumnFile(fileNameRoot+"np_shannon")));
						outputNames.push_back(fileNameRoot+"np_shannon"); outputTypes["np_shannon"].push_back(fileNameRoot+"np_shannon");
					}else if (Estimators[i] == "heip") { 
						cDisplays.push_back(new CollectDisplay(new Heip(), new OneColumnFile(fileNameRoot+"heip")));
						outputNames.push_back(fileNameRoot+"heip"); outputTypes["heip"].push_back(fileNameRoot+"heip");
					}else if (Estimators[i] == "smithwilson") { 
						cDisplays.push_back(new CollectDisplay(new SmithWilson(), new OneColumnFile(fileNameRoot+"smithwilson")));
						outputNames.push_back(fileNameRoot+"smithwilson"); outputTypes["smithwilson"].push_back(fileNameRoot+"smithwilson");
					}else if (Estimators[i] == "simpson") { 
						cDisplays.push_back(new CollectDisplay(new Simpson(), new ThreeColumnFile(fileNameRoot+"simpson")));
						outputNames.push_back(fileNameRoot+"simpson"); outputTypes["simpson"].push_back(fileNameRoot+"simpson");
					}else if (Estimators[i] == "simpsoneven") { 
						cDisplays.push_back(new CollectDisplay(new SimpsonEven(), new OneColumnFile(fileNameRoot+"simpsoneven")));
						outputNames.push_back(fileNameRoot+"simpsoneven"); outputTypes["simpsoneven"].push_back(fileNameRoot+"simpsoneven");
					}else if (Estimators[i] == "invsimpson") { 
						cDisplays.push_back(new CollectDisplay(new InvSimpson(), new ThreeColumnFile(fileNameRoot+"invsimpson")));
						outputNames.push_back(fileNameRoot+"invsimpson"); outputTypes["invsimpson"].push_back(fileNameRoot+"invsimpson");
					}else if (Estimators[i] == "bootstrap") { 
						cDisplays.push_back(new CollectDisplay(new Bootstrap(), new OneColumnFile(fileNameRoot+"bootstrap")));
						outputNames.push_back(fileNameRoot+"bootstrap"); outputTypes["bootstrap"].push_back(fileNameRoot+"bootstrap");
					}else if (Estimators[i] == "geometric") { 
						cDisplays.push_back(new CollectDisplay(new Geom(), new OneColumnFile(fileNameRoot+"geometric")));
						outputNames.push_back(fileNameRoot+"geometric"); outputTypes["geometric"].push_back(fileNameRoot+"geometric");
					}else if (Estimators[i] == "qstat") { 
						cDisplays.push_back(new CollectDisplay(new QStat(), new OneColumnFile(fileNameRoot+"qstat")));
						outputNames.push_back(fileNameRoot+"qstat"); outputTypes["qstat"].push_back(fileNameRoot+"qstat");
					}else if (Estimators[i] == "logseries") { 
						cDisplays.push_back(new CollectDisplay(new LogSD(), new OneColumnFile(fileNameRoot+"logseries")));
						outputNames.push_back(fileNameRoot+"logseries"); outputTypes["logseries"].push_back(fileNameRoot+"logseries");
					}else if (Estimators[i] == "bergerparker") { 
						cDisplays.push_back(new CollectDisplay(new BergerParker(), new OneColumnFile(fileNameRoot+"bergerparker")));
						outputNames.push_back(fileNameRoot+"bergerparker"); outputTypes["bergerparker"].push_back(fileNameRoot+"bergerparker");
					}else if (Estimators[i] == "bstick") { 
						cDisplays.push_back(new CollectDisplay(new BStick(), new ThreeColumnFile(fileNameRoot+"bstick")));
						outputNames.push_back(fileNameRoot+"bstick"); outputTypes["bstick"].push_back(fileNameRoot+"bstick");
					}else if (Estimators[i] == "goodscoverage") { 
						cDisplays.push_back(new CollectDisplay(new GoodsCoverage(), new OneColumnFile(fileNameRoot+"goodscoverage")));
						outputNames.push_back(fileNameRoot+"goodscoverage"); outputTypes["goodscoverage"].push_back(fileNameRoot+"goodscoverage");
					}else if (Estimators[i] == "efron") {
						cDisplays.push_back(new CollectDisplay(new Efron(size), new OneColumnFile(fileNameRoot+"efron")));
						outputNames.push_back(fileNameRoot+"efron"); outputTypes["efron"].push_back(fileNameRoot+"efron");
					}else if (Estimators[i] == "boneh") {
						cDisplays.push_back(new CollectDisplay(new Boneh(size), new OneColumnFile(fileNameRoot+"boneh")));
						outputNames.push_back(fileNameRoot+"boneh"); outputTypes["boneh"].push_back(fileNameRoot+"boneh");
					}else if (Estimators[i] == "solow") {
						cDisplays.push_back(new CollectDisplay(new Solow(size), new OneColumnFile(fileNameRoot+"solow")));
						outputNames.push_back(fileNameRoot+"solow"); outputTypes["solow"].push_back(fileNameRoot+"solow");
					}else if (Estimators[i] == "shen") {
						cDisplays.push_back(new CollectDisplay(new Shen(size, abund), new OneColumnFile(fileNameRoot+"shen")));
						outputNames.push_back(fileNameRoot+"shen"); outputTypes["shen"].push_back(fileNameRoot+"shen");
					}
				}
			}
		
			//if the users entered no valid calculators don't execute command
			if (cDisplays.size() == 0) { return 0; }
			
			input = new InputData(inputFileNames[p], format);
			order = input->getOrderVector();
			string lastLabel = order->getLabel();
			
			//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
			set<string> processedLabels;
			set<string> userLabels = labels;
			
			if (m->control_pressed) {  
				for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}
				for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} outputTypes.clear(); 
				delete input;  
				delete order; 
				m->clearGroups();
				return 0;
			}


			while((order != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
				if (m->control_pressed) { 
					for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}
					for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} outputTypes.clear(); 
					delete input;  
					delete order; 
					m->clearGroups();
					return 0;
				}

				
				if(allLines == 1 || labels.count(order->getLabel()) == 1){
				
					m->mothurOut(order->getLabel()); m->mothurOutEndLine();
					cCurve = new Collect(order, cDisplays);
					cCurve->getCurve(freq);
					delete cCurve;
					
					processedLabels.insert(order->getLabel());
					userLabels.erase(order->getLabel());
					
					
				}
				//you have a label the user want that is smaller than this label and the last label has not already been processed 
				if ((m->anyLabelsToProcess(order->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = order->getLabel();
					
					delete order;
					order = (input->getOrderVector(lastLabel));
					
					m->mothurOut(order->getLabel()); m->mothurOutEndLine();
					cCurve = new Collect(order, cDisplays);
					cCurve->getCurve(freq);
					delete cCurve;
					
					
					processedLabels.insert(order->getLabel());
					userLabels.erase(order->getLabel());
					
					//restore real lastlabel to save below
					order->setLabel(saveLabel);
				}
				
				lastLabel = order->getLabel();	
				
				delete order;		
				order = (input->getOrderVector());
			}
			
			
			if (m->control_pressed) { 
					for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}
					for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} outputTypes.clear(); 
					delete input;  
					m->clearGroups();
					return 0;
			}
				
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
				if (order != NULL) {	delete order;	}
				order = (input->getOrderVector(lastLabel));
				
				m->mothurOut(order->getLabel()); m->mothurOutEndLine();
				
				cCurve = new Collect(order, cDisplays);
				cCurve->getCurve(freq);
				delete cCurve;
				
				if (m->control_pressed) { 
					for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}
					for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} outputTypes.clear(); 
					delete input;  
					delete order;
					m->clearGroups();
					return 0;
				}
				delete order;
			}
			
			for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}
			cDisplays.clear();
			delete input;  
		}
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} return 0; }
				
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "CollectCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
vector<string> CollectCommand::parseSharedFile(string filename) {
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

		return filenames;
	}
	catch(exception& e) {
		m->errorOut(e, "CollectCommand", "parseSharedFile");
		exit(1);
	}
}
//**********************************************************************************************************************

