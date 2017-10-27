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
#include "shannonrange.h"


//**********************************************************************************************************************
vector<string> CollectCommand::setParameters(){	
	try {
		CommandParameter plist("list", "InputTypes", "", "", "LRSS", "LRSS", "none","",false,false,true); parameters.push_back(plist);
		CommandParameter prabund("rabund", "InputTypes", "", "", "LRSS", "LRSS", "none","",false,false,true); parameters.push_back(prabund);
		CommandParameter psabund("sabund", "InputTypes", "", "", "LRSS", "LRSS", "none","",false,false,true); parameters.push_back(psabund);
		CommandParameter pshared("shared", "InputTypes", "", "", "LRSS", "LRSS", "none","",false,false,true); parameters.push_back(pshared);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pfreq("freq", "Number", "", "100", "", "", "","",false,false); parameters.push_back(pfreq);
		CommandParameter pcalc("calc", "Multiple", "sobs-chao-nseqs-coverage-ace-jack-shannon-shannoneven-npshannon-heip-smithwilson-simpson-simpsoneven-invsimpson-bootstrap-geometric-qstat-logseries-bergerparker-bstick-goodscoverage-efron-boneh-solow-shen", "sobs-chao-ace-jack-shannon-npshannon-simpson-shannonrange", "", "", "","",true,false,true); parameters.push_back(pcalc);
		CommandParameter pabund("abund", "Number", "", "10", "", "", "","",false,false); parameters.push_back(pabund);
        CommandParameter palpha("alpha", "Multiple", "0-1-2", "1", "", "", "","",false,false,true); parameters.push_back(palpha);
        CommandParameter psize("size", "Number", "", "0", "", "", "","",false,false); parameters.push_back(psize);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
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
		helpString += "The collect.single command parameters are list, sabund, rabund, shared, label, freq, calc, alpha and abund.  list, sabund, rabund or shared is required unless you have a valid current file. \n";
		helpString += "The collect.single command should be in the following format: \n";
		helpString += "The freq parameter is used indicate when to output your data, by default it is set to 100. But you can set it to a percentage of the number of sequence. For example freq=0.10, means 10%. \n";
		helpString += "collect.single(label=yourLabel, freq=yourFreq, calc=yourEstimators).\n";
		helpString += "Example collect(label=unique-.01-.03, freq=10, calc=sobs-chao-ace-jack).\n";
		helpString += "The default values for freq is 100, and calc are sobs-chao-ace-jack-shannon-npshannon-simpson.\n";
        helpString += "The alpha parameter is used to set the alpha value for the shannonrange calculator.\n";
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
string CollectCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "sobs")             {  pattern =  "[filename],sobs";            }
        else if (type == "chao")        {  pattern =  "[filename],chao";            }
        else if (type == "nseqs")       {  pattern =  "[filename],nseqs";           }
        else if (type == "coverage")    {  pattern =  "[filename],coverage";        }
        else if (type == "ace")         {  pattern =  "[filename],ace";             }
        else if (type == "jack")        {  pattern =  "[filename],jack";            }
        else if (type == "shannon")     {  pattern =  "[filename],shannon";         }
        else if (type == "shannoneven") {  pattern =  "[filename],shannoneven";     }
        else if (type == "shannonrange"){  pattern =  "[filename],shannonrange";    }
        else if (type == "npshannon")   {  pattern =  "[filename],npshannon";       }
        else if (type == "heip")        {  pattern =  "[filename],heip";            }
        else if (type == "smithwilson") {  pattern =  "[filename],smithwilson";     }
        else if (type == "simpson")     {  pattern =  "[filename],simpson";         }
        else if (type == "simpsoneven") {  pattern =  "[filename],simpsoneven";     }
        else if (type == "invsimpson")  {  pattern =  "[filename],invsimpson";      }
        else if (type == "bootstrap")   {  pattern =  "[filename],bootstrap";       }
        else if (type == "geometric")   {  pattern =  "[filename],geometric";       }
        else if (type == "qstat")       {  pattern =  "[filename],qstat";           }
        else if (type == "logseries")   {  pattern =  "[filename],logseries";       }
        else if (type == "bergerparker") {  pattern =  "[filename],bergerparker";   }
        else if (type == "bstick")      {  pattern =  "[filename],bstick";          }
        else if (type == "goodscoverage") {  pattern =  "[filename],goodscoverage"; }
        else if (type == "efron")       {  pattern =  "[filename],efron";           }
        else if (type == "boneh")       {  pattern =  "[filename],boneh";           }
        else if (type == "solow")       {  pattern =  "[filename],solow";           }
        else if (type == "shen")        {  pattern =  "[filename],shen";            }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "CollectCommand", "getOutputPattern");
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
        outputTypes["shannonrange"] = tempOutNames;
		outputTypes["npshannon"] = tempOutNames;
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
			outputTypes["npshannon"] = tempOutNames;
			outputTypes["heip"] = tempOutNames;
			outputTypes["smithwilson"] = tempOutNames;
			outputTypes["simpson"] = tempOutNames;
			outputTypes["simpsoneven"] = tempOutNames;
            outputTypes["shannonrange"] = tempOutNames;
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
            
            temp = validParameter.validFile(parameters, "alpha", false);		if (temp == "not found") { temp = "1"; }
			m->mothurConvert(temp, alpha);
            
            if ((alpha != 0) && (alpha != 1) && (alpha != 2)) { m->mothurOut("[ERROR]: Not a valid alpha value. Valid values are 0, 1 and 2."); m->mothurOutEndLine(); abort=true; }
			
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
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
	
		if ((format != "sharedfile")) { inputFileNames.push_back(inputfile);  }
		else {  inputFileNames = parseSharedFile(sharedfile);  format = "rabund"; }
	
		for (int p = 0; p < inputFileNames.size(); p++) {
			
			if (m->getControl_pressed()) {  outputTypes.clear(); for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	}    return 0; }
			
			if (outputDir == "") { outputDir += m->hasPath(inputFileNames[p]); }
			string fileNameRoot = outputDir + m->getRootName(m->getSimpleName(inputFileNames[p]));
            map<string, string> variables; 
            variables["[filename]"] = fileNameRoot;
			//globaldata->inputFileName = inputFileNames[p];
		
			if (inputFileNames.size() > 1) {
				m->mothurOutEndLine(); m->mothurOut("Processing group " + groups[p]); m->mothurOutEndLine(); m->mothurOutEndLine();
			}
		
			ValidCalculators validCalculator;
			
			for (int i=0; i<Estimators.size(); i++) {
				if (validCalculator.isValidCalculator("single", Estimators[i]) ) { 
					if (Estimators[i] == "sobs") { 
						cDisplays.push_back(new CollectDisplay(new Sobs(), new OneColumnFile(getOutputFileName("sobs", variables))));
						outputNames.push_back(getOutputFileName("sobs", variables)); outputTypes["sobs"].push_back(getOutputFileName("sobs", variables));
					}else if (Estimators[i] == "chao") { 
						cDisplays.push_back(new CollectDisplay(new Chao1(), new ThreeColumnFile(getOutputFileName("chao", variables))));
						outputNames.push_back(getOutputFileName("chao", variables)); outputTypes["chao"].push_back(getOutputFileName("chao", variables));
					}else if (Estimators[i] == "nseqs") { 
						cDisplays.push_back(new CollectDisplay(new NSeqs(), new OneColumnFile(getOutputFileName("nseqs", variables))));
						outputNames.push_back(getOutputFileName("nseqs", variables)); outputTypes["nseqs"].push_back(getOutputFileName("nseqs", variables));
					}else if (Estimators[i] == "coverage") { 
						cDisplays.push_back(new CollectDisplay(new Coverage(), new OneColumnFile(getOutputFileName("coverage", variables))));
						outputNames.push_back(getOutputFileName("coverage", variables)); outputTypes["coverage"].push_back(getOutputFileName("coverage", variables));
					}else if (Estimators[i] == "ace") { 
						cDisplays.push_back(new CollectDisplay(new Ace(abund), new ThreeColumnFile(getOutputFileName("ace", variables))));
						outputNames.push_back(getOutputFileName("ace", variables)); outputTypes["ace"].push_back(getOutputFileName("ace", variables));
					}else if (Estimators[i] == "jack") { 
						cDisplays.push_back(new CollectDisplay(new Jackknife(), new ThreeColumnFile(getOutputFileName("jack", variables))));
						outputNames.push_back(getOutputFileName("jack", variables)); outputTypes["jack"].push_back(getOutputFileName("jack", variables));
					}else if (Estimators[i] == "shannon") { 
						cDisplays.push_back(new CollectDisplay(new Shannon(), new ThreeColumnFile(getOutputFileName("shannon", variables))));
						outputNames.push_back(getOutputFileName("shannon", variables)); outputTypes["shannon"].push_back(getOutputFileName("shannon", variables));
					}else if (Estimators[i] == "shannoneven") { 
						cDisplays.push_back(new CollectDisplay(new ShannonEven(), new OneColumnFile(getOutputFileName("shannoneven", variables))));
						outputNames.push_back(getOutputFileName("shannoneven", variables)); outputTypes["shannoneven"].push_back(getOutputFileName("shannoneven", variables));
                    }else if (Estimators[i] == "shannonrange") {
                            cDisplays.push_back(new CollectDisplay(new RangeShannon(alpha), new ThreeColumnFile(getOutputFileName("shannonrange", variables))));
                            outputNames.push_back(getOutputFileName("shannonrange", variables)); outputTypes["shannoneven"].push_back(getOutputFileName("shannonrange", variables));
					}else if (Estimators[i] == "npshannon") { 
						cDisplays.push_back(new CollectDisplay(new NPShannon(), new OneColumnFile(getOutputFileName("npshannon", variables))));
						outputNames.push_back(getOutputFileName("npshannon", variables)); outputTypes["npshannon"].push_back(getOutputFileName("npshannon", variables));
					}else if (Estimators[i] == "heip") { 
						cDisplays.push_back(new CollectDisplay(new Heip(), new OneColumnFile(getOutputFileName("heip", variables))));
						outputNames.push_back(getOutputFileName("heip", variables)); outputTypes["heip"].push_back(getOutputFileName("heip", variables));
					}else if (Estimators[i] == "smithwilson") { 
						cDisplays.push_back(new CollectDisplay(new SmithWilson(), new OneColumnFile(getOutputFileName("smithwilson", variables))));
						outputNames.push_back(getOutputFileName("smithwilson", variables)); outputTypes["smithwilson"].push_back(getOutputFileName("smithwilson", variables));
					}else if (Estimators[i] == "simpson") { 
						cDisplays.push_back(new CollectDisplay(new Simpson(), new ThreeColumnFile(getOutputFileName("simpson", variables))));
						outputNames.push_back(getOutputFileName("simpson", variables)); outputTypes["simpson"].push_back(getOutputFileName("simpson", variables));
					}else if (Estimators[i] == "simpsoneven") { 
						cDisplays.push_back(new CollectDisplay(new SimpsonEven(), new OneColumnFile(getOutputFileName("simpsoneven", variables))));
						outputNames.push_back(getOutputFileName("simpsoneven", variables)); outputTypes["simpsoneven"].push_back(getOutputFileName("simpsoneven", variables));
					}else if (Estimators[i] == "invsimpson") { 
						cDisplays.push_back(new CollectDisplay(new InvSimpson(), new ThreeColumnFile(getOutputFileName("invsimpson", variables))));
						outputNames.push_back(getOutputFileName("invsimpson", variables)); outputTypes["invsimpson"].push_back(getOutputFileName("invsimpson", variables));
					}else if (Estimators[i] == "bootstrap") { 
						cDisplays.push_back(new CollectDisplay(new Bootstrap(), new OneColumnFile(getOutputFileName("bootstrap", variables))));
						outputNames.push_back(getOutputFileName("bootstrap", variables)); outputTypes["bootstrap"].push_back(getOutputFileName("bootstrap", variables));
					}else if (Estimators[i] == "geometric") { 
						cDisplays.push_back(new CollectDisplay(new Geom(), new OneColumnFile(getOutputFileName("geometric", variables))));
						outputNames.push_back(getOutputFileName("geometric", variables)); outputTypes["geometric"].push_back(getOutputFileName("geometric", variables));
					}else if (Estimators[i] == "qstat") { 
						cDisplays.push_back(new CollectDisplay(new QStat(), new OneColumnFile(getOutputFileName("qstat", variables))));
						outputNames.push_back(getOutputFileName("qstat", variables)); outputTypes["qstat"].push_back(getOutputFileName("qstat", variables));
					}else if (Estimators[i] == "logseries") { 
						cDisplays.push_back(new CollectDisplay(new LogSD(), new OneColumnFile(getOutputFileName("logseries", variables))));
						outputNames.push_back(getOutputFileName("logseries", variables)); outputTypes["logseries"].push_back(getOutputFileName("logseries", variables));
					}else if (Estimators[i] == "bergerparker") { 
						cDisplays.push_back(new CollectDisplay(new BergerParker(), new OneColumnFile(getOutputFileName("bergerparker", variables))));
						outputNames.push_back(getOutputFileName("bergerparker", variables)); outputTypes["bergerparker"].push_back(getOutputFileName("bergerparker", variables));
					}else if (Estimators[i] == "bstick") { 
						cDisplays.push_back(new CollectDisplay(new BStick(), new ThreeColumnFile(getOutputFileName("bstick", variables))));
						outputNames.push_back(getOutputFileName("bstick", variables)); outputTypes["bstick"].push_back(getOutputFileName("bstick", variables));
					}else if (Estimators[i] == "goodscoverage") { 
						cDisplays.push_back(new CollectDisplay(new GoodsCoverage(), new OneColumnFile(getOutputFileName("goodscoverage", variables))));
						outputNames.push_back(getOutputFileName("goodscoverage", variables)); outputTypes["goodscoverage"].push_back(getOutputFileName("goodscoverage", variables));
					}else if (Estimators[i] == "efron") {
						cDisplays.push_back(new CollectDisplay(new Efron(size), new OneColumnFile(getOutputFileName("efron", variables))));
						outputNames.push_back(getOutputFileName("efron", variables)); outputTypes["efron"].push_back(getOutputFileName("efron", variables));
					}else if (Estimators[i] == "boneh") {
						cDisplays.push_back(new CollectDisplay(new Boneh(size), new OneColumnFile(getOutputFileName("boneh", variables))));
						outputNames.push_back(getOutputFileName("boneh", variables)); outputTypes["boneh"].push_back(getOutputFileName("boneh", variables));
					}else if (Estimators[i] == "solow") {
						cDisplays.push_back(new CollectDisplay(new Solow(size), new OneColumnFile(getOutputFileName("solow", variables))));
						outputNames.push_back(getOutputFileName("solow", variables)); outputTypes["solow"].push_back(getOutputFileName("solow", variables));
					}else if (Estimators[i] == "shen") {
						cDisplays.push_back(new CollectDisplay(new Shen(size, abund), new OneColumnFile(getOutputFileName("shen", variables))));
						outputNames.push_back(getOutputFileName("shen", variables)); outputTypes["shen"].push_back(getOutputFileName("shen", variables));
					}
				}
			}
		
			//if the users entered no valid calculators don't execute command
			if (cDisplays.size() == 0) { return 0; }
			
			input = new InputData(inputFileNames[p], format, nullVector);
			order = input->getOrderVector();
			string lastLabel = order->getLabel();
			
			//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
			set<string> processedLabels;
			set<string> userLabels = labels;
			
			if (m->getControl_pressed()) {  
				for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}
				for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} outputTypes.clear(); 
				delete input;  
				delete order; 
				
				return 0;
			}


			while((order != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
				if (m->getControl_pressed()) { 
					for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}
					for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} outputTypes.clear(); 
					delete input;  
					delete order; 
					
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
				if ((m->anyLabelsToProcess(order->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
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
			
			
			if (m->getControl_pressed()) { 
					for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}
					for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} outputTypes.clear(); 
					delete input;  
					
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
			if (needToRun )  {
				if (order != NULL) {	delete order;	}
				order = (input->getOrderVector(lastLabel));
				
				m->mothurOut(order->getLabel()); m->mothurOutEndLine();
				
				cCurve = new Collect(order, cDisplays);
				cCurve->getCurve(freq);
				delete cCurve;
				
				if (m->getControl_pressed()) { 
					for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}
					for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} outputTypes.clear(); 
					delete input;  
					delete order;
					
					return 0;
				}
				delete order;
			}
			
			for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}
			cDisplays.clear();
			delete input;  
		}
		
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} return 0; }
				
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
		
		map<string, string> files;
		map<string, string>::iterator it3;
					
		input = new InputData(filename, "sharedfile", groups);
		SharedRAbundVectors* shared = input->getSharedRAbundVectors();
		
		string sharedFileRoot = m->getRootName(filename);
        groups = shared->getNamesGroups();
		
		//clears file before we start to write to it below
		for (int i=0; i<groups.size(); i++) {
            ofstream temp;
            string group = groups[i];
			m->openOutputFile((sharedFileRoot + group + ".rabund"), temp);
            temp.close();
			filenames.push_back((sharedFileRoot + group + ".rabund"));
            files[group] = (sharedFileRoot + group + ".rabund");
		}
		
		while(shared != NULL) {
            
            vector<SharedRAbundVector*> lookup = shared->getSharedRAbundVectors();
			for (int i = 0; i < lookup.size(); i++) {
                ofstream temp;
                string group = groups[i];
				m->openOutputFileAppend(files[group], temp);
				lookup[i]->getRAbundVector().print(temp);
				temp.close();
			}
		
            for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) { delete lookup[i]; } lookup[i] = NULL; }
			shared = input->getSharedRAbundVectors();
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

