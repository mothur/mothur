/*
 *  sensspeccommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 7/6/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "sensspeccommand.h"
#include "calculator.h"
#include "mcc.hpp"
#include "sensitivity.hpp"
#include "specificity.hpp"
#include "fdr.hpp"
#include "npv.hpp"
#include "ppv.hpp"
#include "f1score.hpp"
#include "tp.hpp"
#include "fp.hpp"
#include "fpfn.hpp"
#include "tptn.hpp"
#include "tn.hpp"
#include "fn.hpp"
#include "accuracy.hpp"

//**********************************************************************************************************************
vector<string> SensSpecCommand::setParameters(){
	try {
		CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none","sensspec",false,true,true); parameters.push_back(plist);
		CommandParameter pphylip("phylip", "InputTypes", "", "", "PhylipColumn", "PhylipColumn", "none","",false,false); parameters.push_back(pphylip);
		CommandParameter pcolumn("column", "InputTypes", "", "", "PhylipColumn", "PhylipColumn", "none","",false,false); parameters.push_back(pcolumn);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount", "none", "none","name",false,false,true); parameters.push_back(pcount);
		CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","name",false,false,true); parameters.push_back(pname);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pcutoff("cutoff", "Number", "", "-1.00", "", "", "","",false,false); parameters.push_back(pcutoff);
		CommandParameter pprecision("precision", "Number", "", "100", "", "", "","",false,false); parameters.push_back(pprecision);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);

		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SensSpecCommand::getHelpString(){
	try {
		string helpString = "";
		helpString += "The sens.spec command determines the quality of the clusters.\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string SensSpecCommand::getOutputPattern(string type) {
    try {
        string pattern = "";

        if (type == "sensspec") {  pattern = "[filename],sensspec"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }

        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SensSpecCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
SensSpecCommand::SensSpecCommand(){
	try {
		abort = true; calledHelp = true;
		setParameters();
		vector<string> tempOutNames;
		outputTypes["sensspec"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "SensSpecCommand");
		exit(1);
	}
}
//***************************************************************************************************************

SensSpecCommand::SensSpecCommand(string option)  {
	try {

		abort = false; calledHelp = false;
		allLines = 1;

		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}

		else {
			string temp;

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
			outputTypes["sensspec"] = tempOutNames;

			//if the user changes the input directory command factory will send this info to us in the output parameter
			string inputDir = validParameter.valid(parameters, "inputdir");
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}

				it = parameters.find("phylip");
				//user has given a template file
				if(it != parameters.end()){
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["phylip"] = inputDir + it->second;		}
				}

				it = parameters.find("column");
				//user has given a template file
				if(it != parameters.end()){
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["column"] = inputDir + it->second;		}
				}

				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["count"] = inputDir + it->second;		}
                }

			}
			//check for required parameters
			listFile = validParameter.validFile(parameters, "list");
			if (listFile == "not found") {
				listFile = current->getListFile();
				if (listFile != "") { m->mothurOut("Using " + listFile + " as input file for the list parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current list file and the list parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (listFile == "not open") { abort = true; }
			else { current->setListFile(listFile); }

			phylipfile = validParameter.validFile(parameters, "phylip");
			if (phylipfile == "not found") { phylipfile = "";  }
			else if (phylipfile == "not open") { abort = true; }
			else { distFile = phylipfile; format = "phylip"; current->setPhylipFile(phylipfile);  }

			columnfile = validParameter.validFile(parameters, "column");
			if (columnfile == "not found") { columnfile = ""; }
			else if (columnfile == "not open") { abort = true; }
			else { distFile = columnfile; format = "column";   current->setColumnFile(columnfile); }

			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not found") { namefile =  "";  }
			else if (namefile == "not open") { namefile = ""; abort = true; }
			else {  current->setNameFile(namefile); }
            
            countfile = validParameter.validFile(parameters, "count");
            if (countfile == "not found") { countfile =  "";  }
            else if (countfile == "not open") { countfile = ""; abort = true; }
            else {  current->setCountFile(countfile); }


			if ((phylipfile == "") && (columnfile == "")) { //is there are current file available for either of these?
				//give priority to column, then phylip
				columnfile = current->getColumnFile();
				if (columnfile != "") {  distFile = columnfile; format = "column";  m->mothurOut("Using " + columnfile + " as input file for the column parameter."); m->mothurOutEndLine(); }
				else {
					phylipfile = current->getPhylipFile();
					if (phylipfile != "") {  distFile = phylipfile; format = "phylip"; m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter."); m->mothurOutEndLine(); }
					else {
						m->mothurOut("No valid current files. You must provide a phylip or column file."); m->mothurOutEndLine();
						abort = true;
					}
				}
			}else if ((phylipfile != "") && (columnfile != "")) { m->mothurOut("When executing a sens.spec command you must enter ONLY ONE of the following: phylip or column."); m->mothurOutEndLine(); abort = true; }

            if (columnfile != "") {
                if ((namefile == "") && (countfile == "")){
                    namefile = current->getNameFile();
                    if (namefile != "") {  m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
                    else {
                        countfile = current->getCountFile();
                        if (countfile != "") {  m->mothurOut("Using " + countfile + " as input file for the count parameter."); m->mothurOutEndLine(); }
                        else {
                            m->mothurOut("You need to provide a namefile or countfile if you are going to use the column format."); m->mothurOutEndLine();
                            abort = true;
                        }	
                    }	
                }
            }

			if ((namefile == "") && (phylipfile != "")) {
                m->mothurOut("[WARNING]: there is no reason to include a name file with a phylip file. Ignoring..."); m->mothurOutEndLine(); abort = false;
            }

			//if the user changes the output directory command factory will send this info to us in the output parameter
			outputDir = validParameter.valid(parameters, "outputdir");
			if (outputDir == "not found"){
				outputDir = "";
				outputDir += util.hasPath(listFile); //if user entered a file with a path then preserve it
			}

			temp = validParameter.valid(parameters, "cutoff");		if (temp == "not found") { temp = "-1.00"; }
			util.mothurConvert(temp, cutoff);

			temp = validParameter.valid(parameters, "precision");	if (temp == "not found") { temp = "100"; }
			util.mothurConvert(temp, precision);

			string label = validParameter.valid(parameters, "label");
			if (label == "not found") { label = ""; }
			else {
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}

            map<string, string> variables;
            variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(listFile));
			sensSpecFileName = getOutputFileName("sensspec",variables);
		}

		m->mothurOutEndLine();
		m->mothurOut("NOTE: sens.spec assumes that only unique sequences were used to generate the distance matrix.");
		m->mothurOutEndLine();
		m->mothurOutEndLine();
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "SensSpecCommand");
		exit(1);
	}
}
//***************************************************************************************************************

int SensSpecCommand::execute(){
	try{
		if (abort) { if (calledHelp) { return 0; }  return 2;	}

        int startTime = time(NULL);

		processListFile();

        if (m->getControl_pressed()) { util.mothurRemove(sensSpecFileName); return 0; }

        m->mothurOut("It took " + toString(time(NULL) - startTime) + " to run sens.spec.\n");

		m->mothurOut("\nOutput File Names: \n");
		m->mothurOut(sensSpecFileName+"\n\n");

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "execute");
		exit(1);
	}
}
//***************************************************************************************************************
int SensSpecCommand::process(ListVector*& list, bool& getCutoff, string& origCutoff){
	try {
        string label = list->getLabel();
        
        if(getCutoff == 1){
            if(label != "unique"){
                origCutoff = label;
                convert(label, cutoff);
                cutoff = util.ceilDist(cutoff, precision);
                origCutoff = toString(util.ceilDist(cutoff, precision));
            }else{ origCutoff = "unique"; cutoff = 0.0000; }
        }
        
        //must read each time because cutoff changes
        string nameOrCount = "";
        string thisNamefile = "";
        map<string, int> counts;
        if (countfile != "") { nameOrCount = "count"; thisNamefile = countfile; CountTable ct; ct.readTable(countfile, false, false); counts = ct.getNameMap(); }
        else if (namefile != "") { nameOrCount = "name"; thisNamefile = namefile; }
        
        string distfile = columnfile;
        if (format == "phylip") { distfile = phylipfile; }
        
        OptiMatrix matrix(distfile, thisNamefile, nameOrCount, format, cutoff, false);
        SensSpecCalc senscalc(matrix, list);
        senscalc.getResults(matrix, truePositives, trueNegatives, falsePositives, falseNegatives);
		outputStatistics(label, origCutoff, list->getNumBins());

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "process");
		exit(1);
	}
}
//***************************************************************************************************************
int SensSpecCommand::processListFile(){
	try{
        setUpOutput();
        
        bool getCutoff = 0;
        string origCutoff = "";
        if(cutoff == -1.00)	{	getCutoff = 1;                                              }
        else 				{	origCutoff = toString(util.ceilDist(cutoff, precision));	}
        
		InputData input(listFile, "list", nullVector);
		ListVector* list = input.getListVector();
		string lastLabel = list->getLabel();

		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;

		while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {

			if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);  }  delete list;  return 0;  }

			if(allLines == 1 || labels.count(list->getLabel()) == 1){
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());

				process(list, getCutoff, origCutoff);
			}

			if ((util.anyLabelsToProcess(list->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {

				string saveLabel = list->getLabel();

				delete list;
				list = input.getListVector(lastLabel);

				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());

				process(list, getCutoff, origCutoff);

				//restore real lastlabel to save below
				list->setLabel(saveLabel);
			}

			lastLabel = list->getLabel();
			delete list;
			list = input.getListVector();
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
			if (list != NULL) {	delete list;	}
			list = input.getListVector(lastLabel);

            process(list, getCutoff, origCutoff);

			delete list;
		}

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "processListFile");
		exit(1);
	}
}
//***************************************************************************************************************
void SensSpecCommand::setUpOutput(){
	try{
		ofstream sensSpecFile;
		util.openOutputFile(sensSpecFileName, sensSpecFile);
        outputNames.push_back(sensSpecFileName); outputTypes["sensspec"].push_back(sensSpecFileName);

		sensSpecFile << "label\tcutoff\tnumotus\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";
        m->mothurOut("label\tcutoff\tnumotus\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");

		sensSpecFile.close();
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "setUpOutput");
		exit(1);
	}
}
//***************************************************************************************************************
void SensSpecCommand::outputStatistics(string label, string cutoff, int numBins){
	try{
		long long tp =  truePositives;
		long long fp =  falsePositives;
		long long tn =  trueNegatives;
		long long fn =  falseNegatives;

        Sensitivity sens;   double sensitivity = sens.getValue(tp, tn, fp, fn);
        Specificity spec;   double specificity = spec.getValue(tp, tn, fp, fn);
        PPV ppv;            double positivePredictiveValue = ppv.getValue(tp, tn, fp, fn);
        NPV npv;            double negativePredictiveValue = npv.getValue(tp, tn, fp, fn);
        FDR fdr;            double falseDiscoveryRate = fdr.getValue(tp, tn, fp, fn);
        Accuracy acc;       double accuracy = acc.getValue(tp, tn, fp, fn);
        MCC mcc;            double matthewsCorrCoef = mcc.getValue(tp, tn, fp, fn);
        F1Score f1;         double f1Score = f1.getValue(tp, tn, fp, fn); 

		ofstream sensSpecFile;
		util.openOutputFileAppend(sensSpecFileName, sensSpecFile);

		sensSpecFile << label << '\t' << cutoff << '\t' << numBins << '\t';
		sensSpecFile << truePositives << '\t' << trueNegatives << '\t' << falsePositives << '\t' << falseNegatives << '\t';
		sensSpecFile << setprecision(4);
		sensSpecFile << sensitivity << '\t' << specificity << '\t' << positivePredictiveValue << '\t' << negativePredictiveValue << '\t';
		sensSpecFile << falseDiscoveryRate << '\t' << accuracy << '\t' << matthewsCorrCoef << '\t' << f1Score << endl;
        
        m->mothurOut(label + "\t" + toString(cutoff) + "\t" + toString(numBins) + "\t"+ toString(truePositives) + "\t" + toString(trueNegatives) + "\t" + toString(falsePositives) + "\t" + toString(falseNegatives) + "\t");
        m->mothurOut(toString(sensitivity) + "\t" + toString(specificity) + "\t" + toString(positivePredictiveValue) + "\t" + toString(negativePredictiveValue) + "\t");
        m->mothurOut(toString(falseDiscoveryRate) + "\t" + toString(accuracy) + "\t" + toString(matthewsCorrCoef) + "\t" + toString(f1Score) + "\n\n");

		sensSpecFile.close();
	}
	catch(exception& e) {
		m->errorOut(e, "SensSpecCommand", "outputStatistics");
		exit(1);
	}
}
//***************************************************************************************************************


