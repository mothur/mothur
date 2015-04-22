/* 
 * File:   kruskalwalliscommand.cpp
 * Author: kiverson
 *
 * Created on June 26, 2012, 11:06 AM
 */

#include "kruskalwalliscommand.h"
#include "linearalgebra.h"

//**********************************************************************************************************************
vector<string> KruskalWallisCommand::setParameters(){
	try {
        CommandParameter pdesign("design", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(pdesign);
        CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none","summary",false,true,true); parameters.push_back(pshared);
        CommandParameter pclass("class", "String", "", "", "", "", "","",false,false); parameters.push_back(pclass);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
        //every command must have inputdir and outputdir.  This allows mothur users to redirect input and output files.
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "KruskalWallisCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string KruskalWallisCommand::getHelpString(){
	try {
		string helpString = "";
		helpString += "The kruskal.wallis command allows you to ....\n";
		helpString += "The kruskal.wallis command parameters are: shared, design, class, label and classes.\n";
		helpString += "The class parameter is used to indicate the which category you would like used for the Kruskal Wallis analysis. If none is provided first category is used.\n";
        helpString += "The label parameter is used to indicate which distances in the shared file you would like to use. labels are separated by dashes.\n";
		helpString += "The kruskal.wallis command should be in the following format: kruskal.wallis(shared=final.an.shared, design=final.design, class=treatment).\n";
        return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "KruskalWallisCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string KruskalWallisCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "kruskall-wallis") {  pattern = "[filename],[distance],kruskall_wallis"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "KruskalWallisCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
KruskalWallisCommand::KruskalWallisCommand(){
	try {
		abort = true; calledHelp = true;
		setParameters();
        vector<string> tempOutNames;
        outputTypes["kruskall-wallis"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "KruskalWallisCommand", "KruskalWallisCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
KruskalWallisCommand::KruskalWallisCommand(string option)  {
	try {
		abort = false; calledHelp = false;
        allLines = 1;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			//valid paramters for this command
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) {
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			vector<string> tempOutNames;
            outputTypes["kruskall-wallis"] = tempOutNames;
            
			//if the user changes the input directory command factory will send this info to us in the output parameter
			string inputDir = validParameter.validFile(parameters, "inputdir", false);
			if (inputDir == "not found"){	inputDir = "";		}
			else {
                
                string path;
				it = parameters.find("design");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["desing"] = inputDir + it->second;		}
				}
				
                it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
            }
            
            //get shared file, it is required
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }
			else if (sharedfile == "not found") {
				//if there is a current shared file, use it
				sharedfile = m->getSharedFile();
				if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current sharedfile and the shared parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setSharedFile(sharedfile); }
            
            //get shared file, it is required
			designfile = validParameter.validFile(parameters, "design", true);
			if (designfile == "not open") { designfile = ""; abort = true; }
			else if (designfile == "not found") {
				//if there is a current shared file, use it
				designfile = m->getDesignFile();
				if (designfile != "") { m->mothurOut("Using " + designfile + " as input file for the design parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current design file and the design parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setDesignFile(designfile); }
            
            //if the user changes the output directory command factory will send this info to us in the output parameter
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){
				outputDir = m->hasPath(sharedfile); //if user entered a file with a path then preserve it
			}
            
            string label = validParameter.validFile(parameters, "label", false);
			if (label == "not found") { label = ""; }
			else {
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
            
            mclass = validParameter.validFile(parameters, "class", false);
			if (mclass == "not found") { mclass = ""; }
            
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "KruskalWallisCommand", "KruskalWallisCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int KruskalWallisCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        DesignMap designMap(designfile);
        
        //if user did not select class use first column
        if (mclass == "") {  mclass = designMap.getDefaultClass(); m->mothurOut("\nYou did not provide a class, using " + mclass +".\n\n"); }
        
        InputData input(sharedfile, "sharedfile");
        vector<SharedRAbundVector*> lookup = input.getSharedRAbundVectors();
        string lastLabel = lookup[0]->getLabel();
        
        //if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
        set<string> processedLabels;
        set<string> userLabels = labels;
        
        
        //as long as you are not at the end of the file or done wih the lines you want
        while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
            
            if (m->control_pressed) { for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  return 0; }
            
            if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){
                
                m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
                
                process(lookup, designMap);
                
                processedLabels.insert(lookup[0]->getLabel());
                userLabels.erase(lookup[0]->getLabel());
            }
            
            if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
                string saveLabel = lookup[0]->getLabel();
                
                for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
                lookup = input.getSharedRAbundVectors(lastLabel);
                m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
                
                process(lookup, designMap);
                
                processedLabels.insert(lookup[0]->getLabel());
                userLabels.erase(lookup[0]->getLabel());
                
                //restore real lastlabel to save below
                lookup[0]->setLabel(saveLabel);
            }
            
            lastLabel = lookup[0]->getLabel();
            //prevent memory leak
            for (int i = 0; i < lookup.size(); i++) {  delete lookup[i]; lookup[i] = NULL; }
            
            if (m->control_pressed) { return 0; }
            
            //get next line to process
            lookup = input.getSharedRAbundVectors();
        }
        
        if (m->control_pressed) {  return 0; }
        
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
            for (int i = 0; i < lookup.size(); i++) { if (lookup[i] != NULL) { delete lookup[i]; } }
            lookup = input.getSharedRAbundVectors(lastLabel);
            
            m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
            process(lookup, designMap);
            
            for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
        }
        
		
        //output files created by command
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
        return 0;
		
    }
	catch(exception& e) {
		m->errorOut(e, "KruskalWallisCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int KruskalWallisCommand::process(vector<SharedRAbundVector*>& lookup, DesignMap& designMap) {
	try {
        map<string, string> variables;
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(sharedfile));
        variables["[distance]"] = lookup[0]->getLabel();
		string outputFileName = getOutputFileName("kruskall-wallis",variables);
        
		ofstream out;
		m->openOutputFile(outputFileName, out);
		outputNames.push_back(outputFileName); outputTypes["kruskall-wallis"].push_back(outputFileName);
        out << "OTULabel\tKW\tPvalue\n";
        
        int numBins = lookup[0]->getNumBins();
        //sanity check to make sure each treatment has a group in the shared file
        set<string> treatments;
        for (int j = 0; j < lookup.size(); j++) {
            string group = lookup[j]->getGroup();
            string treatment = designMap.get(group, mclass); //get value for this group in this category
            treatments.insert(treatment);
        }
        if (treatments.size() < 2) { m->mothurOut("[ERROR]: need at least 2 things to classes to compare, quitting.\n"); m->control_pressed = true; }
        
        LinearAlgebra linear;
        for (int i = 0; i < numBins; i++) {
            if (m->control_pressed) { break; }
            
            vector<spearmanRank> values;
            for (int j = 0; j < lookup.size(); j++) {
                string group = lookup[j]->getGroup();
                string treatment = designMap.get(group, mclass); //get value for this group in this category
                spearmanRank temp(treatment, lookup[j]->getAbundance(i));
                values.push_back(temp);
            }
            
            double pValue = 0.0;
            double H = linear.calcKruskalWallis(values, pValue);
            
            //output H and signifigance
            out << m->currentSharedBinLabels[i] << '\t' << H << '\t' << pValue << endl;
        }
        out.close();
                
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "KruskalWallisCommand", "process");
		exit(1);
	}
}
//**********************************************************************************************************************


