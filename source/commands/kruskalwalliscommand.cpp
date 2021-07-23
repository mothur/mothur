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
		
        abort = false; calledHelp = false; allLines = true;
        
        vector<string> tempOutNames;
        outputTypes["kruskall-wallis"] = tempOutNames;
        
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
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "KruskalWallisCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
KruskalWallisCommand::KruskalWallisCommand(string option)  {
	try {

		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }
			else if (sharedfile == "not found") {
				//if there is a current shared file, use it
				sharedfile = current->getSharedFile();
				if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n");  }
				else { 	m->mothurOut("You have no current sharedfile and the shared parameter is required.\n");  abort = true; }
			}else { current->setSharedFile(sharedfile); }
            
            //get shared file, it is required
			designfile = validParameter.validFile(parameters, "design");
			if (designfile == "not open") { designfile = ""; abort = true; }
			else if (designfile == "not found") {
				//if there is a current shared file, use it
				designfile = current->getDesignFile();
				if (designfile != "") { m->mothurOut("Using " + designfile + " as input file for the design parameter.\n");  }
				else { 	m->mothurOut("You have no current design file and the design parameter is required.\n");  abort = true; }
			}else { current->setDesignFile(designfile); }
            
            
            if (outputdir == ""){
				outputdir = util.hasPath(sharedfile); 
			}
            
            string label = validParameter.valid(parameters, "label");
			if (label == "not found") { label = ""; }
			else {
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				else { allLines = true;  }
			}
            
            mclass = validParameter.valid(parameters, "class");
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
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        DesignMap designMap(designfile); if (m->getControl_pressed()) { return 0; }
        
        //if user did not select class use first column
        if (mclass == "") {  mclass = designMap.getDefaultClass(); m->mothurOut("\nYou did not provide a class, using " + mclass +".\n\n"); }
        
        InputData input(sharedfile, "sharedfile", nullVector);
        set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
        SharedRAbundVectors* lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
        vector<string> currentLabels = lookup->getOTUNames();
        
        while (lookup != NULL) {
            
            if (m->getControl_pressed()) { delete lookup; break; }
            
            vector<SharedRAbundVector*> data = lookup->getSharedRAbundVectors();
            process(data, designMap, currentLabels);
            for (int i = 0; i < data.size(); i++) { delete data[i]; } data.clear();
            delete lookup;
            
            lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
        }
        
        //output files created by command
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
        return 0;
		
    }
	catch(exception& e) {
		m->errorOut(e, "KruskalWallisCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int KruskalWallisCommand::process(vector<SharedRAbundVector*>& lookup, DesignMap& designMap, vector<string> currentLabels) {
	try {
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(sharedfile));
        variables["[distance]"] = lookup[0]->getLabel();
		string outputFileName = getOutputFileName("kruskall-wallis",variables);
        
		ofstream out;
		util.openOutputFile(outputFileName, out);
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
        if (treatments.size() < 2) { m->mothurOut("[ERROR]: need at least 2 things for classes to compare, quitting.\n"); m->setControl_pressed(true); }
        
        LinearAlgebra linear;
        for (int i = 0; i < numBins; i++) {
            if (m->getControl_pressed()) { break; }
            
            vector<spearmanRank> values;
            for (int j = 0; j < lookup.size(); j++) {
                string group = lookup[j]->getGroup();
                string treatment = designMap.get(group, mclass); //get value for this group in this category
                spearmanRank temp(treatment, lookup[j]->get(i));
                values.push_back(temp);
            }
            
            double pValue = 0.0;
            double H = linear.calcKruskalWallis(values, pValue);
            
            //output H and signifigance
            out << currentLabels[i] << '\t' << H << '\t' << pValue << endl;
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


