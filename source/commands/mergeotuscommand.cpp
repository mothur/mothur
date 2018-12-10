//
//  mergeotuscommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 12/10/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "mergeotuscommand.hpp"


//**********************************************************************************************************************
vector<string> MergeOTUsCommand::setParameters(){
    try {
        CommandParameter pconstaxonomy("constaxonomy", "InputTypes", "", "", "none", "none", "none","constaxonomy",false,true, true); parameters.push_back(pconstaxonomy);
        CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none","shared",false,true,true); parameters.push_back(pshared);
        CommandParameter prelabund("relabund", "InputTypes", "", "", "none", "none", "none","relabund",false,true,true); parameters.push_back(prelabund);
        CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none","list",false,true,true); parameters.push_back(plist);
        CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
        CommandParameter ptaxlevel("taxlevel", "Number", "", "-1", "", "", "","",false,false,true); parameters.push_back(ptaxlevel);
        
        
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        vector<string> myArray;
        for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
        return myArray;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeOTUsCommand", "setParameters");
        exit(1);
    }
}
//**********************************************************************************************************************
string MergeOTUsCommand::getHelpString(){
    try {
        string helpString = "";
        helpString += "The merge.otus command parameters are shared, list, relabund, constaxonomy, taxlevel and label.  constaxonomy is a required, unless you have a valid current file.\n";
        helpString += "The taxlevel parameter allows you to specify the taxonomy level you would like to use when merging. Default=maxlevel.\n";
        helpString += "Example merge.otus(shared=yourSharedFile, constaxonomy=yourConsTaxonomyFile).\n";
        
        return helpString;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeOTUsCommand", "getHelpString");
        exit(1);
    }
}
//**********************************************************************************************************************
string MergeOTUsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "shared")           {  pattern = "[filename],merge,[extension]"; }
        else if (type == "list")        {  pattern = "[filename],merge,[extension]"; }
        else if (type == "relabund")    {  pattern = "[filename],merge,[extension]"; }
        else if (type == "constaxonomy") {  pattern = "[filename],merge,[extension]"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeOTUsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
MergeOTUsCommand::MergeOTUsCommand(){
    try {
        abort = true; calledHelp = true;
        setParameters();
        vector<string> tempOutNames;
        outputTypes["shared"] = tempOutNames;
        outputTypes["list"] = tempOutNames;
        outputTypes["relabund"] = tempOutNames;
        outputTypes["constaxonomy"] = tempOutNames;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeOTUsCommand", "MergeOTUsCommand");
        exit(1);
    }
}
//**********************************************************************************************************************

MergeOTUsCommand::MergeOTUsCommand(string option)  {
    try {
        abort = false; calledHelp = false;
        allLines = 1;
        
        //allow user to run help
        if(option == "help") {  help(); abort = true; calledHelp = true; }
        else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        
        else {
            vector<string> myArray = setParameters();
            
            OptionParser parser(option);
            map<string,string> parameters  = parser.getParameters();
            map<string,string>::iterator it;
            
            ValidParameters validParameter;
            
            //check to make sure all parameters are valid for command
            for (it = parameters.begin(); it != parameters.end(); it++) {
                if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
            }
            
            vector<string> tempOutNames;
            outputTypes["shared"] = tempOutNames;
            outputTypes["list"] = tempOutNames;
            outputTypes["relabund"] = tempOutNames;
            outputTypes["constaxonomy"] = tempOutNames;
            
            //if the user changes the input directory command factory will send this info to us in the output parameter
            string inputDir = validParameter.valid(parameters, "inputdir");
            if (inputDir == "not found"){	inputDir = "";		}
            else {
                string path;
                it = parameters.find("shared");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["shared"] = inputDir + it->second;		}
                }
                
                it = parameters.find("constaxonomy");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["constaxonomy"] = inputDir + it->second;		}
                }
                
                it = parameters.find("relabund");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["relabund"] = inputDir + it->second;		}
                }
                
                it = parameters.find("list");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["list"] = inputDir + it->second;		}
                }
            }
            
            sharedfile = validParameter.validFile(parameters, "shared");
            if (sharedfile == "not found") { sharedfile = ""; }
            else if (sharedfile == "not open") { sharedfile = ""; abort = true; }
            else { current->setSharedFile(sharedfile); }
            
            listfile = validParameter.validFile(parameters, "list");
            if (listfile == "not found") { listfile = ""; }
            else if (listfile == "not open") { listfile = ""; abort = true; }
            else { current->setListFile(listfile); }
            
            relabundfile = validParameter.validFile(parameters, "relabund");
            if (relabundfile == "not found") { relabundfile = ""; }
            else if (relabundfile == "not open") { relabundfile = ""; abort = true; }
            else { current->setRelAbundFile(relabundfile); }
            
            constaxfile = validParameter.validFile(parameters, "constaxonomy"); //required
            if (constaxfile == "not found") {
                constaxfile = current->getConsTaxonomyFile();
                if (constaxfile != "") { m->mothurOut("Using " + constaxfile + " as input file for the constaxonomy parameter.\n");  }
                else { 	m->mothurOut("[ERROR]: You have no current constaxonomy file and the constaxonomy parameter is required.\n"); abort = true; }
            }
            else if (constaxfile == "not open") { constaxfile = ""; abort = true; }
            else { current->setConsTaxonomyFile(constaxfile); }
            
            if ((relabundfile == "") && (listfile == "") && (sharedfile == "")) { //no files to merge provided, look for currents
                //is there are current file available for any of these?
                //give priority to shared, then list, then relabund
                //if there is a current shared file, use it
                sharedfile = current->getSharedFile();
                if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n");  }
                else {
                    listfile = current->getListFile();
                    if (listfile != "") { m->mothurOut("Using " + listfile + " as input file for the list parameter.\n");  }
                    else {
                        relabundfile = current->getRelAbundFile();
                        if (relabundfile != "") { m->mothurOut("Using " + relabundfile + " as input file for the rabund parameter.\n");  }
                        else {
                            m->mothurOut("[ERROR]: No valid current files. You must provide a list, relabund or shared file.\n");  abort = true;
                        }
                    }
                }
            }
            
            //if the user changes the output directory command factory will send this info to us in the output parameter
            outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){
                outputDir = "";
                outputDir += util.hasPath(constaxfile); //if user entered a file with a path then preserve it
            }
            
            //check for optional parameter and set defaults
            // ...at some point should added some additional type checking...
            label = validParameter.valid(parameters, "label");
            if (label == "not found") { label = ""; }
            else {
                if(label != "all") {  util.splitAtDash(label, labels);  allLines = 0;  }
                else { allLines = 1;  }
            }
            
            
            
            
            
            
            
            
            
            
            //******************  //get taxlevel parameter *********************************
            
            
            
            
            
            
            
            
            
            
            
            
 
        }
    }
    catch(exception& e) {
        m->errorOut(e, "MergeOTUsCommand", "MergeOTUsCommand");
        exit(1);
    }
}

//**********************************************************************************************************************

MergeOTUsCommand::~MergeOTUsCommand(){}

//**********************************************************************************************************************

int MergeOTUsCommand::execute(){
    try {
        if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); }  return 0; }
        
        //output files created by command
        m->mothurOut("\nOutput File Names: \n");
        for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
        
        string currentName = "";
        itTypes = outputTypes.find("list");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setListFile(currentName); }
        }
        
        itTypes = outputTypes.find("shared");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setSharedFile(currentName); }
        }
        
        itTypes = outputTypes.find("relabund");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setRelAbundFile(currentName); }
        }
        
        //set constaxonomy file as new current constaxonomyfile
        itTypes = outputTypes.find("constaxonomy");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setConsTaxonomyFile(currentName); }
        }
        
        m->mothurOut("\nOutput File Names:\n");
        for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	} m->mothurOutEndLine();

        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeOTUsCommand", "execute");
        exit(1);
    }
}

//**********************************************************************************************************************

