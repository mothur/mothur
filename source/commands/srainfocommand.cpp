//
//  srainfocommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 10/29/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "srainfocommand.hpp"

//**********************************************************************************************************************
vector<string> SRAInfoCommand::setParameters(){
    try {
        CommandParameter psra("sra", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(psra);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        vector<string> myArray;
        for (int i = 0; i < parameters.size(); i++) {    myArray.push_back(parameters[i].name);        }
        return myArray;
    }
    catch(exception& e) {
        m->errorOut(e, "SRAInfoCommand", "setParameters");
        exit(1);
    }
}
//**********************************************************************************************************************
string SRAInfoCommand::getHelpString(){
    try {
        string helpString = "";
        helpString += "The sra.info command reads a sra file and extracts the fastq or sff data.\n";
        helpString += "The sra.info  command parameters are ....\n";
        helpString += "The sra.info command should be in the following format: sra.info(sra=yourSRAFile)\n";
        helpString += "sra.info(sra=SRR000004) \n";
        return helpString;
    }
    catch(exception& e) {
        m->errorOut(e, "SRAInfoCommand", "getHelpString");
        exit(1);
    }
}
//**********************************************************************************************************************
string SRAInfoCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fastq") {  pattern = "[filename],fastq"; }
        else if (type == "sff") {  pattern = "[filename],sff"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SRAInfoCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
SRAInfoCommand::SRAInfoCommand(){
    try {
        abort = true; calledHelp = true;
        setParameters();
        vector<string> tempOutNames;
        outputTypes["summary"] = tempOutNames;
    }
    catch(exception& e) {
        m->errorOut(e, "SRAInfoCommand", "SRAInfoCommand");
        exit(1);
    }
}
//***************************************************************************************************************

SRAInfoCommand::SRAInfoCommand(string option)  {
    try {
        abort = false; calledHelp = false;
        
        //allow user to run help
        if(option == "help") { help(); abort = true; calledHelp = true; }
        else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        
        else {
            vector<string> myArray = setParameters();
            
            OptionParser parser(option);
            map<string,string> parameters = parser.getParameters();
            
            ValidParameters validParameter;
            map<string,string>::iterator it;
            
            //check to make sure all parameters are valid for command
            for (it = parameters.begin(); it != parameters.end(); it++) {
                if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
            }
            
            //if the user changes the input directory command factory will send this info to us in the output parameter
            string inputDir = validParameter.valid(parameters, "inputdir");
            if (inputDir == "not found"){    inputDir = "";        }
            else {
                string path;
                it = parameters.find("sra");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {    parameters["sra"] = inputDir + it->second;        }
                }
            }
            
            //initialize outputTypes
            vector<string> tempOutNames;
            outputTypes["summary"] = tempOutNames;
            
            //check for required parameters
            srafile = validParameter.validFile(parameters, "sra");
            if (srafile == "not open") { srafile = ""; abort = true; }
            else if (srafile == "not found") { m->mothurOut("[ERROR]: The sra parameter is required.\n");  abort = true; } 
            
            outputDir = validParameter.valid(parameters, "outputdir");        if (outputDir == "not found"){
                outputDir = "";
                outputDir += util.hasPath(srafile); //if user entered a file with a path then preserve it
            }
        }
    }
    catch(exception& e) {
        m->errorOut(e, "SRAInfoCommand", "SRAInfoCommand");
        exit(1);
    }
}
//***************************************************************************************************************
int SRAInfoCommand::execute(){
    try{
        
        if (abort) { if (calledHelp) { return 0; }  return 2;    }
        
        
    
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "SRAInfoCommand", "execute");
        exit(1);
    }
}
/**************************************************************************************************/
