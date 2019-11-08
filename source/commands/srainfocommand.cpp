//
//  srainfocommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 10/29/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "srainfocommand.hpp"
#include "systemcommand.h"

//**********************************************************************************************************************
vector<string> SRAInfoCommand::setParameters(){
    try {
        CommandParameter psra("sra", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(psra);
        CommandParameter ptype("type", "Multiple", "fastq-sff", "fastq", "", "", "","",false,false,true); parameters.push_back(ptype);
        CommandParameter pFasterQlocation("fasterq", "String", "", "", "", "", "","",false,false); parameters.push_back(pFasterQlocation);
        CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
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
        helpString += "The sra.info command parameters are ....\n";
        helpString += "The processors parameter allows you to specify how many processors you would like to use. The default is all available. \n";
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
            
            outputType = validParameter.valid(parameters, "type");        if (outputType == "not found"){
                outputType = "fastq";
            }
            if ((outputType == "fastq") || (outputType == "sff")) {}
            else { m->mothurOut("[ERROR]: " + outputType + " is not a valid type option. The sra.info can extract files of type sff or fastq.\n");  abort = true; }
            
            string temp = validParameter.valid(parameters, "processors");    if (temp == "not found"){    temp = current->getProcessors();    }
            processors = current->setProcessors(temp);
            
            vector<string> versionOutputs;
            bool foundTool = false;
            string path = current->getProgramPath();
            string programName = "fasterq_dump"; programName += EXECUTABLE_EXT;
            
            fasterQLocation = validParameter.valid(parameters, "fasterq");
            if (fasterQLocation == "not found") {
                fasterQLocation = "";
                foundTool = util.findTool(programName, fasterQLocation, path, versionOutputs, current->getLocations());
            }else {
                //test to make sure vsearch exists
                ifstream in;
                fasterQLocation = util.getFullPathName(fasterQLocation);
                bool ableToOpen = util.openInputFile(fasterQLocation, in, "no error"); in.close();
                if(!ableToOpen) {
                    m->mothurOut(fasterQLocation + " file does not exist or cannot be opened, ignoring.\n"); fasterQLocation = "";
                    programName = util.getSimpleName(fasterQLocation); fasterQLocation = "";
                    foundTool = util.findTool(programName, fasterQLocation, path, versionOutputs, current->getLocations());
                }
            }
          
            if (foundTool && !abort) { //check fasterq_dump version
                if (versionOutputs.size() >= 3) {
                    string version = versionOutputs[2];
                                                
                    if (version != "2.9.6") {
                        m->mothurOut("[ERROR]: " + programName + " version found = " + version + ". Mothur requires version 2.9.6 which is distributed with mothur's executable or available for download here, https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software\n");  abort = true;
                    }else { m->mothurOut("Using " + programName + " version " + version + ".\n"); }
                }
            }
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: fasterq_dump location using " + fasterQLocation + "\n"); }
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
        
        if (outputType == "fastq")      {  runFastqDump();  }
        //else if (outputType == "sff")   {  runSFFDump();    }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "SRAInfoCommand", "execute");
        exit(1);
    }
}
//***************************************************************************************************************
void SRAInfoCommand::runFastqDump(){
    try{
        
        vector<char*> cPara;
        string fasterQCommand = fasterQLocation;
        fasterQCommand = "\"" + fasterQCommand + "\" " + srafile + " ";
        
        char* tempFasterQ;
        tempFasterQ= new char[fasterQCommand.length()+1];
        *tempFasterQ = '\0';
        strncat(tempFasterQ, fasterQCommand.c_str(), fasterQCommand.length());
        cPara.push_back(tempFasterQ);
        
        //--force - overwrite output files if they exist
        char* force = new char[8];  force[0] = '\0'; strncat(force, "--force", 7);
        cPara.push_back(force);
        
        //-S|--split-files                 write reads into different files
        char* splitFiles = new char[3];     splitFiles[0] = '\0'; strncat(splitFiles, "-S", 2);
        cPara.push_back(splitFiles);
        
        //-3|--split-3                     writes single reads in special file
        char* splitSingleFiles = new char[3];     splitSingleFiles[0] = '\0'; strncat(splitSingleFiles, "-3", 2);
        cPara.push_back(splitSingleFiles);
        
        //--threads=processors
        char* threads = new char[10];  threads[0] = '\0'; strncat(threads, "--threads", 9);
        cPara.push_back(threads);
        string numProcessors = toString(processors);
        char* tempThreads = new char[numProcessors.length()+1];
        *tempThreads = '\0'; strncat(tempThreads, numProcessors.c_str(), numProcessors.length());
        cPara.push_back(tempThreads);
        
        if (outputDir != "") {
            char* outDir = new char[9];  outDir[0] = '\0'; strncat(outDir, "--outdir", 8);
            cPara.push_back(outDir);
            char* tempoutputDir = new char[outputDir.length()+1];
            *tempoutputDir = '\0'; strncat(tempoutputDir, outputDir.c_str(), outputDir.length());
            cPara.push_back(tempoutputDir);
        }
       
        char** fasterQParameters;
        fasterQParameters = new char*[cPara.size()];
        string commandString = "";
        for (int i = 0; i < cPara.size(); i++) {  fasterQParameters[i] = cPara[i];  commandString += toString(cPara[i]) + " "; }
        
#if defined NON_WINDOWS
#else
        commandString = "\"" + commandString + "\"";
#endif
        
        if (m->getDebug()) { m->mothurOut("[DEBUG]: fasterq_dump command = " + commandString + ".\n"); }
        m->mothurOut("[DEBUG]: fasterq_dump command = " + commandString + ".\n");
        
        system(commandString.c_str());
        
        //free memory
        for(int i = 0; i < cPara.size(); i++)  {  delete cPara[i];  }
        delete[] fasterQParameters;
        
        return;
    }
    catch(exception& e) {
        m->errorOut(e, "SRAInfoCommand", "runFastqDump");
        exit(1);
    }
}
/**************************************************************************************************/
