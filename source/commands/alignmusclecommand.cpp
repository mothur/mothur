//
//  alignmusclecommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 2/16/22.
//  Copyright © 2022 Schloss Lab. All rights reserved.
//

#include "alignmusclecommand.hpp"

//**********************************************************************************************************************
vector<string> AlignMuscleCommand::setParameters(){
    try {
        CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta",false,true,true); parameters.push_back(pfasta);
        CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
        CommandParameter pmusclelocation("muscle", "String", "", "", "", "", "","",false,false); parameters.push_back(pmusclelocation);
        CommandParameter pmethod("method", "Multiple", "align-super5", "align", "", "", "","",false,false,true); parameters.push_back(pmethod);
        CommandParameter pperturb("perturb", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pperturb);
        CommandParameter pperm("perm", "Multiple", "none-abc-acb-bca", "none", "", "", "","",false,false,true); parameters.push_back(pperm);
        CommandParameter pstratified("stratified", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pstratified);
        CommandParameter pdiversified("diversified", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pdiversified);
        CommandParameter preplicates("replicates", "Number", "", "4", "", "", "","",false,false); parameters.push_back(preplicates);
        CommandParameter pconsiters("consiters", "Number", "", "2", "", "", "","",false,false); parameters.push_back(pconsiters);
        CommandParameter piters("refineiters", "Number", "", "100", "", "", "","",false,false); parameters.push_back(piters);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
            
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["fasta"] = tempOutNames;

        vector<string> myArray;
        for (int i = 0; i < parameters.size(); i++) {    myArray.push_back(parameters[i].name);        }
        return myArray;
    }
    catch(exception& e) {
        m->errorOut(e, "AlignMuscleCommand", "setParameters");
        exit(1);
    }
}
//**********************************************************************************************************************
string AlignMuscleCommand::getHelpString(){
    try {
        string helpString = "";
        helpString += "The align.muscle command creates multiple alignments of protein sequences.\n";
        helpString += "This command is a wrapper for muscle written by Robert C. Edgar.\n";
        helpString += "The align.muscle command parameters are fasta, method, perturb, perm, stratified, diversified, replicates, consiters, refineiters, processors and muscle.\n";
        helpString += "The fasta parameter allows you to enter the fasta file containing your sequences, and is required, unless you have a valid current fasta file. \n";
        helpString += "The perturb parameter allows you to provide a random number seed for generating HMM perturbations. Default=0. https://drive5.com/muscle5/manual/hmm_perturbations.html\n";
        helpString += "The perm parameter specifies the guide tree permutation. PERM can be none, abc, acb and bca, default=none. https://drive5.com/muscle5/manual/guide_tree_permutations.html.\n";
        helpString += "The stratified parameter allows you to indicate you would like to generate a stratified ensemble. https://drive5.com/muscle5/manual/stratified_ensemble.html\n";
        helpString += "The diversified parameter allows you to indicate you would like to generate a diversified ensemble. https://drive5.com/muscle5/manual/diversified_ensemble.html\n";
        helpString += "The replicates parameter aloows you to indicate the number of replicates, default 4 for -stratified and 100 for -diversified. With -stratified, one replicate is generated for each guide tree permutation, so the total number of replicates is 4×N.\n";
        helpString += "The consiters parameter allows you to indicate the number of consistency iterations. Default 2.\n";
        helpString += "The refineiters parameter allows you to indicate the number of refinement iterations. Default 100.\n";
        helpString += "The muscle parameter allows you to specify the name and location of your muscle executable. By default mothur will look in your path and mothur's executable and mothur tools locations.  You can set the muscle location as follows, muscle=/usr/bin/muscle5.1.\n";
        helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is all available. \n";
        helpString += "The align.muscle command should be in the following format: \n";
        helpString += "align.muscle(fasta=yourFastaFile) \n";
        helpString += "Example: align.muscle(fasta=prot.fasta) \n";
            
        return helpString;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraUchimeCommand", "getHelpString");
        exit(1);
    }
}
//**********************************************************************************************************************
string AlignMuscleCommand::getCommonQuestions(){
    try {
        vector<string> questions, issues, qanswers, ianswers, howtos, hanswers;
        
        string commonQuestions = util.getFormattedHelp(questions, qanswers, issues, ianswers, howtos, hanswers);
        
        return commonQuestions;
    }
    catch(exception& e) {
        m->errorOut(e, "AlignMuscleCommand", "getCommonQuestions");
        exit(1);
    }
}
//**********************************************************************************************************************
string AlignMuscleCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],[tag].fasta"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "AlignMuscleCommand", "getOutputPattern");
        exit(1);
    }
}
//***************************************************************************************************************
AlignMuscleCommand::AlignMuscleCommand(string option) : Command()  {
    try {
        //allow user to run help
        if(option == "help") { help(); abort = true; calledHelp = true; }
        else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
        
        else {
            OptionParser parser(option, setParameters());
            map<string,string> parameters = parser.getParameters();
            
            ValidParameters validParameter;
            fastafile = validParameter.validFile(parameters, "fasta");
            if (fastafile == "not found") {
                fastafile = current->getFastaFile();
                if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter.\n"); }
                else {     m->mothurOut("[ERROR]: You have no current fasta file and the fasta parameter is required.\n");  abort = true; }
            }
            else if (fastafile == "not open") { abort = true; }
            else { current->setFastaFile(fastafile); }
            
            string temp = validParameter.valid(parameters, "processors");    if (temp == "not found"){    temp = current->getProcessors();    }
            processors = current->setProcessors(temp);
            
            perturb = validParameter.valid(parameters, "perturb");    if (perturb == "not found"){  usePerturb = false;  perturb = "0";    }else{  usePerturb = true;  }
            
            method = validParameter.valid(parameters, "method");
            if (method == "not found") {  method = "align";}
            
            if ((method == "align") || (method == "super5")) {}
            else { m->mothurOut("[WARNING]: " + method + " is not a valid method. Options are align or super5, using align.\n"); method = "align"; }
              
            perm = validParameter.valid(parameters, "perm");
            if (perm == "not found") {  perm = "none";  usePerm = false; } else { usePerm = true; }
            
            if ((perm == "none") || (perm == "abc") || (perm == "acb") || (perm == "bca")) {}
            else { m->mothurOut("[WARNING]: " + perm + " is not a valid perm option. Options are none, abc, acb or bca, using none.\n"); perm = "none"; }
            
            temp = validParameter.valid(parameters, "stratified");            if (temp == "not found") { temp = "f"; }
            stratified = util.isTrue(temp);
            if (stratified) { replicates = "4"; }
           
            temp = validParameter.valid(parameters, "diversified");            if (temp == "not found") { temp = "f"; }
            diversified = util.isTrue(temp);
            if (diversified) { replicates = "100"; }
            
            replicates = validParameter.valid(parameters, "replicates");  if (replicates == "not found")            { useReplicates = false;   }    else{ useReplicates = true;         }
         
            consiters = validParameter.valid(parameters, "consiters");   if (consiters == "not found")            { useConsiters = false; consiters = "2";                }    else{ useConsiters = true;            }
            
            refineiters = validParameter.valid(parameters, "refineiters");   if (refineiters == "not found")            { useRefineiters = false; refineiters = "100";                }    else{ useRefineiters = true;            }

                    
            vector<string> versionOutputs;
            bool foundTool = false;
            string programName = "muscle"; programName += EXECUTABLE_EXT;
            
            muscleLocation = validParameter.validFile(parameters, "muscle");
            if (muscleLocation == "not found") {
                muscleLocation = "";
                foundTool = util.findTool(programName, muscleLocation, versionOutputs, current->getLocations());
            }
            else {
                //test to make sure muscle exists
                muscleLocation = util.getFullPathName(muscleLocation);
                ifstream in; foundTool = util.openInputFile(muscleLocation, in, "no error"); in.close();
                if(!foundTool) {
                    m->mothurOut(muscleLocation + " file does not exist or cannot be opened, ignoring.\n"); muscleLocation = "";
                    foundTool = util.findTool(programName, muscleLocation, versionOutputs, current->getLocations());
                }
            }
            
            if (!foundTool) { abort = true; }
            
            muscleLocation = util.getFullPathName(muscleLocation);
            if (m->getDebug()) { m->mothurOut("[DEBUG]: muscle location using " + muscleLocation + "\n"); }
        }
    }
    catch(exception& e) {
        m->errorOut(e, "AlignMuscleCommand", "AlignMuscleCommand");
        exit(1);
    }
}

//***************************************************************************************************************

int AlignMuscleCommand::execute(){
    try{
        
        if (abort) { if (calledHelp) { return 0; }  return 2;    }

        m->mothurOut("Processing sequences from " + fastafile + " ...\n" );
        
        long start = time(NULL);
        
        
        
        if (outputdir == "") { outputdir = util.hasPath(fastafile);  }
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(fastafile));
        //variables["[tag]"] = "denovo";
        //string outputFileName = getOutputFileName("fasta", variables);
                
        
        
        
        
        //set accnos file as new current accnosfile
        string currentName = "";
        itTypes = outputTypes.find("fasta");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
        }
        
        m->mothurOut("\nOutput File Names: \n");
        for (int i = 0; i < outputNames.size(); i++) {    m->mothurOut(outputNames[i]); m->mothurOutEndLine();    }
        m->mothurOutEndLine();
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "AlignMuscleCommand", "execute");
        exit(1);
    }
}
//**********************************************************************************************************************
void AlignMuscleCommand::driver(){
    try {
        string outputFileName = "";
        //to allow for spaces in the path
       /* params->driverOutputFName = "\"" + params->driverOutputFName + "\"";
        params->formattedFastaFilename = "\"" + params->formattedFastaFilename + "\"";
        params->driverAlns = "\"" + params->driverAlns + "\"";
        
        if (params->formattedFastaFilename.length() > 257) {
            params->m->mothurOut("[ERROR]: " + params->formattedFastaFilename + " filename is " + toString(params->formattedFastaFilename.length()) + " long. The uchime program can't handle files with a full path longer than 257 characters, please correct.\n"); params->m->setControl_pressed(true); return 0;
        }else if ((params->driverAlns.length() > 257) && (params->vars->chimealns)) {
            params->m->mothurOut("[ERROR]: " + params->driverAlns + " filename is " + toString(params->driverAlns.length()) + " long. The uchime program can't handle files with a full path longer than 257 characters, please correct.\n"); params->m->setControl_pressed(true); return 0;
        }else if (params->driverOutputFName.length() > 257) {
            params->m->mothurOut("[ERROR]: " + params->driverOutputFName + " filename is " + toString(params->driverOutputFName.length()) + " long. The uchime program can't handle files with a full path longer than 257 characters, please correct input file name.\n"); params->m->setControl_pressed(true); return 0;
        }
        */
        vector<char*> cPara;
        
        string muscleCommand = muscleLocation;
        muscleCommand = "\"" + muscleCommand + "\" ";
        cPara.push_back(util.mothurConvert(muscleCommand));
        
        cPara.push_back(util.mothurConvert("-"+method));
        cPara.push_back(util.mothurConvert(fastafile));
        
        //output filename
        cPara.push_back(util.mothurConvert("-output"));
        cPara.push_back(util.mothurConvert(outputFileName));
        
        if (usePerturb) {
            cPara.push_back(util.mothurConvert("-perturb"));
            cPara.push_back(util.mothurConvert(perturb));
        }
        
        if (usePerm) {
            cPara.push_back(util.mothurConvert("-perm"));
            cPara.push_back(util.mothurConvert(perm));
        }
        
        if (stratified) { cPara.push_back(util.mothurConvert("-stratified")); }
        if (diversified) { cPara.push_back(util.mothurConvert("-diversified")); }

        if (useReplicates) {
            cPara.push_back(util.mothurConvert("-replicates"));
            cPara.push_back(util.mothurConvert(replicates));
        }
        
        if (useConsiters) {
            cPara.push_back(util.mothurConvert("-consiters"));
            cPara.push_back(util.mothurConvert(consiters));
        }
        
        if (useConsiters) {
            cPara.push_back(util.mothurConvert("-consiters"));
            cPara.push_back(util.mothurConvert(consiters));
        }
        
        if (useRefineiters) {
            cPara.push_back(util.mothurConvert("-refineiters"));
            cPara.push_back(util.mothurConvert(refineiters));
        }
        
        char** muscleParameters;
        muscleParameters = new char*[cPara.size()];
        string commandString = "";
        for (int i = 0; i < cPara.size(); i++) {  muscleParameters[i] = cPara[i];  commandString += toString(cPara[i]) + " "; }
        
#if defined NON_WINDOWS
#else
        commandString = "\"" + commandString + "\"";
#endif
        if (m->getDebug()) { m->mothurOut("[DEBUG]: muscle command = " + commandString + ".\n"); }
        system(commandString.c_str());
        
        //free memory
        for(int i = 0; i < cPara.size(); i++)  {  delete cPara[i];  }
        delete[] muscleParameters;
        
        //remove "" from filenames
        //params->driverOutputFName = params->driverOutputFName.substr(1, params->driverOutputFName.length()-2);
        //params->formattedFastaFilename = params->formattedFastaFilename.substr(1, params->formattedFastaFilename.length()-2);
        //params->driverAlns = params->driverAlns.substr(1, params->driverAlns.length()-2);
        
    }
    catch(exception& e) {
        m->errorOut(e, "AlignMuscleCommand", "driver");
        exit(1);
    }
}
//***************************************************************************************************************

