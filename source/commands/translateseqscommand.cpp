//
//  translateseqscommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 11/8/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#include "translateseqscommand.hpp"

//**********************************************************************************************************************
vector<string> TranslateSeqsCommand::setParameters(){
    try {
        CommandParameter pfasta("fasta", "InputTypes", "", "", "FastaReport", "none", "none","summary",false,true,true); parameters.push_back(pfasta);
        CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
        CommandParameter pframes("frames", "Multiple", "1-2-3--1--2--3", "1", "", "", "","",true,false,true); parameters.push_back(pframes);
        CommandParameter pstop("stop", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pstop);
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
        m->errorOut(e, "TranslateSeqsCommand", "setParameters");
        exit(1);
    }
}
//**********************************************************************************************************************
string TranslateSeqsCommand::getHelpString(){
    try {
        string helpString = "";
        helpString += "The translate.seqs command reads a fastafile containing dna and translates it to amino acids, reads a fasta file containing amino acids and translates it to dna.\n";
        helpString += "The translate.seqs command parameters are fasta, frames, stop and processors.\n";
        helpString += "The frames parameter is used to indicate the reading frames you want to use. Options are 1, 2, 3, -1, -2, -3. You can select multiple frames by separating them with '|' characters. For example: frames=1|-1|2.\n";
        helpString += "The frames parameter is used to indicate when to stop the translation. If T, then if the translation hits a stop codon, it stops before that codon. If F, it returns the full translation with a * as the stop codon. Default=t.\n";
        helpString += "The translate.seqs command should be in the following format: \n";
        helpString += "translate.seqs(fasta=yourFastaFile, processors=2) \n";
            
        return helpString;
    }
    catch(exception& e) {
        m->errorOut(e, "TranslateSeqsCommand", "getHelpString");
        exit(1);
    }
}
//**********************************************************************************************************************
string TranslateSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],[tag],fasta"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "TranslateSeqsCommand", "getOutputPattern");
        exit(1);
    }
}
//***************************************************************************************************************

TranslateSeqsCommand::TranslateSeqsCommand(string option) : Command()  {
    try {
        if(option == "help") { help(); abort = true; calledHelp = true; }
        else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
        
        else {
            OptionParser parser(option, setParameters());
            map<string,string> parameters = parser.getParameters();
            
            ValidParameters validParameter;
            fastafile = validParameter.validFile(parameters, "fasta");
            if (fastafile == "not open") { abort = true; }
            else if (fastafile == "not found") {  fastafile = "";  }
            else { current->setFastaFile(fastafile); }
            
            string temp = validParameter.valid(parameters, "processors");    if (temp == "not found"){    temp = current->getProcessors();    }
            processors = current->setProcessors(temp);
            
            temp = validParameter.valid(parameters, "frames");    if (temp == "not found"){    temp = "1";    }
            vector<string> fs; util.splitAtChar(temp, fs, '|');
            for (int i = 0; i < fs.size(); i++) {
                int thisFrame; util.mothurConvert(fs[i], thisFrame);
                frames.push_back(thisFrame);
            }
            
            temp = validParameter.valid(parameters, "stop");        if (temp == "not found") { temp = "T"; }
            stop = util.isTrue(temp);
            
        }
    }
    catch(exception& e) {
        m->errorOut(e, "TranslateSeqsCommand", "TranslateSeqsCommand");
        exit(1);
    }
}
//***************************************************************************************************************

int TranslateSeqsCommand::execute(){
    try{
        
        if (abort) { if (calledHelp) { return 0; }  return 2;    }
        
        long start = time(NULL);
        
        for (int i = 0; i < frames.size(); i++) {
            cout << frames[i] << endl;
        }
        
        string thisOutputDir = outputdir;
        if (outputdir == "") {  thisOutputDir += util.hasPath(fastafile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastafile));
        variables["[tag]"] = "...... TODO ........ ";
        string outputFileName = getOutputFileName("fasta", variables);
        outputTypes["fasta"].push_back(outputFileName);  outputNames.push_back(outputFileName);
    
        
        
        
        
        //if dna,
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "TranslateSeqsCommand", "execute");
        exit(1);
    }
}

//***************************************************************************************************************
