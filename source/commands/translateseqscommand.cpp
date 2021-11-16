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
        //CommandParameter pdna("dna", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pdna);
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
        helpString += "The translate.seqs command reads a fastafile containing dna and translates it to amino acids, .....\n";
        helpString += "The translate.seqs command parameters are fasta, frames, stop, dna and processors.\n";
        helpString += "The frames parameter is used to indicate the reading frames you want to use. Options are 1, 2, 3, -1, -2, -3. You can select multiple frames by separating them with '|' characters. For example: frames=1|-1|2.\n";
        helpString += "The stop parameter is used to indicate when to stop the translation. If T, then if the translation hits a stop codon, it stops before that codon. If F, it returns the full translation with a * as the stop codon. Default=t.\n";
        //helpString += "The dna parameter is used to indicate which type of data is in your fasta file. Default=t, meaning the data is dna. dna=f, means the data in the fasta file is amino acids. \n";
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
        
        if (type == "fasta") {  pattern = "[filename],[tag2],fasta-[filename],[tag],align"; }
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
                if ((thisFrame == 1) || (thisFrame == 2) || (thisFrame == 3) || (thisFrame == -2) || (thisFrame == -1) || (thisFrame == -3)) {
                    frames.push_back(thisFrame); }
                else {
                    m->mothurOut("[WARNING]: " + fs[i] + " is not a valid frame option. Options include 1,2,3,-1,-2,-3. Ignoring " + fs[i] + "\n");
                }
            }
            if (frames.size() == 0) { abort = true; }
            
            temp = validParameter.valid(parameters, "stop");        if (temp == "not found") { temp = "T"; }
            stop = util.isTrue(temp);
            
            //temp = validParameter.valid(parameters, "dna");        if (temp == "not found") { temp = "T"; }
            //dna = util.isTrue(temp);
            
            //force single frames for aa to dna translation
            //if (dna) { frames.clear(); frames.push_back(1); }
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
        
        string thisOutputDir = outputdir;
        if (outputdir == "") {  thisOutputDir += util.hasPath(fastafile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastafile));
        
        vector<linePair> lines; vector<double> positions;
#if defined NON_WINDOWS
        positions = util.divideFile(fastafile, processors);
        for (int i = 0; i < (positions.size()-1); i++) {    lines.push_back(linePair(positions[i], positions[(i+1)]));    }
#else
        long long numFastaSeqs = 0;
        positions = util.setFilePosFasta(filename, numFastaSeqs);
        if (numFastaSeqs < processors) { processors = numFastaSeqs; m->mothurOut("Reducing processors to " + toString(numFastaSeqs) + ".\n");  }
        
        //figure out how many sequences you have to process
        int numSeqsPerProcessor = numFastaSeqs / processors;
        for (int i = 0; i < processors; i++) {
            int startIndex =  i * numSeqsPerProcessor;
            if(i == (processors - 1)){    numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor;     }
            lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
        }
#endif
        
        double numSeqs = 0;
        for (int i = 0; i < frames.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            m->mothurOut("\nTranslating sequences to amino acids using frame " + toString(frames[i]) + ":\n");
            
            variables["[tag2]"] = "aa"+toString(frames[i]);
            string outputFileName = getOutputFileName("fasta", variables);
            outputTypes["fasta"].push_back(outputFileName);  outputNames.push_back(outputFileName);
            
            numSeqs = createProcessesTranslateDNAtoAminoAcids(outputFileName, lines, frames[i]);
        }
        
        m->mothurOut("\nIt took " + toString(time(NULL) - start) + " seconds to translate " + toString(numSeqs) + " sequences.\n");

        //set accnos file as new current accnosfile
        string currentName = "";
        itTypes = outputTypes.find("fasta");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
        }
        
        m->mothurOut("\nOutput File Names:\n");
        for (int i = 0; i < outputNames.size(); i++) {    m->mothurOut(outputNames[i]); m->mothurOutEndLine();    }
        m->mothurOutEndLine();
 
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "TranslateSeqsCommand", "execute");
        exit(1);
    }
}
/**********************************************************************************************************************
void translateToDNADriver(translateSeqsStruct* params) {
    try {
        ifstream inFASTA; params->util.openInputFile(params->inputFilename, inFASTA);
        inFASTA.seekg(params->filePos.start);

        bool done = false;
        long long count = 0;
        while (!done) {
            
            if (params->m->getControl_pressed()) {  break; }
            
            Protein seq(inFASTA); params->util.gobble(inFASTA);
            
            if (seq.getName() != "") {
                string translation; getDNA(seq.getAligned(), translation, params->m);
                params->outputWriter->write('>' + seq.getName() + '\n' + translation + '\n');
            }
            
            #if defined NON_WINDOWS
                unsigned long long pos = inFASTA.tellg();
                if ((pos == -1) || (pos >= params->filePos.end)) { break; }
            #else
                if (count == params->filePos.end) { break; }
            #endif

            //report progress
            if((count) % 1000 == 0){    params->m->mothurOutJustToScreen(toString(count) + "\n");         }

        }
        //report progress
        if((count) % 1000 != 0){    params->m->mothurOutJustToScreen(toString(count) + "\n");         }
        params->numSeqs = count; inFASTA.close();
    }
    catch(exception& e) {
        params->m->errorOut(e, "TranslateSeqsCommand", "translateToDNADriver");
        exit(1);
    }
}*/
//**********************************************************************************************************************
void translateToAminoAcidDriver(translateSeqsStruct* params) {
    try {
        ifstream inFASTA; params->util.openInputFile(params->inputFilename, inFASTA);
        inFASTA.seekg(params->filePos.start);

        bool done = false;
        long long count = 0;
        while (!done) {
            
            if (params->m->getControl_pressed()) {  break; }
            
            Sequence seq(inFASTA); params->util.gobble(inFASTA);
            
            if (seq.getName() != "") {
                Protein prot = seq.getProtein(params->frame);
                prot.printProtein(params->outputWriter);
                count++;
            }
            
            #if defined NON_WINDOWS
                unsigned long long pos = inFASTA.tellg();
                if ((pos == -1) || (pos >= params->filePos.end)) { break; }
            #else
                if (count == params->filePos.end) { break; }
            #endif

            //report progress
            if((count) % 1000 == 0){    params->m->mothurOutJustToScreen(toString(count) + "\n");         }

        }
        //report progress
        if((count) % 1000 != 0){    params->m->mothurOutJustToScreen(toString(count) + "\n");         }
        params->numSeqs = count; inFASTA.close();
    }
    catch(exception& e) {
        params->m->errorOut(e, "TranslateSeqsCommand", "translateToAminoAcidDriver");
        exit(1);
    }
}
//***************************************************************************************************************
//    translateSeqsStruct (linePair fP, OutputWriter* oFName, string fname, bool st, bool dn, int frame) {
double TranslateSeqsCommand::createProcessesTranslateDNAtoAminoAcids(string outputFileName, vector<linePair> lines, int frame) {
    try {
        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<translateSeqsStruct*> data;
        
        auto synchronizedOutputFile = std::make_shared<SynchronizedOutputFile>(outputFileName);
        
        for (int i = 0; i < processors-1; i++) {
            
            OutputWriter* threadOutputWriter = new OutputWriter(synchronizedOutputFile);
            
            translateSeqsStruct* dataBundle = new translateSeqsStruct(lines[i+1], threadOutputWriter, fastafile, stop, dna, frame);
            data.push_back(dataBundle);

            workerThreads.push_back(new std::thread(translateToAminoAcidDriver, dataBundle));
         }
        
        OutputWriter* threadOutputWriter = new OutputWriter(synchronizedOutputFile);
        translateSeqsStruct* dataBundle = new translateSeqsStruct(lines[0], threadOutputWriter, fastafile, stop, dna, frame);
        
        translateToAminoAcidDriver(dataBundle);
        double num = dataBundle->numSeqs;
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            num += data[i]->numSeqs;
            
            delete data[i]->outputWriter;
            delete data[i];
            delete workerThreads[i];
        }
        synchronizedOutputFile->close(); delete threadOutputWriter;  delete dataBundle;
        
        return num;
    }
    catch(exception& e) {
        m->errorOut(e, "TranslateSeqsCommand", "createProcesses");
        exit(1);
    }
}
//***************************************************************************************************************
