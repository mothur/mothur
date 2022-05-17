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
        CommandParameter pamino("amino", "InputTypes", "", "", "FastaReport", "none", "none","summary",false,true,true); parameters.push_back(pamino);
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
        helpString += "The translate.seqs command reads a fastafile containing dna and translates it to amino acids. Alternatively, you can provide an amino file with aligned amino acids and fasta file with unaligned dna sequences and mothur will align the dna to the amino acid. Mothur expects the aligned amino acids and the fasta file to contain the same sequence names.\n";
        helpString += "The translate.seqs command parameters are fasta, amino, frames, stop and processors.\n";
        helpString += "The fasta parameter is used to provide a file containing DNA sequences. It is required.\n";
        helpString += "The amino parameter is used to provide a file related to the fasta file containing amino acid sequences. The amino file should be aligned and mothur will align the DNA reads to the amino acids. Mothur assumes both files are in the same frame.\n";
        helpString += "The frames parameter is used to indicate the reading frames you want to use. Options are 1, 2, 3, -1, -2, -3. You can select multiple frames by separating them with '|' characters. For example: frames=1|-1|2.\n";
        helpString += "The stop parameter is used to indicate when to stop the translation. If T, then if the translation hits a stop codon, it stops before that codon. If F, it returns the full translation with a * as the stop codon. Default=t.\n";
        helpString += "The translate.seqs command should be in the following format: \n";
        helpString += "translate.seqs(fasta=yourFastaFile, processors=2) or translate.seqs(amino=yourProteinSeqs, fasta=yourDNASeqs)\n";
            
        getCommonQuestions();
        
        return helpString;
    }
    catch(exception& e) {
        m->errorOut(e, "TranslateSeqsCommand", "getHelpString");
        exit(1);
    }
}
//**********************************************************************************************************************
string TranslateSeqsCommand::getCommonQuestions(){
    try {
        vector<string> questions, issues, qanswers, ianswers, howtos, hanswers;
        
        string commonQuestions = util.getFormattedHelp(questions, qanswers, issues, ianswers, howtos, hanswers);

        return commonQuestions;
    }
    catch(exception& e) {
        m->errorOut(e, "TranslateSeqsCommand", "getCommonQuestions");
        exit(1);
    }
}
//**********************************************************************************************************************
string TranslateSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],align-[filename],[tag],fasta"; }
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
            if (fastafile == "not open") { fastafile = ""; abort = true; }
            else if (fastafile == "not found") {  fastafile = ""; abort = true; }
            else { current->setFastaFile(fastafile); }
            
            aminofile = validParameter.validFile(parameters, "amino");
            if (aminofile == "not open") { aminofile = ""; abort = true; }
            else if (aminofile == "not found") {  aminofile = "";  }
            
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
        
        if (aminofile != "")    { alignDNAAmino();          }  //if amino file is provided then we are aligning dna to aligned amino acids, assumes same frame
        else                    { translateDNAtoAmino();    }
        
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
//**********************************************************************************************************************
void align(Alignment* alignment, Sequence& seq, Protein& prot, MothurOut* m) {
    try {
        
        int alignmentSize = max(prot.getAligned().size(), (seq.getUnaligned().size() / 3)+1);
        
        if (alignmentSize+1 > alignment->getnRows()) {
            if (m->getDebug()) { m->mothurOut("[DEBUG]: " + seq.getName() + " " + toString(seq.getUnaligned().length()) + " " + toString(alignment->getnRows()) + " \n"); }
            alignment->resize(alignmentSize+2);
        }
        
        alignment->align(seq, prot);
        
        seq.setAligned(alignment->getSeqAAln());
                
        delete alignment;
    }
    catch(exception& e) {
        m->errorOut(e, "TranslateSeqsCommand", "alignDNA");
        exit(1);
    }
}

//**********************************************************************************************************************
//aligns dna to aligned amino acids
void alignAminoDriver(alignAminoStruct* params) {
    try {
        ifstream inFASTA; params->util.openInputFile(params->fastaFilename, inFASTA);
        inFASTA.seekg(params->fastaPos.start);
        
        ifstream inAMINO; params->util.openInputFile(params->aminoFilename, inAMINO);
        inAMINO.seekg(params->aminoPos.start);

        bool done = false; long long count = 0;
        
        while (!done) {
            
            if (params->m->getControl_pressed()) {  break; }
            
            Sequence seq(inFASTA); gobble(inFASTA);
            Protein prot(inAMINO); gobble(inAMINO);
            
            if ((seq.getName() != "") && (prot.getName() != "") && (seq.getName() == prot.getName()))  {
                
                align(params->alignment, seq, prot, params->m); count++;
                
                params->outputWriter->write('>' + seq.getName() + '\n' + seq.getAligned() + '\n');
            }
            
            #if defined NON_WINDOWS
                unsigned long long pos = inFASTA.tellg();
                if ((pos == -1) || (pos >= params->fastaPos.end)) { break; }
                unsigned long long pos2 = inAMINO.tellg();
                if ((pos2 == -1) || (pos2 >= params->aminoPos.end)) { break; }
            #else
                if (count == params->fastaPos.end) { break; }
                if (count == params->aminoPos.end) { break; }
            #endif

            //report progress
            if((count) % 1000 == 0){    params->m->mothurOutJustToScreen(toString(count) + "\n");         }

        }
        //report progress
        if((count) % 1000 != 0){    params->m->mothurOutJustToScreen(toString(count) + "\n");         }
        params->numSeqs = count; inFASTA.close(); inAMINO.close();
    }
    catch(exception& e) {
        params->m->errorOut(e, "TranslateSeqsCommand", "alignAminoDriver");
        exit(1);
    }
}
//**********************************************************************************************************************
void TranslateSeqsCommand::alignDNAAmino() {
    try {
        long start = time(nullptr);
        
        //fills lines and alines. Also sets dnaAligned and aminoAligned
        setLines();
        
        //create output file names
        string thisOutputDir = outputdir; string outputFileName = "";
        map<string, string> variables;
        if (outputdir == "") {  thisOutputDir += util.hasPath(fastafile);  }
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastafile));
        outputFileName = getOutputFileName("fasta", variables);
        outputTypes["fasta"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        
        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<alignAminoStruct*> data;
        
        auto synchronizedOutputFile = std::make_shared<SynchronizedOutputFile>(outputFileName);
        
        for (int i = 0; i < processors-1; i++) {
            
            OutputWriter* threadOutputWriter = new OutputWriter(synchronizedOutputFile);
            
            alignAminoStruct* dataBundle = new alignAminoStruct(lines[i+1], aLines[i+1], threadOutputWriter, fastafile, aminofile, stop);
            data.push_back(dataBundle);

            workerThreads.push_back(new std::thread(alignAminoDriver, dataBundle));
         }
        
        OutputWriter* threadOutputWriter = new OutputWriter(synchronizedOutputFile);
        alignAminoStruct* dataBundle = new alignAminoStruct(lines[0], aLines[0], threadOutputWriter, fastafile, aminofile, stop);
        
        alignAminoDriver(dataBundle);
        double num = dataBundle->numSeqs;
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            num += data[i]->numSeqs;
            
            delete data[i]->outputWriter;
            delete data[i];
            delete workerThreads[i];
        }
        synchronizedOutputFile->close(); delete threadOutputWriter;  delete dataBundle;
                            
        m->mothurOut("\nIt took " + toString(time(nullptr) - start) + " seconds to align " + toString(num) + " sequences.\n");
        
    }
    catch(exception& e) {
        m->errorOut(e, "TranslateSeqsCommand", "alignDNAAmino");
        exit(1);
    }
}
//**********************************************************************************************************************
void TranslateSeqsCommand::translateDNAtoAmino() {
    try {
        long start = time(nullptr);
        
        string thisOutputDir = outputdir;
        if (outputdir == "") {  thisOutputDir += util.hasPath(fastafile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastafile));
        
        vector<double> positions;
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
            
            variables["[tag]"] = "aa"+toString(frames[i]);
            string outputFileName = getOutputFileName("fasta", variables);
            outputTypes["fasta"].push_back(outputFileName);  outputNames.push_back(outputFileName);
            
            numSeqs = createProcessesTranslateDNAtoAminoAcids(outputFileName, lines, frames[i]);
        }
        
        m->mothurOut("\nIt took " + toString(time(nullptr) - start) + " seconds to translate " + toString(numSeqs) + " sequences.\n");
        
    }
    catch(exception& e) {
        m->errorOut(e, "TranslateSeqsCommand", "translateDNAtoAmino");
        exit(1);
    }
}
//**********************************************************************************************************************
void translateToAminoAcidDriver(translateSeqsStruct* params) {
    try {
        ifstream inFASTA; params->util.openInputFile(params->inputFilename, inFASTA);
        inFASTA.seekg(params->filePos.start);

        bool done = false;
        long long count = 0;
        while (!done) {
            
            if (params->m->getControl_pressed()) {  break; }
            
            Sequence seq(inFASTA); gobble(inFASTA);
            
            if (seq.getName() != "") {
                Protein prot = seq.getProtein(params->frame, params->stop);
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
            
            translateSeqsStruct* dataBundle = new translateSeqsStruct(lines[i+1], threadOutputWriter, fastafile, stop, frame);
            data.push_back(dataBundle);

            workerThreads.push_back(new std::thread(translateToAminoAcidDriver, dataBundle));
         }
        
        OutputWriter* threadOutputWriter = new OutputWriter(synchronizedOutputFile);
        translateSeqsStruct* dataBundle = new translateSeqsStruct(lines[0], threadOutputWriter, fastafile, stop, frame);
        
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

bool TranslateSeqsCommand::setLines() {
    try {
        
        vector<double> fastaFilePos;
        vector<double> afileFilePos;
        
#if defined NON_WINDOWS
        //set file positions for fasta file
        fastaFilePos = util.divideFile(fastafile, processors);
        
        //get name of first sequence in each chunk
        map<string, int> firstSeqNames;
        for (int i = 0; i < (fastaFilePos.size()-1); i++) {
            
            ifstream in; util.openInputFile(fastafile, in);
            in.seekg(fastaFilePos[i]);
            
            //adjust start if null strings
            if (i == 0) {  util.zapGremlins(in); gobble(in);  }
            
            Sequence temp(in);
            firstSeqNames[temp.getName()] = i;
                        
            in.close();
        }
        
        if(aminofile != "")    {
            //seach for filePos of each first names in the aminofile and save in afileFilePos
            ifstream inAmino; util.openInputFile(aminofile, inAmino);
            
            string input;
            while(!inAmino.eof()){
                input = util.getline(inAmino);
                
                if (input.length() != 0) {
                    if(input[0] == '>'){ //this is a sequence name line
                        istringstream nameStream(input);
                        
                        string sname = "";  nameStream >> sname;
                        sname = sname.substr(1);
                        
                        util.checkName(sname);
                        
                        map<string, int>::iterator it = firstSeqNames.find(sname);
                        
                        if(it != firstSeqNames.end()) { //this is the start of a new chunk
                            double pos = inAmino.tellg();
                            afileFilePos.push_back(pos - input.length() - 1);
                            firstSeqNames.erase(it);
                        }
                    }
                }
                
                if (firstSeqNames.size() == 0) { break; }
            }
            inAmino.close();
            
            
            if (firstSeqNames.size() != 0) {
                for (map<string, int>::iterator it = firstSeqNames.begin(); it != firstSeqNames.end(); it++) {
                    m->mothurOut(it->first + " is in your fasta file and not in your amino file, aborting.\n"); m->setControl_pressed(true);
                }
                return false;
            }
            
            //get last file position of qfile
            FILE * pFile;
            double size;
            
            //get num bytes in file
            aminofile = util.getFullPathName(aminofile);
            pFile = fopen (aminofile.c_str(),"rb");
            if (pFile==nullptr) perror ("Error opening file");
            else{
                fseek (pFile, 0, SEEK_END);
                size=ftell (pFile);
                fclose (pFile);
            }
            
            afileFilePos.push_back(size);
        }
        
        for (int i = 0; i < (fastaFilePos.size()-1); i++) {
            if (m->getDebug()) { m->mothurOut("[DEBUG]: " + toString(i) +'\t' + toString(fastaFilePos[i]) + '\t' + toString(fastaFilePos[i+1]) + '\n'); }
            lines.push_back(linePair(fastaFilePos[i], fastaFilePos[(i+1)]));
            if (aminofile != "") {  aLines.push_back(linePair(afileFilePos[i], afileFilePos[(i+1)]));  }
        }
        
#else
        
        long long numFastaSeqs = 0;
        fastaFilePos = util.setFilePosFasta(fastafile, numFastaSeqs);
        if (numFastaSeqs < processors) { processors = numFastaSeqs; }
        
        if (aminofile != "") {
            long long numAminoSeqs = 0;
            afileFilePos = util.setFilePosFasta(aminofile, numAminoSeqs);
            
            if (numFastaSeqs != numAminoSeqs) {
                m->mothurOut("[ERROR]: You have " + toString(numFastaSeqs) + " sequences in your fasta file, but " + toString(numAminoSeqs) + " sequences in your amino file, please correct.\n");  m->setControl_pressed(true); return false;
            }
        }
        
        //figure out how many sequences you have to process
        int numSeqsPerProcessor = numFastaSeqs / processors;
        for (int i = 0; i < processors; i++) {
            int startIndex =  i * numSeqsPerProcessor;
            if(i == (processors - 1)){    numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor;     }
            lines.push_back(linePair(fastaFilePos[startIndex], numSeqsPerProcessor));
            if (aminofile != "") {  aLines.push_back(linePair(afileFilePos[startIndex], numSeqsPerProcessor)); }
        }
        
#endif
        
        return true;
    }
    catch(exception& e) {
        m->errorOut(e, "TrimSeqsCommand", "setLines");
        exit(1);
    }
}
//***************************************************************************************************************
