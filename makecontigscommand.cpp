//
//  makecontigscommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/15/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "makecontigscommand.h"

//**********************************************************************************************************************
vector<string> MakeContigsCommand::setParameters(){	
	try {
		CommandParameter pfasta("ffastq", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pfasta);
        CommandParameter prfasta("rfastq", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(prfasta);
		CommandParameter palign("align", "Multiple", "needleman-gotoh", "needleman", "", "", "",false,false); parameters.push_back(palign);
		CommandParameter pmatch("match", "Number", "", "1.0", "", "", "",false,false); parameters.push_back(pmatch);
		CommandParameter pmismatch("mismatch", "Number", "", "-1.0", "", "", "",false,false); parameters.push_back(pmismatch);
		CommandParameter pgapopen("gapopen", "Number", "", "-2.0", "", "", "",false,false); parameters.push_back(pgapopen);
		CommandParameter pgapextend("gapextend", "Number", "", "-1.0", "", "", "",false,false); parameters.push_back(pgapextend);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "",false,false); parameters.push_back(pprocessors);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeContigsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string MakeContigsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The make.contigs command reads a forward fastq file and a reverse fastq file and outputs new fasta and quality files.\n";
		helpString += "The make.contigs command parameters are ffastq, rfastq, align, match, mismatch, gapopen, gapextend and processors.\n";
		helpString += "The ffastq and rfastq parameter is required.\n";
		helpString += "The align parameter allows you to specify the alignment method to use.  Your options are: gotoh, needleman, blast and noalign. The default is needleman.\n";
		helpString += "The match parameter allows you to specify the bonus for having the same base. The default is 1.0.\n";
		helpString += "The mistmatch parameter allows you to specify the penalty for having different bases.  The default is -1.0.\n";
		helpString += "The gapopen parameter allows you to specify the penalty for opening a gap in an alignment. The default is -2.0.\n";
		helpString += "The gapextend parameter allows you to specify the penalty for extending a gap in an alignment.  The default is -1.0.\n";
        helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n";
        helpString += "The make.contigs command should be in the following format: \n";
		helpString += "make.contigs(ffastq=yourForwardFastqFile, rfastq=yourReverseFastqFile, align=yourAlignmentMethod) \n";
		helpString += "Note: No spaces between parameter labels (i.e. ffastq), '=' and parameters (i.e.yourForwardFastqFile).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeContigsCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
MakeContigsCommand::MakeContigsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["qfile"] = tempOutNames;
        outputTypes["mismatch"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeContigsCommand", "MakeContigsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
MakeContigsCommand::MakeContigsCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
        
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters(); 
			
			ValidParameters validParameter("pairwise.seqs");
			map<string, string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["qfile"] = tempOutNames;
            outputTypes["mismatch"] = tempOutNames;
			
            
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else { 
				string path;
                it = parameters.find("ffastq");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["ffastq"] = inputDir + it->second;		}
				}
                
                it = parameters.find("rfastq");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["rfastq"] = inputDir + it->second;		}
				}
            }
            
            ffastqfile = validParameter.validFile(parameters, "ffastq", true);
			if (ffastqfile == "not open") { ffastqfile = ""; abort = true; }	
			else if (ffastqfile == "not found") { ffastqfile = ""; abort=true;  m->mothurOut("The ffastq parameter is required.\n"); }
			
			rfastqfile = validParameter.validFile(parameters, "rfastq", true);
			if (rfastqfile == "not open") { rfastqfile = ""; abort = true; }	
			else if (rfastqfile == "not found") { rfastqfile = ""; abort=true;  m->mothurOut("The rfastq parameter is required.\n"); }
            
            //if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(ffastqfile);		}
			

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.validFile(parameters, "match", false);		if (temp == "not found"){	temp = "1.0";			}
			m->mothurConvert(temp, match);  
			
			temp = validParameter.validFile(parameters, "mismatch", false);		if (temp == "not found"){	temp = "-1.0";			}
			m->mothurConvert(temp, misMatch);  
            if (misMatch > 0) { m->mothurOut("[ERROR]: mismatch must be negative.\n"); abort=true; }
			
			temp = validParameter.validFile(parameters, "gapopen", false);		if (temp == "not found"){	temp = "-2.0";			}
			m->mothurConvert(temp, gapOpen);  
            if (gapOpen > 0) { m->mothurOut("[ERROR]: gapopen must be negative.\n"); abort=true; }
			
			temp = validParameter.validFile(parameters, "gapextend", false);	if (temp == "not found"){	temp = "-1.0";			}
			m->mothurConvert(temp, gapExtend); 
            if (gapExtend > 0) { m->mothurOut("[ERROR]: gapextend must be negative.\n"); abort=true; }
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
			
			align = validParameter.validFile(parameters, "align", false);		if (align == "not found"){	align = "needleman";	}
			if ((align != "needleman") && (align != "gotoh")) { m->mothurOut(align + " is not a valid alignment method. Options are needleman or gotoh. I will use needleman."); m->mothurOutEndLine(); align = "needleman"; }
        }
		
	}
	catch(exception& e) {
		m->errorOut(e, "MakeContigsCommand", "MakeContigsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int MakeContigsCommand::execute(){
	try {
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        //read ffastq and rfastq files creating fasta and qual files.
        //this function will create a forward and reverse, fasta and qual files for each processor.
        //files has an entry for each processor. files[i][0] = forwardFasta, files[i][1] = forwardQual, files[i][2] = reverseFasta, files[i][3] = reverseQual
        int numReads = 0;
        int start = time(NULL);
        longestBase = 1000;
        m->mothurOut("Reading fastq data...\n"); 
        vector< vector<string> > files = readFastqFiles(numReads);  
        m->mothurOut("Done.\n");
    
        if (m->control_pressed) { return 0; }
        
        string outFastaFile = outputDir + m->getRootName(m->getSimpleName(ffastqfile)) + "contigs.fasta";
        string outQualFile = outputDir + m->getRootName(m->getSimpleName(ffastqfile)) + "contigs.qual";
        string outMisMatchFile = outputDir + m->getRootName(m->getSimpleName(ffastqfile)) + "contigs.mismatches";
        outputNames.push_back(outFastaFile); outputTypes["fasta"].push_back(outFastaFile);
        outputNames.push_back(outQualFile); outputTypes["qfile"].push_back(outQualFile);
        outputNames.push_back(outMisMatchFile); outputTypes["mismatch"].push_back(outMisMatchFile);
        
        m->mothurOut("Making contigs...\n"); 
        createProcesses(files, outFastaFile, outQualFile, outMisMatchFile);
        m->mothurOut("Done.\n");
        
        //remove temp fasta and qual files
        for (int i = 0; i < processors; i++) { for(int j = 0; j < files[i].size(); j++) { m->mothurRemove(files[i][j]); }  }
        
        if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); }  return 0; }
        
        m->mothurOut("It took " + toString(time(NULL) - start) + " secs to process " + toString(numReads) + " sequences.\n");
        
        string currentFasta = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentFasta = (itTypes->second)[0]; m->setFastaFile(currentFasta); }
		}
        
        string currentQual = "";
		itTypes = outputTypes.find("qfile");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentQual = (itTypes->second)[0]; m->setQualFile(currentQual); }
		}
		
        //output files created by command
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "MakeContigsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int MakeContigsCommand::createProcesses(vector< vector<string> > files, string outputFasta, string outputQual, string outputMisMatches) {
	try {
		int num = 0;
		vector<int> processIDS;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		int process = 0;
		
		//loop through and create all the processes you want
		while (process != processors-1) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driver(files[process], outputFasta + toString(getpid()) + ".temp", outputQual + toString(getpid()) + ".temp", outputMisMatches + toString(getpid()) + ".temp");
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = outputFasta + toString(getpid()) + ".num.temp";
				m->openOutputFile(tempFile, out);
				out << num << endl;
				out.close();
				
				exit(0);
			}else { 
				m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
				for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
				exit(0);
			}
		}
		
		//do my part
		num = driver(files[processors-1], outputFasta, outputQual, outputMisMatches);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
        
		for (int i = 0; i < processIDS.size(); i++) {
			ifstream in;
			string tempFile =  outputFasta + toString(processIDS[i]) + ".num.temp";
			m->openInputFile(tempFile, in);
			if (!in.eof()) { int tempNum = 0; in >> tempNum; num += tempNum; }
			in.close(); m->mothurRemove(tempFile);
        }
    #else
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the contigsData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<contigsData*> pDataArray; 
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
		
		//Create processor worker threads.
		for( int i=0; i<processors-1; i++ ){
			string extension = toString(i) + ".temp";
			
			contigsData* tempcontig = new contigsData(files[i], (outputFasta + extension), (outputQual + extension), (outputMisMatches + extension), align, m, match, misMatch, gapOpen, gapExtend, i);
			pDataArray.push_back(tempcontig);
			processIDS.push_back(i);
            
			hThreadArray[i] = CreateThread(NULL, 0, MyContigsThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);   
		}
		
        num = driver(files[processors-1], outputFasta, outputQual, outputMisMatches);		
        
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
			num += pDataArray[i]->count;
            CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
        }
				
    #endif	
        
        for (int i = 0; i < processIDS.size(); i++) {
			m->appendFiles((outputFasta + toString(processIDS[i]) + ".temp"), outputFasta);
			m->mothurRemove((outputFasta + toString(processIDS[i]) + ".temp"));
			
			m->appendFiles((outputQual + toString(processIDS[i]) + ".temp"), outputQual);
			m->mothurRemove((outputQual + toString(processIDS[i]) + ".temp"));
            
            m->appendFiles((outputMisMatches + toString(processIDS[i]) + ".temp"), outputMisMatches);
			m->mothurRemove((outputMisMatches + toString(processIDS[i]) + ".temp"));
		}
		
		return num;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeContigsCommand", "createProcesses");
		exit(1);
	}
}
//**********************************************************************************************************************
int MakeContigsCommand::driver(vector<string> files, string outputFasta, string outputQual, string outputMisMatches){
    try {
        
        Alignment* alignment;
        if(align == "gotoh")			{	alignment = new GotohOverlap(gapOpen, gapExtend, match, misMatch, longestBase);			}
		else if(align == "needleman")	{	alignment = new NeedlemanOverlap(gapOpen, match, misMatch, longestBase);				}
        
        int num = 0;
        string thisffastafile = files[0];
        string thisfqualfile = files[1];
        string thisrfastafile = files[2];
        string thisrqualfile = files[3];
        
        if (m->debug) {  m->mothurOut("[DEBUG]: ffasta = " + thisffastafile + ".\n[DEBUG]: fqual = " + thisfqualfile + ".\n[DEBUG]: rfasta = " + thisrfastafile + ".\n[DEBUG]: rqual = " + thisrqualfile + ".\n"); }
        
        ifstream inFFasta, inRFasta, inFQual, inRQual;
        m->openInputFile(thisffastafile, inFFasta);
        m->openInputFile(thisfqualfile, inFQual);
        m->openInputFile(thisrfastafile, inRFasta);
        m->openInputFile(thisrqualfile, inRQual);
        
        ofstream outFasta, outQual, outMisMatch;
        m->openOutputFile(outputFasta, outFasta);
        m->openOutputFile(outputQual, outQual);
        m->openOutputFile(outputMisMatches, outMisMatch);
        outMisMatch << "Name\tLength\tMisMatches\n";
        
        while ((!inFQual.eof()) && (!inFFasta.eof()) && (!inRFasta.eof()) && (!inRQual.eof())) {
            
            if (m->control_pressed) { break; }
            
            //read seqs and quality info
            Sequence fSeq(inFFasta); m->gobble(inFFasta);
            Sequence rSeq(inRFasta); m->gobble(inRFasta);
            QualityScores fQual(inFQual); m->gobble(inFQual);
            QualityScores rQual(inRQual); m->gobble(inRQual);
            
            //flip the reverse reads
            rSeq.reverseComplement();
            rQual.flipQScores();
            
            //pairwise align
            alignment->align(fSeq.getUnaligned(), rSeq.getUnaligned());
            map<int, int> ABaseMap = alignment->getSeqAAlnBaseMap();
            map<int, int> BBaseMap = alignment->getSeqBAlnBaseMap();
            fSeq.setAligned(alignment->getSeqAAln());
            rSeq.setAligned(alignment->getSeqBAln());
            int length = fSeq.getAligned().length();
        
            //traverse alignments merging into one contiguous seq
            string contig = "";
            vector<int> contigScores; 
            int numMismatches = 0;
            string seq1 = fSeq.getAligned();
            string seq2 = rSeq.getAligned();
           
            vector<int> scores1 = fQual.getQualityScores();
            vector<int> scores2 = rQual.getQualityScores();
            
            for (int i = 0; i < length; i++) {
                if (seq1[i] == seq2[i]) { //match, add base and choose highest score
                    contig += seq1[i];
                    contigScores.push_back(scores1[ABaseMap[i]]);
                    if (scores1[ABaseMap[i]] < scores2[BBaseMap[i]]) { contigScores[i] = scores2[BBaseMap[i]]; }
                }else if (((seq1[i] == '.') || (seq1[i] == '-')) && ((seq2[i] != '-') && (seq2[i] != '.'))) { //seq1 is a gap and seq2 is a base, choose seq2
                    contig += seq2[i];
                    contigScores.push_back(scores2[BBaseMap[i]]);
                }else if (((seq2[i] == '.') || (seq2[i] == '-')) && ((seq1[i] != '-') && (seq1[i] != '.'))) { //seq2 is a gap and seq1 is a base, choose seq1
                    contig += seq1[i];
                    contigScores.push_back(scores1[ABaseMap[i]]);
                }else if (((seq1[i] != '-') && (seq1[i] != '.')) && ((seq2[i] != '-') && (seq2[i] != '.'))) { //both bases choose one with better quality
                    char c = seq1[i];
                    contigScores.push_back(scores1[ABaseMap[i]]);
                    if (scores1[ABaseMap[i]] < scores2[BBaseMap[i]]) { contigScores[i] = scores2[BBaseMap[i]]; c = seq2[i]; }
                    contig += c;
                    numMismatches++;
                }else { //should never get here
                    m->mothurOut("[ERROR]: case I didn't think of seq1 = " + toString(seq1[i]) + " and seq2 = " + toString(seq2[i]) + "\n");
                }
            }
            
            //output
            outFasta << ">" << fSeq.getName() << endl << contig << endl;
            outQual << ">" << fSeq.getName() << endl;
            for (int i = 0; i < contigScores.size(); i++) { outQual << contigScores[i] << ' '; }
            outQual << endl;
            outMisMatch << fSeq.getName() << '\t' << contig.length() << '\t' << numMismatches << endl;
            
            num++;
            
			//report progress
			if((num) % 1000 == 0){	m->mothurOut(toString(num)); m->mothurOutEndLine();		}
		}
        
		//report progress
		if((num) % 1000 != 0){	m->mothurOut(toString(num)); m->mothurOutEndLine();		}
        
        inFFasta.close();
        inFQual.close();
        inRFasta.close();
        inRQual.close();
        outFasta.close();
        outQual.close();
        outMisMatch.close();
        delete alignment;
        
        if (m->control_pressed) { m->mothurRemove(outputQual); m->mothurRemove(outputFasta);  m->mothurRemove(outputMisMatches);}
        
        return num;
    }
	catch(exception& e) {
		m->errorOut(e, "MakeContigsCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
vector< vector<string> > MakeContigsCommand::readFastqFiles(int& count){
    try {
        vector< vector<string> > files;
        
        //maps processors number to file pointer
        map<int, vector<ofstream*> > tempfiles;  //tempfiles[0] = forwardFasta, [1] = forwardQual, [2] = reverseFasta, [3] = reverseQual
        map<int, vector<ofstream*> >::iterator it;
        
        //create files to write to
        for (int i = 0; i < processors; i++) {
            vector<ofstream*> temp;
            ofstream* outFF = new ofstream;     temp.push_back(outFF);
            ofstream* outFQ = new ofstream;     temp.push_back(outFQ);
            ofstream* outRF = new ofstream;     temp.push_back(outRF);
            ofstream* outRQ = new ofstream;     temp.push_back(outRQ);
            tempfiles[i] = temp;
            
            vector<string> names;
            string ffastafilename = outputDir + m->getRootName(m->getSimpleName(ffastqfile)) + toString(i) + "ffasta.temp";
            string rfastafilename = outputDir + m->getRootName(m->getSimpleName(rfastqfile)) + toString(i) + "rfasta.temp";
            string fqualfilename = outputDir + m->getRootName(m->getSimpleName(ffastqfile)) + toString(i) + "fqual.temp";
            string rqualfilename = outputDir + m->getRootName(m->getSimpleName(rfastqfile)) + toString(i) + "rqual.temp";
            names.push_back(ffastafilename); names.push_back(fqualfilename);
            names.push_back(rfastafilename); names.push_back(rqualfilename);
            files.push_back(names);
            
            m->openOutputFile(ffastafilename, *outFF);
            m->openOutputFile(rfastafilename, *outRF);
            m->openOutputFile(fqualfilename, *outFQ);
            m->openOutputFile(rqualfilename, *outRQ);
        }
        
        if (m->control_pressed) {
            //close files, delete ofstreams
            for (it = tempfiles.begin(); it!=tempfiles.end(); it++) { for (int i = 0; i < (it->second).size(); i++) { (*(it->second)[i]).close();  delete (it->second)[i]; } }
            //remove files
            for (int i = 0; i < files.size(); i++) {  
                for(int j = 0; j < files[i].size(); j++) { m->mothurRemove(files[i][j]); }
            }
        }
        
        ifstream inForward;
        m->openInputFile(ffastqfile, inForward);
        
        ifstream inReverse;
        m->openInputFile(rfastqfile, inReverse);
        
        count = 0;
        while ((!inForward.eof()) && (!inReverse.eof())) {
            
            if (m->control_pressed) { for (it = tempfiles.begin(); it!=tempfiles.end(); it++) { for (int i = 0; i < (it->second).size(); i++) { (*(it->second)[i]).close();  delete (it->second)[i]; } } for (int i = 0; i < files.size(); i++) {  for(int j = 0; j < files[i].size(); j++) { m->mothurRemove(files[i][j]); } } inForward.close(); inReverse.close(); return files; }
            
            //get a read from forward and reverse fastq files
            fastqRead fread = readFastq(inForward);
            fastqRead rread = readFastq(inReverse);
            checkReads(fread, rread);
            
            if (m->control_pressed) { for (it = tempfiles.begin(); it!=tempfiles.end(); it++) { for (int i = 0; i < (it->second).size(); i++) { (*(it->second)[i]).close();  delete (it->second)[i]; } } for (int i = 0; i < files.size(); i++) {  for(int j = 0; j < files[i].size(); j++) { m->mothurRemove(files[i][j]); } } inForward.close(); inReverse.close(); return files; }
            
            //if the reads are okay write to output files
            int process = count % processors;
            
            *(tempfiles[process][0]) << ">" << fread.name << endl << fread.sequence << endl;
            *(tempfiles[process][1]) << ">" << fread.name << endl;
            for (int i = 0; i < fread.scores.size(); i++) { *(tempfiles[process][1]) << fread.scores[i] << " "; }
            *(tempfiles[process][1]) << endl;
            *(tempfiles[process][2]) << ">" << rread.name << endl << rread.sequence << endl;
            *(tempfiles[process][3]) << ">" << rread.name << endl;
            for (int i = 0; i < rread.scores.size(); i++) { *(tempfiles[process][3]) << rread.scores[i] << " "; }
            *(tempfiles[process][3]) << endl;
            
            count++;
            
            //report progress
			if((count) % 10000 == 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}
			
		}
		//report progress
		if((count) % 10000 != 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}
		
        
        
        //close files, delete ofstreams
        for (it = tempfiles.begin(); it!=tempfiles.end(); it++) { for (int i = 0; i < (it->second).size(); i++) { (*(it->second)[i]).close();  delete (it->second)[i]; } }
        inForward.close();
        inReverse.close();
        
        //adjust for really large processors or really small files
        if (count < processors) { 
            for (int i = count; i < processors; i++) { for(int j = 0; j < files[i].size(); j++) { m->mothurRemove(files[i][j]); } files[i].clear(); }
            files.resize(count);
            processors = count; 
        }
        
        return files;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "readFastqFiles");
        exit(1);
    }
}
//**********************************************************************************************************************
fastqRead MakeContigsCommand::readFastq(ifstream& in){
    try {
        fastqRead read;
        
        //read sequence name
        string name = m->getline(in); m->gobble(in);
        if (name == "") {  m->mothurOut("[ERROR]: Blank fasta name."); m->mothurOutEndLine(); m->control_pressed = true; return read; }
        else if (name[0] != '@') { m->mothurOut("[ERROR]: reading " + name + " expected a name with @ as a leading character."); m->mothurOutEndLine(); m->control_pressed = true; return read; }
        else { name = name.substr(1); }
        
        //read sequence
        string sequence = m->getline(in); m->gobble(in);
        if (sequence == "") {  m->mothurOut("[ERROR]: missing sequence for " + name); m->mothurOutEndLine(); m->control_pressed = true; return read; }
        
        //read sequence name
        string name2 = m->getline(in); m->gobble(in);
        if (name2 == "") {  m->mothurOut("[ERROR]: Blank quality name."); m->mothurOutEndLine(); m->control_pressed = true; return read; }
        else if (name2[0] != '+') { m->mothurOut("[ERROR]: reading " + name2 + " expected a name with + as a leading character."); m->mothurOutEndLine(); m->control_pressed = true; return read; }
        else { name2 = name2.substr(1);  }
        
        //read quality scores
        string quality = m->getline(in); m->gobble(in);
        if (quality == "") {  m->mothurOut("[ERROR]: missing quality for " + name2); m->mothurOutEndLine(); m->control_pressed = true; return read; }
        
        //sanity check sequence length and number of quality scores match
        if (name2 != "") { if (name != name2) { m->mothurOut("[ERROR]: names do not match. read " + name + " for fasta and " + name2 + " for quality."); m->mothurOutEndLine(); m->control_pressed = true; return read; } }
        if (quality.length() != sequence.length()) { m->mothurOut("[ERROR]: Lengths do not match for sequence " + name + ". Read " + toString(sequence.length()) + " characters for fasta and " + toString(quality.length()) + " characters for quality scores."); m->mothurOutEndLine(); m->control_pressed = true; return read; }
        
        vector<int> qualScores;
		int controlChar = int('@');
		for (int i = 0; i < quality.length(); i++) { 
			int temp = int(quality[i]);
			temp -= controlChar;
			
			qualScores.push_back(temp);
		}

        read.name = name;
        read.sequence = sequence;
        read.scores = qualScores;

        return read;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "readFastq");
        exit(1);
    }
}
//**********************************************************************************************************************
bool MakeContigsCommand::checkReads(fastqRead& forward, fastqRead& reverse){
    try {
        bool good = true;
        
        //fix names
        if ((forward.name.length() > 2) && (reverse.name.length() > 2)) {
            forward.name = forward.name.substr(0, forward.name.length()-2);
            reverse.name = reverse.name.substr(0, reverse.name.length()-2);
        }else { good = false; m->control_pressed = true; }
        
        //do names match
        if (forward.name != reverse.name) {
            m->mothurOut("[ERROR]: read " + forward.name + " from " + ffastqfile + ", but read " + reverse.name + " from " + rfastqfile + ".\n");
            good = false; m->control_pressed = true;
        }
        
        //do sequence lengths match
        if (forward.sequence.length() != reverse.sequence.length()) {
            m->mothurOut("[ERROR]: For sequence " + forward.name + " I read a sequence of length " + toString(forward.sequence.length()) + " from " + ffastqfile + ", but read a sequence of length " + toString(reverse.sequence.length()) + " from " + rfastqfile + ".\n");
            good = false; m->control_pressed = true;
        }
        
        //do number of qual scores match 
        if (forward.scores.size() != reverse.scores.size()) {
            m->mothurOut("[ERROR]: For sequence " + forward.name + " I read " + toString(forward.scores.size()) + " quality scores from " + ffastqfile + ", but read  " + toString(reverse.scores.size()) + " quality scores from " + rfastqfile + ".\n");
            good = false; m->control_pressed = true;
        }

        return good;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "readFastq");
        exit(1);
    }
}
//**********************************************************************************************************************




