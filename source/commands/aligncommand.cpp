/*
 *  aligncommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/15/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 *	This version of nast does everything I think that the greengenes nast server does and then some.  I have added the 
 *	feature of allowing users to define their database, kmer size for searching, alignment penalty values and alignment 
 *	method.  This latter feature is perhaps most significant.  nastPlus enables a user to use either a Needleman-Wunsch 
 *	(non-affine gap penalty) or Gotoh (affine gap penalty) pairwise alignment algorithm.  This is significant because it
 *	allows for a global alignment and not the local alignment provided by bLAst.  Furthermore, it has the potential to
 *	provide a better alignment because of the banding method employed by blast (I'm not sure about this).
 *
 */

#include "aligncommand.h"

//**********************************************************************************************************************
vector<string> AlignCommand::setParameters(){	
	try {
		CommandParameter ptemplate("reference", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(ptemplate);
		CommandParameter pcandidate("fasta", "InputTypes", "", "", "none", "none", "none","fasta-alignreport-accnos",false,true,true); parameters.push_back(pcandidate);
		CommandParameter psearch("search", "Multiple", "kmer-blast-suffix", "kmer", "", "", "","",false,false,true); parameters.push_back(psearch);
		CommandParameter pksize("ksize", "Number", "", "8", "", "", "","",false,false); parameters.push_back(pksize);
		CommandParameter pmatch("match", "Number", "", "1.0", "", "", "","",false,false); parameters.push_back(pmatch);
		CommandParameter palign("align", "Multiple", "needleman-gotoh-blast-noalign", "needleman", "", "", "","",false,false,true); parameters.push_back(palign);
		CommandParameter pmismatch("mismatch", "Number", "", "-1.0", "", "", "","",false,false); parameters.push_back(pmismatch);
		CommandParameter pgapopen("gapopen", "Number", "", "-5.0", "", "", "","",false,false); parameters.push_back(pgapopen);
		CommandParameter pgapextend("gapextend", "Number", "", "-2.0", "", "", "","",false,false); parameters.push_back(pgapextend);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pflip("flip", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pflip);
		CommandParameter pthreshold("threshold", "Number", "", "0.50", "", "", "","",false,false); parameters.push_back(pthreshold);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);

		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string AlignCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The align.seqs command reads a file containing sequences and creates an alignment file and a report file.";
		helpString += "The align.seqs command parameters are reference, fasta, search, ksize, align, match, mismatch, gapopen, gapextend and processors.";
		helpString += "The reference and fasta parameters are required. You may leave fasta blank if you have a valid fasta file. You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amzon.fasta.";
		helpString += "The search parameter allows you to specify the method to find most similar template.  Your options are: suffix, kmer and blast. The default is kmer.";
		helpString += "The align parameter allows you to specify the alignment method to use.  Your options are: gotoh, needleman, blast and noalign. The default is needleman.";
		helpString += "The ksize parameter allows you to specify the kmer size for finding most similar template to candidate.  The default is 8.";
		helpString += "The match parameter allows you to specify the bonus for having the same base. The default is 1.0.";
		helpString += "The mistmatch parameter allows you to specify the penalty for having different bases.  The default is -1.0.";
		helpString += "The gapopen parameter allows you to specify the penalty for opening a gap in an alignment. The default is -5.0.";
		helpString += "The gapextend parameter allows you to specify the penalty for extending a gap in an alignment.  The default is -2.0.";
		helpString += "The flip parameter is used to specify whether or not you want mothur to try the reverse complement if a sequence falls below the threshold.  The default is false.";
		helpString += "The threshold is used to specify a cutoff at which an alignment is deemed 'bad' and the reverse complement may be tried. The default threshold is 0.50, meaning 50% of the bases are removed in the alignment.";
		helpString += "If the flip parameter is set to true the reverse complement of the sequence is aligned and the better alignment is reported. Default=t";
		helpString += "The default for the threshold parameter is 0.50, meaning at least 50% of the bases must remain or the sequence is reported as potentially reversed.";
		helpString += "The align.seqs command should be in the following format:";
		helpString += "align.seqs(reference=yourTemplateFile, fasta=yourCandidateFile, align=yourAlignmentMethod, search=yourSearchmethod, ksize=yourKmerSize, match=yourMatchBonus, mismatch=yourMismatchpenalty, gapopen=yourGapopenPenalty, gapextend=yourGapExtendPenalty)";
		helpString += "Example align.seqs(candidate=candidate.fasta, template=core.filtered, align=kmer, search=gotoh, ksize=8, match=2.0, mismatch=3.0, gapopen=-2.0, gapextend=-1.0)";
		helpString += "Note: No spaces between parameter labels (i.e. candidate), '=' and parameters (i.e.yourFastaFile).";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string AlignCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],align"; } //makes file like: amazon.align
        else if (type == "alignreport") {  pattern = "[filename],align.report"; }
        else if (type == "accnos") {  pattern = "[filename],flip.accnos"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "AlignCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
AlignCommand::AlignCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["alignreport"] = tempOutNames;
		outputTypes["accnos"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "AlignCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
AlignCommand::AlignCommand(string option)  {
	try {
		abort = false; calledHelp = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true;}
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters(); 
			
			ValidParameters validParameter("align.seqs");
			map<string, string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["alignreport"] = tempOutNames;
			outputTypes["accnos"] = tempOutNames;
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			

			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;

				it = parameters.find("reference");

				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["reference"] = inputDir + it->second;		}
				}
			}

			candidateFileName = validParameter.validFile(parameters, "fasta", false);
			if (candidateFileName == "not found") { 
				//if there is a current fasta file, use it
				string filename = m->getFastaFile(); 
				if (filename != "") { candidateFileNames.push_back(filename); m->mothurOut("Using " + filename + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the candidate parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { 
				m->splitAtDash(candidateFileName, candidateFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < candidateFileNames.size(); i++) {
					//candidateFileNames[i] = m->getFullPathName(candidateFileNames[i]);
					
					bool ignore = false;
					if (candidateFileNames[i] == "current") { 
						candidateFileNames[i] = m->getFastaFile(); 
						if (candidateFileNames[i] != "") {  m->mothurOut("Using " + candidateFileNames[i] + " as input file for the fasta parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current fastafile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							candidateFileNames.erase(candidateFileNames.begin()+i);
							i--;
						}
					}
					
					if (!ignore) {
					
						if (inputDir != "") {
							string path = m->hasPath(candidateFileNames[i]);
							//if the user has not given a path then, add inputdir. else leave path alone.
							if (path == "") {	candidateFileNames[i] = inputDir + candidateFileNames[i];		}
						}
		
						int ableToOpen;
						ifstream in;
						ableToOpen = m->openInputFile(candidateFileNames[i], in, "noerror");
						in.close();	
						
						//if you can't open it, try default location
						if (ableToOpen == 1) {
							if (m->getDefaultPath() != "") { //default path is set
								string tryPath = m->getDefaultPath() + m->getSimpleName(candidateFileNames[i]);
								m->mothurOut("Unable to open " + candidateFileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								candidateFileNames[i] = tryPath;
							}
						}
						
						//if you can't open it, try output location
						if (ableToOpen == 1) {
							if (m->getOutputDir() != "") { //default path is set
								string tryPath = m->getOutputDir() + m->getSimpleName(candidateFileNames[i]);
								m->mothurOut("Unable to open " + candidateFileNames[i] + ". Trying output directory " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								candidateFileNames[i] = tryPath;
							}
						}
						
										

						if (ableToOpen == 1) { 
							m->mothurOut("Unable to open " + candidateFileNames[i] + ". It will be disregarded."); m->mothurOutEndLine(); 
							//erase from file list
							candidateFileNames.erase(candidateFileNames.begin()+i);
							i--;
						}else {
							m->setFastaFile(candidateFileNames[i]);
						}
					}
				}
				
				//make sure there is at least one valid file left
				if (candidateFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
		
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.validFile(parameters, "ksize", false);		if (temp == "not found"){	temp = "8";				}
			m->mothurConvert(temp, kmerSize); 
			
			temp = validParameter.validFile(parameters, "match", false);		if (temp == "not found"){	temp = "1.0";			}
			m->mothurConvert(temp, match);  
			
			temp = validParameter.validFile(parameters, "mismatch", false);		if (temp == "not found"){	temp = "-1.0";			}
			m->mothurConvert(temp, misMatch);  
			
			temp = validParameter.validFile(parameters, "gapopen", false);		if (temp == "not found"){	temp = "-5.0";			}
			m->mothurConvert(temp, gapOpen);  
			
			temp = validParameter.validFile(parameters, "gapextend", false);	if (temp == "not found"){	temp = "-2.0";			}
			m->mothurConvert(temp, gapExtend); 
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors); 
			
			temp = validParameter.validFile(parameters, "flip", false);			if (temp == "not found"){	temp = "t";				}
			flip = m->isTrue(temp);
			
			//this has to go after save so that if the user sets save=t and provides no reference we abort
			templateFileName = validParameter.validFile(parameters, "reference", true);
			if (templateFileName == "not found") { m->mothurOut("[ERROR]: The reference parameter is a required for the align.seqs command, aborting.\n"); abort = true;
			}else if (templateFileName == "not open") { abort = true; }	
			
			temp = validParameter.validFile(parameters, "threshold", false);	if (temp == "not found"){	temp = "0.50";			}
			m->mothurConvert(temp, threshold); 
			
			search = validParameter.validFile(parameters, "search", false);		if (search == "not found"){	search = "kmer";		}
			if ((search != "suffix") && (search != "kmer") && (search != "blast")) { m->mothurOut("invalid search option: choices are kmer, suffix or blast."); m->mothurOutEndLine(); abort=true; }
			
			align = validParameter.validFile(parameters, "align", false);		if (align == "not found"){	align = "needleman";	}
			if ((align != "needleman") && (align != "gotoh") && (align != "blast") && (align != "noalign")) { m->mothurOut("invalid align option: choices are needleman, gotoh, blast or noalign."); m->mothurOutEndLine(); abort=true; }

		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "AlignCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
AlignCommand::~AlignCommand(){	

	if (abort == false) {
		for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
		delete templateDB;
	}
}
//**********************************************************************************************************************

int AlignCommand::execute(){
	try {
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}

		templateDB = new AlignmentDB(templateFileName, search, kmerSize, gapOpen, gapExtend, match, misMatch, m->getRandomNumber(), true);
		
		for (int s = 0; s < candidateFileNames.size(); s++) {
			if (m->getControl_pressed()) { outputTypes.clear(); return 0; }
			
			m->mothurOut("Aligning sequences from " + candidateFileNames[s] + " ..." ); m->mothurOutEndLine();
			
			if (outputDir == "") {  outputDir += m->hasPath(candidateFileNames[s]); }
            map<string, string> variables; variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(candidateFileNames[s]));
			string alignFileName = getOutputFileName("fasta", variables);  
			string reportFileName = getOutputFileName("alignreport", variables);
			string accnosFileName = getOutputFileName("accnos", variables);
            
			bool hasAccnos = true;
            vector<long long> numFlipped;
            numFlipped.push_back(0); //numflipped because reverse was better
            numFlipped.push_back(0); //total number of sequences with over 50% of bases removed
			
			int numFastaSeqs = 0;
			for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();

			vector<unsigned long long> positions; 
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
			positions = m->divideFile(candidateFileNames[s], processors);
			for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(new linePair(positions[i], positions[(i+1)]));	}
		#else
            positions = m->setFilePosFasta(candidateFileNames[s], numFastaSeqs);
            if (numFastaSeqs < processors) { processors = numFastaSeqs; }
            
            //figure out how many sequences you have to process
            int numSeqsPerProcessor = numFastaSeqs / processors;
            for (int i = 0; i < processors; i++) {
                int startIndex =  i * numSeqsPerProcessor;
                if(i == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor; 	}
                lines.push_back(new linePair(positions[startIndex], numSeqsPerProcessor));
            }
		#endif
			
            numFastaSeqs = createProcesses(alignFileName, reportFileName, accnosFileName, candidateFileNames[s], numFlipped);
				
			if (m->getControl_pressed()) { m->mothurRemove(accnosFileName); m->mothurRemove(alignFileName); m->mothurRemove(reportFileName); outputTypes.clear();  return 0; }
			
			//delete accnos file if its blank else report to user
			if (m->isBlank(accnosFileName)) {  m->mothurRemove(accnosFileName);  hasAccnos = false; }
			else { 
				m->mothurOut("[WARNING]: " + toString(numFlipped[1]) + " of your sequences generated alignments that eliminated too many bases, a list is provided in " + accnosFileName + ".");
				if (!flip) {
					m->mothurOut(" If you set the flip parameter to true mothur will try aligning the reverse compliment as well. flip=t");
				}else{  m->mothurOut("[NOTE]: " + toString(numFlipped[0]) + " of your sequences were reversed to produce a better alignment.");  }
				m->mothurOutEndLine();
			}

			outputNames.push_back(alignFileName); outputTypes["fasta"].push_back(alignFileName);
			outputNames.push_back(reportFileName); outputTypes["alignreport"].push_back(reportFileName);
			if (hasAccnos)	{	outputNames.push_back(accnosFileName);	outputTypes["accnos"].push_back(accnosFileName);  }
		}
		
		//set align file as new current fastafile
		string currentFasta = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentFasta = (itTypes->second)[0]; m->setFastaFile(currentFasta); }
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "execute");
		exit(1);
	}
}
struct alignStruct {
    string alignFName;
    string reportFName;
    string accnosFName;
    string inputFilename;
    string alignMethod, search;
    float match, misMatch, gapOpen, gapExtend, threshold;
    bool flip;
    long long numSeqs;
    
    vector<long long> flippedResults;
    linePair filePos;
    
    MothurOut* m;
    AlignmentDB* templateDB;
    mutex* mutex;
    
       alignStruct (linePair fP, string aFName, string reFName, string ac, string fname, MothurOut* mo, string al, float ma, float misMa, float gOpen, float gExtend, float thr, bool fl, AlignmentDB* tB, string se) {
        
        filePos.start = fP.start;
        filePos.end = fP.end;
        alignFName = aFName;
        reportFName = reFName;
        accnosFName = ac;
        inputFilename = fname;
        numSeqs = 0;
        m = mo;
        alignMethod = al;
        match = ma;
        misMatch = misMa;
        gapOpen = gOpen;
        gapExtend = gExtend;
        threshold = thr;
        flip = fl;
        templateDB = tB;
        search = se;
        flippedResults.resize(2, 0);
    }
    
};
//**********************************************************************************************************************
void alignDriver(alignStruct* params) {
	try {
		ofstream alignmentFile;
		params->m->openOutputFile(params->alignFName, alignmentFile);
		
		ofstream accnosFile;
		params->m->openOutputFile(params->accnosFName, accnosFile);
		
		NastReport report(params->reportFName);
		
		ifstream inFASTA;
		params->m->openInputFile(params->inputFilename, inFASTA);

		inFASTA.seekg(params->filePos.start);

		bool done = false;
        
		long long count = 0;
        long long numFlipped_0 = 0;
        long long numFlipped_1 = 0;
		
		//moved this into driver to avoid deep copies in windows paralellized version
		Alignment* alignment;
		int longestBase = params->templateDB->getLongestBase();
        if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: template longest base = "  + toString(longestBase) + " \n"); }
		if(params->alignMethod == "gotoh")			{	alignment = new GotohOverlap(params->gapOpen, params->gapExtend, params->match, params->misMatch, longestBase);			}
		else if(params->alignMethod == "needleman")	{	alignment = new NeedlemanOverlap(params->gapOpen, params->match, params->misMatch, longestBase);				}
		else if(params->alignMethod == "blast")		{	alignment = new BlastAlignment(params->gapOpen, params->gapExtend, params->match, params->misMatch);		}
		else if(params->alignMethod == "noalign")		{	alignment = new NoAlign();													}
		else {
			params->m->mothurOut(params->alignMethod + " is not a valid alignment option. I will run the command using needleman.");
			params->m->mothurOutEndLine();
			alignment = new NeedlemanOverlap(params->gapOpen, params->match, params->misMatch, longestBase);
		}
	
		while (!done) {
			
			if (params->m->getControl_pressed()) {  break; }
			
			Sequence* candidateSeq = new Sequence(inFASTA);  params->m->gobble(inFASTA);
			report.setCandidate(candidateSeq);

			int origNumBases = candidateSeq->getNumBases();
			string originalUnaligned = candidateSeq->getUnaligned();
			int numBasesNeeded = origNumBases * params->threshold;
	
			if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
				if (candidateSeq->getUnaligned().length()+1 > alignment->getnRows()) {
                    if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: " + candidateSeq->getName() + " " + toString(candidateSeq->getUnaligned().length()) + " " + toString(alignment->getnRows()) + " \n"); }
					alignment->resize(candidateSeq->getUnaligned().length()+2);
				}
                
                float searchScore;
				Sequence temp = params->templateDB->findClosestSequence(candidateSeq, searchScore);
				Sequence* templateSeq = new Sequence(temp.getName(), temp.getAligned());
								
				Nast* nast = new Nast(alignment, candidateSeq, templateSeq);
		
				Sequence* copy;
				
				Nast* nast2;
				bool needToDeleteCopy = false;  //this is needed in case you have you enter the ifs below
												//since nast does not make a copy of hte sequence passed, and it is used by the reporter below
												//you can't delete the copy sequence til after you report, but you may choose not to create it in the first place
												//so this bool tells you if you need to delete it
												
				//if there is a possibility that this sequence should be reversed
				if (candidateSeq->getNumBases() < numBasesNeeded) {
					numFlipped_1++;
					string wasBetter =  "";
					//if the user wants you to try the reverse
					if (params->flip) {
				
						//get reverse compliment
						copy = new Sequence(candidateSeq->getName(), originalUnaligned);
						copy->reverseComplement();
                        
                        if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: flipping "  + candidateSeq->getName() + " \n"); }
						
						//rerun alignment
						Sequence temp2 = params->templateDB->findClosestSequence(copy, searchScore);
						Sequence* templateSeq2 = new Sequence(temp2.getName(), temp2.getAligned());
                        
                        if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: closest template "  + temp2.getName() + " \n"); }
						
						nast2 = new Nast(alignment, copy, templateSeq2);
                        
                        if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: completed Nast2 "  + candidateSeq->getName() + " flipped numBases = " + toString(copy->getNumBases()) + " old numbases = " + toString(candidateSeq->getNumBases()) +" \n"); }
			
						//check if any better
						if (copy->getNumBases() > candidateSeq->getNumBases()) {
							candidateSeq->setAligned(copy->getAligned());  //use reverse compliments alignment since its better
                            delete templateSeq;
							templateSeq = templateSeq2;
							delete nast;
							nast = nast2;
							needToDeleteCopy = true;
							wasBetter = "\treverse complement produced a better alignment, so mothur used the reverse complement.";
                            numFlipped_0++;
						}else{  
							wasBetter = "\treverse complement did NOT produce a better alignment so it was not used, please check sequence.";
							delete nast2;
                            delete templateSeq2;
							delete copy;	
						}
                        if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: done.\n"); }
					}
					
					//create accnos file with names
					accnosFile << candidateSeq->getName() << wasBetter << endl;
				}
				
				report.setTemplate(templateSeq);
				report.setSearchParameters(params->search, searchScore);
				report.setAlignmentParameters(params->alignMethod, alignment);
				report.setNastParameters(*nast);
	
				alignmentFile << '>' << candidateSeq->getName() << '\n' << candidateSeq->getAligned() << endl;
				
				report.print();
				delete nast;
                delete templateSeq;
				if (needToDeleteCopy) {   delete copy;   }
                
				count++;
			}
			delete candidateSeq;
			
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				unsigned long long pos = inFASTA.tellg();
				if ((pos == -1) || (pos >= params->filePos.end)) { break; }
			#else
				if (count == params->filePos.end) { break; }
			#endif
			
			//report progress
			if((count) % 100 == 0){	params->m->mothurOutJustToScreen(toString(count) + "\n"); 		}
			
		}
		//report progress
		if((count) % 100 != 0){	params->m->mothurOutJustToScreen(toString(count) + "\n"); 		}
        
        params->numSeqs += count;
        params->flippedResults[0] += numFlipped_0;
        params->flippedResults[1] += numFlipped_1;
        
		delete alignment;
		alignmentFile.close();
		inFASTA.close();
		accnosFile.close();
		
	}
	catch(exception& e) {
		params->m->errorOut(e, "AlignCommand", "driver");
		exit(1);
	}
}
/**************************************************************************************************/
//void alignDriver(linePair* filePos, string alignFName, string reportFName, string accnosFName, string filename, vector<long long>& numFlipped,MothurOut* m, string align, float match, float misMatch, float gapOpen, float gapExtend, float threshold, bool flip, AlignmentDB* templateDB, string search, long long& count) {
long long AlignCommand::createProcesses(string alignFileName, string reportFileName, string accnosFName, string filename, vector<long long>& numFlipped) {
	try {
        //create array of worker threads
        vector<thread*> workerThreads;
        vector<alignStruct*> data;
        mutex mutex;
        
        long long num = 0;
        for (int i = 0; i < numFlipped.size(); i++) { numFlipped[i] = 0; }
        
        time_t start, end;
        time(&start);
        //Lauch worker threads
        for (int i = 0; i < processors-1; ++i) {
            //alignStruct (linePair fP, string aFName, string reFName, string ac, string fname, MothurOut* mo, string al, float ma, float misMa, float gOpen, float gExtend, float thr, bool fl, AlignmentDB* tB, string se)
            alignStruct* dataBundle = new alignStruct(*lines[i+1], (alignFileName + toString(i+1) + ".temp"),
                                                        reportFileName + toString(i+1) + ".temp",
                                                        accnosFName + toString(i+1) + ".temp", filename,
                                                        m, align, match, misMatch, gapOpen, gapExtend, threshold, flip, templateDB, search);
            data.push_back(dataBundle);

            workerThreads.push_back(new thread(alignDriver, dataBundle));
         }
        
        alignStruct* dataBundle = new alignStruct(*lines[0], (alignFileName), reportFileName, accnosFName, filename,
                                                  m, align, match, misMatch, gapOpen, gapExtend, threshold, flip, templateDB, search);
        alignDriver(dataBundle);
        numFlipped[0] = dataBundle->flippedResults[0];
        numFlipped[1] = dataBundle->flippedResults[1];
        num = dataBundle->numSeqs;
        
        for (int i = 0; i < processors-1; ++i) {
            workerThreads[i]->join();
            num += data[i]->numSeqs;
            numFlipped[0] += data[i]->flippedResults[0];
            numFlipped[1] += data[i]->flippedResults[1];
            delete data[i];
            delete workerThreads[i];
        }
        
        
        time(&end);
        m->mothurOut("It took " + toString(difftime(end, start)) + " secs to align " + toString(num) + " sequences.\n\n");
        
        vector<string> nonBlankAccnosFiles;
        if (!(m->isBlank(accnosFName))) { nonBlankAccnosFiles.push_back(accnosFName); }
        else { m->mothurRemove(accnosFName); } //remove so other files can be renamed to it
        
        for (int i = 0; i < processors-1; i++) {
            m->appendFiles((alignFileName + toString(i+1) + ".temp"), alignFileName);
            m->mothurRemove((alignFileName + toString(i+1) + ".temp"));
            
            appendReportFiles((reportFileName + toString(i+1) + ".temp"), reportFileName);
            m->mothurRemove((reportFileName + toString(i+1) + ".temp"));
            
            if (!(m->isBlank(accnosFName + toString(i+1) + ".temp"))) {
                nonBlankAccnosFiles.push_back(accnosFName + toString(i+1) + ".temp");
            }else { m->mothurRemove((accnosFName + toString(i+1) + ".temp"));  }
            
        }
        
        //append accnos files
        if (nonBlankAccnosFiles.size() != 0) {
            rename(nonBlankAccnosFiles[0].c_str(), accnosFName.c_str());
            
            for (int h=1; h < nonBlankAccnosFiles.size(); h++) {
                m->appendFiles(nonBlankAccnosFiles[h], accnosFName);
                m->mothurRemove(nonBlankAccnosFiles[h]);
            }
        }else { //recreate the accnosfile if needed
            ofstream out;
            m->openOutputFile(accnosFName, out);
            out.close();
        }
        
        return num;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "createProcesses");
		exit(1);
	}
}
//**********************************************************************************************************************

void AlignCommand::appendReportFiles(string temp, string filename) {
	try{
		
		ofstream output;
		ifstream input;
		m->openOutputFileAppend(filename, output);
		m->openInputFile(temp, input);

		while (!input.eof())	{	char c = input.get(); if (c == 10 || c == 13){	break;	}	} // get header line
				
        char buffer[4096];        
        while (!input.eof()) {
            input.read(buffer, 4096);
            output.write(buffer, input.gcount());
        }
		
		input.close();
		output.close();
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "appendReportFiles");
		exit(1);
	}
}
//**********************************************************************************************************************
