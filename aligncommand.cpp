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
#include "sequence.hpp"

#include "gotohoverlap.hpp"
#include "needlemanoverlap.hpp"
#include "blastalign.hpp"
#include "noalign.hpp"

#include "nast.hpp"
#include "nastreport.hpp"


//**********************************************************************************************************************
vector<string> AlignCommand::getValidParameters(){	
	try {
		string AlignArray[] =  {"template","candidate","search","ksize","align","match","mismatch","gapopen","gapextend", "processors","flip","threshold","outputdir","inputdir"};
		vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> AlignCommand::getRequiredParameters(){	
	try {
		string AlignArray[] =  {"template","candidate"};
		vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> AlignCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************
AlignCommand::AlignCommand(){	
	try {
		//initialize outputTypes
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
		abort = false;
	
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			
			//valid paramters for this command
			string AlignArray[] =  {"template","candidate","search","ksize","align","match","mismatch","gapopen","gapextend", "processors","flip","threshold","outputdir","inputdir"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
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

				it = parameters.find("template");

				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["template"] = inputDir + it->second;		}
				}
			}

			//check for required parameters
			templateFileName = validParameter.validFile(parameters, "template", true);
			
			if (templateFileName == "not found") { 
				m->mothurOut("template is a required parameter for the align.seqs command."); 
				m->mothurOutEndLine();
				abort = true; 
			}else if (templateFileName == "not open") { abort = true; }	
			
			candidateFileName = validParameter.validFile(parameters, "candidate", false);
			if (candidateFileName == "not found") { m->mothurOut("candidate is a required parameter for the align.seqs command."); m->mothurOutEndLine(); abort = true;  }
			else { 
				m->splitAtDash(candidateFileName, candidateFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < candidateFileNames.size(); i++) {
					//candidateFileNames[i] = m->getFullPathName(candidateFileNames[i]);
					
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
					
					//if you can't open it, try default location
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
					}
					
				}
				
				//make sure there is at least one valid file left
				if (candidateFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
		
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.validFile(parameters, "ksize", false);		if (temp == "not found"){	temp = "8";				}
			convert(temp, kmerSize); 
			
			temp = validParameter.validFile(parameters, "match", false);		if (temp == "not found"){	temp = "1.0";			}
			convert(temp, match);  
			
			temp = validParameter.validFile(parameters, "mismatch", false);		if (temp == "not found"){	temp = "-1.0";			}
			convert(temp, misMatch);  
			
			temp = validParameter.validFile(parameters, "gapopen", false);		if (temp == "not found"){	temp = "-2.0";			}
			convert(temp, gapOpen);  
			
			temp = validParameter.validFile(parameters, "gapextend", false);	if (temp == "not found"){	temp = "-1.0";			}
			convert(temp, gapExtend); 
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = "1";				}
			convert(temp, processors); 
			
			temp = validParameter.validFile(parameters, "flip", false);			if (temp == "not found"){	temp = "f";				}
			flip = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "threshold", false);	if (temp == "not found"){	temp = "0.50";			}
			convert(temp, threshold); 
			
			search = validParameter.validFile(parameters, "search", false);		if (search == "not found"){	search = "kmer";		}
			
			align = validParameter.validFile(parameters, "align", false);		if (align == "not found"){	align = "needleman";	}
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
		delete alignment;
	}
}

//**********************************************************************************************************************

void AlignCommand::help(){
	try {
		m->mothurOut("The align.seqs command reads a file containing sequences and creates an alignment file and a report file.\n");
		m->mothurOut("The align.seqs command parameters are template, candidate, search, ksize, align, match, mismatch, gapopen, gapextend and processors.\n");
		m->mothurOut("The template and candidate parameters are required. You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amzon.fasta \n");
		m->mothurOut("The search parameter allows you to specify the method to find most similar template.  Your options are: suffix, kmer and blast. The default is kmer.\n");
		m->mothurOut("The align parameter allows you to specify the alignment method to use.  Your options are: gotoh, needleman, blast and noalign. The default is needleman.\n");
		m->mothurOut("The ksize parameter allows you to specify the kmer size for finding most similar template to candidate.  The default is 8.\n");
		m->mothurOut("The match parameter allows you to specify the bonus for having the same base. The default is 1.0.\n");
		m->mothurOut("The mistmatch parameter allows you to specify the penalty for having different bases.  The default is -1.0.\n");
		m->mothurOut("The gapopen parameter allows you to specify the penalty for opening a gap in an alignment. The default is -2.0.\n");
		m->mothurOut("The gapextend parameter allows you to specify the penalty for extending a gap in an alignment.  The default is -1.0.\n");
		m->mothurOut("The flip parameter is used to specify whether or not you want mothur to try the reverse complement if a sequence falls below the threshold.  The default is false.\n");
		m->mothurOut("The threshold is used to specify a cutoff at which an alignment is deemed 'bad' and the reverse complement may be tried. The default threshold is 0.50, meaning 50% of the bases are removed in the alignment.\n");
		m->mothurOut("If the flip parameter is set to true the reverse complement of the sequence is aligned and the better alignment is reported.\n");
		m->mothurOut("The default for the threshold parameter is 0.50, meaning at least 50% of the bases must remain or the sequence is reported as potentially reversed.\n");
		m->mothurOut("The align.seqs command should be in the following format: \n");
		m->mothurOut("align.seqs(template=yourTemplateFile, candidate=yourCandidateFile, align=yourAlignmentMethod, search=yourSearchmethod, ksize=yourKmerSize, match=yourMatchBonus, mismatch=yourMismatchpenalty, gapopen=yourGapopenPenalty, gapextend=yourGapExtendPenalty) \n");
		m->mothurOut("Example align.seqs(candidate=candidate.fasta, template=core.filtered, align=kmer, search=gotoh, ksize=8, match=2.0, mismatch=3.0, gapopen=-2.0, gapextend=-1.0)\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. candidate), '=' and parameters (i.e.yourFastaFile).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "help");
		exit(1);
	}
}


//**********************************************************************************************************************

int AlignCommand::execute(){
	try {
		if (abort == true) {	return 0;	}

		templateDB = new AlignmentDB(templateFileName, search, kmerSize, gapOpen, gapExtend, match, misMatch);
		int longestBase = templateDB->getLongestBase();
		
		if(align == "gotoh")			{	alignment = new GotohOverlap(gapOpen, gapExtend, match, misMatch, longestBase);			}
		else if(align == "needleman")	{	alignment = new NeedlemanOverlap(gapOpen, match, misMatch, longestBase);				}
		else if(align == "blast")		{	alignment = new BlastAlignment(gapOpen, gapExtend, match, misMatch);		}
		else if(align == "noalign")		{	alignment = new NoAlign();													}
		else {
			m->mothurOut(align + " is not a valid alignment option. I will run the command using needleman.");
			m->mothurOutEndLine();
			alignment = new NeedlemanOverlap(gapOpen, match, misMatch, longestBase);
		}
		
		for (int s = 0; s < candidateFileNames.size(); s++) {
			if (m->control_pressed) { outputTypes.clear(); return 0; }
			
			m->mothurOut("Aligning sequences from " + candidateFileNames[s] + " ..." ); m->mothurOutEndLine();
			
			if (outputDir == "") {  outputDir += m->hasPath(candidateFileNames[s]); }
			string alignFileName = outputDir + m->getRootName(m->getSimpleName(candidateFileNames[s])) + "align";
			string reportFileName = outputDir + m->getRootName(m->getSimpleName(candidateFileNames[s])) + "align.report";
			string accnosFileName = outputDir + m->getRootName(m->getSimpleName(candidateFileNames[s])) + "flip.accnos";
			bool hasAccnos = true;
			
			int numFastaSeqs = 0;
			for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
			int start = time(NULL);
		
#ifdef USE_MPI	
				int pid, end, numSeqsPerProcessor; 
				int tag = 2001;
				vector<unsigned long int> MPIPos;
				MPIWroteAccnos = false;
			
				MPI_Status status; 
				MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
				MPI_Comm_size(MPI_COMM_WORLD, &processors); 

				MPI_File inMPI;
				MPI_File outMPIAlign;
				MPI_File outMPIReport;
				MPI_File outMPIAccnos;
				
				int outMode=MPI_MODE_CREATE|MPI_MODE_WRONLY; 
				int inMode=MPI_MODE_RDONLY; 
				
				char outAlignFilename[1024];
				strcpy(outAlignFilename, alignFileName.c_str());
				
				char outReportFilename[1024];
				strcpy(outReportFilename, reportFileName.c_str());
				
				char outAccnosFilename[1024];
				strcpy(outAccnosFilename, accnosFileName.c_str());
				
				char inFileName[1024];
				strcpy(inFileName, candidateFileNames[s].c_str());
				
				MPI_File_open(MPI_COMM_WORLD, inFileName, inMode, MPI_INFO_NULL, &inMPI);  //comm, filename, mode, info, filepointer
				MPI_File_open(MPI_COMM_WORLD, outAlignFilename, outMode, MPI_INFO_NULL, &outMPIAlign);
				MPI_File_open(MPI_COMM_WORLD, outReportFilename, outMode, MPI_INFO_NULL, &outMPIReport);
				MPI_File_open(MPI_COMM_WORLD, outAccnosFilename, outMode, MPI_INFO_NULL, &outMPIAccnos);
				
				if (m->control_pressed) { MPI_File_close(&inMPI);  MPI_File_close(&outMPIAlign);  MPI_File_close(&outMPIReport);  MPI_File_close(&outMPIAccnos); outputTypes.clear(); return 0; }
				
				if (pid == 0) { //you are the root process 
					
					MPIPos = m->setFilePosFasta(candidateFileNames[s], numFastaSeqs); //fills MPIPos, returns numSeqs
					
					//send file positions to all processes
					for(int i = 1; i < processors; i++) { 
						MPI_Send(&numFastaSeqs, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
						MPI_Send(&MPIPos[0], (numFastaSeqs+1), MPI_LONG, i, tag, MPI_COMM_WORLD);
					}
					
					//figure out how many sequences you have to align
					numSeqsPerProcessor = numFastaSeqs / processors;
					int startIndex =  pid * numSeqsPerProcessor;
					if(pid == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - pid * numSeqsPerProcessor; 	}
					
					//align your part
					driverMPI(startIndex, numSeqsPerProcessor, inMPI, outMPIAlign, outMPIReport, outMPIAccnos, MPIPos);
					
					if (m->control_pressed) { MPI_File_close(&inMPI);  MPI_File_close(&outMPIAlign);  MPI_File_close(&outMPIReport);  MPI_File_close(&outMPIAccnos); outputTypes.clear(); return 0; }

					for (int i = 1; i < processors; i++) {
						bool tempResult;
						MPI_Recv(&tempResult, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
						if (tempResult != 0) { MPIWroteAccnos = true; }
					}
				}else{ //you are a child process
					MPI_Recv(&numFastaSeqs, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
					MPIPos.resize(numFastaSeqs+1);
					MPI_Recv(&MPIPos[0], (numFastaSeqs+1), MPI_LONG, 0, tag, MPI_COMM_WORLD, &status);

					
					//figure out how many sequences you have to align
					numSeqsPerProcessor = numFastaSeqs / processors;
					int startIndex =  pid * numSeqsPerProcessor;
					if(pid == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - pid * numSeqsPerProcessor; 	}
					
					
					//align your part
					driverMPI(startIndex, numSeqsPerProcessor, inMPI, outMPIAlign, outMPIReport, outMPIAccnos, MPIPos);
					
					if (m->control_pressed) { MPI_File_close(&inMPI);  MPI_File_close(&outMPIAlign);  MPI_File_close(&outMPIReport);  MPI_File_close(&outMPIAccnos); outputTypes.clear(); return 0; }

					MPI_Send(&MPIWroteAccnos, 1, MPI_INT, 0, tag, MPI_COMM_WORLD); 
				}
				
				//close files 
				MPI_File_close(&inMPI);
				MPI_File_close(&outMPIAlign);
				MPI_File_close(&outMPIReport);
				MPI_File_close(&outMPIAccnos);
				
				//delete accnos file if blank
				if (pid == 0) {
					//delete accnos file if its blank else report to user
					if (MPIWroteAccnos) { 
						m->mothurOut("Some of you sequences generated alignments that eliminated too many bases, a list is provided in " + accnosFileName + ".");
						if (!flip) {
							m->mothurOut(" If you set the flip parameter to true mothur will try aligning the reverse compliment as well."); 
						}else{  m->mothurOut(" If the reverse compliment proved to be better it was reported.");  }
						m->mothurOutEndLine();
					}else { 
						//MPI_Info info;
						//MPI_File_delete(outAccnosFilename, info);
						hasAccnos = false;	
						remove(accnosFileName.c_str()); 
					}
				}
				
#else

		vector<unsigned long int> positions = m->divideFile(candidateFileNames[s], processors);
		for (int i = 0; i < (positions.size()-1); i++) {
			lines.push_back(new linePair(positions[i], positions[(i+1)]));
		}	
	#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			if(processors == 1){
				numFastaSeqs = driver(lines[0], alignFileName, reportFileName, accnosFileName, candidateFileNames[s]);
				
				if (m->control_pressed) { remove(accnosFileName.c_str()); remove(alignFileName.c_str()); remove(reportFileName.c_str()); outputTypes.clear(); return 0; }
				
				//delete accnos file if its blank else report to user
				if (m->isBlank(accnosFileName)) {  remove(accnosFileName.c_str());  hasAccnos = false; }
				else { 
					m->mothurOut("Some of you sequences generated alignments that eliminated too many bases, a list is provided in " + accnosFileName + ".");
					if (!flip) {
						m->mothurOut(" If you set the flip parameter to true mothur will try aligning the reverse compliment as well."); 
					}else{  m->mothurOut(" If the reverse compliment proved to be better it was reported.");  }
					m->mothurOutEndLine();
				}
			}else{
				processIDS.resize(0);
				
				numFastaSeqs = createProcesses(alignFileName, reportFileName, accnosFileName, candidateFileNames[s]); 
				
				rename((alignFileName + toString(processIDS[0]) + ".temp").c_str(), alignFileName.c_str());
				rename((reportFileName + toString(processIDS[0]) + ".temp").c_str(), reportFileName.c_str());
				
				//append alignment and report files
				for(int i=1;i<processors;i++){
					appendAlignFiles((alignFileName + toString(processIDS[i]) + ".temp"), alignFileName);
					remove((alignFileName + toString(processIDS[i]) + ".temp").c_str());
					
					appendReportFiles((reportFileName + toString(processIDS[i]) + ".temp"), reportFileName);
					remove((reportFileName + toString(processIDS[i]) + ".temp").c_str());
				}
				
				vector<string> nonBlankAccnosFiles;
				//delete blank accnos files generated with multiple processes
				for(int i=0;i<processors;i++){  
					if (!(m->isBlank(accnosFileName + toString(processIDS[i]) + ".temp"))) {
						nonBlankAccnosFiles.push_back(accnosFileName + toString(processIDS[i]) + ".temp");
					}else { remove((accnosFileName + toString(processIDS[i]) + ".temp").c_str());  }
				}
				
				//append accnos files
				if (nonBlankAccnosFiles.size() != 0) { 
					rename(nonBlankAccnosFiles[0].c_str(), accnosFileName.c_str());
					
					for (int h=1; h < nonBlankAccnosFiles.size(); h++) {
						appendAlignFiles(nonBlankAccnosFiles[h], accnosFileName);
						remove(nonBlankAccnosFiles[h].c_str());
					}
					m->mothurOut("Some of you sequences generated alignments that eliminated too many bases, a list is provided in " + accnosFileName + ".");
					if (!flip) {
						m->mothurOut(" If you set the flip parameter to true mothur will try aligning the reverse compliment as well."); 
					}else{  m->mothurOut(" If the reverse compliment proved to be better it was reported.");  }
					m->mothurOutEndLine();
				}else{ hasAccnos = false;  }
				
				if (m->control_pressed) { remove(accnosFileName.c_str()); remove(alignFileName.c_str()); remove(reportFileName.c_str()); outputTypes.clear(); return 0; }
			}
	#else
			numFastaSeqs = driver(lines[0], alignFileName, reportFileName, accnosFileName, candidateFileNames[s]);
			
			if (m->control_pressed) { remove(accnosFileName.c_str()); remove(alignFileName.c_str()); remove(reportFileName.c_str()); outputTypes.clear();  return 0; }
			
			//delete accnos file if its blank else report to user
			if (m->isBlank(accnosFileName)) {  remove(accnosFileName.c_str());  hasAccnos = false; }
			else { 
				m->mothurOut("Some of you sequences generated alignments that eliminated too many bases, a list is provided in " + accnosFileName + ".");
				if (!flip) {
					m->mothurOut(" If you set the flip parameter to true mothur will try aligning the reverse compliment as well."); 
				}else{  m->mothurOut(" If the reverse compliment proved to be better it was reported.");  }
				m->mothurOutEndLine();
			}
			
	#endif

#endif		


		#ifdef USE_MPI
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should output to screen
		#endif

			outputNames.push_back(alignFileName); outputTypes["fasta"].push_back(alignFileName);
			outputNames.push_back(reportFileName); outputTypes["alignreport"].push_back(reportFileName);
			if (hasAccnos)	{	outputNames.push_back(accnosFileName);	outputTypes["accnos"].push_back(accnosFileName);  }
			
		#ifdef USE_MPI
			}
		#endif

			m->mothurOut("It took " + toString(time(NULL) - start) + " secs to align " + toString(numFastaSeqs) + " sequences.");
			m->mothurOutEndLine();
			m->mothurOutEndLine();
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

//**********************************************************************************************************************

int AlignCommand::driver(linePair* filePos, string alignFName, string reportFName, string accnosFName, string filename){
	try {
		ofstream alignmentFile;
		m->openOutputFile(alignFName, alignmentFile);
		
		ofstream accnosFile;
		m->openOutputFile(accnosFName, accnosFile);
		
		NastReport report(reportFName);
		
		ifstream inFASTA;
		m->openInputFile(filename, inFASTA);

		inFASTA.seekg(filePos->start);

		bool done = false;
		int count = 0;
	
		while (!done) {
			
			if (m->control_pressed) {  return 0; }
			
			Sequence* candidateSeq = new Sequence(inFASTA);  m->gobble(inFASTA);
			
			int origNumBases = candidateSeq->getNumBases();
			string originalUnaligned = candidateSeq->getUnaligned();
			int numBasesNeeded = origNumBases * threshold;
	
			if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
				if (candidateSeq->getUnaligned().length() > alignment->getnRows()) {
					alignment->resize(candidateSeq->getUnaligned().length()+1);
				}
								
				Sequence temp = templateDB->findClosestSequence(candidateSeq);
				Sequence* templateSeq = &temp;
			
				float searchScore = templateDB->getSearchScore();
								
				Nast* nast = new Nast(alignment, candidateSeq, templateSeq);
		
				Sequence* copy;
				
				Nast* nast2;
				bool needToDeleteCopy = false;  //this is needed in case you have you enter the ifs below
												//since nast does not make a copy of hte sequence passed, and it is used by the reporter below
												//you can't delete the copy sequence til after you report, but you may choose not to create it in the first place
												//so this bool tells you if you need to delete it
												
				//if there is a possibility that this sequence should be reversed
				if (candidateSeq->getNumBases() < numBasesNeeded) {
					
					string wasBetter =  "";
					//if the user wants you to try the reverse
					if (flip) {
				
						//get reverse compliment
						copy = new Sequence(candidateSeq->getName(), originalUnaligned);
						copy->reverseComplement();
						
						//rerun alignment
						Sequence temp2 = templateDB->findClosestSequence(copy);
						Sequence* templateSeq2 = &temp2;
						
						searchScore = templateDB->getSearchScore();
						
						nast2 = new Nast(alignment, copy, templateSeq2);
			
						//check if any better
						if (copy->getNumBases() > candidateSeq->getNumBases()) {
							candidateSeq->setAligned(copy->getAligned());  //use reverse compliments alignment since its better
							templateSeq = templateSeq2; 
							delete nast;
							nast = nast2;
							needToDeleteCopy = true;
							wasBetter = "\treverse complement produced a better alignment, so mothur used the reverse complement.";
						}else{  
							wasBetter = "\treverse complement did NOT produce a better alignment so it was not used, please check sequence.";
							delete nast2;
							delete copy;	
						}
					}
					
					//create accnos file with names
					accnosFile << candidateSeq->getName() << wasBetter << endl;
				}
				
				report.setCandidate(candidateSeq);
				report.setTemplate(templateSeq);
				report.setSearchParameters(search, searchScore);
				report.setAlignmentParameters(align, alignment);
				report.setNastParameters(*nast);
	
				alignmentFile << '>' << candidateSeq->getName() << '\n' << candidateSeq->getAligned() << endl;
				
				report.print();
				delete nast;
				if (needToDeleteCopy) {   delete copy;   }
				
				count++;
			}
			delete candidateSeq;
			
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				unsigned long int pos = inFASTA.tellg();
				if ((pos == -1) || (pos >= filePos->end)) { break; }
			#else
				if (inFASTA.eof()) { break; }
			#endif
			
			//report progress
			if((count) % 100 == 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}
			
		}
		//report progress
		if((count) % 100 != 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}
		
		alignmentFile.close();
		inFASTA.close();
		accnosFile.close();
		
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
#ifdef USE_MPI
int AlignCommand::driverMPI(int start, int num, MPI_File& inMPI, MPI_File& alignFile, MPI_File& reportFile, MPI_File& accnosFile, vector<unsigned long int>& MPIPos){
	try {
		string outputString = "";
		MPI_Status statusReport; 
		MPI_Status statusAlign; 
		MPI_Status statusAccnos; 
		MPI_Status status; 
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
	
		NastReport report;
		
		if (pid == 0) {
			outputString = report.getHeaders();
			int length = outputString.length();
            
			char* buf = new char[length];
			memcpy(buf, outputString.c_str(), length);
		
			MPI_File_write_shared(reportFile, buf, length, MPI_CHAR, &statusReport);

            delete buf;
		}
		
		for(int i=0;i<num;i++){
		
			if (m->control_pressed) {  return 0; }

			//read next sequence
			int length = MPIPos[start+i+1] - MPIPos[start+i];

			char* buf4 = new char[length];
			//memcpy(buf4, outputString.c_str(), length);

			MPI_File_read_at(inMPI, MPIPos[start+i], buf4, length, MPI_CHAR, &status);
			
			string tempBuf = buf4;

			delete buf4;

			if (tempBuf.length() > length) { tempBuf = tempBuf.substr(0, length);  }
	
			istringstream iss (tempBuf,istringstream::in);

			Sequence* candidateSeq = new Sequence(iss);  
	
			int origNumBases = candidateSeq->getNumBases();
			string originalUnaligned = candidateSeq->getUnaligned();
			int numBasesNeeded = origNumBases * threshold;
	
			if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
				if (candidateSeq->getUnaligned().length() > alignment->getnRows()) {
					alignment->resize(candidateSeq->getUnaligned().length()+1);
				}
								
				Sequence temp = templateDB->findClosestSequence(candidateSeq);
				Sequence* templateSeq = &temp;
				
				float searchScore = templateDB->getSearchScore();
								
				Nast* nast = new Nast(alignment, candidateSeq, templateSeq);
				Sequence* copy;
				
				Nast* nast2;
				bool needToDeleteCopy = false;  //this is needed in case you have you enter the ifs below
												//since nast does not make a copy of hte sequence passed, and it is used by the reporter below
												//you can't delete the copy sequence til after you report, but you may choose not to create it in the first place
												//so this bool tells you if you need to delete it
												
				//if there is a possibility that this sequence should be reversed
				if (candidateSeq->getNumBases() < numBasesNeeded) {
					
					string wasBetter = "";
					//if the user wants you to try the reverse
					if (flip) {
						//get reverse compliment
						copy = new Sequence(candidateSeq->getName(), originalUnaligned);
						copy->reverseComplement();
						
						//rerun alignment
						Sequence temp2 = templateDB->findClosestSequence(copy);
						Sequence* templateSeq2 = &temp2;
						
						searchScore = templateDB->getSearchScore();
						
						nast2 = new Nast(alignment, copy, templateSeq2);
			
						//check if any better
						if (copy->getNumBases() > candidateSeq->getNumBases()) {
							candidateSeq->setAligned(copy->getAligned());  //use reverse compliments alignment since its better
							templateSeq = templateSeq2; 
							delete nast;
							nast = nast2;
							needToDeleteCopy = true;
							wasBetter = "\treverse complement produced a better alignment, so mothur used the reverse complement.";
						}else{  
							wasBetter = "\treverse complement did NOT produce a better alignment, please check sequence.";
							delete nast2;
							delete copy;	
						}
					}
					
					//create accnos file with names
					outputString = candidateSeq->getName() + wasBetter + "\n";
					
					//send results to parent
					int length = outputString.length();

					char* buf = new char[length];
					memcpy(buf, outputString.c_str(), length);
				
					MPI_File_write_shared(accnosFile, buf, length, MPI_CHAR, &statusAccnos);
					delete buf;
					MPIWroteAccnos = true;
				}
				
				report.setCandidate(candidateSeq);
				report.setTemplate(templateSeq);
				report.setSearchParameters(search, searchScore);
				report.setAlignmentParameters(align, alignment);
				report.setNastParameters(*nast);
	
				outputString =  ">" + candidateSeq->getName() + "\n" + candidateSeq->getAligned() + "\n";
				
				//send results to parent
				int length = outputString.length();
				char* buf2 = new char[length];
				memcpy(buf2, outputString.c_str(), length);
				
				MPI_File_write_shared(alignFile, buf2, length, MPI_CHAR, &statusAlign);
				
				delete buf2;

				outputString = report.getReport();
				
				//send results to parent
				length = outputString.length();
				char* buf3 = new char[length];
				memcpy(buf3, outputString.c_str(), length);
				
				MPI_File_write_shared(reportFile, buf3, length, MPI_CHAR, &statusReport);
				
				delete buf3;
				delete nast;
				if (needToDeleteCopy) {   delete copy;   }
			}
			delete candidateSeq;
			
			//report progress
			if((i+1) % 100 == 0){	cout << (toString(i+1)) << endl;		}
		}
		//report progress
		if((num) % 100 != 0){	cout << (toString(num)) << endl;		}
		
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "driverMPI");
		exit(1);
	}
}
#endif
/**************************************************************************************************/

int AlignCommand::createProcesses(string alignFileName, string reportFileName, string accnosFName, string filename) {
	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 0;
		int num = 0;
		//		processIDS.resize(0);
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driver(lines[process], alignFileName + toString(getpid()) + ".temp", reportFileName + toString(getpid()) + ".temp", accnosFName + toString(getpid()) + ".temp", filename);
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = alignFileName + toString(getpid()) + ".num.temp";
				m->openOutputFile(tempFile, out);
				out << num << endl;
				out.close();
				
				exit(0);
			}else { m->mothurOut("unable to spawn the necessary processes."); m->mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		for (int i = 0; i < processIDS.size(); i++) {
			ifstream in;
			string tempFile =  alignFileName + toString(processIDS[i]) + ".num.temp";
			m->openInputFile(tempFile, in);
			if (!in.eof()) { int tempNum = 0; in >> tempNum; num += tempNum; }
			in.close(); remove(tempFile.c_str());
		}
		
		return num;
#endif		
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "createProcesses");
		exit(1);
	}
}
/**************************************************************************************************/

void AlignCommand::appendAlignFiles(string temp, string filename) {
	try{
		
		ofstream output;
		ifstream input;
		m->openOutputFileAppend(filename, output);
		m->openInputFile(temp, input);
		
		while(char c = input.get()){
			if(input.eof())		{	break;			}
			else				{	output << c;	}
		}
		
		input.close();
		output.close();
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCommand", "appendAlignFiles");
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
				
		while(char c = input.get()){
			if(input.eof())		{	break;			}
			else				{	output << c;	}
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
