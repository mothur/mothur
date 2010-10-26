/*
 *  pairwiseseqscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 10/20/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "pairwiseseqscommand.h"
#include "sequence.hpp"

#include "gotohoverlap.hpp"
#include "needlemanoverlap.hpp"
#include "blastalign.hpp"
#include "noalign.hpp"
#include "nast.hpp"

#include "ignoregaps.h"
#include "eachgapdist.h"
#include "eachgapignore.h"
#include "onegapdist.h"
#include "onegapignore.h"


//**********************************************************************************************************************
vector<string> PairwiseSeqsCommand::getValidParameters(){	
	try {
		string AlignArray[] =  {"fasta","align","match","mismatch","gapopen","gapextend", "processors","calc","compress","cutoff","countends","output","outputdir","inputdir"};
		vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PairwiseSeqsCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> PairwiseSeqsCommand::getRequiredParameters(){	
	try {
		string AlignArray[] =  {"fasta"};
		vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PairwiseSeqsCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> PairwiseSeqsCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PairwiseSeqsCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************
PairwiseSeqsCommand::PairwiseSeqsCommand(){	
	try {
		abort = true;
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["phylip"] = tempOutNames;
		outputTypes["column"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "PairwiseSeqsCommand", "PairwiseSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
PairwiseSeqsCommand::PairwiseSeqsCommand(string option)  {
	try {
		abort = false;
	
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			
			//valid paramters for this command
			string AlignArray[] =  {"fasta","align","match","mismatch","gapopen","gapextend", "processors","cutoff","compress","calc","countends","output","outputdir","inputdir"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
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
			outputTypes["phylip"] = tempOutNames;
			outputTypes["column"] = tempOutNames;
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			

			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			
			if (inputDir == "not found"){	inputDir = "";		}

			fastaFileName = validParameter.validFile(parameters, "fasta", false);
			if (fastaFileName == "not found") { m->mothurOut("fasta is a required parameter for the pairwise.seqs command."); m->mothurOutEndLine(); abort = true;  }
			else { 
				m->splitAtDash(fastaFileName, fastaFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastaFileNames.size(); i++) {
					
					if (inputDir != "") {
						string path = m->hasPath(fastaFileNames[i]);
						//if the user has not given a path then, add inputdir. else leave path alone.
						if (path == "") {	fastaFileNames[i] = inputDir + fastaFileNames[i];		}
					}
	
					int ableToOpen;
					ifstream in;

					ableToOpen = m->openInputFile(fastaFileNames[i], in, "noerror");
				
					//if you can't open it, try default location
					if (ableToOpen == 1) {
						if (m->getDefaultPath() != "") { //default path is set
							string tryPath = m->getDefaultPath() + m->getSimpleName(fastaFileNames[i]);
							m->mothurOut("Unable to open " + fastaFileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
							ifstream in2;
							ableToOpen = m->openInputFile(tryPath, in2, "noerror");
							in2.close();
							fastaFileNames[i] = tryPath;
						}
					}
					
					//if you can't open it, try output location
					if (ableToOpen == 1) {
						if (m->getOutputDir() != "") { //default path is set
							string tryPath = m->getOutputDir() + m->getSimpleName(fastaFileNames[i]);
							m->mothurOut("Unable to open " + fastaFileNames[i] + ". Trying output directory " + tryPath); m->mothurOutEndLine();
							ifstream in2;
							ableToOpen = m->openInputFile(tryPath, in2, "noerror");
							in2.close();
							fastaFileNames[i] = tryPath;
						}
					}
					
					in.close();					

					if (ableToOpen == 1) { 
						m->mothurOut("Unable to open " + fastaFileNames[i] + ". It will be disregarded."); m->mothurOutEndLine(); 
						//erase from file list
						fastaFileNames.erase(fastaFileNames.begin()+i);
						i--;
					}
				}
				
				//make sure there is at least one valid file left
				if (fastaFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
		
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
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
			
			temp = validParameter.validFile(parameters, "cutoff", false);		if(temp == "not found"){	temp = "1.0"; }
			convert(temp, cutoff); 
			
			temp = validParameter.validFile(parameters, "countends", false);	if(temp == "not found"){	temp = "T";	}
			countends = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "compress", false);		if(temp == "not found"){  temp = "F"; }
			compress = m->isTrue(temp); 
			
			align = validParameter.validFile(parameters, "align", false);		if (align == "not found"){	align = "needleman";	}
			
			output = validParameter.validFile(parameters, "output", false);		if(output == "not found"){	output = "column"; }
			if ((output != "column") && (output != "lt") && (output != "square")) { m->mothurOut(output + " is not a valid output form. Options are column, lt and square. I will use column."); m->mothurOutEndLine(); output = "column"; }
			
			calc = validParameter.validFile(parameters, "calc", false);			
			if (calc == "not found") { calc = "onegap";  }
			else { 
				 if (calc == "default")  {  calc = "onegap";  }
			}
			m->splitAtDash(calc, Estimators);
			
			ValidCalculators validCalculator;
			if (countends) {
				for (int i=0; i<Estimators.size(); i++) {
					if (validCalculator.isValidCalculator("distance", Estimators[i]) == true) { 
						if (Estimators[i] == "nogaps")			{	distCalculator = new ignoreGaps();	}
						else if (Estimators[i] == "eachgap")	{	distCalculator = new eachGapDist();	}
						else if (Estimators[i] == "onegap")		{	distCalculator = new oneGapDist();	}
					}
				}
			}else {
				for (int i=0; i<Estimators.size(); i++) {
					if (validCalculator.isValidCalculator("distance", Estimators[i]) == true) { 
						if (Estimators[i] == "nogaps")		{	distCalculator = new ignoreGaps();					}
						else if (Estimators[i] == "eachgap"){	distCalculator = new eachGapIgnoreTermGapDist();	}
						else if (Estimators[i] == "onegap")	{	distCalculator = new oneGapIgnoreTermGapDist();		}
					}
				}
			}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "PairwiseSeqsCommand", "PairwiseSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
PairwiseSeqsCommand::~PairwiseSeqsCommand(){}
//**********************************************************************************************************************

void PairwiseSeqsCommand::help(){
	try {
		m->mothurOut("The pairwise.seqs command reads a fasta file and creates distance matrix.\n");
		m->mothurOut("The pairwise.seqs command parameters are fasta, align, match, mismatch, gapopen, gapextend, calc, output, cutoff and processors.\n");
		m->mothurOut("The fasta parameter is required. You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amzon.fasta \n");
		m->mothurOut("The align parameter allows you to specify the alignment method to use.  Your options are: gotoh, needleman, blast and noalign. The default is needleman.\n");
		m->mothurOut("The match parameter allows you to specify the bonus for having the same base. The default is 1.0.\n");
		m->mothurOut("The mistmatch parameter allows you to specify the penalty for having different bases.  The default is -1.0.\n");
		m->mothurOut("The gapopen parameter allows you to specify the penalty for opening a gap in an alignment. The default is -2.0.\n");
		m->mothurOut("The gapextend parameter allows you to specify the penalty for extending a gap in an alignment.  The default is -1.0.\n");
		m->mothurOut("The calc parameter allows you to specify the method of calculating the distances.  Your options are: nogaps, onegap or eachgap. The default is onegap.\n");
		m->mothurOut("The countends parameter allows you to specify whether to include terminal gaps in distance.  Your options are: T or F. The default is T.\n");
		m->mothurOut("The cutoff parameter allows you to specify maximum distance to keep. The default is 1.0.\n");
		m->mothurOut("The output parameter allows you to specify format of your distance matrix. Options are column, lt, and square. The default is column.\n");
		m->mothurOut("The compress parameter allows you to indicate that you want the resulting distance file compressed.  The default is false.\n");
		m->mothurOut("The pairwise.seqs command should be in the following format: \n");
		m->mothurOut("pairwise.seqs(fasta=yourfastaFile, align=yourAlignmentMethod) \n");
		m->mothurOut("Example pairwise.seqs(fasta=candidate.fasta, align=blast)\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "PairwiseSeqsCommand", "help");
		exit(1);
	}
}


//**********************************************************************************************************************

int PairwiseSeqsCommand::execute(){
	try {
		if (abort == true) {	return 0;	}
		
		int longestBase = 2000; //will need to update this in driver if we find sequences with more bases.  hardcoded so we don't have the pre-read user fasta file.
		
		if(align == "gotoh")			{	alignment = new GotohOverlap(gapOpen, gapExtend, match, misMatch, longestBase);			}
		else if(align == "needleman")	{	alignment = new NeedlemanOverlap(gapOpen, match, misMatch, longestBase);				}
		else if(align == "blast")		{	alignment = new BlastAlignment(gapOpen, gapExtend, match, misMatch);		}
		else if(align == "noalign")		{	alignment = new NoAlign();													}
		else {
			m->mothurOut(align + " is not a valid alignment option. I will run the command using needleman.");
			m->mothurOutEndLine();
			alignment = new NeedlemanOverlap(gapOpen, match, misMatch, longestBase);
		}
		
		cutoff += 0.005;
		
		for (int s = 0; s < fastaFileNames.size(); s++) {
			if (m->control_pressed) { outputTypes.clear(); return 0; }
			
			m->mothurOut("Processing sequences from " + fastaFileNames[s] + " ..." ); m->mothurOutEndLine();
			
			if (outputDir == "") {  outputDir += m->hasPath(fastaFileNames[s]); }
			
			ifstream inFASTA;
			m->openInputFile(fastaFileNames[s], inFASTA);
			alignDB = SequenceDB(inFASTA); 
			inFASTA.close();
			
			int numSeqs = alignDB.getNumSeqs();
			int startTime = time(NULL);
			string outputFile = "";
				
			if (output == "lt") { //does the user want lower triangle phylip formatted file 
				outputFile = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + "phylip.dist";
				remove(outputFile.c_str()); outputTypes["phylip"].push_back(outputFile);
			}else if (output == "column") { //user wants column format
				outputFile = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + "dist";
				outputTypes["column"].push_back(outputFile);
				remove(outputFile.c_str());
			}else { //assume square
				outputFile = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + "square.dist";
				remove(outputFile.c_str());
				outputTypes["phylip"].push_back(outputFile);
			}
			
			#ifdef USE_MPI
		
			int pid, start, end; 
			int tag = 2001;
					
			MPI_Status status; 
			MPI_Comm_size(MPI_COMM_WORLD, &processors); //set processors to the number of mpi processes running
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
			
			//each process gets where it should start and stop in the file
			if (output != "square") {
				start = int (sqrt(float(pid)/float(processors)) * numSeqs);
				end = int (sqrt(float(pid+1)/float(processors)) * numSeqs);
			}else{
				start = int ((float(pid)/float(processors)) * numSeqs);
				end = int ((float(pid+1)/float(processors)) * numSeqs);
			}
			
			if (output == "column") {
				MPI_File outMPI;
				int amode=MPI_MODE_CREATE|MPI_MODE_WRONLY; 

			char filename[1024];
				strcpy(filename, outputFile.c_str());
				
				MPI_File_open(MPI_COMM_WORLD, filename, amode, MPI_INFO_NULL, &outMPI);

				if (pid == 0) { //you are the root process 
				
					//do your part
					string outputMyPart;
					
					driverMPI(start, end, outMPI, cutoff); 
					
					if (m->control_pressed) { outputTypes.clear(); MPI_File_close(&outMPI);  remove(outputFile.c_str()); delete distCalculator;  return 0; }
				
					//wait on chidren
					for(int i = 1; i < processors; i++) { 
						if (m->control_pressed) { outputTypes.clear();  MPI_File_close(&outMPI);   remove(outputFile.c_str()); delete distCalculator;  return 0; }
						
						char buf[5];
						MPI_Recv(buf, 5, MPI_CHAR, i, tag, MPI_COMM_WORLD, &status); 
					}
				}else { //you are a child process
					//do your part
					driverMPI(start, end, outMPI, cutoff); 
					
					if (m->control_pressed) { outputTypes.clear();  MPI_File_close(&outMPI);  remove(outputFile.c_str()); delete distCalculator;  return 0; }
				
					char buf[5];
					strcpy(buf, "done"); 
					//tell parent you are done.
					MPI_Send(buf, 5, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
				}
				
				MPI_File_close(&outMPI);
				
			}else { //lower triangle format
				if (pid == 0) { //you are the root process 
				
					//do your part
					string outputMyPart;
					unsigned long int mySize;
					
					if (output != "square"){ driverMPI(start, end, outputFile, mySize); }
					else { driverMPI(start, end, outputFile, mySize, output); }
		
					if (m->control_pressed) {  outputTypes.clear();   remove(outputFile.c_str()); delete distCalculator;  return 0; }
					
					int amode=MPI_MODE_APPEND|MPI_MODE_WRONLY|MPI_MODE_CREATE; //
					MPI_File outMPI;
					MPI_File inMPI;

					char filename[1024];
					strcpy(filename, outputFile.c_str());

					MPI_File_open(MPI_COMM_SELF, filename, amode, MPI_INFO_NULL, &outMPI);

					//wait on chidren
					for(int b = 1; b < processors; b++) { 
						unsigned long int fileSize;
						
						if (m->control_pressed) { outputTypes.clear();  MPI_File_close(&outMPI);  remove(outputFile.c_str());  delete distCalculator;  return 0; }
						
						MPI_Recv(&fileSize, 1, MPI_LONG, b, tag, MPI_COMM_WORLD, &status); 
						
						string outTemp = outputFile + toString(b) + ".temp";

						char* buf = new char[outTemp.length()];
						memcpy(buf, outTemp.c_str(), outTemp.length());
						
						MPI_File_open(MPI_COMM_SELF, buf, MPI_MODE_DELETE_ON_CLOSE|MPI_MODE_RDONLY, MPI_INFO_NULL, &inMPI);
						delete buf;

						int count = 0;
						while (count < fileSize) { 
							char buf2[1];
							MPI_File_read(inMPI, buf2, 1, MPI_CHAR, &status);
							MPI_File_write(outMPI, buf2, 1, MPI_CHAR, &status);
							count += 1;
						}
						
						MPI_File_close(&inMPI); //deleted on close
					}
					
					MPI_File_close(&outMPI);
				}else { //you are a child process
					//do your part
					unsigned long int size;
					if (output != "square"){ driverMPI(start, end, (outputFile + toString(pid) + ".temp"), size); }
					else { driverMPI(start, end, (outputFile + toString(pid) + ".temp"), size, output); }
					
					if (m->control_pressed) { delete distCalculator;  return 0; }
				
					//tell parent you are done.
					MPI_Send(&size, 1, MPI_LONG, 0, tag, MPI_COMM_WORLD);
				}
			}
			MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
	#else		
					
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			//if you don't need to fork anything
			if(processors == 1){
				if (output != "square") {  driver(0, numSeqs, outputFile, cutoff); }
				else { driver(0, numSeqs, outputFile, "square");  }
			}else{ //you have multiple processors
				
				for (int i = 0; i < processors; i++) {
					lines.push_back(new linePair());
					if (output != "square") {
						lines[i]->start = int (sqrt(float(i)/float(processors)) * numSeqs);
						lines[i]->end = int (sqrt(float(i+1)/float(processors)) * numSeqs);
					}else{
						lines[i]->start = int ((float(i)/float(processors)) * numSeqs);
						lines[i]->end = int ((float(i+1)/float(processors)) * numSeqs);
					}
				}

				createProcesses(outputFile); 
			
				map<int, int>::iterator it = processIDS.begin();
				rename((outputFile + toString(it->second) + ".temp").c_str(), outputFile.c_str());
				it++;
				
				//append and remove temp files
				for (; it != processIDS.end(); it++) {
					m->appendFiles((outputFile + toString(it->second) + ".temp"), outputFile);
					remove((outputFile + toString(it->second) + ".temp").c_str());
				}
			}
		#else
			//ifstream inFASTA;
			if (output != "square") {  driver(0, numSeqs, outputFile, cutoff); }
			else { driver(0, numSeqs, outputFile, "square");  }
		#endif
		
	#endif
			if (m->control_pressed) { outputTypes.clear();  delete distCalculator; remove(outputFile.c_str()); return 0; }
			
			#ifdef USE_MPI
				MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
						
				if (pid == 0) { //only one process should output to screen
			#endif
			
			ifstream fileHandle;
			fileHandle.open(outputFile.c_str());
			if(fileHandle) {
				m->gobble(fileHandle);
				if (fileHandle.eof()) { m->mothurOut(outputFile + " is blank. This can result if there are no distances below your cutoff.");  m->mothurOutEndLine(); }
			}
			
			if (compress) {
				m->mothurOut("Compressing..."); m->mothurOutEndLine();
				m->mothurOut("(Replacing " + outputFile + " with " + outputFile + ".gz)"); m->mothurOutEndLine();
				system(("gzip -v " + outputFile).c_str());
				outputNames.push_back(outputFile + ".gz");
			}else { outputNames.push_back(outputFile); }
			
			#ifdef USE_MPI
				}
			#endif
			
			m->mothurOut("It took " + toString(time(NULL) - startTime) + " to calculate the distances for " + toString(numSeqs) + " sequences."); m->mothurOutEndLine();
			
			if (m->control_pressed) { outputTypes.clear();  delete distCalculator; remove(outputFile.c_str()); return 0; }
		}
			
		delete distCalculator;
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "PairwiseSeqsCommand", "execute");
		exit(1);
	}
}

/**************************************************************************************************/
void PairwiseSeqsCommand::createProcesses(string filename) {
	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 0;
		processIDS.clear();
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS[lines[process]->end] = pid;  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				if (output != "square") {  driver(lines[process]->start, lines[process]->end, filename + toString(getpid()) + ".temp", cutoff); }
				else { driver(lines[process]->start, lines[process]->end, filename + toString(getpid()) + ".temp", "square"); }
				exit(0);
			}else { m->mothurOut("unable to spawn the necessary processes."); m->mothurOutEndLine(); exit(0); }
		}
	
		//force parent to wait until all the processes are done
		for (map<int, int>::iterator it = processIDS.begin(); it != processIDS.end(); it++) { 
			int temp = it->second;
			wait(&temp);
		}
#endif
	}
	catch(exception& e) {
		m->errorOut(e, "PairwiseSeqsCommand", "createProcesses");
		exit(1);
	}
}

/**************************************************************************************************/
/////// need to fix to work with calcs and sequencedb
int PairwiseSeqsCommand::driver(int startLine, int endLine, string dFileName, float cutoff){
	try {

		int startTime = time(NULL);
		
		//column file
		ofstream outFile(dFileName.c_str(), ios::trunc);
		outFile.setf(ios::fixed, ios::showpoint);
		outFile << setprecision(4);
		
		if((output == "lt") && startLine == 0){	outFile << alignDB.getNumSeqs() << endl;	}
		
		for(int i=startLine;i<endLine;i++){
			if(output == "lt")	{	
				string name = alignDB.get(i).getName();
				if (name.length() < 10) { //pad with spaces to make compatible
					while (name.length() < 10) {  name += " ";  }
				}
				outFile << name << '\t';	
			}
			
			for(int j=0;j<i;j++){
				
				if (m->control_pressed) { outFile.close(); return 0;  }
				
				if (alignDB.get(i).getUnaligned().length() > alignment->getnRows()) {
					alignment->resize(alignDB.get(i).getUnaligned().length()+1);
				}
				
				if (alignDB.get(j).getUnaligned().length() > alignment->getnRows()) {
					alignment->resize(alignDB.get(j).getUnaligned().length()+1);
				}
				
				Sequence* seqI = new Sequence(alignDB.get(i).getName(), alignDB.get(i).getAligned());
				Sequence* seqJ = new Sequence(alignDB.get(j).getName(), alignDB.get(j).getAligned());
				
				Nast(alignment, seqI, seqJ);
				
				distCalculator->calcDist(*seqI, *seqJ);
				double dist = distCalculator->getDist();
				
				delete seqI; delete seqJ;
				
				if(dist <= cutoff){
					if (output == "column") { outFile << alignDB.get(i).getName() << ' ' << alignDB.get(j).getName() << ' ' << dist << endl; }
				}
				if (output == "lt") {  outFile << dist << '\t'; }
			}
			
			if (output == "lt") { outFile << endl; }
			
			if(i % 100 == 0){
				m->mothurOut(toString(i) + "\t" + toString(time(NULL) - startTime)); m->mothurOutEndLine();
			}
			
		}
		m->mothurOut(toString(endLine-1) + "\t" + toString(time(NULL) - startTime)); m->mothurOutEndLine();
		
		outFile.close();
		
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "PairwiseSeqsCommand", "driver");
		exit(1);
	}
}
/**************************************************************************************************/
/////// need to fix to work with calcs and sequencedb
int PairwiseSeqsCommand::driver(int startLine, int endLine, string dFileName, string square){
	try {

		int startTime = time(NULL);
		
		//column file
		ofstream outFile(dFileName.c_str(), ios::trunc);
		outFile.setf(ios::fixed, ios::showpoint);
		outFile << setprecision(4);
		
		if(startLine == 0){	outFile << alignDB.getNumSeqs() << endl;	}
		
		for(int i=startLine;i<endLine;i++){
				
			string name = alignDB.get(i).getName();
			//pad with spaces to make compatible
			if (name.length() < 10) { while (name.length() < 10) {  name += " ";  } }
				
			outFile << name << '\t';	
			
			for(int j=0;j<alignDB.getNumSeqs();j++){
				
				if (m->control_pressed) { outFile.close(); return 0;  }
				
				if (alignDB.get(i).getUnaligned().length() > alignment->getnRows()) {
					alignment->resize(alignDB.get(i).getUnaligned().length()+1);
				}
				
				if (alignDB.get(j).getUnaligned().length() > alignment->getnRows()) {
					alignment->resize(alignDB.get(j).getUnaligned().length()+1);
				}
				
				Sequence* seqI = new Sequence(alignDB.get(i).getName(), alignDB.get(i).getAligned());
				Sequence* seqJ = new Sequence(alignDB.get(j).getName(), alignDB.get(j).getAligned());
				
				Nast(alignment, seqI, seqJ);
				
				distCalculator->calcDist(*seqI, *seqJ);
				double dist = distCalculator->getDist();
				
				delete seqI; delete seqJ;
								
				outFile << dist << '\t'; 
			}
			
			outFile << endl; 
			
			if(i % 100 == 0){
				m->mothurOut(toString(i) + "\t" + toString(time(NULL) - startTime)); m->mothurOutEndLine();
			}
			
		}
		m->mothurOut(toString(endLine-1) + "\t" + toString(time(NULL) - startTime)); m->mothurOutEndLine();
		
		outFile.close();
		
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "PairwiseSeqsCommand", "driver");
		exit(1);
	}
}
#ifdef USE_MPI
/**************************************************************************************************/
/////// need to fix to work with calcs and sequencedb
int PairwiseSeqsCommand::driverMPI(int startLine, int endLine, MPI_File& outMPI, float cutoff){
	try {
		MPI_Status status;
		int startTime = time(NULL);
		
		string outputString = "";
		
		for(int i=startLine;i<endLine;i++){
	
			for(int j=0;j<i;j++){
				
				if (m->control_pressed) {  return 0;  }
				
				if (alignDB.get(i).getUnaligned().length() > alignment->getnRows()) {
					alignment->resize(alignDB.get(i).getUnaligned().length()+1);
				}
				
				if (alignDB.get(j).getUnaligned().length() > alignment->getnRows()) {
					alignment->resize(alignDB.get(j).getUnaligned().length()+1);
				}
				
				Sequence* seqI = new Sequence(alignDB.get(i).getName(), alignDB.get(i).getAligned());
				Sequence* seqJ = new Sequence(alignDB.get(j).getName(), alignDB.get(j).getAligned());
				
				Nast(alignment, seqI, seqJ);
				
				distCalculator->calcDist(*seqI, *seqJ);
				double dist = distCalculator->getDist();
				
				delete seqI; delete seqJ;
				
				if(dist <= cutoff){
					 outputString += (alignDB.get(i).getName() + ' ' + alignDB.get(j).getName() + ' ' + toString(dist) + '\n'); 
				}
			}
			
			if(i % 100 == 0){
				//m->mothurOut(toString(i) + "\t" + toString(time(NULL) - startTime)); m->mothurOutEndLine();
				cout << i << '\t' << (time(NULL) - startTime) << endl;
			}
			
			 
			//send results to parent
			int length = outputString.length();

			char* buf = new char[length];
			memcpy(buf, outputString.c_str(), length);
			
			MPI_File_write_shared(outMPI, buf, length, MPI_CHAR, &status);
			outputString = "";
			delete buf;
			
		}
		
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "PairwiseSeqsCommand", "driverMPI");
		exit(1);
	}
}
/**************************************************************************************************/
/////// need to fix to work with calcs and sequencedb
int PairwiseSeqsCommand::driverMPI(int startLine, int endLine, string file, unsigned long int& size){
	try {
		MPI_Status status;
		
		MPI_File outMPI;
		int amode=MPI_MODE_CREATE|MPI_MODE_WRONLY; 

		char filename[1024];
		strcpy(filename, file.c_str());

		MPI_File_open(MPI_COMM_SELF, filename, amode, MPI_INFO_NULL, &outMPI);

		int startTime = time(NULL);
		
		string outputString = "";
		size = 0;
		
		if(startLine == 0){	outputString += toString(alignDB.getNumSeqs()) + "\n";	}
		
		for(int i=startLine;i<endLine;i++){
				
			string name = alignDB.get(i).getName();
			if (name.length() < 10) { //pad with spaces to make compatible
				while (name.length() < 10) {  name += " ";  }
			}
			outputString += name + "\t";	
			
			for(int j=0;j<i;j++){
				
				if (m->control_pressed) {  return 0;  }
				
				if (alignDB.get(i).getUnaligned().length() > alignment->getnRows()) {
					alignment->resize(alignDB.get(i).getUnaligned().length()+1);
				}
				
				if (alignDB.get(j).getUnaligned().length() > alignment->getnRows()) {
					alignment->resize(alignDB.get(j).getUnaligned().length()+1);
				}
				
				Sequence* seqI = new Sequence(alignDB.get(i).getName(), alignDB.get(i).getAligned());
				Sequence* seqJ = new Sequence(alignDB.get(j).getName(), alignDB.get(j).getAligned());
				
				Nast(alignment, seqI, seqJ);
				
				distCalculator->calcDist(*seqI, *seqJ);
				double dist = distCalculator->getDist();
				
				delete seqI; delete seqJ;
				
				outputString += toString(dist) + "\t"; 
			}
			
			outputString += "\n"; 
			
			//send results to parent
			int length = outputString.length();
			char* buf = new char[length];
			memcpy(buf, outputString.c_str(), length);
			
			MPI_File_write(outMPI, buf, length, MPI_CHAR, &status);
			size += outputString.length();
			outputString = "";
			delete buf;
		}
		
		MPI_File_close(&outMPI);
		
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "PairwiseSeqsCommand", "driverMPI");
		exit(1);
	}
}
/**************************************************************************************************/
/////// need to fix to work with calcs and sequencedb
int PairwiseSeqsCommand::driverMPI(int startLine, int endLine, string file, unsigned long int& size, string square){
	try {
		MPI_Status status;
		
		MPI_File outMPI;
		int amode=MPI_MODE_CREATE|MPI_MODE_WRONLY; 

		char filename[1024];
		strcpy(filename, file.c_str());

		MPI_File_open(MPI_COMM_SELF, filename, amode, MPI_INFO_NULL, &outMPI);
		
		int startTime = time(NULL);
		
		string outputString = "";
		size = 0;
		
		if(startLine == 0){	outputString += toString(alignDB.getNumSeqs()) + "\n";	}
		
		for(int i=startLine;i<endLine;i++){
				
			string name = alignDB.get(i).getName();
			if (name.length() < 10) { //pad with spaces to make compatible
				while (name.length() < 10) {  name += " ";  }
			}
			outputString += name + "\t";	
			
			for(int j=0;j<alignDB.getNumSeqs();j++){
				
				if (m->control_pressed) {  return 0;  }
				
				if (alignDB.get(i).getUnaligned().length() > alignment->getnRows()) {
					alignment->resize(alignDB.get(i).getUnaligned().length()+1);
				}
				
				if (alignDB.get(j).getUnaligned().length() > alignment->getnRows()) {
					alignment->resize(alignDB.get(j).getUnaligned().length()+1);
				}
				
				Sequence* seqI = new Sequence(alignDB.get(i).getName(), alignDB.get(i).getAligned());
				Sequence* seqJ = new Sequence(alignDB.get(j).getName(), alignDB.get(j).getAligned());
				
				Nast(alignment, seqI, seqJ);
				
				distCalculator->calcDist(*seqI, *seqJ);
				double dist = distCalculator->getDist();
				
				delete seqI; delete seqJ;
				
				outputString += toString(dist) + "\t"; 
			}
			
			outputString += "\n"; 

			//send results to parent
			int length = outputString.length();
			char* buf = new char[length];
			memcpy(buf, outputString.c_str(), length);
			
			MPI_File_write(outMPI, buf, length, MPI_CHAR, &status);
			size += outputString.length();
			outputString = "";
			delete buf;
		}
		
		MPI_File_close(&outMPI);
		
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "PairwiseSeqsCommand", "driverMPI");
		exit(1);
	}
}
#endif
/**************************************************************************************************/

