/*
 *  filterseqscommand.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/4/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "filterseqscommand.h"
#include "sequence.hpp"

/**************************************************************************************/

FilterSeqsCommand::FilterSeqsCommand(string option)  {
	try {
		abort = false;
		filterFileName = "";
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta", "trump", "soft", "hard", "vertical", "outputdir","inputdir", "processors"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("hard");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["hard"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			fasta = validParameter.validFile(parameters, "fasta", false);
			if (fasta == "not found") { m->mothurOut("fasta is a required parameter for the filter.seqs command."); m->mothurOutEndLine(); abort = true;  }
			else { 
				splitAtDash(fasta, fastafileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastafileNames.size(); i++) {
					if (inputDir != "") {
						string path = hasPath(fastafileNames[i]);
						//if the user has not given a path then, add inputdir. else leave path alone.
						if (path == "") {	fastafileNames[i] = inputDir + fastafileNames[i];		}
					}

					int ableToOpen;
					ifstream in;
					ableToOpen = openInputFile(fastafileNames[i], in);
					if (ableToOpen == 1) { 
						m->mothurOut(fastafileNames[i] + " will be disregarded."); m->mothurOutEndLine(); 
						//erase from file list
						fastafileNames.erase(fastafileNames.begin()+i);
						i--;
					}else{  
						string simpleName = getSimpleName(fastafileNames[i]);
						filterFileName += simpleName.substr(0, simpleName.find_first_of('.'));
					}
					in.close();
				}
				
				//make sure there is at least one valid file left
				if (fastafileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			if (!abort) {
				//if the user changes the output directory command factory will send this info to us in the output parameter 
				outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
					outputDir = "";	
					outputDir += hasPath(fastafileNames[0]); //if user entered a file with a path then preserve it	
				}
			}
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			
			string temp;
			temp = validParameter.validFile(parameters, "trump", false);			if (temp == "not found") { temp = "*"; }
			trump = temp[0];
			
			temp = validParameter.validFile(parameters, "soft", false);				if (temp == "not found") { soft = 0; }
			else {  soft = (float)atoi(temp.c_str()) / 100.0;  }
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = "1";				}
			convert(temp, processors); 
			
			hard = validParameter.validFile(parameters, "hard", true);				if (hard == "not found") { hard = ""; }
			else if (hard == "not open") { abort = true; }	
			
			vertical = validParameter.validFile(parameters, "vertical", false);		if (vertical == "not found") { vertical = "T"; }
			
			numSeqs = 0;
			
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "FilterSeqsCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void FilterSeqsCommand::help(){
	try {
		m->mothurOut("The filter.seqs command reads a file containing sequences and creates a .filter and .filter.fasta file.\n");
		m->mothurOut("The filter.seqs command parameters are fasta, trump, soft, hard and vertical. \n");
		m->mothurOut("The fasta parameter is required. You may enter several fasta files to build the filter from and filter, by separating their names with -'s.\n");
		m->mothurOut("For example: fasta=abrecovery.fasta-amazon.fasta \n");
		m->mothurOut("The trump parameter .... The default is ...\n");
		m->mothurOut("The soft parameter .... The default is ....\n");
		m->mothurOut("The hard parameter .... The default is ....\n");
		m->mothurOut("The vertical parameter .... The default is T.\n");
		m->mothurOut("The filter.seqs command should be in the following format: \n");
		m->mothurOut("filter.seqs(fasta=yourFastaFile, trump=yourTrump, soft=yourSoft, hard=yourHard, vertical=yourVertical) \n");
		m->mothurOut("Example filter.seqs(fasta=abrecovery.fasta, trump=..., soft=..., hard=..., vertical=T).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n\n");
		
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "help");
		exit(1);
	}
}

/**************************************************************************************/

int FilterSeqsCommand::execute() {	
	try {
	
		if (abort == true) { return 0; }
		vector<string> outputNames;
		
		ifstream inFASTA;
		openInputFile(fastafileNames[0], inFASTA);
		
		Sequence testSeq(inFASTA);
		alignmentLength = testSeq.getAlignLength();
		inFASTA.close();
		
		////////////create filter/////////////////
		
		filter = createFilter();
		
		ofstream outFilter;
		
		string filterFile = outputDir + filterFileName + ".filter";
		openOutputFile(filterFile, outFilter);
		outFilter << filter << endl;
		outFilter.close();
		outputNames.push_back(filterFile);
		
		
		////////////run filter/////////////////
		
		numSeqs = 0;
		for (int i = 0; i < fastafileNames.size(); i++) {
			ifstream in;
			openInputFile(fastafileNames[i], in);
			string filteredFasta = outputDir + getRootName(getSimpleName(fastafileNames[i])) + "filter.fasta";
			ofstream outFASTA;
			openOutputFile(filteredFasta, outFASTA);
			outputNames.push_back(filteredFasta);
			
			
			while(!in.eof()){
				if (m->control_pressed) { in.close(); outFASTA.close(); for(int i = 0; i < outputNames.size(); i++) { remove(outputNames[i].c_str()); }  return 0; }
				
				Sequence seq(in);
				if (seq.getName() != "") {
					string align = seq.getAligned();
					string filterSeq = "";
					
					for(int j=0;j<alignmentLength;j++){
						if(filter[j] == '1'){
							filterSeq += align[j];
						}
					}
					
					outFASTA << '>' << seq.getName() << endl << filterSeq << endl;
					numSeqs++;
				}
				gobble(in);
			}
			outFASTA.close();
			in.close();
		}
		
		int filteredLength = 0;
		for(int i=0;i<alignmentLength;i++){
			if(filter[i] == '1'){	filteredLength++;	}
		}
		
		if (m->control_pressed) {  for(int i = 0; i < outputNames.size(); i++) { remove(outputNames[i].c_str()); }  return 0; }

		
		m->mothurOutEndLine();
		m->mothurOut("Length of filtered alignment: " + toString(filteredLength)); m->mothurOutEndLine();
		m->mothurOut("Number of columns removed: " + toString((alignmentLength-filteredLength))); m->mothurOutEndLine();
		m->mothurOut("Length of the original alignment: " + toString(alignmentLength)); m->mothurOutEndLine();
		m->mothurOut("Number of sequences used to construct filter: " + toString(numSeqs)); m->mothurOutEndLine();
		
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for(int i = 0; i < outputNames.size(); i++) {  m->mothurOut(outputNames[i]); m->mothurOutEndLine();	 }
		m->mothurOutEndLine();

		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************/
string FilterSeqsCommand::createFilter() {	
	try {
		string filterString = "";
		
		Filters F;
		
		if (soft != 0)			{  F.setSoft(soft);		}
		if (trump != '*')		{  F.setTrump(trump);	}
		
		F.setLength(alignmentLength);
		
		if(soft != 0 || isTrue(vertical)){
			F.initialize();
		}
		
		if(hard.compare("") != 0)	{	F.doHard(hard);		}
		else						{	F.setFilter(string(alignmentLength, '1'));	}
		
		numSeqs = 0;
		
		if(trump != '*' || isTrue(vertical) || soft != 0){
			for (int s = 0; s < fastafileNames.size(); s++) {
			
				for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
				
				#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
					if(processors == 1){
						ifstream inFASTA;
						openInputFile(fastafileNames[s], inFASTA);
						int numFastaSeqs=count(istreambuf_iterator<char>(inFASTA),istreambuf_iterator<char>(), '>');
						inFASTA.close();
						
						numSeqs += numFastaSeqs;
				
						lines.push_back(new linePair(0, numFastaSeqs));
				
						driverCreateFilter(F, fastafileNames[s], lines[0]);
					}else{
						vector<int> positions;
				
						ifstream inFASTA;
						openInputFile(fastafileNames[s], inFASTA);
				
						string input;
						while(!inFASTA.eof()){
							input = getline(inFASTA);
							if (input.length() != 0) {
								if(input[0] == '>'){	long int pos = inFASTA.tellg(); positions.push_back(pos - input.length() - 1);	}
							}
						}
						inFASTA.close();
				
						int numFastaSeqs = positions.size();
						
						numSeqs += numFastaSeqs;
				
						int numSeqsPerProcessor = numFastaSeqs / processors;
				
						for (int i = 0; i < processors; i++) {
							long int startPos = positions[ i * numSeqsPerProcessor ];
							if(i == processors - 1){
								numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor;
							}
							lines.push_back(new linePair(startPos, numSeqsPerProcessor));
						}
				
						createProcessesCreateFilter(F, fastafileNames[s]); 
					}
				#else
					ifstream inFASTA;
					openInputFile(fastafileNames[s], inFASTA);
					int numFastaSeqs=count(istreambuf_iterator<char>(inFASTA),istreambuf_iterator<char>(), '>');
					inFASTA.close();
						
					numSeqs += numFastaSeqs;
				
					lines.push_back(new linePair(0, numFastaSeqs));
				
					driverCreateFilter(F, lines[0], fastafileNames[s]);
				#endif
			
			
			}
		}
		
		F.setNumSeqs(numSeqs);
				
		if(isTrue(vertical) == 1)	{	F.doVertical();	}
		if(soft != 0)				{	F.doSoft();		}
		
		filterString = F.getFilter();

		return filterString;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "createFilter");
		exit(1);
	}
}
/**************************************************************************************/
int FilterSeqsCommand::driverCreateFilter(Filters& F, string filename, linePair* line) {	
	try {
		
		ifstream in;
		openInputFile(filename, in);
				
		in.seekg(line->start);
		
		for(int i=0;i<line->numSeqs;i++){
				
			if (m->control_pressed) { in.close(); return 1; }
					
			Sequence seq(in);
			if (seq.getName() != "") {
					if(trump != '*'){	F.doTrump(seq);	}
					if(isTrue(vertical) || soft != 0){	F.getFreqs(seq);	}
					cout.flush();
			}
		}
				
		in.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "driverCreateFilter");
		exit(1);
	}
}
/**************************************************************************************************/

int FilterSeqsCommand::createProcessesCreateFilter(Filters& F, string filename) {
	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 0;
		int exitCommand = 1;
		vector<int> processIDS;
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = vfork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				driverCreateFilter(F, filename, lines[process]);
				exit(0);
			}else { m->mothurOut("unable to spawn the necessary processes."); m->mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		return exitCommand;
#endif		
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "createProcessesCreateFilter");
		exit(1);
	}
}
/**************************************************************************************/
