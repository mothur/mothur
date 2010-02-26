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
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta", "trump", "soft", "hard", "vertical", "outputdir","inputdir"};
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
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not found") { m->mothurOut("fasta is a required parameter for the filter.seqs command."); m->mothurOutEndLine(); abort = true; }
			else if (fastafile == "not open") { abort = true; }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += hasPath(fastafile); //if user entered a file with a path then preserve it	
			}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			
			string temp;
			temp = validParameter.validFile(parameters, "trump", false);			if (temp == "not found") { temp = "*"; }
			trump = temp[0];
			
			temp = validParameter.validFile(parameters, "soft", false);				if (temp == "not found") { soft = 0; }
			else {  soft = (float)atoi(temp.c_str()) / 100.0;  }
			
			hard = validParameter.validFile(parameters, "hard", true);				if (hard == "not found") { hard = ""; }
			else if (hard == "not open") { abort = true; }	
			
			vertical = validParameter.validFile(parameters, "vertical", false);		if (vertical == "not found") { vertical = "T"; }
			
			numSeqs = 0;
			
			if (abort == false) {
			
				if (soft != 0)			{  F.setSoft(soft);		}
				if (trump != '*')		{  F.setTrump(trump);	}
			
			}
						
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
		m->mothurOut("The fasta parameter is required.\n");
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
		
		ifstream inFASTA;
		openInputFile(fastafile, inFASTA);
		
		Sequence testSeq(inFASTA);
		alignmentLength = testSeq.getAlignLength();
		inFASTA.seekg(0);
		
		F.setLength(alignmentLength);
		
		if(soft != 0 || isTrue(vertical)){
			F.initialize();
		}
		
		if(hard.compare("") != 0)	{	F.doHard(hard);		}
		else						{	F.setFilter(string(alignmentLength, '1'));	}

		if(trump != '*' || isTrue(vertical) || soft != 0){
			while(!inFASTA.eof()){	//read through and create the filter...
				Sequence seq(inFASTA);
				if (seq.getName() != "") {
					if(trump != '*'){	F.doTrump(seq);	}
					if(isTrue(vertical) || soft != 0){	F.getFreqs(seq);	}
					numSeqs++;
					cout.flush();
				}
			}
		
		}
		inFASTA.close();
		F.setNumSeqs(numSeqs);
		
		
		if(isTrue(vertical) == 1)	{	F.doVertical();	}
		if(soft != 0)				{	F.doSoft();		}
		
		filter = F.getFilter();

		ofstream outFilter;
		string filterFile = outputDir + getRootName(getSimpleName(fastafile)) + "filter";
		openOutputFile(filterFile, outFilter);
		outFilter << filter << endl;
		outFilter.close();
		
		ifstream inFasta2;
		openInputFile(fastafile, inFasta2);
		string filteredFasta = outputDir + getRootName(getSimpleName(fastafile)) + "filter.fasta";
		ofstream outFASTA;
		openOutputFile(filteredFasta, outFASTA);

		numSeqs = 0;
		while(!inFasta2.eof()){
			Sequence seq(inFasta2);
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
			gobble(inFasta2);
		}
		outFASTA.close();
		inFasta2.close();
		
		
		int filteredLength = 0;
		for(int i=0;i<alignmentLength;i++){
			if(filter[i] == '1'){	filteredLength++;	}
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Length of filtered alignment: " + toString(filteredLength)); m->mothurOutEndLine();
		m->mothurOut("Number of columns removed: " + toString((alignmentLength-filteredLength))); m->mothurOutEndLine();
		m->mothurOut("Length of the original alignment: " + toString(alignmentLength)); m->mothurOutEndLine();
		m->mothurOut("Number of sequences used to construct filter: " + toString(numSeqs)); m->mothurOutEndLine();
		
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOut(filterFile); m->mothurOutEndLine();	
		m->mothurOut(filteredFasta); m->mothurOutEndLine();
		m->mothurOutEndLine();

		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "FilterSeqsCommand", "execute");
		exit(1);
	}
}

/**************************************************************************************/
