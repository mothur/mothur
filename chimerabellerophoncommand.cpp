/*
 *  chimerabellerophoncommand.cpp
 *  Mothur
 *
 *  Created by westcott on 4/1/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "chimerabellerophoncommand.h"
#include "bellerophon.h"

//**********************************************************************************************************************
vector<string> ChimeraBellerophonCommand::getValidParameters(){	
	try {
		string AlignArray[] =  {"fasta","filter","correction","processors","window","increment","outputdir","inputdir"};
		vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraBellerophonCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> ChimeraBellerophonCommand::getRequiredParameters(){	
	try {
		string AlignArray[] =  {"fasta"};
		vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraBellerophonCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> ChimeraBellerophonCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraBellerophonCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************
ChimeraBellerophonCommand::ChimeraBellerophonCommand(){	
	try {
		abort = true;
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["chimera"] = tempOutNames;
		outputTypes["accnos"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraBellerophonCommand", "ChimeraBellerophonCommand");
		exit(1);
	}
}
//***************************************************************************************************************
ChimeraBellerophonCommand::ChimeraBellerophonCommand(string option)  {
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta","filter","correction","processors","window","increment","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter("chimera.bellerophon");
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["chimera"] = tempOutNames;
			outputTypes["accnos"] = tempOutNames;
		
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			
			fastafile = validParameter.validFile(parameters, "fasta", false);
			if (fastafile == "not found") { fastafile = ""; m->mothurOut("fasta is a required parameter for the chimera.bellerophon command."); m->mothurOutEndLine(); abort = true;  }
			else { 
				m->splitAtDash(fastafile, fastaFileNames);
				
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
					
					//if you can't open it, try default location
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
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}

			string temp;
			temp = validParameter.validFile(parameters, "filter", false);			if (temp == "not found") { temp = "F"; }
			filter = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "correction", false);		if (temp == "not found") { temp = "T"; }
			correction = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "processors", false);		if (temp == "not found") { temp = "1"; }
			convert(temp, processors);
			
			temp = validParameter.validFile(parameters, "window", false);			if (temp == "not found") { temp = "0"; }
			convert(temp, window);
			
			temp = validParameter.validFile(parameters, "increment", false);		if (temp == "not found") { temp = "25"; }
			convert(temp, increment);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraBellerophonCommand", "ChimeraBellerophonCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void ChimeraBellerophonCommand::help(){
	try {
		m->mothurOut("The chimera.bellerophon command reads a fastafile and creates list of potentially chimeric sequences.\n");
		m->mothurOut("The chimera.bellerophon command parameters are fasta, filter, correction, processors, window, increment. The fasta parameter is required.\n");
		m->mothurOut("The fasta parameter is required.  You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amzon.fasta \n");
		m->mothurOut("The filter parameter allows you to specify if you would like to apply a vertical and 50% soft filter, default=false. \n");
		m->mothurOut("The correction parameter allows you to put more emphasis on the distance between highly similar sequences and less emphasis on the differences between remote homologs, default=true.\n");
		m->mothurOut("The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n");
		#ifdef USE_MPI
		m->mothurOut("When using MPI, the processors parameter is set to the number of MPI processes running. \n");
		#endif
		m->mothurOut("The window parameter allows you to specify the window size for searching for chimeras, default is 1/4 sequence length. \n");
		m->mothurOut("The increment parameter allows you to specify how far you move each window while finding chimeric sequences, default is 25.\n");
		m->mothurOut("chimera.bellerophon(fasta=yourFastaFile, filter=yourFilter, correction=yourCorrection, processors=yourProcessors) \n");
		m->mothurOut("Example: chimera.bellerophon(fasta=AD.align, filter=True, correction=true, window=200) \n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n\n");	
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraBellerophonCommand", "help");
		exit(1);
	}
}

//***************************************************************************************************************

ChimeraBellerophonCommand::~ChimeraBellerophonCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************

int ChimeraBellerophonCommand::execute(){
	try{
		
		if (abort == true) { return 0; }
		
		for (int i = 0; i < fastaFileNames.size(); i++) {
			
			m->mothurOut("Checking sequences from " + fastaFileNames[i] + " ..." ); m->mothurOutEndLine();
			
			int start = time(NULL);	
			
			chimera = new Bellerophon(fastaFileNames[i], filter, correction, window, increment, processors, outputDir);	
			
			if (outputDir == "") { outputDir = m->hasPath(fastaFileNames[i]);  }//if user entered a file with a path then preserve it		
			string outputFileName = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[i])) +  "bellerophon.chimeras";
			string accnosFileName = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[i])) + "bellerophon.accnos";
			
			chimera->getChimeras();
			
			if (m->control_pressed) { delete chimera; for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());	} outputTypes.clear(); return 0;	}
			
		#ifdef USE_MPI
			MPI_File outMPI;
			MPI_File outMPIAccnos;
			
			int outMode=MPI_MODE_CREATE|MPI_MODE_WRONLY; 
			
			char outFilename[1024];
			strcpy(outFilename, accnosFileName.c_str());

			char FileName[1024];
			strcpy(FileName, outputFileName.c_str());

			MPI_File_open(MPI_COMM_WORLD, FileName, outMode, MPI_INFO_NULL, &outMPI);  //comm, filename, mode, info, filepointer
			MPI_File_open(MPI_COMM_WORLD, outFilename, outMode, MPI_INFO_NULL, &outMPIAccnos);
	
			numSeqs = chimera->print(outMPI, outMPIAccnos);
			
			MPI_File_close(&outMPI);
			MPI_File_close(&outMPIAccnos);

		#else
		
			ofstream out;
			m->openOutputFile(outputFileName, out);
			
			ofstream out2;
			m->openOutputFile(accnosFileName, out2);
			
			numSeqs = chimera->print(out, out2);
			out.close();
			out2.close(); 
			
		#endif
			
			if (m->control_pressed) { remove(accnosFileName.c_str()); remove(outputFileName.c_str()); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());	} outputTypes.clear(); delete chimera;	return 0;	}
			
			m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences.");	m->mothurOutEndLine(); m->mothurOutEndLine();
			
			outputNames.push_back(outputFileName);  outputTypes["chimera"].push_back(outputFileName);
			outputNames.push_back(accnosFileName);  outputTypes["accnos"].push_back(accnosFileName);
			
			delete chimera;
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
				
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraBellerophonCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

