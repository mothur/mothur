/*
 *  clearcutcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 5/11/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "clearcutcommand.h"
#ifdef __cplusplus
extern "C" {
#endif
#include "clearcut.h"
#ifdef __cplusplus
}
#endif

//**********************************************************************************************************************
vector<string> ClearcutCommand::setParameters(){	
	try {
		CommandParameter pphylip("phylip", "InputTypes", "", "", "FastaPhylip", "FastaPhylip", "none",false,false); parameters.push_back(pphylip);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "FastaPhylip", "FastaPhylip", "none",false,false); parameters.push_back(pfasta);
		CommandParameter pverbose("verbose", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pverbose);
		CommandParameter pquiet("quiet", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pquiet);
		CommandParameter pversion("version", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pversion);
		CommandParameter pseed("seed", "String", "", "", "*", "", "",false,false); parameters.push_back(pseed);
		CommandParameter pnorandom("norandom", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pnorandom);
		CommandParameter pshuffle("shuffle", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pshuffle);
		CommandParameter pneighbor("neighbor", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(pneighbor);
		CommandParameter pexpblen("expblen", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pexpblen);
		CommandParameter pexpdist("expdist", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pexpdist);
		CommandParameter pDNA("DNA", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pDNA);
		CommandParameter pprotein("protein", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pprotein);
		CommandParameter pjukes("jukes", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pjukes);
		CommandParameter pkimura("kimura", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pkimura);
		CommandParameter pstdout("stdout", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pstdout);
		CommandParameter pntrees("ntrees", "Number", "", "1", "", "", "",false,false); parameters.push_back(pntrees);
		CommandParameter pmatrixout("matrixout", "String", "", "", "", "", "",false,false); parameters.push_back(pmatrixout);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClearcutCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ClearcutCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The clearcut command interfaces mothur with the clearcut program written by Initiative for Bioinformatics and Evolutionary Studies (IBEST) at the University of Idaho.\n";
		helpString += "For more information about clearcut refer to http://bioinformatics.hungry.com/clearcut/ \n";
		helpString += "The clearcut command parameters are phylip, fasta, version, verbose, quiet, seed, norandom, shuffle, neighbor, expblen, expdist, ntrees, matrixout, stdout, kimura, jukes, protein, DNA. \n";
		helpString += "The phylip parameter allows you to enter your phylip formatted distance matrix. \n";
		helpString += "The fasta parameter allows you to enter your aligned fasta file, if you enter a fastafile you specify if the sequences are DNA or protein using the DNA or protein parameters. \n";
		
		helpString += "The version parameter prints out the version of clearcut you are using, default=F. \n";
		helpString += "The verbose parameter prints out more output from clearcut, default=F. \n";
		helpString += "The quiet parameter turns on silent operation mode, default=F. \n";
		helpString += "The seed parameter allows you to explicitly set the PRNG seed to a specific value. \n";
		helpString += "The norandom parameter allows you to attempt joins deterministically, default=F. \n";
		helpString += "The shuffle parameter allows you to randomly shuffle the distance matrix, default=F. \n";
		helpString += "The neighbor parameter allows you to use traditional Neighbor-Joining algorithm, default=T. \n";
		
		helpString += "The DNA parameter allows you to indicate your fasta file contains DNA sequences, default=F. \n";
		helpString += "The protein parameter allows you to indicate your fasta file contains protein sequences, default=F. \n";
		
		helpString += "The stdout parameter outputs your tree to STDOUT, default=F. \n";
		helpString += "The matrixout parameter allows you to specify a filename to output a distance matrix to. \n";
		helpString += "The ntrees parameter allows you to specify the number of output trees, default=1. \n";
		helpString += "The expblen parameter allows you to use exponential notation for branch lengths, default=F. \n";
		helpString += "The expdist parameter allows you to use exponential notation for distance outputs, default=F. \n";
		
		helpString += "The clearcut command should be in the following format: \n";
		helpString += "clearcut(phylip=yourDistanceFile) \n";
		helpString += "Example: clearcut(phylip=abrecovery.phylip.dist) \n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ClearcutCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
ClearcutCommand::ClearcutCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["tree"] = tempOutNames;
		outputTypes["matrixout"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ClearcutCommand", "ClearcutCommand");
		exit(1);
	}
}
/**************************************************************************************/
ClearcutCommand::ClearcutCommand(string option)  {	
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["tree"] = tempOutNames;
			outputTypes["matrixout"] = tempOutNames;

			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("phylip");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["phylip"] = inputDir + it->second;		}
				}
			}

			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { fastafile = ""; abort = true; }
			else if (fastafile == "not found") { fastafile = ""; }	
			else { inputFile = fastafile;  m->setFastaFile(fastafile); }
			
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not open") { phylipfile = ""; abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; }
			else { inputFile = phylipfile;  m->setPhylipFile(phylipfile); }
				
			if ((phylipfile == "") && (fastafile == "")) {  
				//is there are current file available for either of these?
				//give priority to phylip, then fasta
				phylipfile = m->getPhylipFile(); 
				if (phylipfile != "") {  inputFile = phylipfile; m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter."); m->mothurOutEndLine(); }
				else { 
					fastafile = m->getFastaFile(); 
					if (fastafile != "") { inputFile = fastafile;  m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("No valid current files. You must provide a phylip or fasta file before you can use the clearcut command."); m->mothurOutEndLine(); 
						abort = true;
					}
				}
			}
			if ((phylipfile != "") && (fastafile != "")) {  m->mothurOut("You must provide either a phylip formatted distance matrix or an aligned fasta file, not BOTH."); m->mothurOutEndLine(); abort=true; }

			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(inputFile);	}
			
			string temp;
			temp = validParameter.validFile(parameters, "version", false);		if (temp == "not found"){	temp = "F";			}
			version = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "verbose", false);		if (temp == "not found"){	temp = "F";			}
			verbose = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "quiet", false);		if (temp == "not found"){	temp = "F";			}
			quiet = m->isTrue(temp); 
			
			seed = validParameter.validFile(parameters, "seed", false);			if (seed == "not found"){	seed = "*";			}
			
			temp = validParameter.validFile(parameters, "norandom", false);		if (temp == "not found"){	temp = "F";			}
			norandom = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "shuffle", false);		if (temp == "not found"){	temp = "F";			}
			shuffle = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "neighbor", false);		if (temp == "not found"){	temp = "T";			}
			neighbor = m->isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "DNA", false);			if (temp == "not found"){	temp = "F";			}
			DNA = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "protein", false);		if (temp == "not found"){	temp = "F";			}
			protein = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "jukes", false);		if (temp == "not found"){	temp = "F";			}
			jukes = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "kimura", false);		if (temp == "not found"){	temp = "F";			}
			kimura = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "stdout", false);		if (temp == "not found"){	temp = "F";			}
			stdoutWanted = m->isTrue(temp); 
			
			matrixout = validParameter.validFile(parameters, "matrixout", false);	if (matrixout == "not found"){	matrixout = "";		}
			
			ntrees = validParameter.validFile(parameters, "ntrees", false);		if (ntrees == "not found"){	ntrees = "1";		}
			
			temp = validParameter.validFile(parameters, "expblen", false);		if (temp == "not found"){	temp = "F";			}
			expblen = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "expdist", false);		if (temp == "not found"){	temp = "F";			}
			expdist = m->isTrue(temp);
			
			if ((fastafile != "") && ((!DNA) && (!protein))) { m->mothurOut("You must specify the type of sequences you are using: DNA or protein"); m->mothurOutEndLine(); abort=true; }
		}

	}
	catch(exception& e) {
		m->errorOut(e, "ClearcutCommand", "ClearcutCommand");
		exit(1);
	}
}
/**************************************************************************************/
int ClearcutCommand::execute() {	
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//prepare filename
		string outputName = outputDir + m->getRootName(m->getSimpleName(inputFile)) + "tre";
		outputNames.push_back(outputName); outputTypes["tree"].push_back(outputName);
		
		vector<char*> cPara;
		
		char* tempClearcut = new char[8];  
		strcpy(tempClearcut, "clearcut"); 
		cPara.push_back(tempClearcut);
				
		//you gave us a distance matrix
		if (phylipfile != "") {  char* temp = new char[10];  strcpy(temp, "--distance");  cPara.push_back(temp);	}
		
		//you gave us a fastafile
		if (fastafile != "") { char* temp = new char[11];  strcpy(temp, "--alignment");  cPara.push_back(temp); 	}
		
		if (version)			{  char* temp = new char[9];  strcpy(temp, "--version");  cPara.push_back(temp);	}
		if (verbose)			{  char* temp = new char[9];  strcpy(temp, "--verbose");  cPara.push_back(temp);	}
		if (quiet)				{  char* temp = new char[7];  strcpy(temp, "--quiet");  cPara.push_back(temp);	}
		if (seed != "*")		{  
			string tempSeed = "--seed=" + seed;
			char* temp = new char[tempSeed.length()];
			strcpy(temp, tempSeed.c_str());
			cPara.push_back(temp);
		}
		if (norandom)			{  char* temp = new char[10];  strcpy(temp, "--norandom");  cPara.push_back(temp);	}
		if (shuffle)			{  char* temp = new char[9];  strcpy(temp, "--shuffle");  cPara.push_back(temp);	}
		if (neighbor)			{  char* temp = new char[10];  strcpy(temp, "--neighbor");  cPara.push_back(temp);	}
		
		string tempIn = "--in=" + inputFile;  
		char* tempI = new char[tempIn.length()];
		strcpy(tempI, tempIn.c_str());
		cPara.push_back(tempI);
		
		if (stdoutWanted)		{  char* temp = new char[8];  strcpy(temp, "--stdout");  cPara.push_back(temp);	}
		else{  
			string tempOut = "--out=" + outputName;  
			
			char* temp = new char[tempOut.length()];
			strcpy(temp, tempOut.c_str());
			cPara.push_back(temp);
		}
			
		if (DNA)				{  char* temp = new char[5];  strcpy(temp, "--DNA");  cPara.push_back(temp);		}
		if (protein)			{  char* temp = new char[9];  strcpy(temp, "--protein");  cPara.push_back(temp);	}
		if (jukes)				{  char* temp = new char[7];  strcpy(temp, "--jukes");  cPara.push_back(temp);		}
		if (kimura)				{ char* temp = new char[8];  strcpy(temp, "--kimura");  cPara.push_back(temp);		}
		if (matrixout != "")	{  
			string tempMatrix =  "--matrixout=" + outputDir + matrixout; 
			char* temp = new char[tempMatrix.length()];
			strcpy(temp, tempMatrix.c_str());
			cPara.push_back(temp);
			outputNames.push_back((outputDir + matrixout));
			outputTypes["matrixout"].push_back((outputDir + matrixout));
		}

		if (ntrees != "1")		{  
			string tempNtrees = "--ntrees=" + ntrees; 
			char* temp = new char[tempNtrees.length()];
			strcpy(temp, tempNtrees.c_str());
			cPara.push_back(temp);
		}

		if (expblen)			{ char* temp = new char[9];  strcpy(temp, "--expblen");  cPara.push_back(temp); 	}
		if (expdist)			{ char* temp = new char[9];  strcpy(temp, "--expdist");  cPara.push_back(temp);	}
		
		char** clearcutParameters;
		clearcutParameters = new char*[cPara.size()];
		for (int i = 0; i < cPara.size(); i++) {  clearcutParameters[i] = cPara[i];  }
		int numArgs = cPara.size();
		
		clearcut_main(numArgs, clearcutParameters); 
		
		//free memory
		for(int i = 0; i < cPara.size(); i++)  {  delete[] cPara[i];  }
		delete[] clearcutParameters; 
		
		if (!stdoutWanted) {	
			
			//set first tree file as new current treefile
			string currentTree = "";
			itTypes = outputTypes.find("tree");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentTree = (itTypes->second)[0]; m->setTreeFile(currentTree); }
			}
			
			m->mothurOutEndLine();
			m->mothurOut("Output File Names: "); m->mothurOutEndLine();
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
			m->mothurOutEndLine();
		}

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClearcutCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************/




