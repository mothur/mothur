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
vector<string> ClearcutCommand::getValidParameters(){	
	try {
		string AlignArray[] =  {"fasta","phylip","version","verbose","quiet","seed","norandom","shuffle","neighbor","expblen",
								"expdist","ntrees","matrixout","stdout","kimura","jukes","protein","DNA","outputdir","inputdir"};
		vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClearcutCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
ClearcutCommand::ClearcutCommand(){	
	try {
		abort = true; calledHelp = true; 
		vector<string> tempOutNames;
		outputTypes["tree"] = tempOutNames;
		outputTypes["matrixout"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ClearcutCommand", "ClearcutCommand");
		exit(1);
	}
}//**********************************************************************************************************************
vector<string> ClearcutCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"fasta","phylip","or"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClearcutCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> ClearcutCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClearcutCommand", "getRequiredFiles");
		exit(1);
	}
}
/**************************************************************************************/
ClearcutCommand::ClearcutCommand(string option)  {	
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta","phylip","version","verbose","quiet","seed","norandom","shuffle","neighbor","expblen",
								"expdist","ntrees","matrixout","stdout","kimura","jukes","protein","DNA","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
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
			else { inputFile = fastafile;  }
			
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not open") { phylipfile = ""; abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; }
			else { inputFile = phylipfile;  }
				
			if ((phylipfile == "") && (fastafile == "")) {  m->mothurOut("You must provide either a phylip formatted distance matrix or an aligned fasta file."); m->mothurOutEndLine(); abort=true; }
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
//**********************************************************************************************************************

void ClearcutCommand::help(){
	try {
		m->mothurOut("The clearcut command interfaces mothur with the clearcut program written by Initiative for Bioinformatics and Evolutionary Studies (IBEST) at the University of Idaho.\n");
		m->mothurOut("For more information about clearcut refer to http://bioinformatics.hungry.com/clearcut/ \n");
		m->mothurOut("The clearcut command parameters are phylip, fasta, version, verbose, quiet, seed, norandom, shuffle, neighbor, expblen, expdist, ntrees, matrixout, stdout, kimura, jukes, protein, DNA. \n");
		m->mothurOut("The phylip parameter allows you to enter your phylip formatted distance matrix. \n");
		m->mothurOut("The fasta parameter allows you to enter your aligned fasta file, if you enter a fastafile you specify if the sequences are DNA or protein using the DNA or protein parameters. \n");
		
		m->mothurOut("The version parameter prints out the version of clearcut you are using, default=F. \n");
		m->mothurOut("The verbose parameter prints out more output from clearcut, default=F. \n");
		m->mothurOut("The quiet parameter turns on silent operation mode, default=F. \n");
		m->mothurOut("The seed parameter allows you to explicitly set the PRNG seed to a specific value. \n");
		m->mothurOut("The norandom parameter allows you to attempt joins deterministically, default=F. \n");
		m->mothurOut("The shuffle parameter allows you to randomly shuffle the distance matrix, default=F. \n");
		m->mothurOut("The neighbor parameter allows you to use traditional Neighbor-Joining algorithm, default=T. \n");
		
		m->mothurOut("The DNA parameter allows you to indicate your fasta file contains DNA sequences, default=F. \n");
		m->mothurOut("The protein parameter allows you to indicate your fasta file contains protein sequences, default=F. \n");
		
		m->mothurOut("The stdout parameter outputs your tree to STDOUT, default=F. \n");
		m->mothurOut("The matrixout parameter allows you to specify a filename to output a distance matrix to. \n");
		m->mothurOut("The ntrees parameter allows you to specify the number of output trees, default=1. \n");
		m->mothurOut("The expblen parameter allows you to use exponential notation for branch lengths, default=F. \n");
		m->mothurOut("The expdist parameter allows you to use exponential notation for distance outputs, default=F. \n");

		m->mothurOut("The clearcut command should be in the following format: \n");
		m->mothurOut("clearcut(phylip=yourDistanceFile) \n");
		m->mothurOut("Example: clearcut(phylip=abrecovery.phylip.dist) \n");	
		
	}
	catch(exception& e) {
		m->errorOut(e, "ClearcutCommand", "help");
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




