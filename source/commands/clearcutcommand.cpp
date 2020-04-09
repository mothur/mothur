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
		CommandParameter pphylip("phylip", "InputTypes", "", "", "FastaPhylip", "FastaPhylip", "none","tree",false,false,true); parameters.push_back(pphylip);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "FastaPhylip", "FastaPhylip", "none","tree",false,false,true); parameters.push_back(pfasta);
		CommandParameter pverbose("verbose", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pverbose);
		CommandParameter pquiet("quiet", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pquiet);
		CommandParameter pversion("version", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pversion);
		CommandParameter prseed("rseed", "String", "", "", "*", "", "","",false,false); parameters.push_back(prseed);
		CommandParameter pnorandom("norandom", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pnorandom);
		CommandParameter pshuffle("shuffle", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pshuffle);
		CommandParameter pneighbor("neighbor", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pneighbor);
		CommandParameter pexpblen("expblen", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pexpblen);
		CommandParameter pexpdist("expdist", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pexpdist);
		CommandParameter pDNA("DNA", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pDNA);
		CommandParameter pprotein("protein", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pprotein);
		CommandParameter pjukes("jukes", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pjukes);
		CommandParameter pkimura("kimura", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pkimura);
		CommandParameter pstdout("stdout", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pstdout);
		CommandParameter pntrees("ntrees", "Number", "", "1", "", "", "","",false,false); parameters.push_back(pntrees);
		CommandParameter pmatrixout("matrixout", "String", "", "", "", "", "","",false,false); parameters.push_back(pmatrixout);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["tree"] = tempOutNames;
        outputTypes["matrixout"] = tempOutNames;
		
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
		helpString += "The rseed parameter allows you to explicitly set the PRNG seed to a specific value. \n";
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
string ClearcutCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "matrixout") {  pattern = "[filename],"; } 
        else if (type == "tree") {  pattern = "[filename],tre"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ClearcutCommand", "getOutputPattern");
        exit(1);
    }
}
/**************************************************************************************/
ClearcutCommand::ClearcutCommand(string option)  {	
	try {
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not open") { fastafile = ""; abort = true; }
			else if (fastafile == "not found") { fastafile = ""; }	
			else { inputFile = fastafile;  current->setFastaFile(fastafile); }
			
			phylipfile = validParameter.validFile(parameters, "phylip");
			if (phylipfile == "not open") { phylipfile = ""; abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; }
			else { inputFile = phylipfile;  current->setPhylipFile(phylipfile); }
				
			if ((phylipfile == "") && (fastafile == "")) {  
				//is there are current file available for either of these?
				//give priority to phylip, then fasta
				phylipfile = current->getPhylipFile(); 
				if (phylipfile != "") {  inputFile = phylipfile; m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter.\n");  }
				else { 
					fastafile = current->getFastaFile(); 
					if (fastafile != "") { inputFile = fastafile;  m->mothurOut("Using " + fastafile + " as input file for the fasta parameter.\n"); }
					else { 
						m->mothurOut("No valid current files. You must provide a phylip or fasta file before you can use the clearcut command.\n");
						abort = true;
					}
				}
			}
			if ((phylipfile != "") && (fastafile != "")) {  m->mothurOut("You must provide either a phylip formatted distance matrix or an aligned fasta file, not BOTH.\n");  abort=true; }

			
			 
			if (outputdir == ""){	outputdir = util.hasPath(inputFile);	}
			
			string temp;
			temp = validParameter.valid(parameters, "version");		if (temp == "not found"){	temp = "F";			}
			version = util.isTrue(temp);
			
			temp = validParameter.valid(parameters, "verbose");		if (temp == "not found"){	temp = "F";			}
			verbose = util.isTrue(temp); 
			
			temp = validParameter.valid(parameters, "quiet");		if (temp == "not found"){	temp = "F";			}
			quiet = util.isTrue(temp); 
			
			seed = validParameter.valid(parameters, "rseed");			if (seed == "not found"){	seed = "*";			}
			
			temp = validParameter.valid(parameters, "norandom");		if (temp == "not found"){	temp = "F";			}
			norandom = util.isTrue(temp); 
			
			temp = validParameter.valid(parameters, "shuffle");		if (temp == "not found"){	temp = "F";			}
			shuffle = util.isTrue(temp); 
			
			temp = validParameter.valid(parameters, "neighbor");		if (temp == "not found"){	temp = "T";			}
			neighbor = util.isTrue(temp); 
			
			temp = validParameter.valid(parameters, "DNA");			if (temp == "not found"){	temp = "F";			}
			DNA = util.isTrue(temp);
			
			temp = validParameter.valid(parameters, "protein");		if (temp == "not found"){	temp = "F";			}
			protein = util.isTrue(temp);
			
			temp = validParameter.valid(parameters, "jukes");		if (temp == "not found"){	temp = "F";			}
			jukes = util.isTrue(temp);
			
			temp = validParameter.valid(parameters, "kimura");		if (temp == "not found"){	temp = "F";			}
			kimura = util.isTrue(temp);
			
			temp = validParameter.valid(parameters, "stdout");		if (temp == "not found"){	temp = "F";			}
			stdoutWanted = util.isTrue(temp); 
			
			matrixout = validParameter.valid(parameters, "matrixout");	if (matrixout == "not found"){	matrixout = "";		}
			
			ntrees = validParameter.valid(parameters, "ntrees");		if (ntrees == "not found"){	ntrees = "1";		}
			
			temp = validParameter.valid(parameters, "expblen");		if (temp == "not found"){	temp = "F";			}
			expblen = util.isTrue(temp);
			
			temp = validParameter.valid(parameters, "expdist");		if (temp == "not found"){	temp = "F";			}
			expdist = util.isTrue(temp);
			
			if ((fastafile != "") && ((!DNA) && (!protein))) { m->mothurOut("You must specify the type of sequences you are using: DNA or protein.\n"); abort=true; }
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
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		//prepare filename
        map<string, string> variables; 
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(inputFile));
		string outputName = getOutputFileName("tree", variables);
		outputNames.push_back(outputName); outputTypes["tree"].push_back(outputName);
		
        int numArgs = 4; //clearcut, in, out and fastafile or phylipfile
        if (version) { numArgs++; } if (verbose) { numArgs++; } if (quiet) { numArgs++; } if (seed != "*")	{ numArgs++; } if (norandom) { numArgs++; }
        if (shuffle) { numArgs++; } if (neighbor) { numArgs++; } if (stdoutWanted) { numArgs++; } if (DNA)	{ numArgs++; } if (protein) { numArgs++; }
        if (jukes) { numArgs++; } if (kimura) { numArgs++; } if (matrixout != "") { numArgs++; } if (ntrees != "1")	{ numArgs++; } if (expblen) { numArgs++; } if (expdist) { numArgs++; }
        
        
        char** clearcutParameters;
		clearcutParameters = new char*[numArgs];
        
        clearcutParameters[0] = new char[9];
		*clearcutParameters[0] = '\0'; strncat(clearcutParameters[0], "clearcut", 8);
				
		//you gave us a distance matrix
        if (phylipfile != "") {  clearcutParameters[1] = new char[11];  *clearcutParameters[1] = '\0'; strncat(clearcutParameters[1], "--distance", 10); 	}
		
		//you gave us a fastafile
		if (fastafile != "") { clearcutParameters[1] = new char[12];  *clearcutParameters[1] = '\0'; strncat(clearcutParameters[1], "--alignment", 11);  	}
		
        int parameterCount = 2;
		if (version)			{  clearcutParameters[parameterCount] = new char[10];  *clearcutParameters[parameterCount] = '\0'; strncat(clearcutParameters[parameterCount], "--version", 9);  	parameterCount++; }
		if (verbose)			{  clearcutParameters[parameterCount] = new char[10];  *clearcutParameters[parameterCount] = '\0'; strncat(clearcutParameters[parameterCount], "--verbose", 9);  parameterCount++;	}
		if (quiet)				{  clearcutParameters[parameterCount] = new char[8];  *clearcutParameters[parameterCount] = '\0'; strncat(clearcutParameters[parameterCount], "--quiet", 7);  parameterCount++;	}
		if (seed != "*")		{  
			string tempSeed = "--seed=" + seed;
			clearcutParameters[parameterCount] = new char[tempSeed.length()+1];
			*clearcutParameters[parameterCount] = '\0'; strncat(clearcutParameters[parameterCount], tempSeed.c_str(), tempSeed.length());
			parameterCount++;
		}
		if (norandom)			{  clearcutParameters[parameterCount] = new char[11];  *clearcutParameters[parameterCount] = '\0'; strncat(clearcutParameters[parameterCount], "--norandom", 10);  parameterCount++;	}
		if (shuffle)			{  clearcutParameters[parameterCount] = new char[10];  *clearcutParameters[parameterCount] = '\0'; strncat(clearcutParameters[parameterCount], "--shuffle", 9);  parameterCount++;	}
		if (neighbor)			{  clearcutParameters[parameterCount] = new char[11];  *clearcutParameters[parameterCount] = '\0'; strncat(clearcutParameters[parameterCount], "--neighbor", 10);  parameterCount++;	}
		
		string tempIn = "--in=" + inputFile;  
        clearcutParameters[parameterCount] = new char[tempIn.length()+1];
		*clearcutParameters[parameterCount] = '\0'; strncat(clearcutParameters[parameterCount], tempIn.c_str(), tempIn.length());
        parameterCount++;
		
		if (stdoutWanted)		{  clearcutParameters[parameterCount] = new char[9];  *clearcutParameters[parameterCount] = '\0'; strncat(clearcutParameters[parameterCount], "--stdout", 8);  parameterCount++;	}
		else{  
            string tempOut = "--out=" + outputName;
			clearcutParameters[parameterCount] = new char[tempOut.length()+1];
			*clearcutParameters[parameterCount] = '\0'; strncat(clearcutParameters[parameterCount], tempOut.c_str(), tempOut.length());
			parameterCount++;
            
		}
			
		if (DNA)				{  clearcutParameters[parameterCount] = new char[6];  *clearcutParameters[parameterCount] = '\0'; strncat(clearcutParameters[parameterCount], "--DNA", 5);  parameterCount++;		}
		if (protein)			{  clearcutParameters[parameterCount] = new char[10];  *clearcutParameters[parameterCount] = '\0'; strncat(clearcutParameters[parameterCount], "--protein", 9);  parameterCount++;	}
		if (jukes)				{  clearcutParameters[parameterCount] = new char[8];  *clearcutParameters[parameterCount] = '\0'; strncat(clearcutParameters[parameterCount], "--jukes", 7);  parameterCount++;	}
		if (kimura)				{ clearcutParameters[parameterCount] = new char[9];  *clearcutParameters[parameterCount] = '\0'; strncat(clearcutParameters[parameterCount], "--kimura", 8);  parameterCount++;		}
		if (matrixout != "")	{  
			string tempMatrix =  "--matrixout=" + outputdir + matrixout;
			clearcutParameters[parameterCount] = new char[tempMatrix.length()+1];
			*clearcutParameters[parameterCount] = '\0'; strncat(clearcutParameters[parameterCount], tempMatrix.c_str(), tempMatrix.length());
			parameterCount++;
			outputNames.push_back((outputdir + matrixout));
			outputTypes["matrixout"].push_back((outputdir + matrixout));
		}

		if (ntrees != "1")		{  
			string tempNtrees = "--ntrees=" + ntrees; 
			clearcutParameters[parameterCount] = new char[tempNtrees.length()+1];
			*clearcutParameters[parameterCount] = '\0'; strncat(clearcutParameters[parameterCount], tempNtrees.c_str(), tempNtrees.length());
			parameterCount++;
		}

		if (expblen)			{ clearcutParameters[parameterCount] = new char[10];  *clearcutParameters[parameterCount] = '\0'; strncat(clearcutParameters[parameterCount], "--expblen", 9);  parameterCount++; 	}
		if (expdist)			{ clearcutParameters[parameterCount] = new char[10];  *clearcutParameters[parameterCount] = '\0'; strncat(clearcutParameters[parameterCount], "--expdist", 9);  parameterCount++;	}
        
        errno = 0;
		clearcut_main(numArgs, clearcutParameters); 
		
		//free memory
		for(int i = 0; i < numArgs; i++)  {  delete[] clearcutParameters[i];  }
		delete[] clearcutParameters; 
		
		if (!stdoutWanted) {	
			
			//set first tree file as new current treefile
			string currentTree = "";
			itTypes = outputTypes.find("tree");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentTree = (itTypes->second)[0]; current->setTreeFile(currentTree); }
			}
			
			m->mothurOut("\nOutput File Names:\n");
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




