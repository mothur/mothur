/*
 *  clearcutcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 5/11/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "clearcutcommand.h"

/**************************************************************************************/
ClearcutCommand::ClearcutCommand(string option)  {	
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
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
				
				it = parameters.find("phylip");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
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
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}
			
			string temp;
			temp = validParameter.validFile(parameters, "version", false);		if (temp == "not found"){	temp = "F";			}
			version = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "verbose", false);		if (temp == "not found"){	temp = "F";			}
			verbose = isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "quiet", false);		if (temp == "not found"){	temp = "F";			}
			quiet = isTrue(temp); 
			
			seed = validParameter.validFile(parameters, "seed", false);			if (seed == "not found"){	seed = "*";			}
			
			temp = validParameter.validFile(parameters, "norandom", false);		if (temp == "not found"){	temp = "F";			}
			norandom = isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "shuffle", false);		if (temp == "not found"){	temp = "F";			}
			shuffle = isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "neighbor", false);		if (temp == "not found"){	temp = "F";			}
			neighbor = isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "DNA", false);			if (temp == "not found"){	temp = "F";			}
			DNA = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "protein", false);		if (temp == "not found"){	temp = "F";			}
			protein = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "jukes", false);		if (temp == "not found"){	temp = "F";			}
			jukes = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "kimura", false);		if (temp == "not found"){	temp = "F";			}
			kimura = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "stdout", false);		if (temp == "not found"){	temp = "F";			}
			stdoutWanted = isTrue(temp); 
			
			matrixout = validParameter.validFile(parameters, "matrixout", false);	if (matrixout == "not found"){	matrixout = "";		}
			
			ntrees = validParameter.validFile(parameters, "ntrees", false);		if (ntrees == "not found"){	ntrees = "1";		}
			
			temp = validParameter.validFile(parameters, "expblen", false);		if (temp == "not found"){	temp = "F";			}
			expblen = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "expdist", false);		if (temp == "not found"){	temp = "F";			}
			expdist = isTrue(temp);
			
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
		m->mothurOut("The clearcut executable must be in a folder called clearcut in the same folder as your mothur executable, similar to mothur's requirements for using blast. \n");
		m->mothurOut("The clearcut command parameters are phylip, fasta, version, verbose, quiet, seed, norandom, shuffle, neighbor, expblen, expdist, ntrees, matrixout, stdout, kimura, jukes, protein, DNA. \n");
		m->mothurOut("The phylip parameter allows you to enter your phylip formatted distance matrix. \n");
		m->mothurOut("The fasta parameter allows you to enter your aligned fasta file, if you enter a fastafile you specify if the sequences are DNA or protein using the DNA or protein parameters. \n");
		
		m->mothurOut("The version parameter prints out the version of clearcut you are using, default=F. \n");
		m->mothurOut("The verbose parameter prints out more output from clearcut, default=F. \n");
		m->mothurOut("The quiet parameter turns on silent operation mode, default=F. \n");
		m->mothurOut("The seed parameter allows you to explicitly set the PRNG seed to a specific value. \n");
		m->mothurOut("The norandom parameter allows you to attempt joins deterministically, default=F. \n");
		m->mothurOut("The shuffle parameter allows you to randomly shuffle the distance matrix, default=F. \n");
		m->mothurOut("The neighbor parameter allows you to use traditional Neighbor-Joining algorithm, default=F. \n");
		
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
		
		if (abort == true) { return 0; }
				
		//prepare filename
		string outputName = outputDir + getRootName(getSimpleName(inputFile)) + "tre";
		
		//get location of clearcut
		GlobalData* globaldata = GlobalData::getInstance();
		string path = globaldata->argv;
		path = path.substr(0, (path.find_last_of('m')));
		
		string clearcutCommand = "";
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			clearcutCommand = path + "clearcut/clearcut ";
		#else
			clearcutCommand = path + "clearcut\\clearcut ";
		#endif
		
		//you gave us a distance matrix
		if (phylipfile != "") { clearcutCommand += "--distance "; 	}
		
		//you gave us a fastafile
		if (fastafile != "") { clearcutCommand += "--alignment "; 	}
		
		if (version)			{  clearcutCommand += "--version ";		}
		if (verbose)			{  clearcutCommand += "--verbose ";		}
		if (quiet)				{  clearcutCommand += "--quiet ";		}
		if (seed != "*")		{  clearcutCommand += "--seed=" + seed + " "; }
		if (norandom)			{  clearcutCommand += "--norandom ";	}
		if (shuffle)			{  clearcutCommand += "--shuffle ";		}
		if (neighbor)			{  clearcutCommand += "--neighbor ";	}
		
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			clearcutCommand += "--in=" + inputFile + " "; 
		#else
			clearcutCommand += "--in=\"" + inputFile + "\" "; 
		#endif
		
		if (stdoutWanted)		{  clearcutCommand += "--stdout ";		}
		else{  
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				clearcutCommand += "--out=" + outputName + " "; }
			#else
				clearcutCommand += "--out=\"" + outputName + "\" "; }
			#endif
		
		if (DNA)				{  clearcutCommand += "--DNA ";			}
		if (protein)			{  clearcutCommand += "--protein ";		}
		if (jukes)				{  clearcutCommand += "--jukes ";		}
		if (kimura)				{  clearcutCommand += "--kimura ";		}
		if (matrixout != "")	{  clearcutCommand += "--matrixout=" + matrixout + " ";  }
		if (ntrees != "1")		{  clearcutCommand += "--ntrees=" + ntrees + " "; }
		if (expblen)			{  clearcutCommand += "--expblen ";		}
		if (expdist)			{  clearcutCommand += "--expdist ";		}
	
		//run clearcut
		system(clearcutCommand.c_str());
		
		if (!stdoutWanted) {	
			m->mothurOutEndLine();
			m->mothurOut("Output File Names: "); m->mothurOutEndLine();
			m->mothurOut(outputName); m->mothurOutEndLine();
			if (matrixout != "")	{  m->mothurOut(matrixout); m->mothurOutEndLine();  }
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




