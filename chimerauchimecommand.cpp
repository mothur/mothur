/*
 *  chimerauchimecommand.cpp
 *  Mothur
 *
 *  Created by westcott on 5/13/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "chimerauchimecommand.h"
#include "deconvolutecommand.h"
#include "uc.h"
#include "sequence.hpp"


//**********************************************************************************************************************
vector<string> ChimeraUchimeCommand::setParameters(){	
	try {
		CommandParameter ptemplate("reference", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(ptemplate);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pname);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "",false,false); parameters.push_back(pprocessors);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraUchimeCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraUchimeCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The chimera.uchime command reads a fastafile and referencefile and outputs potentially chimeric sequences.\n";
		helpString += "This command is a wrapper for uchime written by Robert C. Edgar.\n";
		helpString += "The chimera.uchime command parameters are fasta, name, reference and processors.\n";
		helpString += "The fasta parameter allows you to enter the fasta file containing your potentially chimeric sequences, and is required, unless you have a valid current fasta file. \n";
		helpString += "The name parameter allows you to provide a name file, if you are using template=self. \n";
		helpString += "You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amazon.fasta \n";
		helpString += "The reference parameter allows you to enter a reference file containing known non-chimeric sequences, and is required. You may also set template=self, in this case the abundant sequences will be used as potential parents. \n";
		helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n";
#ifdef USE_MPI
		helpString += "When using MPI, the processors parameter is set to the number of MPI processes running. \n";
#endif
		helpString += "The chimera.uchime command should be in the following format: \n";
		helpString += "chimera.uchime(fasta=yourFastaFile, reference=yourTemplate) \n";
		helpString += "Example: chimera.uchime(fasta=AD.align, reference=silva.gold.align) \n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraUchimeCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
ChimeraUchimeCommand::ChimeraUchimeCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
		vector<string> tempOutNames;
		outputTypes["chimera"] = tempOutNames;
		outputTypes["accnos"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraUchimeCommand", "ChimeraUchimeCommand");
		exit(1);
	}
}
//***************************************************************************************************************
ChimeraUchimeCommand::ChimeraUchimeCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter("chimera.uchime");
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			vector<string> tempOutNames;
			outputTypes["chimera"] = tempOutNames;
			outputTypes["accnos"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", false);
			if (fastafile == "not found") { 				
				//if there is a current fasta file, use it
				string filename = m->getFastaFile(); 
				if (filename != "") { fastaFileNames.push_back(filename); m->mothurOut("Using " + filename + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { 
				m->splitAtDash(fastafile, fastaFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastaFileNames.size(); i++) {
					
					bool ignore = false;
					if (fastaFileNames[i] == "current") { 
						fastaFileNames[i] = m->getFastaFile(); 
						if (fastaFileNames[i] != "") {  m->mothurOut("Using " + fastaFileNames[i] + " as input file for the fasta parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current fastafile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							fastaFileNames.erase(fastaFileNames.begin()+i);
							i--;
						}
					}
					
					if (!ignore) {
						
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
				}
				
				//make sure there is at least one valid file left
				if (fastaFileNames.size() == 0) { m->mothurOut("[ERROR]: no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			
			//check for required parameters
			bool hasName = true;
			namefile = validParameter.validFile(parameters, "name", false);
			if (namefile == "not found") { namefile = "";  hasName = false; }
			else { 
				m->splitAtDash(namefile, nameFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < nameFileNames.size(); i++) {
					
					bool ignore = false;
					if (nameFileNames[i] == "current") { 
						nameFileNames[i] = m->getNameFile(); 
						if (nameFileNames[i] != "") {  m->mothurOut("Using " + nameFileNames[i] + " as input file for the name parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current namefile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							nameFileNames.erase(nameFileNames.begin()+i);
							i--;
						}
					}
					
					if (!ignore) {
						
						if (inputDir != "") {
							string path = m->hasPath(nameFileNames[i]);
							//if the user has not given a path then, add inputdir. else leave path alone.
							if (path == "") {	nameFileNames[i] = inputDir + nameFileNames[i];		}
						}
						
						int ableToOpen;
						ifstream in;
						
						ableToOpen = m->openInputFile(nameFileNames[i], in, "noerror");
						
						//if you can't open it, try default location
						if (ableToOpen == 1) {
							if (m->getDefaultPath() != "") { //default path is set
								string tryPath = m->getDefaultPath() + m->getSimpleName(nameFileNames[i]);
								m->mothurOut("Unable to open " + nameFileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								nameFileNames[i] = tryPath;
							}
						}
						
						if (ableToOpen == 1) {
							if (m->getOutputDir() != "") { //default path is set
								string tryPath = m->getOutputDir() + m->getSimpleName(nameFileNames[i]);
								m->mothurOut("Unable to open " + nameFileNames[i] + ". Trying output directory " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								nameFileNames[i] = tryPath;
							}
						}
						
						in.close();
						
						if (ableToOpen == 1) { 
							m->mothurOut("Unable to open " + nameFileNames[i] + ". It will be disregarded."); m->mothurOutEndLine(); 
							//erase from file list
							nameFileNames.erase(nameFileNames.begin()+i);
							i--;
						}
					}
				}
				
				//make sure there is at least one valid file left
				if (nameFileNames.size() == 0) { m->mothurOut("[ERROR]: no valid name files."); m->mothurOutEndLine(); abort = true; }
			}
			
			if (hasName && (nameFileNames.size() != fastaFileNames.size())) { m->mothurOut("[ERROR]: The number of namefiles does not match the number of fastafiles, please correct."); m->mothurOutEndLine(); abort=true; }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}
			
			
			string path;
			it = parameters.find("reference");
			//user has given a template file
			if(it != parameters.end()){ 
				if (it->second == "self") { templatefile = "self"; }
				else {
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["reference"] = inputDir + it->second;		}
					
					templatefile = validParameter.validFile(parameters, "reference", true);
					if (templatefile == "not open") { abort = true; }
					else if (templatefile == "not found") { templatefile = "";  m->mothurOut("reference is a required parameter for the chimera.slayer command."); m->mothurOutEndLine(); abort = true;  }	
				}
			}
			
			string temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			convert(temp, processors);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayerCommand", "ChimeraSlayerCommand");
		exit(1);
	}
}
//***************************************************************************************************************

int ChimeraUchimeCommand::execute(){
	try{
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		for (int s = 0; s < fastaFileNames.size(); s++) {
			
			m->mothurOut("Checking sequences from " + fastaFileNames[s] + " ..." ); m->mothurOutEndLine();
			
			int start = time(NULL);	
			string nameFile = "";
			
			if (templatefile == "self") { //you want to run slayer with a refernce template
				
				#ifdef USE_MPI	
					int pid; 
					MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
					if (pid == 0) { //you are the root process 
				#endif	
				
				if (processors != 1) { m->mothurOut("When using template=self, mothur can only use 1 processor, continuing."); m->mothurOutEndLine(); processors = 1; }
				if (nameFileNames.size() != 0) { //you provided a namefile and we don't need to create one
					nameFile = nameFileNames[s];
				}else {
					m->mothurOutEndLine(); m->mothurOut("No namesfile given, running unique.seqs command to generate one."); m->mothurOutEndLine(); m->mothurOutEndLine();
					
					//use unique.seqs to create new name and fastafile
					string inputString = "fasta=" + fastaFileNames[s];
					m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
					m->mothurOut("Running command: unique.seqs(" + inputString + ")"); m->mothurOutEndLine(); 
					
					Command* uniqueCommand = new DeconvoluteCommand(inputString);
					uniqueCommand->execute();
					
					map<string, vector<string> > filenames = uniqueCommand->getOutputFiles();
					
					delete uniqueCommand;
					
					m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
					
					nameFile = filenames["name"][0];
					fastaFileNames[s] = filenames["fasta"][0];
				}
				
				//create input file for uchime
				//read through fastafile and store info
				map<string, string> seqs;
				ifstream in;
				m->openInputFile(fastaFileNames[s], in);
				
				while (!in.eof()) {
					
					if (m->control_pressed) { in.close(); for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	}  return 0; }
					
					Sequence seq(in); m->gobble(in);
					seqs[seq.getName()] = seq.getAligned();
				}
				in.close();
				
				//read namefile
				vector<seqPriorityNode> nameMapCount;
				int error = m->readNames(nameFile, nameMapCount, seqs);
				
				if (m->control_pressed) { for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	}  return 0; }
				
				if (error == 1) { for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	}  return 0; }
				if (seqs.size() != nameMapCount.size()) { m->mothurOut( "The number of sequences in your fastafile does not match the number of sequences in your namefile, aborting."); m->mothurOutEndLine(); for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	}  return 0; }
				
				sort(nameMapCount.begin(), nameMapCount.end(), compareSeqPriorityNodes);
				
				string newFasta = fastaFileNames[s] + ".temp";
				ofstream out;
				m->openOutputFile(newFasta, out);
				
				//print new file in order of
				for (int i = 0; i < nameMapCount.size(); i++) {
					out << ">" << nameMapCount[i].name  << "/ab=" << nameMapCount[i].numIdentical << "/" << endl << nameMapCount[i].seq << endl;
				}
				out.close();
				
				fastaFileNames[s] = newFasta;
						
				#ifdef USE_MPI	
					}
				#endif
				if (m->control_pressed) {  for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	}  return 0;	}				
			}
			
			if (outputDir == "") { outputDir = m->hasPath(fastaFileNames[s]);  }//if user entered a file with a path then preserve it				
			string outputFileName = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + "slayer.chimera";
			string accnosFileName = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s]))  + "slayer.accnos";
			
			if (m->control_pressed) {  for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	}  return 0;	}
			
			int numSeqs = 0;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			if(processors == 1){ numSeqs = driver(outputFileName, fastaFileNames[s], accnosFileName); }
			else{	numSeqs = createProcesses(outputFileName, fastaFileNames[s], accnosFileName); }
#else
			numSeqs = driver(outputFileName, fastaFileNames[s], accnosFileName);
#endif
			if (m->control_pressed) { for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	} return 0; }

			
			outputNames.push_back(outputFileName); outputTypes["chimera"].push_back(outputFileName);
			outputNames.push_back(accnosFileName); outputTypes["accnos"].push_back(accnosFileName);
			
			m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences.");	m->mothurOutEndLine();
		}
		
		//set accnos file as new current accnosfile
		string current = "";
		itTypes = outputTypes.find("accnos");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setAccnosFile(current); }
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}	
		m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraUchimeCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int ChimeraUchimeCommand::driver(string outputFName, string filename, string accnos){
	try {
		
		vector<char*> cPara;
		
		char* tempUchime = new char[8];  
		strcpy(tempUchime, "./uchime "); 
		cPara.push_back(tempUchime);
		
		char* tempIn = new char[7];  
		strcpy(tempIn, "--input"); 
		cPara.push_back(tempIn);
		char* temp = new char[filename.length()];
		strcpy(temp, filename.c_str());
		cPara.push_back(temp);
		
		//are you using a reference file
		if (templatefile != "self") {
						
			//add reference file
			char* tempRef = new char[4]; 
			strcpy(tempRef, "--db"); 
			cPara.push_back(tempRef);  
			char* tempR = new char[templatefile.length()];
			strcpy(tempR, templatefile.c_str());
			cPara.push_back(tempR);
		}
		
		char* tempO = new char[11]; 
		strcpy(tempO, "--uchimeout"); 
		cPara.push_back(tempO);
		char* tempout = new char[outputFName.length()];
		strcpy(tempout, outputFName.c_str());
		cPara.push_back(tempout);
		
		char** uchimeParameters;
		uchimeParameters = new char*[cPara.size()];
		for (int i = 0; i < cPara.size(); i++) {  uchimeParameters[i] = cPara[i];  } 
		int numArgs = cPara.size();
		
		uchime_main(numArgs, uchimeParameters); 
		
		//free memory
		for(int i = 0; i < cPara.size(); i++)  {  delete[] cPara[i];  }
		delete[] uchimeParameters; 
		
		//create accnos file from uchime results
		ifstream in; 
		m->openInputFile(outputFName, in);
		
		ofstream out;
		m->openOutputFile(accnos, out);
		
		int num = 0;
		while(!in.eof()) {
			
			if (m->control_pressed) { break; }
			
			string name = "";
			string chimeraFlag = "";
			in >> chimeraFlag >> name;
			
			//fix name if needed
			if (templatefile != "self") { 
				name = name.substr(0, name.length()-1); //rip off last /
				name = name.substr(0, name.find_last_of('/'));
			}
			
			for (int i = 0; i < 15; i++) {  in >> chimeraFlag; }
			m->gobble(in);
			
			if (chimeraFlag == "Y") {  out << name << endl; }
			num++;
		}
		in.close();
		out.close();
		
		return num;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraUchimeCommand", "driver");
		exit(1);
	}
}
/**************************************************************************************************/

int ChimeraUchimeCommand::createProcesses(string outputFileName, string filename, string accnos) {
	try {
		
		processIDS.clear();
		int process = 1;
		int num = 0;
		
		//break up file into multiple files
		vector<string> files;
		m->divideFile(filename, processors, files);
		
		if (m->control_pressed) {  return 0;  }
		
#ifdef USE_MPI	
		int pid, numSeqsPerProcessor; 
		int tag = 2001;
		
		MPI_Status status; 
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
		MPI_Comm_size(MPI_COMM_WORLD, &processors); 
				
		if (pid == 0) { //you are the root process 
			num = driver(outputFileName, files[0], accnos);
			
			if (templatefile != "self") {
				//wait on chidren
				for(int j = 1; j < processors; j++) { 
					int temp;
					MPI_Recv(&temp, 1, MPI_INT, j, tag, MPI_COMM_WORLD, &status);
					num += temp;
					
					m->appendFiles((outputFileName + toString(j) + ".temp"), outputFileName);
					remove((outputFileName + toString(j) + ".temp").c_str());
					
					m->appendFiles((accnos + toString(j) + ".temp"), accnos);
					remove((accnos + toString(j) + ".temp").c_str());
				}
			}
		}else{ //you are a child process
			if (templatefile != "self") { //if template=self we can only use 1 processor
				num = driver(outputFileName+toString(pid) + ".temp", files[pid], accnos+toString(pid) + ".temp");	
				
				//send numSeqs to parent
				MPI_Send(&num, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
			}
		}

		MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
#else
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driver(outputFileName + toString(getpid()) + ".temp", files[process], accnos + toString(getpid()) + ".temp");
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = outputFileName + toString(getpid()) + ".num.temp";
				m->openOutputFile(tempFile, out);
				out << num << endl;
				out.close();
				
				exit(0);
			}else { 
				m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
				for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
				exit(0);
			}
		}
		
		//do my part
		num = driver(outputFileName, files[0], accnos);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		for (int i = 0; i < processIDS.size(); i++) {
			ifstream in;
			string tempFile =  outputFileName + toString(processIDS[i]) + ".num.temp";
			m->openInputFile(tempFile, in);
			if (!in.eof()) { int tempNum = 0; in >> tempNum; num += tempNum; }
			in.close(); remove(tempFile.c_str());
		}
		
		
		//append output files
		for(int i=0;i<processIDS[i];i++){
			m->appendFiles((outputFileName + toString(processIDS[i]) + ".temp"), outputFileName);
			remove((outputFileName + toString(processIDS[i]) + ".temp").c_str());
			
			m->appendFiles((accnos + toString(processIDS[i]) + ".temp"), accnos);
			remove((accnos + toString(processIDS[i]) + ".temp").c_str());
		}
#endif		
		//get rid of the file pieces.
		for (int i = 0; i < files.size(); i++) { remove(files[i].c_str()); }
		
		return num;	
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraUchimeCommand", "createProcesses");
		exit(1);
	}
}

/**************************************************************************************************/

