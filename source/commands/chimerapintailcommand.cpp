/*
 *  chimerapintailcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 4/1/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "chimerapintailcommand.h"
#include "pintail.h"


//**********************************************************************************************************************
vector<string> ChimeraPintailCommand::setParameters(){	
	try {
		CommandParameter ptemplate("reference", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(ptemplate);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","chimera-accnos",false,true,true); parameters.push_back(pfasta);
		CommandParameter pconservation("conservation", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pconservation);
		CommandParameter pquantile("quantile", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pquantile);
		CommandParameter pfilter("filter", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pfilter);
		CommandParameter pwindow("window", "Number", "", "0", "", "", "","","",false,false); parameters.push_back(pwindow);
		CommandParameter pincrement("increment", "Number", "", "25", "", "", "","",false,false); parameters.push_back(pincrement);
		CommandParameter pmask("mask", "String", "", "", "", "", "","",false,false); parameters.push_back(pmask);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);

		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPintailCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraPintailCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The chimera.pintail command reads a fastafile and referencefile and outputs potentially chimeric sequences.\n";
		helpString += "This command was created using the algorithms described in the 'At Least 1 in 20 16S rRNA Sequence Records Currently Held in the Public Repositories is Estimated To Contain Substantial Anomalies' paper by Kevin E. Ashelford 1, Nadia A. Chuzhanova 3, John C. Fry 1, Antonia J. Jones 2 and Andrew J. Weightman 1.\n";
		helpString += "The chimera.pintail command parameters are fasta, reference, filter, mask, processors, window, increment, conservation and quantile.\n";
		helpString += "The fasta parameter allows you to enter the fasta file containing your potentially chimeric sequences, and is required unless you have a valid current fasta file. \n";
		helpString += "You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amzon.fasta \n";
		helpString += "The reference parameter allows you to enter a reference file containing known non-chimeric sequences, and is required. \n";
		helpString += "The filter parameter allows you to specify if you would like to apply a vertical and 50% soft filter. \n";
		helpString += "The mask parameter allows you to specify a file containing one sequence you wish to use as a mask for the your sequences, by default no mask is applied.  You can apply an ecoli mask by typing, mask=default. \n";
		helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n";
		helpString += "The window parameter allows you to specify the window size for searching for chimeras, default=300. \n";
		helpString += "The increment parameter allows you to specify how far you move each window while finding chimeric sequences, default=25.\n";
		helpString += "The conservation parameter allows you to enter a frequency file containing the highest bases frequency at each place in the alignment.\n";
		helpString += "The quantile parameter allows you to enter a file containing quantiles for a template files sequences, if you use the filter the quantile file generated becomes unique to the fasta file you used.\n";
		helpString += "The chimera.pintail command should be in the following format: \n";
		helpString += "chimera.pintail(fasta=yourFastaFile, reference=yourTemplate) \n";
		helpString += "Example: chimera.pintail(fasta=AD.align, reference=silva.bacteria.fasta) \n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPintailCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraPintailCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "chimera") {  pattern = "[filename],[tag],pintail.chimeras-[filename],pintail.chimeras"; } 
        else if (type == "accnos") {  pattern = "[filename],[tag],pintail.accnos-[filename],pintail.accnos"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraPintailCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ChimeraPintailCommand::ChimeraPintailCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
		vector<string> tempOutNames;
		outputTypes["chimera"] = tempOutNames;
		outputTypes["accnos"] = tempOutNames;
        
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPintailCommand", "ChimeraPintailCommand");
		exit(1);
	}
}
//***************************************************************************************************************
ChimeraPintailCommand::ChimeraPintailCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter("chimera.pintail");
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			vector<string> tempOutNames;
			outputTypes["chimera"] = tempOutNames;
			outputTypes["accnos"] = tempOutNames;
        
		
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("reference");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["reference"] = inputDir + it->second;		}
				}
				
				it = parameters.find("conservation");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["conservation"] = inputDir + it->second;		}
				}
				
				it = parameters.find("quantile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["quantile"] = inputDir + it->second;		}
				}
			}

			
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
		
						bool ableToOpen;
						ifstream in;
						
						ableToOpen = m->openInputFile(fastaFileNames[i], in, "noerror");
					
						//if you can't open it, try default location
						if (!ableToOpen) {
							if (m->getDefaultPath() != "") { //default path is set
								string tryPath = m->getDefaultPath() + m->getSimpleName(fastaFileNames[i]);
								m->mothurOut("Unable to open " + fastaFileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								fastaFileNames[i] = tryPath;
							}
						}
						
						if (!ableToOpen) {
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
						
						if (!ableToOpen) { 
							m->mothurOut("Unable to open " + fastaFileNames[i] + ". It will be disregarded."); m->mothurOutEndLine(); 
							//erase from file list
							fastaFileNames.erase(fastaFileNames.begin()+i);
							i--;
						}else {
							m->setFastaFile(fastaFileNames[i]);
						}
					}
				}
				
				//make sure there is at least one valid file left
				if (fastaFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			string temp;
			temp = validParameter.validFile(parameters, "filter", false);			if (temp == "not found") { temp = "F"; }
			filter = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
			
			temp = validParameter.validFile(parameters, "window", false);			if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, window);
			
			temp = validParameter.validFile(parameters, "increment", false);		if (temp == "not found") { temp = "25"; }
			m->mothurConvert(temp, increment);
			
			
			//this has to go after save so that if the user sets save=t and provides no reference we abort
			templatefile = validParameter.validFile(parameters, "reference", true);
			if (templatefile == "not found") { m->mothurOut("[ERROR]: The reference parameter is a required, aborting.\n"); abort = true;
			}else if (templatefile == "not open") { abort = true; }
			
			
			maskfile = validParameter.validFile(parameters, "mask", false);
			if (maskfile == "not found") { maskfile = "";  }	
			else if (maskfile != "default")  { 
				if (inputDir != "") {
					string path = m->hasPath(maskfile);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	maskfile = inputDir + maskfile;		}
				}

				ifstream in;
				bool	ableToOpen = m->openInputFile(maskfile, in, "no error");
				if (!ableToOpen) { 
					if (m->getDefaultPath() != "") { //default path is set
							string tryPath = m->getDefaultPath() + m->getSimpleName(maskfile);
							m->mothurOut("Unable to open " + maskfile + ". Trying default " + tryPath); m->mothurOutEndLine();
							ifstream in2;
							ableToOpen = m->openInputFile(tryPath, in2, "noerror");
							in2.close();
							maskfile = tryPath;
					}
				}
				
				if (!ableToOpen) {
						if (m->getOutputDir() != "") { //default path is set
							string tryPath = m->getOutputDir() + m->getSimpleName(maskfile);
							m->mothurOut("Unable to open " + maskfile + ". Trying output directory " + tryPath); m->mothurOutEndLine();
							ifstream in2;
							ableToOpen = m->openInputFile(tryPath, in2, "noerror");
							in2.close();
							maskfile = tryPath;
						}
				}
				
				in.close();
					
				if (!ableToOpen) { 
						m->mothurOut("Unable to open " + maskfile + "."); m->mothurOutEndLine(); 
						abort = true;
				}
			}

			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}
		
			consfile = validParameter.validFile(parameters, "conservation", true);
			if (consfile == "not open") { abort = true; }
			else if (consfile == "not found") { 
				consfile = "";  
				//check for consfile
				string tempConsFile = m->getRootName(inputDir + m->getSimpleName(templatefile)) + "freq";
				ifstream FileTest(tempConsFile.c_str());
				if(FileTest){	
					bool GoodFile = m->checkReleaseVersion(FileTest, m->getVersion());
					if (GoodFile) {  
						m->mothurOut("I found " + tempConsFile + " in your input file directory. I will use it to save time."); m->mothurOutEndLine();  consfile = tempConsFile;  FileTest.close();	
					}
				}else {
					string tempConsFile = m->getDefaultPath() + m->getRootName(m->getSimpleName(templatefile)) + "freq";
					ifstream FileTest2(tempConsFile.c_str());
					if(FileTest2){	
						bool GoodFile = m->checkReleaseVersion(FileTest2, m->getVersion());
						if (GoodFile) {  
							m->mothurOut("I found " + tempConsFile + " in your input file directory. I will use it to save time."); m->mothurOutEndLine();  consfile = tempConsFile;  FileTest2.close();	
						}
					}
				}
			}	
			
			quanfile = validParameter.validFile(parameters, "quantile", true);
			if (quanfile == "not open") { abort = true; }
			else if (quanfile == "not found") { quanfile = ""; }
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPintailCommand", "ChimeraPintailCommand");
		exit(1);
	}
}
//***************************************************************************************************************

int ChimeraPintailCommand::execute(){
	try{
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		for (int s = 0; s < fastaFileNames.size(); s++) {
				
			m->mothurOut("Checking sequences from " + fastaFileNames[s] + " ..." ); m->mothurOutEndLine();

			int start = time(NULL);	
			
			//set user options
			if (maskfile == "default") { m->mothurOut("I am using the default 236627 EU009184.1 Shigella dysenteriae str. FBD013."); m->mothurOutEndLine();  }
			
			//check for quantile to save the time
			string baseName = templatefile;
			
			string tempQuan = "";
			if ((!filter) && (maskfile == "")) {
				tempQuan = inputDir + m->getRootName(m->getSimpleName(baseName)) + "pintail.quan";
			}else if ((!filter) && (maskfile != "")) { 
				tempQuan = inputDir + m->getRootName(m->getSimpleName(baseName)) + "pintail.masked.quan";
			}else if ((filter) && (maskfile != "")) { 
				tempQuan = inputDir + m->getRootName(m->getSimpleName(baseName)) + "pintail.filtered." + m->getSimpleName(m->getRootName(fastaFileNames[s])) + "masked.quan";
			}else if ((filter) && (maskfile == "")) { 
				tempQuan = inputDir + m->getRootName(m->getSimpleName(baseName)) + "pintail.filtered." + m->getSimpleName(m->getRootName(fastaFileNames[s])) + "quan";
			}
			
			ifstream FileTest(tempQuan.c_str());
			if(FileTest){	
				bool GoodFile = m->checkReleaseVersion(FileTest, m->getVersion());
				if (GoodFile) {  
					m->mothurOut("I found " + tempQuan + " in your input file directory. I will use it to save time."); m->mothurOutEndLine();  quanfile = tempQuan;  FileTest.close();	
				}
			}else {
				string tryPath = m->getDefaultPath();
				string tempQuan = "";
				if ((!filter) && (maskfile == "")) {
					tempQuan = tryPath + m->getRootName(m->getSimpleName(baseName)) + "pintail.quan";
				}else if ((!filter) && (maskfile != "")) { 
					tempQuan = tryPath + m->getRootName(m->getSimpleName(baseName)) + "pintail.masked.quan";
				}else if ((filter) && (maskfile != "")) { 
					tempQuan = tryPath + m->getRootName(m->getSimpleName(baseName)) + "pintail.filtered." + m->getSimpleName(m->getRootName(fastaFileNames[s])) + "masked.quan";
				}else if ((filter) && (maskfile == "")) { 
					tempQuan = tryPath + m->getRootName(m->getSimpleName(baseName)) + "pintail.filtered." + m->getSimpleName(m->getRootName(fastaFileNames[s])) + "quan";
				}
				
				ifstream FileTest2(tempQuan.c_str());
				if(FileTest2){	
					bool GoodFile = m->checkReleaseVersion(FileTest2, m->getVersion());
					if (GoodFile) {  
						m->mothurOut("I found " + tempQuan + " in your input file directory. I will use it to save time."); m->mothurOutEndLine();  quanfile = tempQuan;  FileTest2.close();	
					}
				}
			}
			chimera = new Pintail(fastaFileNames[s], templatefile, filter, processors, maskfile, consfile, quanfile, window, increment, outputDir);
			
			if (outputDir == "") { outputDir = m->hasPath(fastaFileNames[s]);  }//if user entered a file with a path then preserve it
			string outputFileName, accnosFileName;
            map<string, string> variables;
            variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastaFileNames[s]));
			if (maskfile != "") { variables["[tag]"] = m->getSimpleName(m->getRootName(maskfile)); }
            outputFileName = getOutputFileName("chimera", variables);
            accnosFileName = getOutputFileName("accnos", variables);
			
			
			if (m->getControl_pressed()) { delete chimera; for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	}  return 0;	}
			
			if (chimera->getUnaligned()) { 
				m->mothurOut("Your template sequences are different lengths, please correct."); m->mothurOutEndLine(); 
				delete chimera;
				return 0; 
			}
			templateSeqsLength = chimera->getLength();
		
								
			//break up file
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				vector<unsigned long long> positions = m->divideFile(fastaFileNames[s], processors);
			
				for (int i = 0; i < (positions.size()-1); i++) {
					lines.push_back(new linePair(positions[i], positions[(i+1)]));
				}	
			
				if(processors == 1){
		
					numSeqs = driver(lines[0], outputFileName, fastaFileNames[s], accnosFileName);
					
					if (m->getControl_pressed()) { outputTypes.clear(); m->mothurRemove(outputFileName); m->mothurRemove(accnosFileName); for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	} for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear(); delete chimera; return 0; }
					
				}else{
					processIDS.resize(0);
					
					numSeqs = createProcesses(outputFileName, fastaFileNames[s], accnosFileName); 
				
					rename((outputFileName + toString(processIDS[0]) + ".temp").c_str(), outputFileName.c_str());
					rename((accnosFileName + toString(processIDS[0]) + ".temp").c_str(), accnosFileName.c_str());
						
					//append output files
					for(int i=1;i<processors;i++){
						m->appendFiles((outputFileName + toString(processIDS[i]) + ".temp"), outputFileName);
						m->mothurRemove((outputFileName + toString(processIDS[i]) + ".temp"));
					}
					
					//append output files
					for(int i=1;i<processors;i++){
						m->appendFiles((accnosFileName + toString(processIDS[i]) + ".temp"), accnosFileName);
						m->mothurRemove((accnosFileName + toString(processIDS[i]) + ".temp"));
					}
										
					if (m->getControl_pressed()) { 
						m->mothurRemove(outputFileName); 
						m->mothurRemove(accnosFileName);
						for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	} outputTypes.clear();
						for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
						delete chimera;
						return 0;
					}
				}

			#else
				lines.push_back(new linePair(0, 1000));
				numSeqs = driver(lines[0], outputFileName, fastaFileNames[s], accnosFileName);
				
				if (m->getControl_pressed()) { outputTypes.clear(); m->mothurRemove(outputFileName); m->mothurRemove(accnosFileName); for (int j = 0; j < outputNames.size(); j++) {	m->mothurRemove(outputNames[j]);	} for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear(); delete chimera; return 0; }
			#endif
			
			delete chimera;
			for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
			
			outputNames.push_back(outputFileName); outputTypes["chimera"].push_back(outputFileName);
			outputNames.push_back(accnosFileName); outputTypes["accnos"].push_back(accnosFileName);
			
			m->mothurOutEndLine();
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
		m->errorOut(e, "ChimeraPintailCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int ChimeraPintailCommand::driver(linePair* filePos, string outputFName, string filename, string accnos){
	try {
		ofstream out;
		m->openOutputFile(outputFName, out);
		
		ofstream out2;
		m->openOutputFile(accnos, out2);
		
		ifstream inFASTA;
		m->openInputFile(filename, inFASTA);

		inFASTA.seekg(filePos->start);

		bool done = false;
		int count = 0;
	
		while (!done) {
				
			if (m->getControl_pressed()) {	return 1;	}
		
			Sequence* candidateSeq = new Sequence(inFASTA);  m->gobble(inFASTA);
				
			if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
				
				if (candidateSeq->getAligned().length() != templateSeqsLength)  {  //chimeracheck does not require seqs to be aligned
					m->mothurOut(candidateSeq->getName() + " is not the same length as the template sequences. Skipping."); m->mothurOutEndLine();
				}else{
					//find chimeras
					chimera->getChimeras(candidateSeq);
					
					if (m->getControl_pressed()) {	delete candidateSeq; return 1;	}
		
					//print results
					chimera->print(out, out2);
				}
				count++;
			}
			delete candidateSeq;
			
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				unsigned long long pos = inFASTA.tellg();
				if ((pos == -1) || (pos >= filePos->end)) { break; }
			#else
				if (inFASTA.eof()) { break; }
			#endif
			
			//report progress
			if((count) % 100 == 0){	m->mothurOutJustToScreen("Processing sequence: " + toString(count) + "\n"); 		}
		}
		//report progress
		if((count) % 100 != 0){	m->mothurOutJustToScreen("Processing sequence: " + toString(count) + "\n"); 		}
		
		out.close();
		out2.close();
		inFASTA.close();
				
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPintailCommand", "driver");
		exit(1);
	}
}
/**************************************************************************************************/

int ChimeraPintailCommand::createProcesses(string outputFileName, string filename, string accnos) {
	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		int process = 0;
		int num = 0;
        bool recalc = false;
		
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driver(lines[process], outputFileName + toString(m->mothurGetpid(process)) + ".temp", filename, accnos + toString(m->mothurGetpid(process)) + ".temp");
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = outputFileName + toString(m->mothurGetpid(process)) + ".num.temp";
				m->openOutputFile(tempFile, out);
				out << num << endl;
				out.close();

				exit(0);
            }else {
                m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(process) + "\n"); processors = process;
                for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                //wait to die
                for (int i=0;i<processIDS.size();i++) {
                    int temp = processIDS[i];
                    wait(&temp);
                }
                m->setControl_pressed(false);
                recalc = true;
                break;
			}
		}

        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->setControl_pressed(false);
					 processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");

            for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
            
            vector<unsigned long long> positions = m->divideFile(filename, processors);
            for (int i = 0; i < (positions.size()-1); i++) {
                lines.push_back(new linePair(positions[i], positions[(i+1)]));
            }
            num = 0;
            processIDS.resize(0);
            process = 0;
            
            while (process != processors) {
                pid_t pid = fork();
                
                if (pid > 0) {
                    processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                    process++;
                }else if (pid == 0){
                    num = driver(lines[process], outputFileName + toString(m->mothurGetpid(process)) + ".temp", filename, accnos + toString(m->mothurGetpid(process)) + ".temp");
                    
                    //pass numSeqs to parent
                    ofstream out;
                    string tempFile = outputFileName + toString(m->mothurGetpid(process)) + ".num.temp";
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
        }

		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		for (int i = 0; i < processIDS.size(); i++) {
			ifstream in;
			string tempFile =  outputFileName + toString(processIDS[i]) + ".num.temp";
			m->openInputFile(tempFile, in);
			if (!in.eof()) { int tempNum = 0; in >> tempNum; num += tempNum; }
			in.close(); m->mothurRemove(tempFile);
		}
		
		return num;
#endif		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPintailCommand", "createProcesses");
		exit(1);
	}
}

/**************************************************************************************************/


