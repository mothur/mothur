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
			inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("reference");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["reference"] = inputDir + it->second;		}
				}
				
				it = parameters.find("conservation");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["conservation"] = inputDir + it->second;		}
				}
				
				it = parameters.find("quantile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["quantile"] = inputDir + it->second;		}
				}
			}

			
			//check for required parameters
			fastafile = validParameter.valid(parameters, "fasta");
			if (fastafile == "not found") { 				
				//if there is a current fasta file, use it
				string filename = current->getFastaFile(); 
				if (filename != "") { fastaFileNames.push_back(filename); m->mothurOut("Using " + filename + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { 
				util.splitAtDash(fastafile, fastaFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastaFileNames.size(); i++) {
					
					bool ignore = false;
					if (fastaFileNames[i] == "current") { 
						fastaFileNames[i] = current->getFastaFile(); 
						if (fastaFileNames[i] != "") {  m->mothurOut("Using " + fastaFileNames[i] + " as input file for the fasta parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current fastafile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							fastaFileNames.erase(fastaFileNames.begin()+i);
							i--;
						}
					}
					
                    if (!ignore) {
                        if (util.checkLocations(fastaFileNames[i], current->getLocations())) { current->setFastaFile(fastaFileNames[i]); }
                        else { fastaFileNames.erase(fastaFileNames.begin()+i); i--; } //erase from file list
                    }
                }
				
				//make sure there is at least one valid file left
				if (fastaFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			string temp;
			temp = validParameter.valid(parameters, "filter");			if (temp == "not found") { temp = "F"; }
			filter = util.isTrue(temp);
						
			temp = validParameter.valid(parameters, "window");			if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, window);
			
			temp = validParameter.valid(parameters, "increment");		if (temp == "not found") { temp = "25"; }
			util.mothurConvert(temp, increment);
			
			
			//this has to go after save so that if the user sets save=t and provides no reference we abort
			templatefile = validParameter.validFile(parameters, "reference");
			if (templatefile == "not found") { m->mothurOut("[ERROR]: The reference parameter is a required, aborting.\n"); abort = true;
			}else if (templatefile == "not open") { abort = true; }
			
			
			maskfile = validParameter.valid(parameters, "mask");
			if (maskfile == "not found") { maskfile = "";  }	
			else if (maskfile != "default")  {
                if (util.checkLocations(maskfile, current->getLocations())) {  }
                else { m->mothurOut("Unable to open " + maskfile + "."); m->mothurOutEndLine(); abort = true; } //erase from file list
			}

			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";	}
		
			consfile = validParameter.validFile(parameters, "conservation");
			if (consfile == "not open") { abort = true; }
			else if (consfile == "not found") { 
				consfile = "";  
				//check for consfile
				string tempConsFile = util.getRootName(inputDir + util.getSimpleName(templatefile)) + "freq";
				ifstream FileTest(tempConsFile.c_str());
				if(FileTest){	
					bool GoodFile = util.checkReleaseVersion(FileTest, current->getVersion());
					if (GoodFile) {  
						m->mothurOut("I found " + tempConsFile + " in your input file directory. I will use it to save time."); m->mothurOutEndLine();  consfile = tempConsFile;  FileTest.close();	
					}
				}else {
					string tempConsFile = current->getDefaultPath() + util.getRootName(util.getSimpleName(templatefile)) + "freq";
					ifstream FileTest2(tempConsFile.c_str());
					if(FileTest2){	
						bool GoodFile = util.checkReleaseVersion(FileTest2, current->getVersion());
						if (GoodFile) {  
							m->mothurOut("I found " + tempConsFile + " in your input file directory. I will use it to save time."); m->mothurOutEndLine();  consfile = tempConsFile;  FileTest2.close();	
						}
					}
				}
			}	
			
			quanfile = validParameter.validFile(parameters, "quantile");
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
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		for (int s = 0; s < fastaFileNames.size(); s++) {
				
			m->mothurOut("Checking sequences from " + fastaFileNames[s] + " ..." ); m->mothurOutEndLine();

			long start = time(NULL);	
			
			//set user options
			if (maskfile == "default") { m->mothurOut("I am using the default 236627 EU009184.1 Shigella dysenteriae str. FBD013."); m->mothurOutEndLine();  }
			
			//check for quantile to save the time
			string baseName = templatefile;
			
			string tempQuan = "";
			if ((!filter) && (maskfile == "")) {
				tempQuan = inputDir + util.getRootName(util.getSimpleName(baseName)) + "pintail.quan";
			}else if ((!filter) && (maskfile != "")) { 
				tempQuan = inputDir + util.getRootName(util.getSimpleName(baseName)) + "pintail.masked.quan";
			}else if ((filter) && (maskfile != "")) { 
				tempQuan = inputDir + util.getRootName(util.getSimpleName(baseName)) + "pintail.filtered." + util.getSimpleName(util.getRootName(fastaFileNames[s])) + "masked.quan";
			}else if ((filter) && (maskfile == "")) { 
				tempQuan = inputDir + util.getRootName(util.getSimpleName(baseName)) + "pintail.filtered." + util.getSimpleName(util.getRootName(fastaFileNames[s])) + "quan";
			}
			
			ifstream FileTest(tempQuan.c_str());
			if(FileTest){	
				bool GoodFile = util.checkReleaseVersion(FileTest, current->getVersion());
				if (GoodFile) {  
					m->mothurOut("I found " + tempQuan + " in your input file directory. I will use it to save time."); m->mothurOutEndLine();  quanfile = tempQuan;  FileTest.close();	
				}
			}else {
				string tryPath = current->getDefaultPath();
				string tempQuan = "";
				if ((!filter) && (maskfile == "")) {
					tempQuan = tryPath + util.getRootName(util.getSimpleName(baseName)) + "pintail.quan";
				}else if ((!filter) && (maskfile != "")) { 
					tempQuan = tryPath + util.getRootName(util.getSimpleName(baseName)) + "pintail.masked.quan";
				}else if ((filter) && (maskfile != "")) { 
					tempQuan = tryPath + util.getRootName(util.getSimpleName(baseName)) + "pintail.filtered." + util.getSimpleName(util.getRootName(fastaFileNames[s])) + "masked.quan";
				}else if ((filter) && (maskfile == "")) { 
					tempQuan = tryPath + util.getRootName(util.getSimpleName(baseName)) + "pintail.filtered." + util.getSimpleName(util.getRootName(fastaFileNames[s])) + "quan";
				}
				
				ifstream FileTest2(tempQuan.c_str());
				if(FileTest2){	
					bool GoodFile = util.checkReleaseVersion(FileTest2, current->getVersion());
					if (GoodFile) {  
						m->mothurOut("I found " + tempQuan + " in your input file directory. I will use it to save time."); m->mothurOutEndLine();  quanfile = tempQuan;  FileTest2.close();	
					}
				}
			}
			chimera = new Pintail(fastaFileNames[s], templatefile, filter, maskfile, consfile, quanfile, window, increment, outputDir, current->getVersion());
			
			if (outputDir == "") { outputDir = util.hasPath(fastaFileNames[s]);  }//if user entered a file with a path then preserve it
			string outputFileName, accnosFileName;
            map<string, string> variables;
            variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastaFileNames[s]));
			if (maskfile != "") { variables["[tag]"] = util.getSimpleName(util.getRootName(maskfile)); }
            outputFileName = getOutputFileName("chimera", variables);
            accnosFileName = getOutputFileName("accnos", variables);
			
			
			if (m->getControl_pressed()) { delete chimera; for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0;	}
			
			if (chimera->getUnaligned()) { 
				m->mothurOut("Your template sequences are different lengths, please correct."); m->mothurOutEndLine(); 
				delete chimera;
				return 0; 
			}
			templateSeqsLength = chimera->getLength();
		
            numSeqs = driver(outputFileName, fastaFileNames[s], accnosFileName);
            
            if (m->getControl_pressed()) { outputTypes.clear(); util.mothurRemove(outputFileName); util.mothurRemove(accnosFileName); for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  delete chimera; return 0; }
								
			delete chimera;
			
			outputNames.push_back(outputFileName); outputTypes["chimera"].push_back(outputFileName);
			outputNames.push_back(accnosFileName); outputTypes["accnos"].push_back(accnosFileName);
			
			m->mothurOutEndLine();
			m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences.");	m->mothurOutEndLine();
		}
		
		//set accnos file as new current accnosfile
		string currentName = "";
		itTypes = outputTypes.find("accnos");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setAccnosFile(currentName); }
		}
		
		m->mothurOut("\nOutput File Names: \n"); 
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

int ChimeraPintailCommand::driver(string outputFName, string filename, string accnos){
	try {
		ofstream out;
		util.openOutputFile(outputFName, out);
		
		ofstream out2;
		util.openOutputFile(accnos, out2);
		
		ifstream inFASTA;
		util.openInputFile(filename, inFASTA);

		int count = 0;
	
		while (!inFASTA.eof()) {
				
			if (m->getControl_pressed()) {	break;	}
		
			Sequence* candidateSeq = new Sequence(inFASTA);  util.gobble(inFASTA);
				
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


