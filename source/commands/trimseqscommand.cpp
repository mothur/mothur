/*
 *  trimseqscommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 6/6/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "trimseqscommand.h"
#include "needlemanoverlap.hpp"
#include "trimoligos.h"


//**********************************************************************************************************************
vector<string> TrimSeqsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta",false,true,true); parameters.push_back(pfasta);
		CommandParameter poligos("oligos", "InputTypes", "", "", "none", "none", "none","group",false,false,true); parameters.push_back(poligos);
		CommandParameter pqfile("qfile", "InputTypes", "", "", "none", "none", "none","qfile",false,false,true); parameters.push_back(pqfile);
		CommandParameter pname("name", "InputTypes", "", "", "namecount", "none", "none","name",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "namecount", "none", "none","count",false,false,true); parameters.push_back(pcount);
		CommandParameter pflip("flip", "Boolean", "", "F", "", "", "","",false,false,true); parameters.push_back(pflip);
        CommandParameter preorient("checkorient", "Boolean", "", "F", "", "", "","",false,false,true); parameters.push_back(preorient);
		CommandParameter pmaxambig("maxambig", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pmaxambig);
		CommandParameter pmaxhomop("maxhomop", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pmaxhomop);
		CommandParameter pminlength("minlength", "Number", "", "1", "", "", "","",false,false); parameters.push_back(pminlength);
		CommandParameter pmaxlength("maxlength", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pmaxlength);
		CommandParameter ppdiffs("pdiffs", "Number", "", "0", "", "", "","",false,false,true); parameters.push_back(ppdiffs);
		CommandParameter pbdiffs("bdiffs", "Number", "", "0", "", "", "","",false,false,true); parameters.push_back(pbdiffs);
        CommandParameter pldiffs("ldiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pldiffs);
		CommandParameter psdiffs("sdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(psdiffs);
        CommandParameter ptdiffs("tdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(ptdiffs);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pallfiles("allfiles", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pallfiles);
		CommandParameter pkeepforward("keepforward", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pkeepforward);
        CommandParameter plogtransform("logtransform", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(plogtransform);
		CommandParameter pqtrim("qtrim", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pqtrim);
		CommandParameter pqthreshold("qthreshold", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pqthreshold);
		CommandParameter pqaverage("qaverage", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pqaverage);
		CommandParameter prollaverage("rollaverage", "Number", "", "0", "", "", "","",false,false); parameters.push_back(prollaverage);
		CommandParameter pqwindowaverage("qwindowaverage", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pqwindowaverage);
		CommandParameter pqstepsize("qstepsize", "Number", "", "1", "", "", "","",false,false); parameters.push_back(pqstepsize);
		CommandParameter pqwindowsize("qwindowsize", "Number", "", "50", "", "", "","",false,false); parameters.push_back(pqwindowsize);
		CommandParameter pkeepfirst("keepfirst", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pkeepfirst);
		CommandParameter premovelast("removelast", "Number", "", "0", "", "", "","",false,false); parameters.push_back(premovelast);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
			
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string TrimSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The trim.seqs command reads a fastaFile and creates 2 new fasta files, .trim.fasta and scrap.fasta, as well as group files if you provide and oligos file.\n";
		helpString += "The .trim.fasta contains sequences that meet your requirements, and the .scrap.fasta contains those which don't.\n";
		helpString += "The trim.seqs command parameters are fasta, name, count, flip, checkorient, oligos, maxambig, maxhomop, minlength, maxlength, qfile, qthreshold, qaverage, diffs, qtrim, keepfirst, removelast, logtransform and allfiles.\n";
		helpString += "The fasta parameter is required.\n";
		helpString += "The flip parameter will output the reverse compliment of your trimmed sequence. The default is false.\n";
        helpString += "The checkorient parameter will check the reverse compliment of the sequence if the barcodes and primers cannot be found in the forward. The default is false.\n";
		helpString += "The oligos parameter allows you to provide an oligos file.\n";
		helpString += "The name parameter allows you to provide a names file with your fasta file.\n";
        helpString += "The count parameter allows you to provide a count file with your fasta file.\n";
		helpString += "The maxambig parameter allows you to set the maximum number of ambiguous bases allowed. The default is -1.\n";
		helpString += "The maxhomop parameter allows you to set a maximum homopolymer length. \n";
		helpString += "The minlength parameter allows you to set and minimum sequence length. \n";
		helpString += "The maxlength parameter allows you to set and maximum sequence length. \n";
		helpString += "The tdiffs parameter is used to specify the total number of differences allowed in the sequence. The default is pdiffs + bdiffs + sdiffs + ldiffs.\n";
		helpString += "The bdiffs parameter is used to specify the number of differences allowed in the barcode. The default is 0.\n";
		helpString += "The pdiffs parameter is used to specify the number of differences allowed in the primer. The default is 0.\n";
        helpString += "The ldiffs parameter is used to specify the number of differences allowed in the linker. The default is 0.\n";
		helpString += "The sdiffs parameter is used to specify the number of differences allowed in the spacer. The default is 0.\n";
		helpString += "The qfile parameter allows you to provide a quality file.\n";
		helpString += "The qthreshold parameter allows you to set a minimum quality score allowed. \n";
		helpString += "The qaverage parameter allows you to set a minimum average quality score allowed. \n";
		helpString += "The qwindowsize parameter allows you to set a number of bases in a window. Default=50.\n";
		helpString += "The qwindowaverage parameter allows you to set a minimum average quality score allowed over a window. \n";
		helpString += "The rollaverage parameter allows you to set a minimum rolling average quality score allowed over a window. \n";
		helpString += "The qstepsize parameter allows you to set a number of bases to move the window over. Default=1.\n";
        helpString += "The logtransform parameter allows you to indicate you want the averages for the qwindowaverage, rollaverage and qaverage to be calculated using a logtransform. Default=F.\n";
		helpString += "The allfiles parameter will create separate group and fasta file for each grouping. The default is F.\n";
		helpString += "The keepforward parameter allows you to indicate whether you want the forward primer removed or not. The default is F, meaning remove the forward primer.\n";
		helpString += "The qtrim parameter will trim sequence from the point that they fall below the qthreshold and put it in the .trim file if set to true. The default is T.\n";
		helpString += "The keepfirst parameter trims the sequence to the first keepfirst number of bases after the barcode or primers are removed, before the sequence is checked to see if it meets the other requirements. \n";
		helpString += "The removelast removes the last removelast number of bases after the barcode or primers are removed, before the sequence is checked to see if it meets the other requirements.\n";
		helpString += "The trim.seqs command should be in the following format: \n";
		helpString += "trim.seqs(fasta=yourFastaFile, flip=yourFlip, oligos=yourOligos, maxambig=yourMaxambig,  \n";
		helpString += "maxhomop=yourMaxhomop, minlength=youMinlength, maxlength=yourMaxlength)  \n";	
		helpString += "Example trim.seqs(fasta=abrecovery.fasta, flip=..., oligos=..., maxambig=..., maxhomop=..., minlength=..., maxlength=...).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		helpString += "For more details please check out the wiki http://www.mothur.org/wiki/Trim.seqs .\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string TrimSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "qfile") {  pattern = "[filename],[tag],qual"; } 
        else if (type == "fasta") {  pattern = "[filename],[tag],fasta"; } 
        else if (type == "group") {  pattern = "[filename],groups"; }
        else if (type == "name") {  pattern = "[filename],[tag],names"; }
        else if (type == "count") {  pattern = "[filename],[tag],count_table-[filename],count_table"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "TrimSeqsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************

TrimSeqsCommand::TrimSeqsCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["qfile"] = tempOutNames;
		outputTypes["group"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "TrimSeqsCommand");
		exit(1);
	}
}
//***************************************************************************************************************

TrimSeqsCommand::TrimSeqsCommand(string option)  {
	try {
		
		abort = false; calledHelp = false;   
		comboStarts = 0;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["qfile"] = tempOutNames;
			outputTypes["group"] = tempOutNames;
			outputTypes["name"] = tempOutNames;
            outputTypes["count"] = tempOutNames;
			
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
				
				it = parameters.find("oligos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["oligos"] = inputDir + it->second;		}
				}
				
				it = parameters.find("qfile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["qfile"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
				
			}

			
			//check for required parameters
			fastaFile = validParameter.validFile(parameters, "fasta", true);
			if (fastaFile == "not found") { 				
				fastaFile = m->getFastaFile(); 
				if (fastaFile != "") { m->mothurOut("Using " + fastaFile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else if (fastaFile == "not open") { abort = true; }	
			else { m->setFastaFile(fastaFile); }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(fastaFile); //if user entered a file with a path then preserve it	
			}
		
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.validFile(parameters, "flip", false);
			if (temp == "not found")    {	flip = 0;	}
			else {  flip = m->isTrue(temp);		}
		
			temp = validParameter.validFile(parameters, "oligos", true);
			if (temp == "not found"){	oligoFile = "";		}
			else if(temp == "not open"){	abort = true;	} 
			else					{	oligoFile = temp; m->setOligosFile(oligoFile);		}
			
			
			temp = validParameter.validFile(parameters, "maxambig", false);		if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, maxAmbig);  

			temp = validParameter.validFile(parameters, "maxhomop", false);		if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, maxHomoP);  

			temp = validParameter.validFile(parameters, "minlength", false);	if (temp == "not found") { temp = "1"; }
			m->mothurConvert(temp, minLength); 
			
			temp = validParameter.validFile(parameters, "maxlength", false);	if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, maxLength);
			
			temp = validParameter.validFile(parameters, "bdiffs", false);		if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, bdiffs);
			
			temp = validParameter.validFile(parameters, "pdiffs", false);		if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, pdiffs);
            
            temp = validParameter.validFile(parameters, "ldiffs", false);		if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, ldiffs);
            
            temp = validParameter.validFile(parameters, "sdiffs", false);		if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, sdiffs);
			
			temp = validParameter.validFile(parameters, "tdiffs", false);		if (temp == "not found") { int tempTotal = pdiffs + bdiffs + ldiffs + sdiffs;  temp = toString(tempTotal); }
			m->mothurConvert(temp, tdiffs);
			
			if(tdiffs == 0){	tdiffs = bdiffs + pdiffs + ldiffs + sdiffs;	}
			
			temp = validParameter.validFile(parameters, "qfile", true);	
			if (temp == "not found")	{	qFileName = "";		}
			else if(temp == "not open")	{	abort = true;		}
			else						{	qFileName = temp;	m->setQualFile(qFileName); }
			
			temp = validParameter.validFile(parameters, "name", true);	
			if (temp == "not found")	{	nameFile = "";		}
			else if(temp == "not open")	{	nameFile = "";	abort = true;		}
			else						{	nameFile = temp;	m->setNameFile(nameFile); }
            
            countfile = validParameter.validFile(parameters, "count", true);
			if (countfile == "not open") { abort = true; countfile = ""; }	
			else if (countfile == "not found") { countfile = ""; }
			else { m->setCountTableFile(countfile); }
			
            if ((countfile != "") && (nameFile != "")) { m->mothurOut("You must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
			
			temp = validParameter.validFile(parameters, "qthreshold", false);	if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, qThreshold);
			
			temp = validParameter.validFile(parameters, "qtrim", false);		if (temp == "not found") { temp = "t"; }
			qtrim = m->isTrue(temp);

			temp = validParameter.validFile(parameters, "rollaverage", false);	if (temp == "not found") { temp = "0"; }
			convert(temp, qRollAverage);

			temp = validParameter.validFile(parameters, "qwindowaverage", false);if (temp == "not found") { temp = "0"; }
			convert(temp, qWindowAverage);

			temp = validParameter.validFile(parameters, "qwindowsize", false);	if (temp == "not found") { temp = "50"; }
			convert(temp, qWindowSize);

			temp = validParameter.validFile(parameters, "qstepsize", false);	if (temp == "not found") { temp = "1"; }
			convert(temp, qWindowStep);

			temp = validParameter.validFile(parameters, "qaverage", false);		if (temp == "not found") { temp = "0"; }
			convert(temp, qAverage);

			temp = validParameter.validFile(parameters, "keepfirst", false);	if (temp == "not found") { temp = "0"; }
			convert(temp, keepFirst);

			temp = validParameter.validFile(parameters, "removelast", false);	if (temp == "not found") { temp = "0"; }
			convert(temp, removeLast);
			
			temp = validParameter.validFile(parameters, "allfiles", false);		if (temp == "not found") { temp = "F"; }
			allFiles = m->isTrue(temp);
            
            temp = validParameter.validFile(parameters, "keepforward", false);		if (temp == "not found") { temp = "F"; }
			keepforward = m->isTrue(temp);
            
            temp = validParameter.validFile(parameters, "logtransform", false);		if (temp == "not found") { temp = "F"; }
			logtransform = m->isTrue(temp);
            
            temp = validParameter.validFile(parameters, "checkorient", false);		if (temp == "not found") { temp = "F"; }
			reorient = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors); 
			
			
			if(allFiles && (oligoFile == "")){
				m->mothurOut("You selected allfiles, but didn't enter an oligos.  Ignoring the allfiles request."); m->mothurOutEndLine();
			}
			if((qAverage != 0 && qThreshold != 0) && qFileName == ""){
				m->mothurOut("You didn't provide a quality file name, quality criteria will be ignored."); m->mothurOutEndLine();
				qAverage=0;
				qThreshold=0;
			}
			if(!flip && oligoFile=="" && !maxLength && !minLength && (maxAmbig==-1) && !maxHomoP && qFileName == ""){		
				m->mothurOut("You didn't set any options... quiting command."); m->mothurOutEndLine();
				abort = true;
			}
			
            if (countfile == "") {
                if (nameFile == "") {
                    vector<string> files; files.push_back(fastaFile);
                    parser.getNameFile(files);
                }
            }
		}

	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "TrimSeqsCommand");
		exit(1);
	}
}
//***************************************************************************************************************

int TrimSeqsCommand::execute(){
	try{
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
        pairedOligos = false;
		numFPrimers = 0;  //this needs to be initialized
		numRPrimers = 0;
        numSpacers = 0;
        numLinkers = 0;
		createGroup = false;
		vector<vector<string> > fastaFileNames;
		vector<vector<string> > qualFileNames;
		vector<vector<string> > nameFileNames;
		
        map<string, string> variables; 
		variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastaFile));
        variables["[tag]"] = "trim";
		string trimSeqFile = getOutputFileName("fasta",variables);
        string trimQualFile = getOutputFileName("qfile",variables);
		outputNames.push_back(trimSeqFile); outputTypes["fasta"].push_back(trimSeqFile);
        
        variables["[tag]"] = "scrap";
		string scrapSeqFile = getOutputFileName("fasta",variables);
        string scrapQualFile = getOutputFileName("qfile",variables);
		outputNames.push_back(scrapSeqFile); outputTypes["fasta"].push_back(scrapSeqFile);
		
		if (qFileName != "") {
			outputNames.push_back(trimQualFile);
			outputNames.push_back(scrapQualFile);
			outputTypes["qfile"].push_back(trimQualFile);
			outputTypes["qfile"].push_back(scrapQualFile); 
		}
		
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(nameFile));
        variables["[tag]"] = "trim";
		string trimNameFile = getOutputFileName("name",variables);
        variables["[tag]"] = "scrap";
		string scrapNameFile = getOutputFileName("name",variables);
		
		if (nameFile != "") {
			m->readNames(nameFile, nameMap);
			outputNames.push_back(trimNameFile);
			outputNames.push_back(scrapNameFile);
			outputTypes["name"].push_back(trimNameFile);
			outputTypes["name"].push_back(scrapNameFile); 
		}
        
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(countfile));
        variables["[tag]"] = "trim";
        string trimCountFile = getOutputFileName("count",variables);
        variables["[tag]"] = "scrap";
		string scrapCountFile = getOutputFileName("count",variables);
		
		if (countfile != "") {
            CountTable ct;
            ct.readTable(countfile, true, false);
            nameCount = ct.getNameMap();
			outputNames.push_back(trimCountFile);
			outputNames.push_back(scrapCountFile);
			outputTypes["count"].push_back(trimCountFile);
			outputTypes["count"].push_back(scrapCountFile); 
		}

		
		if (m->control_pressed) { return 0; }
		
		string outputGroupFileName;
		if(oligoFile != ""){
			createGroup = getOligos(fastaFileNames, qualFileNames, nameFileNames);
			if ((createGroup) && (countfile == "")){
                map<string, string> myvariables; 
                myvariables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastaFile));
				outputGroupFileName = getOutputFileName("group",myvariables);
				outputNames.push_back(outputGroupFileName); outputTypes["group"].push_back(outputGroupFileName);
			}
		}
        
        if (m->control_pressed) { return 0; }
            
        //fills lines and qlines
		setLines(fastaFile, qFileName);
		
        if(processors == 1){
            driverCreateTrim(fastaFile, qFileName, trimSeqFile, scrapSeqFile, trimQualFile, scrapQualFile, trimNameFile, scrapNameFile, trimCountFile, scrapCountFile, outputGroupFileName, fastaFileNames, qualFileNames, nameFileNames, lines[0], qLines[0]);
        }else{
            createProcessesCreateTrim(fastaFile, qFileName, trimSeqFile, scrapSeqFile, trimQualFile, scrapQualFile, trimNameFile, scrapNameFile, trimCountFile, scrapCountFile, outputGroupFileName, fastaFileNames, qualFileNames, nameFileNames); 
        }	
		
		
		if (m->control_pressed) {  return 0; }			
 	
		if(allFiles){
			map<string, string> uniqueFastaNames;// so we don't add the same groupfile multiple times
			map<string, string>::iterator it;
			set<string> namesToRemove;
			for(int i=0;i<fastaFileNames.size();i++){
				for(int j=0;j<fastaFileNames[0].size();j++){
					if (fastaFileNames[i][j] != "") {
						if (namesToRemove.count(fastaFileNames[i][j]) == 0) {
							if(m->isBlank(fastaFileNames[i][j])){
								m->mothurRemove(fastaFileNames[i][j]);
								namesToRemove.insert(fastaFileNames[i][j]);
							
								if(qFileName != ""){
									m->mothurRemove(qualFileNames[i][j]);
									namesToRemove.insert(qualFileNames[i][j]);
								}
								
								if(nameFile != ""){
									m->mothurRemove(nameFileNames[i][j]);
									namesToRemove.insert(nameFileNames[i][j]);
								}
							}else{	
								it = uniqueFastaNames.find(fastaFileNames[i][j]);
								if (it == uniqueFastaNames.end()) {	
									uniqueFastaNames[fastaFileNames[i][j]] = barcodeNameVector[i];	
								}	
							}
						}
					}
				}
			}
			
			//remove names for outputFileNames, just cleans up the output
			vector<string> outputNames2;
			for(int i = 0; i < outputNames.size(); i++) { if (namesToRemove.count(outputNames[i]) == 0) { outputNames2.push_back(outputNames[i]); } }
			outputNames = outputNames2;
			
            for (it = uniqueFastaNames.begin(); it != uniqueFastaNames.end(); it++) {
                ifstream in;
                m->openInputFile(it->first, in);
                
                ofstream out;
                map<string, string> myvariables; 
                myvariables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(it->first));
                string thisGroupName = "";
                if (countfile == "") { thisGroupName = getOutputFileName("group",myvariables); outputNames.push_back(thisGroupName); outputTypes["group"].push_back(thisGroupName); }
                else {  thisGroupName = getOutputFileName("count",myvariables); outputNames.push_back(thisGroupName); outputTypes["count"].push_back(thisGroupName);  }
                m->openOutputFile(thisGroupName, out);
                
                if (countfile != "") {  out << "Representative_Sequence\ttotal\t" << it->second << endl;  }
                
                while (!in.eof()){
                    if (m->control_pressed) { break; }
                    
                    Sequence currSeq(in); m->gobble(in);
                    if (countfile == "") {  
                        out << currSeq.getName() << '\t' << it->second << endl;  
                        
                        if (nameFile != "") {
                            map<string, string>::iterator itName = nameMap.find(currSeq.getName());
                            if (itName != nameMap.end()) { 
                                vector<string> thisSeqsNames; 
                                m->splitAtChar(itName->second, thisSeqsNames, ',');
                                for (int k = 1; k < thisSeqsNames.size(); k++) { //start at 1 to skip self
                                    out << thisSeqsNames[k] << '\t' << it->second << endl;
                                }
                            }else { m->mothurOut("[ERROR]: " + currSeq.getName() + " is not in your namefile, please correct."); m->mothurOutEndLine(); }							
                        }
                    }else { 
                        map<string, int>::iterator itTotalReps = nameCount.find(currSeq.getName());
                        if (itTotalReps != nameCount.end()) { out << currSeq.getName() << '\t' << itTotalReps->second << '\t' << itTotalReps->second << endl; }
                        else { m->mothurOut("[ERROR]: " + currSeq.getName() + " is not in your count file, please correct."); m->mothurOutEndLine(); }
                    }
                }
                in.close();
                out.close();
            }
            
            if (countfile != "") { //create countfile with group info included
                CountTable* ct = new CountTable();
                ct->readTable(trimCountFile, true, false);
                map<string, int> justTrimmedNames = ct->getNameMap();
                delete ct;
                
                CountTable newCt;
                for (map<string, int>::iterator itCount = groupCounts.begin(); itCount != groupCounts.end(); itCount++) { newCt.addGroup(itCount->first); }
                vector<int> tempCounts; tempCounts.resize(groupCounts.size(), 0);
                for (map<string, int>::iterator itNames = justTrimmedNames.begin(); itNames != justTrimmedNames.end(); itNames++) {
                    newCt.push_back(itNames->first, tempCounts); //add it to the table with no abundance so we can set the groups abundance
                    map<string, string>::iterator it2 = groupMap.find(itNames->first);
                    if (it2 != groupMap.end()) { newCt.setAbund(itNames->first, it2->second, itNames->second); }
                    else { m->mothurOut("[ERROR]: missing group info for " + itNames->first + "."); m->mothurOutEndLine(); m->control_pressed = true; }
                }
                newCt.printTable(trimCountFile);
            }
		}
		
		if (m->control_pressed) {	for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0;	}

		//output group counts
		m->mothurOutEndLine();
		int total = 0;
		if (groupCounts.size() != 0) {  m->mothurOut("Group count: \n");  }
		for (map<string, int>::iterator it = groupCounts.begin(); it != groupCounts.end(); it++) {
			 total += it->second; m->mothurOut(it->first + "\t" + toString(it->second)); m->mothurOutEndLine(); 
		}
		if (total != 0) { m->mothurOut("Total of all groups is " + toString(total)); m->mothurOutEndLine(); }
		
		if (m->control_pressed) {	for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0;	}

		//set fasta file as new current fastafile
		string current = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
		}
		
		itTypes = outputTypes.find("name");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setNameFile(current); }
		}
		
		itTypes = outputTypes.find("qfile");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setQualFile(current); }
		}
		
		itTypes = outputTypes.find("group");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setGroupFile(current); }
		}
        
        itTypes = outputTypes.find("count");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setCountTableFile(current); }
		}

		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;	
			
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "execute");
		exit(1);
	}
}
		
/**************************************************************************************/
int TrimSeqsCommand::driverCreateTrim(string filename, string qFileName, string trimFileName, string scrapFileName, string trimQFileName, string scrapQFileName, string trimNFileName, string scrapNFileName, string trimCFileName, string scrapCFileName, string groupFileName, vector<vector<string> > fastaFileNames, vector<vector<string> > qualFileNames, vector<vector<string> > nameFileNames, linePair line, linePair qline) {	
		
	try {
		
		ofstream trimFASTAFile;
		m->openOutputFile(trimFileName, trimFASTAFile);
		
		ofstream scrapFASTAFile;
		m->openOutputFile(scrapFileName, scrapFASTAFile);
		
		ofstream trimQualFile;
		ofstream scrapQualFile;
		if(qFileName != ""){
			m->openOutputFile(trimQFileName, trimQualFile);
			m->openOutputFile(scrapQFileName, scrapQualFile);
		}
		
		ofstream trimNameFile;
		ofstream scrapNameFile;
		if(nameFile != ""){
			m->openOutputFile(trimNFileName, trimNameFile);
			m->openOutputFile(scrapNFileName, scrapNameFile);
		}
		
        ofstream trimCountFile;
		ofstream scrapCountFile;
		if(countfile != ""){
			m->openOutputFile(trimCFileName, trimCountFile);
			m->openOutputFile(scrapCFileName, scrapCountFile);
            if (line.start == 0) { trimCountFile << "Representative_Sequence\ttotal" << endl; scrapCountFile << "Representative_Sequence\ttotal" << endl; }
		}
		
		ofstream outGroupsFile;
		if ((createGroup) && (countfile == "")){	m->openOutputFile(groupFileName, outGroupsFile);   }
		if(allFiles){
			for (int i = 0; i < fastaFileNames.size(); i++) { //clears old file
				for (int j = 0; j < fastaFileNames[i].size(); j++) { //clears old file
					if (fastaFileNames[i][j] != "") {
						ofstream temp;
						m->openOutputFile(fastaFileNames[i][j], temp);			temp.close();
						if(qFileName != ""){
							m->openOutputFile(qualFileNames[i][j], temp);			temp.close();
						}
						
						if(nameFile != ""){
							m->openOutputFile(nameFileNames[i][j], temp);			temp.close();
						}
					}
				}
			}
		}
		
		ifstream inFASTA;
		m->openInputFile(filename, inFASTA);
		inFASTA.seekg(line.start);
		
		ifstream qFile;
		if(qFileName != "")	{
			m->openInputFile(qFileName, qFile);
			qFile.seekg(qline.start);  
		}
		
		int count = 0;
		bool moreSeqs = 1;
        int numBarcodes = barcodes.size();
		TrimOligos* trimOligos = NULL;
        if (pairedOligos)   {   trimOligos = new TrimOligos(pdiffs, bdiffs, 0, 0, pairedPrimers, pairedBarcodes, false); numBarcodes = pairedBarcodes.size(); }
        else                {   trimOligos = new TrimOligos(pdiffs, bdiffs, ldiffs, sdiffs, primers, barcodes, revPrimer, linker, spacer);  }
        
        TrimOligos* rtrimOligos = NULL;
        if (reorient) {
            //create reoriented primer and barcode pairs
            map<int, oligosPair> rpairedPrimers, rpairedBarcodes;
            for (map<int, oligosPair>::iterator it = pairedPrimers.begin(); it != pairedPrimers.end(); it++) {
                  oligosPair tempPair(reverseOligo((it->second).reverse), (reverseOligo((it->second).forward))); //reversePrimer, rc ForwardPrimer
                rpairedPrimers[it->first] = tempPair;
                //cout  << reverseOligo((it->second).reverse) << '\t' << (reverseOligo((it->second).forward)) << '\t' << primerNameVector[it->first] << endl;
            }
            for (map<int, oligosPair>::iterator it = pairedBarcodes.begin(); it != pairedBarcodes.end(); it++) {
                 oligosPair tempPair(reverseOligo((it->second).reverse), (reverseOligo((it->second).forward))); //reverseBarcode, rc ForwardBarcode
                rpairedBarcodes[it->first] = tempPair;
                 //cout  << reverseOligo((it->second).reverse) << '\t' << (reverseOligo((it->second).forward)) << '\t' << barcodeNameVector[it->first] << endl;
            }
            int index = rpairedBarcodes.size();
            for (map<string, int>::iterator it = barcodes.begin(); it != barcodes.end(); it++) {
                oligosPair tempPair("", reverseOligo((it->first))); //reverseBarcode, rc ForwardBarcode
                rpairedBarcodes[index] = tempPair; index++;
                //cout  << reverseOligo((it->second).reverse) << '\t' << (reverseOligo((it->second).forward)) << '\t' << barcodeNameVector[it->first] << endl;
            }
            
            index = rpairedPrimers.size();
            for (map<string, int>::iterator it = primers.begin(); it != primers.end(); it++) {
                oligosPair tempPair("", reverseOligo((it->first))); //reverseBarcode, rc ForwardBarcode
                rpairedPrimers[index] = tempPair; index++;
                //cout  << reverseOligo((it->second).reverse) << '\t' << (reverseOligo((it->second).forward)) << '\t' << primerNameVector[it->first] << endl;
            }

            rtrimOligos = new TrimOligos(pdiffs, bdiffs, 0, 0, rpairedPrimers, rpairedBarcodes, false); numBarcodes = rpairedBarcodes.size();
        }
        
		while (moreSeqs) {
				
            int obsBDiffs = 0;
            int obsPDiffs = 0;
            
            
			if (m->control_pressed) {
                delete trimOligos; if (reorient) { delete rtrimOligos; }
				inFASTA.close(); trimFASTAFile.close(); scrapFASTAFile.close();
				if ((createGroup) && (countfile == "")) {	 outGroupsFile.close();   }
                if(qFileName != "")	{	qFile.close();	scrapQualFile.close(); trimQualFile.close();	}
                if(nameFile != "")	{	scrapNameFile.close(); trimNameFile.close();	}
                if(countfile != "")	{	scrapCountFile.close(); trimCountFile.close();	}
				for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0;
			}
			
			int success = 1;
			string trashCode = "";
            string commentString = "";
			int currentSeqsDiffs = 0;

			Sequence currSeq(inFASTA); m->gobble(inFASTA);
			//cout << currSeq.getName() << '\t' << currSeq.getUnaligned() << endl;
            Sequence savedSeq(currSeq.getName(), currSeq.getAligned());
            
			QualityScores currQual; QualityScores savedQual;
			if(qFileName != ""){
				currQual = QualityScores(qFile);  m->gobble(qFile);
                savedQual.setName(currQual.getName()); savedQual.setScores(currQual.getScores());
                //cout << currQual.getName() << endl;
			}
			  
			string origSeq = currSeq.getUnaligned();
			if (origSeq != "") {
				
				int barcodeIndex = 0;
				int primerIndex = 0;
				
                
//                cout << currSeq.getName() << '\t'; cout.flush();
                
                if(numLinkers != 0){
					success = trimOligos->stripLinker(currSeq, currQual);
					if(success > ldiffs)		{	trashCode += 'k';	}
					else{ currentSeqsDiffs += success;  }

				}
                
				if(numBarcodes != 0){
					vector<int> results = trimOligos->stripBarcode(currSeq, currQual, barcodeIndex);
                    if (pairedOligos) {
                        success = results[0] + results[2];
                        commentString += "fbdiffs=" + toString(results[0]) + "(" + trimOligos->getCodeValue(results[1], bdiffs) + "), rbdiffs=" + toString(results[2]) + "(" + trimOligos->getCodeValue(results[3], bdiffs) + ") ";
                    }
                    else {
                        success = results[0];
                        commentString += "bdiffs=" + toString(results[0]) + "(" + trimOligos->getCodeValue(results[1], bdiffs) + ") ";
                    }
					if(success > bdiffs)		{  trashCode += 'b'; }
					else{ currentSeqsDiffs += success;  }
				}
                obsBDiffs = success;
                
//                cout << success << '\t'; cout.flush();
                
				//cout << currSeq.getName() << '\t' << currSeq.getUnaligned() << endl;
                if(numSpacers != 0){
					success = trimOligos->stripSpacer(currSeq, currQual);
					if(success > sdiffs)		{	trashCode += 's';	}
					else{ currentSeqsDiffs += success;  }

				}
                
				if(numFPrimers != 0){
					vector<int> results = trimOligos->stripForward(currSeq, currQual, primerIndex, keepforward);
                    if (pairedOligos) {
                        success = results[0] + results[2];
                        commentString += "fpdiffs=" + toString(results[0]) + "(" + trimOligos->getCodeValue(results[1], pdiffs) + "), rpdiffs=" + toString(results[2]) + "(" + trimOligos->getCodeValue(results[3], pdiffs) + ") ";
                    }
                    else {
                        success = results[0];
                        commentString += "fpdiffs=" + toString(results[0]) + "(" + trimOligos->getCodeValue(results[1], pdiffs) + ") ";
                    }
					if(success > pdiffs)		{  trashCode += 'f';  }
					else{ currentSeqsDiffs += success;  }
				}
                obsPDiffs = success;
                
//                cout << success << '\t'; cout.flush();
				
//                cout << currentSeqsDiffs << endl;
                
				if(numRPrimers != 0){
					vector<int> results =  trimOligos->stripReverse(currSeq, currQual);
                    success = results[0];
                    commentString += "rpdiffs=" + toString(results[0]) + "(" + trimOligos->getCodeValue(results[1], pdiffs) + ") ";
                    if(success > pdiffs)		{	trashCode += 'r';	}
                    else{ currentSeqsDiffs += success;  }
				}
                
                if (currentSeqsDiffs > tdiffs)	{	trashCode += 't';   }
                
                if (reorient && (trashCode != "")) { //if you failed and want to check the reverse
                    int thisSuccess = 0;
                    string thisTrashCode = "";
                    string thiscommentString = "";
                    int thisCurrentSeqsDiffs = 0;
                    
                    int thisBarcodeIndex = 0;
                    int thisPrimerIndex = 0;
                    //cout << currSeq.getName() << '\t' << savedSeq.getUnaligned() << endl;
                    if(numBarcodes != 0){
                        vector<int> results =  rtrimOligos->stripBarcode(savedSeq, savedQual, thisBarcodeIndex);
                        if (pairedOligos) {
                            thisSuccess = results[0] + results[2];
                            thiscommentString += "fbdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], bdiffs) + "), rbdiffs=" + toString(results[2]) + "(" + rtrimOligos->getCodeValue(results[3], bdiffs) + ") ";
                        }
                        else {
                            thisSuccess = results[0];
                            thiscommentString += "bdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], bdiffs) + ") ";
                        }
                        if(thisSuccess > bdiffs)		{ thisTrashCode += "b"; }
                        else{ thisCurrentSeqsDiffs += thisSuccess;  }
                    }
                    
                    int revBDiffs = thisSuccess;
                    
//                    cout << thisSuccess << '\t'; cout.flush();

                    
                    //cout << currSeq.getName() << '\t' << savedSeq.getUnaligned() << endl;
                    if(numFPrimers != 0){
                        vector<int> results =  rtrimOligos->stripForward(savedSeq, savedQual, thisPrimerIndex, keepforward);
                        if (pairedOligos) {
                            thisSuccess = results[0] + results[2];
                            thiscommentString += "fpdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], pdiffs) + "), rpdiffs=" + toString(results[2]) + "(" + rtrimOligos->getCodeValue(results[3], pdiffs) + ") ";
                        }
                        else {
                            thisSuccess = results[0];
                            thiscommentString += "pdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], pdiffs) + ") ";
                        }
                        if(thisSuccess > pdiffs)		{ thisTrashCode += "f"; }
                        else{ thisCurrentSeqsDiffs += thisSuccess;  }
                    }
                    int revPDiffs = thisSuccess;

//                    cout << thisSuccess << '\t'; cout.flush();

                    if (thisCurrentSeqsDiffs > tdiffs)	{	thisTrashCode += 't';   }
                    
//                    cout << thisCurrentSeqsDiffs << endl;

                    
                    if (thisTrashCode == "") {
                        obsPDiffs = revPDiffs;
                        obsBDiffs = revBDiffs;
                        
                        trashCode = thisTrashCode;
                        success = thisSuccess;
                        currentSeqsDiffs = thisCurrentSeqsDiffs;
                        barcodeIndex = thisBarcodeIndex;
                        commentString = thiscommentString;
                        primerIndex = thisPrimerIndex;
                        savedSeq.reverseComplement();
                        currSeq.setAligned(savedSeq.getAligned());
                        if(qFileName != ""){
                            savedQual.flipQScores();
                            currQual.setScores(savedQual.getScores());
                        }
                    }else { trashCode += "(" + thisTrashCode + ")";  }
                }
                
				if(keepFirst != 0){
					success = keepFirstTrim(currSeq, currQual);
				}
				
				if(removeLast != 0){
					success = removeLastTrim(currSeq, currQual);
					if(!success)				{	trashCode += 'l';	}
				}

				
				if(qFileName != ""){
					int origLength = currSeq.getNumBases();
					
					if(qThreshold != 0)			{	success = currQual.stripQualThreshold(currSeq, qThreshold);			}
					else if(qAverage != 0)		{	success = currQual.cullQualAverage(currSeq, qAverage, logtransform);				}
					else if(qRollAverage != 0)	{	success = currQual.stripQualRollingAverage(currSeq, qRollAverage, logtransform);	}
					else if(qWindowAverage != 0){	success = currQual.stripQualWindowAverage(currSeq, qWindowStep, qWindowSize, qWindowAverage, logtransform);	}
					else						{	success = 1;				}
					
					//you don't want to trim, if it fails above then scrap it
					if ((!qtrim) && (origLength != currSeq.getNumBases())) {  success = 0; }
					
					if(!success)				{	trashCode += 'q';	}
				}				
		
				if(minLength > 0 || maxLength > 0){
					success = cullLength(currSeq);
					if(!success)				{	trashCode += 'l';	}
				}
				if(maxHomoP > 0){
					success = cullHomoP(currSeq);
					if(!success)				{	trashCode += 'h';	}
				}
				if(maxAmbig != -1){
					success = cullAmbigs(currSeq);
					if(!success)				{	trashCode += 'n';	}
				}
				
				if(flip){		// should go last			
					currSeq.reverseComplement();
					if(qFileName != ""){
						currQual.flipQScores();	
					}
				}
				
                if (m->debug) { m->mothurOut("[DEBUG]: " + currSeq.getName() + ", trashcode= " + trashCode); if (trashCode.length() != 0) { m->mothurOutEndLine(); } }

                string seqComment = currSeq.getComment();
                currSeq.setComment("\t" + commentString + "\t" + seqComment);
                
				if(trashCode.length() == 0){
                    string thisGroup = "";
                    if (createGroup) {
						if(numBarcodes != 0){
							thisGroup = barcodeNameVector[barcodeIndex];
							if (numFPrimers != 0) {
								if (primerNameVector[primerIndex] != "") { 
									if(thisGroup != "") {
										thisGroup += "." + primerNameVector[primerIndex]; 
									}else {
										thisGroup = primerNameVector[primerIndex]; 
									}
								} 
							}
                        }
                    }
                    
                    int pos = thisGroup.find("ignore");
                    if (pos == string::npos) {
                        currSeq.setAligned(currSeq.getUnaligned());
                        currSeq.printSequence(trimFASTAFile);
                        
                        if(qFileName != ""){
                            currQual.printQScores(trimQualFile);
                        }
                        
                        
                        if(nameFile != ""){
                            map<string, string>::iterator itName = nameMap.find(currSeq.getName());
                            if (itName != nameMap.end()) {  trimNameFile << itName->first << '\t' << itName->second << endl; }
                            else { m->mothurOut("[ERROR]: " + currSeq.getName() + " is not in your namefile, please correct."); m->mothurOutEndLine(); }
                        }
                        
                        int numRedundants = 0;
                        if (countfile != "") {
                            map<string, int>::iterator itCount = nameCount.find(currSeq.getName());
                            if (itCount != nameCount.end()) { 
                                trimCountFile << itCount->first << '\t' << itCount->second << endl;
                                numRedundants = itCount->second-1;
                            }else { m->mothurOut("[ERROR]: " + currSeq.getName() + " is not in your count file, please correct."); m->mothurOutEndLine(); }
                        }
                        
                        if (createGroup) {
                            if(numBarcodes != 0){
                                                                
                                if (m->debug) { m->mothurOut(", group= " + thisGroup + "\n"); }
                                
                                if (countfile == "") { outGroupsFile << currSeq.getName() << '\t' << thisGroup << endl; }
                                else {   groupMap[currSeq.getName()] = thisGroup; }
                                
                                if (nameFile != "") {
                                    map<string, string>::iterator itName = nameMap.find(currSeq.getName());
                                    if (itName != nameMap.end()) { 
                                        vector<string> thisSeqsNames; 
                                        m->splitAtChar(itName->second, thisSeqsNames, ',');
                                        numRedundants = thisSeqsNames.size()-1; //we already include ourselves below
                                        for (int k = 1; k < thisSeqsNames.size(); k++) { //start at 1 to skip self
                                            outGroupsFile << thisSeqsNames[k] << '\t' << thisGroup << endl;
                                        }
                                    }else { m->mothurOut("[ERROR]: " + currSeq.getName() + " is not in your namefile, please correct."); m->mothurOutEndLine(); }							
                                }
                                
                                map<string, int>::iterator it = groupCounts.find(thisGroup);
                                if (it == groupCounts.end()) {	groupCounts[thisGroup] = 1 + numRedundants; }
                                else { groupCounts[it->first] += (1 + numRedundants); }
								
                            }
                        }
                        
                        if(allFiles){
                            ofstream output;
                            m->openOutputFileAppend(fastaFileNames[barcodeIndex][primerIndex], output);
                            currSeq.printSequence(output);
                            output.close();
                            
                            if(qFileName != ""){
                                m->openOutputFileAppend(qualFileNames[barcodeIndex][primerIndex], output);
                                currQual.printQScores(output);
                                output.close();							
                            }
                            
                            if(nameFile != ""){
                                map<string, string>::iterator itName = nameMap.find(currSeq.getName());
                                if (itName != nameMap.end()) { 
                                    m->openOutputFileAppend(nameFileNames[barcodeIndex][primerIndex], output);
                                    output << itName->first << '\t' << itName->second << endl; 
                                    output.close();
                                }else { m->mothurOut("[ERROR]: " + currSeq.getName() + " is not in your namefile, please correct."); m->mothurOutEndLine(); }
                            }
                        }
                    }
				}
				else{
					if(nameFile != ""){ //needs to be before the currSeq name is changed
						map<string, string>::iterator itName = nameMap.find(currSeq.getName());
						if (itName != nameMap.end()) {  scrapNameFile << itName->first << '\t' << itName->second << endl; }
						else { m->mothurOut("[ERROR]: " + currSeq.getName() + " is not in your namefile, please correct."); m->mothurOutEndLine(); }
					}
                    if (countfile != "") {
                        map<string, int>::iterator itCount = nameCount.find(currSeq.getName());
                        if (itCount != nameCount.end()) { 
                            trimCountFile << itCount->first << '\t' << itCount->second << endl;
                        }else { m->mothurOut("[ERROR]: " + currSeq.getName() + " is not in your count file, please correct."); m->mothurOutEndLine(); }
                    }
                    
					currSeq.setName(currSeq.getName() + '|' + trashCode);
					currSeq.setUnaligned(origSeq);
					currSeq.setAligned(origSeq);
					currSeq.printSequence(scrapFASTAFile);
					if(qFileName != ""){
						currQual.printQScores(scrapQualFile);
					}
				}
				count++;
			}
			
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				unsigned long long pos = inFASTA.tellg();
				if ((pos == -1) || (pos >= line.end)) { break; }
			
			#else
				if (inFASTA.eof()) { break; }
			#endif
			
			//report progress
			if((count) % 1000 == 0){	m->mothurOutJustToScreen(toString(count)+"\n"); 		}
			
		}
		//report progress
		if((count) % 1000 != 0){	m->mothurOutJustToScreen(toString(count)+"\n");		}
		
		delete trimOligos;
        if (reorient) { delete rtrimOligos; }
		inFASTA.close();
		trimFASTAFile.close();
		scrapFASTAFile.close();
		if (createGroup) {	 outGroupsFile.close();   }
		if(qFileName != "")	{	qFile.close();	scrapQualFile.close(); trimQualFile.close();	}
		if(nameFile != "")	{	scrapNameFile.close(); trimNameFile.close();	}
        if(countfile != "")	{	scrapCountFile.close(); trimCountFile.close();	}
		
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "driverCreateTrim");
		exit(1);
	}
}

/**************************************************************************************************/

int TrimSeqsCommand::createProcessesCreateTrim(string filename, string qFileName, string trimFASTAFileName, string scrapFASTAFileName, string trimQualFileName, string scrapQualFileName, string trimNameFileName, string scrapNameFileName, string trimCountFileName, string scrapCountFileName, string groupFile, vector<vector<string> > fastaFileNames, vector<vector<string> > qualFileNames, vector<vector<string> > nameFileNames) {
	try {
        
        int process = 1;
		int exitCommand = 1;
		processIDS.clear();
        bool recalc = false;
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				
				vector<vector<string> > tempFASTAFileNames = fastaFileNames;
				vector<vector<string> > tempPrimerQualFileNames = qualFileNames;
				vector<vector<string> > tempNameFileNames = nameFileNames;

				if(allFiles){
					ofstream temp;

					for(int i=0;i<tempFASTAFileNames.size();i++){
						for(int j=0;j<tempFASTAFileNames[i].size();j++){
							if (tempFASTAFileNames[i][j] != "") {
								tempFASTAFileNames[i][j] += toString(getpid()) + ".temp";
								m->openOutputFile(tempFASTAFileNames[i][j], temp);			temp.close();

								if(qFileName != ""){
									tempPrimerQualFileNames[i][j] += toString(getpid()) + ".temp";
									m->openOutputFile(tempPrimerQualFileNames[i][j], temp);		temp.close();
								}
								if(nameFile != ""){
									tempNameFileNames[i][j] += toString(getpid()) + ".temp";
									m->openOutputFile(tempNameFileNames[i][j], temp);		temp.close();
								}
							}
						}
					}
				}
							
				driverCreateTrim(filename,
								 qFileName,
								 (trimFASTAFileName + toString(getpid()) + ".temp"),
								 (scrapFASTAFileName + toString(getpid()) + ".temp"),
								 (trimQualFileName + toString(getpid()) + ".temp"),
								 (scrapQualFileName + toString(getpid()) + ".temp"),
								 (trimNameFileName + toString(getpid()) + ".temp"),
								 (scrapNameFileName + toString(getpid()) + ".temp"),
                                 (trimCountFileName + toString(getpid()) + ".temp"),
								 (scrapCountFileName + toString(getpid()) + ".temp"),
								 (groupFile + toString(getpid()) + ".temp"),
								 tempFASTAFileNames,
								 tempPrimerQualFileNames,
								 tempNameFileNames,
								 lines[process],
								 qLines[process]);
                
                if (m->debug) { m->mothurOut("[DEBUG]: " + toString(lines[process].start) + '\t' + toString(qLines[process].start) + '\t' + toString(getpid()) + '\n'); }
				
				//pass groupCounts to parent
				if(createGroup){
					ofstream out;
					string tempFile = filename + toString(getpid()) + ".num.temp";
					m->openOutputFile(tempFile, out);
					
					out << groupCounts.size() << endl;
					
					for (map<string, int>::iterator it = groupCounts.begin(); it != groupCounts.end(); it++) {
						out << it->first << '\t' << it->second << endl;
					}
                    
                    out << groupMap.size() << endl;
                    for (map<string, string>::iterator it = groupMap.begin(); it != groupMap.end(); it++) {
						out << it->first << '\t' << it->second << endl;
					}
					out.close();
				}
				exit(0);
			}else { 
                m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(process) + "\n"); processors = process;
                for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                //wait to die
                for (int i=0;i<processIDS.size();i++) {
                    int temp = processIDS[i];
                    wait(&temp);
                }
                m->control_pressed = false;
                for (int i=0;i<processIDS.size();i++) {
                    m->mothurRemove(trimFASTAFileName + (toString(processIDS[i]) + ".temp"));
                    m->mothurRemove(scrapFASTAFileName + (toString(processIDS[i]) + ".temp"));
                    m->mothurRemove(trimQualFileName + (toString(processIDS[i]) + ".temp"));
                    m->mothurRemove(scrapQualFileName + (toString(processIDS[i]) + ".temp"));
                    m->mothurRemove(trimNameFileName + (toString(processIDS[i]) + ".temp"));
                    m->mothurRemove(scrapNameFileName + (toString(processIDS[i]) + ".temp"));
                    m->mothurRemove(trimCountFileName + (toString(processIDS[i]) + ".temp"));
                    m->mothurRemove(scrapCountFileName + (toString(processIDS[i]) + ".temp"));
                    m->mothurRemove(groupFile + (toString(processIDS[i]) + ".temp"));
                    if (createGroup) {
                        string tempFile = filename + (toString(processIDS[i])) + ".num.temp";
                        m->mothurRemove(tempFile);
                    }
                    if(allFiles){
                        for(int i=0;i<fastaFileNames.size();i++){
                            for(int j=0;j<fastaFileNames[0].size();j++){
                                if (fastaFileNames[i][j] != "") {
                                    string tempFile = fastaFileNames[i][j] +(toString(processIDS[i])) + ".temp";
                                    m->mothurRemove(tempFile);
                                    
                                    if(qFileName != ""){
                                        string tempFile = qualFileNames[i][j] +(toString(processIDS[i])) + ".temp";
                                        m->mothurRemove(tempFile);
                                    }
                                    if(nameFile != ""){
                                        string tempFile = nameFileNames[i][j] +(toString(processIDS[i])) + ".temp";
                                        m->mothurRemove(tempFile);
                                    }
                                }
                            }
                        }
                    }
                }
                recalc = true;
                break;
			}
		}
		
        
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->control_pressed = false;  for (int i=0;i<processIDS.size();i++) {}processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            
            //redo file divide
            lines.clear();
            setLines(fastaFile, qFileName);
            
            exitCommand = 1;
            processIDS.resize(0);
            process = 1;
            
            while (process != processors) {
                int pid = fork();
                
                if (pid > 0) {
                    processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                    process++;
                }else if (pid == 0){
                    
                    vector<vector<string> > tempFASTAFileNames = fastaFileNames;
                    vector<vector<string> > tempPrimerQualFileNames = qualFileNames;
                    vector<vector<string> > tempNameFileNames = nameFileNames;
                    
                    if(allFiles){
                        ofstream temp;
                        
                        for(int i=0;i<tempFASTAFileNames.size();i++){
                            for(int j=0;j<tempFASTAFileNames[i].size();j++){
                                if (tempFASTAFileNames[i][j] != "") {
                                    tempFASTAFileNames[i][j] += toString(getpid()) + ".temp";
                                    m->openOutputFile(tempFASTAFileNames[i][j], temp);			temp.close();
                                    
                                    if(qFileName != ""){
                                        tempPrimerQualFileNames[i][j] += toString(getpid()) + ".temp";
                                        m->openOutputFile(tempPrimerQualFileNames[i][j], temp);		temp.close();
                                    }
                                    if(nameFile != ""){
                                        tempNameFileNames[i][j] += toString(getpid()) + ".temp";
                                        m->openOutputFile(tempNameFileNames[i][j], temp);		temp.close();
                                    }
                                }
                            }
                        }
                    }
                    
                    driverCreateTrim(filename,
                                     qFileName,
                                     (trimFASTAFileName + toString(getpid()) + ".temp"),
                                     (scrapFASTAFileName + toString(getpid()) + ".temp"),
                                     (trimQualFileName + toString(getpid()) + ".temp"),
                                     (scrapQualFileName + toString(getpid()) + ".temp"),
                                     (trimNameFileName + toString(getpid()) + ".temp"),
                                     (scrapNameFileName + toString(getpid()) + ".temp"),
                                     (trimCountFileName + toString(getpid()) + ".temp"),
                                     (scrapCountFileName + toString(getpid()) + ".temp"),
                                     (groupFile + toString(getpid()) + ".temp"),
                                     tempFASTAFileNames,
                                     tempPrimerQualFileNames,
                                     tempNameFileNames,
                                     lines[process],
                                     qLines[process]);
                    
                    if (m->debug) { m->mothurOut("[DEBUG]: " + toString(lines[process].start) + '\t' + toString(qLines[process].start) + '\t' + toString(getpid()) + '\n'); }
                    
                    //pass groupCounts to parent
                    if(createGroup){
                        ofstream out;
                        string tempFile = filename + toString(getpid()) + ".num.temp";
                        m->openOutputFile(tempFile, out);
                        
                        out << groupCounts.size() << endl;
                        
                        for (map<string, int>::iterator it = groupCounts.begin(); it != groupCounts.end(); it++) {
                            out << it->first << '\t' << it->second << endl;
                        }
                        
                        out << groupMap.size() << endl;
                        for (map<string, string>::iterator it = groupMap.begin(); it != groupMap.end(); it++) {
                            out << it->first << '\t' << it->second << endl;
                        }
                        out.close();
                    }
                    exit(0);
                }else { 
                    m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
                    for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                    exit(0);
                }
            }            
        }

        
		//parent do my part
		ofstream temp;
		m->openOutputFile(trimFASTAFileName, temp);		temp.close();
		m->openOutputFile(scrapFASTAFileName, temp);	temp.close();
		if(qFileName != ""){
			m->openOutputFile(trimQualFileName, temp);		temp.close();
			m->openOutputFile(scrapQualFileName, temp);		temp.close();
		}
		if (nameFile != "") {
			m->openOutputFile(trimNameFileName, temp);		temp.close();
			m->openOutputFile(scrapNameFileName, temp);		temp.close();
		}
        if (countfile != "") {
			m->openOutputFile(trimCountFileName, temp);		temp.close();
			m->openOutputFile(scrapCountFileName, temp);		temp.close();
		}

		driverCreateTrim(filename, qFileName, trimFASTAFileName, scrapFASTAFileName, trimQualFileName, scrapQualFileName, trimNameFileName, scrapNameFileName, trimCountFileName, scrapCountFileName, groupFile, fastaFileNames, qualFileNames, nameFileNames, lines[0], qLines[0]);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
#else
        //////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the trimData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<trimData*> pDataArray; 
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
		
		//Create processor worker threads.
		for( int h=0; h<processors-1; h++){
			
            string extension = "";
			if (h != 0) { extension = toString(h) + ".temp"; processIDS.push_back(h); }
            vector<vector<string> > tempFASTAFileNames = fastaFileNames;
            vector<vector<string> > tempPrimerQualFileNames = qualFileNames;
            vector<vector<string> > tempNameFileNames = nameFileNames;
            
            if(allFiles){
                ofstream temp;
                
                for(int i=0;i<tempFASTAFileNames.size();i++){
                    for(int j=0;j<tempFASTAFileNames[i].size();j++){
                        if (tempFASTAFileNames[i][j] != "") {
                            tempFASTAFileNames[i][j] += extension;
                            m->openOutputFile(tempFASTAFileNames[i][j], temp);			temp.close();
                            
                            if(qFileName != ""){
                                tempPrimerQualFileNames[i][j] += extension;
                                m->openOutputFile(tempPrimerQualFileNames[i][j], temp);		temp.close();
                            }
                            if(nameFile != ""){
                                tempNameFileNames[i][j] += extension;
                                m->openOutputFile(tempNameFileNames[i][j], temp);		temp.close();
                            }
                        }
                    }
                }
            }

            
			trimData* tempTrim = new trimData(filename,
                                              qFileName, nameFile, countfile,
                                              (trimFASTAFileName+extension),
                                              (scrapFASTAFileName+extension),
                                              (trimQualFileName+extension),
                                              (scrapQualFileName+extension),
                                              (trimNameFileName+extension),
                                              (scrapNameFileName+extension),
                                              (trimCountFileName+extension),
                                              (scrapCountFileName+extension),
                                              (groupFile+extension),
                                              tempFASTAFileNames,
                                              tempPrimerQualFileNames,
                                              tempNameFileNames,
                                              lines[h].start, lines[h].end, qLines[h].start, qLines[h].end, m,
                                              pdiffs, bdiffs, ldiffs, sdiffs, tdiffs, primers, barcodes, revPrimer, linker, spacer, pairedBarcodes, pairedPrimers, pairedOligos,
                                             primerNameVector, barcodeNameVector, createGroup, allFiles, keepforward, keepFirst, removeLast,
                                              qWindowStep, qWindowSize, qWindowAverage, qtrim, qThreshold, qAverage, qRollAverage, logtransform, 
                                             minLength, maxAmbig, maxHomoP, maxLength, flip, reorient, nameMap, nameCount);
			pDataArray.push_back(tempTrim);
            
			hThreadArray[h] = CreateThread(NULL, 0, MyTrimThreadFunction, pDataArray[h], 0, &dwThreadIdArray[h]);   
		}
        
        //parent do my part
		ofstream temp;
		m->openOutputFile(trimFASTAFileName, temp);		temp.close();
		m->openOutputFile(scrapFASTAFileName, temp);	temp.close();
		if(qFileName != ""){
			m->openOutputFile(trimQualFileName, temp);		temp.close();
			m->openOutputFile(scrapQualFileName, temp);		temp.close();
		}
		if (nameFile != "") {
			m->openOutputFile(trimNameFileName, temp);		temp.close();
			m->openOutputFile(scrapNameFileName, temp);		temp.close();
		}
        vector<vector<string> > tempFASTAFileNames = fastaFileNames;
        vector<vector<string> > tempPrimerQualFileNames = qualFileNames;
        vector<vector<string> > tempNameFileNames = nameFileNames;
        if(allFiles){
            ofstream temp;
            string extension = toString(processors-1) + ".temp";
            for(int i=0;i<tempFASTAFileNames.size();i++){
                for(int j=0;j<tempFASTAFileNames[i].size();j++){
                    if (tempFASTAFileNames[i][j] != "") {
                        tempFASTAFileNames[i][j] += extension;
                        m->openOutputFile(tempFASTAFileNames[i][j], temp);			temp.close();
                        
                        if(qFileName != ""){
                            tempPrimerQualFileNames[i][j] += extension;
                            m->openOutputFile(tempPrimerQualFileNames[i][j], temp);		temp.close();
                        }
                        if(nameFile != ""){
                            tempNameFileNames[i][j] += extension;
                            m->openOutputFile(tempNameFileNames[i][j], temp);		temp.close();
                        }
                    }
                }
            }
        }
        
		driverCreateTrim(filename, qFileName, (trimFASTAFileName + toString(processors-1) + ".temp"), (scrapFASTAFileName + toString(processors-1) + ".temp"), (trimQualFileName + toString(processors-1) + ".temp"), (scrapQualFileName + toString(processors-1) + ".temp"), (trimNameFileName + toString(processors-1) + ".temp"), (scrapNameFileName + toString(processors-1) + ".temp"), (trimCountFileName + toString(processors-1) + ".temp"), (scrapCountFileName + toString(processors-1) + ".temp"), (groupFile + toString(processors-1) + ".temp"), tempFASTAFileNames, tempPrimerQualFileNames, tempNameFileNames, lines[processors-1], qLines[processors-1]);
        processIDS.push_back(processors-1);

        
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
            if (pDataArray[i]->count != pDataArray[i]->lineEnd) {
                m->mothurOut("[ERROR]: process " + toString(i) + " only processed " + toString(pDataArray[i]->count) + " of " + toString(pDataArray[i]->lineEnd) + " sequences assigned to it, quitting. \n"); m->control_pressed = true;
            }
			for (map<string, int>::iterator it = pDataArray[i]->groupCounts.begin(); it != pDataArray[i]->groupCounts.end(); it++) {
                map<string, int>::iterator it2 = groupCounts.find(it->first);
                if (it2 == groupCounts.end()) {	groupCounts[it->first] = it->second; }
                else { groupCounts[it->first] += it->second; }
            }
            for (map<string, string>::iterator it = pDataArray[i]->groupMap.begin(); it != pDataArray[i]->groupMap.end(); it++) {
                map<string, string>::iterator it2 = groupMap.find(it->first);
                if (it2 == groupMap.end()) {	groupMap[it->first] = it->second; }
                else { m->mothurOut("[ERROR]: " + it->first + " is in your fasta file more than once. Sequence names must be unique. please correct.\n");  }
            }
            CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
        
#endif		
        
        
        //append files
		for(int i=0;i<processIDS.size();i++){
			
			m->mothurOut("Appending files from process " + toString(processIDS[i])); m->mothurOutEndLine();
			
			m->appendFiles((trimFASTAFileName + toString(processIDS[i]) + ".temp"), trimFASTAFileName);
			m->mothurRemove((trimFASTAFileName + toString(processIDS[i]) + ".temp"));
			m->appendFiles((scrapFASTAFileName + toString(processIDS[i]) + ".temp"), scrapFASTAFileName);
			m->mothurRemove((scrapFASTAFileName + toString(processIDS[i]) + ".temp"));
			
			if(qFileName != ""){
				m->appendFiles((trimQualFileName + toString(processIDS[i]) + ".temp"), trimQualFileName);
				m->mothurRemove((trimQualFileName + toString(processIDS[i]) + ".temp"));
				m->appendFiles((scrapQualFileName + toString(processIDS[i]) + ".temp"), scrapQualFileName);
				m->mothurRemove((scrapQualFileName + toString(processIDS[i]) + ".temp"));
			}
			
			if(nameFile != ""){
				m->appendFiles((trimNameFileName + toString(processIDS[i]) + ".temp"), trimNameFileName);
				m->mothurRemove((trimNameFileName + toString(processIDS[i]) + ".temp"));
				m->appendFiles((scrapNameFileName + toString(processIDS[i]) + ".temp"), scrapNameFileName);
				m->mothurRemove((scrapNameFileName + toString(processIDS[i]) + ".temp"));
			}
            
            if(countfile != ""){
				m->appendFiles((trimCountFileName + toString(processIDS[i]) + ".temp"), trimCountFileName);
				m->mothurRemove((trimCountFileName + toString(processIDS[i]) + ".temp"));
				m->appendFiles((scrapCountFileName + toString(processIDS[i]) + ".temp"), scrapCountFileName);
				m->mothurRemove((scrapCountFileName + toString(processIDS[i]) + ".temp"));
			}
			
			if((createGroup)&&(countfile == "")){
				m->appendFiles((groupFile + toString(processIDS[i]) + ".temp"), groupFile);
				m->mothurRemove((groupFile + toString(processIDS[i]) + ".temp"));
			}
			
			
			if(allFiles){
				for(int j=0;j<fastaFileNames.size();j++){
					for(int k=0;k<fastaFileNames[j].size();k++){
						if (fastaFileNames[j][k] != "") {
							m->appendFiles((fastaFileNames[j][k] + toString(processIDS[i]) + ".temp"), fastaFileNames[j][k]);
							m->mothurRemove((fastaFileNames[j][k] + toString(processIDS[i]) + ".temp"));
							
							if(qFileName != ""){
								m->appendFiles((qualFileNames[j][k] + toString(processIDS[i]) + ".temp"), qualFileNames[j][k]);
								m->mothurRemove((qualFileNames[j][k] + toString(processIDS[i]) + ".temp"));
							}
							
							if(nameFile != ""){
								m->appendFiles((nameFileNames[j][k] + toString(processIDS[i]) + ".temp"), nameFileNames[j][k]);
								m->mothurRemove((nameFileNames[j][k] + toString(processIDS[i]) + ".temp"));
							}
						}
					}
				}
			}
			
            #if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
			if(createGroup){
				ifstream in;
				string tempFile =  filename + toString(processIDS[i]) + ".num.temp";
				m->openInputFile(tempFile, in);
				int tempNum;
				string group;
				
				in >> tempNum; m->gobble(in);
				
				if (tempNum != 0) {
					for (int i = 0; i < tempNum; i++) { 
                        int groupNum;
						in >> group >> groupNum; m->gobble(in);
                        
						map<string, int>::iterator it = groupCounts.find(group);
						if (it == groupCounts.end()) {	groupCounts[group] = groupNum; }
						else { groupCounts[it->first] += groupNum; }
					}
				}
                in >> tempNum; m->gobble(in);
                if (tempNum != 0) {
					for (int i = 0; i < tempNum; i++) { 
                        string group, seqName;
						in >> seqName >> group; m->gobble(in);
                        
						map<string, string>::iterator it = groupMap.find(seqName);
						if (it == groupMap.end()) {	groupMap[seqName] = group; }
						else { m->mothurOut("[ERROR]: " + seqName + " is in your fasta file more than once. Sequence names must be unique. please correct.\n");  }
					}
				}
                
				in.close(); m->mothurRemove(tempFile);
			}
            #endif
		}

        return exitCommand;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "createProcessesCreateTrim");
		exit(1);
	}
}

/**************************************************************************************************/

int TrimSeqsCommand::setLines(string filename, string qfilename) {
	try {
        
        vector<unsigned long long> fastaFilePos;
		vector<unsigned long long> qfileFilePos;
		
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		//set file positions for fasta file
		fastaFilePos = m->divideFile(filename, processors);
		
		//get name of first sequence in each chunk
		map<string, int> firstSeqNames;
		for (int i = 0; i < (fastaFilePos.size()-1); i++) {
			ifstream in;
			m->openInputFile(filename, in);
			in.seekg(fastaFilePos[i]);
            
            //adjust start if null strings
            if (i == 0) {  m->zapGremlins(in); m->gobble(in);  }
		
			Sequence temp(in); 
			firstSeqNames[temp.getName()] = i;
		
			in.close();
		}
		
		if(qfilename != "")	{
            //seach for filePos of each first name in the qfile and save in qfileFilePos
            ifstream inQual;
            m->openInputFile(qfilename, inQual);
            
            string input;
            while(!inQual.eof()){	
                input = m->getline(inQual);
                
                if (input.length() != 0) {
                    if(input[0] == '>'){ //this is a sequence name line
                        istringstream nameStream(input);
                        
                        string sname = "";  nameStream >> sname;
                        sname = sname.substr(1);
                        
                        m->checkName(sname);
                        
                        map<string, int>::iterator it = firstSeqNames.find(sname);
                        
                        if(it != firstSeqNames.end()) { //this is the start of a new chunk
                            unsigned long long pos = inQual.tellg(); 
                            qfileFilePos.push_back(pos - input.length() - 1);	
                            firstSeqNames.erase(it);
                        }
                    }
                }
                
                if (firstSeqNames.size() == 0) { break; }
            }
            inQual.close();
            
            
            if (firstSeqNames.size() != 0) { 
                for (map<string, int>::iterator it = firstSeqNames.begin(); it != firstSeqNames.end(); it++) {
                    m->mothurOut(it->first + " is in your fasta file and not in your quality file, not using quality file."); m->mothurOutEndLine();
                }
                qFileName = "";
                return processors;
            }
            
            //get last file position of qfile
            FILE * pFile;
            unsigned long long size;
            
            //get num bytes in file
            pFile = fopen (qfilename.c_str(),"rb");
            if (pFile==NULL) perror ("Error opening file");
            else{
                fseek (pFile, 0, SEEK_END);
                size=ftell (pFile);
                fclose (pFile);
            }
            
            qfileFilePos.push_back(size);
        }
        
        for (int i = 0; i < (fastaFilePos.size()-1); i++) {
            if (m->debug) { m->mothurOut("[DEBUG]: " + toString(i) +'\t' + toString(fastaFilePos[i]) + '\t' + toString(fastaFilePos[i+1]) + '\n'); }
			lines.push_back(linePair(fastaFilePos[i], fastaFilePos[(i+1)]));
			if (qfilename != "") {  qLines.push_back(linePair(qfileFilePos[i], qfileFilePos[(i+1)]));  }
		}	
		if(qfilename == "")	{	qLines = lines;	} //files with duds
		
		return processors;
		
		#else
            
        if (processors == 1) { //save time
			//fastaFilePos.push_back(0); qfileFilePos.push_back(0);
			//fastaFilePos.push_back(1000); qfileFilePos.push_back(1000);
            lines.push_back(linePair(0, 1000));
            if (qfilename != "") {  qLines.push_back(linePair(0, 1000)); }
        }else{
            int numFastaSeqs = 0;
            fastaFilePos = m->setFilePosFasta(filename, numFastaSeqs); 
            if (fastaFilePos.size() < processors) { processors = fastaFilePos.size(); }
        
            if (qfilename != "") { 
                int numQualSeqs = 0;
                qfileFilePos = m->setFilePosFasta(qfilename, numQualSeqs); 
                
                if (numFastaSeqs != numQualSeqs) {
                    m->mothurOut("[ERROR]: You have " + toString(numFastaSeqs) + " sequences in your fasta file, but " + toString(numQualSeqs) + " sequences in your quality file."); m->mothurOutEndLine(); m->control_pressed = true; 
                }
            }
        
            //figure out how many sequences you have to process
            int numSeqsPerProcessor = numFastaSeqs / processors;
            for (int i = 0; i < processors; i++) {
                int startIndex =  i * numSeqsPerProcessor;
                if(i == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor; 	}
                lines.push_back(linePair(fastaFilePos[startIndex], numSeqsPerProcessor));
                if (qfilename != "") {  qLines.push_back(linePair(qfileFilePos[startIndex], numSeqsPerProcessor)); }
            }
        }
            if(qfilename == "")	{	qLines = lines;	} //files with duds
			return 1;
		
		#endif
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "setLines");
		exit(1);
	}
}

//***************************************************************************************************************

bool TrimSeqsCommand::getOligos(vector<vector<string> >& fastaFileNames, vector<vector<string> >& qualFileNames, vector<vector<string> >& nameFileNames){
	try {
		ifstream inOligos;
		m->openInputFile(oligoFile, inOligos);
		
		ofstream test;
		
		string type, oligo, roligo, group;
        bool hasPrimer = false; bool hasPairedBarcodes = false;

		int indexPrimer = 0;
		int indexBarcode = 0;
        int indexPairedPrimer = 0;
		int indexPairedBarcode = 0;
        set<string> uniquePrimers;
        set<string> uniqueBarcodes;
		
		while(!inOligos.eof()){

			inOligos >> type; 
            
		 	if (m->debug) { m->mothurOut("[DEBUG]: reading type - " + type + ".\n"); }	
            
			if(type[0] == '#'){
				while (!inOligos.eof())	{	char c = inOligos.get();  if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
				m->gobble(inOligos);
			}
			else{
				m->gobble(inOligos);
				//make type case insensitive
				for(int i=0;i<type.length();i++){	type[i] = toupper(type[i]);  }
				
				inOligos >> oligo;
                
                if (m->debug) { m->mothurOut("[DEBUG]: reading - " + oligo + ".\n"); }
				
				for(int i=0;i<oligo.length();i++){
					oligo[i] = toupper(oligo[i]);
					if(oligo[i] == 'U')	{	oligo[i] = 'T';	}
				}
				
				if(type == "FORWARD"){
					group = "";
					
					// get rest of line in case there is a primer name
					while (!inOligos.eof())	{	
						char c = inOligos.get(); 
						if (c == 10 || c == 13 || c == -1){	break;	}
						else if (c == 32 || c == 9){;} //space or tab
						else { 	group += c;  }
					} 
					
					//check for repeat barcodes
					map<string, int>::iterator itPrime = primers.find(oligo);
					if (itPrime != primers.end()) { m->mothurOut("primer " + oligo + " is in your oligos file already."); m->mothurOutEndLine();  }
					
                    if (m->debug) {  if (group != "") { m->mothurOut("[DEBUG]: reading group " + group + ".\n"); }else{ m->mothurOut("[DEBUG]: no group for primer " + oligo + ".\n"); }  }
                    
					primers[oligo]=indexPrimer; indexPrimer++;		
					primerNameVector.push_back(group);
				}
                else if (type == "PRIMER"){
                    m->gobble(inOligos);
					
                    inOligos >> roligo;
                    
                    for(int i=0;i<roligo.length();i++){
                        roligo[i] = toupper(roligo[i]);
                        if(roligo[i] == 'U')	{	roligo[i] = 'T';	}
                    }
                    roligo = reverseOligo(roligo);
                    
                    group = "";
                    
					// get rest of line in case there is a primer name
					while (!inOligos.eof())	{
						char c = inOligos.get();
						if (c == 10 || c == 13 || c == -1){	break;	}
						else if (c == 32 || c == 9){;} //space or tab
						else { 	group += c;  }
					}
                    
                    oligosPair newPrimer(oligo, roligo);
                    
                     if (m->debug) { m->mothurOut("[DEBUG]: primer pair " + newPrimer.forward + " " + newPrimer.reverse + ", and group = " + group + ".\n"); }
					
					//check for repeat barcodes
                    string tempPair = oligo+roligo;
                    if (uniquePrimers.count(tempPair) != 0) { m->mothurOut("primer pair " + newPrimer.forward + " " + newPrimer.reverse + " is in your oligos file already."); m->mothurOutEndLine();  }
                    else { uniquePrimers.insert(tempPair); }
					
                    if (m->debug) {  if (group != "") { m->mothurOut("[DEBUG]: reading group " + group + ".\n"); }else{ m->mothurOut("[DEBUG]: no group for primer pair " + newPrimer.forward + " " + newPrimer.reverse + ".\n"); }  }
                    
					pairedPrimers[indexPairedPrimer]=newPrimer; indexPairedPrimer++;
					primerNameVector.push_back(group);
                    hasPrimer = true;
                }
				else if(type == "REVERSE"){
					//Sequence oligoRC("reverse", oligo);
					//oligoRC.reverseComplement();
                    string oligoRC = reverseOligo(oligo);
					revPrimer.push_back(oligoRC);
				}
				else if(type == "BARCODE"){
					inOligos >> group;
                    
                    //barcode lines can look like   BARCODE   atgcatgc   groupName  - for 454 seqs
                    //or                            BARCODE   atgcatgc   atgcatgc    groupName  - for illumina data that has forward and reverse info
                    
                    string temp = "";
                    while (!inOligos.eof())	{
						char c = inOligos.get();
						if (c == 10 || c == 13 || c == -1){	break;	}
						else if (c == 32 || c == 9){;} //space or tab
						else { 	temp += c;  }
					}
					
                    //then this is illumina data with 4 columns
                    if (temp != "") {
                        hasPairedBarcodes = true;
                        string reverseBarcode = group; //reverseOligo(group); //reverse barcode
                        group = temp;
                        
                        for(int i=0;i<reverseBarcode.length();i++){
                            reverseBarcode[i] = toupper(reverseBarcode[i]);
                            if(reverseBarcode[i] == 'U')	{	reverseBarcode[i] = 'T';	}
                        }
                        
                        reverseBarcode = reverseOligo(reverseBarcode);
                        oligosPair newPair(oligo, reverseBarcode);
                        
                        if (m->debug) { m->mothurOut("[DEBUG]: barcode pair " + newPair.forward + " " + newPair.reverse + ", and group = " + group + ".\n"); }
                        
                        //check for repeat barcodes
                        string tempPair = oligo+reverseBarcode;
                        if (uniqueBarcodes.count(tempPair) != 0) { m->mothurOut("barcode pair " + newPair.forward + " " + newPair.reverse +  " is in your oligos file already, disregarding."); m->mothurOutEndLine();  }
                        else { uniqueBarcodes.insert(tempPair); }
                        
                        pairedBarcodes[indexPairedBarcode]=newPair; indexPairedBarcode++;
                        barcodeNameVector.push_back(group);
                    }else {				
                        //check for repeat barcodes
                        map<string, int>::iterator itBar = barcodes.find(oligo);
                        if (itBar != barcodes.end()) { m->mothurOut("barcode " + oligo + " is in your oligos file already."); m->mothurOutEndLine();  }
                        
                        barcodes[oligo]=indexBarcode; indexBarcode++;
                        barcodeNameVector.push_back(group);
                    }
				}else if(type == "LINKER"){
					linker.push_back(oligo);
				}else if(type == "SPACER"){
					spacer.push_back(oligo);
				}
				else{	m->mothurOut("[WARNING]: " + type + " is not recognized as a valid type. Choices are forward, reverse, and barcode. Ignoring " + oligo + "."); m->mothurOutEndLine(); }
			}
			m->gobble(inOligos);
		}	
		inOligos.close();
		
        if (hasPairedBarcodes || hasPrimer) {
            pairedOligos = true;
            if ((primers.size() != 0) || (barcodes.size() != 0) || (linker.size() != 0) || (spacer.size() != 0) || (revPrimer.size() != 0)) { m->control_pressed = true;  m->mothurOut("[ERROR]: cannot mix paired primers and barcodes with non paired or linkers and spacers, quitting."); m->mothurOutEndLine();  return 0; }
        }
        
		if(barcodeNameVector.size() == 0 && primerNameVector[0] == ""){	allFiles = 0;	}
        
		//add in potential combos
		if(barcodeNameVector.size() == 0){
			barcodes[""] = 0;
			barcodeNameVector.push_back("");			
		}
		
		if(primerNameVector.size() == 0){
			primers[""] = 0;
			primerNameVector.push_back("");			
		}
		
		fastaFileNames.resize(barcodeNameVector.size());
		for(int i=0;i<fastaFileNames.size();i++){
			fastaFileNames[i].assign(primerNameVector.size(), "");
		}
		if(qFileName != "")	{	qualFileNames = fastaFileNames;	}
		if(nameFile != "")	{	nameFileNames = fastaFileNames;	}
		
		if(allFiles){
			set<string> uniqueNames; //used to cleanup outputFileNames
            if (pairedOligos) {
                for(map<int, oligosPair>::iterator itBar = pairedBarcodes.begin();itBar != pairedBarcodes.end();itBar++){
                    for(map<int, oligosPair>::iterator itPrimer = pairedPrimers.begin();itPrimer != pairedPrimers.end(); itPrimer++){
                        
                        string primerName = primerNameVector[itPrimer->first];
                        string barcodeName = barcodeNameVector[itBar->first];
                        
                        if ((primerName == "ignore") || (barcodeName == "ignore")) { } //do nothing
                        else {
                            string comboGroupName = "";
                            string fastaFileName = "";
                            string qualFileName = "";
                            string nameFileName = "";
                            string countFileName = "";
                            
                            if(primerName == ""){
                                comboGroupName = barcodeNameVector[itBar->first];
                            }
                            else{
                                if(barcodeName == ""){
                                    comboGroupName = primerNameVector[itPrimer->first];
                                }
                                else{
                                    comboGroupName = barcodeNameVector[itBar->first] + "." + primerNameVector[itPrimer->first];
                                }
                            }
                            
                            
                            ofstream temp;
                            map<string, string> variables;
                            variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastaFile));
                            variables["[tag]"] = comboGroupName;
                            fastaFileName = getOutputFileName("fasta", variables);
                            if (uniqueNames.count(fastaFileName) == 0) {
                                outputNames.push_back(fastaFileName);
                                outputTypes["fasta"].push_back(fastaFileName);
                                uniqueNames.insert(fastaFileName);
                            }
                            
                            fastaFileNames[itBar->first][itPrimer->first] = fastaFileName;
                            m->openOutputFile(fastaFileName, temp);		temp.close();
                            
                            if(qFileName != ""){
                                variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(qFileName));
                                qualFileName = getOutputFileName("qfile", variables);
                                if (uniqueNames.count(qualFileName) == 0) {
                                    outputNames.push_back(qualFileName);
                                    outputTypes["qfile"].push_back(qualFileName);
                                }
                                
                                qualFileNames[itBar->first][itPrimer->first] = qualFileName;
                                m->openOutputFile(qualFileName, temp);		temp.close();
                            }
                            
                            if(nameFile != ""){
                                variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(nameFile));
                                nameFileName = getOutputFileName("name", variables);
                                if (uniqueNames.count(nameFileName) == 0) {
                                    outputNames.push_back(nameFileName);
                                    outputTypes["name"].push_back(nameFileName);
                                }
                                
                                nameFileNames[itBar->first][itPrimer->first] = nameFileName;
                                m->openOutputFile(nameFileName, temp);		temp.close();
                            }
                        }
                    }
                }
            }else {
                for(map<string, int>::iterator itBar = barcodes.begin();itBar != barcodes.end();itBar++){
                    for(map<string, int>::iterator itPrimer = primers.begin();itPrimer != primers.end(); itPrimer++){
                        
                        string primerName = primerNameVector[itPrimer->second];
                        string barcodeName = barcodeNameVector[itBar->second];
                        
                        if ((primerName == "ignore") || (barcodeName == "ignore")) { } //do nothing 
                        else {
                            string comboGroupName = "";
                            string fastaFileName = "";
                            string qualFileName = "";
                            string nameFileName = "";
                            string countFileName = "";
                            
                            if(primerName == ""){
                                comboGroupName = barcodeNameVector[itBar->second];
                            }
                            else{
                                if(barcodeName == ""){
                                    comboGroupName = primerNameVector[itPrimer->second];
                                }
                                else{
                                    comboGroupName = barcodeNameVector[itBar->second] + "." + primerNameVector[itPrimer->second];
                                }
                            }
                            
                            
                            ofstream temp;
                            map<string, string> variables; 
                            variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastaFile));
                            variables["[tag]"] = comboGroupName;
                            fastaFileName = getOutputFileName("fasta", variables);
                            if (uniqueNames.count(fastaFileName) == 0) {
                                outputNames.push_back(fastaFileName);
                                outputTypes["fasta"].push_back(fastaFileName);
                                uniqueNames.insert(fastaFileName);
                            }
                            
                            fastaFileNames[itBar->second][itPrimer->second] = fastaFileName;
                            m->openOutputFile(fastaFileName, temp);		temp.close();
                            
                            if(qFileName != ""){
                                variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(qFileName));
                                qualFileName = getOutputFileName("qfile", variables);
                                if (uniqueNames.count(qualFileName) == 0) {
                                    outputNames.push_back(qualFileName);
                                    outputTypes["qfile"].push_back(qualFileName);
                                }
                                
                                qualFileNames[itBar->second][itPrimer->second] = qualFileName;
                                m->openOutputFile(qualFileName, temp);		temp.close();
                            }
                            
                            if(nameFile != ""){
                                variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(nameFile));
                                nameFileName = getOutputFileName("name", variables);
                                if (uniqueNames.count(nameFileName) == 0) {
                                    outputNames.push_back(nameFileName);
                                    outputTypes["name"].push_back(nameFileName);
                                }
                                
                                nameFileNames[itBar->second][itPrimer->second] = nameFileName;
                                m->openOutputFile(nameFileName, temp);		temp.close();
                            }
                        }
                    }
                }
            }
		}
		numFPrimers = primers.size();
        if (pairedOligos) { numFPrimers  = pairedPrimers.size(); }
		numRPrimers = revPrimer.size();
        numLinkers = linker.size();
        numSpacers = spacer.size();
		
		bool allBlank = true;
		for (int i = 0; i < barcodeNameVector.size(); i++) {
			if (barcodeNameVector[i] != "") {
				allBlank = false;
				break;
			}
		}
		for (int i = 0; i < primerNameVector.size(); i++) {
			if (primerNameVector[i] != "") {
				allBlank = false;
				break;
			}
		}

		if (allBlank) {
			m->mothurOut("[WARNING]: your oligos file does not contain any group names.  mothur will not create a groupfile."); m->mothurOutEndLine();
			allFiles = false;
			return false;
		}
		
		return true;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "getOligos");
		exit(1);
	}
}
//***************************************************************************************************************

bool TrimSeqsCommand::keepFirstTrim(Sequence& sequence, QualityScores& qscores){
	try {
		bool success = 1;
		if(qscores.getName() != ""){
			qscores.trimQScores(-1, keepFirst);
		}

//        sequence.printSequence(cout);cout << endl;
        
		sequence.trim(keepFirst);
        
//        sequence.printSequence(cout);cout << endl << endl;;

		return success;
	}
	catch(exception& e) {
		m->errorOut(e, "keepFirstTrim", "countDiffs");
		exit(1);
	}
	
}	

//***************************************************************************************************************

bool TrimSeqsCommand::removeLastTrim(Sequence& sequence, QualityScores& qscores){
	try {
		bool success = 0;
		
		int length = sequence.getNumBases() - removeLast;
		
		if(length > 0){
			if(qscores.getName() != ""){
				qscores.trimQScores(-1, length);
			}
			sequence.trim(length);
			success = 1;
		}
		else{
			success = 0;
		}

		return success;
	}
	catch(exception& e) {
		m->errorOut(e, "removeLastTrim", "countDiffs");
		exit(1);
	}
	
}	

//***************************************************************************************************************

bool TrimSeqsCommand::cullLength(Sequence& seq){
	try {
	
		int length = seq.getNumBases();
		bool success = 0;	//guilty until proven innocent
		
		if(length >= minLength && maxLength == 0)			{	success = 1;	}
		else if(length >= minLength && length <= maxLength)	{	success = 1;	}
		else												{	success = 0;	}
		
		return success;
	
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "cullLength");
		exit(1);
	}
	
}

//***************************************************************************************************************

bool TrimSeqsCommand::cullHomoP(Sequence& seq){
	try {
		int longHomoP = seq.getLongHomoPolymer();
		bool success = 0;	//guilty until proven innocent
		
		if(longHomoP <= maxHomoP){	success = 1;	}
		else					{	success = 0;	}
		
		return success;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "cullHomoP");
		exit(1);
	}
	
}
//********************************************************************/
string TrimSeqsCommand::reverseOligo(string oligo){
	try {
        string reverse = "";
        
        for(int i=oligo.length()-1;i>=0;i--){
            
            if(oligo[i] == 'A')		{	reverse += 'T';	}
            else if(oligo[i] == 'T'){	reverse += 'A';	}
            else if(oligo[i] == 'U'){	reverse += 'A';	}
            
            else if(oligo[i] == 'G'){	reverse += 'C';	}
            else if(oligo[i] == 'C'){	reverse += 'G';	}
            
            else if(oligo[i] == 'R'){	reverse += 'Y';	}
            else if(oligo[i] == 'Y'){	reverse += 'R';	}
            
            else if(oligo[i] == 'M'){	reverse += 'K';	}
            else if(oligo[i] == 'K'){	reverse += 'M';	}
            
            else if(oligo[i] == 'W'){	reverse += 'W';	}
            else if(oligo[i] == 'S'){	reverse += 'S';	}
            
            else if(oligo[i] == 'B'){	reverse += 'V';	}
            else if(oligo[i] == 'V'){	reverse += 'B';	}
            
            else if(oligo[i] == 'D'){	reverse += 'H';	}
            else if(oligo[i] == 'H'){	reverse += 'D';	}
            
            else						{	reverse += 'N';	}
        }
        
        
        return reverse;
    }
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "reverseOligo");
		exit(1);
	}
}

//***************************************************************************************************************

bool TrimSeqsCommand::cullAmbigs(Sequence& seq){
	try {
		int numNs = seq.getAmbigBases();
		bool success = 0;	//guilty until proven innocent
		
		if(numNs <= maxAmbig)	{	success = 1;	}
		else					{	success = 0;	}
		
		return success;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "cullAmbigs");
		exit(1);
	}
	
}
//***************************************************************************************************************

