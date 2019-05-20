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
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
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
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["qfile"] = tempOutNames;
			outputTypes["group"] = tempOutNames;
			outputTypes["name"] = tempOutNames;
            outputTypes["count"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("oligos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["oligos"] = inputDir + it->second;		}
				}
				
				it = parameters.find("qfile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["qfile"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
				
			}

			
			//check for required parameters
			fastaFile = validParameter.validFile(parameters, "fasta");
			if (fastaFile == "not found") { 				
				fastaFile = current->getFastaFile(); 
				if (fastaFile != "") { m->mothurOut("Using " + fastaFile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else if (fastaFile == "not open") { abort = true; }	
			else { current->setFastaFile(fastaFile); }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += util.hasPath(fastaFile); //if user entered a file with a path then preserve it	
			}
		
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp = validParameter.valid(parameters, "flip");
			if (temp == "not found")    {	flip = 0;	}
			else {  flip = util.isTrue(temp);		}
		
			temp = validParameter.validFile(parameters, "oligos");
			if (temp == "not found"){	oligoFile = "";		}
			else if(temp == "not open"){	abort = true;	} 
			else					{	oligoFile = temp; current->setOligosFile(oligoFile);		}
			
			
			temp = validParameter.valid(parameters, "maxambig");		if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, maxAmbig);  

			temp = validParameter.valid(parameters, "maxhomop");		if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, maxHomoP);  

			temp = validParameter.valid(parameters, "minlength");	if (temp == "not found") { temp = "1"; }
			util.mothurConvert(temp, minLength); 
			
			temp = validParameter.valid(parameters, "maxlength");	if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, maxLength);
			
			temp = validParameter.valid(parameters, "bdiffs");		if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, bdiffs);
			
			temp = validParameter.valid(parameters, "pdiffs");		if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, pdiffs);
            
            temp = validParameter.valid(parameters, "ldiffs");		if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, ldiffs);
            
            temp = validParameter.valid(parameters, "sdiffs");		if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, sdiffs);
			
			temp = validParameter.valid(parameters, "tdiffs");		if (temp == "not found") { int tempTotal = pdiffs + bdiffs + ldiffs + sdiffs;  temp = toString(tempTotal); }
			util.mothurConvert(temp, tdiffs);
			
			if(tdiffs == 0){	tdiffs = bdiffs + pdiffs + ldiffs + sdiffs;	}
			
			temp = validParameter.validFile(parameters, "qfile");	
			if (temp == "not found")	{	qFileName = "";		}
			else if(temp == "not open")	{	abort = true;		}
			else						{	qFileName = temp;	current->setQualFile(qFileName); }
			
			temp = validParameter.validFile(parameters, "name");	
			if (temp == "not found")	{	nameFile = "";		}
			else if(temp == "not open")	{	nameFile = "";	abort = true;		}
			else						{	nameFile = temp;	current->setNameFile(nameFile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { abort = true; countfile = ""; }	
			else if (countfile == "not found") { countfile = ""; }
			else { current->setCountFile(countfile); }
			
            if ((countfile != "") && (nameFile != "")) { m->mothurOut("You must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
			
			temp = validParameter.valid(parameters, "qthreshold");	if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, qThreshold);
			
			temp = validParameter.valid(parameters, "qtrim");		if (temp == "not found") { temp = "t"; }
			qtrim = util.isTrue(temp);

			temp = validParameter.valid(parameters, "rollaverage");	if (temp == "not found") { temp = "0"; }
			convert(temp, qRollAverage);

			temp = validParameter.valid(parameters, "qwindowaverage");if (temp == "not found") { temp = "0"; }
			convert(temp, qWindowAverage);

			temp = validParameter.valid(parameters, "qwindowsize");	if (temp == "not found") { temp = "50"; }
			convert(temp, qWindowSize);

			temp = validParameter.valid(parameters, "qstepsize");	if (temp == "not found") { temp = "1"; }
			convert(temp, qWindowStep);

			temp = validParameter.valid(parameters, "qaverage");		if (temp == "not found") { temp = "0"; }
			convert(temp, qAverage);

			temp = validParameter.valid(parameters, "keepfirst");	if (temp == "not found") { temp = "0"; }
			convert(temp, keepFirst);

			temp = validParameter.valid(parameters, "removelast");	if (temp == "not found") { temp = "0"; }
			convert(temp, removeLast);
			
			temp = validParameter.valid(parameters, "allfiles");		if (temp == "not found") { temp = "F"; }
			allFiles = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "keepforward");		if (temp == "not found") { temp = "F"; }
			keepforward = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "logtransform");		if (temp == "not found") { temp = "F"; }
			logtransform = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "checkorient");		if (temp == "not found") { temp = "F"; }
			reorient = util.isTrue(temp);
			
			temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
			
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
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
                }
            }
		}
        
        pairedOligos = false;
        createGroup = false;

	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "TrimSeqsCommand");
		exit(1);
	}
}
//***************************************************************************************************************

int TrimSeqsCommand::execute(){
	try{
	
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
        map<string, string> variables; 
		variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastaFile));
        variables["[tag]"] = "trim";
		string trimSeqFile = getOutputFileName("fasta",variables);
        string trimQualFile = getOutputFileName("qfile",variables);
		outputNames.push_back(trimSeqFile); outputTypes["fasta"].push_back(trimSeqFile);
        
        variables["[tag]"] = "scrap";
		string scrapSeqFile = getOutputFileName("fasta",variables);
        string scrapQualFile = getOutputFileName("qfile",variables);
		outputNames.push_back(scrapSeqFile); outputTypes["fasta"].push_back(scrapSeqFile);
		
		if (qFileName != "") {
			outputNames.push_back(trimQualFile); outputNames.push_back(scrapQualFile);
			outputTypes["qfile"].push_back(trimQualFile); outputTypes["qfile"].push_back(scrapQualFile);
		}
		
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(nameFile));
        variables["[tag]"] = "trim";
		string trimNameFile = getOutputFileName("name",variables);
        variables["[tag]"] = "scrap";
		string scrapNameFile = getOutputFileName("name",variables);
		
		if (nameFile != "") {
			util.readNames(nameFile, nameMap);
			outputNames.push_back(trimNameFile); outputNames.push_back(scrapNameFile);
			outputTypes["name"].push_back(trimNameFile); outputTypes["name"].push_back(scrapNameFile);
		}
        
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(countfile));
        variables["[tag]"] = "trim";
        string trimCountFile = getOutputFileName("count",variables);
        variables["[tag]"] = "scrap";
		string scrapCountFile = getOutputFileName("count",variables);
		
		if (countfile != "") {
            CountTable ct;
            ct.readTable(countfile, true, false);
            nameCount = ct.getNameMap();
			outputNames.push_back(trimCountFile); outputNames.push_back(scrapCountFile);
			outputTypes["count"].push_back(trimCountFile); outputTypes["count"].push_back(scrapCountFile);
		}

		if (m->getControl_pressed()) { return 0; }
		
        int startTime = time(NULL);
        
        long long numSeqs = createProcessesCreateTrim(fastaFile, qFileName, trimSeqFile, scrapSeqFile, trimQualFile, scrapQualFile, trimNameFile, scrapNameFile, trimCountFile, scrapCountFile);
        
        m->mothurOut("It took " + toString(time(NULL) - startTime) + " secs to trim " + toString(numSeqs) + " sequences.\n");
		
		if (m->getControl_pressed()) {	for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0;	}

		//output group counts
		m->mothurOutEndLine();
		int total = 0;
		if (groupCounts.size() != 0) {  m->mothurOut("Group count: \n");  }
		for (map<string, int>::iterator it = groupCounts.begin(); it != groupCounts.end(); it++) {
			 total += it->second; m->mothurOut(it->first + "\t" + toString(it->second)); m->mothurOutEndLine(); 
		}
		if (total != 0) { m->mothurOut("Total of all groups is " + toString(total)); m->mothurOutEndLine(); }
		
		if (m->getControl_pressed()) {	for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0;	}

		//set fasta file as new current fastafile
		string currentName = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
		}
		
		itTypes = outputTypes.find("name");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setNameFile(currentName); }
		}
		
		itTypes = outputTypes.find("qfile");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setQualFile(currentName); }
		}
		
		itTypes = outputTypes.find("group");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setGroupFile(currentName); }
		}
        
        itTypes = outputTypes.find("count");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
		}

		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
		
		return 0;	
			
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "execute");
		exit(1);
	}
}
//***************************************************************************************************************

bool keepFirstTrim(Sequence& sequence, QualityScores& qscores, int keepFirst){
    bool success = 1;
    if(qscores.getName() != ""){ qscores.trimQScores(-1, keepFirst); }
        
    sequence.trim(keepFirst);
        
    return success;
}
//***************************************************************************************************************

bool removeLastTrim(Sequence& sequence, QualityScores& qscores, int removeLast){
    bool success = 0;
        
    int length = sequence.getNumBases() - removeLast;
        
    if(length > 0){
        if(qscores.getName() != ""){ qscores.trimQScores(-1, length); }
        sequence.trim(length);
        success = 1;
    }
    else{ success = 0; }
        
    return success;
}
//***************************************************************************************************************
bool cullLength(Sequence& seq, int minLength, int maxLength){
    int length = seq.getNumBases();
    bool success = 0;	//guilty until proven innocent
        
    if(length >= minLength && maxLength == 0)			{	success = 1;	}
    else if(length >= minLength && length <= maxLength)	{	success = 1;	}
    else												{	success = 0;	}
        
    return success;
}
//***************************************************************************************************************
bool cullHomoP(Sequence& seq, int maxHomoP){
    int longHomoP = seq.getLongHomoPolymer();
    bool success = 0;	//guilty until proven innocent
        
    if(longHomoP <= maxHomoP){	success = 1;	}
    else					{	success = 0;	}
        
    return success;
}
//***************************************************************************************************************
bool cullAmbigs(Sequence& seq, int maxAmbig){
    int numNs = seq.getAmbigBases();
    bool success = 0;	//guilty until proven innocent
        
    if(numNs <= maxAmbig)	{	success = 1;	}
    else					{	success = 0;	}
        
    return success;
}
/**************************************************************************************************/
struct trimData {
    unsigned long long start, end;
    MothurOut* m;
    Oligos oligos;
    string filename, qFileName;
    OutputWriter* trimFileName;
    OutputWriter* scrapFileName;
    OutputWriter* trimQFileName;
    OutputWriter* scrapQFileName;
    set<string> badNames;
    unsigned long long lineStart, lineEnd, qlineStart, qlineEnd;
    bool flip, allFiles, qtrim, keepforward, createGroup, pairedOligos, reorient, logtransform;
    int maxAmbig, maxHomoP, minLength, maxLength, tdiffs, bdiffs, pdiffs, ldiffs, sdiffs;
    int qWindowSize, qWindowStep, keepFirst, removeLast, count;
    double qRollAverage, qThreshold, qWindowAverage, qAverage;
    vector<string> revPrimer;
    map<string, int> barcodes;
    map<string, int> primers;
    map<string, int> nameCount;
    vector<string>  linker;
    vector<string>  spacer;
    map<string, int> combos;
    map<string, int> groupCounts;
    map<string, string> nameMap;
    map<string, string> groupMap;
    map<int, oligosPair> pairedBarcodes;
    map<int, oligosPair> pairedPrimers;
    Utils util;
    
    trimData(){}
    trimData(string fn, string qn, OutputWriter* tn, OutputWriter* sn, OutputWriter* tqn, OutputWriter* sqn, unsigned long long lstart, unsigned long long lend, unsigned long long qstart, unsigned long long qend, map<string, string> nm, map<string, int> ncount) {
        filename = fn;
        qFileName = qn;
        trimFileName = tn;
        scrapFileName = sn;
        trimQFileName = tqn;
        scrapQFileName = sqn;
        lineStart = lstart;
        lineEnd = lend;
        qlineStart = qstart;
        qlineEnd = qend;
        m = MothurOut::getInstance();
        nameCount = ncount;
        nameMap = nm;

    }
    void setOligosOptions(Oligos olig, int pd, int bd, int ld, int sd, int td, bool po, bool cGroup, bool aFiles, bool keepF, int keepfi, int removeL,
                          int WindowStep, int WindowSize, int WindowAverage, bool trim, double Threshold, double Average, double RollAverage, bool lt,
                          int minL, int maxA, int maxH, int maxL, bool fli, bool reo) {
        oligos = olig;
        pdiffs = pd;
        bdiffs = bd;
        ldiffs = ld;
        sdiffs = sd;
        tdiffs = td;
        barcodes = oligos.getBarcodes();
        pairedPrimers = oligos.getPairedPrimers();
        pairedBarcodes = oligos.getPairedBarcodes();
        pairedOligos = po;
        primers = oligos.getPrimers();
        revPrimer = oligos.getReversePrimers();
        linker = oligos.getLinkers();
        spacer = oligos.getSpacers();
        
        createGroup = cGroup;
        allFiles = aFiles;
        keepforward = keepF;
        keepFirst = keepfi;
        removeLast = removeL;
        qWindowStep = WindowStep;
        qWindowSize = WindowSize;
        qWindowAverage = WindowAverage;
        qtrim = trim;
        qThreshold = Threshold;
        qAverage = Average;
        qRollAverage = RollAverage;
        logtransform = lt;
        minLength = minL;
        maxAmbig = maxA;
        maxHomoP = maxH;
        maxLength = maxL;
        flip = fli;
        reorient = reo;
        count = 0;
    }
};
/**************************************************************************************/
//string filename, string qFileName, string trimFileName, string scrapFileName, string trimQFileName, string scrapQFileName, string trimNFileName, string scrapNFileName, string trimCFileName, string scrapCFileName, string groupFileName, vector<vector<string> > fastaFileNames, vector<vector<string> > qualFileNames, vector<vector<string> > nameFileNames, linePair line, linePair qline
int driverTrim(trimData* params) {
		
	try {
        int numFPrimers, numRPrimers, numLinkers, numSpacers;
        numFPrimers = params->primers.size();
        numRPrimers = params->revPrimer.size();
        numLinkers = params->linker.size();
        numSpacers = params->spacer.size();
		
		ifstream inFASTA;
		params->util.openInputFile(params->filename, inFASTA);
		inFASTA.seekg(params->lineStart);
		
		ifstream qFile;
		if(params->qFileName != "")	{
			params->util.openInputFile(params->qFileName, qFile);
			qFile.seekg(params->qlineStart);
		}
		
		bool moreSeqs = 1;
        int numBarcodes = params->barcodes.size();
		TrimOligos* trimOligos = NULL;
        if (params->pairedOligos)   {   trimOligos = new TrimOligos(params->pdiffs, params->bdiffs, 0, 0, params->pairedPrimers, params->pairedBarcodes, false); numBarcodes = params->pairedBarcodes.size(); }
        else                {   trimOligos = new TrimOligos(params->pdiffs, params->bdiffs, params->ldiffs, params->sdiffs, params->primers, params->barcodes, params->revPrimer, params->linker, params->spacer);  }
        
        TrimOligos* rtrimOligos = NULL;
        if (params->reorient) {
            //create reoriented primer and barcode pairs
            map<int, oligosPair> rpairedPrimers, rpairedBarcodes;
            for (map<int, oligosPair>::iterator it = params->pairedPrimers.begin(); it != params->pairedPrimers.end(); it++) {
                  oligosPair tempPair(params->util.reverseOligo((it->second).reverse), (params->util.reverseOligo((it->second).forward))); //reversePrimer, rc ForwardPrimer
                rpairedPrimers[it->first] = tempPair;
            }
            for (map<int, oligosPair>::iterator it = params->pairedBarcodes.begin(); it != params->pairedBarcodes.end(); it++) {
                 oligosPair tempPair(params->util.reverseOligo((it->second).reverse), (params->util.reverseOligo((it->second).forward))); //reverseBarcode, rc ForwardBarcode
                rpairedBarcodes[it->first] = tempPair;
            }
            int index = rpairedBarcodes.size();
            for (map<string, int>::iterator it = params->barcodes.begin(); it != params->barcodes.end(); it++) {
                oligosPair tempPair("", params->util.reverseOligo((it->first))); //reverseBarcode, rc ForwardBarcode
                rpairedBarcodes[index] = tempPair; index++;
            }
            
            index = rpairedPrimers.size();
            for (map<string, int>::iterator it = params->primers.begin(); it != params->primers.end(); it++) {
                oligosPair tempPair("", params->util.reverseOligo((it->first))); //reverseBarcode, rc ForwardBarcode
                rpairedPrimers[index] = tempPair; index++;
            }

            rtrimOligos = new TrimOligos(params->pdiffs, params->bdiffs, 0, 0, rpairedPrimers, rpairedBarcodes, false); numBarcodes = rpairedBarcodes.size();
        }
        
        if(numBarcodes == 0){ params->createGroup = false; }
        
		while (moreSeqs) {
				
            int obsBDiffs = 0;
            int obsPDiffs = 0;
            
			if (params->m->getControl_pressed()) { break; }
			
			int success = 1;
			string trashCode = "";
            string commentString = "";
			int currentSeqsDiffs = 0;

			Sequence currSeq(inFASTA); params->util.gobble(inFASTA);
            Sequence savedSeq(currSeq.getName(), currSeq.getAligned());
            
			QualityScores currQual; QualityScores savedQual;
			if(params->qFileName != ""){
				currQual = QualityScores(qFile);  params->util.gobble(qFile);
                savedQual.setName(currQual.getName()); savedQual.setScores(currQual.getScores());
			}
			  
			string origSeq = currSeq.getUnaligned();
			if (origSeq != "") {
				
				int barcodeIndex = 0;
				int primerIndex = 0;
               
                if(numLinkers != 0){
					success = trimOligos->stripLinker(currSeq, currQual);
					if(success > params->ldiffs)		{	trashCode += 'k';	}
					else{ currentSeqsDiffs += success;  }

				}
                
				if(numBarcodes != 0){
					vector<int> results = trimOligos->stripBarcode(currSeq, currQual, barcodeIndex);
                    if (params->pairedOligos) {
                        success = results[0] + results[2];
                        commentString += "fbdiffs=" + toString(results[0]) + "(" + trimOligos->getCodeValue(results[1], params->bdiffs) + "), rbdiffs=" + toString(results[2]) + "(" + trimOligos->getCodeValue(results[3], params->bdiffs) + ") ";
                    }
                    else {
                        success = results[0];
                        commentString += "bdiffs=" + toString(results[0]) + "(" + trimOligos->getCodeValue(results[1], params->bdiffs) + ") ";
                    }
					if(success > params->bdiffs)		{  trashCode += 'b'; }
					else{ currentSeqsDiffs += success;  }
				}
                obsBDiffs = success;
                
               if(numSpacers != 0){
					success = trimOligos->stripSpacer(currSeq, currQual);
					if(success > params->sdiffs)		{	trashCode += 's';	}
					else{ currentSeqsDiffs += success;  }

				}
                
				if(numFPrimers != 0){
					vector<int> results = trimOligos->stripForward(currSeq, currQual, primerIndex, params->keepforward);
                    if (params->pairedOligos) {
                        success = results[0] + results[2];
                        commentString += "fpdiffs=" + toString(results[0]) + "(" + trimOligos->getCodeValue(results[1], params->pdiffs) + "), rpdiffs=" + toString(results[2]) + "(" + trimOligos->getCodeValue(results[3], params->pdiffs) + ") ";
                    }
                    else {
                        success = results[0];
                        commentString += "fpdiffs=" + toString(results[0]) + "(" + trimOligos->getCodeValue(results[1], params->pdiffs) + ") ";
                    }
					if(success > params->pdiffs)		{  trashCode += 'f';  }
					else{ currentSeqsDiffs += success;  }
				}
                obsPDiffs = success;
                
				if(numRPrimers != 0){
					vector<int> results =  trimOligos->stripReverse(currSeq, currQual);
                    success = results[0];
                    commentString += "rpdiffs=" + toString(results[0]) + "(" + trimOligos->getCodeValue(results[1], params->pdiffs) + ") ";
                    if(success > params->pdiffs)		{	trashCode += 'r';	}
                    else{ currentSeqsDiffs += success;  }
				}
                
                if (currentSeqsDiffs > params->tdiffs)	{	trashCode += 't';   }
                
                if (params->reorient && (trashCode != "")) { //if you failed and want to check the reverse
                    int thisSuccess = 0;
                    string thisTrashCode = "";
                    string thiscommentString = "";
                    int thisCurrentSeqsDiffs = 0;
                    
                    int thisBarcodeIndex = 0;
                    int thisPrimerIndex = 0;
                    
                    if(numBarcodes != 0){
                        vector<int> results =  rtrimOligos->stripBarcode(savedSeq, savedQual, thisBarcodeIndex);
                        if (params->pairedOligos) {
                            thisSuccess = results[0] + results[2];
                            thiscommentString += "fbdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], params->bdiffs) + "), rbdiffs=" + toString(results[2]) + "(" + rtrimOligos->getCodeValue(results[3], params->bdiffs) + ") ";
                        }
                        else {
                            thisSuccess = results[0];
                            thiscommentString += "bdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], params->bdiffs) + ") ";
                        }
                        if(thisSuccess > params->bdiffs)		{ thisTrashCode += "b"; }
                        else{ thisCurrentSeqsDiffs += thisSuccess;  }
                    }
                    
                    int revBDiffs = thisSuccess;
                    
                    if(numFPrimers != 0){
                        vector<int> results =  rtrimOligos->stripForward(savedSeq, savedQual, thisPrimerIndex, params->keepforward);
                        if (params->pairedOligos) {
                            thisSuccess = results[0] + results[2];
                            thiscommentString += "fpdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], params->pdiffs) + "), rpdiffs=" + toString(results[2]) + "(" + rtrimOligos->getCodeValue(results[3], params->pdiffs) + ") ";
                        }
                        else {
                            thisSuccess = results[0];
                            thiscommentString += "pdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], params->pdiffs) + ") ";
                        }
                        if(thisSuccess > params->pdiffs)		{ thisTrashCode += "f"; }
                        else{ thisCurrentSeqsDiffs += thisSuccess;  }
                    }
                    int revPDiffs = thisSuccess;

                    if (thisCurrentSeqsDiffs > params->tdiffs)	{	thisTrashCode += 't';   }
                    
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
                        if(params->qFileName != ""){
                            savedQual.flipQScores();
                            currQual.setScores(savedQual.getScores());
                        }
                    }else { trashCode += "(" + thisTrashCode + ")";  }
                }
                
				if(params->keepFirst != 0){ success = keepFirstTrim(currSeq, currQual, params->keepFirst); }
				
				if(params->removeLast != 0){
					success = removeLastTrim(currSeq, currQual, params->removeLast);
					if(!success)				{	trashCode += 'l';	}
				}
				
				if(params->qFileName != ""){
					int origLength = currSeq.getNumBases();
					
					if(!params->util.isEqual(params->qThreshold,0))			{	success = currQual.stripQualThreshold(currSeq, params->qThreshold);			}
					else if(!params->util.isEqual(params->qAverage, 0))		{	success = currQual.cullQualAverage(currSeq, params->qAverage, params->logtransform);				}
					else if(!params->util.isEqual(params->qRollAverage, 0))	{	success = currQual.stripQualRollingAverage(currSeq, params->qRollAverage, params->logtransform);	}
					else if(!params->util.isEqual(params->qWindowAverage, 0)){	success = currQual.stripQualWindowAverage(currSeq, params->qWindowStep, params->qWindowSize, params->qWindowAverage, params->logtransform);	}
					else						{	success = 1;				}
					
					//you don't want to trim, if it fails above then scrap it
					if ((!params->qtrim) && (origLength != currSeq.getNumBases())) {  success = 0; }
					
					if(!success)				{	trashCode += 'q';	}
				}				
		
				if(params->minLength > 0 || params->maxLength > 0){
					success = cullLength(currSeq, params->minLength, params->maxLength);
					if(!success)				{	trashCode += 'l';	}
				}
				if(params->maxHomoP > 0){
					success = cullHomoP(currSeq, params->maxHomoP);
					if(!success)				{	trashCode += 'h';	}
				}
				if(params->maxAmbig != -1){
					success = cullAmbigs(currSeq, params->maxAmbig);
					if(!success)				{	trashCode += 'n';	}
				}
				
				if(params->flip){		// should go last
					currSeq.reverseComplement();
					if(params->qFileName != ""){ currQual.flipQScores();	 }
				}
				
                if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: " + currSeq.getName() + ", trashcode= " + trashCode + "\n");  }

                string seqComment = currSeq.getComment();
                currSeq.setComment("\t" + commentString + "\t" + seqComment);
                
                if(trashCode.length() == 0){
                    string thisGroup = "";
                    if (params->createGroup) { thisGroup = params->oligos.getGroupName(barcodeIndex, primerIndex); }
                    
                    int pos = thisGroup.find("ignore");
                    if (pos == string::npos) {
                        currSeq.setAligned(currSeq.getUnaligned());
                        currSeq.printSequence(params->trimFileName);
                        
                        if(params->qFileName != ""){ currQual.printQScores(params->trimQFileName); }
                        
                        if (params->createGroup) {
                            if (params->m->getDebug()) { params->m->mothurOut(", group= " + thisGroup + "\n"); }
                            params->groupMap[currSeq.getName()] = thisGroup;
                            
                            map<string, int>::iterator it = params->groupCounts.find(thisGroup);
                            if (it == params->groupCounts.end()) {	params->groupCounts[thisGroup] = 1; }
                            else { params->groupCounts[it->first] += 1; }
                        }
                    }
                }
                else{
                    params->badNames.insert(currSeq.getName());
                    currSeq.setName(currSeq.getName() + " | " + trashCode);
                    currSeq.setUnaligned(origSeq);
                    currSeq.setAligned(origSeq);
                    currSeq.printSequence(params->scrapFileName);
                    if(params->qFileName != ""){ currQual.printQScores(params->scrapQFileName); }
                }
				params->count++;
			}
			
			#if defined NON_WINDOWS
				unsigned long long pos = inFASTA.tellg();
				if ((pos == -1) || (pos >= params->lineEnd)) { break; }
			#else
				if ((params->count == params->lineEnd) || (inFASTA.eof())) { break; }
			#endif
			
			//report progress
			if((params->count) % 1000 == 0){	params->m->mothurOutJustToScreen(toString(params->count)+"\n"); 		}
			
		}
		//report progress
		if((params->count) % 1000 != 0){	params->m->mothurOutJustToScreen(toString(params->count)+"\n");		}
		
		delete trimOligos;
        if (params->reorient) { delete rtrimOligos; }
		inFASTA.close();
		if(params->qFileName != "")	{	qFile.close();		}
        
        return params->count;
    }
	catch(exception& e) {
		params->m->errorOut(e, "TrimSeqsCommand", "driverTrim");
		exit(1);
	}
}

/**************************************************************************************************/
long long TrimSeqsCommand::createProcessesCreateTrim(string filename, string qFileName, string trimFASTAFileName, string scrapFASTAFileName, string trimQualFileName, string scrapQualFileName, string trimNameFileName, string scrapNameFileName, string trimCountFileName, string scrapCountFileName) {
	try {
        string groupFile;
        Oligos oligos;
        if(oligoFile != ""){
            oligos.read(oligoFile);
            
            if (m->getControl_pressed()) { return 0; } //error in reading oligos
            
            if (oligos.hasPairedBarcodes()) { pairedOligos = true; }
            else { pairedOligos = false; }
            
            vector<string> groupNames = oligos.getGroupNames();
            if (groupNames.size() == 0) { allFiles = 0;   }
            else { createGroup = true; }

            if ((createGroup) && (countfile == "")){
                map<string, string> myvariables;
                myvariables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastaFile));
                groupFile = getOutputFileName("group",myvariables);
                outputNames.push_back(groupFile); outputTypes["group"].push_back(groupFile);
            }
        }
        
        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<trimData*> data;
        
        //fills lines and qlines
        setLines(filename, qFileName);
        
        auto synchronizedOutputFastaTrimFile = std::make_shared<SynchronizedOutputFile>(trimFASTAFileName);
        auto synchronizedOutputFastaScrapFile = std::make_shared<SynchronizedOutputFile>(scrapFASTAFileName);
        auto synchronizedOutputQTrimFile = std::make_shared<SynchronizedOutputFile>(trimQualFileName);
        auto synchronizedOutputQScrapFile = std::make_shared<SynchronizedOutputFile>(scrapQualFileName);
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            OutputWriter* threadFastaTrimWriter = new OutputWriter(synchronizedOutputFastaTrimFile);
            OutputWriter* threadFastaScrapWriter = new OutputWriter(synchronizedOutputFastaScrapFile);
            OutputWriter* threadQTrimWriter = NULL;
            OutputWriter* threadQScrapWriter = NULL;
            if (qFileName != "") {
                threadQTrimWriter = new OutputWriter(synchronizedOutputQTrimFile);
                threadQScrapWriter = new OutputWriter(synchronizedOutputQScrapFile);
            }
            //string fn, string qn, OutputWriter* tn, OutputWriter* sn, OutputWriter* tqn, OutputWriter* sqn, unsigned long long lstart, unsigned long long lend, unsigned long long qstart, unsigned long long qend, map<string, string> nm, map<string, int> ncoun
            trimData* dataBundle = new trimData(filename, qFileName, threadFastaTrimWriter, threadFastaScrapWriter, threadQTrimWriter, threadQScrapWriter, lines[i+1].start, lines[i+1].end, qLines[i+1].start, qLines[i+1].end, nameMap, nameCount);
            dataBundle->setOligosOptions(oligos, pdiffs, bdiffs, ldiffs, sdiffs, tdiffs, pairedOligos, createGroup, allFiles, keepforward, keepFirst, removeLast,
                                           qWindowStep, qWindowSize, qWindowAverage, qtrim, qThreshold, qAverage, qRollAverage, logtransform,
                                           minLength, maxAmbig, maxHomoP, maxLength, flip, reorient);
            data.push_back(dataBundle);
            
            workerThreads.push_back(new std::thread(driverTrim, dataBundle));
        }
        
        OutputWriter* threadFastaTrimWriter = new OutputWriter(synchronizedOutputFastaTrimFile);
        OutputWriter* threadFastaScrapWriter = new OutputWriter(synchronizedOutputFastaScrapFile);
        OutputWriter* threadQTrimWriter = NULL;
        OutputWriter* threadQScrapWriter = NULL;
        if (qFileName != "") {
            threadQTrimWriter = new OutputWriter(synchronizedOutputQTrimFile);
            threadQScrapWriter = new OutputWriter(synchronizedOutputQScrapFile);
        }
        //string fn, string qn, OutputWriter* tn, OutputWriter* sn, OutputWriter* tqn, OutputWriter* sqn, unsigned long long lstart, unsigned long long lend, unsigned long long qstart, unsigned long long qend, map<string, string> nm, map<string, int> ncoun
        trimData* dataBundle = new trimData(filename, qFileName, threadFastaTrimWriter, threadFastaScrapWriter, threadQTrimWriter, threadQScrapWriter, lines[0].start, lines[0].end, qLines[0].start, qLines[0].end, nameMap, nameCount);
        dataBundle->setOligosOptions(oligos, pdiffs, bdiffs, ldiffs, sdiffs, tdiffs, pairedOligos,
                                    createGroup, allFiles, keepforward, keepFirst, removeLast,
                                     qWindowStep, qWindowSize, qWindowAverage, qtrim, qThreshold, qAverage, qRollAverage, logtransform,
                                     minLength, maxAmbig, maxHomoP, maxLength, flip, reorient);

        driverTrim(dataBundle);
        long long num = dataBundle->count;
        
        set<string> badNames = dataBundle->badNames;
        groupCounts = dataBundle->groupCounts;
        groupMap = dataBundle->groupMap;

        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            num += data[i]->count;
            
            delete data[i]->trimFileName;
            delete data[i]->scrapFileName;
            if (qFileName != "") {
                delete data[i]->trimQFileName;
                delete data[i]->scrapQFileName;
            }
            badNames.insert(data[i]->badNames.begin(), data[i]->badNames.end());
            groupMap.insert(data[i]->groupMap.begin(), data[i]->groupMap.end());
            //merge counts
            for (map<string, int>::iterator it = data[i]->groupCounts.begin(); it != data[i]->groupCounts.end(); it++) {
                map<string, int>::iterator itMine = groupCounts.find(it->first);
                if (itMine != groupCounts.end()) { itMine->second += it->second; }
                else { groupCounts[it->first] = it->second; }
            }
            
            delete data[i];
            delete workerThreads[i];
        }
        synchronizedOutputFastaTrimFile->close(); synchronizedOutputQScrapFile->close();
        delete threadFastaTrimWriter;
        delete threadFastaScrapWriter;
        if (qFileName != "") {
            synchronizedOutputQTrimFile->close(); synchronizedOutputFastaTrimFile->close();
            delete threadQTrimWriter; delete threadQScrapWriter;
        }
        delete dataBundle;
        
        processNamesCountFiles(trimFASTAFileName, badNames, groupMap, trimNameFileName, scrapNameFileName, trimCountFileName, scrapCountFileName, groupFile);
        
        return num;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "createProcessesCreateTrim");
		exit(1);
	}
}
/**************************************************************************************************/
int TrimSeqsCommand::processNamesCountFiles(string trimFasta, set<string> badNames, map<string, string> groupMap, string trimNameFileName, string scrapNameFileName, string trimCountFileName, string scrapCountFileName, string groupFile) {
    try {
        
        if(nameFile != ""){
            ifstream inName;
            util.openInputFile(nameFile, inName);
            
            ofstream outTrimName;
            util.openOutputFile(trimNameFileName, outTrimName);
            
            ofstream outScrapName;
            util.openOutputFile(scrapNameFileName, outScrapName);
            
            string repName, redundNames;
            while (!inName.eof()) {
                if (m->getControl_pressed()) { break; }
                
                inName >> repName; util.gobble(inName);
                inName >> redundNames; util.gobble(inName);
                
                if (badNames.count(repName) == 0) { //good seq, put in trim
                    outTrimName << repName << '\t' << redundNames << endl;
                    int numRedunds = util.getNumNames(redundNames);
                    if (numRedunds != 1) { //update groupCounts and add to groupMap
                        map<string, string>::iterator itGroup = groupMap.find(repName);
                        if (itGroup != groupMap.end()) {
                            vector<string> seqNames; util.splitAtComma(redundNames, seqNames);
                            for (int i = 1; i < seqNames.size(); i++) { groupMap[seqNames[i]] = itGroup->second; }
                            groupCounts[itGroup->second] += (numRedunds - 1);
                        }
                    }
                }else { //bad seq -> scrap file
                    outScrapName << repName << '\t' << redundNames << endl;
                }
            }
            outScrapName.close();
            outTrimName.close();
        }
        
        //print group or count file.
        if (groupFile != "") {
            ofstream outGroup;
            util.openOutputFile(groupFile, outGroup);
            for (map<string, string>::iterator itGroup = groupMap.begin(); itGroup != groupMap.end(); itGroup++) {
                outGroup << itGroup->first << '\t' << itGroup->second << endl;
            }
            outGroup.close();
        }else {
            map<string, int> justTrimmedNames;
            map<string, int> justScrappedNames;
            for (map<string, int>::iterator itCount = nameCount.begin(); itCount != nameCount.end(); itCount++) {
                if (badNames.count(itCount->first) != 0) { justTrimmedNames[itCount->first] = itCount->second; }
                else { justScrappedNames[itCount->first] = itCount->second; }
            }
            
            CountTable newCt; //trimmed file
            for (map<string, int>::iterator itCount = groupCounts.begin(); itCount != groupCounts.end(); itCount++) { newCt.addGroup(itCount->first); }
            vector<int> tempCounts; tempCounts.resize(groupCounts.size(), 0);
            for (map<string, int>::iterator itNames = justTrimmedNames.begin(); itNames != justTrimmedNames.end(); itNames++) {
                newCt.push_back(itNames->first, tempCounts); //add it to the table with no abundance so we can set the groups abundance
                map<string, string>::iterator it2 = groupMap.find(itNames->first);
                if (it2 != groupMap.end()) { newCt.setAbund(itNames->first, it2->second, itNames->second); }
                else { m->mothurOut("[ERROR]: missing group info for " + itNames->first + ".\n");  m->setControl_pressed(true); }
            }
            newCt.printTable(trimCountFileName);
            
            CountTable newScrapCt; //scrap file
            newScrapCt.addGroup("scrap");
            vector<int> myCounts; tempCounts.resize(1, 0);
            for (map<string, int>::iterator itNames = justScrappedNames.begin(); itNames != justScrappedNames.end(); itNames++) {
                newScrapCt.push_back(itNames->first, myCounts); //add it to the table with no abundance so we can set the groups abundance
                newScrapCt.setAbund(itNames->first, "scrap", itNames->second);
            }
            newScrapCt.printTable(scrapCountFileName);
        }

        if(allFiles){
            //run split.groups command
            //use unique.seqs to create new name and fastafile
            string inputString = "fasta=" + trimFasta;
            if (nameFile != "")     { inputString += ", name=" + trimNameFileName;     }
            if (countfile != "")    { inputString += ", count=" + trimCountFileName;   }
            else {  inputString += ", group=" + groupFile; }
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Generating allfiles... Running command: split.groups(" + inputString + ")\n");
            current->setMothurCalling(true);
            
            Command* splitCommand = new SplitGroupCommand(inputString);
            splitCommand->execute();
            
            map<string, vector<string> > filenames = splitCommand->getOutputFiles();
            
            delete splitCommand;
            current->setMothurCalling(false);
            m->mothurOut("/******************************************/\n");
        }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "TrimSeqsCommand", "processNamesCountFiles");
        exit(1);
    }
}
/**************************************************************************************************/

int TrimSeqsCommand::setLines(string filename, string qfilename) {
	try {
        
        vector<unsigned long long> fastaFilePos;
		vector<unsigned long long> qfileFilePos;
		
		#if defined NON_WINDOWS
		//set file positions for fasta file
		fastaFilePos = util.divideFile(filename, processors);
		
		//get name of first sequence in each chunk
		map<string, int> firstSeqNames;
		for (int i = 0; i < (fastaFilePos.size()-1); i++) {
			ifstream in;
			util.openInputFile(filename, in);
			in.seekg(fastaFilePos[i]);
            
            //adjust start if null strings
            if (i == 0) {  util.zapGremlins(in); util.gobble(in);  }
		
			Sequence temp(in); 
			firstSeqNames[temp.getName()] = i;
		
			in.close();
		}
		
		if(qfilename != "")	{
            //seach for filePos of each first name in the qfile and save in qfileFilePos
            ifstream inQual;
            util.openInputFile(qfilename, inQual);
            
            string input;
            while(!inQual.eof()){	
                input = util.getline(inQual);
                
                if (input.length() != 0) {
                    if(input[0] == '>'){ //this is a sequence name line
                        istringstream nameStream(input);
                        
                        string sname = "";  nameStream >> sname;
                        sname = sname.substr(1);
                        
                        util.checkName(sname);
                        
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
            qfilename = util.getFullPathName(qfilename);
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
            if (m->getDebug()) { m->mothurOut("[DEBUG]: " + toString(i) +'\t' + toString(fastaFilePos[i]) + '\t' + toString(fastaFilePos[i+1]) + '\n'); }
			lines.push_back(linePair(fastaFilePos[i], fastaFilePos[(i+1)]));
			if (qfilename != "") {  qLines.push_back(linePair(qfileFilePos[i], qfileFilePos[(i+1)]));  }
		}	
		if(qfilename == "")	{	qLines = lines;	} //files with duds
		
		return processors;
		
		#else
        
            long long numFastaSeqs = 0;
            fastaFilePos = util.setFilePosFasta(filename, numFastaSeqs); 
            if (numFastaSeqs < processors) { processors = numFastaSeqs; }
        
            if (qfilename != "") { 
                long long numQualSeqs = 0;
                qfileFilePos = util.setFilePosFasta(qfilename, numQualSeqs); 
                
                if (numFastaSeqs != numQualSeqs) {
                    m->mothurOut("[ERROR]: You have " + toString(numFastaSeqs) + " sequences in your fasta file, but " + toString(numQualSeqs) + " sequences in your quality file."); m->mothurOutEndLine(); m->setControl_pressed(true); 
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

