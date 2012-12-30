//
//  makecontigscommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/15/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "makecontigscommand.h"

//**********************************************************************************************************************
vector<string> MakeContigsCommand::setParameters(){	
	try {
		CommandParameter pfastq("ffastq", "InputTypes", "", "", "FastaFastqFile", "FastaFastqFile", "fastqGroup","fasta-qfile",false,false,true); parameters.push_back(pfastq);
        CommandParameter prfastq("rfastq", "InputTypes", "", "", "none", "none", "fastqGroup","fasta-qfile",false,false,true); parameters.push_back(prfastq);
        CommandParameter pfasta("ffasta", "InputTypes", "", "", "FastaFastqFile", "FastaFastqFile", "fastaGroup","fasta",false,false,true); parameters.push_back(pfasta);
        CommandParameter prfasta("rfasta", "InputTypes", "", "", "none", "none", "none","fastaGroup",false,false,true); parameters.push_back(prfasta);
        CommandParameter pfqual("fqfile", "InputTypes", "", "", "none", "none", "qfileGroup","qfile",false,false,true); parameters.push_back(pfqual);
        CommandParameter prqual("rqfile", "InputTypes", "", "", "none", "none", "qfileGroup","qfile",false,false,true); parameters.push_back(prqual);
        CommandParameter pfile("file", "InputTypes", "", "", "FastaFastqFile", "FastaFastqFile", "none","fasta-qfile",false,false,true); parameters.push_back(pfile);
        CommandParameter poligos("oligos", "InputTypes", "", "", "none", "none", "none","group",false,false,true); parameters.push_back(poligos);
		CommandParameter ppdiffs("pdiffs", "Number", "", "0", "", "", "","",false,false,true); parameters.push_back(ppdiffs);
		CommandParameter pbdiffs("bdiffs", "Number", "", "0", "", "", "","",false,false,true); parameters.push_back(pbdiffs);
//        CommandParameter pldiffs("ldiffs", "Number", "", "0", "", "", "",false,false); parameters.push_back(pldiffs);
//		CommandParameter psdiffs("sdiffs", "Number", "", "0", "", "", "",false,false); parameters.push_back(psdiffs);
        CommandParameter ptdiffs("tdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(ptdiffs);

        CommandParameter palign("align", "Multiple", "needleman-gotoh", "needleman", "", "", "","",false,false); parameters.push_back(palign);
        CommandParameter pallfiles("allfiles", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pallfiles);
		CommandParameter pmatch("match", "Number", "", "1.0", "", "", "","",false,false); parameters.push_back(pmatch);
		CommandParameter pmismatch("mismatch", "Number", "", "-1.0", "", "", "","",false,false); parameters.push_back(pmismatch);
		CommandParameter pgapopen("gapopen", "Number", "", "-2.0", "", "", "","",false,false); parameters.push_back(pgapopen);
		CommandParameter pgapextend("gapextend", "Number", "", "-1.0", "", "", "","",false,false); parameters.push_back(pgapextend);
        CommandParameter pthreshold("threshold", "Number", "", "40", "", "", "","",false,false); parameters.push_back(pthreshold);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeContigsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string MakeContigsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The make.contigs command reads a file, forward fastq file and a reverse fastq file or forward fasta and reverse fasta files and outputs new fasta.  It will also provide new quality files if the fastq or file parameter is used.\n";
        helpString += "If an oligos file is provided barcodes and primers will be trimmed, and a group file will be created.\n";
		helpString += "The make.contigs command parameters are ffastq, rfastq, oligos, tdiffs, bdiffs, ldiffs, sdiffs, pdiffs, align, match, mismatch, gapopen, gapextend, allfiles and processors.\n";
		helpString += "The ffastq and rfastq, file, or ffasta and rfasta parameters are required.\n";
        helpString += "The file parameter is 2 column file containing the forward fastq files in the first column and their matching reverse fastq files in the second column.  Mothur will process each pair and create a combined fasta and qual file with all the sequences.\n";
        helpString += "The ffastq and rfastq parameters are used to provide a forward fastq and reverse fastq file to process.  If you provide one, you must provide the other.\n";
        helpString += "The ffasta and rfasta parameters are used to provide a forward fasta and reverse fasta file to process.  If you provide one, you must provide the other.\n";
        helpString += "The fqfile and rqfile parameters are used to provide a forward quality and reverse quality files to process with the ffasta and rfasta parameters.  If you provide one, you must provide the other.\n";
		helpString += "The align parameter allows you to specify the alignment method to use.  Your options are: gotoh and needleman. The default is needleman.\n";
        helpString += "The tdiffs parameter is used to specify the total number of differences allowed in the sequence. The default is pdiffs + bdiffs + sdiffs + ldiffs.\n";
		helpString += "The bdiffs parameter is used to specify the number of differences allowed in the barcode. The default is 0.\n";
		helpString += "The pdiffs parameter is used to specify the number of differences allowed in the primer. The default is 0.\n";
        helpString += "The ldiffs parameter is used to specify the number of differences allowed in the linker. The default is 0.\n";
		helpString += "The sdiffs parameter is used to specify the number of differences allowed in the spacer. The default is 0.\n";
		helpString += "The match parameter allows you to specify the bonus for having the same base. The default is 1.0.\n";
		helpString += "The mistmatch parameter allows you to specify the penalty for having different bases.  The default is -1.0.\n";
		helpString += "The gapopen parameter allows you to specify the penalty for opening a gap in an alignment. The default is -2.0.\n";
		helpString += "The gapextend parameter allows you to specify the penalty for extending a gap in an alignment.  The default is -1.0.\n";
        helpString += "The threshold parameter allows you to set a quality scores threshold. In the case where we are trying to decide whether to keep a base or remove it because the base is compared to a gap in the other fragment, if the base has a quality score below the threshold we eliminate it. Default=40.\n";
        helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n";
        helpString += "The allfiles parameter will create separate group and fasta file for each grouping. The default is F.\n";
        helpString += "The make.contigs command should be in the following format: \n";
		helpString += "make.contigs(ffastq=yourForwardFastqFile, rfastq=yourReverseFastqFile, align=yourAlignmentMethod) \n";
		helpString += "Note: No spaces between parameter labels (i.e. ffastq), '=' and parameters (i.e.yourForwardFastqFile).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeContigsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string MakeContigsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],[tag],contigs.fasta"; } 
        else if (type == "qfile") {  pattern = "[filename],[tag],contigs.qual"; } 
        else if (type == "group") {  pattern = "[filename],[tag],contigs.groups"; }
        else if (type == "mismatch") {  pattern = "[filename],[tag],contigs.mismatch"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
MakeContigsCommand::MakeContigsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["qfile"] = tempOutNames;
        outputTypes["group"] = tempOutNames;
        outputTypes["mismatch"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeContigsCommand", "MakeContigsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
MakeContigsCommand::MakeContigsCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
        
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters(); 
			
			ValidParameters validParameter("pairwise.seqs");
			map<string, string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["qfile"] = tempOutNames;
            outputTypes["mismatch"] = tempOutNames;
            outputTypes["group"] = tempOutNames;
			
            
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else { 
				string path;
                it = parameters.find("ffastq");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["ffastq"] = inputDir + it->second;		}
				}
                
                it = parameters.find("rfastq");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["rfastq"] = inputDir + it->second;		}
				}
                
                it = parameters.find("ffasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["ffasta"] = inputDir + it->second;		}
				}
                
                it = parameters.find("rfasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["rfasta"] = inputDir + it->second;		}
				}
                
                it = parameters.find("fqfile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fqfile"] = inputDir + it->second;		}
				}
                
                it = parameters.find("rqfile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["rqfile"] = inputDir + it->second;		}
				}
                
                it = parameters.find("file");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["file"] = inputDir + it->second;		}
				}
                
                it = parameters.find("oligos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["oligos"] = inputDir + it->second;		}
				}
            }
            
            ffastqfile = validParameter.validFile(parameters, "ffastq", true);
			if (ffastqfile == "not open") {  abort = true; }	
			else if (ffastqfile == "not found") { ffastqfile = ""; }
			
			rfastqfile = validParameter.validFile(parameters, "rfastq", true);
			if (rfastqfile == "not open") {  abort = true; }	
			else if (rfastqfile == "not found") { rfastqfile = "";  }
            
            ffastafile = validParameter.validFile(parameters, "ffasta", true);
			if (ffastafile == "not open") {  abort = true; }	
			else if (ffastafile == "not found") { ffastafile = ""; }
			
			rfastafile = validParameter.validFile(parameters, "rfasta", true);
			if (rfastafile == "not open") {  abort = true; }	
			else if (rfastafile == "not found") { rfastafile = "";  }
            
            fqualfile = validParameter.validFile(parameters, "fqfile", true);
			if (fqualfile == "not open") {  abort = true; }	
			else if (fqualfile == "not found") { fqualfile = ""; }
			
			rqualfile = validParameter.validFile(parameters, "rqfile", true);
			if (rqualfile == "not open") {  abort = true; }	
			else if (rqualfile == "not found") { rqualfile = "";  }
            
            file = validParameter.validFile(parameters, "file", true);
			if (file == "not open") {  abort = true; }	
			else if (file == "not found") { file = "";  }
            
            //provide at least
            if ((file == "") && (ffastafile == "") && (ffastqfile == "")) { abort = true; m->mothurOut("[ERROR]: The file, ffastq and rfastq or ffasta and rfasta parameters are required.\n"); }
            if ((file != "") && ((ffastafile != "") || (ffastqfile != ""))) { abort = true; m->mothurOut("[ERROR]: The file, ffastq and rfastq or ffasta and rfasta parameters are required.\n"); }
            if ((ffastqfile != "") && (rfastqfile == "")) {  abort = true; m->mothurOut("[ERROR]: If you provide use the ffastq, you must provide a rfastq file.\n"); }
            if ((ffastqfile == "") && (rfastqfile != "")) {  abort = true; m->mothurOut("[ERROR]: If you provide use the rfastq, you must provide a ffastq file.\n"); }
            if ((ffastafile != "") && (rfastafile == "")) {  abort = true; m->mothurOut("[ERROR]: If you provide use the ffasta, you must provide a rfasta file.\n"); }
            if ((ffastafile == "") && (rfastafile != "")) {  abort = true; m->mothurOut("[ERROR]: If you provide use the rfasta, you must provide a ffasta file.\n"); }
            if ((fqualfile != "") && (rqualfile == "")) {  abort = true; m->mothurOut("[ERROR]: If you provide use the fqfile, you must provide a rqfile file.\n"); }
            if ((fqualfile == "") && (rqualfile != "")) {  abort = true; m->mothurOut("[ERROR]: If you provide use the rqfile, you must provide a fqfile file.\n"); }
            if (((fqualfile != "") || (rqualfile != "")) && ((ffastafile == "") || (rfastafile == ""))) {
                abort = true; m->mothurOut("[ERROR]: If you provide use the rqfile or fqfile file, you must provide the ffasta and rfasta parameters.\n");
            }
            
            oligosfile = validParameter.validFile(parameters, "oligos", true);
			if (oligosfile == "not found")      {	oligosfile = "";	}
			else if(oligosfile == "not open")   {	abort = true;       } 
			else {	 m->setOligosFile(oligosfile);		}
            
            //if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
                 outputDir = ""; 
            }
			

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.validFile(parameters, "match", false);		if (temp == "not found"){	temp = "1.0";			}
			m->mothurConvert(temp, match);  
			
			temp = validParameter.validFile(parameters, "mismatch", false);		if (temp == "not found"){	temp = "-1.0";			}
			m->mothurConvert(temp, misMatch);  
            if (misMatch > 0) { m->mothurOut("[ERROR]: mismatch must be negative.\n"); abort=true; }
			
			temp = validParameter.validFile(parameters, "gapopen", false);		if (temp == "not found"){	temp = "-2.0";			}
			m->mothurConvert(temp, gapOpen);  
            if (gapOpen > 0) { m->mothurOut("[ERROR]: gapopen must be negative.\n"); abort=true; }
			
			temp = validParameter.validFile(parameters, "gapextend", false);	if (temp == "not found"){	temp = "-1.0";			}
			m->mothurConvert(temp, gapExtend); 
            if (gapExtend > 0) { m->mothurOut("[ERROR]: gapextend must be negative.\n"); abort=true; }
			
            temp = validParameter.validFile(parameters, "threshold", false);	if (temp == "not found"){	temp = "40";			}
			m->mothurConvert(temp, threshold); 
            if ((threshold < 0) || (threshold > 40)) { m->mothurOut("[ERROR]: threshold must be between 0 and 40.\n"); abort=true; }

			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
            
            temp = validParameter.validFile(parameters, "bdiffs", false);		if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, bdiffs);
			
			temp = validParameter.validFile(parameters, "pdiffs", false);		if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, pdiffs);
            
  //          temp = validParameter.validFile(parameters, "ldiffs", false);		if (temp == "not found") { temp = "0"; }
//			m->mothurConvert(temp, ldiffs);
            ldiffs = 0;
            
 //           temp = validParameter.validFile(parameters, "sdiffs", false);		if (temp == "not found") { temp = "0"; }
 //           m->mothurConvert(temp, sdiffs);
            sdiffs = 0;
			
			temp = validParameter.validFile(parameters, "tdiffs", false);		if (temp == "not found") { int tempTotal = pdiffs + bdiffs;  temp = toString(tempTotal); }
			m->mothurConvert(temp, tdiffs);
			
			if(tdiffs == 0){	tdiffs = bdiffs + pdiffs;	}  //+ ldiffs + sdiffs;

            temp = validParameter.validFile(parameters, "allfiles", false);		if (temp == "not found") { temp = "F"; }
			allFiles = m->isTrue(temp);
			
			align = validParameter.validFile(parameters, "align", false);		if (align == "not found"){	align = "needleman";	}
			if ((align != "needleman") && (align != "gotoh")) { m->mothurOut(align + " is not a valid alignment method. Options are needleman or gotoh. I will use needleman."); m->mothurOutEndLine(); align = "needleman"; }
        }
		
	}
	catch(exception& e) {
		m->errorOut(e, "MakeContigsCommand", "MakeContigsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int MakeContigsCommand::execute(){
	try {
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        //read ffastq and rfastq files creating fasta and qual files.
        //this function will create a forward and reverse, fasta and qual files for each processor.
        //files has an entry for each processor. files[i][0] = forwardFasta, files[i][1] = forwardQual, files[i][2] = reverseFasta, files[i][3] = reverseQual.  filesToProcess is for each filepair in the file parameter file.  for ffastq and rfastq this will be size 1.
        unsigned long int numReads = 0;
        int start = time(NULL);
        longestBase = 1000;
        m->mothurOut("Reading fastq data...\n"); 
        vector < vector< vector<string> > > filesToProcess = preProcessData(numReads);
        m->mothurOut("Done.\n");
       
        if (m->control_pressed) { return 0; }
        
        map<string, string> cvars;
        string compOutputDir = outputDir;
        if (outputDir == "") { compOutputDir = m->hasPath(file); }
        cvars["[filename]"] = compOutputDir + m->getRootName(m->getSimpleName(file));
        cvars["[tag]"] = "";
        string compositeGroupFile = getOutputFileName("group",cvars);
        cvars["[tag]"] = "trim";
        string compositeFastaFile = getOutputFileName("fasta",cvars);
        string compositeQualFile = getOutputFileName("qfile",cvars);
        cvars["[tag]"] = "scrap";
        string compositeScrapFastaFile = getOutputFileName("fasta",cvars);
        string compositeScrapQualFile = getOutputFileName("qfile",cvars);
        cvars["[tag]"] = "";
        string compositeMisMatchFile = getOutputFileName("mismatch",cvars);
        
        if (filesToProcess.size() > 1) { //clear files for append below
            ofstream outCTFasta, outCTQual, outCSFasta, outCSQual, outCMisMatch;
            m->openOutputFile(compositeFastaFile, outCTFasta); outCTFasta.close();
            m->openOutputFile(compositeScrapFastaFile, outCSFasta); outCSFasta.close();
            m->openOutputFile(compositeMisMatchFile, outCMisMatch); outCMisMatch.close();
            m->openOutputFile(compositeQualFile, outCTQual); outCTQual.close();
            m->openOutputFile(compositeScrapQualFile, outCSQual); outCSQual.close();
            outputNames.push_back(compositeFastaFile); outputTypes["fasta"].push_back(compositeFastaFile);
            outputNames.push_back(compositeQualFile); outputTypes["qfile"].push_back(compositeQualFile);
            outputNames.push_back(compositeMisMatchFile); outputTypes["mismatch"].push_back(compositeMisMatchFile);
            outputNames.push_back(compositeScrapFastaFile); outputTypes["fasta"].push_back(compositeScrapFastaFile);
            outputNames.push_back(compositeScrapQualFile); outputTypes["qfile"].push_back(compositeScrapQualFile);
        }
        
        for (int l = 0; l < filesToProcess.size(); l++) {
            
            m->mothurOut("\n>>>>>\tProcessing " + filesToProcess[l][0][0] + " (file " + toString(l+1) + " of " + toString(filesToProcess.size()) + ")\t<<<<<\n");
            
            vector<vector<string> > fastaFileNames;
            vector<vector<string> > qualFileNames;
            createGroup = false;
            string outputGroupFileName;
            map<string, string> variables; 
            string thisOutputDir = outputDir;
            if (outputDir == "") {  thisOutputDir = m->hasPath(filesToProcess[l][0][0]); }
            variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(filesToProcess[l][0][0]));
            variables["[tag]"] = "";
            if(oligosfile != ""){
                createGroup = getOligos(fastaFileNames, qualFileNames, variables["[filename]"]);
                if (createGroup) { 
                    outputGroupFileName = getOutputFileName("group",variables);
                    outputNames.push_back(outputGroupFileName); outputTypes["group"].push_back(outputGroupFileName);
                }
            }
            
            variables["[tag]"] = "trim";
            string outFastaFile = getOutputFileName("fasta",variables);
            string outQualFile = getOutputFileName("qfile",variables);
            variables["[tag]"] = "scrap";
            string outScrapFastaFile = getOutputFileName("fasta",variables);
            string outScrapQualFile = getOutputFileName("qfile",variables);
            variables["[tag]"] = "";
            string outMisMatchFile = getOutputFileName("mismatch",variables);
            outputNames.push_back(outFastaFile); outputTypes["fasta"].push_back(outFastaFile);
            outputNames.push_back(outScrapFastaFile); outputTypes["fasta"].push_back(outScrapFastaFile);
            if (filesToProcess[l][0][1] != "") {
                outputNames.push_back(outQualFile); outputTypes["qfile"].push_back(outQualFile);
                outputNames.push_back(outScrapQualFile); outputTypes["qfile"].push_back(outScrapQualFile);
            }
            outputNames.push_back(outMisMatchFile); outputTypes["mismatch"].push_back(outMisMatchFile);
            
            m->mothurOut("Making contigs...\n"); 
            createProcesses(filesToProcess[l], outFastaFile, outQualFile, outScrapFastaFile, outScrapQualFile, outMisMatchFile, fastaFileNames, qualFileNames);
            m->mothurOut("Done.\n");
            
            //remove temp fasta and qual files
            for (int i = 0; i < processors; i++) { for(int j = 0; j < filesToProcess[l][i].size(); j++) { m->mothurRemove(filesToProcess[l][i][j]); }  }
            
            if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); }  return 0; }
            
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
                                    
                                    if (filesToProcess[l][0][1] != "") {
                                        m->mothurRemove(qualFileNames[i][j]);
                                        namesToRemove.insert(qualFileNames[i][j]);
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
                    string thisGroupName = thisOutputDir + m->getRootName(m->getSimpleName(it->first));
                    thisGroupName += getOutputFileName("group",variables); outputNames.push_back(thisGroupName); outputTypes["group"].push_back(thisGroupName); 
                    m->openOutputFile(thisGroupName, out);
                    
                    while (!in.eof()){
                        if (m->control_pressed) { break; }
                        
                        Sequence currSeq(in); m->gobble(in);
                        out << currSeq.getName() << '\t' << it->second << endl;  
                    }
                    in.close();
                    out.close();
                }
            }
            
            if (createGroup) {
                ofstream outGroup;
                m->openOutputFile(outputGroupFileName, outGroup);
                for (map<string, string>::iterator itGroup = groupMap.begin(); itGroup != groupMap.end(); itGroup++) {
                    outGroup << itGroup->first << '\t' << itGroup->second << endl;
                }
                outGroup.close();
            }
            
            if (filesToProcess.size() > 1) { //merge into large combo files
                if (createGroup) {  
                    if (l == 0) { 
                        ofstream outCGroup;
                        m->openOutputFile(compositeGroupFile, outCGroup); outCGroup.close();
                        outputNames.push_back(compositeGroupFile); outputTypes["group"].push_back(compositeGroupFile);
                    }
                    m->appendFiles(outputGroupFileName, compositeGroupFile);  
                }
                m->appendFiles(outMisMatchFile, compositeMisMatchFile);
                m->appendFiles(outFastaFile, compositeFastaFile);
                m->appendFiles(outQualFile, compositeQualFile);
                m->appendFiles(outScrapFastaFile, compositeScrapFastaFile);
                m->appendFiles(outScrapQualFile, compositeScrapQualFile);
            }
        }
        m->mothurOut("It took " + toString(time(NULL) - start) + " secs to process " + toString(numReads) + " sequences.\n");
        
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
        
        string currentFasta = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentFasta = (itTypes->second)[0]; m->setFastaFile(currentFasta); }
		}
        
        string currentQual = "";
		itTypes = outputTypes.find("qfile");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentQual = (itTypes->second)[0]; m->setQualFile(currentQual); }
		}
        
        string currentGroup = "";
		itTypes = outputTypes.find("group");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentGroup = (itTypes->second)[0]; m->setGroupFile(currentGroup); }
		}
		
        //output files created by command
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "MakeContigsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
vector< vector< vector<string> > > MakeContigsCommand::preProcessData(unsigned long int& numReads) {
	try {
        vector< vector< vector<string> > > filesToProcess;
        
        if (ffastqfile != "") { //reading one file
            vector< vector<string> > files = readFastqFiles(numReads, ffastqfile, rfastqfile); 
            //adjust for really large processors or really small files
            if (numReads == 0) {  m->mothurOut("[ERROR]: no good reads.\n"); m->control_pressed = true; }
            if (numReads < processors) { 
                for (int i = numReads; i < processors; i++) { for(int j = 0; j < files[i].size(); j++) { m->mothurRemove(files[i][j]); } files[i].clear(); }
                files.resize(numReads);
                processors = numReads; 
            }
            filesToProcess.push_back(files);
        }else if (file != "") { //reading multiple files
            //return only valid pairs
            vector< vector<string> > filePairsToProcess = readFileNames(file);
            
            if (m->control_pressed) { return filesToProcess; }
            
            if (filePairsToProcess.size() != 0) {
                for (int i = 0; i < filePairsToProcess.size(); i++) {
                    
                    if (m->control_pressed) { for (int l = 0; l < filesToProcess.size(); l++) { for (int k = 0; k < filesToProcess[l].size(); k++) { for(int j = 0; j < filesToProcess[l][k].size(); j++) { m->mothurRemove(filesToProcess[l][k][j]); } filesToProcess[l][k].clear(); } return filesToProcess; } }
                    
                    unsigned long int thisFilesReads;
                    vector< vector<string> > files = readFastqFiles(thisFilesReads, filePairsToProcess[i][0], filePairsToProcess[i][1]); 
                    
                    //adjust for really large processors or really small files
                    if (thisFilesReads < processors) { 
                        m->mothurOut("[ERROR]: " + filePairsToProcess[i][0] + " has less than " + toString(processors) + " good reads, skipping\n"); 
                        for (int k = 0; k < files.size(); k++) { for(int j = 0; j < files[k].size(); j++) { m->mothurRemove(files[k][j]); } files[k].clear(); }
                    }else {
                        filesToProcess.push_back(files);
                        numReads += thisFilesReads;
                    }
                }
                //all files are bad
                if (numReads == 0) {  m->control_pressed = true; }
            }
        }else if (ffastafile != "") {
            vector< vector<string> > files = readFastaFiles(numReads, ffastafile, rfastafile);
            //adjust for really large processors or really small files
            if (numReads == 0) {  m->mothurOut("[ERROR]: no good reads.\n"); m->control_pressed = true; }
            if (numReads < processors) { 
                for (int i = numReads; i < processors; i++) { for(int j = 0; j < files[i].size(); j++) { m->mothurRemove(files[i][j]); } files[i].clear(); }
                files.resize(numReads);
                processors = numReads; 
            }
            filesToProcess.push_back(files);
        }else { m->control_pressed = true; } //should not get here
        
        return filesToProcess;
    }
	catch(exception& e) {
		m->errorOut(e, "MakeContigsCommand", "preProcessData");
		exit(1);
	}
}
//**********************************************************************************************************************
int MakeContigsCommand::createProcesses(vector< vector<string> > files, string outputFasta, string outputQual, string outputScrapFasta, string outputScrapQual, string outputMisMatches, vector<vector<string> > fastaFileNames, vector<vector<string> > qualFileNames) {
	try {
		int num = 0;
		vector<int> processIDS;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		int process = 0;
		
		//loop through and create all the processes you want
		while (process != processors-1) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
                vector<vector<string> > tempFASTAFileNames = fastaFileNames;
				vector<vector<string> > tempPrimerQualFileNames = qualFileNames;
                
				if(allFiles){
					ofstream temp;
                    
					for(int i=0;i<tempFASTAFileNames.size();i++){
						for(int j=0;j<tempFASTAFileNames[i].size();j++){
							if (tempFASTAFileNames[i][j] != "") {
								tempFASTAFileNames[i][j] += toString(getpid()) + ".temp";
								m->openOutputFile(tempFASTAFileNames[i][j], temp);			temp.close();
                                
                                if (files[processors-1][1] != "") {
                                    tempPrimerQualFileNames[i][j] += toString(getpid()) + ".temp";
                                    m->openOutputFile(tempPrimerQualFileNames[i][j], temp);		temp.close();
                                }
							}
						}
					}
				}

				num = driver(files[process], 
                             outputFasta + toString(getpid()) + ".temp", 
                             outputQual + toString(getpid()) + ".temp", 
                             outputScrapFasta + toString(getpid()) + ".temp", 
                             outputScrapQual + toString(getpid()) + ".temp",
                             outputMisMatches + toString(getpid()) + ".temp",
                             tempFASTAFileNames,
                             tempPrimerQualFileNames);
				
				//pass groupCounts to parent
                ofstream out;
                string tempFile = toString(getpid()) + ".num.temp";
                m->openOutputFile(tempFile, out);
                out << num << endl;
				if(createGroup){
					out << groupCounts.size() << endl;
					
					for (map<string, int>::iterator it = groupCounts.begin(); it != groupCounts.end(); it++) {
						out << it->first << '\t' << it->second << endl;
					}
                    
                    out << groupMap.size() << endl;
                    for (map<string, string>::iterator it = groupMap.begin(); it != groupMap.end(); it++) {
						out << it->first << '\t' << it->second << endl;
					}
				}
                out.close();
				
				exit(0);
			}else { 
				m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
				for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
				exit(0);
			}
		}
		
        ofstream temp;
		m->openOutputFile(outputFasta, temp);		temp.close();
        m->openOutputFile(outputScrapFasta, temp);		temp.close();
        if (files[processors-1][1] != "") {
            m->openOutputFile(outputScrapQual, temp);		temp.close();
            m->openOutputFile(outputQual, temp);	temp.close();
        }
        
		//do my part
		num = driver(files[processors-1], outputFasta, outputQual, outputScrapFasta, outputScrapQual, outputMisMatches, fastaFileNames, qualFileNames);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
        
		for (int i = 0; i < processIDS.size(); i++) {
            ifstream in;
            string tempFile = toString(processIDS[i]) + ".num.temp";
            m->openInputFile(tempFile, in);
            int tempNum;
            in >> tempNum; num += tempNum; m->gobble(in);
            
			if(createGroup){
				string group;
				in >> tempNum; m->gobble(in);
				
				if (tempNum != 0) {
					for (int j = 0; j < tempNum; j++) { 
                        int groupNum;
						in >> group >> groupNum; m->gobble(in);
                        
						map<string, int>::iterator it = groupCounts.find(group);
						if (it == groupCounts.end()) {	groupCounts[group] = groupNum; }
						else { groupCounts[it->first] += groupNum; }
					}
				}
                in >> tempNum; m->gobble(in);
                if (tempNum != 0) {
					for (int j = 0; j < tempNum; j++) { 
                        string group, seqName;
						in >> seqName >> group; m->gobble(in);
                        
						map<string, string>::iterator it = groupMap.find(seqName);
						if (it == groupMap.end()) {	groupMap[seqName] = group; }
						else { m->mothurOut("[ERROR]: " + seqName + " is in your fasta file more than once. Sequence names must be unique. please correct.\n");  }
					}
				}
			}
            in.close(); m->mothurRemove(tempFile);
        }
    #else
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the contigsData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<contigsData*> pDataArray; 
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
		
		//Create processor worker threads.
		for( int h=0; h<processors-1; h++ ){
			string extension = "";
			if (h != 0) { extension = toString(h) + ".temp"; processIDS.push_back(h); }
            vector<vector<string> > tempFASTAFileNames = fastaFileNames;
            vector<vector<string> > tempPrimerQualFileNames = qualFileNames;
            
            if(allFiles){
                ofstream temp;
                
                for(int i=0;i<tempFASTAFileNames.size();i++){
                    for(int j=0;j<tempFASTAFileNames[i].size();j++){
                        if (tempFASTAFileNames[i][j] != "") {
                            tempFASTAFileNames[i][j] += extension;
                            m->openOutputFile(tempFASTAFileNames[i][j], temp);			temp.close();
                            
                            if (files[processors-1][1] != "") {
                                tempPrimerQualFileNames[i][j] += extension;
                                m->openOutputFile(tempPrimerQualFileNames[i][j], temp);		temp.close();
                            }
                        }
                    }
                }
            }

			          
			contigsData* tempcontig = new contigsData(files[h], (outputFasta + extension), (outputQual + extension), (outputScrapFasta + extension), (outputScrapQual + extension),(outputMisMatches + extension), align, m, match, misMatch, gapOpen, gapExtend, threshold, barcodes, primers, tempFASTAFileNames, tempPrimerQualFileNames, barcodeNameVector, primerNameVector, pdiffs, bdiffs, tdiffs, createGroup, allFiles, h);
			pDataArray.push_back(tempcontig);
            
			hThreadArray[h] = CreateThread(NULL, 0, MyContigsThreadFunction, pDataArray[h], 0, &dwThreadIdArray[h]);   
		}
        
        vector<vector<string> > tempFASTAFileNames = fastaFileNames;
        vector<vector<string> > tempPrimerQualFileNames = qualFileNames;

        if(allFiles){
            ofstream temp;
            string extension = toString(processors-1) + ".temp";
            
            for(int i=0;i<tempFASTAFileNames.size();i++){
                for(int j=0;j<tempFASTAFileNames[i].size();j++){
                    if (tempFASTAFileNames[i][j] != "") {
                        tempFASTAFileNames[i][j] += extension;
                        m->openOutputFile(tempFASTAFileNames[i][j], temp);			temp.close();
                        
                        if (files[processors-1][1] != "") {
                            tempPrimerQualFileNames[i][j] += extension;
                            m->openOutputFile(tempPrimerQualFileNames[i][j], temp);		temp.close();
                        }
                    }
                }
            }
        }

		//parent do my part
		ofstream temp;
		m->openOutputFile(outputFasta, temp);		temp.close();
        m->openOutputFile(outputScrapFasta, temp);		temp.close();
        if (files[processors-1][1] != "") {
            m->openOutputFile(outputScrapQual, temp);		temp.close();
            m->openOutputFile(outputQual, temp);	temp.close();
		}
        
        //do my part
        processIDS.push_back(processors-1);
		num = driver(files[processors-1], (outputFasta+ toString(processors-1) + ".temp"), (outputQual+ toString(processors-1) + ".temp"), (outputScrapFasta+ toString(processors-1) + ".temp"), (outputScrapQual+ toString(processors-1) + ".temp"), (outputMisMatches+ toString(processors-1) + ".temp"), tempFASTAFileNames, tempPrimerQualFileNames);	
        
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
			num += pDataArray[i]->count;
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
        
        for (int i = 0; i < processIDS.size(); i++) {
			m->appendFiles((outputFasta + toString(processIDS[i]) + ".temp"), outputFasta);
			m->mothurRemove((outputFasta + toString(processIDS[i]) + ".temp"));
			
			m->appendFiles((outputScrapFasta + toString(processIDS[i]) + ".temp"), outputScrapFasta);
			m->mothurRemove((outputScrapFasta + toString(processIDS[i]) + ".temp"));
			
            if (files[processors-1][1] != "") {
                m->appendFiles((outputScrapQual + toString(processIDS[i]) + ".temp"), outputScrapQual);
                m->mothurRemove((outputScrapQual + toString(processIDS[i]) + ".temp"));
                
                m->appendFiles((outputQual + toString(processIDS[i]) + ".temp"), outputQual);
                m->mothurRemove((outputQual + toString(processIDS[i]) + ".temp"));
            }
            
            m->appendFiles((outputMisMatches + toString(processIDS[i]) + ".temp"), outputMisMatches);
			m->mothurRemove((outputMisMatches + toString(processIDS[i]) + ".temp"));
            
            if(allFiles){
				for(int j=0;j<fastaFileNames.size();j++){
					for(int k=0;k<fastaFileNames[j].size();k++){
						if (fastaFileNames[j][k] != "") {
							m->appendFiles((fastaFileNames[j][k] + toString(processIDS[i]) + ".temp"), fastaFileNames[j][k]);
							m->mothurRemove((fastaFileNames[j][k] + toString(processIDS[i]) + ".temp"));
							
                            if (files[processors-1][1] != "") {
                                m->appendFiles((qualFileNames[j][k] + toString(processIDS[i]) + ".temp"), qualFileNames[j][k]);
                                m->mothurRemove((qualFileNames[j][k] + toString(processIDS[i]) + ".temp"));
                            }
						}
					}
				}
			}
		}
		
		return num;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeContigsCommand", "createProcesses");
		exit(1);
	}
}
//**********************************************************************************************************************
int MakeContigsCommand::driver(vector<string> files, string outputFasta, string outputQual, string outputScrapFasta, string outputScrapQual, string outputMisMatches, vector<vector<string> > fastaFileNames, vector<vector<string> > qualFileNames){
    try {
        
        Alignment* alignment;
        if(align == "gotoh")			{	alignment = new GotohOverlap(gapOpen, gapExtend, match, misMatch, longestBase);			}
		else if(align == "needleman")	{	alignment = new NeedlemanOverlap(gapOpen, match, misMatch, longestBase);				}
        
        int num = 0;
        string thisffastafile = files[0];
        string thisfqualfile = files[1];
        string thisrfastafile = files[2];
        string thisrqualfile = files[3];
        
        if (m->debug) {  m->mothurOut("[DEBUG]: ffasta = " + thisffastafile + ".\n[DEBUG]: fqual = " + thisfqualfile + ".\n[DEBUG]: rfasta = " + thisrfastafile + ".\n[DEBUG]: rqual = " + thisrqualfile + ".\n"); }
        
        ifstream inFFasta, inRFasta, inFQual, inRQual;
        ofstream outFasta, outQual, outMisMatch, outScrapFasta, outScrapQual;
        m->openInputFile(thisffastafile, inFFasta);
        m->openInputFile(thisrfastafile, inRFasta);
        if (thisfqualfile != "") {
            m->openInputFile(thisfqualfile, inFQual);
            m->openInputFile(thisrqualfile, inRQual);
            m->openOutputFile(outputScrapQual, outScrapQual);
            m->openOutputFile(outputQual, outQual);
        }
        m->openOutputFile(outputFasta, outFasta);
        m->openOutputFile(outputScrapFasta, outScrapFasta);
        m->openOutputFile(outputMisMatches, outMisMatch);
        outMisMatch << "Name\tLength\tMisMatches\n";
        
        TrimOligos trimOligos(pdiffs, bdiffs, 0, 0, primers, barcodes);
        
        while ((!inFFasta.eof()) && (!inRFasta.eof())) {
            
            if (m->control_pressed) { break; }
            
            int success = 1;
            string trashCode = "";
            int currentSeqsDiffs = 0;

            //read seqs and quality info
            Sequence fSeq(inFFasta); m->gobble(inFFasta);
            Sequence rSeq(inRFasta); m->gobble(inRFasta);
            QualityScores* fQual = NULL; QualityScores* rQual = NULL;
            if (thisfqualfile != "") {
                fQual = new QualityScores(inFQual); m->gobble(inFQual);
                rQual = new QualityScores(inRQual); m->gobble(inRQual);
            }
            
            int barcodeIndex = 0;
            int primerIndex = 0;
            
            if(barcodes.size() != 0){
                if (thisfqualfile != "") {
                    success = trimOligos.stripBarcode(fSeq, rSeq, *fQual, *rQual, barcodeIndex);
                }else {
                    success = trimOligos.stripBarcode(fSeq, rSeq, barcodeIndex);
                }
                if(success > bdiffs)		{	trashCode += 'b';	}
                else{ currentSeqsDiffs += success;  }
            }
            
            if(primers.size() != 0){
                if (thisfqualfile != "") {
                    success = trimOligos.stripForward(fSeq, rSeq, *fQual, *rQual, primerIndex);
                }else {
                    success = trimOligos.stripForward(fSeq, rSeq, primerIndex);
                }
                if(success > pdiffs)		{	trashCode += 'f';	}
                else{ currentSeqsDiffs += success;  }
            }
            
            if (currentSeqsDiffs > tdiffs)	{	trashCode += 't';   }
            
            //flip the reverse reads
            rSeq.reverseComplement();
            if (thisfqualfile != "") { rQual->flipQScores(); }

            //pairwise align
            alignment->align(fSeq.getUnaligned(), rSeq.getUnaligned());
            map<int, int> ABaseMap = alignment->getSeqAAlnBaseMap();
            map<int, int> BBaseMap = alignment->getSeqBAlnBaseMap();
            fSeq.setAligned(alignment->getSeqAAln());
            rSeq.setAligned(alignment->getSeqBAln());
            int length = fSeq.getAligned().length();
            
            //traverse alignments merging into one contiguous seq
            string contig = "";
            vector<int> contigScores; 
            int numMismatches = 0;
            string seq1 = fSeq.getAligned();
            string seq2 = rSeq.getAligned();
            vector<int> scores1, scores2; 
            if (thisfqualfile != "") {
                scores1 = fQual->getQualityScores();
                scores2 = rQual->getQualityScores();
                delete fQual; delete rQual;
            }
            
            // if (num < 5) {  cout << fSeq.getStartPos() << '\t' << fSeq.getEndPos() << '\t' << rSeq.getStartPos() << '\t' << rSeq.getEndPos() << endl; }
            int overlapStart = fSeq.getStartPos();
            int seq2Start = rSeq.getStartPos();
            //bigger of the 2 starting positions is the location of the overlapping start
            if (overlapStart < seq2Start) { //seq2 starts later so take from 0 to seq2Start from seq1
                overlapStart = seq2Start; 
                for (int i = 0; i < overlapStart; i++) {
                    contig += seq1[i];
                    if (thisfqualfile != "") { contigScores.push_back(scores1[ABaseMap[i]]); }
                }
            }else { //seq1 starts later so take from 0 to overlapStart from seq2
                for (int i = 0; i < overlapStart; i++) {
                    contig += seq2[i];
                    if (thisfqualfile != "") { contigScores.push_back(scores2[BBaseMap[i]]); }
                }
            }
            
            int seq1End = fSeq.getEndPos();
            int seq2End = rSeq.getEndPos();
            int overlapEnd = seq1End;
            if (seq2End < overlapEnd) { overlapEnd = seq2End; }  //smallest end position is where overlapping ends
            
            for (int i = overlapStart; i < overlapEnd; i++) {
                if (seq1[i] == seq2[i]) { //match, add base and choose highest score
                    contig += seq1[i];
                    if (thisfqualfile != "") { 
                        contigScores.push_back(scores1[ABaseMap[i]]); 
                        if (scores1[ABaseMap[i]] < scores2[BBaseMap[i]]) { contigScores[contigScores.size()-1] = scores2[BBaseMap[i]]; }
                    }
                }else if (((seq1[i] == '.') || (seq1[i] == '-')) && ((seq2[i] != '-') && (seq2[i] != '.'))) { //seq1 is a gap and seq2 is a base, choose seq2, unless quality score for base is below threshold. In that case eliminate base
                    if (thisfqualfile != "") {
                        if (scores2[BBaseMap[i]] < threshold) { } //
                        else {
                            contig += seq2[i];
                            contigScores.push_back(scores2[BBaseMap[i]]);
                        }
                    }else { contig += seq2[i]; } //with no quality info, then we keep it?
                }else if (((seq2[i] == '.') || (seq2[i] == '-')) && ((seq1[i] != '-') && (seq1[i] != '.'))) { //seq2 is a gap and seq1 is a base, choose seq1, unless quality score for base is below threshold. In that case eliminate base
                    if (thisfqualfile != "") {
                        if (scores1[ABaseMap[i]] < threshold) { } //
                        else {
                            contig += seq1[i];
                            contigScores.push_back(scores1[ABaseMap[i]]);
                        }
                    }else { contig += seq1[i]; } //with no quality info, then we keep it?
                }else if (((seq1[i] != '-') && (seq1[i] != '.')) && ((seq2[i] != '-') && (seq2[i] != '.'))) { //both bases choose one with better quality
                    if (thisfqualfile != "") {
                        char c = seq1[i];
                        contigScores.push_back(scores1[ABaseMap[i]]);
                        if (scores1[ABaseMap[i]] < scores2[BBaseMap[i]]) { contigScores[contigScores.size()-1] = scores2[BBaseMap[i]]; c = seq2[i]; }
                        contig += c;
                        numMismatches++;
                    }else { numMismatches++; } //cant decide, so eliminate and mark as mismatch
                }else { //should never get here
                    m->mothurOut("[ERROR]: case I didn't think of seq1 = " + toString(seq1[i]) + " and seq2 = " + toString(seq2[i]) + "\n");
                }
            }
            
            if (seq1End < seq2End) { //seq1 ends before seq2 so take from overlap to length from seq2
                for (int i = overlapEnd; i < length; i++) {
                    contig += seq2[i];
                    if (thisfqualfile != "") { contigScores.push_back(scores2[BBaseMap[i]]); }
                }
            }else { //seq2 ends before seq1 so take from overlap to length from seq1
                for (int i = overlapEnd; i < length; i++) {
                    contig += seq1[i];
                    if (thisfqualfile != "") { contigScores.push_back(scores1[ABaseMap[i]]); }
                }
                
            }
            
            if(trashCode.length() == 0){
                bool ignore = false;
                
                if (m->debug) { m->mothurOut(fSeq.getName()); }
                
                if (createGroup) {
                    if(barcodes.size() != 0){
                        string thisGroup = barcodeNameVector[barcodeIndex];
                        if (primers.size() != 0) { 
                            if (primerNameVector[primerIndex] != "") { 
                                if(thisGroup != "") {
                                    thisGroup += "." + primerNameVector[primerIndex]; 
                                }else {
                                    thisGroup = primerNameVector[primerIndex]; 
                                }
                            } 
                        }
                        
                        if (m->debug) { m->mothurOut(", group= " + thisGroup + "\n"); }
                        
                        int pos = thisGroup.find("ignore");
                        if (pos == string::npos) {
                            groupMap[fSeq.getName()] = thisGroup; 
                        
                            map<string, int>::iterator it = groupCounts.find(thisGroup);
                            if (it == groupCounts.end()) {	groupCounts[thisGroup] = 1; }
                            else { groupCounts[it->first] ++; }
                        }else { ignore = true; }
                        
                    }
                }
                if (m->debug) { m->mothurOut("\n"); }
                
                if(allFiles && !ignore){
                    ofstream output;
                    m->openOutputFileAppend(fastaFileNames[barcodeIndex][primerIndex], output);
                    output << ">" << fSeq.getName() << endl << contig << endl;
                    output.close();
                    
                    if (thisfqualfile != "") {
                        m->openOutputFileAppend(qualFileNames[barcodeIndex][primerIndex], output);
                        output << ">" << fSeq.getName() << endl;
                        for (int i = 0; i < contigScores.size(); i++) { output << contigScores[i] << ' '; }
                        output << endl;
                        output.close();	
                    }
                }
                
                //output
                outFasta << ">" << fSeq.getName() << endl << contig << endl;
                if (thisfqualfile != "") {
                    outQual << ">" << fSeq.getName() << endl;
                    for (int i = 0; i < contigScores.size(); i++) { outQual << contigScores[i] << ' '; }
                    outQual << endl;
                }
                outMisMatch << fSeq.getName() << '\t' << contig.length() << '\t' << numMismatches << endl;
            }else {
                //output
                outScrapFasta << ">" << fSeq.getName() << " | " << trashCode << endl << contig << endl;
                if (thisfqualfile != "") {
                    outScrapQual << ">" << fSeq.getName() << " | " << trashCode << endl;
                    for (int i = 0; i < contigScores.size(); i++) { outScrapQual << contigScores[i] << ' '; }
                    outScrapQual << endl;
                }
            }
            num++;
            
			//report progress
			if((num) % 1000 == 0){	m->mothurOut(toString(num)); m->mothurOutEndLine();		}
		}
        
		//report progress
		if((num) % 1000 != 0){	m->mothurOut(toString(num)); m->mothurOutEndLine();		}
        
        inFFasta.close();
        inRFasta.close();
        outFasta.close();
        outScrapFasta.close();
        outMisMatch.close();
        if (thisfqualfile != "") {
            inFQual.close();
            inRQual.close();
            outQual.close();
            outScrapQual.close();
        }
        delete alignment;
        
        if (m->control_pressed) {  m->mothurRemove(outputFasta); m->mothurRemove(outputScrapFasta);m->mothurRemove(outputMisMatches); if (thisfqualfile != "") { m->mothurRemove(outputQual); m->mothurRemove(outputScrapQual); } }
    
        return num;
    }
	catch(exception& e) {
		m->errorOut(e, "MakeContigsCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
vector< vector<string> > MakeContigsCommand::readFastqFiles(unsigned long int& count, string ffastq, string rfastq){
    try {
        vector< vector<string> > files;
        //maps processors number to file pointer
        map<int, vector<ofstream*> > tempfiles;  //tempfiles[0] = forwardFasta, [1] = forwardQual, [2] = reverseFasta, [3] = reverseQual
        map<int, vector<ofstream*> >::iterator it;
        
        //create files to write to
        for (int i = 0; i < processors; i++) {
            vector<ofstream*> temp;
            ofstream* outFF = new ofstream;     temp.push_back(outFF);
            ofstream* outFQ = new ofstream;     temp.push_back(outFQ);
            ofstream* outRF = new ofstream;     temp.push_back(outRF);
            ofstream* outRQ = new ofstream;     temp.push_back(outRQ);
            tempfiles[i] = temp;
            
            vector<string> names;
            string thisOutputDir = outputDir;
            if (outputDir == "") { thisOutputDir = m->hasPath(ffastq); }
            string ffastafilename = thisOutputDir + m->getRootName(m->getSimpleName(ffastq)) + toString(i) + "ffastatemp";
            string rfastafilename = thisOutputDir + m->getRootName(m->getSimpleName(rfastq)) + toString(i) + "rfastatemp";
            string fqualfilename = thisOutputDir + m->getRootName(m->getSimpleName(ffastq)) + toString(i) + "fqualtemp";
            string rqualfilename = thisOutputDir + m->getRootName(m->getSimpleName(rfastq)) + toString(i) + "rqualtemp";
            names.push_back(ffastafilename); names.push_back(fqualfilename);
            names.push_back(rfastafilename); names.push_back(rqualfilename);
            files.push_back(names);
            
            m->openOutputFile(ffastafilename, *outFF);
            m->openOutputFile(rfastafilename, *outRF);
            m->openOutputFile(fqualfilename, *outFQ);
            m->openOutputFile(rqualfilename, *outRQ);
        }
        
        if (m->control_pressed) {
            //close files, delete ofstreams
            for (it = tempfiles.begin(); it!=tempfiles.end(); it++) { for (int i = 0; i < (it->second).size(); i++) { (*(it->second)[i]).close();  delete (it->second)[i]; } }
            //remove files
            for (int i = 0; i < files.size(); i++) {  
                for(int j = 0; j < files[i].size(); j++) { m->mothurRemove(files[i][j]); }
            }
        }
        
        ifstream inForward;
        m->openInputFile(ffastq, inForward);
        
        ifstream inReverse;
        m->openInputFile(rfastq, inReverse);
        
        count = 0;
        map<string, fastqRead> uniques;
        map<string, fastqRead>::iterator itUniques;
        while ((!inForward.eof()) || (!inReverse.eof())) {
            
            if (m->control_pressed) { for (it = tempfiles.begin(); it!=tempfiles.end(); it++) { for (int i = 0; i < (it->second).size(); i++) { (*(it->second)[i]).close();  delete (it->second)[i]; } } for (int i = 0; i < files.size(); i++) {  for(int j = 0; j < files[i].size(); j++) { m->mothurRemove(files[i][j]); } } inForward.close(); inReverse.close(); return files; }
            
            //get a read from forward and reverse fastq files
            bool ignoref, ignorer;
            fastqRead thisFread, thisRread;
            if (!inForward.eof()) {  thisFread = readFastq(inForward, ignoref); }
            else { ignoref = true; }
            if (!inReverse.eof()) { thisRread = readFastq(inReverse, ignorer);  }
            else { ignorer = true; }
            
            vector<pairFastqRead> reads = getReads(ignoref, ignorer, thisFread, thisRread, uniques);
           
            for (int i = 0; i < reads.size(); i++) {
                fastqRead fread = reads[i].forward;
                fastqRead rread = reads[i].reverse;
                
                if (m->debug) { m->mothurOut(toString(count) + '\t' + fread.name + '\t' + rread.name + '\n'); }
               
                if (checkReads(fread, rread, ffastq, rfastq)) {
                    if (m->control_pressed) { for (it = tempfiles.begin(); it!=tempfiles.end(); it++) { for (int i = 0; i < (it->second).size(); i++) { (*(it->second)[i]).close();  delete (it->second)[i]; } } for (int i = 0; i < files.size(); i++) {  for(int j = 0; j < files[i].size(); j++) { m->mothurRemove(files[i][j]); } } inForward.close(); inReverse.close(); return files; }
                    
                    //if the reads are okay write to output files
                    int process = count % processors;
                    
                    *(tempfiles[process][0]) << ">" << fread.name << endl << fread.sequence << endl;
                    *(tempfiles[process][1]) << ">" << fread.name << endl;
                    for (int i = 0; i < fread.scores.size(); i++) { *(tempfiles[process][1]) << fread.scores[i] << " "; }
                    *(tempfiles[process][1]) << endl;
                    *(tempfiles[process][2]) << ">" << rread.name << endl << rread.sequence << endl;
                    *(tempfiles[process][3]) << ">" << rread.name << endl;
                    for (int i = 0; i < rread.scores.size(); i++) { *(tempfiles[process][3]) << rread.scores[i] << " "; }
                    *(tempfiles[process][3]) << endl;
                    
                    count++;
                    
                    //report progress
                    if((count) % 10000 == 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}
                }
            }
		}
		//report progress
		if((count) % 10000 != 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}
        
        if (uniques.size() != 0) {
            for (itUniques = uniques.begin(); itUniques != uniques.end(); itUniques++) {
                m->mothurOut("[WARNING]: did not find paired read for " + itUniques->first + ", ignoring.\n");
            }
            m->mothurOutEndLine();
        }
        
        //close files, delete ofstreams
        for (it = tempfiles.begin(); it!=tempfiles.end(); it++) { for (int i = 0; i < (it->second).size(); i++) { (*(it->second)[i]).close();  delete (it->second)[i]; } }
        inForward.close();
        inReverse.close();
        
        return files;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "readFastqFiles");
        exit(1);
    }
}
//**********************************************************************************************************************
vector< vector<string> > MakeContigsCommand::readFastaFiles(unsigned long int& count, string ffasta, string rfasta){
    try {
        vector< vector<string> > files;
        //maps processors number to file pointer
        map<int, vector<ofstream*> > tempfiles;  //tempfiles[0] = forwardFasta, [1] = forwardQual, [2] = reverseFasta, [3] = reverseQual
        map<int, vector<ofstream*> >::iterator it;
        
        //create files to write to
        for (int i = 0; i < processors; i++) {
            vector<ofstream*> temp;
            ofstream* outFF = new ofstream;     temp.push_back(outFF);
            ofstream* outFQ = new ofstream;     temp.push_back(outFQ);
            ofstream* outRF = new ofstream;     temp.push_back(outRF);
            ofstream* outRQ = new ofstream;     temp.push_back(outRQ);
            tempfiles[i] = temp;
            
            vector<string> names;
            string thisOutputDir = outputDir;
            if (outputDir == "") { thisOutputDir = m->hasPath(ffasta); }
            string ffastafilename = thisOutputDir + m->getRootName(m->getSimpleName(ffasta)) + toString(i) + "ffastatemp";
            string rfastafilename = thisOutputDir + m->getRootName(m->getSimpleName(rfasta)) + toString(i) + "rfastatemp";
            string fqualfilename = "";
            if (fqualfile != "") { fqualfilename = thisOutputDir + m->getRootName(m->getSimpleName(fqualfile)) + toString(i) + "fqual.temp";  m->openOutputFile(fqualfilename, *outFQ); }
            string rqualfilename = "";
            if (rqualfile != "") { rqualfilename = thisOutputDir + m->getRootName(m->getSimpleName(rqualfile)) + toString(i) + "rqual.temp"; m->openOutputFile(rqualfilename, *outRQ); }
            names.push_back(ffastafilename); names.push_back(fqualfilename);
            names.push_back(rfastafilename); names.push_back(rqualfilename);
            files.push_back(names);
            
            m->openOutputFile(ffastafilename, *outFF);
            m->openOutputFile(rfastafilename, *outRF);
        }
        
        if (m->control_pressed) {
            //close files, delete ofstreams
            for (it = tempfiles.begin(); it!=tempfiles.end(); it++) { for (int i = 0; i < (it->second).size(); i++) { (*(it->second)[i]).close();  delete (it->second)[i]; } }
            //remove files
            for (int i = 0; i < files.size(); i++) {  
                for(int j = 0; j < files[i].size(); j++) { m->mothurRemove(files[i][j]); }
            }
        }
        
        ifstream inForwardFasta;
        m->openInputFile(ffasta, inForwardFasta);
        
        ifstream inReverseFasta;
        m->openInputFile(rfasta, inReverseFasta);
        
        ifstream inForwardQual, inReverseQual;
        if (fqualfile != "") { m->openInputFile(fqualfile, inForwardQual); m->openInputFile(rqualfile, inReverseQual); }
        
        count = 0;
        map<string, fastqRead> uniques;
        map<string, fastqRead>::iterator itUniques;
        while ((!inForwardFasta.eof()) || (!inReverseFasta.eof())) {
            
            if (m->control_pressed) { for (it = tempfiles.begin(); it!=tempfiles.end(); it++) { for (int i = 0; i < (it->second).size(); i++) { (*(it->second)[i]).close();  delete (it->second)[i]; } } for (int i = 0; i < files.size(); i++) {  for(int j = 0; j < files[i].size(); j++) { m->mothurRemove(files[i][j]); } } inReverseFasta.close(); inForwardFasta.close(); if (fqualfile != "") { inReverseQual.close(); inReverseQual.close(); } return files; }
            
            //get a reads from forward and reverse fasta files
            bool ignoref, ignorer;
            fastqRead thisFread, thisRread;
            if (!inForwardFasta.eof()) {  
                ignoref = false; 
                Sequence temp(inForwardFasta); m->gobble(inForwardFasta);
                thisFread.name = temp.getName();
                thisFread.sequence = temp.getUnaligned();
            }else { ignoref = true; }
            if (!inReverseFasta.eof()) {  
                ignorer = false; 
                Sequence temp(inReverseFasta); m->gobble(inReverseFasta);
                thisRread.name = temp.getName();
                thisRread.sequence = temp.getUnaligned();  
            }else { ignorer = true; }
            
            //get qual reads if given
            if (fqualfile != "") {
                if (!inForwardQual.eof() && !ignoref) {  
                    QualityScores temp(inForwardQual); m->gobble(inForwardQual);
                    //if forward files dont match ignore read
                    if (thisFread.name != temp.getName()) { ignoref = true; } 
                    else { thisFread.scores = temp.getQualityScores(); }
                }else { ignoref = true; }
                if (!inReverseQual.eof() && !ignorer) {  
                    QualityScores temp(inReverseQual); m->gobble(inReverseQual);
                    //if reverse files dont match ignore read
                    if (thisRread.name != temp.getName()) { ignorer = true; } 
                    else { thisRread.scores = temp.getQualityScores(); }
                }else { ignorer = true; }
            }
            
            vector<pairFastqRead> reads = getReads(ignoref, ignorer, thisFread, thisRread, uniques);
            
            for (int i = 0; i < reads.size(); i++) {
                fastqRead fread = reads[i].forward;
                fastqRead rread = reads[i].reverse;
                
                if (m->debug) { m->mothurOut(toString(count) + '\t' + fread.name + '\t' + rread.name + '\n'); }
                
                if (checkReads(fread, rread, ffasta, rfasta)) {
                    if (m->control_pressed) { for (it = tempfiles.begin(); it!=tempfiles.end(); it++) { for (int i = 0; i < (it->second).size(); i++) { (*(it->second)[i]).close();  delete (it->second)[i]; } } for (int i = 0; i < files.size(); i++) {  for(int j = 0; j < files[i].size(); j++) { m->mothurRemove(files[i][j]); } } inReverseFasta.close(); inForwardFasta.close(); if (fqualfile != "") { inReverseQual.close(); inReverseQual.close(); } return files; }
                    
                    //if the reads are okay write to output files
                    int process = count % processors;
                    
                    *(tempfiles[process][0]) << ">" << fread.name << endl << fread.sequence << endl;
                    *(tempfiles[process][2]) << ">" << rread.name << endl << rread.sequence << endl;
                    if (fqualfile != "") { //if you have quality info, print it
                        *(tempfiles[process][1]) << ">" << fread.name << endl;
                        for (int i = 0; i < fread.scores.size(); i++) { *(tempfiles[process][1]) << fread.scores[i] << " "; }
                        *(tempfiles[process][1]) << endl;
                        *(tempfiles[process][3]) << ">" << rread.name << endl;
                        for (int i = 0; i < rread.scores.size(); i++) { *(tempfiles[process][3]) << rread.scores[i] << " "; }
                        *(tempfiles[process][3]) << endl;
                    }
                    count++;
                    
                    //report progress
                    if((count) % 10000 == 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}
                }
            }
		}
		//report progress
		if((count) % 10000 != 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}
        
        if (uniques.size() != 0) {
            for (itUniques = uniques.begin(); itUniques != uniques.end(); itUniques++) {
                m->mothurOut("[WARNING]: did not find paired read for " + itUniques->first + ", ignoring.\n");
            }
            m->mothurOutEndLine();
        }
        
        //close files, delete ofstreams
        for (it = tempfiles.begin(); it!=tempfiles.end(); it++) { for (int i = 0; i < (it->second).size(); i++) { (*(it->second)[i]).close();  delete (it->second)[i]; } }
        inReverseFasta.close(); 
        inForwardFasta.close(); 
        if (fqualfile != "") { inReverseQual.close(); inReverseQual.close(); }
        
        return files;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "readFastaFiles");
        exit(1);
    }
}
//**********************************************************************************************************************
vector<pairFastqRead> MakeContigsCommand::getReads(bool ignoref, bool ignorer, fastqRead forward, fastqRead reverse, map<string, fastqRead>& uniques){
    try {
        vector<pairFastqRead> reads;
        map<string, fastqRead>::iterator itUniques;
            
        if (!ignoref && !ignorer) {
            if (forward.name == reverse.name) { 
                pairFastqRead temp(forward, reverse);
                reads.push_back(temp);
            }else {
                bool match = false;
                //if no match are the names only different by 1 and 2?
                string tempFRead = forward.name.substr(0, forward.name.length()-1);
                string tempRRead = reverse.name.substr(0, reverse.name.length()-1);
                if (tempFRead == tempRRead) {
                    if ((forward.name[forward.name.length()-1] == '1') && (reverse.name[reverse.name.length()-1] == '2')) {
                        forward.name = tempFRead;
                        reverse.name = tempRRead;
                        pairFastqRead temp(forward, reverse);
                        reads.push_back(temp);
                        match = true;
                    }
                }
                
                if (!match) {
                    //look for forward pair
                    itUniques = uniques.find(forward.name);
                    if (itUniques != uniques.end()) {  //we have the pair for this read
                        pairFastqRead temp(forward, itUniques->second);
                        reads.push_back(temp);
                        uniques.erase(itUniques);
                    }else { //save this read for later
                        uniques[forward.name] = forward;
                    }
                    
                    //look for reverse pair
                    itUniques = uniques.find(reverse.name);
                    if (itUniques != uniques.end()) {  //we have the pair for this read
                        pairFastqRead temp(itUniques->second, reverse);
                        reads.push_back(temp);
                        uniques.erase(itUniques);
                    }else { //save this read for later
                        uniques[reverse.name] = reverse;
                    }
                }
                                
            }
        }else if (!ignoref && ignorer) { //ignore reverse keep forward
            //look for forward pair
            itUniques = uniques.find(forward.name);
            if (itUniques != uniques.end()) {  //we have the pair for this read
                pairFastqRead temp(forward, itUniques->second);
                reads.push_back(temp);
                uniques.erase(itUniques);
            }else { //save this read for later
                uniques[forward.name] = forward;
            }

        }else if (ignoref && !ignorer) { //ignore forward keep reverse
            //look for reverse pair
            itUniques = uniques.find(reverse.name);
            if (itUniques != uniques.end()) {  //we have the pair for this read
                pairFastqRead temp(itUniques->second, reverse);
                reads.push_back(temp);
                uniques.erase(itUniques);
            }else { //save this read for later
                uniques[reverse.name] = reverse;
            }
        }//else ignore both and do nothing
        
        return reads;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "readFastqFiles");
        exit(1);
    }
}
//**********************************************************************************************************************
fastqRead MakeContigsCommand::readFastq(ifstream& in, bool& ignore){
    try {
        fastqRead read;
        
        ignore = false;
        
        //read sequence name
        string line = m->getline(in); m->gobble(in);
        vector<string> pieces = m->splitWhiteSpace(line);
        string name = "";  if (pieces.size() != 0) { name = pieces[0]; }
        if (name == "") {  m->mothurOut("[WARNING]: Blank fasta name, ignoring read."); m->mothurOutEndLine(); ignore=true;  }
        else if (name[0] != '@') { m->mothurOut("[WARNING]: reading " + name + " expected a name with @ as a leading character, ignoring read."); m->mothurOutEndLine(); ignore=true; }
        else { name = name.substr(1); }
        
        //read sequence
        string sequence = m->getline(in); m->gobble(in);
        if (sequence == "") {  m->mothurOut("[WARNING]: missing sequence for " + name + ", ignoring."); ignore=true; }
        
        //read sequence name
        line = m->getline(in); m->gobble(in);
        pieces = m->splitWhiteSpace(line);
        string name2 = "";  if (pieces.size() != 0) { name2 = pieces[0]; }
        if (name2 == "") {  m->mothurOut("[WARNING]: expected a name with + as a leading character, ignoring."); ignore=true; }
        else if (name2[0] != '+') { m->mothurOut("[WARNING]: reading " + name2 + " expected a name with + as a leading character, ignoring."); ignore=true; }
        else { name2 = name2.substr(1); if (name2 == "") { name2 = name; } }
        
        //read quality scores
        string quality = m->getline(in); m->gobble(in);
        if (quality == "") {  m->mothurOut("[WARNING]: missing quality for " + name2 + ", ignoring."); ignore=true; }
         
        //sanity check sequence length and number of quality scores match
        if (name2 != "") { if (name != name2) { m->mothurOut("[WARNING]: names do not match. read " + name + " for fasta and " + name2 + " for quality, ignoring."); ignore=true; } }
        if (quality.length() != sequence.length()) { m->mothurOut("[WARNING]: Lengths do not match for sequence " + name + ". Read " + toString(sequence.length()) + " characters for fasta and " + toString(quality.length()) + " characters for quality scores, ignoring read."); ignore=true; }
        
        vector<int> qualScores;
		int controlChar = int('!');
		for (int i = 0; i < quality.length(); i++) { 
			int temp = int(quality[i]);
			temp -= controlChar;
			
			qualScores.push_back(temp);
		}
    
        read.name = name;
        read.sequence = sequence;
        read.scores = qualScores;

        return read;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "readFastq");
        exit(1);
    }
}
//**********************************************************************************************************************
bool MakeContigsCommand::checkReads(fastqRead& forward, fastqRead& reverse, string ffile, string rfile){
    try {
        bool good = true;
        
        //do sequence lengths match
        if (forward.sequence.length() != reverse.sequence.length()) {
            m->mothurOut("[WARNING]: For sequence " + forward.name + " I read a sequence of length " + toString(forward.sequence.length()) + " from " + ffile + ", but read a sequence of length " + toString(reverse.sequence.length()) + " from " + rfile + ", ignoring.\n");
            good = false; 
        }
        
        //do number of qual scores match 
        if (forward.scores.size() != reverse.scores.size()) {
            m->mothurOut("[WARNING]: For sequence " + forward.name + " I read " + toString(forward.scores.size()) + " quality scores from " + ffile + ", but read  " + toString(reverse.scores.size()) + " quality scores from " + rfile + ", ignoring.\n");
            good = false; 
        }

        return good;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "checkReads");
        exit(1);
    }
}
//***************************************************************************************************************
vector< vector<string> > MakeContigsCommand::readFileNames(string filename){
	try {
        vector< vector<string> > files;
        string forward, reverse;
        
        ifstream in;
        m->openInputFile(filename, in);
        
        while(!in.eof()) {
            
            if (m->control_pressed) { return files; }
            
            in >> forward; m->gobble(in);
            in >> reverse; m->gobble(in);
            
            //check to make sure both are able to be opened
            ifstream in2;
            int openForward = m->openInputFile(forward, in2, "noerror");
            
            //if you can't open it, try default location
            if (openForward == 1) {
                if (m->getDefaultPath() != "") { //default path is set
                    string tryPath = m->getDefaultPath() + m->getSimpleName(forward);
                    m->mothurOut("Unable to open " + forward + ". Trying default " + tryPath); m->mothurOutEndLine();
                    ifstream in3;
                    openForward = m->openInputFile(tryPath, in3, "noerror");
                    in3.close();
                    forward = tryPath;
                }
            }
            
            //if you can't open it, try output location
            if (openForward == 1) {
                if (m->getOutputDir() != "") { //default path is set
                    string tryPath = m->getOutputDir() + m->getSimpleName(forward);
                    m->mothurOut("Unable to open " + forward + ". Trying output directory " + tryPath); m->mothurOutEndLine();
                    ifstream in4;
                    openForward = m->openInputFile(tryPath, in4, "noerror");
                    forward = tryPath;
                    in4.close();
                }
            }
            
            if (openForward == 1) { //can't find it
                m->mothurOut("[WARNING]: can't find " + forward + ", ignoring pair.\n"); 
            }else{  in2.close();  }
            
            ifstream in3;
            int openReverse = m->openInputFile(reverse, in3, "noerror");
            
            //if you can't open it, try default location
            if (openReverse == 1) {
                if (m->getDefaultPath() != "") { //default path is set
                    string tryPath = m->getDefaultPath() + m->getSimpleName(reverse);
                    m->mothurOut("Unable to open " + reverse + ". Trying default " + tryPath); m->mothurOutEndLine();
                    ifstream in3;
                    openReverse = m->openInputFile(tryPath, in3, "noerror");
                    in3.close();
                    reverse = tryPath;
                }
            }
            
            //if you can't open it, try output location
            if (openReverse == 1) {
                if (m->getOutputDir() != "") { //default path is set
                    string tryPath = m->getOutputDir() + m->getSimpleName(reverse);
                    m->mothurOut("Unable to open " + reverse + ". Trying output directory " + tryPath); m->mothurOutEndLine();
                    ifstream in4;
                    openReverse = m->openInputFile(tryPath, in4, "noerror");
                    reverse = tryPath;
                    in4.close();
                }
            }
            
            if (openReverse == 1) { //can't find it
                m->mothurOut("[WARNING]: can't find " + reverse + ", ignoring pair.\n"); 
            }else{  in3.close();  }
            
            if ((openForward != 1) && (openReverse != 1)) { //good pair
                vector<string> pair;
                pair.push_back(forward);
                pair.push_back(reverse);
                files.push_back(pair);
            }
            
        }
        in.close();
        
        return files;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "checkReads");
        exit(1);
    }
}
//***************************************************************************************************************
//illumina data requires paired forward and reverse data
//BARCODE   atgcatgc   atgcatgc    groupName 
//PRIMER   atgcatgc   atgcatgc    groupName  
//PRIMER   atgcatgc   atgcatgc  
bool MakeContigsCommand::getOligos(vector<vector<string> >& fastaFileNames, vector<vector<string> >& qualFileNames, string rootname){
	try {
		ifstream in;
		m->openInputFile(oligosfile, in);
		
		ofstream test;
		
		string type, foligo, roligo, group;
        
		int indexPrimer = 0;
		int indexBarcode = 0;
        set<string> uniquePrimers;
        set<string> uniqueBarcodes;
		
		while(!in.eof()){
            
			in >> type; 
            cout << type << endl;
		 	if (m->debug) { m->mothurOut("[DEBUG]: reading type - " + type + ".\n"); }	
            
			if(type[0] == '#'){
				while (!in.eof())	{	char c = in.get();  if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
				m->gobble(in);
			}
			else{
				m->gobble(in);
				//make type case insensitive
				for(int i=0;i<type.length();i++){	type[i] = toupper(type[i]);  }
				
				in >> foligo;
                
                if (m->debug) { m->mothurOut("[DEBUG]: reading - " + foligo + ".\n"); }
				
				for(int i=0;i<foligo.length();i++){
					foligo[i] = toupper(foligo[i]);
					if(foligo[i] == 'U')	{	foligo[i] = 'T';	}
				}
				
				if(type == "FORWARD"){
					m->gobble(in);
					
                    in >> roligo;
                    
                    for(int i=0;i<roligo.length();i++){
                        roligo[i] = toupper(roligo[i]);
                        if(roligo[i] == 'U')	{	roligo[i] = 'T';	}
                    }
                    //roligo = reverseOligo(roligo);
                    
                    group = "";
                    
					// get rest of line in case there is a primer name
					while (!in.eof())	{	
						char c = in.get(); 
						if (c == 10 || c == 13){	break;	}
						else if (c == 32 || c == 9){;} //space or tab
						else { 	group += c;  }
					} 
                    
                    oligosPair newPrimer(foligo, roligo);
					
					//check for repeat barcodes
                    string tempPair = foligo+roligo;
                    if (uniquePrimers.count(tempPair) != 0) { m->mothurOut("primer pair " + newPrimer.forward + " " + newPrimer.reverse + " is in your oligos file already."); m->mothurOutEndLine();  }
                    else { uniquePrimers.insert(tempPair); }
					
                    if (m->debug) {  if (group != "") { m->mothurOut("[DEBUG]: reading group " + group + ".\n"); }else{ m->mothurOut("[DEBUG]: no group for primer pair " + newPrimer.forward + " " + newPrimer.reverse + ".\n"); }  }
                    
					primers[indexPrimer]=newPrimer; indexPrimer++;		
					primerNameVector.push_back(group);
				}else if(type == "BARCODE"){
					m->gobble(in);
					
                    in >> roligo;
                    
                    for(int i=0;i<roligo.length();i++){
                        roligo[i] = toupper(roligo[i]);
                        if(roligo[i] == 'U')	{	roligo[i] = 'T';	}
                    }
                    //roligo = reverseOligo(roligo);
                    
                    oligosPair newPair(foligo, roligo);
                    
                    group = "";
                    while (!in.eof())	{	
						char c = in.get(); 
						if (c == 10 || c == 13){	break;	}
						else if (c == 32 || c == 9){;} //space or tab
						else { 	group += c;  }
					} 
					
                    if (m->debug) { m->mothurOut("[DEBUG]: barcode pair " + newPair.forward + " " + newPair.reverse + ", and group = " + group + ".\n"); }
                        
                    //check for repeat barcodes
                    string tempPair = foligo+roligo;
                    if (uniqueBarcodes.count(tempPair) != 0) { m->mothurOut("barcode pair " + newPair.forward + " " + newPair.reverse +  " is in your oligos file already, disregarding."); m->mothurOutEndLine();  }
                    else { uniqueBarcodes.insert(tempPair); }
                        
                    barcodes[indexBarcode]=newPair; indexBarcode++;
					barcodeNameVector.push_back(group);
                    cout << group << endl;
				}else if(type == "LINKER"){
					linker.push_back(foligo);
                    m->mothurOut("[WARNING]: make.contigs is not setup to remove linkers, ignoring.\n");
				}else if(type == "SPACER"){
					spacer.push_back(foligo);
                    m->mothurOut("[WARNING]: make.contigs is not setup to remove spacers, ignoring.\n");
				}
				else{	m->mothurOut("[WARNING]: " + type + " is not recognized as a valid type. Choices are primer, barcode, linker and spacer. Ignoring " + foligo + "."); m->mothurOutEndLine(); }
			}
			m->gobble(in);
		}	
		in.close();
		
		if(barcodeNameVector.size() == 0 && primerNameVector[0] == ""){	allFiles = 0;	}
		
		//add in potential combos
		if(barcodeNameVector.size() == 0){
            oligosPair temp("", "");
			barcodes[0] = temp;
			barcodeNameVector.push_back("");			
		}
		
		if(primerNameVector.size() == 0){
            oligosPair temp("", "");
			primers[0] = temp;
			primerNameVector.push_back("");			
		}
		
		fastaFileNames.resize(barcodeNameVector.size());
		for(int i=0;i<fastaFileNames.size();i++){
			fastaFileNames[i].assign(primerNameVector.size(), "");
		}
		qualFileNames = fastaFileNames;	
		
		if(allFiles){
			set<string> uniqueNames; //used to cleanup outputFileNames
			for(map<int, oligosPair>::iterator itBar = barcodes.begin();itBar != barcodes.end();itBar++){
				for(map<int, oligosPair>::iterator itPrimer = primers.begin();itPrimer != primers.end(); itPrimer++){
					
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
                        fastaFileName = rootname + comboGroupName + ".fasta";
                        if (uniqueNames.count(fastaFileName) == 0) {
                            outputNames.push_back(fastaFileName);
                            outputTypes["fasta"].push_back(fastaFileName);
                            uniqueNames.insert(fastaFileName);
                        }
                        
                        fastaFileNames[itBar->first][itPrimer->first] = fastaFileName;
                        m->openOutputFile(fastaFileName, temp);		temp.close();
                        
                        if ((fqualfile != "") || (ffastqfile != "") || (file != "")) {
                            qualFileName = rootname + ".qual";
                            if (uniqueNames.count(qualFileName) == 0) {
                                outputNames.push_back(qualFileName);
                                outputTypes["qfile"].push_back(qualFileName);
                            }
                            
                            qualFileNames[itBar->first][itPrimer->first] = qualFileName;
                            m->openOutputFile(qualFileName, temp);		temp.close();
                        }
                    }
				}
			}
		}
		
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
		m->errorOut(e, "MakeContigsCommand", "getOligos");
		exit(1);
	}
}
//********************************************************************/
string MakeContigsCommand::reverseOligo(string oligo){
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
		m->errorOut(e, "MakeContigsCommand", "reverseOligo");
		exit(1);
	}
}
//**********************************************************************************************************************




