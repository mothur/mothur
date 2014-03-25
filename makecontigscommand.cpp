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
        CommandParameter pfqual("fqfile", "InputTypes", "", "", "none", "none", "qfileGroup","",false,false,true); parameters.push_back(pfqual);
        CommandParameter prqual("rqfile", "InputTypes", "", "", "none", "none", "qfileGroup","",false,false,true); parameters.push_back(prqual);
        CommandParameter pfile("file", "InputTypes", "", "", "FastaFastqFile", "FastaFastqFile", "none","fasta-qfile",false,false,true); parameters.push_back(pfile);
        CommandParameter poligos("oligos", "InputTypes", "", "", "none", "none", "none","group",false,false,true); parameters.push_back(poligos);
        CommandParameter pfindex("findex", "InputTypes", "", "", "none", "none", "none","",false,false,true); parameters.push_back(pfindex);
        CommandParameter prindex("rindex", "InputTypes", "", "", "none", "none", "none","",false,false,true); parameters.push_back(prindex);
		CommandParameter ppdiffs("pdiffs", "Number", "", "0", "", "", "","",false,false,true); parameters.push_back(ppdiffs);
		CommandParameter pbdiffs("bdiffs", "Number", "", "0", "", "", "","",false,false,true); parameters.push_back(pbdiffs);
        CommandParameter ptdiffs("tdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(ptdiffs);

        CommandParameter palign("align", "Multiple", "needleman-gotoh", "needleman", "", "", "","",false,false); parameters.push_back(palign);
        CommandParameter pallfiles("allfiles", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pallfiles);
        CommandParameter ptrimoverlap("trimoverlap", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(ptrimoverlap);
		CommandParameter pmatch("match", "Number", "", "1.0", "", "", "","",false,false); parameters.push_back(pmatch);
		CommandParameter pmismatch("mismatch", "Number", "", "-1.0", "", "", "","",false,false); parameters.push_back(pmismatch);
		CommandParameter pgapopen("gapopen", "Number", "", "-2.0", "", "", "","",false,false); parameters.push_back(pgapopen);
		CommandParameter pgapextend("gapextend", "Number", "", "-1.0", "", "", "","",false,false); parameters.push_back(pgapextend);
        CommandParameter pthreshold("insert", "Number", "", "20", "", "", "","",false,false); parameters.push_back(pthreshold);
        CommandParameter pdeltaq("deltaq", "Number", "", "6", "", "", "","",false,false); parameters.push_back(pdeltaq);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
        CommandParameter pformat("format", "Multiple", "sanger-illumina-solexa-illumina1.8+", "illumina1.8+", "", "", "","",false,false,true); parameters.push_back(pformat);
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
		helpString += "The make.contigs command reads a file, forward fastq file and a reverse fastq file or forward fasta and reverse fasta files and outputs new fasta. \n";
        helpString += "If an oligos file is provided barcodes and primers will be trimmed, and a group file will be created.\n";
        helpString += "If a forward index or reverse index file is provided barcodes be trimmed, and a group file will be created. The oligos parameter is required if an index file is given.\n";
		helpString += "The make.contigs command parameters are file, ffastq, rfastq, ffasta, rfasta, fqfile, rqfile, oligos, findex, rindex, format, tdiffs, bdiffs, pdiffs, align, match, mismatch, gapopen, gapextend, insert, deltaq, allfiles and processors.\n";
		helpString += "The ffastq and rfastq, file, or ffasta and rfasta parameters are required.\n";
        helpString += "The file parameter is 2, 3 or 4 column file containing the forward fastq files in the first column and their matching reverse fastq files in the second column, or a groupName then forward fastq file and reverse fastq file, or forward fastq file then reverse fastq then forward index and reverse index file.  If you only have one index file add 'none' for the other one.  Mothur will process each pair and create a combined fasta and report file with all the sequences.\n";
        helpString += "The ffastq and rfastq parameters are used to provide a forward fastq and reverse fastq file to process.  If you provide one, you must provide the other.\n";
        helpString += "The ffasta and rfasta parameters are used to provide a forward fasta and reverse fasta file to process.  If you provide one, you must provide the other.\n";
        helpString += "The fqfile and rqfile parameters are used to provide a forward quality and reverse quality files to process with the ffasta and rfasta parameters.  If you provide one, you must provide the other.\n";
		helpString += "The format parameter is used to indicate whether your sequences are sanger, solexa, illumina1.8+ or illumina, default=illumina1.8+.\n";
        helpString += "The findex and rindex parameters are used to provide a forward index and reverse index files to process.  \n";
        helpString += "The align parameter allows you to specify the alignment method to use.  Your options are: gotoh and needleman. The default is needleman.\n";
        helpString += "The tdiffs parameter is used to specify the total number of differences allowed in the sequence. The default is pdiffs + bdiffs + sdiffs + ldiffs.\n";
		helpString += "The bdiffs parameter is used to specify the number of differences allowed in the barcode. The default is 0.\n";
		helpString += "The pdiffs parameter is used to specify the number of differences allowed in the primer. The default is 0.\n";
        //helpString += "The ldiffs parameter is used to specify the number of differences allowed in the linker. The default is 0.\n";
		//helpString += "The sdiffs parameter is used to specify the number of differences allowed in the spacer. The default is 0.\n";
		helpString += "The match parameter allows you to specify the bonus for having the same base. The default is 1.0.\n";
		helpString += "The mistmatch parameter allows you to specify the penalty for having different bases.  The default is -1.0.\n";
        helpString += "The deltaq parameter allows you to specify the delta allowed between quality scores of a mismatched base.  For example in the overlap, if deltaq=5 and in the alignment seqA, pos 200 has a quality score of 30 and the same position in seqB has a quality score of 20, you take the base from seqA (30-20 >= 5).  If the quality score in seqB is 28 then the base in the consensus will be an N (30-28<5) The default is 6.\n";
		helpString += "The gapopen parameter allows you to specify the penalty for opening a gap in an alignment. The default is -2.0.\n";
		helpString += "The gapextend parameter allows you to specify the penalty for extending a gap in an alignment.  The default is -1.0.\n";
        helpString += "The insert parameter allows you to set a quality scores threshold. In the case where we are trying to decide whether to keep a base or remove it because the base is compared to a gap in the other fragment, if the base has a quality score equal to or below the threshold we eliminate it. Default=20.\n";
        helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n";
        helpString += "The allfiles parameter will create separate group and fasta file for each grouping. The default is F.\n";

        helpString += "The trimoverlap parameter allows you to trim the sequences to only the overlapping section. The default is F.\n";
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
        else if (type == "group") {  pattern = "[filename],[tag],contigs.groups"; }
        else if (type == "report") {  pattern = "[filename],[tag],contigs.report"; }
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
        outputTypes["group"] = tempOutNames;
        outputTypes["report"] = tempOutNames;
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
        createFileGroup = false; createOligosGroup = false;
        
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
            outputTypes["report"] = tempOutNames;
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
                
                it = parameters.find("findex");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["findex"] = inputDir + it->second;		}
				}
                
                it = parameters.find("rindex");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["rindex"] = inputDir + it->second;		}
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
            
            findexfile = validParameter.validFile(parameters, "findex", true);
			if (findexfile == "not found")      {	findexfile = "";	}
			else if(findexfile == "not open")   {	abort = true;       }
            
            rindexfile = validParameter.validFile(parameters, "rindex", true);
			if (rindexfile == "not found")      {	rindexfile = "";	}
			else if(rindexfile == "not open")   {	abort = true;       }
        
            if ((rindexfile != "") || (findexfile != "")) {
				if (oligosfile == ""){
					oligosfile = m->getOligosFile();
					if (oligosfile != "") {  m->mothurOut("Using " + oligosfile + " as input file for the oligos parameter.\n");  }
					else {
                        m->mothurOut("You need to provide an oligos file if you are going to use an index file.\n"); abort = true;
					}
				}
                
                //can only use an index file with the fastq parameters not fasta and qual
                if ((ffastafile != "") || (rfastafile != "")) {
                    m->mothurOut("[ERROR]: You can only use an index file with the fastq parameters or the file option.\n"); abort = true;
                }
			}
            
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
			
            temp = validParameter.validFile(parameters, "insert", false);	if (temp == "not found"){	temp = "20";			}
			m->mothurConvert(temp, insert); 
            if ((insert < 0) || (insert > 40)) { m->mothurOut("[ERROR]: insert must be between 0 and 40.\n"); abort=true; }

            temp = validParameter.validFile(parameters, "deltaq", false);	if (temp == "not found"){	temp = "6";			}
			m->mothurConvert(temp, deltaq);
            
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
            
            
            temp = validParameter.validFile(parameters, "trimoverlap", false);		if (temp == "not found") { temp = "F"; }
			trimOverlap = m->isTrue(temp);
			
			align = validParameter.validFile(parameters, "align", false);		if (align == "not found"){	align = "needleman";	}
			if ((align != "needleman") && (align != "gotoh")) { m->mothurOut(align + " is not a valid alignment method. Options are needleman or gotoh. I will use needleman."); m->mothurOutEndLine(); align = "needleman"; }
            
            format = validParameter.validFile(parameters, "format", false);		if (format == "not found"){	format = "illumina1.8+";	}
            
            if ((format != "sanger") && (format != "illumina") && (format != "illumina1.8+") && (format != "solexa"))  { 
				m->mothurOut(format + " is not a valid format. Your format choices are sanger, solexa, illumina1.8+ and illumina, aborting." ); m->mothurOutEndLine();
				abort=true;
			}
            
            //fill convert table - goes from solexa to sanger. Used fq_all2std.pl as a reference.
            for (int i = -64; i < 65; i++) { 
                char temp = (char) ((int)(33 + 10*log(1+pow(10,(i/10.0)))/log(10)+0.499));
                convertTable.push_back(temp);
            }
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
        cvars["[tag]"] = "scrap";
        string compositeScrapFastaFile = getOutputFileName("fasta",cvars);
        cvars["[tag]"] = "";
        string compositeMisMatchFile = getOutputFileName("report",cvars);
        
        if (filesToProcess.size() > 1) { //clear files for append below
            ofstream outCTFasta, outCTQual, outCSFasta, outCSQual, outCMisMatch;
            m->openOutputFile(compositeFastaFile, outCTFasta); outCTFasta.close();
            m->openOutputFile(compositeScrapFastaFile, outCSFasta); outCSFasta.close();
            m->openOutputFile(compositeMisMatchFile, outCMisMatch); outCMisMatch.close();
            outputNames.push_back(compositeFastaFile); outputTypes["fasta"].push_back(compositeFastaFile);
            outputNames.push_back(compositeMisMatchFile); outputTypes["report"].push_back(compositeMisMatchFile);
            outputNames.push_back(compositeScrapFastaFile); outputTypes["fasta"].push_back(compositeScrapFastaFile);
        }
        
        map<string, int> totalGroupCounts;
        
        for (int l = 0; l < filesToProcess.size(); l++) {
            
            m->mothurOut("\n>>>>>\tProcessing " + filesToProcess[l][0][0] + " (file " + toString(l+1) + " of " + toString(filesToProcess.size()) + ")\t<<<<<\n");
            
            groupCounts.clear();
            groupMap.clear();
            vector<vector<string> > fastaFileNames;
            createOligosGroup = false;
            string outputGroupFileName;
            map<string, string> variables; 
            string thisOutputDir = outputDir;
            if (outputDir == "") {  thisOutputDir = m->hasPath(filesToProcess[l][0][0]); }
            variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(filesToProcess[l][0][0]));
            variables["[tag]"] = "";
            if(oligosfile != ""){  createOligosGroup = getOligos(fastaFileNames, variables["[filename]"]);  }
            if (createOligosGroup || createFileGroup) {
                outputGroupFileName = getOutputFileName("group",variables);
            }
            
            //give group in file file precedence
            if (createFileGroup) {  createOligosGroup = false; }
            
            variables["[tag]"] = "trim";
            string outFastaFile = getOutputFileName("fasta",variables);
            variables["[tag]"] = "scrap";
            string outScrapFastaFile = getOutputFileName("fasta",variables);
            variables["[tag]"] = "";
            string outMisMatchFile = getOutputFileName("report",variables);
                        
            m->mothurOut("Making contigs...\n"); 
            createProcesses(filesToProcess[l], outFastaFile, outScrapFastaFile, outMisMatchFile, fastaFileNames, l);
            
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
                    out.close();
                    in.close();
                }
            }
            
            if (createFileGroup || createOligosGroup) {
                ofstream outGroup;
                m->openOutputFile(outputGroupFileName, outGroup);
                for (map<string, string>::iterator itGroup = groupMap.begin(); itGroup != groupMap.end(); itGroup++) {
                    outGroup << itGroup->first << '\t' << itGroup->second << endl;
                }
                outGroup.close();
            }
            
            if (filesToProcess.size() > 1) { //merge into large combo files
                if (createFileGroup || createOligosGroup) {
                    if (l == 0) {
                        ofstream outCGroup;
                        m->openOutputFile(compositeGroupFile, outCGroup); outCGroup.close();
                        outputNames.push_back(compositeGroupFile); outputTypes["group"].push_back(compositeGroupFile);
                    }
                    m->appendFiles(outputGroupFileName, compositeGroupFile);
                    if (!allFiles) { m->mothurRemove(outputGroupFileName);  }
                    else { outputNames.push_back(outputGroupFileName); outputTypes["group"].push_back(outputGroupFileName); }
                    
                    for (map<string, int>::iterator itGroups = groupCounts.begin(); itGroups != groupCounts.end(); itGroups++) {
                        map<string, int>::iterator itTemp = totalGroupCounts.find(itGroups->first);
                        if (itTemp == totalGroupCounts.end()) { totalGroupCounts[itGroups->first] = itGroups->second; } //new group create it in totalGroups
                        else { itTemp->second += itGroups->second; } //existing group, update total
                    }
                }
                if (l == 0) {  m->appendFiles(outMisMatchFile, compositeMisMatchFile);  }
                else {  m->appendFilesWithoutHeaders(outMisMatchFile, compositeMisMatchFile);  }
                m->appendFiles(outFastaFile, compositeFastaFile);
                m->appendFiles(outScrapFastaFile, compositeScrapFastaFile);
                if (!allFiles) {
                    m->mothurRemove(outMisMatchFile);
                    m->mothurRemove(outFastaFile);
                    m->mothurRemove(outScrapFastaFile);
                }else {
                    outputNames.push_back(outFastaFile); outputTypes["fasta"].push_back(outFastaFile);
                    outputNames.push_back(outScrapFastaFile); outputTypes["fasta"].push_back(outScrapFastaFile);
                    outputNames.push_back(outMisMatchFile); outputTypes["report"].push_back(outMisMatchFile);
                }
            }else {
                totalGroupCounts = groupCounts;
                outputNames.push_back(outFastaFile); outputTypes["fasta"].push_back(outFastaFile);
                outputNames.push_back(outScrapFastaFile); outputTypes["fasta"].push_back(outScrapFastaFile);
                outputNames.push_back(outMisMatchFile); outputTypes["report"].push_back(outMisMatchFile);
                if (createFileGroup || createOligosGroup) {
                     outputNames.push_back(outputGroupFileName); outputTypes["group"].push_back(outputGroupFileName); 
                }
            }
            m->mothurOut("Done.\n");
        }
        
        m->mothurOut("It took " + toString(time(NULL) - start) + " secs to process " + toString(numReads) + " sequences.\n");
        
        if (m->control_pressed) {	for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0;	}
        
		//output group counts
		m->mothurOutEndLine();
		int total = 0;
		if (totalGroupCounts.size() != 0) {  m->mothurOut("Group count: \n");  }
		for (map<string, int>::iterator it = totalGroupCounts.begin(); it != totalGroupCounts.end(); it++) {
            total += it->second; m->mothurOut(it->first + "\t" + toString(it->second)); m->mothurOutEndLine(); 
		}
		if (total != 0) { m->mothurOut("Total of all groups is " + toString(total)); m->mothurOutEndLine(); }
		
		if (m->control_pressed) {	for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0;	}
        
        string currentFasta = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentFasta = (itTypes->second)[0]; m->setFastaFile(currentFasta); }
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
            vector< vector<string> > files = readFastqFiles(numReads, ffastqfile, rfastqfile, findexfile, rindexfile);
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
                    vector< vector<string> > files = readFastqFiles(thisFilesReads, filePairsToProcess[i][0], filePairsToProcess[i][1], filePairsToProcess[i][2], filePairsToProcess[i][3]);
                    
                    //adjust for really large processors or really small files
                    if (thisFilesReads < processors) { 
                        m->mothurOut("[ERROR]: " + filePairsToProcess[i][0] + " has less than " + toString(processors) + " good reads, skipping\n"); 
                        for (int k = 0; k < files.size(); k++) { for(int j = 0; j < files[k].size(); j++) { m->mothurRemove(files[k][j]); } files[k].clear(); }
                        //remove from file2Group if necassary
                        map<int, string> cFile2Group;
                        for (map<int, string>::iterator it = file2Group.begin(); it != file2Group.end(); it++) {
                            if ((it->first) < i) { cFile2Group[it->first] = it->second; }
                            else if ((it->first) == i) { } //do nothing, we removed files for i
                            else { cFile2Group[(it->first-1)] = it->second; } //adjust files because i was removed
                        }
                        file2Group = cFile2Group;
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
int MakeContigsCommand::createProcesses(vector< vector<string> > files, string outputFasta, string outputScrapFasta, string outputMisMatches, vector<vector<string> > fastaFileNames, int index) {
	try {
		int num = 0;
		vector<int> processIDS;
        string group = "";
        map<int, string>::iterator it = file2Group.find(index);
        if (it != file2Group.end()) { group = it->second; }
        
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
                
				if(allFiles){
					ofstream temp;
                    
					for(int i=0;i<tempFASTAFileNames.size();i++){
						for(int j=0;j<tempFASTAFileNames[i].size();j++){
							if (tempFASTAFileNames[i][j] != "") {
								tempFASTAFileNames[i][j] += m->mothurGetpid(process) + ".temp";
								m->openOutputFile(tempFASTAFileNames[i][j], temp);			temp.close();
							}
						}
					}
				}
                
				num = driver(files[process], 
                             outputFasta + m->mothurGetpid(process) + ".temp",
                             outputScrapFasta + m->mothurGetpid(process) + ".temp",
                             outputMisMatches + m->mothurGetpid(process) + ".temp",
                             tempFASTAFileNames, process, group);
				
				//pass groupCounts to parent
                ofstream out;
                string tempFile = m->mothurGetpid(process) + ".num.temp";
                m->openOutputFile(tempFile, out);
                out << num << endl;
				if (createFileGroup || createOligosGroup) {
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
                
		//do my part
		num = driver(files[processors-1], outputFasta, outputScrapFasta,  outputMisMatches, fastaFileNames, processors-1, group);
		
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
            
			if (createFileGroup || createOligosGroup) {
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
                        
            if(allFiles){
                ofstream temp;
                
                for(int i=0;i<tempFASTAFileNames.size();i++){
                    for(int j=0;j<tempFASTAFileNames[i].size();j++){
                        if (tempFASTAFileNames[i][j] != "") {
                            tempFASTAFileNames[i][j] += extension;
                            m->openOutputFile(tempFASTAFileNames[i][j], temp);			temp.close();
                        }
                    }
                }
            }
			         
			contigsData* tempcontig = new contigsData(group, files[h], (outputFasta + extension), (outputScrapFasta + extension), (outputMisMatches + extension), align, m, match, misMatch, gapOpen, gapExtend, insert, deltaq, barcodes, primers, tempFASTAFileNames, barcodeNameVector, primerNameVector, pdiffs, bdiffs, tdiffs, createOligosGroup, createFileGroup, allFiles, trimOverlap, h);
			pDataArray.push_back(tempcontig);
            
			hThreadArray[h] = CreateThread(NULL, 0, MyContigsThreadFunction, pDataArray[h], 0, &dwThreadIdArray[h]);   
		}
        
        vector<vector<string> > tempFASTAFileNames = fastaFileNames;

        if(allFiles){
            ofstream temp;
            string extension = toString(processors-1) + ".temp";
            
            for(int i=0;i<tempFASTAFileNames.size();i++){
                for(int j=0;j<tempFASTAFileNames[i].size();j++){
                    if (tempFASTAFileNames[i][j] != "") {
                        tempFASTAFileNames[i][j] += extension;
                        m->openOutputFile(tempFASTAFileNames[i][j], temp);			temp.close();
                    }
                }
            }
        }

		//parent do my part
		ofstream temp;
		m->openOutputFile(outputFasta, temp);		temp.close();
        m->openOutputFile(outputScrapFasta, temp);		temp.close();
        
        //do my part
        processIDS.push_back(processors-1);
		num = driver(files[processors-1], (outputFasta+ toString(processors-1) + ".temp"),  (outputScrapFasta+ toString(processors-1) + ".temp"),  (outputMisMatches+ toString(processors-1) + ".temp"), tempFASTAFileNames, processors-1, group);
        
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
			num += pDataArray[i]->count;
            if (!pDataArray[i]->done) {
                m->mothurOut("[ERROR]: process " + toString(i) + " only processed " + toString(pDataArray[i]->count) + " of sequences assigned to it, quitting. \n"); m->control_pressed = true; 
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
        
        for (int i = 0; i < processIDS.size(); i++) {
			m->appendFiles((outputFasta + toString(processIDS[i]) + ".temp"), outputFasta);
			m->mothurRemove((outputFasta + toString(processIDS[i]) + ".temp"));
			
			m->appendFiles((outputScrapFasta + toString(processIDS[i]) + ".temp"), outputScrapFasta);
			m->mothurRemove((outputScrapFasta + toString(processIDS[i]) + ".temp"));
            
            m->appendFilesWithoutHeaders((outputMisMatches + toString(processIDS[i]) + ".temp"), outputMisMatches);
			m->mothurRemove((outputMisMatches + toString(processIDS[i]) + ".temp"));
            
            if(allFiles){
				for(int j=0;j<fastaFileNames.size();j++){
					for(int k=0;k<fastaFileNames[j].size();k++){
						if (fastaFileNames[j][k] != "") {
							m->appendFiles((fastaFileNames[j][k] + toString(processIDS[i]) + ".temp"), fastaFileNames[j][k]);
							m->mothurRemove((fastaFileNames[j][k] + toString(processIDS[i]) + ".temp"));
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
int MakeContigsCommand::driver(vector<string> files, string outputFasta, string outputScrapFasta, string outputMisMatches, vector<vector<string> > fastaFileNames, int process, string group){
    try {
        
        Alignment* alignment;
        if(align == "gotoh")			{	alignment = new GotohOverlap(gapOpen, gapExtend, match, misMatch, longestBase);			}
		else if(align == "needleman")	{	alignment = new NeedlemanOverlap(gapOpen, match, misMatch, longestBase);				}
        
        int num = 0;
        string thisffastafile = files[0];
        string thisfqualfile = files[1];
        string thisrfastafile = files[2];
        string thisrqualfile = files[3];
        string thisfindexfile = files[4];
        string thisrindexfile = files[5];
        
        if (m->debug) {  m->mothurOut("[DEBUG]: ffasta = " + thisffastafile + ".\n[DEBUG]: fqual = " + thisfqualfile + ".\n[DEBUG]: rfasta = " + thisrfastafile + ".\n[DEBUG]: rqual = " + thisrqualfile + ".\n[DEBUG]: findex = " + thisfindexfile + ".\n[DEBUG]: rindex = " + thisrindexfile + ".\n"); }
        
        ifstream inFFasta, inRFasta, inFQual, inRQual, inFIndex, inRIndex;
        ofstream outFasta, outMisMatch, outScrapFasta;
        m->openInputFile(thisffastafile, inFFasta);
        m->openInputFile(thisrfastafile, inRFasta);
        if (thisfqualfile != "") {
            m->openInputFile(thisfqualfile, inFQual);
            m->openInputFile(thisrqualfile, inRQual);
        }
        
        if (thisfindexfile != "") { m->openInputFile(thisfindexfile, inFIndex);  }
        if (thisrindexfile != "") { m->openInputFile(thisrindexfile, inRIndex);  }
        
        m->openOutputFile(outputFasta, outFasta);
        m->openOutputFile(outputScrapFasta, outScrapFasta);
        m->openOutputFile(outputMisMatches, outMisMatch);
        outMisMatch << "Name\tLength\tOverlap_Length\tOverlap_Start\tOverlap_End\tMisMatches\tNum_Ns\n";  
        
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
            Sequence findexBarcode("findex", "NONE");  Sequence rindexBarcode("rindex", "NONE");
            if (thisfindexfile != "") {
                Sequence temp(inFIndex); m->gobble(inFIndex);
                findexBarcode.setAligned(temp.getAligned());
            }
            
            if (thisrindexfile != "") {
                Sequence temp(inRIndex); m->gobble(inRIndex);
                rindexBarcode.setAligned(temp.getAligned());
            }
            
            int barcodeIndex = 0;
            int primerIndex = 0;
                        
            if(barcodes.size() != 0){
                if (thisfqualfile != "") {
                    if ((thisfindexfile != "") || (thisrindexfile != "")) {
                        success = trimOligos.stripBarcode(findexBarcode, rindexBarcode, *fQual, *rQual, barcodeIndex);
                    }else {
                        success = trimOligos.stripBarcode(fSeq, rSeq, *fQual, *rQual, barcodeIndex);
                    }
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
                for (int i = 0; i < overlapStart; i++) { contig += seq1[i];  }
            }else { //seq1 starts later so take from 0 to overlapStart from seq2
                for (int i = 0; i < overlapStart; i++) {  contig += seq2[i]; }
            }
            
            int seq1End = fSeq.getEndPos();
            int seq2End = rSeq.getEndPos();
            int overlapEnd = seq1End;
            if (seq2End < overlapEnd) { overlapEnd = seq2End; }  //smallest end position is where overlapping ends
            
            int oStart = contig.length();
            //cout << fSeq.getAligned()  << endl; cout << rSeq.getAligned() << endl;
            for (int i = overlapStart; i < overlapEnd; i++) {
                //cout << seq1[i] << ' ' << seq2[i] << ' ' << scores1[ABaseMap[i]] << ' ' << scores2[BBaseMap[i]] << endl;
                if (seq1[i] == seq2[i]) { //match, add base and choose highest score
                    contig += seq1[i];
                }else if (((seq1[i] == '.') || (seq1[i] == '-')) && ((seq2[i] != '-') && (seq2[i] != '.'))) { //seq1 is a gap and seq2 is a base, choose seq2, unless quality score for base is below insert. In that case eliminate base
                    if (thisfqualfile != "") {
                        if (scores2[BBaseMap[i]] <= insert) { } //
                        else { contig += seq2[i];  }
                    }else { contig += seq2[i]; } //with no quality info, then we keep it?
                }else if (((seq2[i] == '.') || (seq2[i] == '-')) && ((seq1[i] != '-') && (seq1[i] != '.'))) { //seq2 is a gap and seq1 is a base, choose seq1, unless quality score for base is below insert. In that case eliminate base
                    if (thisfqualfile != "") {
                        if (scores1[ABaseMap[i]] <= insert) { } //
                        else { contig += seq1[i];  }
                    }else { contig += seq1[i]; } //with no quality info, then we keep it?
                }else if (((seq1[i] != '-') && (seq1[i] != '.')) && ((seq2[i] != '-') && (seq2[i] != '.'))) { //both bases choose one with better quality
                    if (thisfqualfile != "") {
                        if (abs(scores1[ABaseMap[i]] - scores2[BBaseMap[i]]) >= deltaq) { //is the difference in qual scores >= deltaq, if yes choose base with higher score
                            char c = seq1[i];
                            if (scores1[ABaseMap[i]] < scores2[BBaseMap[i]]) { c = seq2[i]; }
                            contig += c;
                        }else { //if no, base becomes n
                            contig += 'N';
                        }
                        numMismatches++;
                    }else { numMismatches++; } //cant decide, so eliminate and mark as mismatch
                }else { //should never get here
                    m->mothurOut("[ERROR]: case I didn't think of seq1 = " + toString(seq1[i]) + " and seq2 = " + toString(seq2[i]) + "\n");
                }
            }
            int oend = contig.length();
            if (seq1End < seq2End) { //seq1 ends before seq2 so take from overlap to length from seq2
                for (int i = overlapEnd; i < length; i++) { contig += seq2[i];  }
            }else { //seq2 ends before seq1 so take from overlap to length from seq1
                for (int i = overlapEnd; i < length; i++) {  contig += seq1[i]; }
            }
            //cout << contig << endl;
            //exit(1);
            if (trimOverlap) { contig = contig.substr(overlapStart-1, oend-oStart);  if (contig.length() == 0) { trashCode += "l"; } }
            
            if(trashCode.length() == 0){
                bool ignore = false;
                
                if (m->debug) { m->mothurOut(fSeq.getName()); }
                
                if (createOligosGroup) {
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
                }else if (createFileGroup) {
                    int pos = group.find("ignore");
                    if (pos == string::npos) {
                        groupMap[fSeq.getName()] = group;
                        
                        map<string, int>::iterator it = groupCounts.find(group);
                        if (it == groupCounts.end()) {	groupCounts[group] = 1; }
                        else { groupCounts[it->first] ++; }
                    }else { ignore = true; }
                }
                if (m->debug) { m->mothurOut("\n"); }
                
                if(allFiles && !ignore){
                    ofstream output;
                    m->openOutputFileAppend(fastaFileNames[barcodeIndex][primerIndex], output);
                    output << ">" << fSeq.getName() << endl << contig << endl;
                    output.close();
                }
                
                //output
                outFasta << ">" << fSeq.getName() << endl << contig << endl;
                int numNs = 0;
                for (int i = 0; i < contig.length(); i++) { if (contig[i] == 'N') { numNs++; }  }
                outMisMatch << fSeq.getName() << '\t' << contig.length() << '\t' << (oend-oStart) << '\t' << oStart << '\t' << oend << '\t' << numMismatches << '\t' << numNs << endl;
            }else {
                //output
                outScrapFasta << ">" << fSeq.getName() << " | " << trashCode << endl << contig << endl;
            }
            num++;
            
			//report progress
			if((num) % 1000 == 0){	m->mothurOutJustToScreen(toString(num)); m->mothurOutEndLine();		}
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
        }
        delete alignment;
        
        if (m->control_pressed) {  m->mothurRemove(outputFasta); m->mothurRemove(outputScrapFasta);m->mothurRemove(outputMisMatches);  }
    
        return num;
    }
	catch(exception& e) {
		m->errorOut(e, "MakeContigsCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
vector< vector<string> > MakeContigsCommand::readFastqFiles(unsigned long int& count, string ffastq, string rfastq, string findex, string rindex){
    try {
        vector< vector<string> > files;
        //maps processors number to file pointer
        map<int, vector<ofstream*> > tempfiles;  //tempfiles[0] = forwardFasta, [1] = forwardQual, [2] = reverseFasta, [3] = reverseQual, tempfiles[4] = forwardIndex, [4] = forwardReverse
        map<int, vector<ofstream*> >::iterator it;
        
        //create files to write to
        for (int i = 0; i < processors; i++) {
            vector<ofstream*> temp;
            ofstream* outFF = new ofstream;     temp.push_back(outFF);
            ofstream* outFQ = new ofstream;     temp.push_back(outFQ);
            ofstream* outRF = new ofstream;     temp.push_back(outRF);
            ofstream* outRQ = new ofstream;     temp.push_back(outRQ);
            ofstream* outFI = new ofstream;     temp.push_back(outFI);  
            ofstream* outRI = new ofstream;     temp.push_back(outRI);  
            tempfiles[i] = temp;
            
            vector<string> names;
            string thisOutputDir = outputDir;
            if (outputDir == "") { thisOutputDir = m->hasPath(ffastq); }
            string ffastafilename = thisOutputDir + m->getRootName(m->getSimpleName(ffastq)) + toString(i) + "ffastatemp";
            string rfastafilename = thisOutputDir + m->getRootName(m->getSimpleName(rfastq)) + toString(i) + "rfastatemp";
            string fqualfilename = thisOutputDir + m->getRootName(m->getSimpleName(ffastq)) + toString(i) + "fqualtemp";
            string rqualfilename = thisOutputDir + m->getRootName(m->getSimpleName(rfastq)) + toString(i) + "rqualtemp";
            string findexfilename = ""; string rindexfilename = "";
            noneOk = false; //flag to oligos file read that its okay to allow for non paired barcodes
            if (findex != "") {  findexfilename = thisOutputDir + m->getRootName(m->getSimpleName(findex)) + toString(i) + "findextemp"; m->openOutputFile(findexfilename, *outFI);  noneOk = true; }
            if (rindex != "") {  rindexfilename = thisOutputDir + m->getRootName(m->getSimpleName(rindex)) + toString(i) + "rindextemp";  m->openOutputFile(rindexfilename, *outRI); noneOk = true; }
            names.push_back(ffastafilename); names.push_back(fqualfilename);
            names.push_back(rfastafilename); names.push_back(rqualfilename);
            names.push_back(findexfilename); names.push_back(rindexfilename);
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
                for(int j = 0; j < files[i].size(); j++) { if (files[i][j] != "") {  m->mothurRemove(files[i][j]); }  }
            }
        }
        
        ifstream inForward;
        m->openInputFile(ffastq, inForward);
        
        ifstream inReverse;
        m->openInputFile(rfastq, inReverse);
        
        ifstream infIndex, inrIndex;
        bool findexIsGood = false;
        bool rindexIsGood = false;
        if (findex != "") { m->openInputFile(findex, infIndex); findexIsGood = true; }
        if (rindex != "") { m->openInputFile(rindex, inrIndex); rindexIsGood = true; }
        
        count = 0;
        map<string, fastqRead> uniques;
        map<string, fastqRead> iUniques;
        map<string, pairFastqRead> pairUniques;
        map<string, fastqRead>::iterator itUniques;
        while ((!inForward.eof()) || (!inReverse.eof()) || (findexIsGood) || (rindexIsGood)) {
            
            if (m->control_pressed) { for (it = tempfiles.begin(); it!=tempfiles.end(); it++) { for (int i = 0; i < (it->second).size(); i++) { (*(it->second)[i]).close();  delete (it->second)[i]; } } for (int i = 0; i < files.size(); i++) {  for(int j = 0; j < files[i].size(); j++) { if (files[i][j] != "") { m->mothurRemove(files[i][j]); } } } inForward.close(); inReverse.close(); if (findex != "") { infIndex.close(); } if (findex != "") { inrIndex.close(); } return files; }
            
            //get a read from forward and reverse fastq files
            bool ignoref, ignorer, ignorefi, ignoreri;
            fastqRead thisFread, thisRread, thisFIread, thisRIread;
            if (!inForward.eof()) {  thisFread = readFastq(inForward, ignoref); }
            else { ignoref = true; }
            if (!inReverse.eof()) { thisRread = readFastq(inReverse, ignorer);  }
            else { ignorer = true; }
            if (findexIsGood) { thisFIread = readFastq(infIndex, ignorefi); if (infIndex.eof()) { findexIsGood = false; } }
            else { ignorefi = true; }
            if (rindexIsGood) { thisRIread = readFastq(inrIndex, ignoreri);  if (inrIndex.eof()) { rindexIsGood = false; } }
            else { ignoreri = true; }
            
            bool allowOne = false;
            if ((findex == "") || (rindex == "")) { allowOne = true; }
            vector<pairFastqRead> frReads = getReads(ignoref, ignorer, thisFread, thisRread, uniques, false);
            vector<pairFastqRead> friReads = getReads(ignorefi, ignoreri, thisFIread, thisRIread, iUniques, allowOne);
            
            //add in index info if provided
            vector<pairFastqRead> reads = frReads;
            if ((findex != "") || (rindex != "")) {  reads = mergeReads(frReads, friReads, pairUniques); }
            
            for (int i = 0; i < reads.size(); i++) {
                fastqRead fread = reads[i].forward;
                fastqRead rread = reads[i].reverse;
                fastqRead firead = reads[i].findex;
                fastqRead riread = reads[i].rindex;
                
                if (m->debug) { m->mothurOut(toString(count) + '\t' + fread.name + '\t' + rread.name + '\n'); if (findex != "") { m->mothurOut(toString(count) + '\t' + firead.name + '\n'); } if (rindex != "") { m->mothurOut(toString(count) + '\t' + riread.name + '\n'); } }
               
                //if (checkReads(fread, rread, ffastq, rfastq)) {
                    if (m->control_pressed) { for (it = tempfiles.begin(); it!=tempfiles.end(); it++) { for (int i = 0; i < (it->second).size(); i++) { (*(it->second)[i]).close();  delete (it->second)[i]; } } for (int i = 0; i < files.size(); i++) {  for(int j = 0; j < files[i].size(); j++) { if (files[i][j] != "") { m->mothurRemove(files[i][j]); } } } inForward.close(); inReverse.close(); if (findex != "") { infIndex.close(); } if (findex != "") { inrIndex.close(); } return files; }
                    
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
                    if (findex != "") {  *(tempfiles[process][4]) << ">" << firead.name << endl << firead.sequence << endl; }
                    if (rindex != "") {  *(tempfiles[process][5]) << ">" << riread.name << endl << riread.sequence << endl; }
                
                    count++;
                    
                    //report progress
                    if((count) % 10000 == 0){	m->mothurOutJustToScreen(toString(count)); m->mothurOutEndLine();		}
                //}
            }
		}
		//report progress
		if((count) % 10000 != 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}
        
        if (uniques.size() != 0) {
            for (itUniques = uniques.begin(); itUniques != uniques.end(); itUniques++) {
                if (m->control_pressed) { break; }
                m->mothurOut("[WARNING]: did not find paired read for " + itUniques->first + ", ignoring.\n");
            }
            for (map<string, pairFastqRead>:: iterator it = pairUniques.begin(); it != pairUniques.end(); it++) {
                if (m->control_pressed) { break; }
                m->mothurOut("[WARNING]: did not find paired read for " + (it->first).substr(1) + ", ignoring.\n");
            }
            m->mothurOutEndLine();
        }
        
        //close files, delete ofstreams
        for (it = tempfiles.begin(); it!=tempfiles.end(); it++) { for (int i = 0; i < (it->second).size(); i++) { (*(it->second)[i]).close();  delete (it->second)[i]; } }
        inForward.close();
        inReverse.close();
        if (findex != "") { infIndex.close(); }
        if (rindex != "") { inrIndex.close(); }
       
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
            string findexfilename = ""; string rindexfilename = "";
            names.push_back(ffastafilename); names.push_back(fqualfilename);
            names.push_back(rfastafilename); names.push_back(rqualfilename);
            names.push_back(findexfilename); names.push_back(rindexfilename);
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
            
            vector<pairFastqRead> reads = getReads(ignoref, ignorer, thisFread, thisRread, uniques, false);
            
            for (int i = 0; i < reads.size(); i++) {
                fastqRead fread = reads[i].forward;
                fastqRead rread = reads[i].reverse;
                
                if (m->debug) { m->mothurOut(toString(count) + '\t' + fread.name + '\t' + rread.name + '\n'); }
                
               // if (checkReads(fread, rread, ffasta, rfasta)) {
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
                //}
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
vector<pairFastqRead> MakeContigsCommand::getReads(bool ignoref, bool ignorer, fastqRead forward, fastqRead reverse, map<string, fastqRead>& uniques, bool allowOne){
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
            if (allowOne) {
                fastqRead dummy;
                pairFastqRead temp(forward, dummy);
                reads.push_back(temp);
            }else {
                //look for forward pair
                itUniques = uniques.find(forward.name);
                if (itUniques != uniques.end()) {  //we have the pair for this read
                    pairFastqRead temp(forward, itUniques->second);
                    reads.push_back(temp);
                    uniques.erase(itUniques);
                }else { //save this read for later
                    uniques[forward.name] = forward;
                }
            }
        }else if (ignoref && !ignorer) { //ignore forward keep reverse
            if (allowOne) {
                fastqRead dummy;
                pairFastqRead temp(dummy, reverse);
                reads.push_back(temp);
            }else {
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
        }//else ignore both and do nothing
        
        return reads;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "readFastqFiles");
        exit(1);
    }
}
//**********************************************************************************************************************
//look through the reads from the forward and reverse files and try to find matching reads from index files.
vector<pairFastqRead> MakeContigsCommand::mergeReads(vector<pairFastqRead> thisReads, vector<pairFastqRead> indexes, map<string, pairFastqRead>& uniques){
    try {
        vector<pairFastqRead> reads;
        map<string, pairFastqRead>::iterator itUniques;
        
        set<int> foundIndexes;
        for (int i = 0; i < thisReads.size(); i++) {
            bool found = false;
            for (int j = 0; j < indexes.size(); j++) {
                
                //incase only one index
                string indexName = indexes[j].forward.name;
                if (indexName == "") { indexName = indexes[j].reverse.name; }
                
                if (thisReads[i].forward.name == indexName){
                    thisReads[i].findex = indexes[j].forward;
                    thisReads[i].rindex = indexes[j].reverse;
                    reads.push_back(thisReads[i]);
                    found = true;
                    foundIndexes.insert(j);
                }
            }
            
            if (!found) {
                //look for forward pair
                itUniques = uniques.find('i'+thisReads[i].forward.name);
                if (itUniques != uniques.end()) {  //we have the pair for this read
                    thisReads[i].findex = itUniques->second.forward;
                    thisReads[i].rindex = itUniques->second.reverse;
                    reads.push_back(thisReads[i]);
                    uniques.erase(itUniques);
                }else { //save this read for later
                    uniques['r'+thisReads[i].forward.name] = thisReads[i];
                }
            }
        }
        
        if (foundIndexes.size() != indexes.size()) { //if we didnt match all the indexes look for them in uniques
            for (int j = 0; j < indexes.size(); j++) {
                if (foundIndexes.count(j) == 0) { //we didnt find this one
                    //incase only one index
                    string indexName = indexes[j].forward.name;
                    if (indexName == "") { indexName = indexes[j].reverse.name; }
                    
                    //look for forward pair
                    itUniques = uniques.find('r'+indexName);
                    if (itUniques != uniques.end()) {  //we have the pair for this read
                        pairFastqRead temp(itUniques->second.forward, itUniques->second.reverse, indexes[j].forward, indexes[j].reverse);
                        reads.push_back(temp);
                        uniques.erase(itUniques);
                    }else { //save this read for later
                        uniques['i'+indexName] = indexes[j];
                    }
                }
            }
        }
       
        
        return reads;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "mergeReads");
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
        
        vector<int> qualScores = convertQual(quality);
        
        m->checkName(name);
        read.name = name;
        read.sequence = sequence;
        read.scores = qualScores;
        
        if (m->debug) { m->mothurOut("[DEBUG]: " + read.name + " " + read.sequence + " " + quality + "\n"); }

        return read;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "readFastq");
        exit(1);
    }
}
/**********************************************************************************************************************
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
}*/
//***************************************************************************************************************
//lines can be 2, 3, or 4 columns
// forward.fastq reverse.fastq -> 2 column
// groupName forward.fastq reverse.fastq -> 3 column
// forward.fastq reverse.fastq forward.index.fastq  reverse.index.fastq  -> 4 column
// forward.fastq reverse.fastq none  reverse.index.fastq  -> 4 column
// forward.fastq reverse.fastq forward.index.fastq  none  -> 4 column
vector< vector<string> > MakeContigsCommand::readFileNames(string filename){
	try {
        vector< vector<string> > files;
        string forward, reverse, findex, rindex;
        
        ifstream in;
        m->openInputFile(filename, in);
        
        while(!in.eof()) {
            
            if (m->control_pressed) { return files; }
            
            string line = m->getline(in);  m->gobble(in);
            vector<string> pieces = m->splitWhiteSpace(line);
            
            string group = "";
            if (pieces.size() == 2) {
                forward = pieces[0];
                reverse = pieces[1];
                group = "";
                findex = "";
                rindex = "";
            }else if (pieces.size() == 3) {
                group = pieces[0];
                forward = pieces[1];
                reverse = pieces[2];
                findex = "";
                rindex = "";
                createFileGroup = true;
            }else if (pieces.size() == 4) {
                forward = pieces[0];
                reverse = pieces[1];
                findex = pieces[2];
                rindex = pieces[3];
                if ((findex == "none") || (findex == "NONE")){ findex = ""; }
                if ((rindex == "none") || (rindex == "NONE")){ rindex = ""; }
            }else {
                m->mothurOut("[ERROR]: file lines can be 2, 3, or 4 columns. The forward fastq files in the first column and their matching reverse fastq files in the second column, or a groupName then forward fastq file and reverse fastq file, or forward fastq file then reverse fastq then forward index and reverse index file.  If you only have one index file add 'none' for the other one. \n"); m->control_pressed = true;
            }
            
            if (m->debug) { m->mothurOut("[DEBUG]: group = " + group + ", forward = " + forward + ", reverse = " + reverse + ", forwardIndex = " + findex + ", reverseIndex = " + rindex + ".\n"); }
            
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
            
            int openFindex = 0;
            if (findex != "") {
                ifstream in4;
                openFindex = m->openInputFile(findex, in4, "noerror"); in4.close();
                
                //if you can't open it, try default location
                if (openFindex == 1) {
                    if (m->getDefaultPath() != "") { //default path is set
                        string tryPath = m->getDefaultPath() + m->getSimpleName(findex);
                        m->mothurOut("Unable to open " + findex + ". Trying default " + tryPath); m->mothurOutEndLine();
                        ifstream in5;
                        openFindex = m->openInputFile(tryPath, in5, "noerror");
                        in5.close();
                        findex = tryPath;
                    }
                }
                
                //if you can't open it, try output location
                if (openFindex == 1) {
                    if (m->getOutputDir() != "") { //default path is set
                        string tryPath = m->getOutputDir() + m->getSimpleName(findex);
                        m->mothurOut("Unable to open " + findex + ". Trying output directory " + tryPath); m->mothurOutEndLine();
                        ifstream in6;
                        openFindex = m->openInputFile(tryPath, in6, "noerror");
                        findex = tryPath;
                        in6.close();
                    }
                }
                
                if (openFindex == 1) { //can't find it
                    m->mothurOut("[WARNING]: can't find " + findex + ", ignoring pair.\n");
                }
            }
            
            int openRindex = 0;
            if (rindex != "") {
                ifstream in7;
                openRindex = m->openInputFile(rindex, in7, "noerror"); in7.close();
                
                //if you can't open it, try default location
                if (openRindex == 1) {
                    if (m->getDefaultPath() != "") { //default path is set
                        string tryPath = m->getDefaultPath() + m->getSimpleName(rindex);
                        m->mothurOut("Unable to open " + rindex + ". Trying default " + tryPath); m->mothurOutEndLine();
                        ifstream in8;
                        openRindex = m->openInputFile(tryPath, in8, "noerror");
                        in8.close();
                        rindex = tryPath;
                    }
                }
                
                //if you can't open it, try output location
                if (openRindex == 1) {
                    if (m->getOutputDir() != "") { //default path is set
                        string tryPath = m->getOutputDir() + m->getSimpleName(rindex);
                        m->mothurOut("Unable to open " + rindex + ". Trying output directory " + tryPath); m->mothurOutEndLine();
                        ifstream in9;
                        openRindex = m->openInputFile(tryPath, in9, "noerror");
                        rindex = tryPath;
                        in9.close();
                    }
                }
                
                if (openRindex == 1) { //can't find it
                    m->mothurOut("[WARNING]: can't find " + rindex + ", ignoring pair.\n");
                }
            }

            
            if ((openForward != 1) && (openReverse != 1) && (openFindex != 1) && (openRindex != 1)) { //good pair
                file2Group[files.size()] = group;
                vector<string> pair;
                pair.push_back(forward);
                pair.push_back(reverse);
                pair.push_back(findex);
                pair.push_back(rindex);
                if (((findex != "") || (rindex != "")) && (oligosfile == "")) { m->mothurOut("[ERROR]: You need to provide an oligos file if you are going to use an index file.\n"); m->control_pressed = true;  }
                files.push_back(pair);
            }
        }
        in.close();
        
        return files;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "readFileNames");
        exit(1);
    }
}
//***************************************************************************************************************
//illumina data requires paired forward and reverse data
//BARCODE   atgcatgc   atgcatgc    groupName 
//PRIMER   atgcatgc   atgcatgc    groupName  
//PRIMER   atgcatgc   atgcatgc  
bool MakeContigsCommand::getOligos(vector<vector<string> >& fastaFileNames, string rootname){
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
				
				if(type == "PRIMER"){
					m->gobble(in);
					
                    in >> roligo;
                    
                    for(int i=0;i<roligo.length();i++){
                        roligo[i] = toupper(roligo[i]);
                        if(roligo[i] == 'U')	{	roligo[i] = 'T';	}
                    }
                    //roligo = reverseOligo(roligo);
                    
                    if (m->debug) { m->mothurOut("[DEBUG]: reading - " + roligo + ".\n"); }
                    
                    group = "";
                    
					// get rest of line in case there is a primer name
					while (!in.eof())	{	
						char c = in.get(); 
						if (c == 10 || c == 13 || c == -1){	break;	}
						else if (c == 32 || c == 9){;} //space or tab
						else { 	group += c;  }
					} 
                    
                    oligosPair newPrimer(foligo, roligo);
                    
                    if (m->debug) { m->mothurOut("[DEBUG]: primer pair " + newPrimer.forward + " " + newPrimer.reverse + ", and group = " + group + ".\n"); }
                    
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
                    
                    if ((foligo == "NONE") || (roligo == "NONE")) { if (!noneOk) { m->mothurOut("[ERROR]: barcodes must be paired unless you are using an index file.\n"); m->control_pressed = true; } }
                    
                    group = "";
                    while (!in.eof())	{	
						char c = in.get(); 
						if (c == 10 || c == 13 || c == -1){	break;	}
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
vector<int> MakeContigsCommand::convertQual(string qual) {
	try {
		vector<int> qualScores;
        bool negativeScores = false;
		
		for (int i = 0; i < qual.length(); i++) { 
            
            int temp = 0;
            temp = int(qual[i]);
            if (format == "illumina") {
                temp -= 64; //char '@'
            }else if (format == "illumina1.8+") {
                    temp -= int('!'); //char '!'
            }else if (format == "solexa") {
                temp = int(convertTable[temp]); //convert to sanger
                temp -= int('!'); //char '!'
            }else {
                temp -= int('!'); //char '!'
            }
            
            if (temp < -5) { negativeScores = true; }
			qualScores.push_back(temp);
		}
		
        if (negativeScores) { m->mothurOut("[ERROR]: finding negative quality scores, do you have the right format selected? http://en.wikipedia.org/wiki/FASTQ_format#Encoding \n");  m->control_pressed = true;  }
        
		return qualScores;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeContigsCommand", "convertQual");
		exit(1);
	}
}

//**********************************************************************************************************************




