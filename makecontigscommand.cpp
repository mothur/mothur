//
//  makecontigscommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/15/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "makecontigscommand.h"

//********************************************************************************************************************
//sorts biggest to smallest
inline bool compareFileSizes(vector<string> left, vector<string> right){
    
    FILE * pFile;
    long leftsize = 0;
    
    //get num bytes in file
    string filename = left[0];
    pFile = fopen (filename.c_str(),"rb");
    string error = "Error opening " + filename;
    if (pFile==NULL) perror (error.c_str());
    else{
        fseek (pFile, 0, SEEK_END);
        leftsize=ftell (pFile);
        fclose (pFile);
    }
    
    FILE * pFile2;
    long rightsize = 0;
    
    //get num bytes in file
    filename = right[0];
    pFile2 = fopen (filename.c_str(),"rb");
    error = "Error opening " + filename;
    if (pFile2==NULL) perror (error.c_str());
    else{
        fseek (pFile2, 0, SEEK_END);
        rightsize=ftell (pFile2);
        fclose (pFile2);
    }
    
    return (leftsize > rightsize);
}
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
        CommandParameter preorient("checkorient", "Boolean", "", "F", "", "", "","",false,false,true); parameters.push_back(preorient);
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
        helpString += "The checkorient parameter will check look for the reverse compliment of the barcode or primer in the sequence. If found the sequence is flipped. The default is false.\n";
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
        else if (type == "qfile") {  pattern = "[filename],[tag],contigs.qual"; }
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
        outputTypes["qfile"] = tempOutNames;
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
            outputTypes["qfile"] = tempOutNames;
            outputTypes["report"] = tempOutNames;
            outputTypes["group"] = tempOutNames;
			
            
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			inputDir = validParameter.validFile(parameters, "inputdir", false);		
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
            
            temp = validParameter.validFile(parameters, "checkorient", false);		if (temp == "not found") { temp = "F"; }
			reorient = m->isTrue(temp);
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
        
        unsigned long long numReads = 0;
        map<string, int> totalGroupCounts;
        int start = time(NULL);
        longestBase = 1000;
        
        if (file != "") {
            numReads = processMultipleFileOption(totalGroupCounts);
        }else if ((ffastqfile != "") || (ffastafile != "")) {
            numReads = processSingleFileOption(totalGroupCounts);
        }else {  return 0; }
        
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
unsigned long long MakeContigsCommand::processSingleFileOption(map<string, int>& totalGroupCounts) {
    try {
        bool hasQual = false;
        unsigned long long numReads = 0;
        string inputFile = "";
        vector<string> fileInputs;
        vector<string> qualOrIndexInputs;
        vector<linePair> lines;
        vector<linePair> qLines;
        delim = '>';
        map<string, string> variables;
        string thisOutputDir = outputDir;
        
        if (ffastafile != "") {
            inputFile = ffastafile;
            if (outputDir == "") {  thisOutputDir = m->hasPath(inputFile); }
            fileInputs.push_back(ffastafile); fileInputs.push_back(rfastafile);
            
            if (fqualfile != "") {
                hasQual = true;
                qualOrIndexInputs.push_back(fqualfile); qualOrIndexInputs.push_back(rqualfile);
                variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(fqualfile));
                variables["[tag]"] = "trim";
                outQualFile = getOutputFileName("qfile",variables);
                variables["[tag]"] = "scrap";
                outScrapQualFile = getOutputFileName("qfile",variables);
            }else {
                outQualFile = ""; outScrapQualFile = "";
            }
            
            variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(inputFile));
            delim = '>';
        }else { //ffastqfile
            hasQual = true;
            inputFile = ffastqfile;
            if (outputDir == "") {  thisOutputDir = m->hasPath(inputFile); }
            variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(inputFile));
            variables["[tag]"] = "trim";
            outQualFile = getOutputFileName("qfile",variables);
            variables["[tag]"] = "scrap";
            outScrapQualFile = getOutputFileName("qfile",variables);
            
            fileInputs.push_back(ffastqfile); fileInputs.push_back(rfastqfile);
            if ((findexfile != "") || (rindexfile != "")){
                qualOrIndexInputs.push_back("NONE"); qualOrIndexInputs.push_back("NONE");
                if (findexfile != "") { qualOrIndexInputs[0] = findexfile; }
                if (rindexfile != "") { qualOrIndexInputs[1] = rindexfile; }
            }
            delim = '@';
        }
        variables["[tag]"] = "trim";
        outFastaFile = getOutputFileName("fasta",variables);
        variables["[tag]"] = "scrap";
        outScrapFastaFile = getOutputFileName("fasta",variables);
        variables["[tag]"] = "";
        outMisMatchFile = getOutputFileName("report",variables);
        
        //divides the files so that the processors can share the workload.
        setLines(fileInputs, qualOrIndexInputs, lines, qLines, delim);
        
        vector<vector<string> > fastaFileNames, qualFileNames;
        map<string, string> uniqueFastaNames;// so we don't add the same groupfile multiple times
        createOligosGroup = false;
        oligos = new Oligos();
        numBarcodes = 0; numFPrimers= 0; numLinkers= 0; numSpacers = 0; numRPrimers = 0;
        
        if(oligosfile != "")                        {       createOligosGroup = getOligos(fastaFileNames, qualFileNames, variables["[filename]"], uniqueFastaNames);    }
        if (createOligosGroup || createFileGroup)   {       outputGroupFileName = getOutputFileName("group",variables);                                                 }
        
        //give group in file file precedence
        if (createFileGroup) {  createOligosGroup = false; }
        
        m->mothurOut("Making contigs...\n");
        numReads = createProcesses(fileInputs, qualOrIndexInputs, outFastaFile, outScrapFastaFile, outQualFile, outScrapQualFile, outMisMatchFile, fastaFileNames, qualFileNames, lines, qLines, "");
        
        if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); }  delete oligos; return 0; }
        
        if(allFiles){
            // so we don't add the same groupfile multiple times
            map<string, string>::iterator it;
            set<string> namesToRemove;
            for(int i=0;i<fastaFileNames.size();i++){
                for(int j=0;j<fastaFileNames[0].size();j++){
                    if (fastaFileNames[i][j] != "") {
                        if (namesToRemove.count(fastaFileNames[i][j]) == 0) {
                            if(m->isBlank(fastaFileNames[i][j])){
                                m->mothurRemove(fastaFileNames[i][j]);
                                namesToRemove.insert(fastaFileNames[i][j]);
                                uniqueFastaNames.erase(fastaFileNames[i][j]); //remove from list for group file print
                                
                                m->mothurRemove(qualFileNames[i][j]);
                                namesToRemove.insert(qualFileNames[i][j]);
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
                variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(it->first));
                string thisGroupName = getOutputFileName("group",variables); outputNames.push_back(thisGroupName); outputTypes["group"].push_back(thisGroupName);
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
        
        if (file == "") {
            totalGroupCounts = groupCounts;
            outputNames.push_back(outFastaFile); outputTypes["fasta"].push_back(outFastaFile);
            outputNames.push_back(outScrapFastaFile); outputTypes["fasta"].push_back(outScrapFastaFile);
            if (hasQual) {
                outputNames.push_back(outQualFile); outputTypes["qfile"].push_back(outQualFile);
                outputNames.push_back(outScrapQualFile); outputTypes["qfile"].push_back(outScrapQualFile);
            }
            outputNames.push_back(outMisMatchFile); outputTypes["report"].push_back(outMisMatchFile);
            if (createFileGroup || createOligosGroup) {
                outputNames.push_back(outputGroupFileName); outputTypes["group"].push_back(outputGroupFileName);
            }
        }
        m->mothurOut("Done.\n");
        delete oligos;
        
        return numReads;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "processSingleFileOption");
        exit(1);
    }
}
//**********************************************************************************************************************
unsigned long long MakeContigsCommand::processMultipleFileOption(map<string, int>& totalGroupCounts) {
    try {
        unsigned long long numReads = 0;
        
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
        string compositeQualFile = getOutputFileName("qfile",cvars);
        cvars["[tag]"] = "scrap";
        string compositeScrapQualFile = getOutputFileName("qfile",cvars);
        cvars["[tag]"] = "";
        string compositeMisMatchFile = getOutputFileName("report",cvars);
        
        ofstream outCTFasta, outCTQual, outCSFasta, outCSQual, outCMisMatch;
        m->openOutputFile(compositeFastaFile, outCTFasta); outCTFasta.close();
        m->openOutputFile(compositeQualFile, outCTQual); outCTQual.close();
        m->openOutputFile(compositeScrapFastaFile, outCSFasta); outCSFasta.close();
        m->openOutputFile(compositeScrapQualFile, outCSQual); outCSQual.close();
        m->openOutputFile(compositeMisMatchFile, outCMisMatch); outCMisMatch.close();
        outputNames.push_back(compositeFastaFile); outputTypes["fasta"].push_back(compositeFastaFile);
        outputNames.push_back(compositeQualFile); outputTypes["qfile"].push_back(compositeQualFile);
        outputNames.push_back(compositeMisMatchFile); outputTypes["report"].push_back(compositeMisMatchFile);
        outputNames.push_back(compositeScrapFastaFile); outputTypes["fasta"].push_back(compositeScrapFastaFile);
        outputNames.push_back(compositeScrapQualFile); outputTypes["qfile"].push_back(compositeScrapQualFile);
        
        //read file
        vector< vector<string> > fileInputs = readFileNames(file);
        
        for (int l = 0; l < fileInputs.size(); l++) {
            
            m->mothurOut("\n>>>>>\tProcessing file pair " + fileInputs[l][0] + " - " + fileInputs[l][1] + " (files " + toString(l+1) + " of " + toString(fileInputs.size()) + ")\t<<<<<\n");
            
            ffastqfile = fileInputs[l][0];
            rfastqfile = fileInputs[l][1];
            findexfile = fileInputs[l][2];
            rindexfile = fileInputs[l][3];
            string group = file2Group[l];
            groupCounts.clear();
            groupMap.clear();
            
            //run file as if it was a single
            numReads += processSingleFileOption(groupCounts);
            
            //append to combo files
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
            m->appendFiles(outQualFile, compositeQualFile);
            m->appendFiles(outScrapQualFile, compositeScrapQualFile);
            if (!allFiles) {
                m->mothurRemove(outMisMatchFile);
                m->mothurRemove(outFastaFile);
                m->mothurRemove(outScrapFastaFile);
                m->mothurRemove(outQualFile);
                m->mothurRemove(outScrapQualFile);
            }else {
                outputNames.push_back(outFastaFile); outputTypes["fasta"].push_back(outFastaFile);
                outputNames.push_back(outScrapFastaFile); outputTypes["fasta"].push_back(outScrapFastaFile);
                outputNames.push_back(outQualFile); outputTypes["qfile"].push_back(outQualFile);
                outputNames.push_back(outScrapQualFile); outputTypes["qfile"].push_back(outScrapQualFile);
                outputNames.push_back(outMisMatchFile); outputTypes["report"].push_back(outMisMatchFile);
            }
        }
        
         return numReads;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "processMultipleFileOption");
        exit(1);
    }
}
//**********************************************************************************************************************
//fileInputs[0] = forward Fasta or Forward Fastq, fileInputs[1] = reverse Fasta or reverse Fastq. if qualOrIndexFiles.size() != 0, then qualOrIndexFiles[0] = forward qual or Forward index, qualOrIndexFiles[1] = reverse qual or reverse index.
//lines[0] - ffasta, lines[1] - rfasta) - processor1
//lines[2] - ffasta, lines[3] - rfasta) - processor2
//lines[4] - ffasta, lines[5] - rfasta) - processor3
//...
//qlines[0] - fqual or findex, qlines[1] - rqual or rindex) - processor1
//qlines[2] - fqual or findex, qlines[3] - rqual or rindex) - processor2
//qlines[4] - fqual or findex, qlines[5] - rqual or rindex) - processor3
//...
//if using index files and only have 1 then the file name = NONE, and entries are duds. Copies of other index file.
//if no index files are given, then qualOrIndexFiles.size() == 0.
unsigned long long MakeContigsCommand::createProcesses(vector<string> fileInputs, vector<string> qualOrIndexFiles, string outputFasta, string outputScrapFasta, string outputQual, string outputScrapQual, string outputMisMatches, vector<vector<string> > fastaFileNames, vector<vector<string> > qualFileNames, vector<linePair> lines, vector<linePair> qLines, string group) {
	try {
		int num = 0;
		vector<int> processIDS;
        
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		int process = 1;
		
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
                vector<vector<string> > tempFASTAFileNames = fastaFileNames;
                vector<vector<string> > tempQUALFileNames = qualFileNames;
                
				if(allFiles){
					ofstream temp;
                    
					for(int i=0;i<tempFASTAFileNames.size();i++){
						for(int j=0;j<tempFASTAFileNames[i].size();j++){
							if (tempFASTAFileNames[i][j] != "") {
								tempFASTAFileNames[i][j] += m->mothurGetpid(process) + ".temp";
								m->openOutputFile(tempFASTAFileNames[i][j], temp);			temp.close();
							}
                            if (tempQUALFileNames[i][j] != "") {
                                tempQUALFileNames[i][j] += m->mothurGetpid(process) + ".temp";
                                m->openOutputFile(tempQUALFileNames[i][j], temp);			temp.close();
                            }
						}
					}
				}
                
                int spot = process*2;
				num = driver(fileInputs, qualOrIndexFiles,
                             outputFasta + m->mothurGetpid(process) + ".temp",
                             outputScrapFasta + m->mothurGetpid(process) + ".temp",
                             outputQual + m->mothurGetpid(process) + ".temp",
                             outputScrapQual + m->mothurGetpid(process) + ".temp",
                             outputMisMatches + m->mothurGetpid(process) + ".temp",
                             tempFASTAFileNames, tempQUALFileNames, lines[spot], lines[spot+1], qLines[spot], qLines[spot+1], group);
				
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
        m->openOutputFile(outputQual, temp);		temp.close();
        m->openOutputFile(outputScrapQual, temp);		temp.close();
    
		//do my part
        int spot = 0;
		num = driver(fileInputs, qualOrIndexFiles, outputFasta, outputScrapFasta,  outputQual, outputScrapQual, outputMisMatches, fastaFileNames, qualFileNames, lines[spot], lines[spot+1], qLines[spot], qLines[spot+1], group);
		
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
            vector<vector<string> > tempQUALFileNames = qualFileNames;
                        
            if(allFiles){
                ofstream temp;
                
                for(int i=0;i<tempFASTAFileNames.size();i++){
                    for(int j=0;j<tempFASTAFileNames[i].size();j++){
                        if (tempFASTAFileNames[i][j] != "") {
                            tempFASTAFileNames[i][j] += extension;
                            m->openOutputFile(tempFASTAFileNames[i][j], temp);			temp.close();
                        }
                        if (tempQUALFileNames[i][j] != "") {
                            tempQUALFileNames[i][j] += extension;
                            m->openOutputFile(tempQUALFileNames[i][j], temp);			temp.close();
                        }
                    }
                }
            }
			
            int spot = (h)*2;
			contigsData* tempcontig = new contigsData(format, delim, group, fileInputs, qualOrIndexFiles, (outputFasta + extension), (outputScrapFasta + extension), (outputQual + extension), (outputScrapQual + extension), (outputMisMatches + extension), align, m, match, misMatch, gapOpen, gapExtend, insert, deltaq, tempFASTAFileNames, tempQUALFileNames, oligosfile, reorient, pdiffs, bdiffs, tdiffs, createOligosGroup, createFileGroup, allFiles, trimOverlap, lines[spot], lines[spot+1], qLines[spot], qLines[spot+1], h);
			pDataArray.push_back(tempcontig);
            
			hThreadArray[h] = CreateThread(NULL, 0, MyContigsThreadFunction, pDataArray[h], 0, &dwThreadIdArray[h]);   
		}
        
        vector<vector<string> > tempFASTAFileNames = fastaFileNames;
        vector<vector<string> > tempQUALFileNames = qualFileNames;

        if(allFiles){
            ofstream temp;
            string extension = toString(processors-1) + ".temp";
            
            for(int i=0;i<tempFASTAFileNames.size();i++){
                for(int j=0;j<tempFASTAFileNames[i].size();j++){
                    if (tempFASTAFileNames[i][j] != "") {
                        tempFASTAFileNames[i][j] += extension;
                        m->openOutputFile(tempFASTAFileNames[i][j], temp);			temp.close();
                    }
                    if (tempQUALFileNames[i][j] != "") {
                        tempQUALFileNames[i][j] += extension;
                        m->openOutputFile(tempQUALFileNames[i][j], temp);			temp.close();
                    }
                }
            }
        }

		//parent do my part
		ofstream temp, temp2, temp3, temp4;
		m->openOutputFile(outputFasta, temp);		temp.close();
        m->openOutputFile(outputScrapFasta, temp2);		temp2.close();
        m->openOutputFile(outputQual, temp3);		temp3.close();
        m->openOutputFile(outputScrapQual, temp4);		temp4.close();
        
        //do my part
        int spot = (processors-1)*2;
        processIDS.push_back(processors-1);
        num = driver(fileInputs, qualOrIndexFiles, (outputFasta+ toString(processors-1) + ".temp"),  (outputScrapFasta+ toString(processors-1) + ".temp"),  (outputQual+ toString(processors-1) + ".temp"),  (outputScrapQual+ toString(processors-1) + ".temp"), (outputMisMatches+ toString(processors-1) + ".temp"), tempFASTAFileNames, tempQUALFileNames, lines[spot], lines[spot+1], qLines[spot], qLines[spot+1], group);

        
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
            
            m->appendFiles((outputQual + toString(processIDS[i]) + ".temp"), outputQual);
            m->mothurRemove((outputQual + toString(processIDS[i]) + ".temp"));
            
            m->appendFiles((outputScrapQual + toString(processIDS[i]) + ".temp"), outputScrapQual);
            m->mothurRemove((outputScrapQual + toString(processIDS[i]) + ".temp"));
            
            m->appendFilesWithoutHeaders((outputMisMatches + toString(processIDS[i]) + ".temp"), outputMisMatches);
			m->mothurRemove((outputMisMatches + toString(processIDS[i]) + ".temp"));
            
            if(allFiles){
				for(int j=0;j<fastaFileNames.size();j++){
					for(int k=0;k<fastaFileNames[j].size();k++){
						if (fastaFileNames[j][k] != "") {
							m->appendFiles((fastaFileNames[j][k] + toString(processIDS[i]) + ".temp"), fastaFileNames[j][k]);
							m->mothurRemove((fastaFileNames[j][k] + toString(processIDS[i]) + ".temp"));
						}
                        if (qualFileNames[j][k] != "") {
                            m->appendFiles((qualFileNames[j][k] + toString(processIDS[i]) + ".temp"), qualFileNames[j][k]);
                            m->mothurRemove((qualFileNames[j][k] + toString(processIDS[i]) + ".temp"));
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
int MakeContigsCommand::driver(vector<string> inputFiles, vector<string> qualOrIndexFiles, string outputFasta, string outputScrapFasta, string outputQual, string outputScrapQual,  string outputMisMatches, vector<vector<string> > fastaFileNames, vector<vector<string> > qualFileNames, linePair linesInput, linePair linesInputReverse, linePair qlinesInput, linePair qlinesInputReverse, string group){
    try {
        
        Alignment* alignment;
        if(align == "gotoh")			{	alignment = new GotohOverlap(gapOpen, gapExtend, match, misMatch, longestBase);			}
		else if(align == "needleman")	{	alignment = new NeedlemanOverlap(gapOpen, match, misMatch, longestBase);				}
        
        int num = 0;
        string thisfqualindexfile, thisrqualindexfile, thisffastafile, thisrfastafile;
        thisfqualindexfile = ""; thisrqualindexfile = "";
        thisffastafile = inputFiles[0]; thisrfastafile = inputFiles[1];
        if (qualOrIndexFiles.size() != 0) {
            thisfqualindexfile = qualOrIndexFiles[0];
            thisrqualindexfile = qualOrIndexFiles[1];
        }
        
        if (m->debug) {  m->mothurOut("[DEBUG]: ffasta = " + thisffastafile + ".\n[DEBUG]: rfasta = " + thisrfastafile + ".\n[DEBUG]: fqualindex = " + thisfqualindexfile + ".\n[DEBUG]: rqualindex = " + thisfqualindexfile + ".\n"); }
        
        ifstream inFFasta, inRFasta, inFQualIndex, inRQualIndex;
        ofstream outFasta, outMisMatch, outScrapFasta, outQual, outScrapQual;
        m->openInputFile(thisffastafile, inFFasta);
        m->openInputFile(thisrfastafile, inRFasta);
        
        inFFasta.seekg(linesInput.start);
        inRFasta.seekg(linesInputReverse.start);
        
        if (thisfqualindexfile != "") {
            if (thisfqualindexfile != "NONE") {
                m->openInputFile(thisfqualindexfile, inFQualIndex);
                inFQualIndex.seekg(qlinesInput.start);
            }
            else {  thisfqualindexfile = ""; }
            if (thisrqualindexfile != "NONE") {
                m->openInputFile(thisrqualindexfile, inRQualIndex);
                inRQualIndex.seekg(qlinesInputReverse.start);
            }
            else { thisrqualindexfile = ""; }
        }
        
        m->openOutputFile(outputFasta, outFasta);
        m->openOutputFile(outputScrapFasta, outScrapFasta);
        m->openOutputFile(outputMisMatches, outMisMatch);
        bool hasQuality = false;
        outMisMatch << "Name\tLength\tOverlap_Length\tOverlap_Start\tOverlap_End\tMisMatches\tNum_Ns\n";
        if (delim == '@') { //fastq files so make an output quality
            m->openOutputFile(outputQual, outQual);
            m->openOutputFile(outputScrapQual, outScrapQual);
            hasQuality = true;
        }else if ((delim == '>') && (qualOrIndexFiles.size() != 0)) { //fasta and qual files
            m->openOutputFile(outputQual, outQual);
            m->openOutputFile(outputScrapQual, outScrapQual);
            hasQuality = true;
        }
        
        TrimOligos trimOligos(pdiffs, bdiffs, 0, 0, oligos->getPairedPrimers(), oligos->getPairedBarcodes());
        
        TrimOligos* rtrimOligos = NULL;
        if (reorient) {  rtrimOligos = new TrimOligos(pdiffs, bdiffs, 0, 0, oligos->getReorientedPairedPrimers(), oligos->getReorientedPairedBarcodes()); numBarcodes = oligos->getReorientedPairedBarcodes().size();    }
        
        while ((!inFFasta.eof()) && (!inRFasta.eof())) {
            
            if (m->control_pressed) { break; }
            
            int success = 1;
            string trashCode = "";
            string commentString = "";
            int currentSeqsDiffs = 0;
            bool hasIndex = false;
            
            bool ignore; ignore = false;
            Sequence fSeq, rSeq;
            QualityScores* fQual = NULL; QualityScores* rQual = NULL;
            QualityScores* savedFQual = NULL; QualityScores* savedRQual = NULL;
            Sequence findexBarcode("findex", "NONE");  Sequence rindexBarcode("rindex", "NONE");
            if (delim == '@') { //fastq files
                bool tignore;
                FastqRead fread(inFFasta, tignore, format); m->gobble(inFFasta);
                FastqRead rread(inRFasta, ignore, format); m->gobble(inRFasta);
                if (tignore) { ignore=true; }
                fSeq.setName(fread.getName()); fSeq.setAligned(fread.getSeq());
                rSeq.setName(rread.getName()); rSeq.setAligned(rread.getSeq());
                fQual = new QualityScores(fread.getName(), fread.getScores());
                rQual = new QualityScores(rread.getName(), rread.getScores());
                if (thisfqualindexfile != "") { //forward index file
                    FastqRead firead(inFQualIndex, tignore, format); m->gobble(inFQualIndex);
                    if (tignore) { ignore=true; }
                    findexBarcode.setAligned(firead.getSeq());
                    if (firead.getName() != fread.getName()) { m->mothurOut("[WARNING]: name mismatch in forward index file. Ignoring, " + fread.getName() + ".\n"); ignore = true; }
                    hasIndex = true;
                }
                if (thisrqualindexfile != "") { //reverse index file
                    FastqRead riread(inRQualIndex, tignore, format); m->gobble(inRQualIndex);
                    if (tignore) { ignore=true; }
                    rindexBarcode.setAligned(riread.getSeq());
                    if (riread.getName() != fread.getName()) { m->mothurOut("[WARNING]: name mismatch in reverse index file. Ignoring, " + fread.getName() + ".\n"); ignore = true; }
                    hasIndex = true;
                }
                if (fread.getName() != rread.getName()) { m->mothurOut("[WARNING]: name mismatch in forward and reverse fastq file. Ignoring, " + fread.getName() + ".\n"); ignore = true; }
            }else { //reading fasta and maybe qual
                Sequence tfSeq(inFFasta); m->gobble(inFFasta);
                Sequence trSeq(inRFasta); m->gobble(inRFasta);
                fSeq.setName(tfSeq.getName()); fSeq.setAligned(tfSeq.getAligned());
                rSeq.setName(trSeq.getName()); rSeq.setAligned(trSeq.getAligned());
                if (thisfqualindexfile != "") {
                    fQual = new QualityScores(inFQualIndex); m->gobble(inFQualIndex);
                    rQual = new QualityScores(inRQualIndex); m->gobble(inRQualIndex);
                    savedFQual = new QualityScores(fQual->getName(), fQual->getQualityScores());
                    savedRQual = new QualityScores(rQual->getName(), rQual->getQualityScores());
                    if (fQual->getName() != tfSeq.getName()) { m->mothurOut("[WARNING]: name mismatch in forward quality file. Ignoring, " + tfSeq.getName() + ".\n"); ignore = true; }
                    if (rQual->getName() != trSeq.getName()) { m->mothurOut("[WARNING]: name mismatch in reverse quality file. Ignoring, " + trSeq.getName() + ".\n"); ignore = true; }
                }
                if (tfSeq.getName() != trSeq.getName()) { m->mothurOut("[WARNING]: name mismatch in forward and reverse fasta file. Ignoring, " + tfSeq.getName() + ".\n"); ignore = true; }
            }
            
            if (!ignore) {
                int barcodeIndex = 0;
                int primerIndex = 0;
                Sequence savedFSeq(fSeq.getName(), fSeq.getAligned());  Sequence savedRSeq(rSeq.getName(), rSeq.getAligned());
                Sequence savedFindex(findexBarcode.getName(), findexBarcode.getAligned()); Sequence savedRIndex(rindexBarcode.getName(), rindexBarcode.getAligned());
                
                if(numBarcodes != 0){
                    vector<int> results;
                    if (hasQuality) {
                        if (hasIndex) {
                            results = trimOligos.stripBarcode(findexBarcode, rindexBarcode, *fQual, *rQual, barcodeIndex);
                        }else {
                            results = trimOligos.stripBarcode(fSeq, rSeq, *fQual, *rQual, barcodeIndex);
                        }
                    }else {
                        results = trimOligos.stripBarcode(fSeq, rSeq, barcodeIndex);
                    }
                    success = results[0] + results[2];
                    commentString += "fbdiffs=" + toString(results[0]) + "(" + trimOligos.getCodeValue(results[1], bdiffs) + "), rbdiffs=" + toString(results[2]) + "(" + trimOligos.getCodeValue(results[3], bdiffs) + ") ";
                    if(success > bdiffs)		{	trashCode += 'b';	}
                    else{ currentSeqsDiffs += success;  }
                }
                
                if(numFPrimers != 0){
                    vector<int> results;
                    if (hasQuality) {
                        results = trimOligos.stripForward(fSeq, rSeq, *fQual, *rQual, primerIndex);
                    }else {
                        results = trimOligos.stripForward(fSeq, rSeq, primerIndex);
                    }
                    success = results[0] + results[2];
                    commentString += "fpdiffs=" + toString(results[0]) + "(" + trimOligos.getCodeValue(results[1], pdiffs) + "), rpdiffs=" + toString(results[2]) + "(" + trimOligos.getCodeValue(results[3], pdiffs) + ") ";
                    if(success > pdiffs)		{	trashCode += 'f';	}
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
                    
                    if(numBarcodes != 0){
                        vector<int> results;
                        if (hasQuality) {
                            if (hasIndex) {
                                results = rtrimOligos->stripBarcode(savedFindex, savedRIndex, *savedFQual, *savedRQual, thisBarcodeIndex);
                            }else {
                                results = rtrimOligos->stripBarcode(savedFSeq, savedRSeq, *savedFQual, *savedRQual, thisBarcodeIndex);
                            }
                        }else {
                            results = rtrimOligos->stripBarcode(savedFSeq, savedRSeq, thisBarcodeIndex);
                        }
                        thisSuccess = results[0] + results[2];
                        thiscommentString += "fbdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], bdiffs) + "), rbdiffs=" + toString(results[2]) + "(" + rtrimOligos->getCodeValue(results[3], bdiffs) + ") ";
                        if(thisSuccess > bdiffs)		{	thisTrashCode += 'b';	}
                        else{ thisCurrentSeqsDiffs += thisSuccess;  }
                    }
                    
                    if(numFPrimers != 0){
                        vector<int> results;
                        if (hasQuality) {
                            results = rtrimOligos->stripForward(savedFSeq, savedRSeq, *savedFQual, *savedRQual, thisPrimerIndex);
                        }else {
                            results = rtrimOligos->stripForward(savedFSeq, savedRSeq, thisPrimerIndex);
                        }
                        thisSuccess = results[0] + results[2];
                        thiscommentString += "fpdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], pdiffs) + "), rpdiffs=" + toString(results[2]) + "(" + rtrimOligos->getCodeValue(results[3], pdiffs) + ") ";
                        if(thisSuccess > pdiffs)		{	thisTrashCode += 'f';	}
                        else{ thisCurrentSeqsDiffs += thisSuccess;  }
                    }
                    
                    if (thisCurrentSeqsDiffs > tdiffs)	{	thisTrashCode += 't';   }
                    
                    if (thisTrashCode == "") {
                        trashCode = thisTrashCode;
                        success = thisSuccess;
                        currentSeqsDiffs = thisCurrentSeqsDiffs;
                        commentString = thiscommentString;
                        barcodeIndex = thisBarcodeIndex;
                        primerIndex = thisPrimerIndex;
                        savedFSeq.reverseComplement();
                        savedRSeq.reverseComplement();
                        fSeq.setAligned(savedFSeq.getAligned());
                        rSeq.setAligned(savedRSeq.getAligned());
                        if(hasQuality){
                            savedFQual->flipQScores(); savedRQual->flipQScores();
                            fQual->setScores(savedFQual->getScores()); rQual->setScores(savedRQual->getScores());
                        }
                    }else { trashCode += "(" + thisTrashCode + ")";  }
                }
                
                
                //flip the reverse reads
                rSeq.reverseComplement();
                if (hasQuality) { rQual->flipQScores(); }
                
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
                vector<int> scores1, scores2, contigScores;
                if (hasQuality) {
                    scores1 = fQual->getQualityScores();
                    scores2 = rQual->getQualityScores();
                    delete fQual; delete rQual;  delete savedFQual; delete savedRQual;
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
                        if (hasQuality) {
                            if (scores2[BBaseMap[i]] <= insert) { } //
                            else { contig += seq2[i];  }
                        }else { contig += seq2[i]; } //with no quality info, then we keep it?
                    }else if (((seq2[i] == '.') || (seq2[i] == '-')) && ((seq1[i] != '-') && (seq1[i] != '.'))) { //seq2 is a gap and seq1 is a base, choose seq1, unless quality score for base is below insert. In that case eliminate base
                        if (hasQuality) {
                            if (scores1[ABaseMap[i]] <= insert) { } //
                            else { contig += seq1[i];  }
                        }else { contig += seq1[i]; } //with no quality info, then we keep it?
                    }else if (((seq1[i] != '-') && (seq1[i] != '.')) && ((seq2[i] != '-') && (seq2[i] != '.'))) { //both bases choose one with better quality
                        if (hasQuality) {
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
                        string thisGroup = oligos->getGroupName(barcodeIndex, primerIndex);
                        if (m->debug) { m->mothurOut(", group= " + thisGroup + "\n"); }
                        
                        int pos = thisGroup.find("ignore");
                        if (pos == string::npos) {
                            groupMap[fSeq.getName()] = thisGroup;
                            
                            map<string, int>::iterator it = groupCounts.find(thisGroup);
                            if (it == groupCounts.end()) {	groupCounts[thisGroup] = 1; }
                            else { groupCounts[it->first] ++; }
                        }else { ignore = true; }
                    }else if (createFileGroup) { //for 3 column file option
                        int pos = group.find("ignore");
                        if (pos == string::npos) {
                            groupMap[fSeq.getName()] = group;
                            
                            map<string, int>::iterator it = groupCounts.find(group);
                            if (it == groupCounts.end()) {	groupCounts[group] = 1; }
                            else { groupCounts[it->first] ++; }
                        }else { ignore = true; }
                    }
                    if (m->debug) { m->mothurOut("\n"); }
                    
                    if(!ignore){
                        //output
                        outFasta << ">" << fSeq.getName() << '\t' << commentString << endl << contig << endl;
                        outQual << ">" << fSeq.getName() << '\t' << commentString << endl;
                        for (int i = 0; i < contigScores.size(); i++) { outQual << contigScores[i] << " "; }  outQual << endl;
                        
                        int numNs = 0;
                        for (int i = 0; i < contig.length(); i++) { if (contig[i] == 'N') { numNs++; }  }
                        outMisMatch << fSeq.getName() << '\t' << contig.length() << '\t' << (oend-oStart) << '\t' << oStart << '\t' << oend << '\t' << numMismatches << '\t' << numNs << endl;
                        
                        if (allFiles) {
                            ofstream output;
                            m->openOutputFileAppend(fastaFileNames[barcodeIndex][primerIndex], output);
                            output << ">" << fSeq.getName() << '\t' << commentString << endl << contig << endl;
                            output.close();
                            
                            ofstream output2;
                            m->openOutputFileAppend(qualFileNames[barcodeIndex][primerIndex], output2);
                            output2 << ">" << fSeq.getName() << '\t' << commentString << endl;
                            for (int i = 0; i < contigScores.size(); i++) { output2 << contigScores[i] << " "; }  output2 << endl;
                            output2.close();
                        }
                    }
                }else {
                    //output
                    outScrapFasta << ">" << fSeq.getName() << " | " << trashCode << '\t' << commentString << endl << contig << endl;
                    outScrapQual << ">" << fSeq.getName() << '\t' << commentString << endl;
                    for (int i = 0; i < contigScores.size(); i++) { outScrapQual << contigScores[i] << " "; }  outScrapQual << endl;
                }
            }
            num++;
            
            #if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
                unsigned long long pos = inFFasta.tellg();
                if ((pos == -1) || (pos >= linesInput.end)) { break; }
            #else
                if (inFFasta.eof()) { break; }
            #endif
            
			//report progress
            if((num) % 1000 == 0){	m->mothurOutJustToScreen(toString(num)); m->mothurOutEndLine();	 }
		}
        
		//report progress
		if((num) % 1000 != 0){	m->mothurOut(toString(num)); m->mothurOutEndLine();		}
        
        inFFasta.close();
        inRFasta.close();
        outFasta.close();
        outScrapFasta.close();
        outMisMatch.close();
        if (delim == '@') {
            if (thisfqualindexfile != "") { inFQualIndex.close(); }
            if (thisrqualindexfile != "") { inRQualIndex.close(); }
            outQual.close();
            outScrapQual.close();
        }else{
            if (hasQuality) {
                inFQualIndex.close();
                inRQualIndex.close();
                outQual.close();
                outScrapQual.close();
            }
        }
        delete alignment;
        if (reorient) { delete rtrimOligos; }
        
        if (m->control_pressed) {  m->mothurRemove(outputFasta); m->mothurRemove(outputScrapFasta);m->mothurRemove(outputMisMatches);  }
    
        return num;
    }
	catch(exception& e) {
		m->errorOut(e, "MakeContigsCommand", "driver");
		exit(1);
	}
}
/**************************************************************************************************/

int MakeContigsCommand::setLines(vector<string> fasta, vector<string> qual, vector<linePair>& lines, vector<linePair>& qLines, char delim) {
    try {
        lines.clear();
        qLines.clear();
        vector<unsigned long long> fastaFilePos;
        vector<unsigned long long> qfileFilePos;
        vector<unsigned long long> temp;
        
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        //set file positions for fasta file
        fastaFilePos = m->divideFile(fasta[0], processors, delim);
        
        //get name of first sequence in each chunk
        map<string, int> firstSeqNames;
        for (int i = 0; i < (fastaFilePos.size()-1); i++) {
            ifstream in;
            m->openInputFile(fasta[0], in);
            in.seekg(fastaFilePos[i]);
            string name = "";
            
            if (delim == '>') {
                Sequence temp(in);
                name = temp.getName();
            }else {
                string line = m->getline(in); m->gobble(in);
                vector<string> pieces = m->splitWhiteSpace(line);
                name = pieces[0];
                name = name.substr(1);
                m->checkName(name);
            }
            firstSeqNames[name] = i;
            
            in.close();
        }
        
        map<string, int> copy;
        if (qual.size() != 0) { copy = firstSeqNames; }
        
        //look for match in reverse file
        ifstream in2;
        m->openInputFile(fasta[1], in2);
        
        string input;
        while(!in2.eof()){
            input = m->getline(in2); m->gobble(in2);
            
            if (input.length() != 0) {
                if(input[0] == delim){ //this is a name line
                    vector<string> pieces = m->splitWhiteSpace(input);
                    string name = pieces[0];
                    name = name.substr(1);
                    m->checkName(name);
                    
                    map<string, int>::iterator it = firstSeqNames.find(name);
                    
                    if(it != firstSeqNames.end()) { //this is the start of a new chunk
                        unsigned long long pos = in2.tellg();
                        qfileFilePos.push_back(pos - input.length() - 1);
                        firstSeqNames.erase(it);
                    }
                }
            }
            
            if (firstSeqNames.size() == 0) { break; }
        }
        in2.close();
        
        //get last file position of reverse fasta[1]
        FILE * pFile;
        unsigned long long size;
        
        //get num bytes in file
        pFile = fopen (fasta[1].c_str(),"rb");
        if (pFile==NULL) perror ("Error opening file");
        else{
            fseek (pFile, 0, SEEK_END);
            size=ftell (pFile);
            fclose (pFile);
        }
        qfileFilePos.push_back(size);

        
        if (firstSeqNames.size() != 0) {
            for (map<string, int>::iterator it = firstSeqNames.begin(); it != firstSeqNames.end(); it++) {
                if (delim == '>') {
                    m->mothurOut(it->first + " is in your forward fasta file and not in your reverse file, please remove it using the remove.seqs command before proceeding."); m->mothurOutEndLine();
                }else {
                    m->mothurOut(it->first + " is in your forward fastq file and not in your reverse file, please remove it using the remove.seqs command before proceeding."); m->mothurOutEndLine();
                }
            }
            m->control_pressed = true;
            return processors;
        }

        //fill lines with paired forward and reverse fasta lines
        for (int i = 0; i < (fastaFilePos.size()-1); i++) {
            if (m->debug) { m->mothurOut("[DEBUG]: forward " + toString(i) +'\t' + toString(fastaFilePos[i]) + '\t' + toString(fastaFilePos[i+1]) + '\n'); }
            lines.push_back(linePair(fastaFilePos[i], fastaFilePos[(i+1)]));
            if (m->debug) { m->mothurOut("[DEBUG]: reverse " + toString(i) +'\t' + toString(qfileFilePos[i]) + '\t' + toString(qfileFilePos[i+1]) + '\n'); }
            lines.push_back(linePair(qfileFilePos[i], qfileFilePos[(i+1)]));
        }
        
        qfileFilePos.clear();
        
        if (qual.size() != 0) {
            firstSeqNames = copy;
            
            if (qual[0] != "NONE") {
                //seach for filePos of each first name in the qfile and save in qfileFilePos
                ifstream inQual;
                m->openInputFile(qual[0], inQual);
                
                string input;
                while(!inQual.eof()){
                    input = m->getline(inQual); m->gobble(inQual);
                    
                    if (input.length() != 0) {
                        if(input[0] == delim){ //this is a sequence name line
                            vector<string> pieces = m->splitWhiteSpace(input);
                            string name = pieces[0];
                            name = name.substr(1);
                            m->checkName(name);
                            
                            map<string, int>::iterator it = firstSeqNames.find(name);
                            
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
                
                //get last file position of reverse qual[0]
                FILE * pFile;
                unsigned long long size;
                
                //get num bytes in file
                pFile = fopen (qual[0].c_str(),"rb");
                if (pFile==NULL) perror ("Error opening file");
                else{
                    fseek (pFile, 0, SEEK_END);
                    size=ftell (pFile);
                    fclose (pFile);
                }
                qfileFilePos.push_back(size);
                
                
                if (firstSeqNames.size() != 0) {
                    for (map<string, int>::iterator it = firstSeqNames.begin(); it != firstSeqNames.end(); it++) {
                        if (delim == '>') {
                            m->mothurOut(it->first + " is in your forward fasta file and reverse fasta file, but not your forward qfile, please remove it using the remove.seqs command before proceeding."); m->mothurOutEndLine();
                        }else {
                            m->mothurOut(it->first + " is in your forward fastq file and reverse fastq file, but not your forward index, please remove it using the remove.seqs command before proceeding."); m->mothurOutEndLine();
                        }
                    }
                    m->control_pressed = true;
                    return processors;
                }
            }
            firstSeqNames = copy;
            
            if (qual[1] != "NONE") {
                ifstream inQual2;
                m->openInputFile(qual[1], inQual2);
                
                while(!inQual2.eof()){
                    input = m->getline(inQual2); m->gobble(inQual2);
                    
                    if (input.length() != 0) {
                        if(input[0] == delim){ //this is a sequence name line
                            vector<string> pieces = m->splitWhiteSpace(input);
                            string name = pieces[0];
                            name = name.substr(1);
                            
                            m->checkName(name);
                            
                            map<string, int>::iterator it = firstSeqNames.find(name);
                            
                            if(it != firstSeqNames.end()) { //this is the start of a new chunk
                                unsigned long long pos = inQual2.tellg();
                                temp.push_back(pos - input.length() - 1);
                                firstSeqNames.erase(it);
                            }
                        }
                    }
                    
                    if (firstSeqNames.size() == 0) { break; }
                }
                inQual2.close();
                
                //get last file position of reverse qual[1]
                FILE * pFile2;
                
                //get num bytes in file
                pFile2 = fopen (qual[1].c_str(),"rb");
                if (pFile2==NULL) perror ("Error opening file");
                else{
                    fseek (pFile2, 0, SEEK_END);
                    size=ftell (pFile2);
                    fclose (pFile2);
                }
                temp.push_back(size);
                
                
                if (firstSeqNames.size() != 0) {
                    for (map<string, int>::iterator it = firstSeqNames.begin(); it != firstSeqNames.end(); it++) {
                        if (delim == '>') {
                            m->mothurOut(it->first + " is in your forward fasta file, reverse fasta file, and forward qfile but not your reverse qfile, please remove it using the remove.seqs command before proceeding."); m->mothurOutEndLine();
                        }else {
                            if (qual[0] != "NONE") {
                                m->mothurOut(it->first + " is in your forward fastq file, reverse fastq file, and forward index but not your reverse index, please remove it using the remove.seqs command before proceeding."); m->mothurOutEndLine();
                            }else {
                                m->mothurOut(it->first + " is in your forward fastq file, reverse fastq file, but not your reverse index, please remove it using the remove.seqs command before proceeding."); m->mothurOutEndLine();
                            }
                        }
                    }
                    m->control_pressed = true;
                    return processors;
                }
            }
            
            if (qual[0] == "NONE") { qfileFilePos = temp; } //fill with duds, if both were NONE then qual.size() == 0
            if (qual[1] == "NONE") { temp = qfileFilePos; } //fill with duds, if both were NONE then qual.size() == 0
            
            
            //fill lines with paired forward and reverse fasta lines
            for (int i = 0; i < (fastaFilePos.size()-1); i++) {
                if (m->debug) { m->mothurOut("[DEBUG]: forward " + toString(i) +'\t' + toString(qfileFilePos[i]) + '\t' + toString(qfileFilePos[i+1]) + '\n'); }
                qLines.push_back(linePair(qfileFilePos[i], qfileFilePos[(i+1)]));
                if (m->debug) { m->mothurOut("[DEBUG]: reverse " + toString(i) +'\t' + toString(temp[i]) + '\t' + toString(temp[i+1]) + '\n'); }
                qLines.push_back(linePair(temp[i], temp[(i+1)]));
            }
        }else {  qLines = lines;	} //files with duds
        
        
        return processors;
        
#else
        
        if (processors == 1) { //save time
            //fastaFilePos.push_back(0); qfileFilePos.push_back(0);
            //fastaFilePos.push_back(1000); qfileFilePos.push_back(1000);
            lines.push_back(linePair(0, 1000)); lines.push_back(linePair(0, 1000)); //fasta[0], fasta[1] - forward and reverse
            if (qual.size() != 0) {  qLines.push_back(linePair(0, 1000)); qLines.push_back(linePair(0, 1000)); } //qual[0], qual[1] - forward and reverse
        }else{
            unsigned long long numFastaSeqs = 0;
            fastaFilePos = m->setFilePosFasta(fasta[0], numFastaSeqs, delim); //forward
            if (fastaFilePos.size() < processors) { processors = fastaFilePos.size(); }
            
            unsigned long long numRFastaSeqs = 0;
            qfileFilePos = m->setFilePosFasta(fasta[1], numRFastaSeqs, delim); //reverse
            
            if (numFastaSeqs != numRFastaSeqs) {
                if (delim == '>') {
                    m->mothurOut("[ERROR]: You have " + toString(numFastaSeqs) + " sequences in your forward fasta file, but " + toString(numRFastaSeqs) + " sequences in your reverse fasta file. Please use the list.seqs and get.seqs commands to make the files match before proceeding."); m->mothurOutEndLine(); m->control_pressed = true; return processors;
                }else {
                    m->mothurOut("[ERROR]: You have " + toString(numFastaSeqs) + " sequences in your forward fastq file, but " + toString(numRFastaSeqs) + " sequences in your reverse fastq file. Please use the list.seqs and get.seqs commands to make the files match before proceeding."); m->mothurOutEndLine(); m->control_pressed = true; return processors;
                }
            }
            
            //figure out how many sequences you have to process
            unsigned long long numSeqsPerProcessor = numFastaSeqs / processors;
            for (int i = 0; i < processors; i++) {
                unsigned long long startIndex =  i * numSeqsPerProcessor;
                if(i == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor; 	}
                lines.push_back(linePair(fastaFilePos[startIndex], numSeqsPerProcessor)); //forward
                lines.push_back(linePair(qfileFilePos[startIndex], numSeqsPerProcessor)); //reverse
            }
            
            
            if (qual.size() != 0) {
                long long numFQualSeqs = 0;
                long long numRQualSeqs = 0;
                fastaFilePos.clear();
                qfileFilePos.clear();
                
                if (qual[0] != "NONE") {  fastaFilePos = m->setFilePosFasta(qual[0], numFQualSeqs, '>');  } //forward index or qual file
                if (qual[1] != "NONE") {  qfileFilePos = m->setFilePosFasta(qual[1], numRQualSeqs, '>');  }//reverse index or qual file
                
                if (qual[0] == "NONE") { fastaFilePos = qfileFilePos; numFQualSeqs = numRQualSeqs; } //fill with duds, if both were NONE then qual.size() == 0
                if (qual[1] == "NONE") { qfileFilePos = fastaFilePos; numRQualSeqs = numFQualSeqs; } //fill with duds, if both were NONE then qual.size() == 0

                
                if ((numFQualSeqs != numRQualSeqs) || (numFQualSeqs != numFastaSeqs)){
                    if (delim == '>') {
                        m->mothurOut("[ERROR]: You have " + toString(numFastaSeqs) + " sequences in your forward fasta file, " + toString(numRFastaSeqs) + " sequences in your reverse fasta file, " + toString(numFQualSeqs) + " sequences in your forward qual file, " + toString(numRQualSeqs) + " sequences in your reverse qual file. Please use the list.seqs and get.seqs commands to make the files match before proceeding."); m->mothurOutEndLine(); m->control_pressed = true; return processors;
                    }else {
                        if (qual[0] != "NONE") {
                            m->mothurOut("[ERROR]: You have " + toString(numFastaSeqs) + " sequences in your forward fastq file, " + toString(numRFastaSeqs) + " sequences in your reverse fastq file and " + toString(numRQualSeqs) + " sequences in your reverse index file. Please use the list.seqs and get.seqs commands to make the files match before proceeding."); m->mothurOutEndLine(); m->control_pressed = true; return processors;
                        }else if (qual[1] != "NONE") {
                            m->mothurOut("[ERROR]: You have " + toString(numFastaSeqs) + " sequences in your forward fastq file, " + toString(numRFastaSeqs) + " sequences in your reverse fastq file and " + toString(numFQualSeqs) + " sequences in your forward index file. Please use the list.seqs and get.seqs commands to make the files match before proceeding."); m->mothurOutEndLine(); m->control_pressed = true; return processors;
                        }else {
                            m->mothurOut("[ERROR]: You have " + toString(numFastaSeqs) + " sequences in your forward fastq file, " + toString(numRFastaSeqs) + " sequences in your reverse fastq file, " + toString(numFQualSeqs) + " sequences in your forward index file, " + toString(numRQualSeqs) + " sequences in your reverse index file. Please use the list.seqs and get.seqs commands to make the files match before proceeding."); m->mothurOutEndLine(); m->control_pressed = true; return processors;
                        }
                    }
                }
                
                //figure out how many sequences you have to process
                unsigned long long numSeqsPerProcessor = numFQualSeqs / processors;
                for (int i = 0; i < processors; i++) {
                    unsigned long long startIndex =  i * numSeqsPerProcessor;
                    if(i == (processors - 1)){	numSeqsPerProcessor = numFQualSeqs - i * numSeqsPerProcessor; 	}
                    qLines.push_back(linePair(fastaFilePos[startIndex], numSeqsPerProcessor)); //forward
                    qLines.push_back(linePair(qfileFilePos[startIndex], numSeqsPerProcessor)); //reverse
                }
  
            }else { qLines = lines;	} //files with duds
        }
        if(qfilename == "")	{	qLines = lines;	} //files with duds
        return 1;
        
#endif
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "setLines");
        exit(1);
    }
}

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
                if ((findex == "none") || (findex == "NONE")){ findex = "NONE"; }
                if ((rindex == "none") || (rindex == "NONE")){ rindex = "NONE"; }
            }else {
                m->mothurOut("[ERROR]: file lines can be 2, 3, or 4 columns. The forward fastq files in the first column and their matching reverse fastq files in the second column, or a groupName then forward fastq file and reverse fastq file, or forward fastq file then reverse fastq then forward index and reverse index file.  If you only have one index file add 'none' for the other one. \n"); m->control_pressed = true;
            }
            
            if (m->debug) { m->mothurOut("[DEBUG]: group = " + group + ", forward = " + forward + ", reverse = " + reverse + ", forwardIndex = " + findex + ", reverseIndex = " + rindex + ".\n"); }
            
            if (inputDir != "") {
                string path = m->hasPath(forward);
                if (path == "") {  forward = inputDir + forward;  }
                
                path = m->hasPath(reverse);
                if (path == "") {  reverse = inputDir + reverse;  }
                
                if (findex != "") {
                    path = m->hasPath(findex);
                    if (path == "") {  findex = inputDir + findex;  }
                }
                
                if (rindex != "") {
                    path = m->hasPath(rindex);
                    if (path == "") {  rindex = inputDir + rindex;  }
                }
            }
            
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
            if ((findex != "") && (findex != "NONE")){
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
            if ((rindex != "") && (rindex != "NONE")) {
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
bool MakeContigsCommand::getOligos(vector<vector<string> >& fastaFileNames, vector<vector<string> >& qualFileNames, string rootname, map<string, string>& fastaFile2Group){
	try {
        if (m->debug) { m->mothurOut("[DEBUG]: oligosfile = " + oligosfile + "\n"); }
        
        bool allBlank = false;
        oligos->read(oligosfile, false);
        
        if (m->control_pressed) { return false; } //error in reading oligos
        
        if (oligos->hasPairedBarcodes()) {
            numFPrimers = oligos->getPairedPrimers().size();
            numBarcodes = oligos->getPairedBarcodes().size();
        }else {
            m->mothurOut("[ERROR]: make.contigs requires paired barcodes and primers. You can set one end to NONE if you are using an index file.\n"); m->control_pressed = true;
        }
    
        if (m->control_pressed) { return false; }
    
        numLinkers = oligos->getLinkers().size();
        numSpacers = oligos->getSpacers().size();
        numRPrimers = oligos->getReversePrimers().size();
        if (numLinkers != 0) { m->mothurOut("[WARNING]: make.contigs is not setup to remove linkers, ignoring.\n"); }
        if (numSpacers != 0) { m->mothurOut("[WARNING]: make.contigs is not setup to remove spacers, ignoring.\n"); }
       
        vector<string> groupNames = oligos->getGroupNames();
        if (groupNames.size() == 0) { allFiles = 0; allBlank = true;  }
        
        
        fastaFileNames.resize(oligos->getBarcodeNames().size());
		for(int i=0;i<fastaFileNames.size();i++){
            for(int j=0;j<oligos->getPrimerNames().size();j++){  fastaFileNames[i].push_back(""); }
		}
        
        qualFileNames = fastaFileNames;
        
        if (allFiles) {
            set<string> uniqueNames; //used to cleanup outputFileNames
            map<int, oligosPair> barcodes = oligos->getPairedBarcodes();
            map<int, oligosPair> primers = oligos->getPairedPrimers();
            for(map<int, oligosPair>::iterator itBar = barcodes.begin();itBar != barcodes.end();itBar++){
                for(map<int, oligosPair>::iterator itPrimer = primers.begin();itPrimer != primers.end(); itPrimer++){
                    
                    string primerName = oligos->getPrimerName(itPrimer->first);
                    string barcodeName = oligos->getBarcodeName(itBar->first);
                    
                    if ((primerName == "ignore") || (barcodeName == "ignore")) { } //do nothing
                    else if ((primerName == "") && (barcodeName == "")) { } //do nothing
                    else {
                        string comboGroupName = "";
                        string fastaFileName = "";
                        string qualFileName = "";
                        
                        if(primerName == ""){
                            comboGroupName = barcodeName;
                        }else{
                            if(barcodeName == ""){
                                comboGroupName = primerName;
                            }
                            else{
                                comboGroupName = barcodeName + "." + primerName;
                            }
                        }
                        
                        
                        ofstream temp, temp2;
                        map<string, string> variables;
                        variables["[filename]"] = rootname;
                        variables["[tag]"] = comboGroupName;
                        fastaFileName = getOutputFileName("fasta", variables);
                        qualFileName = getOutputFileName("qfile", variables);
                        if (uniqueNames.count(fastaFileName) == 0) {
                            outputNames.push_back(fastaFileName);
                            outputTypes["fasta"].push_back(fastaFileName);
                            uniqueNames.insert(fastaFileName);
                            fastaFile2Group[fastaFileName] = comboGroupName;
                            
                            outputNames.push_back(qualFileName);
                            outputTypes["qfile"].push_back(qualFileName);
                            uniqueNames.insert(qualFileName);

                        }
                        
                        fastaFileNames[itBar->first][itPrimer->first] = fastaFileName;
                        m->openOutputFile(fastaFileName, temp);		temp.close();
                        //cout << fastaFileName << endl;
                        
                        qualFileNames[itBar->first][itPrimer->first] = qualFileName;
                        m->openOutputFile(qualFileName, temp2);		temp2.close();
                    }
                }
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
//**********************************************************************************************************************

