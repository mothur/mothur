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
        CommandParameter palign("align", "Multiple", "needleman-gotoh-kmer", "needleman", "", "", "","",false,false); parameters.push_back(palign);
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
        CommandParameter pksize("ksize", "Number", "", "8", "", "", "","",false,false); parameters.push_back(pksize);
        CommandParameter pkmer("kmer", "Number", "", "2", "", "", "","",false,false); parameters.push_back(pkmer);
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
        helpString += "The align parameter allows you to specify the alignment method to use.  Your options are: kmer, gotoh and needleman. The default is needleman.\n";
        helpString += "The ksize parameter allows you to set the kmer size if you are doing align=kmer. Default=8.\n";
        helpString += "The kmer parameter allows you to set the number of sequence locations for a particular k-mer. When attempting to align the sequences, mothur will store the location of every k-mer in a table. If the same k-mer is present multiple times, only the first ones will be stored until the table is full; when this occurs, an FML error is emitted. If the sequences are highly repetitive, lost positions can prevent good alignments; this can be alleviated by increasing this amount. This should be small (no more than 10; the default is 2), or the k-mer table will be extremely large, using a large amount of RAM per thread. Try increasing the value until FML errors go away.. Default=2.\n";
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
            
            temp = validParameter.validFile(parameters, "ksize", false);	if (temp == "not found"){	temp = "7";			}
            m->mothurConvert(temp, kmerSize);
            
            temp = validParameter.validFile(parameters, "kmer", false);	if (temp == "not found"){	temp = "2";			}
            m->mothurConvert(temp, numKmers);

            
            temp = validParameter.validFile(parameters, "trimoverlap", false);		if (temp == "not found") { temp = "F"; }
			trimOverlap = m->isTrue(temp);
			
			align = validParameter.validFile(parameters, "align", false);		if (align == "not found"){	align = "needleman";	}
			if ((align != "needleman") && (align != "gotoh") && (align != "kmer")) { m->mothurOut(align + " is not a valid alignment method. Options are kmer, needleman or gotoh. I will use needleman."); m->mothurOutEndLine(); align = "needleman"; }
            
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
			contigsData* tempcontig = new contigsData(format, delim, group, fileInputs, qualOrIndexFiles, (outputFasta + extension), (outputScrapFasta + extension), (outputQual + extension), (outputScrapQual + extension), (outputMisMatches + extension), align, m, match, misMatch, gapOpen, gapExtend, insert, deltaq, tempFASTAFileNames, tempQUALFileNames, oligosfile, reorient, pdiffs, bdiffs, tdiffs, kmerSize, createOligosGroup, createFileGroup, allFiles, trimOverlap, lines[spot].start, lines[spot].end, lines[spot+1].start, lines[spot+1].end, qLines[spot].start, qLines[spot].end, qLines[spot+1].start, qLines[spot+1].end, h);
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
        const double qual_match_simple_bayesian[][47] = {
            { -1.09861, -1.32887, -1.55913, -1.78939, -2.01965, -2.2499, -2.48016, -2.71042, -2.94068, -3.17094, -3.4012, -3.63146, -3.86171, -4.09197, -4.32223, -4.55249, -4.78275, -5.01301, -5.24327, -5.47352, -5.70378, -5.93404, -6.1643, -6.39456, -6.62482, -6.85508, -7.08533, -7.31559, -7.54585, -7.77611, -8.00637, -8.23663, -8.46688, -8.69714, -8.9274, -9.15766, -9.38792, -9.61818, -9.84844, -10.0787, -10.309, -10.5392, -10.7695, -10.9997, -11.23, -11.4602, -11.6905},
            { -1.32887, -1.37587, -1.41484, -1.44692, -1.47315, -1.49449, -1.51178, -1.52572, -1.53694, -1.54593, -1.55314, -1.5589, -1.5635, -1.56717, -1.5701, -1.57243, -1.57428, -1.57576, -1.57693, -1.57786, -1.5786, -1.57919, -1.57966, -1.58003, -1.58033, -1.58057, -1.58075, -1.5809, -1.58102, -1.58111, -1.58119, -1.58125, -1.58129, -1.58133, -1.58136, -1.58138, -1.5814, -1.58142, -1.58143, -1.58144, -1.58145, -1.58145, -1.58146, -1.58146, -1.58146, -1.58146, -1.58147},
            { -1.55913, -1.41484, -1.31343, -1.23963, -1.18465, -1.14303, -1.11117, -1.08657, -1.06744, -1.05251, -1.0408, -1.0316, -1.02436, -1.01863, -1.01411, -1.01054, -1.00771, -1.00546, -1.00368, -1.00227, -1.00115, -1.00027, -0.99956, -0.999001, -0.998557, -0.998204, -0.997924, -0.997702, -0.997525, -0.997385, -0.997273, -0.997185, -0.997114, -0.997059, -0.997014, -0.996979, -0.996951, -0.996929, -0.996911, -0.996897, -0.996886, -0.996877, -0.99687, -0.996865, -0.99686, -0.996857, -0.996854},
            { -1.78939, -1.44692, -1.23963, -1.10098, -1.0031, -0.931648, -0.878319, -0.837896, -0.806912, -0.782967, -0.764347, -0.7498, -0.738394, -0.729426, -0.722359, -0.71678, -0.712372, -0.708883, -0.706121, -0.703933, -0.702197, -0.700821, -0.69973, -0.698863, -0.698176, -0.69763, -0.697196, -0.696852, -0.696579, -0.696362, -0.69619, -0.696053, -0.695944, -0.695858, -0.695789, -0.695735, -0.695692, -0.695657, -0.69563, -0.695608, -0.695591, -0.695577, -0.695566, -0.695558, -0.695551, -0.695546, -0.695541},
            { -2.01965, -1.47315, -1.18465, -1.0031, -0.879224, -0.790712, -0.725593, -0.676729, -0.639547, -0.610968, -0.588834, -0.571596, -0.558111, -0.547528, -0.539201, -0.532636, -0.527451, -0.523352, -0.520107, -0.517538, -0.515502, -0.513887, -0.512606, -0.51159, -0.510784, -0.510144, -0.509636, -0.509232, -0.508912, -0.508658, -0.508456, -0.508295, -0.508168, -0.508067, -0.507986, -0.507922, -0.507872, -0.507831, -0.507799, -0.507774, -0.507754, -0.507738, -0.507725, -0.507715, -0.507707, -0.507701, -0.507695},
            { -2.2499, -1.49449, -1.14303, -0.931648, -0.790712, -0.691393, -0.618979, -0.564976, -0.524066, -0.492723, -0.468507, -0.449682, -0.434976, -0.423448, -0.414384, -0.407243, -0.401606, -0.397151, -0.393627, -0.390836, -0.388625, -0.386872, -0.385482, -0.384379, -0.383503, -0.382809, -0.382257, -0.38182, -0.381472, -0.381196, -0.380977, -0.380803, -0.380664, -0.380554, -0.380467, -0.380398, -0.380343, -0.380299, -0.380264, -0.380237, -0.380215, -0.380198, -0.380184, -0.380173, -0.380164, -0.380157, -0.380152},
            { -2.48016, -1.51178, -1.11117, -0.878319, -0.725593, -0.618979, -0.541714, -0.48433, -0.440984, -0.407844, -0.382281, -0.362431, -0.34694, -0.334804, -0.325268, -0.317757, -0.311831, -0.307149, -0.303445, -0.300513, -0.29819, -0.296348, -0.294888, -0.29373, -0.29281, -0.292081, -0.291502, -0.291042, -0.290677, -0.290387, -0.290157, -0.289974, -0.289829, -0.289713, -0.289622, -0.289549, -0.289491, -0.289445, -0.289409, -0.28938, -0.289357, -0.289339, -0.289324, -0.289313, -0.289304, -0.289296, -0.28929},
            { -2.71042, -1.52572, -1.08657, -0.837896, -0.676729, -0.564976, -0.48433, -0.424604, -0.379581, -0.345208, -0.318723, -0.298173, -0.282146, -0.269595, -0.259737, -0.251976, -0.245853, -0.241016, -0.23719, -0.234162, -0.231763, -0.229861, -0.228354, -0.227158, -0.226208, -0.225455, -0.224857, -0.224383, -0.224006, -0.223707, -0.223469, -0.22328, -0.22313, -0.223011, -0.222917, -0.222842, -0.222782, -0.222734, -0.222697, -0.222667, -0.222643, -0.222624, -0.222609, -0.222597, -0.222588, -0.222581, -0.222575},
            { -2.94068, -1.53694, -1.06744, -0.806912, -0.639547, -0.524066, -0.440984, -0.379581, -0.333359, -0.298107, -0.270966, -0.249919, -0.233512, -0.220668, -0.210582, -0.202642, -0.19638, -0.191434, -0.187522, -0.184426, -0.181973, -0.180029, -0.178488, -0.177265, -0.176295, -0.175525, -0.174914, -0.174428, -0.174043, -0.173737, -0.173494, -0.173301, -0.173148, -0.173026, -0.17293, -0.172853, -0.172792, -0.172744, -0.172705, -0.172675, -0.17265, -0.172631, -0.172616, -0.172604, -0.172594, -0.172586, -0.17258},
            { -3.17094, -1.54593, -1.05251, -0.782967, -0.610968, -0.492723, -0.407844, -0.345208, -0.298107, -0.262213, -0.234592, -0.213183, -0.196498, -0.18344, -0.173188, -0.165119, -0.158755, -0.153729, -0.149755, -0.146609, -0.144117, -0.142143, -0.140577, -0.139335, -0.138349, -0.137567, -0.136946, -0.136453, -0.136062, -0.135751, -0.135504, -0.135308, -0.135153, -0.135029, -0.134931, -0.134853, -0.134791, -0.134742, -0.134703, -0.134672, -0.134647, -0.134628, -0.134612, -0.1346, -0.13459, -0.134582, -0.134576},
            { -3.4012, -1.55314, -1.0408, -0.764347, -0.588834, -0.468507, -0.382281, -0.318723, -0.270966, -0.234592, -0.206614, -0.184935, -0.168044, -0.154827, -0.144451, -0.136285, -0.129846, -0.124761, -0.12074, -0.117558, -0.115037, -0.113039, -0.111455, -0.110198, -0.109202, -0.10841, -0.107782, -0.107284, -0.106888, -0.106574, -0.106324, -0.106126, -0.105968, -0.105843, -0.105744, -0.105665, -0.105602, -0.105553, -0.105513, -0.105482, -0.105457, -0.105437, -0.105421, -0.105409, -0.105399, -0.105391, -0.105385},
            { -3.63146, -1.5589, -1.0316, -0.7498, -0.571596, -0.449682, -0.362431, -0.298173, -0.249919, -0.213183, -0.184935, -0.163052, -0.146004, -0.132667, -0.122198, -0.11396, -0.107464, -0.102334, -0.0982781, -0.0950678, -0.0925252, -0.09051, -0.0889123, -0.0876449, -0.0866394, -0.0858414, -0.0852079, -0.0847051, -0.0843058, -0.0839888, -0.083737, -0.0835371, -0.0833783, -0.0832522, -0.083152, -0.0830725, -0.0830093, -0.0829591, -0.0829192, -0.0828876, -0.0828624, -0.0828425, -0.0828266, -0.082814, -0.082804, -0.082796, -0.0827897},
            { -3.86171, -1.5635, -1.02436, -0.738394, -0.558111, -0.434976, -0.34694, -0.282146, -0.233512, -0.196498, -0.168044, -0.146004, -0.128838, -0.115409, -0.104869, -0.096575, -0.0900357, -0.0848716, -0.0807886, -0.0775572, -0.0749978, -0.0729694, -0.0713612, -0.0700856, -0.0690735, -0.0682703, -0.0676327, -0.0671265, -0.0667247, -0.0664056, -0.0661522, -0.065951, -0.0657912, -0.0656642, -0.0655634, -0.0654833, -0.0654198, -0.0653692, -0.0653291, -0.0652972, -0.0652719, -0.0652518, -0.0652359, -0.0652232, -0.0652131, -0.0652051, -0.0651987},
            { -4.09197, -1.56717, -1.01863, -0.729426, -0.547528, -0.423448, -0.334804, -0.269595, -0.220668, -0.18344, -0.154827, -0.132667, -0.115409, -0.101909, -0.0913142, -0.0829777, -0.0764049, -0.0712146, -0.0671109, -0.0638632, -0.061291, -0.0592525, -0.0576362, -0.0563542, -0.055337, -0.0545298, -0.053889, -0.0533804, -0.0529765, -0.0526558, -0.0524012, -0.0521989, -0.0520383, -0.0519108, -0.0518095, -0.051729, -0.0516651, -0.0516143, -0.051574, -0.051542, -0.0515165, -0.0514963, -0.0514803, -0.0514675, -0.0514574, -0.0514493, -0.051443},
            { -4.32223, -1.5701, -1.01411, -0.722359, -0.539201, -0.414384, -0.325268, -0.259737, -0.210582, -0.173188, -0.144451, -0.122198, -0.104869, -0.0913142, -0.0806768, -0.0723072, -0.0657085, -0.0604979, -0.0563782, -0.0531178, -0.0505356, -0.0484892, -0.0468667, -0.0455797, -0.0445586, -0.0437483, -0.0431051, -0.0425945, -0.0421891, -0.0418671, -0.0416115, -0.0414085, -0.0412473, -0.0411192, -0.0410175, -0.0409368, -0.0408726, -0.0408216, -0.0407812, -0.040749, -0.0407235, -0.0407032, -0.0406871, -0.0406743, -0.0406641, -0.040656, -0.0406496},
            { -4.55249, -1.57243, -1.01054, -0.71678, -0.532636, -0.407243, -0.317757, -0.251976, -0.202642, -0.165119, -0.136285, -0.11396, -0.096575, -0.0829777, -0.0723072, -0.0639118, -0.0572929, -0.0520664, -0.0479342, -0.044664, -0.042074, -0.0400214, -0.038394, -0.0371032, -0.0360791, -0.0352663, -0.0346212, -0.0341091, -0.0337024, -0.0333796, -0.0331232, -0.0329196, -0.0327579, -0.0326294, -0.0325274, -0.0324464, -0.0323821, -0.0323309, -0.0322904, -0.0322581, -0.0322325, -0.0322121, -0.032196, -0.0321831, -0.032173, -0.0321649, -0.0321584},
            { -4.78275, -1.57428, -1.00771, -0.712372, -0.527451, -0.401606, -0.311831, -0.245853, -0.19638, -0.158755, -0.129846, -0.107464, -0.0900357, -0.0764049, -0.0657085, -0.0572929, -0.0506582, -0.0454193, -0.0412773, -0.0379994, -0.0354033, -0.033346, -0.0317148, -0.0304209, -0.0293944, -0.0285798, -0.0279331, -0.0274198, -0.0270122, -0.0266886, -0.0264316, -0.0262275, -0.0260655, -0.0259367, -0.0258345, -0.0257533, -0.0256888, -0.0256376, -0.0255969, -0.0255645, -0.0255389, -0.0255185, -0.0255023, -0.0254894, -0.0254792, -0.0254711, -0.0254646},
            { -5.01301, -1.57576, -1.00546, -0.708883, -0.523352, -0.397151, -0.307149, -0.241016, -0.191434, -0.153729, -0.124761, -0.102334, -0.0848716, -0.0712146, -0.0604979, -0.0520664, -0.0454193, -0.0401706, -0.036021, -0.032737, -0.0301362, -0.028075, -0.0264408, -0.0251447, -0.0241163, -0.0233001, -0.0226523, -0.0221381, -0.0217297, -0.0214055, -0.0211481, -0.0209436, -0.0207812, -0.0206523, -0.0205498, -0.0204685, -0.0204039, -0.0203526, -0.0203118, -0.0202794, -0.0202537, -0.0202333, -0.020217, -0.0202041, -0.0201939, -0.0201858, -0.0201793},
            { -5.24327, -1.57693, -1.00368, -0.706121, -0.520107, -0.393627, -0.303445, -0.23719, -0.187522, -0.149755, -0.12074, -0.0982781, -0.0807886, -0.0671109, -0.0563782, -0.0479342, -0.0412773, -0.036021, -0.0318653, -0.0285766, -0.025972, -0.0239079, -0.0222713, -0.0209733, -0.0199434, -0.0191261, -0.0184774, -0.0179624, -0.0175535, -0.0172288, -0.016971, -0.0167662, -0.0166036, -0.0164745, -0.0163719, -0.0162904, -0.0162257, -0.0161743, -0.0161335, -0.0161011, -0.0160753, -0.0160549, -0.0160386, -0.0160257, -0.0160155, -0.0160073, -0.0160009},
            { -5.47352, -1.57786, -1.00227, -0.703933, -0.517538, -0.390836, -0.300513, -0.234162, -0.184426, -0.146609, -0.117558, -0.0950678, -0.0775572, -0.0638632, -0.0531178, -0.044664, -0.0379994, -0.032737, -0.0285766, -0.0252842, -0.0226766, -0.0206101, -0.0189717, -0.0176722, -0.0166412, -0.015823, -0.0151735, -0.0146579, -0.0142486, -0.0139235, -0.0136654, -0.0134604, -0.0132976, -0.0131684, -0.0130657, -0.0129841, -0.0129193, -0.0128679, -0.012827, -0.0127945, -0.0127688, -0.0127483, -0.012732, -0.0127191, -0.0127088, -0.0127007, -0.0126942},
            { -5.70378, -1.5786, -1.00115, -0.702197, -0.515502, -0.388625, -0.29819, -0.231763, -0.181973, -0.144117, -0.115037, -0.0925252, -0.0749978, -0.061291, -0.0505356, -0.042074, -0.0354033, -0.0301362, -0.025972, -0.0226766, -0.0200667, -0.0179984, -0.0163585, -0.0150578, -0.0140259, -0.0132069, -0.0125569, -0.0120409, -0.0116311, -0.0113058, -0.0110475, -0.0108423, -0.0106794, -0.01055, -0.0104472, -0.0103655, -0.0103007, -0.0102492, -0.0102083, -0.0101758, -0.01015, -0.0101295, -0.0101132, -0.0101003, -0.01009, -0.0100819, -0.0100754},
            { -5.93404, -1.57919, -1.00027, -0.700821, -0.513887, -0.386872, -0.296348, -0.229861, -0.180029, -0.142143, -0.113039, -0.09051, -0.0729694, -0.0592525, -0.0484892, -0.0400214, -0.033346, -0.028075, -0.0239079, -0.0206101, -0.0179984, -0.0159286, -0.0142876, -0.012986, -0.0119533, -0.0111338, -0.0104833, -0.00996692, -0.00955691, -0.00923135, -0.00897283, -0.00876752, -0.00860447, -0.00847497, -0.00837212, -0.00829043, -0.00822555, -0.00817401, -0.00813308, -0.00810056, -0.00807474, -0.00805422, -0.00803793, -0.00802498, -0.0080147, -0.00800654, -0.00800005},
            { -6.1643, -1.57966, -0.99956, -0.69973, -0.512606, -0.385482, -0.294888, -0.228354, -0.178488, -0.140577, -0.111455, -0.0889123, -0.0713612, -0.0576362, -0.0468667, -0.038394, -0.0317148, -0.0264408, -0.0222713, -0.0189717, -0.0163585, -0.0142876, -0.0126457, -0.0113434, -0.0103101, -0.00949014, -0.00883928, -0.00832259, -0.00791235, -0.00758661, -0.00732794, -0.00712252, -0.00695938, -0.00682981, -0.00672691, -0.00664517, -0.00658025, -0.00652869, -0.00648773, -0.0064552, -0.00642936, -0.00640883, -0.00639253, -0.00637958, -0.00636929, -0.00636112, -0.00635463},
            { -6.39456, -1.58003, -0.999001, -0.698863, -0.51159, -0.384379, -0.29373, -0.227158, -0.177265, -0.139335, -0.110198, -0.0876449, -0.0700856, -0.0563542, -0.0455797, -0.0371032, -0.0304209, -0.0251447, -0.0209733, -0.0176722, -0.0150578, -0.012986, -0.0113434, -0.0100405, -0.00900678, -0.00818644, -0.00753529, -0.00701837, -0.00660796, -0.00628208, -0.00602329, -0.00581778, -0.00565457, -0.00552494, -0.00542199, -0.00534022, -0.00527527, -0.00522368, -0.00518271, -0.00515016, -0.00512431, -0.00510378, -0.00508747, -0.00507451, -0.00506422, -0.00505604, -0.00504955},
            { -6.62482, -1.58033, -0.998557, -0.698176, -0.510784, -0.383503, -0.29281, -0.226208, -0.176295, -0.138349, -0.109202, -0.0866394, -0.0690735, -0.055337, -0.0445586, -0.0360791, -0.0293944, -0.0241163, -0.0199434, -0.0166412, -0.0140259, -0.0119533, -0.0103101, -0.00900678, -0.00797271, -0.00715208, -0.00650071, -0.00598361, -0.00557305, -0.00524706, -0.00498818, -0.0047826, -0.00461933, -0.00448966, -0.00438667, -0.00430487, -0.0042399, -0.0041883, -0.00414731, -0.00411475, -0.00408889, -0.00406835, -0.00405203, -0.00403907, -0.00402878, -0.0040206, -0.0040141},
            { -6.85508, -1.58057, -0.998204, -0.69763, -0.510144, -0.382809, -0.292081, -0.225455, -0.175525, -0.137567, -0.10841, -0.0858414, -0.0682703, -0.0545298, -0.0437483, -0.0352663, -0.0285798, -0.0233001, -0.0191261, -0.015823, -0.0132069, -0.0111338, -0.00949014, -0.00818644, -0.00715208, -0.00633122, -0.00567967, -0.00516243, -0.00475176, -0.00442567, -0.00416673, -0.00396109, -0.00379778, -0.00366807, -0.00356505, -0.00348323, -0.00341824, -0.00336662, -0.00332562, -0.00329306, -0.00326719, -0.00324664, -0.00323032, -0.00321736, -0.00320706, -0.00319888, -0.00319238},
            { -7.08533, -1.58075, -0.997924, -0.697196, -0.509636, -0.382257, -0.291502, -0.224857, -0.174914, -0.136946, -0.107782, -0.0852079, -0.0676327, -0.053889, -0.0431051, -0.0346212, -0.0279331, -0.0226523, -0.0184774, -0.0151735, -0.0125569, -0.0104833, -0.00883928, -0.00753529, -0.00650071, -0.00567967, -0.00502798, -0.00451062, -0.00409986, -0.00377371, -0.00351471, -0.00330902, -0.00314567, -0.00301594, -0.0029129, -0.00283106, -0.00276606, -0.00271443, -0.00267342, -0.00264084, -0.00261497, -0.00259442, -0.00257809, -0.00256512, -0.00255482, -0.00254664, -0.00254014},
            { -7.31559, -1.5809, -0.997702, -0.696852, -0.509232, -0.38182, -0.291042, -0.224383, -0.174428, -0.136453, -0.107284, -0.0847051, -0.0671265, -0.0533804, -0.0425945, -0.0341091, -0.0274198, -0.0221381, -0.0179624, -0.0146579, -0.0120409, -0.00996692, -0.00832259, -0.00701837, -0.00598361, -0.00516243, -0.00451062, -0.00399318, -0.00358235, -0.00325613, -0.00299709, -0.00279137, -0.00262799, -0.00249823, -0.00239518, -0.00231332, -0.00224831, -0.00219667, -0.00215565, -0.00212307, -0.00209719, -0.00207664, -0.00206031, -0.00204734, -0.00203704, -0.00202886, -0.00202236},
            { -7.54585, -1.58102, -0.997525, -0.696579, -0.508912, -0.381472, -0.290677, -0.224006, -0.174043, -0.136062, -0.106888, -0.0843058, -0.0667247, -0.0529765, -0.0421891, -0.0337024, -0.0270122, -0.0217297, -0.0175535, -0.0142486, -0.0116311, -0.00955691, -0.00791235, -0.00660796, -0.00557305, -0.00475176, -0.00409986, -0.00358235, -0.00317146, -0.0028452, -0.00258612, -0.00238037, -0.00221697, -0.0020872, -0.00198413, -0.00190226, -0.00183724, -0.00178559, -0.00174457, -0.00171198, -0.0016861, -0.00166554, -0.00164921, -0.00163624, -0.00162594, -0.00161776, -0.00161126},
            { -7.77611, -1.58111, -0.997385, -0.696362, -0.508658, -0.381196, -0.290387, -0.223707, -0.173737, -0.135751, -0.106574, -0.0839888, -0.0664056, -0.0526558, -0.0418671, -0.0333796, -0.0266886, -0.0214055, -0.0172288, -0.0139235, -0.0113058, -0.00923135, -0.00758661, -0.00628208, -0.00524706, -0.00442567, -0.00377371, -0.00325613, -0.0028452, -0.00251891, -0.0022598, -0.00205403, -0.00189061, -0.00176082, -0.00165774, -0.00157586, -0.00151083, -0.00145918, -0.00141815, -0.00138557, -0.00135968, -0.00133912, -0.00132279, -0.00130982, -0.00129951, -0.00129133, -0.00128483},
            { -8.00637, -1.58119, -0.997273, -0.69619, -0.508456, -0.380977, -0.290157, -0.223469, -0.173494, -0.135504, -0.106324, -0.083737, -0.0661522, -0.0524012, -0.0416115, -0.0331232, -0.0264316, -0.0211481, -0.016971, -0.0136654, -0.0110475, -0.00897283, -0.00732794, -0.00602329, -0.00498818, -0.00416673, -0.00351471, -0.00299709, -0.00258612, -0.0022598, -0.00200067, -0.00179488, -0.00163145, -0.00150165, -0.00139855, -0.00131667, -0.00125164, -0.00119998, -0.00115895, -0.00112636, -0.00110047, -0.00107991, -0.00106358, -0.0010506, -0.0010403, -0.00103211, -0.00102561},
            { -8.23663, -1.58125, -0.997185, -0.696053, -0.508295, -0.380803, -0.289974, -0.22328, -0.173301, -0.135308, -0.106126, -0.0835371, -0.065951, -0.0521989, -0.0414085, -0.0329196, -0.0262275, -0.0209436, -0.0167662, -0.0134604, -0.0108423, -0.00876752, -0.00712252, -0.00581778, -0.0047826, -0.00396109, -0.00330902, -0.00279137, -0.00238037, -0.00205403, -0.00179488, -0.00158908, -0.00142563, -0.00129582, -0.00119272, -0.00111084, -0.0010458, -0.000994137, -0.000953104, -0.000920511, -0.000894622, -0.000874059, -0.000857725, -0.000844751, -0.000834445, -0.000826259, -0.000819756},
            { -8.46688, -1.58129, -0.997114, -0.695944, -0.508168, -0.380664, -0.289829, -0.22313, -0.173148, -0.135153, -0.105968, -0.0833783, -0.0657912, -0.0520383, -0.0412473, -0.0327579, -0.0260655, -0.0207812, -0.0166036, -0.0132976, -0.0106794, -0.00860447, -0.00695938, -0.00565457, -0.00461933, -0.00379778, -0.00314567, -0.00262799, -0.00221697, -0.00189061, -0.00163145, -0.00142563, -0.00126218, -0.00113236, -0.00102926, -0.000947368, -0.000882324, -0.000830661, -0.000789625, -0.00075703, -0.00073114, -0.000710576, -0.000694241, -0.000681266, -0.00067096, -0.000662773, -0.00065627},
            { -8.69714, -1.58133, -0.997059, -0.695858, -0.508067, -0.380554, -0.289713, -0.223011, -0.173026, -0.135029, -0.105843, -0.0832522, -0.0656642, -0.0519108, -0.0411192, -0.0326294, -0.0259367, -0.0206523, -0.0164745, -0.0131684, -0.01055, -0.00847497, -0.00682981, -0.00552494, -0.00448966, -0.00366807, -0.00301594, -0.00249823, -0.0020872, -0.00176082, -0.00150165, -0.00129582, -0.00113236, -0.00100254, -0.000899433, -0.000817538, -0.000752491, -0.000700826, -0.000659788, -0.000627192, -0.000601301, -0.000580736, -0.0005644, -0.000551424, -0.000541118, -0.000532931, -0.000526428},
            { -8.9274, -1.58136, -0.997014, -0.695789, -0.507986, -0.380467, -0.289622, -0.222917, -0.17293, -0.134931, -0.105744, -0.083152, -0.0655634, -0.0518095, -0.0410175, -0.0325274, -0.0258345, -0.0205498, -0.0163719, -0.0130657, -0.0104472, -0.00837212, -0.00672691, -0.00542199, -0.00438667, -0.00356505, -0.0029129, -0.00239518, -0.00198413, -0.00165774, -0.00139855, -0.00119272, -0.00102926, -0.000899433, -0.00079632, -0.000714422, -0.000649373, -0.000597706, -0.000556667, -0.00052407, -0.000498178, -0.000477612, -0.000461276, -0.0004483, -0.000437993, -0.000429806, -0.000423302},
            { -9.15766, -1.58138, -0.996979, -0.695735, -0.507922, -0.380398, -0.289549, -0.222842, -0.172853, -0.134853, -0.105665, -0.0830725, -0.0654833, -0.051729, -0.0409368, -0.0324464, -0.0257533, -0.0204685, -0.0162904, -0.0129841, -0.0103655, -0.00829043, -0.00664517, -0.00534022, -0.00430487, -0.00348323, -0.00283106, -0.00231332, -0.00190226, -0.00157586, -0.00131667, -0.00111084, -0.000947368, -0.000817538, -0.000714422, -0.000632522, -0.000567471, -0.000515803, -0.000474763, -0.000442165, -0.000416272, -0.000395705, -0.000379369, -0.000366392, -0.000356085, -0.000347898, -0.000341394},
            { -9.38792, -1.5814, -0.996951, -0.695692, -0.507872, -0.380343, -0.289491, -0.222782, -0.172792, -0.134791, -0.105602, -0.0830093, -0.0654198, -0.0516651, -0.0408726, -0.0323821, -0.0256888, -0.0204039, -0.0162257, -0.0129193, -0.0103007, -0.00822555, -0.00658025, -0.00527527, -0.0042399, -0.00341824, -0.00276606, -0.00224831, -0.00183724, -0.00151083, -0.00125164, -0.0010458, -0.000882324, -0.000752491, -0.000649373, -0.000567471, -0.000502419, -0.00045075, -0.000409709, -0.00037711, -0.000351217, -0.00033065, -0.000314313, -0.000301336, -0.000291028, -0.000282841, -0.000276337},
            { -9.61818, -1.58142, -0.996929, -0.695657, -0.507831, -0.380299, -0.289445, -0.222734, -0.172744, -0.134742, -0.105553, -0.0829591, -0.0653692, -0.0516143, -0.0408216, -0.0323309, -0.0256376, -0.0203526, -0.0161743, -0.0128679, -0.0102492, -0.00817401, -0.00652869, -0.00522368, -0.0041883, -0.00336662, -0.00271443, -0.00219667, -0.00178559, -0.00145918, -0.00119998, -0.000994137, -0.000830661, -0.000700826, -0.000597706, -0.000515803, -0.00045075, -0.000399079, -0.000358037, -0.000325438, -0.000299544, -0.000278977, -0.00026264, -0.000249663, -0.000239355, -0.000231167, -0.000224664},
            { -9.84844, -1.58143, -0.996911, -0.69563, -0.507799, -0.380264, -0.289409, -0.222697, -0.172705, -0.134703, -0.105513, -0.0829192, -0.0653291, -0.051574, -0.0407812, -0.0322904, -0.0255969, -0.0203118, -0.0161335, -0.012827, -0.0102083, -0.00813308, -0.00648773, -0.00518271, -0.00414731, -0.00332562, -0.00267342, -0.00215565, -0.00174457, -0.00141815, -0.00115895, -0.000953104, -0.000789625, -0.000659788, -0.000556667, -0.000474763, -0.000409709, -0.000358037, -0.000316995, -0.000284396, -0.000258502, -0.000237934, -0.000221596, -0.000208619, -0.000198311, -0.000190123, -0.00018362},
            { -10.0787, -1.58144, -0.996897, -0.695608, -0.507774, -0.380237, -0.28938, -0.222667, -0.172675, -0.134672, -0.105482, -0.0828876, -0.0652972, -0.051542, -0.040749, -0.0322581, -0.0255645, -0.0202794, -0.0161011, -0.0127945, -0.0101758, -0.00810056, -0.0064552, -0.00515016, -0.00411475, -0.00329306, -0.00264084, -0.00212307, -0.00171198, -0.00138557, -0.00112636, -0.000920511, -0.00075703, -0.000627192, -0.00052407, -0.000442165, -0.00037711, -0.000325438, -0.000284396, -0.000251796, -0.000225901, -0.000205333, -0.000188996, -0.000176018, -0.00016571, -0.000157522, -0.000151019},
            { -10.309, -1.58145, -0.996886, -0.695591, -0.507754, -0.380215, -0.289357, -0.222643, -0.17265, -0.134647, -0.105457, -0.0828624, -0.0652719, -0.0515165, -0.0407235, -0.0322325, -0.0255389, -0.0202537, -0.0160753, -0.0127688, -0.01015, -0.00807474, -0.00642936, -0.00512431, -0.00408889, -0.00326719, -0.00261497, -0.00209719, -0.0016861, -0.00135968, -0.00110047, -0.000894622, -0.00073114, -0.000601301, -0.000498178, -0.000416272, -0.000351217, -0.000299544, -0.000258502, -0.000225901, -0.000200007, -0.000179438, -0.000163101, -0.000150123, -0.000139815, -0.000131627, -0.000125123},
            { -10.5392, -1.58145, -0.996877, -0.695577, -0.507738, -0.380198, -0.289339, -0.222624, -0.172631, -0.134628, -0.105437, -0.0828425, -0.0652518, -0.0514963, -0.0407032, -0.0322121, -0.0255185, -0.0202333, -0.0160549, -0.0127483, -0.0101295, -0.00805422, -0.00640883, -0.00510378, -0.00406835, -0.00324664, -0.00259442, -0.00207664, -0.00166554, -0.00133912, -0.00107991, -0.000874059, -0.000710576, -0.000580736, -0.000477612, -0.000395705, -0.00033065, -0.000278977, -0.000237934, -0.000205333, -0.000179438, -0.00015887, -0.000142532, -0.000129555, -0.000119246, -0.000111058, -0.000104554},
            { -10.7695, -1.58146, -0.99687, -0.695566, -0.507725, -0.380184, -0.289324, -0.222609, -0.172616, -0.134612, -0.105421, -0.0828266, -0.0652359, -0.0514803, -0.0406871, -0.032196, -0.0255023, -0.020217, -0.0160386, -0.012732, -0.0101132, -0.00803793, -0.00639253, -0.00508747, -0.00405203, -0.00323032, -0.00257809, -0.00206031, -0.00164921, -0.00132279, -0.00106358, -0.000857725, -0.000694241, -0.0005644, -0.000461276, -0.000379369, -0.000314313, -0.00026264, -0.000221596, -0.000188996, -0.000163101, -0.000142532, -0.000126194, -0.000113217, -0.000102908, -9.47203e-05, -8.82164e-05},
            { -10.9997, -1.58146, -0.996865, -0.695558, -0.507715, -0.380173, -0.289313, -0.222597, -0.172604, -0.1346, -0.105409, -0.082814, -0.0652232, -0.0514675, -0.0406743, -0.0321831, -0.0254894, -0.0202041, -0.0160257, -0.0127191, -0.0101003, -0.00802498, -0.00637958, -0.00507451, -0.00403907, -0.00321736, -0.00256512, -0.00204734, -0.00163624, -0.00130982, -0.0010506, -0.000844751, -0.000681266, -0.000551424, -0.0004483, -0.000366392, -0.000301336, -0.000249663, -0.000208619, -0.000176018, -0.000150123, -0.000129555, -0.000113217, -0.000100239, -8.99308e-05, -8.17427e-05, -7.52387e-05},
            { -11.23, -1.58146, -0.99686, -0.695551, -0.507707, -0.380164, -0.289304, -0.222588, -0.172594, -0.13459, -0.105399, -0.082804, -0.0652131, -0.0514574, -0.0406641, -0.032173, -0.0254792, -0.0201939, -0.0160155, -0.0127088, -0.01009, -0.0080147, -0.00636929, -0.00506422, -0.00402878, -0.00320706, -0.00255482, -0.00203704, -0.00162594, -0.00129951, -0.0010403, -0.000834445, -0.00067096, -0.000541118, -0.000437993, -0.000356085, -0.000291028, -0.000239355, -0.000198311, -0.00016571, -0.000139815, -0.000119246, -0.000102908, -8.99308e-05, -7.96225e-05, -7.14344e-05, -6.49304e-05},
            { -11.4602, -1.58146, -0.996857, -0.695546, -0.507701, -0.380157, -0.289296, -0.222581, -0.172586, -0.134582, -0.105391, -0.082796, -0.0652051, -0.0514493, -0.040656, -0.0321649, -0.0254711, -0.0201858, -0.0160073, -0.0127007, -0.0100819, -0.00800654, -0.00636112, -0.00505604, -0.0040206, -0.00319888, -0.00254664, -0.00202886, -0.00161776, -0.00129133, -0.00103211, -0.000826259, -0.000662773, -0.000532931, -0.000429806, -0.000347898, -0.000282841, -0.000231167, -0.000190123, -0.000157522, -0.000131627, -0.000111058, -9.47203e-05, -8.17427e-05, -7.14344e-05, -6.32462e-05, -5.67422e-05},
            { -11.6905, -1.58147, -0.996854, -0.695541, -0.507695, -0.380152, -0.28929, -0.222575, -0.17258, -0.134576, -0.105385, -0.0827897, -0.0651987, -0.051443, -0.0406496, -0.0321584, -0.0254646, -0.0201793, -0.0160009, -0.0126942, -0.0100754, -0.00800005, -0.00635463, -0.00504955, -0.0040141, -0.00319238, -0.00254014, -0.00202236, -0.00161126, -0.00128483, -0.00102561, -0.000819756, -0.00065627, -0.000526428, -0.000423302, -0.000341394, -0.000276337, -0.000224664, -0.00018362, -0.000151019, -0.000125123, -0.000104554, -8.82164e-05, -7.52387e-05, -6.49304e-05, -5.67422e-05, -5.02381e-05}};
        const double qual_mismatch_simple_bayesian[][47] = {
            { -1.50408, -1.40619, -1.33474, -1.28141, -1.24099, -1.21, -1.18606, -1.16744, -1.15289, -1.14148, -1.13251, -1.12545, -1.11987, -1.11546, -1.11197, -1.10921, -1.10702, -1.10529, -1.10391, -1.10282, -1.10195, -1.10126, -1.10072, -1.10028, -1.09994, -1.09967, -1.09945, -1.09928, -1.09914, -1.09903, -1.09895, -1.09888, -1.09882, -1.09878, -1.09874, -1.09872, -1.0987, -1.09868, -1.09867, -1.09865, -1.09865, -1.09864, -1.09863, -1.09863, -1.09863, -1.09862, -1.09862},
            { -1.40619, -1.38979, -1.37696, -1.36688, -1.35894, -1.35268, -1.34774, -1.34383, -1.34073, -1.33828, -1.33634, -1.3348, -1.33358, -1.33261, -1.33184, -1.33123, -1.33074, -1.33036, -1.33005, -1.32981, -1.32962, -1.32946, -1.32934, -1.32924, -1.32917, -1.32911, -1.32906, -1.32902, -1.32899, -1.32896, -1.32895, -1.32893, -1.32892, -1.32891, -1.3289, -1.32889, -1.32889, -1.32889, -1.32888, -1.32888, -1.32888, -1.32888, -1.32888, -1.32887, -1.32887, -1.32887, -1.32887},
            { -1.33474, -1.37696, -1.41181, -1.44039, -1.46368, -1.48258, -1.49786, -1.51016, -1.52003, -1.52795, -1.53428, -1.53934, -1.54338, -1.5466, -1.54916, -1.55121, -1.55283, -1.55412, -1.55515, -1.55597, -1.55662, -1.55713, -1.55754, -1.55787, -1.55813, -1.55833, -1.5585, -1.55863, -1.55873, -1.55881, -1.55888, -1.55893, -1.55897, -1.559, -1.55903, -1.55905, -1.55907, -1.55908, -1.55909, -1.5591, -1.5591, -1.55911, -1.55911, -1.55912, -1.55912, -1.55912, -1.55912},
            { -1.28141, -1.36688, -1.44039, -1.50289, -1.55549, -1.59933, -1.63558, -1.66534, -1.68963, -1.70935, -1.72529, -1.73814, -1.74847, -1.75675, -1.76338, -1.76867, -1.7729, -1.77627, -1.77895, -1.78109, -1.78279, -1.78414, -1.78522, -1.78608, -1.78676, -1.7873, -1.78773, -1.78807, -1.78834, -1.78855, -1.78873, -1.78886, -1.78897, -1.78906, -1.78912, -1.78918, -1.78922, -1.78926, -1.78928, -1.7893, -1.78932, -1.78934, -1.78935, -1.78935, -1.78936, -1.78937, -1.78937},
            { -1.24099, -1.35894, -1.46368, -1.55549, -1.63493, -1.70287, -1.76033, -1.80845, -1.8484, -1.8813, -1.90823, -1.93016, -1.94792, -1.96226, -1.97379, -1.98305, -1.99047, -1.9964, -2.00114, -2.00492, -2.00793, -2.01033, -2.01224, -2.01376, -2.01497, -2.01593, -2.01669, -2.0173, -2.01778, -2.01816, -2.01847, -2.01871, -2.0189, -2.01906, -2.01918, -2.01927, -2.01935, -2.01941, -2.01946, -2.0195, -2.01953, -2.01955, -2.01957, -2.01959, -2.0196, -2.01961, -2.01962},
            { -1.21, -1.35268, -1.48258, -1.59933, -1.70287, -1.79352, -1.87187, -1.93881, -1.99536, -2.04269, -2.08194, -2.11426, -2.14069, -2.1622, -2.17962, -2.19368, -2.20499, -2.21406, -2.22133, -2.22714, -2.23178, -2.23548, -2.23843, -2.24078, -2.24265, -2.24414, -2.24532, -2.24626, -2.24701, -2.2476, -2.24808, -2.24845, -2.24875, -2.24899, -2.24918, -2.24933, -2.24945, -2.24954, -2.24962, -2.24967, -2.24972, -2.24976, -2.24979, -2.24981, -2.24983, -2.24985, -2.24986},
            { -1.18606, -1.34774, -1.49786, -1.63558, -1.76033, -1.87187, -1.97029, -2.05601, -2.12976, -2.19248, -2.24527, -2.28928, -2.32567, -2.35556, -2.37995, -2.39976, -2.41577, -2.42868, -2.43906, -2.44737, -2.45403, -2.45935, -2.4636, -2.46698, -2.46968, -2.47183, -2.47353, -2.47489, -2.47598, -2.47684, -2.47752, -2.47806, -2.47849, -2.47884, -2.47911, -2.47933, -2.4795, -2.47964, -2.47974, -2.47983, -2.4799, -2.47995, -2.48, -2.48003, -2.48006, -2.48008, -2.4801},
            { -1.16744, -1.34383, -1.51016, -1.66534, -1.80845, -1.93881, -2.05601, -2.16001, -2.25109, -2.32986, -2.39718, -2.45408, -2.5017, -2.54122, -2.57376, -2.60038, -2.62204, -2.63959, -2.65376, -2.66515, -2.6743, -2.68162, -2.68748, -2.69215, -2.69588, -2.69886, -2.70122, -2.70311, -2.70461, -2.7058, -2.70675, -2.7075, -2.7081, -2.70858, -2.70896, -2.70926, -2.7095, -2.70969, -2.70984, -2.70996, -2.71005, -2.71013, -2.71019, -2.71024, -2.71028, -2.71031, -2.71033},
            { -1.15289, -1.34073, -1.52003, -1.68963, -1.8484, -1.99536, -2.12976, -2.25109, -2.3592, -2.45427, -2.5368, -2.60759, -2.66762, -2.71801, -2.75994, -2.79454, -2.8229, -2.84602, -2.86477, -2.87992, -2.89212, -2.90191, -2.90977, -2.91605, -2.92106, -2.92507, -2.92826, -2.9308, -2.93282, -2.93444, -2.93572, -2.93674, -2.93755, -2.93819, -2.9387, -2.93911, -2.93943, -2.93969, -2.93989, -2.94005, -2.94018, -2.94029, -2.94037, -2.94043, -2.94048, -2.94052, -2.94056},
            { -1.14148, -1.33828, -1.52795, -1.70935, -1.8813, -2.04269, -2.19248, -2.32986, -2.45427, -2.56545, -2.66352, -2.74891, -2.82235, -2.8848, -2.93733, -2.98112, -3.01733, -3.04705, -3.07131, -3.09101, -3.10693, -3.11977, -3.13008, -3.13835, -3.14496, -3.15025, -3.15447, -3.15784, -3.16052, -3.16265, -3.16435, -3.1657, -3.16678, -3.16763, -3.16831, -3.16885, -3.16928, -3.16962, -3.16989, -3.17011, -3.17028, -3.17041, -3.17052, -3.17061, -3.17068, -3.17073, -3.17077},
            { -1.13251, -1.33634, -1.53428, -1.72529, -1.90823, -2.08194, -2.24527, -2.39718, -2.5368, -2.66352, -2.77704, -2.87741, -2.96499, -3.04048, -3.10478, -3.15899, -3.20424, -3.2417, -3.27249, -3.29764, -3.31808, -3.33462, -3.34796, -3.35868, -3.36728, -3.37416, -3.37966, -3.38405, -3.38756, -3.39035, -3.39257, -3.39434, -3.39574, -3.39686, -3.39775, -3.39846, -3.39902, -3.39947, -3.39982, -3.40011, -3.40033, -3.40051, -3.40065, -3.40076, -3.40085, -3.40092, -3.40098},
            { -1.12545, -1.3348, -1.53934, -1.73814, -1.93016, -2.11426, -2.28928, -2.45408, -2.60759, -2.74891, -2.87741, -2.99272, -3.09485, -3.18412, -3.2612, -3.32696, -3.38246, -3.42885, -3.4673, -3.49893, -3.52479, -3.54582, -3.56284, -3.57658, -3.58762, -3.59648, -3.60357, -3.60925, -3.61377, -3.61738, -3.62026, -3.62255, -3.62438, -3.62583, -3.62698, -3.6279, -3.62863, -3.62921, -3.62967, -3.63004, -3.63033, -3.63056, -3.63075, -3.63089, -3.63101, -3.6311, -3.63117},
            { -1.11987, -1.33358, -1.54338, -1.74847, -1.94792, -2.14069, -2.32567, -2.5017, -2.66762, -2.82235, -2.96499, -3.09485, -3.21154, -3.31504, -3.40563, -3.48395, -3.55084, -3.60736, -3.65465, -3.69388, -3.72617, -3.75259, -3.77408, -3.79149, -3.80553, -3.81683, -3.8259, -3.83316, -3.83897, -3.84361, -3.8473, -3.85025, -3.8526, -3.85447, -3.85595, -3.85713, -3.85807, -3.85882, -3.85942, -3.85989, -3.86026, -3.86056, -3.8608, -3.86099, -3.86114, -3.86126, -3.86135},
            { -1.11546, -1.33261, -1.5466, -1.75675, -1.96226, -2.1622, -2.35556, -2.54122, -2.71801, -2.8848, -3.04048, -3.18412, -3.31504, -3.43281, -3.53737, -3.629, -3.70828, -3.77607, -3.83339, -3.88139, -3.92122, -3.95404, -3.9809, -4.00276, -4.02047, -4.03476, -4.04626, -4.0555, -4.06289, -4.0688, -4.07352, -4.07729, -4.08029, -4.08268, -4.08459, -4.0861, -4.08731, -4.08826, -4.08903, -4.08963, -4.09011, -4.0905, -4.0908, -4.09104, -4.09123, -4.09138, -4.09151},
            { -1.11197, -1.33184, -1.54916, -1.76338, -1.97379, -2.17962, -2.37995, -2.57376, -2.75994, -2.93733, -3.10478, -3.2612, -3.40563, -3.53737, -3.65598, -3.76138, -3.85381, -3.93386, -4.00234, -4.0603, -4.10885, -4.14917, -4.1824, -4.20961, -4.23176, -4.24971, -4.2642, -4.27586, -4.28523, -4.29273, -4.29872, -4.30351, -4.30734, -4.31038, -4.31281, -4.31474, -4.31627, -4.3175, -4.31847, -4.31924, -4.31986, -4.32034, -4.32073, -4.32104, -4.32128, -4.32148, -4.32163},
            { -1.10921, -1.33123, -1.55121, -1.76867, -1.98305, -2.19368, -2.39976, -2.60038, -2.79454, -2.98112, -3.15899, -3.32696, -3.48395, -3.629, -3.76138, -3.88065, -3.9867, -4.07977, -4.16041, -4.22945, -4.2879, -4.3369, -4.3776, -4.41116, -4.43864, -4.46102, -4.47916, -4.49381, -4.5056, -4.51507, -4.52265, -4.52872, -4.53356, -4.53742, -4.5405, -4.54296, -4.54491, -4.54646, -4.5477, -4.54868, -4.54947, -4.55009, -4.55058, -4.55097, -4.55128, -4.55153, -4.55173},
            { -1.10702, -1.33074, -1.55283, -1.7729, -1.99047, -2.20499, -2.41577, -2.62204, -2.8229, -3.01733, -3.20424, -3.38246, -3.55084, -3.70828, -3.85381, -3.9867, -4.10649, -4.21306, -4.30662, -4.38774, -4.45721, -4.51606, -4.5654, -4.60641, -4.64022, -4.66792, -4.69049, -4.70878, -4.72355, -4.73544, -4.74499, -4.75264, -4.75876, -4.76365, -4.76755, -4.77065, -4.77313, -4.7751, -4.77667, -4.77792, -4.77891, -4.7797, -4.78032, -4.78082, -4.78122, -4.78153, -4.78178},
            { -1.10529, -1.33036, -1.55412, -1.77627, -1.9964, -2.21406, -2.42868, -2.63959, -2.84602, -3.04705, -3.2417, -3.42885, -3.60736, -3.77607, -3.93386, -4.07977, -4.21306, -4.33325, -4.44022, -4.53419, -4.61567, -4.68549, -4.74465, -4.79427, -4.83552, -4.86954, -4.89741, -4.92012, -4.93853, -4.9534, -4.96537, -4.97499, -4.98269, -4.98885, -4.99377, -4.9977, -5.00083, -5.00332, -5.0053, -5.00688, -5.00814, -5.00914, -5.00993, -5.01056, -5.01107, -5.01147, -5.01178},
            { -1.10391, -1.33005, -1.55515, -1.77895, -2.00114, -2.22133, -2.43906, -2.65376, -2.86477, -3.07131, -3.27249, -3.4673, -3.65465, -3.83339, -4.00234, -4.16041, -4.30662, -4.44022, -4.56074, -4.66803, -4.76231, -4.84409, -4.91418, -4.97359, -5.02342, -5.06486, -5.09904, -5.12706, -5.14988, -5.16839, -5.18334, -5.19537, -5.20504, -5.21278, -5.21897, -5.22392, -5.22787, -5.23102, -5.23352, -5.23552, -5.23711, -5.23837, -5.23938, -5.24017, -5.24081, -5.24131, -5.24172},
            { -1.10282, -1.32981, -1.55597, -1.78109, -2.00492, -2.22714, -2.44737, -2.66515, -2.87992, -3.09101, -3.29764, -3.49893, -3.69388, -3.88139, -4.0603, -4.22945, -4.38774, -4.53419, -4.66803, -4.78881, -4.89635, -4.99087, -5.07289, -5.1432, -5.2028, -5.25281, -5.29439, -5.32871, -5.35683, -5.37974, -5.39832, -5.41334, -5.42542, -5.43513, -5.44291, -5.44913, -5.4541, -5.45806, -5.46122, -5.46374, -5.46574, -5.46734, -5.46861, -5.46962, -5.47042, -5.47106, -5.47156},
            { -1.10195, -1.32962, -1.55662, -1.78279, -2.00793, -2.23178, -2.45403, -2.6743, -2.89212, -3.10693, -3.31808, -3.52479, -3.72617, -3.92122, -4.10885, -4.2879, -4.45721, -4.61567, -4.76231, -4.89635, -5.01732, -5.12507, -5.21979, -5.30199, -5.37247, -5.43222, -5.48237, -5.52408, -5.55849, -5.5867, -5.60969, -5.62833, -5.64339, -5.65552, -5.66525, -5.67306, -5.6793, -5.68429, -5.68827, -5.69144, -5.69396, -5.69598, -5.69758, -5.69885, -5.69986, -5.70067, -5.70131},
            { -1.10126, -1.32946, -1.55713, -1.78414, -2.01033, -2.23548, -2.45935, -2.68162, -2.90191, -3.11977, -3.33462, -3.54582, -3.75259, -3.95404, -4.14917, -4.3369, -4.51606, -4.68549, -4.84409, -4.99087, -5.12507, -5.2462, -5.35411, -5.44898, -5.53133, -5.60194, -5.66182, -5.71208, -5.75388, -5.78837, -5.81665, -5.83969, -5.85838, -5.87348, -5.88564, -5.89541, -5.90323, -5.90949, -5.91449, -5.91848, -5.92166, -5.9242, -5.92621, -5.92782, -5.92909, -5.93011, -5.93092},
            { -1.10072, -1.32934, -1.55754, -1.78522, -2.01224, -2.23843, -2.4636, -2.68748, -2.90977, -3.13008, -3.34796, -3.56284, -3.77408, -3.9809, -4.1824, -4.3776, -4.5654, -4.74465, -4.91418, -5.07289, -5.21979, -5.35411, -5.47537, -5.5834, -5.67839, -5.76086, -5.83158, -5.89155, -5.9419, -5.98377, -6.01833, -6.04666, -6.06975, -6.08848, -6.10361, -6.1158, -6.12558, -6.13342, -6.1397, -6.14471, -6.14871, -6.15189, -6.15443, -6.15645, -6.15806, -6.15934, -6.16036},
            { -1.10028, -1.32924, -1.55787, -1.78608, -2.01376, -2.24078, -2.46698, -2.69215, -2.91605, -3.13835, -3.35868, -3.57658, -3.79149, -4.00276, -4.20961, -4.41116, -4.60641, -4.79427, -4.97359, -5.1432, -5.30199, -5.44898, -5.5834, -5.70476, -5.81289, -5.90798, -5.99054, -6.06134, -6.12139, -6.17181, -6.21374, -6.24836, -6.27673, -6.29986, -6.31861, -6.33377, -6.34597, -6.35578, -6.36363, -6.36991, -6.37493, -6.37894, -6.38213, -6.38467, -6.3867, -6.38831, -6.38959},
            { -1.09994, -1.32917, -1.55813, -1.78676, -2.01497, -2.24265, -2.46968, -2.69588, -2.92106, -3.14496, -3.36728, -3.58762, -3.80553, -4.02047, -4.23176, -4.43864, -4.64022, -4.83552, -5.02342, -5.2028, -5.37247, -5.53133, -5.67839, -5.81289, -5.93433, -6.04254, -6.1377, -6.22033, -6.29121, -6.35132, -6.40179, -6.44377, -6.47843, -6.50683, -6.52999, -6.54877, -6.56395, -6.57617, -6.58598, -6.59385, -6.60014, -6.60516, -6.60917, -6.61237, -6.61492, -6.61695, -6.61856},
            { -1.09967, -1.32911, -1.55833, -1.7873, -2.01593, -2.24414, -2.47183, -2.69886, -2.92507, -3.15025, -3.37416, -3.59648, -3.81683, -4.03476, -4.24971, -4.46102, -4.66792, -4.86954, -5.06486, -5.25281, -5.43222, -5.60194, -5.76086, -5.90798, -6.04254, -6.16404, -6.27231, -6.36754, -6.45023, -6.52116, -6.58132, -6.63183, -6.67385, -6.70854, -6.73697, -6.76015, -6.77895, -6.79414, -6.80637, -6.8162, -6.82407, -6.83037, -6.8354, -6.83942, -6.84262, -6.84517, -6.8472},
            { -1.09945, -1.32906, -1.5585, -1.78773, -2.01669, -2.24532, -2.47353, -2.70122, -2.92826, -3.15447, -3.37966, -3.60357, -3.8259, -4.04626, -4.2642, -4.47916, -4.69049, -4.89741, -5.09904, -5.29439, -5.48237, -5.66182, -5.83158, -5.99054, -6.1377, -6.27231, -6.39386, -6.50219, -6.59746, -6.6802, -6.75117, -6.81137, -6.86191, -6.90396, -6.93867, -6.96713, -6.99033, -7.00914, -7.02435, -7.03659, -7.04642, -7.0543, -7.06061, -7.06564, -7.06966, -7.07286, -7.07542},
            { -1.09928, -1.32902, -1.55863, -1.78807, -2.0173, -2.24626, -2.47489, -2.70311, -2.9308, -3.15784, -3.38405, -3.60925, -3.83316, -4.0555, -4.27586, -4.49381, -4.70878, -4.92012, -5.12706, -5.32871, -5.52408, -5.71208, -5.89155, -6.06134, -6.22033, -6.36754, -6.50219, -6.62378, -6.73214, -6.82745, -6.91022, -6.98123, -7.04146, -7.09203, -7.13411, -7.16884, -7.19731, -7.22052, -7.23935, -7.25456, -7.26682, -7.27666, -7.28454, -7.29085, -7.29589, -7.29991, -7.30311},
            { -1.09914, -1.32899, -1.55873, -1.78834, -2.01778, -2.24701, -2.47598, -2.70461, -2.93282, -3.16052, -3.38756, -3.61377, -3.83897, -4.06289, -4.28523, -4.5056, -4.72355, -4.93853, -5.14988, -5.35683, -5.55849, -5.75388, -5.9419, -6.12139, -6.29121, -6.45023, -6.59746, -6.73214, -6.85376, -6.96216, -7.0575, -7.1403, -7.21133, -7.27159, -7.32218, -7.36428, -7.39902, -7.42751, -7.45073, -7.46957, -7.48479, -7.49705, -7.50689, -7.51478, -7.52109, -7.52614, -7.53016},
            { -1.09903, -1.32896, -1.55881, -1.78855, -2.01816, -2.2476, -2.47684, -2.7058, -2.93444, -3.16265, -3.39035, -3.61738, -3.84361, -4.0688, -4.29273, -4.51507, -4.73544, -4.9534, -5.16839, -5.37974, -5.5867, -5.78837, -5.98377, -6.17181, -6.35132, -6.52116, -6.6802, -6.82745, -6.96216, -7.0838, -7.19222, -7.28759, -7.37041, -7.44147, -7.50174, -7.55235, -7.59446, -7.62922, -7.65772, -7.68095, -7.6998, -7.71502, -7.72729, -7.73713, -7.74503, -7.75134, -7.75639},
            { -1.09895, -1.32895, -1.55888, -1.78873, -2.01847, -2.24808, -2.47752, -2.70675, -2.93572, -3.16435, -3.39257, -3.62026, -3.8473, -4.07352, -4.29872, -4.52265, -4.74499, -4.96537, -5.18334, -5.39832, -5.60969, -5.81665, -6.01833, -6.21374, -6.40179, -6.58132, -6.75117, -6.91022, -7.0575, -7.19222, -7.31389, -7.42233, -7.51772, -7.60056, -7.67163, -7.73192, -7.78254, -7.82466, -7.85943, -7.88794, -7.91118, -7.93003, -7.94526, -7.95753, -7.96738, -7.97528, -7.98159},
            { -1.09888, -1.32893, -1.55893, -1.78886, -2.01871, -2.24845, -2.47806, -2.7075, -2.93674, -3.1657, -3.39434, -3.62255, -3.85025, -4.07729, -4.30351, -4.52872, -4.75264, -4.97499, -5.19537, -5.41334, -5.62833, -5.83969, -6.04666, -6.24836, -6.44377, -6.63183, -6.81137, -6.98123, -7.1403, -7.28759, -7.42233, -7.54401, -7.65246, -7.74787, -7.83072, -7.90181, -7.96211, -8.01274, -8.05488, -8.08965, -8.11817, -8.14141, -8.16027, -8.1755, -8.18777, -8.19763, -8.20553},
            { -1.09882, -1.32892, -1.55897, -1.78897, -2.0189, -2.24875, -2.47849, -2.7081, -2.93755, -3.16678, -3.39574, -3.62438, -3.8526, -4.08029, -4.30734, -4.53356, -4.75876, -4.98269, -5.20504, -5.42542, -5.64339, -5.85838, -6.06975, -6.27673, -6.47843, -6.67385, -6.86191, -7.04146, -7.21133, -7.37041, -7.51772, -7.65246, -7.77416, -7.88263, -7.97804, -8.06091, -8.132, -8.19232, -8.24296, -8.2851, -8.31988, -8.3484, -8.37165, -8.39051, -8.40575, -8.41802, -8.42788},
            { -1.09878, -1.32891, -1.559, -1.78906, -2.01906, -2.24899, -2.47884, -2.70858, -2.93819, -3.16763, -3.39686, -3.62583, -3.85447, -4.08268, -4.31038, -4.53742, -4.76365, -4.98885, -5.21278, -5.43513, -5.65552, -5.87348, -6.08848, -6.29986, -6.50683, -6.70854, -6.90396, -7.09203, -7.27159, -7.44147, -7.60056, -7.74787, -7.88263, -8.00433, -8.11281, -8.20823, -8.29111, -8.36221, -8.42253, -8.47318, -8.51533, -8.55012, -8.57864, -8.60189, -8.62076, -8.636, -8.64827},
            { -1.09874, -1.3289, -1.55903, -1.78912, -2.01918, -2.24918, -2.47911, -2.70896, -2.9387, -3.16831, -3.39775, -3.62698, -3.85595, -4.08459, -4.31281, -4.5405, -4.76755, -4.99377, -5.21897, -5.44291, -5.66525, -5.88564, -6.10361, -6.31861, -6.52999, -6.73697, -6.93867, -7.13411, -7.32218, -7.50174, -7.67163, -7.83072, -7.97804, -8.11281, -8.23452, -8.34301, -8.43844, -8.52132, -8.59243, -8.65276, -8.70341, -8.74556, -8.78036, -8.80888, -8.83214, -8.851, -8.86625},
            { -1.09872, -1.32889, -1.55905, -1.78918, -2.01927, -2.24933, -2.47933, -2.70926, -2.93911, -3.16885, -3.39846, -3.6279, -3.85713, -4.0861, -4.31474, -4.54296, -4.77065, -4.9977, -5.22392, -5.44913, -5.67306, -5.89541, -6.1158, -6.33377, -6.54877, -6.76015, -6.96713, -7.16884, -7.36428, -7.55235, -7.73192, -7.90181, -8.06091, -8.20823, -8.34301, -8.46472, -8.57322, -8.66866, -8.75154, -8.82266, -8.88299, -8.93365, -8.9758, -9.0106, -9.03913, -9.06239, -9.08126},
            { -1.0987, -1.32889, -1.55907, -1.78922, -2.01935, -2.24945, -2.4795, -2.7095, -2.93943, -3.16928, -3.39902, -3.62863, -3.85807, -4.08731, -4.31627, -4.54491, -4.77313, -5.00083, -5.22787, -5.4541, -5.6793, -5.90323, -6.12558, -6.34597, -6.56395, -6.77895, -6.99033, -7.19731, -7.39902, -7.59446, -7.78254, -7.96211, -8.132, -8.29111, -8.43844, -8.57322, -8.69494, -8.80344, -8.89888, -8.98177, -9.05289, -9.11323, -9.16389, -9.20605, -9.24085, -9.26938, -9.29264},
            { -1.09868, -1.32889, -1.55908, -1.78926, -2.01941, -2.24954, -2.47964, -2.70969, -2.93969, -3.16962, -3.39947, -3.62921, -3.85882, -4.08826, -4.3175, -4.54646, -4.7751, -5.00332, -5.23102, -5.45806, -5.68429, -5.90949, -6.13342, -6.35578, -6.57617, -6.79414, -7.00914, -7.22052, -7.42751, -7.62922, -7.82466, -8.01274, -8.19232, -8.36221, -8.52132, -8.66866, -8.80344, -8.92516, -9.03366, -9.12911, -9.21201, -9.28313, -9.34347, -9.39414, -9.43629, -9.4711, -9.49963},
            { -1.09867, -1.32888, -1.55909, -1.78928, -2.01946, -2.24962, -2.47974, -2.70984, -2.93989, -3.16989, -3.39982, -3.62967, -3.85942, -4.08903, -4.31847, -4.5477, -4.77667, -5.0053, -5.23352, -5.46122, -5.68827, -5.91449, -6.1397, -6.36363, -6.58598, -6.80637, -7.02435, -7.23935, -7.45073, -7.65772, -7.85943, -8.05488, -8.24296, -8.42253, -8.59243, -8.75154, -8.89888, -9.03366, -9.15539, -9.2639, -9.35935, -9.44225, -9.51338, -9.57372, -9.62438, -9.66654, -9.70135},
            { -1.09865, -1.32888, -1.5591, -1.7893, -2.0195, -2.24967, -2.47983, -2.70996, -2.94005, -3.17011, -3.40011, -3.63004, -3.85989, -4.08963, -4.31924, -4.54868, -4.77792, -5.00688, -5.23552, -5.46374, -5.69144, -5.91848, -6.14471, -6.36991, -6.59385, -6.8162, -7.03659, -7.25456, -7.46957, -7.68095, -7.88794, -8.08965, -8.2851, -8.47318, -8.65276, -8.82266, -8.98177, -9.12911, -9.2639, -9.38563, -9.49414, -9.58959, -9.67249, -9.74362, -9.80396, -9.85463, -9.8968},
            { -1.09865, -1.32888, -1.5591, -1.78932, -2.01953, -2.24972, -2.4799, -2.71005, -2.94018, -3.17028, -3.40033, -3.63033, -3.86026, -4.09011, -4.31986, -4.54947, -4.77891, -5.00814, -5.23711, -5.46574, -5.69396, -5.92166, -6.14871, -6.37493, -6.60014, -6.82407, -7.04642, -7.26682, -7.48479, -7.6998, -7.91118, -8.11817, -8.31988, -8.51533, -8.70341, -8.88299, -9.05289, -9.21201, -9.35935, -9.49414, -9.61587, -9.72438, -9.81984, -9.90274, -9.97387, -10.0342, -10.0849},
            { -1.09864, -1.32888, -1.55911, -1.78934, -2.01955, -2.24976, -2.47995, -2.71013, -2.94029, -3.17041, -3.40051, -3.63056, -3.86056, -4.0905, -4.32034, -4.55009, -4.7797, -5.00914, -5.23837, -5.46734, -5.69598, -5.9242, -6.15189, -6.37894, -6.60516, -6.83037, -7.0543, -7.27666, -7.49705, -7.71502, -7.93003, -8.14141, -8.3484, -8.55012, -8.74556, -8.93365, -9.11323, -9.28313, -9.44225, -9.58959, -9.72438, -9.84612, -9.95463, -10.0501, -10.133, -10.2041, -10.2645},
            { -1.09863, -1.32888, -1.55911, -1.78935, -2.01957, -2.24979, -2.48, -2.71019, -2.94037, -3.17052, -3.40065, -3.63075, -3.8608, -4.0908, -4.32073, -4.55058, -4.78032, -5.00993, -5.23938, -5.46861, -5.69758, -5.92621, -6.15443, -6.38213, -6.60917, -6.8354, -7.06061, -7.28454, -7.50689, -7.72729, -7.94526, -8.16027, -8.37165, -8.57864, -8.78036, -8.9758, -9.16389, -9.34347, -9.51338, -9.67249, -9.81984, -9.95463, -10.0764, -10.1849, -10.2803, -10.3632, -10.4344},
            { -1.09863, -1.32887, -1.55912, -1.78935, -2.01959, -2.24981, -2.48003, -2.71024, -2.94043, -3.17061, -3.40076, -3.63089, -3.86099, -4.09104, -4.32104, -4.55097, -4.78082, -5.01056, -5.24017, -5.46962, -5.69885, -5.92782, -6.15645, -6.38467, -6.61237, -6.83942, -7.06564, -7.29085, -7.51478, -7.73713, -7.95753, -8.1755, -8.39051, -8.60189, -8.80888, -9.0106, -9.20605, -9.39414, -9.57372, -9.74362, -9.90274, -10.0501, -10.1849, -10.3066, -10.4151, -10.5106, -10.5935},
            { -1.09863, -1.32887, -1.55912, -1.78936, -2.0196, -2.24983, -2.48006, -2.71028, -2.94048, -3.17068, -3.40085, -3.63101, -3.86114, -4.09123, -4.32128, -4.55128, -4.78122, -5.01107, -5.24081, -5.47042, -5.69986, -5.92909, -6.15806, -6.3867, -6.61492, -6.84262, -7.06966, -7.29589, -7.52109, -7.74503, -7.96738, -8.18777, -8.40575, -8.62076, -8.83214, -9.03913, -9.24085, -9.43629, -9.62438, -9.80396, -9.97387, -10.133, -10.2803, -10.4151, -10.5369, -10.6454, -10.7408},
            { -1.09862, -1.32887, -1.55912, -1.78937, -2.01961, -2.24985, -2.48008, -2.71031, -2.94052, -3.17073, -3.40092, -3.6311, -3.86126, -4.09138, -4.32148, -4.55153, -4.78153, -5.01147, -5.24131, -5.47106, -5.70067, -5.93011, -6.15934, -6.38831, -6.61695, -6.84517, -7.07286, -7.29991, -7.52614, -7.75134, -7.97528, -8.19763, -8.41802, -8.636, -8.851, -9.06239, -9.26938, -9.4711, -9.66654, -9.85463, -10.0342, -10.2041, -10.3632, -10.5106, -10.6454, -10.7671, -10.8756}, 
            { -1.09862, -1.32887, -1.55912, -1.78937, -2.01962, -2.24986, -2.4801, -2.71033, -2.94056, -3.17077, -3.40098, -3.63117, -3.86135, -4.09151, -4.32163, -4.55173, -4.78178, -5.01178, -5.24172, -5.47156, -5.70131, -5.93092, -6.16036, -6.38959, -6.61856, -6.8472, -7.07542, -7.30311, -7.53016, -7.75639, -7.98159, -8.20553, -8.42788, -8.64827, -8.86625, -9.08126, -9.29264, -9.49963, -9.70135, -9.8968, -10.0849, -10.2645, -10.4344, -10.5935, -10.7408, -10.8756, -10.9974}};
    
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
        
        Alignment* alignment;
        if(align == "gotoh")			{	alignment = new GotohOverlap(gapOpen, gapExtend, match, misMatch, longestBase);			}
        else if(align == "needleman")	{	alignment = new NeedlemanOverlap(gapOpen, match, misMatch, longestBase);				}
        else if(align == "kmer")        {
            if (hasQuality) { alignment = new KmerAlign(kmerSize, numKmers);                                            }
            else { m->mothurOut("[ERROR]: kmerAlign requires quality data, aborting.\n"); m->control_pressed = true;    }
        }

        
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
                savedFQual = new QualityScores(fQual->getName(), fQual->getQualityScores());
                savedRQual = new QualityScores(rQual->getName(), rQual->getQualityScores());
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
                cout << fSeq.getUnaligned() << endl << rSeq.getUnaligned() << endl;
                fQual->printQScores(cout); rQual->printQScores(cout);
                
                rSeq.reverseComplement();
                if (hasQuality) { rQual->flipQScores(); }
                cout << fSeq.getUnaligned() << endl << rSeq.getUnaligned() << endl;
                fQual->printQScores(cout); rQual->printQScores(cout);
                
                //pairwise align
                if(align != "kmer")        {  alignment->align(fSeq.getUnaligned(), rSeq.getUnaligned());  }
                else {  alignment->align(fSeq.getUnaligned(), rSeq.getUnaligned(), fQual->getQualityScores(), rQual->getQualityScores());  }
                    
                map<int, int> ABaseMap = alignment->getSeqAAlnBaseMap();
                map<int, int> BBaseMap = alignment->getSeqBAlnBaseMap();
                fSeq.setAligned(alignment->getSeqAAln());
                rSeq.setAligned(alignment->getSeqBAln());
                int length = fSeq.getAligned().length();
                cout << fSeq.getAligned() << endl << rSeq.getAligned() << endl;
                
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
                
                //if (num == 5) {  cout << fSeq.getStartPos() << '\t' << fSeq.getEndPos() << '\t' << rSeq.getStartPos() << '\t' << rSeq.getEndPos() << endl; exit(1); }
                int overlapStart = fSeq.getStartPos()-1;
                int seq2Start = rSeq.getStartPos()-1;
                
                //bigger of the 2 starting positions is the location of the overlapping start
                if (overlapStart < seq2Start) { //seq2 starts later so take from 0 to seq2Start from seq1
                    overlapStart = seq2Start;
                    for (int i = 0; i < overlapStart; i++) { contig += seq1[i];  if (((seq1[i] != '-') && (seq1[i] != '.'))) { contigScores.push_back(scores1[ABaseMap[i]]); } }
                }else { //seq1 starts later so take from 0 to overlapStart from seq2
                    for (int i = 0; i < overlapStart; i++) {  contig += seq2[i]; if (((seq2[i] != '-') && (seq2[i] != '.'))) {  contigScores.push_back(scores2[BBaseMap[i]]); }  }
                }
                
                int seq1End = fSeq.getEndPos();
                int seq2End = rSeq.getEndPos();
                int overlapEnd = seq1End;
                if (seq2End < overlapEnd) { overlapEnd = seq2End; }  //smallest end position is where overlapping ends
                
                int oStart = contig.length();
                //cout << fSeq.getAligned()  << endl; cout << rSeq.getAligned() << endl;
                
                int firstForward = 0; int seq2FirstForward = 0; int lastReverse = seq1.length(); int seq2lastReverse = seq2.length(); bool firstChooseSeq1 = false; bool lastChooseSeq1 = false;
                for (int i = 0; i < seq1.length(); i++) { if ((seq1[i] != '.') && (seq1[i] != '-')) { if (scores1[ABaseMap[i]] == 2) { firstForward++; }else { break; } } }
                for (int i = 0; i < seq2.length(); i++) { if ((seq2[i] != '.') && (seq2[i] != '-')) { if (scores2[BBaseMap[i]] == 2) { seq2FirstForward++; }else { break; } } }
                if (seq2FirstForward > firstForward) { firstForward = seq2FirstForward; firstChooseSeq1 = true; }
                for (int i = seq1.length()-1; i >= 0; i--) { if ((seq1[i] != '.') && (seq1[i] != '-')) { if (scores1[ABaseMap[i]] == 2) { lastReverse--; }else { break; } } }
                for (int i = seq2.length()-1; i >= 0; i--) { if ((seq2[i] != '.') && (seq2[i] != '-')) { if (scores2[BBaseMap[i]] == 2) { seq2lastReverse--; }else { break; } } }
                if (lastReverse > seq2lastReverse) { lastReverse = seq2lastReverse; lastChooseSeq1 = true; }

                
                cout << firstForward << '\t' << lastReverse << endl;
                
                for (int i = overlapStart; i < overlapEnd; i++) {
                    //cout << seq1[i] << ' ' << seq2[i] << ' ' << scores1[ABaseMap[i]] << ' ' << scores2[BBaseMap[i]] << endl;
                    if (seq1[i] == seq2[i]) {
                        contig += seq1[i];
                        if (hasQuality) {
                            contigScores.push_back(convertProb(qual_match_simple_bayesian[PHREDCLAMP(scores1[ABaseMap[i]])][PHREDCLAMP(scores2[BBaseMap[i]])]));
                        }
                    }else if (((seq1[i] == '.') || (seq1[i] == '-')) && ((seq2[i] != '-') && (seq2[i] != '.'))) { //seq1 is a gap and seq2 is a base, choose seq2, unless quality score for base is below insert. In that case eliminate base
                        if (hasQuality) {
                            if (scores2[BBaseMap[i]] <= insert) { } //
                            else {
                                contig += seq2[i];
                                contigScores.push_back(scores2[BBaseMap[i]]);
                            }
                        }else {  contig += seq2[i]; } //with no quality info, then we keep it?
                    }else if (((seq2[i] == '.') || (seq2[i] == '-')) && ((seq1[i] != '-') && (seq1[i] != '.'))) { //seq2 is a gap and seq1 is a base, choose seq1, unless quality score for base is below insert. In that case eliminate base
                        if (hasQuality) {
                            if (scores1[ABaseMap[i]] <= insert) { } //eliminate base
                            else {
                                contig += seq1[i];
                                contigScores.push_back(scores1[ABaseMap[i]]);
                            }
                        }else { contig += seq1[i]; } //with no quality info, then we keep it?
                    }else if (((seq1[i] != '-') && (seq1[i] != '.')) && ((seq2[i] != '-') && (seq2[i] != '.'))) { //both bases choose one with better quality
                        if (hasQuality) {
                            if (abs(scores1[ABaseMap[i]] - scores2[BBaseMap[i]]) >= deltaq) { //is the difference in qual scores >= deltaq, if yes choose base with higher score
                                char c = seq1[i];
                                if (scores1[ABaseMap[i]] < scores2[BBaseMap[i]]) { c = seq2[i]; }
                                contig += c;
                                if ((i >= firstForward) && (i <= lastReverse)) { //in unmasked section
                                    contigScores.push_back(convertProb(qual_mismatch_simple_bayesian[PHREDCLAMP(scores1[ABaseMap[i]])][PHREDCLAMP(scores2[BBaseMap[i]])]));
                                }else if (i < firstForward) {
                                    if (firstChooseSeq1) { contigScores.push_back(scores1[ABaseMap[i]]); }
                                    else { contigScores.push_back(scores2[BBaseMap[i]]); }
                                }else if ((i > lastReverse)) {
                                    if (lastChooseSeq1) { contigScores.push_back(scores1[ABaseMap[i]]);   }
                                    else { contigScores.push_back(scores2[BBaseMap[i]]); }
                                }else { contigScores.push_back(2); } //N
                            }else { //if no, base becomes n
                                contig += 'N'; contigScores.push_back(2);
                            }
                            numMismatches++;
                        }else { numMismatches++; } //cant decide, so eliminate and mark as mismatch
                    }else { //should never get here
                        m->mothurOut("[ERROR]: case I didn't think of seq1 = " + toString(seq1[i]) + " and seq2 = " + toString(seq2[i]) + "\n");
                    }
                    printf("Overlap seq: %i, %i, %i, %c, %i\n", i, scores1[ABaseMap[i]], scores2[BBaseMap[i]], contig[contig.length()-1], contigScores[contigScores.size()-1]);
                }
                int oend = contig.length();
                if (seq1End < seq2End) { //seq1 ends before seq2 so take from overlap to length from seq2
                    for (int i = overlapEnd; i < length; i++) { contig += seq2[i];  if (((seq2[i] != '-') && (seq2[i] != '.'))) {  contigScores.push_back(scores2[BBaseMap[i]]); }
                        printf("Reverse seq: %i, %c, %i\n", i, contig[contig.length()-1], contigScores[contigScores.size()-1]);
                    }
                }else { //seq2 ends before seq1 so take from overlap to length from seq1
                    for (int i = overlapEnd; i < length; i++) {  contig += seq1[i]; if (((seq1[i] != '-') && (seq1[i] != '.'))) { contigScores.push_back(scores1[ABaseMap[i]]); }
                    printf("Reverse seq: %i, %c, %i\n", i, contig[contig.length()-1], contigScores[contigScores.size()-1]);
                    }
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
                    outScrapQual << ">" << fSeq.getName() << " | " << trashCode << '\t' << commentString << endl;
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
            long long numFastaSeqs = 0;
            fastaFilePos = m->setFilePosFasta(fasta[0], numFastaSeqs, delim); //forward
            if (fastaFilePos.size() < processors) { processors = fastaFilePos.size(); }
            
            long long numRFastaSeqs = 0;
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
                
                if (qual[0] != "NONE") {  fastaFilePos = m->setFilePosFasta(qual[0], numFQualSeqs, delim);  } //forward index or qual file
                if (qual[1] != "NONE") {  qfileFilePos = m->setFilePosFasta(qual[1], numRQualSeqs, delim);  }//reverse index or qual file
                
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
        if(qual.size() != 0)	{	qLines = lines;	} //files with duds
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
//***************************************************************************************************************
/**
 * Convert the probability to a quality score.
 */
int MakeContigsCommand::convertProb(double qProb){
    try {
        
        int lower = 0;
        int upper = 46;
        
        if (qProb < qual_score[0])  { return 1; }
        
        while (lower < upper) {
            int mid = lower + (upper - lower) / 2;
            if (qual_score[mid] == qProb) {
                return mid;
            }
            if (mid == lower) {
                return lower;
            } else if (qual_score[mid] > qProb) {
                upper = mid;
            } else if (qual_score[mid] < qProb) {
                lower = mid + 1;
            }
        }
        
        return lower;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "convertProb");
        exit(1);
    }
}
//**********************************************************************************************************************

