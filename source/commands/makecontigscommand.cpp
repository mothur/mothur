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
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
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
            
            temp = validParameter.validFile(parameters, "ksize", false);	if (temp == "not found"){	temp = "8";			}
            m->mothurConvert(temp, kmerSize);
            
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
            
            qual_score.resize(47);
            qual_score[0] = -2; qual_score[1] = -1.58147; qual_score[2] = -0.996843; qual_score[3] = -0.695524; qual_score[4] = -0.507676; qual_score[5] = -0.38013; qual_score[6] = -0.289268; qual_score[7] = -0.222552; qual_score[8] = -0.172557; qual_score[9] = -0.134552; qual_score[10] = -0.105361; qual_score[11] = -0.0827653; qual_score[12] = -0.0651742; qual_score[13] = -0.0514183; qual_score[14] = -0.0406248; qual_score[15] = -0.0321336; qual_score[16] = -0.0254397; qual_score[17] = -0.0201544; qual_score[18] = -0.0159759; qual_score[19] = -0.0126692; qual_score[20] = -0.0100503; qual_score[21] = -0.007975; qual_score[22] = -0.00632956; qual_score[23] = -0.00502447; qual_score[24] = -0.00398902; qual_score[25] = -0.00316729; qual_score[26] = -0.00251505; qual_score[27] = -0.00199726; qual_score[28] = -0.00158615; qual_score[29] = -0.00125972; qual_score[30] = -0.0010005; qual_score[31] = -0.000794644; qual_score[32] = -0.000631156; qual_score[33] = -0.000501313; qual_score[34] = -0.000398186; qual_score[35] = -0.000316278; qual_score[36] = -0.00025122; qual_score[37] = -0.000199546; qual_score[38] = -0.000158502; qual_score[39] = -0.0001259; qual_score[40] = -0.000100005; qual_score[41] = -7.9436e-05; qual_score[42] = -6.30977e-05; qual_score[43] = -5.012e-05; qual_score[44] = -3.98115e-05; qual_score[45] = -3.16233e-05; qual_score[46] = -2.51192e-05;
            
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
        longestBase = 1000; group = "";
   
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
		if (total != 0) { m->mothurOut("\nTotal of all groups is " + toString(total)); m->mothurOutEndLine(); }
		
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
        numReads = createProcesses(fileInputs, qualOrIndexInputs, outFastaFile, outScrapFastaFile, outQualFile, outScrapQualFile, outMisMatchFile, fastaFileNames, qualFileNames, lines, qLines, group);
        
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
        
        if (gz) {
            numReads = createProcessesGroups(fileInputs, compositeGroupFile, compositeFastaFile, compositeScrapFastaFile, compositeQualFile, compositeScrapQualFile, compositeMisMatchFile, totalGroupCounts);
        }else {
            for (int l = 0; l < fileInputs.size(); l++) {
                int startTime = time(NULL);
                
                if (m->control_pressed) { break; }
                
                m->mothurOut("\n>>>>>\tProcessing file pair " + fileInputs[l][0] + " - " + fileInputs[l][1] + " (files " + toString(l+1) + " of " + toString(fileInputs.size()) + ")\t<<<<<\n");
                
                ffastqfile = fileInputs[l][0];
                rfastqfile = fileInputs[l][1];
                findexfile = fileInputs[l][2];
                rindexfile = fileInputs[l][3];
                group = file2Group[l];
                
                groupCounts.clear();
                groupMap.clear();
                
                //run file as if it was a single
                int thisNumReads = processSingleFileOption(groupCounts);
                numReads += thisNumReads;
                
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
                
                m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - startTime) + " secs to assemble " + toString(thisNumReads) + " reads.\n");	m->mothurOutEndLine();
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
unsigned long long MakeContigsCommand::createProcessesGroups(vector< vector<string> > fileInputs, string compositeGroupFile, string compositeFastaFile, string compositeScrapFastaFile, string compositeQualFile, string compositeScrapQualFile, string compositeMisMatchFile, map<string, int>& totalGroupCounts) {
    try {
        unsigned long long num = 0;
        vector<int> processIDS;
        bool recalc = false;
        
        vector<linePair> startEndIndexes;
        
        //divide files between processors
        int remainingPairs = fileInputs.size();
        int startIndex = 0;
        for (int remainingProcessors = processors; remainingProcessors > 0; remainingProcessors--) {
            int numPairs = remainingPairs; //case for last processor
            if (remainingProcessors != 1) { numPairs = ceil(remainingPairs / remainingProcessors); }
            startEndIndexes.push_back(linePair(startIndex, (startIndex+numPairs))); //startIndex, endIndex
            startIndex = startIndex + numPairs;
            remainingPairs = remainingPairs - numPairs;
        }
        
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        int process = 1;
        
        //loop through and create all the processes you want
        while (process != processors) {
            pid_t pid = fork();
            
            if (pid > 0) {
                processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                process++;
            }else if (pid == 0){
                
                num = driverGroups(fileInputs, startEndIndexes[process].start, startEndIndexes[process].end,
                                   compositeGroupFile + m->mothurGetpid(process) + ".temp",
                                   compositeFastaFile + m->mothurGetpid(process) + ".temp",
                                   compositeScrapFastaFile + m->mothurGetpid(process) + ".temp",
                                   compositeQualFile + m->mothurGetpid(process) + ".temp",
                                   compositeScrapQualFile + m->mothurGetpid(process) + ".temp",
                                   compositeMisMatchFile + m->mothurGetpid(process) + ".temp", totalGroupCounts);
                
                //pass groupCounts to parent
                ofstream out;
                string tempFile = m->mothurGetpid(process) + ".num.temp";
                m->openOutputFile(tempFile, out);
                out << num << endl;
                if (createFileGroup || createOligosGroup) {
                    out << totalGroupCounts.size() << endl;
                    
                    for (map<string, int>::iterator it = totalGroupCounts.begin(); it != totalGroupCounts.end(); it++) {
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
                m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(process) + "\n"); processors = process;
                for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                //wait to die
                for (int i=0;i<processIDS.size();i++) {
                    int temp = processIDS[i];
                    wait(&temp);
                }
                m->control_pressed = false;
                for (int i=0;i<processIDS.size();i++) {
                    m->mothurRemove((toString(processIDS[i]) + ".num.temp"));
                }
                recalc = true;
                break;
            }
        }
        
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->control_pressed = false;  for (int i=0;i<processIDS.size();i++) {m->mothurRemove((toString(processIDS[i]) + ".num.temp"));}processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            
            vector<linePair> startEndIndexes;
            
            //divide files between processors
            int remainingPairs = fileInputs.size();
            int startIndex = 0;
            for (int remainingProcessors = processors; remainingProcessors > 0; remainingProcessors--) {
                int numPairs = remainingPairs; //case for last processor
                if (remainingProcessors != 1) { numPairs = ceil(remainingPairs / remainingProcessors); }
                startEndIndexes.push_back(linePair(startIndex, (startIndex+numPairs))); //startIndex, endIndex
                startIndex = startIndex + numPairs;
                remainingPairs = remainingPairs - numPairs;
            }
            
            num = 0;
            processIDS.resize(0);
            process = 1;
            
            //loop through and create all the processes you want
            while (process != processors) {
                pid_t pid = fork();
                
                if (pid > 0) {
                    processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                    process++;
                }else if (pid == 0){
                    
                    num = driverGroups(fileInputs, startEndIndexes[process].start, startEndIndexes[process].end,
                                       compositeGroupFile + m->mothurGetpid(process) + ".temp",
                                       compositeFastaFile + m->mothurGetpid(process) + ".temp",
                                       compositeScrapFastaFile + m->mothurGetpid(process) + ".temp",
                                       compositeQualFile + m->mothurGetpid(process) + ".temp",
                                       compositeScrapQualFile + m->mothurGetpid(process) + ".temp",
                                       compositeMisMatchFile + m->mothurGetpid(process) + ".temp", totalGroupCounts);
                    
                    //pass groupCounts to parent
                    ofstream out;
                    string tempFile = m->mothurGetpid(process) + ".num.temp";
                    m->openOutputFile(tempFile, out);
                    out << num << endl;
                    if (createFileGroup || createOligosGroup) {
                        out << totalGroupCounts.size() << endl;
                        
                        for (map<string, int>::iterator it = totalGroupCounts.begin(); it != totalGroupCounts.end(); it++) {
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
        }
        
        num = driverGroups(fileInputs, startEndIndexes[0].start, startEndIndexes[0].end, compositeGroupFile, compositeFastaFile, compositeScrapFastaFile, compositeQualFile, compositeScrapQualFile, compositeMisMatchFile, totalGroupCounts);
        
        
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
                        
                        map<string, int>::iterator it = totalGroupCounts.find(group);
                        if (it == totalGroupCounts.end()) {	totalGroupCounts[group] = groupNum; }
                        else { totalGroupCounts[it->first] += groupNum; }
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
            
            m->appendFiles((compositeGroupFile + toString(processIDS[i]) + ".temp"), compositeGroupFile);
            m->mothurRemove((compositeGroupFile + toString(processIDS[i]) + ".temp"));
            
            m->appendFiles((compositeFastaFile + toString(processIDS[i]) + ".temp"), compositeFastaFile);
            m->mothurRemove((compositeFastaFile + toString(processIDS[i]) + ".temp"));
            
            m->appendFiles((compositeScrapFastaFile + toString(processIDS[i]) + ".temp"), compositeScrapFastaFile);
            m->mothurRemove((compositeScrapFastaFile + toString(processIDS[i]) + ".temp"));
            
            m->appendFiles((compositeQualFile + toString(processIDS[i]) + ".temp"), compositeQualFile);
            m->mothurRemove((compositeQualFile + toString(processIDS[i]) + ".temp"));
            
            m->appendFiles((compositeScrapQualFile + toString(processIDS[i]) + ".temp"), compositeScrapQualFile);
            m->mothurRemove((compositeScrapQualFile + toString(processIDS[i]) + ".temp"));
            
            m->appendFiles((compositeMisMatchFile + toString(processIDS[i]) + ".temp"), compositeMisMatchFile);
            m->mothurRemove((compositeMisMatchFile + toString(processIDS[i]) + ".temp"));
        }

        
#else
      processors = 1;
      num = driverGroups(fileInputs, 0, fileInputs.size(), compositeGroupFile, compositeFastaFile, compositeScrapFastaFile, compositeQualFile, compositeScrapQualFile, compositeMisMatchFile, totalGroupCounts);
#endif
        return num;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "createProcessesGroups");
        exit(1);
    }
}
//**********************************************************************************************************************
unsigned long long MakeContigsCommand::driverGroups(vector< vector<string> > fileInputs, int start, int end, string compositeGroupFile, string compositeFastaFile, string compositeScrapFastaFile, string compositeQualFile, string compositeScrapQualFile, string compositeMisMatchFile, map<string, int>& totalGroupCounts) {
    try {
        unsigned long long numReads = 0;
        delim = '@';
        
        for (int l = start; l < end; l++) {
            int startTime = time(NULL);
            
            if (m->control_pressed) { break; }
            
            m->mothurOut("\n>>>>>\tProcessing file pair " + fileInputs[l][0] + " - " + fileInputs[l][1] + " (files " + toString(l+1) + " of " + toString(fileInputs.size()) + ")\t<<<<<\n");
            
            ffastqfile = fileInputs[l][0];
            rfastqfile = fileInputs[l][1];
            findexfile = fileInputs[l][2];
            rindexfile = fileInputs[l][3];
            group = file2Group[l];
            
            groupCounts.clear();
            groupMap.clear();
            
            vector<string> thisFileInputs;
            vector<string> thisQualOrIndexInputs;
            vector<linePair> thisLines;
            vector<linePair> thisQLines;
            map<string, string> variables;
            string thisOutputDir = outputDir;
            
            string inputFile = ffastqfile;
            if (outputDir == "") {  thisOutputDir = m->hasPath(inputFile); }
            variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(inputFile));
            variables["[tag]"] = "trim";
            outQualFile = getOutputFileName("qfile",variables);
            variables["[tag]"] = "scrap";
            outScrapQualFile = getOutputFileName("qfile",variables);
            
            thisFileInputs.push_back(ffastqfile); thisFileInputs.push_back(rfastqfile);
            if ((findexfile != "") || (rindexfile != "")){
                thisQualOrIndexInputs.push_back("NONE"); thisQualOrIndexInputs.push_back("NONE");
                if (findexfile != "") { thisQualOrIndexInputs[0] = findexfile; }
                if (rindexfile != "") { thisQualOrIndexInputs[1] = rindexfile; }
            }
            
            variables["[tag]"] = "trim";
            outFastaFile = getOutputFileName("fasta",variables);
            variables["[tag]"] = "scrap";
            outScrapFastaFile = getOutputFileName("fasta",variables);
            variables["[tag]"] = "";
            outMisMatchFile = getOutputFileName("report",variables);
            
            //fake out lines - we are just going to check for end of file. Work is divided by number of files per processor.
            thisLines.push_back(linePair(0, 1000)); thisLines.push_back(linePair(0, 1000)); //fasta[0], fasta[1] - forward and reverse
            thisQLines.push_back(linePair(0, 1000)); thisQLines.push_back(linePair(0, 1000));  //qual[0], qual[1] - forward and reverse
            
            vector<vector<string> > fastaFileNames, qualFileNames;
            map<string, string> uniqueFastaNames;// so we don't add the same groupfile multiple times
            createOligosGroup = false;
            oligos = new Oligos();
            numBarcodes = 0; numFPrimers= 0; numLinkers= 0; numSpacers = 0; numRPrimers = 0;
            
            if(oligosfile != "")                        {       createOligosGroup = getOligos(fastaFileNames, qualFileNames, variables["[filename]"], uniqueFastaNames);    }
            if (createOligosGroup || createFileGroup)   {       outputGroupFileName = getOutputFileName("group",variables);                                                 }
            
            //give group in file file precedence
            if (createFileGroup) {  createOligosGroup = false; }
         
            ofstream temp;
            m->openOutputFile(outFastaFile, temp);		temp.close();
            m->openOutputFile(outScrapFastaFile, temp);		temp.close();
            m->openOutputFile(outQualFile, temp);		temp.close();
            m->openOutputFile(outScrapQualFile, temp);		temp.close();
            
            m->mothurOut("Making contigs...\n");
            
            unsigned long long thisNumReads = driver(thisFileInputs, thisQualOrIndexInputs, outFastaFile, outScrapFastaFile,  outQualFile, outScrapQualFile, outMisMatchFile, fastaFileNames, qualFileNames, thisLines[0], thisLines[1], thisQLines[0], thisQLines[1], group);
            
            numReads += thisNumReads;
            
            m->mothurOut("Done.\n");
            
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
            
            //append to combo files
            if (createFileGroup || createOligosGroup) {
                ofstream outCGroup;
                if (l == 0) {   m->openOutputFile(compositeGroupFile, outCGroup);  outputNames.push_back(compositeGroupFile); outputTypes["group"].push_back(compositeGroupFile);         }
                else {          m->openOutputFileAppend(compositeGroupFile, outCGroup);     }
                
                if (!allFiles) {
                    m->mothurRemove(outputGroupFileName);
                }else {
                    ofstream outGroup;
                    m->openOutputFile(outputGroupFileName, outGroup);
                    
                    for (map<string, string>::iterator itGroup = groupMap.begin(); itGroup != groupMap.end(); itGroup++) {
                        outCGroup << itGroup->first << '\t' << itGroup->second << endl;
                        outGroup << itGroup->first << '\t' << itGroup->second << endl;
                    }
                    outGroup.close();
                }
                outCGroup.close();
                
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

            m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - startTime) + " secs to assemble " + toString(thisNumReads) + " reads.\n");	m->mothurOutEndLine();
        }
        
        return numReads;

    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "driverGroups");
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
        bool recalc = false;
       
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
                m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(process) + "\n"); processors = process;
                for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                //wait to die
                for (int i=0;i<processIDS.size();i++) {
                    int temp = processIDS[i];
                    wait(&temp);
                }
                m->control_pressed = false;
                for (int i=0;i<processIDS.size();i++) {
                    m->mothurRemove((toString(processIDS[i]) + ".num.temp"));
                }
                recalc = true;
                break;
			}
		}
		
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->control_pressed = false;  for (int i=0;i<processIDS.size();i++) {m->mothurRemove((toString(processIDS[i]) + ".num.temp"));}processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            
            //redo file divide
            lines.clear();
            setLines(fileInputs, qualOrIndexFiles, lines, qLines, delim);
            
            num = 0;
            processIDS.resize(0);
            process = 1;
            
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
        }

        
        ofstream temp;
		m->openOutputFile(outputFasta, temp);		temp.close();
        m->openOutputFile(outputScrapFasta, temp);		temp.close();
        if (outputQual != "") {
            m->openOutputFile(outputQual, temp);		temp.close();
            m->openOutputFile(outputScrapQual, temp);		temp.close();
        }
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
			contigsData* tempcontig = new contigsData(format, delim, group, fileInputs, qualOrIndexFiles, (outputFasta + extension), (outputScrapFasta + extension), (outputQual + extension), (outputScrapQual + extension), (outputMisMatches + extension), align, m, match, misMatch, gapOpen, gapExtend, insert, deltaq, tempFASTAFileNames, tempQUALFileNames, oligosfile, reorient, pdiffs, bdiffs, tdiffs, kmerSize, createOligosGroup, createFileGroup, allFiles, trimOverlap, gz, lines[spot].start, lines[spot].end, lines[spot+1].start, lines[spot+1].end, qLines[spot].start, qLines[spot].end, qLines[spot+1].start, qLines[spot+1].end, h);
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
        if (outputQual != "") {
            m->openOutputFile(outputQual, temp3);		temp3.close();
            m->openOutputFile(outputScrapQual, temp4);		temp4.close();
        }
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
            
            if (outputQual != "") {
                m->appendFiles((outputQual + toString(processIDS[i]) + ".temp"), outputQual);
                m->mothurRemove((outputQual + toString(processIDS[i]) + ".temp"));
            
                m->appendFiles((outputScrapQual + toString(processIDS[i]) + ".temp"), outputScrapQual);
                m->mothurRemove((outputScrapQual + toString(processIDS[i]) + ".temp"));
            }
            
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
unsigned long long MakeContigsCommand::driver(vector<string> inputFiles, vector<string> qualOrIndexFiles, string outputFasta, string outputScrapFasta, string outputQual, string outputScrapQual,  string outputMisMatches, vector<vector<string> > fastaFileNames, vector<vector<string> > qualFileNames, linePair linesInput, linePair linesInputReverse, linePair qlinesInput, linePair qlinesInputReverse, string group){
    try {
        
        vector< vector<double> > qual_match_simple_bayesian;
        qual_match_simple_bayesian.resize(47);
        for (int i = 0; i < qual_match_simple_bayesian.size(); i++) { qual_match_simple_bayesian[i].resize(47);  }
        
        vector< vector<double> > qual_mismatch_simple_bayesian;
        qual_mismatch_simple_bayesian.resize(47);
        for (int i = 0; i < qual_mismatch_simple_bayesian.size(); i++) { qual_mismatch_simple_bayesian[i].resize(47);  }
        
        loadQmatchValues(qual_match_simple_bayesian, qual_mismatch_simple_bayesian);
        
        unsigned long long num = 0;
        string thisfqualindexfile, thisrqualindexfile, thisffastafile, thisrfastafile;
        thisfqualindexfile = ""; thisrqualindexfile = "";
        thisffastafile = inputFiles[0]; thisrfastafile = inputFiles[1];
        if (qualOrIndexFiles.size() != 0) {
            thisfqualindexfile = qualOrIndexFiles[0];
            thisrqualindexfile = qualOrIndexFiles[1];
        }
        
        if (m->debug) {  m->mothurOut("[DEBUG]: ffasta = " + thisffastafile + ".\n[DEBUG]: rfasta = " + thisrfastafile + ".\n[DEBUG]: fqualindex = " + thisfqualindexfile + ".\n[DEBUG]: rqualindex = " + thisfqualindexfile + ".\n"); }
        
        ifstream inFFasta, inRFasta, inFQualIndex, inRQualIndex;
        boost::iostreams::filtering_istream inFF, inRF, inFQ, inRQ;
        if (!gz) { //plain text files
            m->openInputFile(thisffastafile, inFFasta);
            m->openInputFile(thisrfastafile, inRFasta);
            
            inFFasta.seekg(linesInput.start);
            inRFasta.seekg(linesInputReverse.start);
        }else { //compressed files - no need to seekg because compressed files divide workload differently
            m->openInputFileBinary(thisffastafile, inFFasta, inFF);
            m->openInputFileBinary(thisrfastafile, inRFasta, inRF);
        }
        
        ofstream outFasta, outMisMatch, outScrapFasta, outQual, outScrapQual;
        if (thisfqualindexfile != "") {
            if (thisfqualindexfile != "NONE") {
                if (!gz) { //plain text files
                    m->openInputFile(thisfqualindexfile, inFQualIndex);
                    inFQualIndex.seekg(qlinesInput.start);
                }else { m->openInputFileBinary(thisfqualindexfile, inFQualIndex, inFQ); } //compressed files - no need to seekg because compressed files divide workload differently
            }
            else {  thisfqualindexfile = ""; }
            if (thisrqualindexfile != "NONE") {
                if (!gz) { //plain text files
                    m->openInputFile(thisrqualindexfile, inRQualIndex);
                    inRQualIndex.seekg(qlinesInputReverse.start);
                }else { m->openInputFileBinary(thisrqualindexfile, inRQualIndex, inRQ); } //compressed files - no need to seekg because compressed files divide workload differently
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
        
        if (m->debug) { if (hasQuality) { m->mothurOut("[DEBUG]: hasQuality = true\n");  } else { m->mothurOut("[DEBUG]: hasQuality = false\n"); } }
        
        TrimOligos trimOligos(pdiffs, bdiffs, 0, 0, oligos->getPairedPrimers(), oligos->getPairedBarcodes());

        TrimOligos* rtrimOligos = NULL;
        if (reorient) {  rtrimOligos = new TrimOligos(pdiffs, bdiffs, 0, 0, oligos->getReorientedPairedPrimers(), oligos->getReorientedPairedBarcodes()); numBarcodes = oligos->getReorientedPairedBarcodes().size();    }
        
        Alignment* alignment;
        if(align == "gotoh")			{	alignment = new GotohOverlap(gapOpen, gapExtend, match, misMatch, longestBase);			}
        else if(align == "needleman")	{	alignment = new NeedlemanOverlap(gapOpen, match, misMatch, longestBase);				}
        else if(align == "kmer")        {   alignment = new KmerAlign(kmerSize);                                                    }

        bool good = true;
        while (good) {
            
            if (m->control_pressed) { break; }
            
            int success = 1;
            string trashCode = "";
            string commentString = "";
            int currentSeqsDiffs = 0;
            
            bool hasIndex = false;
            bool ignore = false;
            Sequence fSeq, rSeq;
            QualityScores* fQual = NULL; QualityScores* rQual = NULL;
            QualityScores* savedFQual = NULL; QualityScores* savedRQual = NULL;
            Sequence findexBarcode("findex", "NONE");  Sequence rindexBarcode("rindex", "NONE");
            
            //read from input files
            if (gz) {   ignore = read(fSeq, rSeq, fQual, rQual, savedFQual, savedRQual, findexBarcode, rindexBarcode, delim, hasIndex, inFF, inRF, inFQ, inRQ, thisfqualindexfile, thisrqualindexfile); }
            else    {   ignore = read(fSeq, rSeq, fQual, rQual, savedFQual, savedRQual, findexBarcode, rindexBarcode, delim, hasIndex, inFFasta, inRFasta, inFQualIndex, inRQualIndex, thisfqualindexfile, thisrqualindexfile); }
            
            //remove primers and barcodes if neccessary
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
                
                //assemble reads
                string contig = "";
                int oend, oStart;
                int numMismatches = 0;
                vector<int> contigScores = assembleFragments(qual_match_simple_bayesian, qual_mismatch_simple_bayesian, fSeq, rSeq, fQual, rQual, savedFQual, savedRQual, hasQuality, alignment, contig, trashCode, oend, oStart, numMismatches);
                
                //prints results to outputs files
                if(trashCode.length() == 0){
                    bool ignore = false;
                    
                    if (m->debug) { m->mothurOut("[DEBUG]: " + fSeq.getName()); }
                    
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
                        if (hasQuality) {
                            outQual << ">" << fSeq.getName() << '\t' << commentString << endl;
                            for (int i = 0; i < contigScores.size(); i++) { outQual << contigScores[i] << " "; }  outQual << endl;
                        }
                        int numNs = 0;
                        for (int i = 0; i < contig.length(); i++) { if (contig[i] == 'N') { numNs++; }  }
                        outMisMatch << fSeq.getName() << '\t' << contig.length() << '\t' << (oend-oStart) << '\t' << oStart << '\t' << oend << '\t' << numMismatches << '\t' << numNs << endl;
                        
                        if (allFiles) {
                            ofstream output;
                            m->openOutputFileAppend(fastaFileNames[barcodeIndex][primerIndex], output);
                            output << ">" << fSeq.getName() << '\t' << commentString << endl << contig << endl;
                            output.close();
                            
                            if (hasQuality) {
                                ofstream output2;
                                m->openOutputFileAppend(qualFileNames[barcodeIndex][primerIndex], output2);
                                output2 << ">" << fSeq.getName() << '\t' << commentString << endl;
                                for (int i = 0; i < contigScores.size(); i++) { output2 << contigScores[i] << " "; }  output2 << endl;
                                output2.close();
                            }
                        }
                    }
                }else {
                    //output
                    outScrapFasta << ">" << fSeq.getName() << " | " << trashCode << '\t' << commentString << endl << contig << endl;
                    if (hasQuality) {
                        outScrapQual << ">" << fSeq.getName() << " | " << trashCode << '\t' << commentString << endl;
                        for (int i = 0; i < contigScores.size(); i++) { outScrapQual << contigScores[i] << " "; }  outScrapQual << endl;
                    }
                }
            }
            num++;
            
            #if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
                if (!gz) {
                    unsigned long long pos = inFFasta.tellg();
                    if ((pos == -1) || (pos >= linesInput.end)) { good = false; break; }
                }else {
                    if (inFF.eof() || inRF.eof()) { good = false; break; }
                }
            #else
                if (!gz) {
                    if ((inFFasta.eof()) || (inRFasta.eof())) { good = false; break; }
                }else {
                    if (inFF.eof() || inRF.eof()) { good = false; break; }
                }
            #endif
            
			//report progress
            if((num) % 1000 == 0){	m->mothurOutJustToScreen(toString(num)); m->mothurOutEndLine();	 }
		}
        
		//report progress
		if((num) % 1000 != 0){	m->mothurOut(toString(num)); m->mothurOutEndLine();		}
        
        //close files
        inFFasta.close();
        inRFasta.close();
        if (gz) { inFF.pop(); inRF.pop(); }
        outFasta.close();
        outScrapFasta.close();
        outMisMatch.close();
        if (delim == '@') {
            if (thisfqualindexfile != "") { inFQualIndex.close(); if (gz) { inFQ.pop(); } }
            if (thisrqualindexfile != "") { inRQualIndex.close(); if (gz) { inRQ.pop(); } }
            outQual.close();
            outScrapQual.close();
        }else{
            if (hasQuality) {
                inFQualIndex.close();
                inRQualIndex.close();
                if (gz) { inFQ.pop(); inRQ.pop(); }
                outQual.close();
                outScrapQual.close();
            }
        }
        
        //cleanup memory
        delete alignment;
        if (reorient) { delete rtrimOligos; }
        
        if (m->control_pressed) {
            m->mothurRemove(outputFasta); m->mothurRemove(outputScrapFasta); m->mothurRemove(outputMisMatches);
            if (hasQuality) { m->mothurRemove(outputQual); m->mothurRemove(outputScrapQual); }
        }
    
        return num;
    }
	catch(exception& e) {
		m->errorOut(e, "MakeContigsCommand", "driver");
		exit(1);
	}
}
/**************************************************************************************************/
//vector<int> contigScores = assembleFragments(qual_match_simple_bayesian, qual_mismatch_simple_bayesian, fSeq, rSeq, alignment, contig);
vector<int> MakeContigsCommand::assembleFragments(vector< vector<double> >&qual_match_simple_bayesian, vector< vector<double> >& qual_mismatch_simple_bayesian, Sequence& fSeq, Sequence& rSeq, QualityScores*& fQual, QualityScores*& rQual, QualityScores*& savedFQual, QualityScores*& savedRQual, bool hasQuality, Alignment*& alignment, string& contig, string& trashCode, int& oend, int& oStart, int& numMismatches) {
    try {
        vector<int> contigScores;
        
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
        string seq1 = fSeq.getAligned();
        string seq2 = rSeq.getAligned();
        vector<int> scores1, scores2;
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
            for (int i = 0; i < overlapStart; i++) { contig += seq1[i];  if (hasQuality) { if (((seq1[i] != '-') && (seq1[i] != '.'))) { contigScores.push_back(scores1[ABaseMap[i]]); } } }
        }else { //seq1 starts later so take from 0 to overlapStart from seq2
            for (int i = 0; i < overlapStart; i++) {  contig += seq2[i]; if (hasQuality) { if (((seq2[i] != '-') && (seq2[i] != '.'))) {  contigScores.push_back(scores2[BBaseMap[i]]); }  } }
        }
        
        int seq1End = fSeq.getEndPos();
        int seq2End = rSeq.getEndPos();
        int overlapEnd = seq1End;
        if (seq2End < overlapEnd) { overlapEnd = seq2End; }  //smallest end position is where overlapping ends
        
        oStart = contig.length();
        //cout << fSeq.getAligned()  << endl; cout << rSeq.getAligned() << endl;
        
        int firstForward = 0; int seq2FirstForward = 0; int lastReverse = seq1.length(); int seq2lastReverse = seq2.length(); bool firstChooseSeq1 = false; bool lastChooseSeq1 = false;
        if (hasQuality) {
            for (int i = 0; i < seq1.length(); i++) { if ((seq1[i] != '.') && (seq1[i] != '-')) { if (scores1[ABaseMap[i]] == 2) { firstForward++; }else { break; } } }
            for (int i = 0; i < seq2.length(); i++) { if ((seq2[i] != '.') && (seq2[i] != '-')) { if (scores2[BBaseMap[i]] == 2) { seq2FirstForward++; }else { break; } } }
            if (seq2FirstForward > firstForward) { firstForward = seq2FirstForward; firstChooseSeq1 = true; }
            for (int i = seq1.length()-1; i >= 0; i--) { if ((seq1[i] != '.') && (seq1[i] != '-')) { if (scores1[ABaseMap[i]] == 2) { lastReverse--; }else { break; } } }
            for (int i = seq2.length()-1; i >= 0; i--) { if ((seq2[i] != '.') && (seq2[i] != '-')) { if (scores2[BBaseMap[i]] == 2) { seq2lastReverse--; }else { break; } } }
            if (lastReverse > seq2lastReverse) { lastReverse = seq2lastReverse; lastChooseSeq1 = true; }
        }
        
        //cout << firstForward << '\t' << lastReverse << endl;
        
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
            //printf("Overlap seq: %i, %i, %i, %c, %i\n", i, scores1[ABaseMap[i]], scores2[BBaseMap[i]], contig[contig.length()-1], contigScores[contigScores.size()-1]);
        }
        oend = contig.length();
        if (seq1End < seq2End) { //seq1 ends before seq2 so take from overlap to length from seq2
            for (int i = overlapEnd; i < length; i++) { contig += seq2[i];  if (hasQuality) { if (((seq2[i] != '-') && (seq2[i] != '.'))) {  contigScores.push_back(scores2[BBaseMap[i]]); } }
            }
        }else { //seq2 ends before seq1 so take from overlap to length from seq1
            for (int i = overlapEnd; i < length; i++) {  contig += seq1[i];  if (hasQuality) { if (((seq1[i] != '-') && (seq1[i] != '.'))) { contigScores.push_back(scores1[ABaseMap[i]]); } }
            }
        }
        //cout << contig << endl;
        
        if (trimOverlap) {
            contig = contig.substr(overlapStart, oend-oStart);
            if (contig.length() == 0) { trashCode += "l"; }
            if (hasQuality) {
                vector<int> newContigScores;
                for (int i = overlapStart; i < oend; i++)  { newContigScores.push_back(contigScores[i]);  }
                contigScores = newContigScores;
            }
        }

        return contigScores;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "assembleFragments");
        exit(1);
    }
}
/**************************************************************************************************/
//ignore = read(fSeq, rSeq, fQual, rQual, savedFQual, savedRQual, findexBarcode, rindexBarcode, delim, hasIndex, inFF, inRF, inFQ, inRQ);
bool MakeContigsCommand::read(Sequence& fSeq, Sequence& rSeq, QualityScores*& fQual, QualityScores*& rQual, QualityScores*& savedFQual, QualityScores*& savedRQual, Sequence& findexBarcode, Sequence& rindexBarcode, char delim, bool& hasIndex, boost::iostreams::filtering_istream& inFF, boost::iostreams::filtering_istream& inRF, boost::iostreams::filtering_istream& inFQ, boost::iostreams::filtering_istream& inRQ, string thisfqualindexfile, string thisrqualindexfile) {
    try {
        bool ignore = false;
        
        if (delim == '@') { //fastq files
            bool tignore;
            FastqRead fread(inFF, tignore, format);
            FastqRead rread(inRF, ignore, format);
            if (!checkName(fread, rread)) {
                FastqRead f2read(inFF, tignore, format);
                if (!checkName(f2read, rread)) {
                    FastqRead r2read(inRF, ignore, format);
                    if (!checkName(fread, r2read)) {
                        m->mothurOut("[WARNING]: name mismatch in forward and reverse fastq file. Ignoring, " + fread.getName() + ".\n"); ignore = true;
                    }else { rread = r2read; }
                }else { fread = f2read; }
            }
            if (tignore) { ignore=true; }
            fSeq.setName(fread.getName()); fSeq.setAligned(fread.getSeq());
            rSeq.setName(rread.getName()); rSeq.setAligned(rread.getSeq());
            fQual = new QualityScores(fread.getName(), fread.getScores());
            rQual = new QualityScores(rread.getName(), rread.getScores());
            savedFQual = new QualityScores(fQual->getName(), fQual->getQualityScores());
            savedRQual = new QualityScores(rQual->getName(), rQual->getQualityScores());
            if (thisfqualindexfile != "") { //forward index file
                FastqRead firead(inFQ, tignore, format);
                if (tignore) { ignore=true; }
                findexBarcode.setAligned(firead.getSeq());
                if (!checkName(fread, firead)) {
                    FastqRead f2iread(inFQ, tignore, format);
                    if (tignore) { ignore=true; }
                    if (!checkName(fread, f2iread)) {
                        m->mothurOut("[WARNING]: name mismatch in forward index file. Ignoring, " + fread.getName() + ".\n"); ignore = true;
                    }else { firead = f2iread; findexBarcode.setAligned(firead.getSeq()); }
                }
                hasIndex = true;
            }
            if (thisrqualindexfile != "") { //reverse index file
                FastqRead riread(inRQ, tignore, format);
                if (tignore) { ignore=true; }
                rindexBarcode.setAligned(riread.getSeq());
                if (!checkName(fread, riread)) {
                    FastqRead r2iread(inRQ, tignore, format); m->gobble(inRQ);
                    if (tignore) { ignore=true; }
                    if (!checkName(fread, r2iread)) {
                        m->mothurOut("[WARNING]: name mismatch in forward index file. Ignoring, " + fread.getName() + ".\n"); ignore = true;
                    }else { riread = r2iread; rindexBarcode.setAligned(riread.getSeq()); }
                }
                hasIndex = true;
            }
        }else { //reading fasta and maybe qual
            Sequence tfSeq(inFF);
            Sequence trSeq(inRF);
            if (!checkName(tfSeq, trSeq)) {
                Sequence t2fSeq(inFF);
                if (!checkName(t2fSeq, trSeq)) {
                    Sequence t2rSeq(inRF);
                    if (!checkName(tfSeq, t2rSeq)) {
                        m->mothurOut("[WARNING]: name mismatch in forward and reverse fasta file. Ignoring, " + tfSeq.getName() + ".\n"); ignore = true;
                    }else { trSeq = t2fSeq; }
                }else { tfSeq = t2fSeq; }
            }
            fSeq.setName(tfSeq.getName()); fSeq.setAligned(tfSeq.getAligned());
            rSeq.setName(trSeq.getName()); rSeq.setAligned(trSeq.getAligned());
            if (thisfqualindexfile != "") {
                fQual = new QualityScores(inFQ); m->gobble(inFQ);
                rQual = new QualityScores(inRQ); m->gobble(inRQ);
                if (!checkName(*fQual, *rQual)) {
                    m->mothurOut("[WARNING]: name mismatch in forward and reverse qual file. Ignoring, " + fQual->getName() + ".\n"); ignore = true;
                }
                savedFQual = new QualityScores(fQual->getName(), fQual->getQualityScores());
                savedRQual = new QualityScores(rQual->getName(), rQual->getQualityScores());
                if (fQual->getName() != tfSeq.getName()) { m->mothurOut("[WARNING]: name mismatch in forward quality file. Ignoring, " + tfSeq.getName() + ".\n"); ignore = true; }
                if (rQual->getName() != trSeq.getName()) { m->mothurOut("[WARNING]: name mismatch in reverse quality file. Ignoring, " + trSeq.getName() + ".\n"); ignore = true; }
            }
            if (tfSeq.getName() != trSeq.getName()) { m->mothurOut("[WARNING]: name mismatch in forward and reverse fasta file. Ignoring, " + tfSeq.getName() + ".\n"); ignore = true; }
        }

        return ignore;
        
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "read");
        exit(1);
    }
}
/**************************************************************************************************/

bool MakeContigsCommand::read(Sequence& fSeq, Sequence& rSeq, QualityScores*& fQual, QualityScores*& rQual, QualityScores*& savedFQual, QualityScores*& savedRQual, Sequence& findexBarcode, Sequence& rindexBarcode, char delim, bool& hasIndex, ifstream& inFFasta, ifstream& inRFasta, ifstream& inFQualIndex, ifstream& inRQualIndex, string thisfqualindexfile, string thisrqualindexfile) {
    try {
        bool ignore = false;
        
        if (delim == '@') { //fastq files
            bool tignore;
            FastqRead fread(inFFasta, tignore, format); m->gobble(inFFasta);
            FastqRead rread(inRFasta, ignore, format); m->gobble(inRFasta);
            if (!checkName(fread, rread)) {
                FastqRead f2read(inFFasta, tignore, format); m->gobble(inFFasta);
                if (!checkName(f2read, rread)) {
                    FastqRead r2read(inRFasta, ignore, format); m->gobble(inRFasta);
                    if (!checkName(fread, r2read)) {
                        m->mothurOut("[WARNING]: name mismatch in forward and reverse fastq file. Ignoring, " + fread.getName() + ".\n"); ignore = true;
                    }else { rread = r2read; }
                }else { fread = f2read; }
            }
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
                if (!checkName(fread, firead)) {
                    FastqRead f2iread(inFQualIndex, tignore, format); m->gobble(inFQualIndex);
                    if (tignore) { ignore=true; }
                    if (!checkName(fread, f2iread)) {
                        m->mothurOut("[WARNING]: name mismatch in forward index file. Ignoring, " + fread.getName() + ".\n"); ignore = true;
                    }else { firead = f2iread; findexBarcode.setAligned(firead.getSeq()); }
                }
                hasIndex = true;
            }
            if (thisrqualindexfile != "") { //reverse index file
                FastqRead riread(inRQualIndex, tignore, format); m->gobble(inRQualIndex);
                if (tignore) { ignore=true; }
                rindexBarcode.setAligned(riread.getSeq());
                if (!checkName(fread, riread)) {
                    FastqRead r2iread(inRQualIndex, tignore, format); m->gobble(inRQualIndex);
                    if (tignore) { ignore=true; }
                    if (!checkName(fread, r2iread)) {
                        m->mothurOut("[WARNING]: name mismatch in forward index file. Ignoring, " + fread.getName() + ".\n"); ignore = true;
                    }else { riread = r2iread; rindexBarcode.setAligned(riread.getSeq()); }
                }
                hasIndex = true;
            }
        }else { //reading fasta and maybe qual
            Sequence tfSeq(inFFasta); m->gobble(inFFasta);
            Sequence trSeq(inRFasta); m->gobble(inRFasta);
            if (!checkName(tfSeq, trSeq)) {
                Sequence t2fSeq(inFFasta); m->gobble(inFFasta);
                if (!checkName(t2fSeq, trSeq)) {
                    Sequence t2rSeq(inRFasta); m->gobble(inRFasta);
                    if (!checkName(tfSeq, t2rSeq)) {
                        m->mothurOut("[WARNING]: name mismatch in forward and reverse fasta file. Ignoring, " + tfSeq.getName() + ".\n"); ignore = true;
                    }else { trSeq = t2fSeq; }
                }else { tfSeq = t2fSeq; }
            }
            fSeq.setName(tfSeq.getName()); fSeq.setAligned(tfSeq.getAligned());
            rSeq.setName(trSeq.getName()); rSeq.setAligned(trSeq.getAligned());
            if (thisfqualindexfile != "") {
                fQual = new QualityScores(inFQualIndex); m->gobble(inFQualIndex);
                rQual = new QualityScores(inRQualIndex); m->gobble(inRQualIndex);
                if (!checkName(*fQual, *rQual)) {
                    m->mothurOut("[WARNING]: name mismatch in forward and reverse qual file. Ignoring, " + fQual->getName() + ".\n"); ignore = true;
                }
                savedFQual = new QualityScores(fQual->getName(), fQual->getQualityScores());
                savedRQual = new QualityScores(rQual->getName(), rQual->getQualityScores());
                if (fQual->getName() != tfSeq.getName()) { m->mothurOut("[WARNING]: name mismatch in forward quality file. Ignoring, " + tfSeq.getName() + ".\n"); ignore = true; }
                if (rQual->getName() != trSeq.getName()) { m->mothurOut("[WARNING]: name mismatch in reverse quality file. Ignoring, " + trSeq.getName() + ".\n"); ignore = true; }
            }
            if (tfSeq.getName() != trSeq.getName()) { m->mothurOut("[WARNING]: name mismatch in forward and reverse fasta file. Ignoring, " + tfSeq.getName() + ".\n"); ignore = true; }
        }
        
        return ignore;
        
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "read");
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
        map<string, int> firstSeqNames; map<string, int> trimmedNames;
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
            trimmedNames[name.substr(0, name.length()-1)];
            in.close();
        }
        
        map<string, int> copy;  map<string, int> tcopy;
        if (qual.size() != 0) { copy = firstSeqNames; tcopy = trimmedNames; }
        
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
                    map<string, int>::iterator itTrimmed = trimmedNames.find(name.substr(0, name.length()-1));
                    
                    if (it != firstSeqNames.end())  { //this is the start of a new chunk
                        unsigned long long pos = in2.tellg();
                        qfileFilePos.push_back(pos - input.length() - 1);
                        firstSeqNames.erase(it);
                    }else if (itTrimmed != trimmedNames.end()) {
                        unsigned long long pos = in2.tellg();
                        qfileFilePos.push_back(pos - input.length() - 1);
                        trimmedNames.erase(itTrimmed);
                    }
                }
            }
            
            if ((firstSeqNames.size() == 0) || (trimmedNames.size() == 0)) { break; }
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

        
        if ((firstSeqNames.size() != 0) && (trimmedNames.size() != 0)){
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
            trimmedNames = tcopy;
            
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
                            map<string, int>::iterator itTrimmed = trimmedNames.find(name.substr(0, name.length()-1));
                            
                            if(it != firstSeqNames.end()) { //this is the start of a new chunk
                                unsigned long long pos = inQual.tellg();
                                qfileFilePos.push_back(pos - input.length() - 1);
                                firstSeqNames.erase(it);
                            }else if (itTrimmed != trimmedNames.end()) {
                                unsigned long long pos = inQual.tellg();
                                qfileFilePos.push_back(pos - input.length() - 1);
                                trimmedNames.erase(itTrimmed);
                            }
                        }
                    }
                    
                    if ((firstSeqNames.size() == 0) || (trimmedNames.size() == 0)) { break; }
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
                
                
                if ((firstSeqNames.size() != 0) && (trimmedNames.size() != 0)){
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
            trimmedNames = tcopy;
            
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
                            map<string, int>::iterator itTrimmed = trimmedNames.find(name.substr(0, name.length()-1));
                            
                            if(it != firstSeqNames.end()) { //this is the start of a new chunk
                                unsigned long long pos = inQual2.tellg();
                                temp.push_back(pos - input.length() - 1);
                                firstSeqNames.erase(it);
                            }else if (itTrimmed != trimmedNames.end()) {
                                unsigned long long pos = inQual2.tellg();
                                qfileFilePos.push_back(pos - input.length() - 1);
                                trimmedNames.erase(itTrimmed);
                            }
                        }
                    }
                    
                    if ((firstSeqNames.size() == 0) || (trimmedNames.size() == 0)) { break; }
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
                
                
                if ((firstSeqNames.size() != 0) && (trimmedNames.size() != 0)){
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
            qLines.push_back(linePair(0, 1000)); qLines.push_back(linePair(0, 1000));  //qual[0], qual[1] - forward and reverse
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
        bool allGZ = true;
        bool allPlainTxt = true;
        
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
                if (m->isGZ(forward)[1]) { allPlainTxt = false;  }
                else {   allGZ = false;  }
                if (m->isGZ(reverse)[1]) { allPlainTxt = false;  }
                else {   allGZ = false;  }
                if ((findex != "") && (findex != "NONE")) {
                    if (m->isGZ(findex)[1]) { allPlainTxt = false;  }
                    else {   allGZ = false;  }
                }
                if ((rindex != "") && (rindex != "NONE")) {
                    if (m->isGZ(rindex)[1]) { allPlainTxt = false;  }
                    else {   allGZ = false;  }
                }
                if (!allGZ && !allPlainTxt) { //mixed bag of files, uh oh...
                    m->mothurOut("[ERROR]: Your files must all be in compressed .gz form or all in plain text form.  Please correct. \n"); m->control_pressed = true;
                }
                pair.push_back(forward);
                pair.push_back(reverse);
                pair.push_back(findex);
                pair.push_back(rindex);
                if (((findex != "") || (rindex != "")) && (oligosfile == "")) { m->mothurOut("[ERROR]: You need to provide an oligos file if you are going to use an index file.\n"); m->control_pressed = true;  }
                files.push_back(pair);
            }
        }
        in.close();
        
        if (allGZ) {
            gz = true;
            #if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
            #else
                processors=1;
            #endif
        }else { gz = false; }
        
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
        
        if (oligos->hasPairedBarcodes() || oligos->hasPairedPrimers()) {
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
 * checks for minor diffs @MS7_15058:1:1101:11899:1633#8/1 @MS7_15058:1:1101:11899:1633#8/2 should match
 */
bool MakeContigsCommand::checkName(FastqRead& forward, FastqRead& reverse){
    try {
        if (forward.getName() == reverse.getName()) {
            return true;
        }else {
            //if no match are the names only different by 1 and 2?
            string tempFRead = forward.getName().substr(0, forward.getName().length()-1);
            string tempRRead = reverse.getName().substr(0, reverse.getName().length()-1);
            if (tempFRead == tempRRead) {
                if ((forward.getName()[forward.getName().length()-1] == '1') && (reverse.getName()[reverse.getName().length()-1] == '2')) {
                    forward.setName(tempFRead);
                    reverse.setName(tempRRead);
                    return true;
                }
            }else {
                //if no match are the names only different by 1 and 2?
                string tempFRead = forward.getName();
                string tempRRead = reverse.getName().substr(0, reverse.getName().length()-1);
                if (tempFRead == tempRRead) {
                    reverse.setName(tempRRead);
                    return true;
                }
            }
        }
        
        return false;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "ckeckName");
        exit(1);
    }
}
//***************************************************************************************************************
/**
 * checks for minor diffs @MS7_15058:1:1101:11899:1633#8/1 @MS7_15058:1:1101:11899:1633#8/2 should match
 */
bool MakeContigsCommand::checkName(Sequence& forward, Sequence& reverse){
    try {
        if (forward.getName() == reverse.getName()) {
            return true;
        }else {
            //if no match are the names only different by 1 and 2?
            string tempFRead = forward.getName().substr(0, forward.getName().length()-1);
            string tempRRead = reverse.getName().substr(0, reverse.getName().length()-1);
            if (tempFRead == tempRRead) {
                if ((forward.getName()[forward.getName().length()-1] == '1') && (reverse.getName()[reverse.getName().length()-1] == '2')) {
                    forward.setName(tempFRead);
                    reverse.setName(tempRRead);
                    return true;
                }
            }else {
                //if no match are the names only different by 1 and 2?
                string tempFRead = forward.getName();
                string tempRRead = reverse.getName().substr(0, reverse.getName().length()-1);
                if (tempFRead == tempRRead) {
                    reverse.setName(tempRRead);
                    return true;
                }
            }
        }
        
        return false;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "ckeckName");
        exit(1);
    }
}
//***************************************************************************************************************
/**
 * checks for minor diffs @MS7_15058:1:1101:11899:1633#8/1 @MS7_15058:1:1101:11899:1633#8/2 should match
 */
bool MakeContigsCommand::checkName(QualityScores& forward, QualityScores& reverse){
    try {
        if (forward.getName() == reverse.getName()) {
            return true;
        }else {
            //if no match are the names only different by 1 and 2?
            string tempFRead = forward.getName().substr(0, forward.getName().length()-1);
            string tempRRead = reverse.getName().substr(0, reverse.getName().length()-1);
            if (tempFRead == tempRRead) {
                if ((forward.getName()[forward.getName().length()-1] == '1') && (reverse.getName()[reverse.getName().length()-1] == '2')) {
                    forward.setName(tempFRead);
                    reverse.setName(tempRRead);
                    return true;
                }
            }else {
                //if no match are the names only different by 1 and 2?
                string tempFRead = forward.getName();
                string tempRRead = reverse.getName().substr(0, reverse.getName().length()-1);
                if (tempFRead == tempRRead) {
                    reverse.setName(tempRRead);
                    return true;
                }
            }
        }
        
        return false;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "ckeckName");
        exit(1);
    }
}

//***************************************************************************************************************
/**
 * checks for minor diffs @MS7_15058:1:1101:11899:1633#8/1 @MS7_15058:1:1101:11899:1633#8/2 should match
 */
bool MakeContigsCommand::checkName(Sequence& forward, QualityScores& reverse){
    try {
        if (forward.getName() == reverse.getName()) {
            return true;
        }else {
            //if no match are the names only different by 1 and 2?
            string tempFRead = forward.getName().substr(0, forward.getName().length()-1);
            string tempRRead = reverse.getName().substr(0, reverse.getName().length()-1);
            if (tempFRead == tempRRead) {
                if ((forward.getName()[forward.getName().length()-1] == '1') && (reverse.getName()[reverse.getName().length()-1] == '2')) {
                    forward.setName(tempFRead);
                    reverse.setName(tempRRead);
                    return true;
                }
            }else {
                //if no match are the names only different by 1 and 2?
                string tempFRead = forward.getName();
                string tempRRead = reverse.getName().substr(0, reverse.getName().length()-1);
                if (tempFRead == tempRRead) {
                    reverse.setName(tempRRead);
                    return true;
                }
            }
        }
        
        return false;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "ckeckName");
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
//***************************************************************************************************************
int MakeContigsCommand::loadQmatchValues(vector< vector<double> >& qual_match_simple_bayesian, vector< vector<double> >& qual_mismatch_simple_bayesian){
    try {
        qual_match_simple_bayesian[0][0] = -1.09861; qual_match_simple_bayesian[0][1] = -1.32887; qual_match_simple_bayesian[0][2] = -1.55913; qual_match_simple_bayesian[0][3] = -1.78939; qual_match_simple_bayesian[0][4] = -2.01965; qual_match_simple_bayesian[0][5] = -2.2499; qual_match_simple_bayesian[0][6] = -2.48016; qual_match_simple_bayesian[0][7] = -2.71042; qual_match_simple_bayesian[0][8] = -2.94068; qual_match_simple_bayesian[0][9] = -3.17094; qual_match_simple_bayesian[0][10] = -3.4012; qual_match_simple_bayesian[0][11] = -3.63146; qual_match_simple_bayesian[0][12] = -3.86171; qual_match_simple_bayesian[0][13] = -4.09197; qual_match_simple_bayesian[0][14] = -4.32223; qual_match_simple_bayesian[0][15] = -4.55249; qual_match_simple_bayesian[0][16] = -4.78275; qual_match_simple_bayesian[0][17] = -5.01301; qual_match_simple_bayesian[0][18] = -5.24327; qual_match_simple_bayesian[0][19] = -5.47352; qual_match_simple_bayesian[0][20] = -5.70378; qual_match_simple_bayesian[0][21] = -5.93404; qual_match_simple_bayesian[0][22] = -6.1643; qual_match_simple_bayesian[0][23] = -6.39456; qual_match_simple_bayesian[0][24] = -6.62482; qual_match_simple_bayesian[0][25] = -6.85508; qual_match_simple_bayesian[0][26] = -7.08533; qual_match_simple_bayesian[0][27] = -7.31559; qual_match_simple_bayesian[0][28] = -7.54585; qual_match_simple_bayesian[0][29] = -7.77611; qual_match_simple_bayesian[0][30] = -8.00637; qual_match_simple_bayesian[0][31] = -8.23663; qual_match_simple_bayesian[0][32] = -8.46688; qual_match_simple_bayesian[0][33] = -8.69714; qual_match_simple_bayesian[0][34] = -8.9274; qual_match_simple_bayesian[0][35] = -9.15766; qual_match_simple_bayesian[0][36] = -9.38792; qual_match_simple_bayesian[0][37] = -9.61818; qual_match_simple_bayesian[0][38] = -9.84844; qual_match_simple_bayesian[0][39] = -10.0787; qual_match_simple_bayesian[0][40] = -10.309; qual_match_simple_bayesian[0][41] = -10.5392; qual_match_simple_bayesian[0][42] = -10.7695; qual_match_simple_bayesian[0][43] = -10.9997; qual_match_simple_bayesian[0][44] = -11.23; qual_match_simple_bayesian[0][45] = -11.4602; qual_match_simple_bayesian[0][46] = -11.6905;
        qual_match_simple_bayesian[1][0] = -1.32887; qual_match_simple_bayesian[1][1] = -1.37587; qual_match_simple_bayesian[1][2] = -1.41484; qual_match_simple_bayesian[1][3] = -1.44692; qual_match_simple_bayesian[1][4] = -1.47315; qual_match_simple_bayesian[1][5] = -1.49449; qual_match_simple_bayesian[1][6] = -1.51178; qual_match_simple_bayesian[1][7] = -1.52572; qual_match_simple_bayesian[1][8] = -1.53694; qual_match_simple_bayesian[1][9] = -1.54593; qual_match_simple_bayesian[1][10] = -1.55314; qual_match_simple_bayesian[1][11] = -1.5589; qual_match_simple_bayesian[1][12] = -1.5635; qual_match_simple_bayesian[1][13] = -1.56717; qual_match_simple_bayesian[1][14] = -1.5701; qual_match_simple_bayesian[1][15] = -1.57243; qual_match_simple_bayesian[1][16] = -1.57428; qual_match_simple_bayesian[1][17] = -1.57576; qual_match_simple_bayesian[1][18] = -1.57693; qual_match_simple_bayesian[1][19] = -1.57786; qual_match_simple_bayesian[1][20] = -1.5786; qual_match_simple_bayesian[1][21] = -1.57919; qual_match_simple_bayesian[1][22] = -1.57966; qual_match_simple_bayesian[1][23] = -1.58003; qual_match_simple_bayesian[1][24] = -1.58033; qual_match_simple_bayesian[1][25] = -1.58057; qual_match_simple_bayesian[1][26] = -1.58075; qual_match_simple_bayesian[1][27] = -1.5809; qual_match_simple_bayesian[1][28] = -1.58102; qual_match_simple_bayesian[1][29] = -1.58111; qual_match_simple_bayesian[1][30] = -1.58119; qual_match_simple_bayesian[1][31] = -1.58125; qual_match_simple_bayesian[1][32] = -1.58129; qual_match_simple_bayesian[1][33] = -1.58133; qual_match_simple_bayesian[1][34] = -1.58136; qual_match_simple_bayesian[1][35] = -1.58138; qual_match_simple_bayesian[1][36] = -1.5814; qual_match_simple_bayesian[1][37] = -1.58142; qual_match_simple_bayesian[1][38] = -1.58143; qual_match_simple_bayesian[1][39] = -1.58144; qual_match_simple_bayesian[1][40] = -1.58145; qual_match_simple_bayesian[1][41] = -1.58145; qual_match_simple_bayesian[1][42] = -1.58146; qual_match_simple_bayesian[1][43] = -1.58146; qual_match_simple_bayesian[1][44] = -1.58146; qual_match_simple_bayesian[1][45] = -1.58146; qual_match_simple_bayesian[1][46] = -1.58147;
        qual_match_simple_bayesian[2][0] = -1.55913; qual_match_simple_bayesian[2][1] = -1.41484; qual_match_simple_bayesian[2][2] = -1.31343; qual_match_simple_bayesian[2][3] = -1.23963; qual_match_simple_bayesian[2][4] = -1.18465; qual_match_simple_bayesian[2][5] = -1.14303; qual_match_simple_bayesian[2][6] = -1.11117; qual_match_simple_bayesian[2][7] = -1.08657; qual_match_simple_bayesian[2][8] = -1.06744; qual_match_simple_bayesian[2][9] = -1.05251; qual_match_simple_bayesian[2][10] = -1.0408; qual_match_simple_bayesian[2][11] = -1.0316; qual_match_simple_bayesian[2][12] = -1.02436; qual_match_simple_bayesian[2][13] = -1.01863; qual_match_simple_bayesian[2][14] = -1.01411; qual_match_simple_bayesian[2][15] = -1.01054; qual_match_simple_bayesian[2][16] = -1.00771; qual_match_simple_bayesian[2][17] = -1.00546; qual_match_simple_bayesian[2][18] = -1.00368; qual_match_simple_bayesian[2][19] = -1.00227; qual_match_simple_bayesian[2][20] = -1.00115; qual_match_simple_bayesian[2][21] = -1.00027; qual_match_simple_bayesian[2][22] = -0.99956; qual_match_simple_bayesian[2][23] = -0.999001; qual_match_simple_bayesian[2][24] = -0.998557; qual_match_simple_bayesian[2][25] = -0.998204; qual_match_simple_bayesian[2][26] = -0.997924; qual_match_simple_bayesian[2][27] = -0.997702; qual_match_simple_bayesian[2][28] = -0.997525; qual_match_simple_bayesian[2][29] = -0.997385; qual_match_simple_bayesian[2][30] = -0.997273; qual_match_simple_bayesian[2][31] = -0.997185; qual_match_simple_bayesian[2][32] = -0.997114; qual_match_simple_bayesian[2][33] = -0.997059; qual_match_simple_bayesian[2][34] = -0.997014; qual_match_simple_bayesian[2][35] = -0.996979; qual_match_simple_bayesian[2][36] = -0.996951; qual_match_simple_bayesian[2][37] = -0.996929; qual_match_simple_bayesian[2][38] = -0.996911; qual_match_simple_bayesian[2][39] = -0.996897; qual_match_simple_bayesian[2][40] = -0.996886; qual_match_simple_bayesian[2][41] = -0.996877; qual_match_simple_bayesian[2][42] = -0.99687; qual_match_simple_bayesian[2][43] = -0.996865; qual_match_simple_bayesian[2][44] = -0.99686; qual_match_simple_bayesian[2][45] = -0.996857; qual_match_simple_bayesian[2][46] = -0.996854;
        qual_match_simple_bayesian[3][0] = -1.78939; qual_match_simple_bayesian[3][1] = -1.44692; qual_match_simple_bayesian[3][2] = -1.23963; qual_match_simple_bayesian[3][3] = -1.10098; qual_match_simple_bayesian[3][4] = -1.0031; qual_match_simple_bayesian[3][5] = -0.931648; qual_match_simple_bayesian[3][6] = -0.878319; qual_match_simple_bayesian[3][7] = -0.837896; qual_match_simple_bayesian[3][8] = -0.806912; qual_match_simple_bayesian[3][9] = -0.782967; qual_match_simple_bayesian[3][10] = -0.764347; qual_match_simple_bayesian[3][11] = -0.7498; qual_match_simple_bayesian[3][12] = -0.738394; qual_match_simple_bayesian[3][13] = -0.729426; qual_match_simple_bayesian[3][14] = -0.722359; qual_match_simple_bayesian[3][15] = -0.71678; qual_match_simple_bayesian[3][16] = -0.712372; qual_match_simple_bayesian[3][17] = -0.708883; qual_match_simple_bayesian[3][18] = -0.706121; qual_match_simple_bayesian[3][19] = -0.703933; qual_match_simple_bayesian[3][20] = -0.702197; qual_match_simple_bayesian[3][21] = -0.700821; qual_match_simple_bayesian[3][22] = -0.69973; qual_match_simple_bayesian[3][23] = -0.698863; qual_match_simple_bayesian[3][24] = -0.698176; qual_match_simple_bayesian[3][25] = -0.69763; qual_match_simple_bayesian[3][26] = -0.697196; qual_match_simple_bayesian[3][27] = -0.696852; qual_match_simple_bayesian[3][28] = -0.696579; qual_match_simple_bayesian[3][29] = -0.696362; qual_match_simple_bayesian[3][30] = -0.69619; qual_match_simple_bayesian[3][31] = -0.696053; qual_match_simple_bayesian[3][32] = -0.695944; qual_match_simple_bayesian[3][33] = -0.695858; qual_match_simple_bayesian[3][34] = -0.695789; qual_match_simple_bayesian[3][35] = -0.695735; qual_match_simple_bayesian[3][36] = -0.695692; qual_match_simple_bayesian[3][37] = -0.695657; qual_match_simple_bayesian[3][38] = -0.69563; qual_match_simple_bayesian[3][39] = -0.695608; qual_match_simple_bayesian[3][40] = -0.695591; qual_match_simple_bayesian[3][41] = -0.695577; qual_match_simple_bayesian[3][42] = -0.695566; qual_match_simple_bayesian[3][43] = -0.695558; qual_match_simple_bayesian[3][44] = -0.695551; qual_match_simple_bayesian[3][45] = -0.695546; qual_match_simple_bayesian[3][46] = -0.695541;
        qual_match_simple_bayesian[4][0] = -2.01965; qual_match_simple_bayesian[4][1] = -1.47315; qual_match_simple_bayesian[4][2] = -1.18465; qual_match_simple_bayesian[4][3] = -1.0031; qual_match_simple_bayesian[4][4] = -0.879224; qual_match_simple_bayesian[4][5] = -0.790712; qual_match_simple_bayesian[4][6] = -0.725593; qual_match_simple_bayesian[4][7] = -0.676729; qual_match_simple_bayesian[4][8] = -0.639547; qual_match_simple_bayesian[4][9] = -0.610968; qual_match_simple_bayesian[4][10] = -0.588834; qual_match_simple_bayesian[4][11] = -0.571596; qual_match_simple_bayesian[4][12] = -0.558111; qual_match_simple_bayesian[4][13] = -0.547528; qual_match_simple_bayesian[4][14] = -0.539201; qual_match_simple_bayesian[4][15] = -0.532636; qual_match_simple_bayesian[4][16] = -0.527451; qual_match_simple_bayesian[4][17] = -0.523352; qual_match_simple_bayesian[4][18] = -0.520107; qual_match_simple_bayesian[4][19] = -0.517538; qual_match_simple_bayesian[4][20] = -0.515502; qual_match_simple_bayesian[4][21] = -0.513887; qual_match_simple_bayesian[4][22] = -0.512606; qual_match_simple_bayesian[4][23] = -0.51159; qual_match_simple_bayesian[4][24] = -0.510784; qual_match_simple_bayesian[4][25] = -0.510144; qual_match_simple_bayesian[4][26] = -0.509636; qual_match_simple_bayesian[4][27] = -0.509232; qual_match_simple_bayesian[4][28] = -0.508912; qual_match_simple_bayesian[4][29] = -0.508658; qual_match_simple_bayesian[4][30] = -0.508456; qual_match_simple_bayesian[4][31] = -0.508295; qual_match_simple_bayesian[4][32] = -0.508168; qual_match_simple_bayesian[4][33] = -0.508067; qual_match_simple_bayesian[4][34] = -0.507986; qual_match_simple_bayesian[4][35] = -0.507922; qual_match_simple_bayesian[4][36] = -0.507872; qual_match_simple_bayesian[4][37] = -0.507831; qual_match_simple_bayesian[4][38] = -0.507799; qual_match_simple_bayesian[4][39] = -0.507774; qual_match_simple_bayesian[4][40] = -0.507754; qual_match_simple_bayesian[4][41] = -0.507738; qual_match_simple_bayesian[4][42] = -0.507725; qual_match_simple_bayesian[4][43] = -0.507715; qual_match_simple_bayesian[4][44] = -0.507707; qual_match_simple_bayesian[4][45] = -0.507701; qual_match_simple_bayesian[4][46] = -0.507695;
        qual_match_simple_bayesian[5][0] = -2.2499; qual_match_simple_bayesian[5][1] = -1.49449; qual_match_simple_bayesian[5][2] = -1.14303; qual_match_simple_bayesian[5][3] = -0.931648; qual_match_simple_bayesian[5][4] = -0.790712; qual_match_simple_bayesian[5][5] = -0.691393; qual_match_simple_bayesian[5][6] = -0.618979; qual_match_simple_bayesian[5][7] = -0.564976; qual_match_simple_bayesian[5][8] = -0.524066; qual_match_simple_bayesian[5][9] = -0.492723; qual_match_simple_bayesian[5][10] = -0.468507; qual_match_simple_bayesian[5][11] = -0.449682; qual_match_simple_bayesian[5][12] = -0.434976; qual_match_simple_bayesian[5][13] = -0.423448; qual_match_simple_bayesian[5][14] = -0.414384; qual_match_simple_bayesian[5][15] = -0.407243; qual_match_simple_bayesian[5][16] = -0.401606; qual_match_simple_bayesian[5][17] = -0.397151; qual_match_simple_bayesian[5][18] = -0.393627; qual_match_simple_bayesian[5][19] = -0.390836; qual_match_simple_bayesian[5][20] = -0.388625; qual_match_simple_bayesian[5][21] = -0.386872; qual_match_simple_bayesian[5][22] = -0.385482; qual_match_simple_bayesian[5][23] = -0.384379; qual_match_simple_bayesian[5][24] = -0.383503; qual_match_simple_bayesian[5][25] = -0.382809; qual_match_simple_bayesian[5][26] = -0.382257; qual_match_simple_bayesian[5][27] = -0.38182; qual_match_simple_bayesian[5][28] = -0.381472; qual_match_simple_bayesian[5][29] = -0.381196; qual_match_simple_bayesian[5][30] = -0.380977; qual_match_simple_bayesian[5][31] = -0.380803; qual_match_simple_bayesian[5][32] = -0.380664; qual_match_simple_bayesian[5][33] = -0.380554; qual_match_simple_bayesian[5][34] = -0.380467; qual_match_simple_bayesian[5][35] = -0.380398; qual_match_simple_bayesian[5][36] = -0.380343; qual_match_simple_bayesian[5][37] = -0.380299; qual_match_simple_bayesian[5][38] = -0.380264; qual_match_simple_bayesian[5][39] = -0.380237; qual_match_simple_bayesian[5][40] = -0.380215; qual_match_simple_bayesian[5][41] = -0.380198; qual_match_simple_bayesian[5][42] = -0.380184; qual_match_simple_bayesian[5][43] = -0.380173; qual_match_simple_bayesian[5][44] = -0.380164; qual_match_simple_bayesian[5][45] = -0.380157; qual_match_simple_bayesian[5][46] = -0.380152;
        qual_match_simple_bayesian[6][0] = -2.48016; qual_match_simple_bayesian[6][1] = -1.51178; qual_match_simple_bayesian[6][2] = -1.11117; qual_match_simple_bayesian[6][3] = -0.878319; qual_match_simple_bayesian[6][4] = -0.725593; qual_match_simple_bayesian[6][5] = -0.618979; qual_match_simple_bayesian[6][6] = -0.541714; qual_match_simple_bayesian[6][7] = -0.48433; qual_match_simple_bayesian[6][8] = -0.440984; qual_match_simple_bayesian[6][9] = -0.407844; qual_match_simple_bayesian[6][10] = -0.382281; qual_match_simple_bayesian[6][11] = -0.362431; qual_match_simple_bayesian[6][12] = -0.34694; qual_match_simple_bayesian[6][13] = -0.334804; qual_match_simple_bayesian[6][14] = -0.325268; qual_match_simple_bayesian[6][15] = -0.317757; qual_match_simple_bayesian[6][16] = -0.311831; qual_match_simple_bayesian[6][17] = -0.307149; qual_match_simple_bayesian[6][18] = -0.303445; qual_match_simple_bayesian[6][19] = -0.300513; qual_match_simple_bayesian[6][20] = -0.29819; qual_match_simple_bayesian[6][21] = -0.296348; qual_match_simple_bayesian[6][22] = -0.294888; qual_match_simple_bayesian[6][23] = -0.29373; qual_match_simple_bayesian[6][24] = -0.29281; qual_match_simple_bayesian[6][25] = -0.292081; qual_match_simple_bayesian[6][26] = -0.291502; qual_match_simple_bayesian[6][27] = -0.291042; qual_match_simple_bayesian[6][28] = -0.290677; qual_match_simple_bayesian[6][29] = -0.290387; qual_match_simple_bayesian[6][30] = -0.290157; qual_match_simple_bayesian[6][31] = -0.289974; qual_match_simple_bayesian[6][32] = -0.289829; qual_match_simple_bayesian[6][33] = -0.289713; qual_match_simple_bayesian[6][34] = -0.289622; qual_match_simple_bayesian[6][35] = -0.289549; qual_match_simple_bayesian[6][36] = -0.289491; qual_match_simple_bayesian[6][37] = -0.289445; qual_match_simple_bayesian[6][38] = -0.289409; qual_match_simple_bayesian[6][39] = -0.28938; qual_match_simple_bayesian[6][40] = -0.289357; qual_match_simple_bayesian[6][41] = -0.289339; qual_match_simple_bayesian[6][42] = -0.289324; qual_match_simple_bayesian[6][43] = -0.289313; qual_match_simple_bayesian[6][44] = -0.289304; qual_match_simple_bayesian[6][45] = -0.289296; qual_match_simple_bayesian[6][46] = -0.28929;
        qual_match_simple_bayesian[7][0] = -2.71042; qual_match_simple_bayesian[7][1] = -1.52572; qual_match_simple_bayesian[7][2] = -1.08657; qual_match_simple_bayesian[7][3] = -0.837896; qual_match_simple_bayesian[7][4] = -0.676729; qual_match_simple_bayesian[7][5] = -0.564976; qual_match_simple_bayesian[7][6] = -0.48433; qual_match_simple_bayesian[7][7] = -0.424604; qual_match_simple_bayesian[7][8] = -0.379581; qual_match_simple_bayesian[7][9] = -0.345208; qual_match_simple_bayesian[7][10] = -0.318723; qual_match_simple_bayesian[7][11] = -0.298173; qual_match_simple_bayesian[7][12] = -0.282146; qual_match_simple_bayesian[7][13] = -0.269595; qual_match_simple_bayesian[7][14] = -0.259737; qual_match_simple_bayesian[7][15] = -0.251976; qual_match_simple_bayesian[7][16] = -0.245853; qual_match_simple_bayesian[7][17] = -0.241016; qual_match_simple_bayesian[7][18] = -0.23719; qual_match_simple_bayesian[7][19] = -0.234162; qual_match_simple_bayesian[7][20] = -0.231763; qual_match_simple_bayesian[7][21] = -0.229861; qual_match_simple_bayesian[7][22] = -0.228354; qual_match_simple_bayesian[7][23] = -0.227158; qual_match_simple_bayesian[7][24] = -0.226208; qual_match_simple_bayesian[7][25] = -0.225455; qual_match_simple_bayesian[7][26] = -0.224857; qual_match_simple_bayesian[7][27] = -0.224383; qual_match_simple_bayesian[7][28] = -0.224006; qual_match_simple_bayesian[7][29] = -0.223707; qual_match_simple_bayesian[7][30] = -0.223469; qual_match_simple_bayesian[7][31] = -0.22328; qual_match_simple_bayesian[7][32] = -0.22313; qual_match_simple_bayesian[7][33] = -0.223011; qual_match_simple_bayesian[7][34] = -0.222917; qual_match_simple_bayesian[7][35] = -0.222842; qual_match_simple_bayesian[7][36] = -0.222782; qual_match_simple_bayesian[7][37] = -0.222734; qual_match_simple_bayesian[7][38] = -0.222697; qual_match_simple_bayesian[7][39] = -0.222667; qual_match_simple_bayesian[7][40] = -0.222643; qual_match_simple_bayesian[7][41] = -0.222624; qual_match_simple_bayesian[7][42] = -0.222609; qual_match_simple_bayesian[7][43] = -0.222597; qual_match_simple_bayesian[7][44] = -0.222588; qual_match_simple_bayesian[7][45] = -0.222581; qual_match_simple_bayesian[7][46] = -0.222575;
        qual_match_simple_bayesian[8][0] = -2.94068; qual_match_simple_bayesian[8][1] = -1.53694; qual_match_simple_bayesian[8][2] = -1.06744; qual_match_simple_bayesian[8][3] = -0.806912; qual_match_simple_bayesian[8][4] = -0.639547; qual_match_simple_bayesian[8][5] = -0.524066; qual_match_simple_bayesian[8][6] = -0.440984; qual_match_simple_bayesian[8][7] = -0.379581; qual_match_simple_bayesian[8][8] = -0.333359; qual_match_simple_bayesian[8][9] = -0.298107; qual_match_simple_bayesian[8][10] = -0.270966; qual_match_simple_bayesian[8][11] = -0.249919; qual_match_simple_bayesian[8][12] = -0.233512; qual_match_simple_bayesian[8][13] = -0.220668; qual_match_simple_bayesian[8][14] = -0.210582; qual_match_simple_bayesian[8][15] = -0.202642; qual_match_simple_bayesian[8][16] = -0.19638; qual_match_simple_bayesian[8][17] = -0.191434; qual_match_simple_bayesian[8][18] = -0.187522; qual_match_simple_bayesian[8][19] = -0.184426; qual_match_simple_bayesian[8][20] = -0.181973; qual_match_simple_bayesian[8][21] = -0.180029; qual_match_simple_bayesian[8][22] = -0.178488; qual_match_simple_bayesian[8][23] = -0.177265; qual_match_simple_bayesian[8][24] = -0.176295; qual_match_simple_bayesian[8][25] = -0.175525; qual_match_simple_bayesian[8][26] = -0.174914; qual_match_simple_bayesian[8][27] = -0.174428; qual_match_simple_bayesian[8][28] = -0.174043; qual_match_simple_bayesian[8][29] = -0.173737; qual_match_simple_bayesian[8][30] = -0.173494; qual_match_simple_bayesian[8][31] = -0.173301; qual_match_simple_bayesian[8][32] = -0.173148; qual_match_simple_bayesian[8][33] = -0.173026; qual_match_simple_bayesian[8][34] = -0.17293; qual_match_simple_bayesian[8][35] = -0.172853; qual_match_simple_bayesian[8][36] = -0.172792; qual_match_simple_bayesian[8][37] = -0.172744; qual_match_simple_bayesian[8][38] = -0.172705; qual_match_simple_bayesian[8][39] = -0.172675; qual_match_simple_bayesian[8][40] = -0.17265; qual_match_simple_bayesian[8][41] = -0.172631; qual_match_simple_bayesian[8][42] = -0.172616; qual_match_simple_bayesian[8][43] = -0.172604; qual_match_simple_bayesian[8][44] = -0.172594; qual_match_simple_bayesian[8][45] = -0.172586; qual_match_simple_bayesian[8][46] = -0.17258;
        qual_match_simple_bayesian[9][0] = -3.17094; qual_match_simple_bayesian[9][1] = -1.54593; qual_match_simple_bayesian[9][2] = -1.05251; qual_match_simple_bayesian[9][3] = -0.782967; qual_match_simple_bayesian[9][4] = -0.610968; qual_match_simple_bayesian[9][5] = -0.492723; qual_match_simple_bayesian[9][6] = -0.407844; qual_match_simple_bayesian[9][7] = -0.345208; qual_match_simple_bayesian[9][8] = -0.298107; qual_match_simple_bayesian[9][9] = -0.262213; qual_match_simple_bayesian[9][10] = -0.234592; qual_match_simple_bayesian[9][11] = -0.213183; qual_match_simple_bayesian[9][12] = -0.196498; qual_match_simple_bayesian[9][13] = -0.18344; qual_match_simple_bayesian[9][14] = -0.173188; qual_match_simple_bayesian[9][15] = -0.165119; qual_match_simple_bayesian[9][16] = -0.158755; qual_match_simple_bayesian[9][17] = -0.153729; qual_match_simple_bayesian[9][18] = -0.149755; qual_match_simple_bayesian[9][19] = -0.146609; qual_match_simple_bayesian[9][20] = -0.144117; qual_match_simple_bayesian[9][21] = -0.142143; qual_match_simple_bayesian[9][22] = -0.140577; qual_match_simple_bayesian[9][23] = -0.139335; qual_match_simple_bayesian[9][24] = -0.138349; qual_match_simple_bayesian[9][25] = -0.137567; qual_match_simple_bayesian[9][26] = -0.136946; qual_match_simple_bayesian[9][27] = -0.136453; qual_match_simple_bayesian[9][28] = -0.136062; qual_match_simple_bayesian[9][29] = -0.135751; qual_match_simple_bayesian[9][30] = -0.135504; qual_match_simple_bayesian[9][31] = -0.135308; qual_match_simple_bayesian[9][32] = -0.135153; qual_match_simple_bayesian[9][33] = -0.135029; qual_match_simple_bayesian[9][34] = -0.134931; qual_match_simple_bayesian[9][35] = -0.134853; qual_match_simple_bayesian[9][36] = -0.134791; qual_match_simple_bayesian[9][37] = -0.134742; qual_match_simple_bayesian[9][38] = -0.134703; qual_match_simple_bayesian[9][39] = -0.134672; qual_match_simple_bayesian[9][40] = -0.134647; qual_match_simple_bayesian[9][41] = -0.134628; qual_match_simple_bayesian[9][42] = -0.134612; qual_match_simple_bayesian[9][43] = -0.1346; qual_match_simple_bayesian[9][44] = -0.13459; qual_match_simple_bayesian[9][45] = -0.134582; qual_match_simple_bayesian[9][46] = -0.134576;
        qual_match_simple_bayesian[10][0] = -3.4012; qual_match_simple_bayesian[10][1] = -1.55314; qual_match_simple_bayesian[10][2] = -1.0408; qual_match_simple_bayesian[10][3] = -0.764347; qual_match_simple_bayesian[10][4] = -0.588834; qual_match_simple_bayesian[10][5] = -0.468507; qual_match_simple_bayesian[10][6] = -0.382281; qual_match_simple_bayesian[10][7] = -0.318723; qual_match_simple_bayesian[10][8] = -0.270966; qual_match_simple_bayesian[10][9] = -0.234592; qual_match_simple_bayesian[10][10] = -0.206614; qual_match_simple_bayesian[10][11] = -0.184935; qual_match_simple_bayesian[10][12] = -0.168044; qual_match_simple_bayesian[10][13] = -0.154827; qual_match_simple_bayesian[10][14] = -0.144451; qual_match_simple_bayesian[10][15] = -0.136285; qual_match_simple_bayesian[10][16] = -0.129846; qual_match_simple_bayesian[10][17] = -0.124761; qual_match_simple_bayesian[10][18] = -0.12074; qual_match_simple_bayesian[10][19] = -0.117558; qual_match_simple_bayesian[10][20] = -0.115037; qual_match_simple_bayesian[10][21] = -0.113039; qual_match_simple_bayesian[10][22] = -0.111455; qual_match_simple_bayesian[10][23] = -0.110198; qual_match_simple_bayesian[10][24] = -0.109202; qual_match_simple_bayesian[10][25] = -0.10841; qual_match_simple_bayesian[10][26] = -0.107782; qual_match_simple_bayesian[10][27] = -0.107284; qual_match_simple_bayesian[10][28] = -0.106888; qual_match_simple_bayesian[10][29] = -0.106574; qual_match_simple_bayesian[10][30] = -0.106324; qual_match_simple_bayesian[10][31] = -0.106126; qual_match_simple_bayesian[10][32] = -0.105968; qual_match_simple_bayesian[10][33] = -0.105843; qual_match_simple_bayesian[10][34] = -0.105744; qual_match_simple_bayesian[10][35] = -0.105665; qual_match_simple_bayesian[10][36] = -0.105602; qual_match_simple_bayesian[10][37] = -0.105553; qual_match_simple_bayesian[10][38] = -0.105513; qual_match_simple_bayesian[10][39] = -0.105482; qual_match_simple_bayesian[10][40] = -0.105457; qual_match_simple_bayesian[10][41] = -0.105437; qual_match_simple_bayesian[10][42] = -0.105421; qual_match_simple_bayesian[10][43] = -0.105409; qual_match_simple_bayesian[10][44] = -0.105399; qual_match_simple_bayesian[10][45] = -0.105391; qual_match_simple_bayesian[10][46] = -0.105385;
        qual_match_simple_bayesian[11][0] = -3.63146; qual_match_simple_bayesian[11][1] = -1.5589; qual_match_simple_bayesian[11][2] = -1.0316; qual_match_simple_bayesian[11][3] = -0.7498; qual_match_simple_bayesian[11][4] = -0.571596; qual_match_simple_bayesian[11][5] = -0.449682; qual_match_simple_bayesian[11][6] = -0.362431; qual_match_simple_bayesian[11][7] = -0.298173; qual_match_simple_bayesian[11][8] = -0.249919; qual_match_simple_bayesian[11][9] = -0.213183; qual_match_simple_bayesian[11][10] = -0.184935; qual_match_simple_bayesian[11][11] = -0.163052; qual_match_simple_bayesian[11][12] = -0.146004; qual_match_simple_bayesian[11][13] = -0.132667; qual_match_simple_bayesian[11][14] = -0.122198; qual_match_simple_bayesian[11][15] = -0.11396; qual_match_simple_bayesian[11][16] = -0.107464; qual_match_simple_bayesian[11][17] = -0.102334; qual_match_simple_bayesian[11][18] = -0.0982781; qual_match_simple_bayesian[11][19] = -0.0950678; qual_match_simple_bayesian[11][20] = -0.0925252; qual_match_simple_bayesian[11][21] = -0.09051; qual_match_simple_bayesian[11][22] = -0.0889123; qual_match_simple_bayesian[11][23] = -0.0876449; qual_match_simple_bayesian[11][24] = -0.0866394; qual_match_simple_bayesian[11][25] = -0.0858414; qual_match_simple_bayesian[11][26] = -0.0852079; qual_match_simple_bayesian[11][27] = -0.0847051; qual_match_simple_bayesian[11][28] = -0.0843058; qual_match_simple_bayesian[11][29] = -0.0839888; qual_match_simple_bayesian[11][30] = -0.083737; qual_match_simple_bayesian[11][31] = -0.0835371; qual_match_simple_bayesian[11][32] = -0.0833783; qual_match_simple_bayesian[11][33] = -0.0832522; qual_match_simple_bayesian[11][34] = -0.083152; qual_match_simple_bayesian[11][35] = -0.0830725; qual_match_simple_bayesian[11][36] = -0.0830093; qual_match_simple_bayesian[11][37] = -0.0829591; qual_match_simple_bayesian[11][38] = -0.0829192; qual_match_simple_bayesian[11][39] = -0.0828876; qual_match_simple_bayesian[11][40] = -0.0828624; qual_match_simple_bayesian[11][41] = -0.0828425; qual_match_simple_bayesian[11][42] = -0.0828266; qual_match_simple_bayesian[11][43] = -0.082814; qual_match_simple_bayesian[11][44] = -0.082804; qual_match_simple_bayesian[11][45] = -0.082796; qual_match_simple_bayesian[11][46] = -0.0827897;
        qual_match_simple_bayesian[12][0] = -3.86171; qual_match_simple_bayesian[12][1] = -1.5635; qual_match_simple_bayesian[12][2] = -1.02436; qual_match_simple_bayesian[12][3] = -0.738394; qual_match_simple_bayesian[12][4] = -0.558111; qual_match_simple_bayesian[12][5] = -0.434976; qual_match_simple_bayesian[12][6] = -0.34694; qual_match_simple_bayesian[12][7] = -0.282146; qual_match_simple_bayesian[12][8] = -0.233512; qual_match_simple_bayesian[12][9] = -0.196498; qual_match_simple_bayesian[12][10] = -0.168044; qual_match_simple_bayesian[12][11] = -0.146004; qual_match_simple_bayesian[12][12] = -0.128838; qual_match_simple_bayesian[12][13] = -0.115409; qual_match_simple_bayesian[12][14] = -0.104869; qual_match_simple_bayesian[12][15] = -0.096575; qual_match_simple_bayesian[12][16] = -0.0900357; qual_match_simple_bayesian[12][17] = -0.0848716; qual_match_simple_bayesian[12][18] = -0.0807886; qual_match_simple_bayesian[12][19] = -0.0775572; qual_match_simple_bayesian[12][20] = -0.0749978; qual_match_simple_bayesian[12][21] = -0.0729694; qual_match_simple_bayesian[12][22] = -0.0713612; qual_match_simple_bayesian[12][23] = -0.0700856; qual_match_simple_bayesian[12][24] = -0.0690735; qual_match_simple_bayesian[12][25] = -0.0682703; qual_match_simple_bayesian[12][26] = -0.0676327; qual_match_simple_bayesian[12][27] = -0.0671265; qual_match_simple_bayesian[12][28] = -0.0667247; qual_match_simple_bayesian[12][29] = -0.0664056; qual_match_simple_bayesian[12][30] = -0.0661522; qual_match_simple_bayesian[12][31] = -0.065951; qual_match_simple_bayesian[12][32] = -0.0657912; qual_match_simple_bayesian[12][33] = -0.0656642; qual_match_simple_bayesian[12][34] = -0.0655634; qual_match_simple_bayesian[12][35] = -0.0654833; qual_match_simple_bayesian[12][36] = -0.0654198; qual_match_simple_bayesian[12][37] = -0.0653692; qual_match_simple_bayesian[12][38] = -0.0653291; qual_match_simple_bayesian[12][39] = -0.0652972; qual_match_simple_bayesian[12][40] = -0.0652719; qual_match_simple_bayesian[12][41] = -0.0652518; qual_match_simple_bayesian[12][42] = -0.0652359; qual_match_simple_bayesian[12][43] = -0.0652232; qual_match_simple_bayesian[12][44] = -0.0652131; qual_match_simple_bayesian[12][45] = -0.0652051; qual_match_simple_bayesian[12][46] = -0.0651987;
        qual_match_simple_bayesian[13][0] = -4.09197; qual_match_simple_bayesian[13][1] = -1.56717; qual_match_simple_bayesian[13][2] = -1.01863; qual_match_simple_bayesian[13][3] = -0.729426; qual_match_simple_bayesian[13][4] = -0.547528; qual_match_simple_bayesian[13][5] = -0.423448; qual_match_simple_bayesian[13][6] = -0.334804; qual_match_simple_bayesian[13][7] = -0.269595; qual_match_simple_bayesian[13][8] = -0.220668; qual_match_simple_bayesian[13][9] = -0.18344; qual_match_simple_bayesian[13][10] = -0.154827; qual_match_simple_bayesian[13][11] = -0.132667; qual_match_simple_bayesian[13][12] = -0.115409; qual_match_simple_bayesian[13][13] = -0.101909; qual_match_simple_bayesian[13][14] = -0.0913142; qual_match_simple_bayesian[13][15] = -0.0829777; qual_match_simple_bayesian[13][16] = -0.0764049; qual_match_simple_bayesian[13][17] = -0.0712146; qual_match_simple_bayesian[13][18] = -0.0671109; qual_match_simple_bayesian[13][19] = -0.0638632; qual_match_simple_bayesian[13][20] = -0.061291; qual_match_simple_bayesian[13][21] = -0.0592525; qual_match_simple_bayesian[13][22] = -0.0576362; qual_match_simple_bayesian[13][23] = -0.0563542; qual_match_simple_bayesian[13][24] = -0.055337; qual_match_simple_bayesian[13][25] = -0.0545298; qual_match_simple_bayesian[13][26] = -0.053889; qual_match_simple_bayesian[13][27] = -0.0533804; qual_match_simple_bayesian[13][28] = -0.0529765; qual_match_simple_bayesian[13][29] = -0.0526558; qual_match_simple_bayesian[13][30] = -0.0524012; qual_match_simple_bayesian[13][31] = -0.0521989; qual_match_simple_bayesian[13][32] = -0.0520383; qual_match_simple_bayesian[13][33] = -0.0519108; qual_match_simple_bayesian[13][34] = -0.0518095; qual_match_simple_bayesian[13][35] = -0.051729; qual_match_simple_bayesian[13][36] = -0.0516651; qual_match_simple_bayesian[13][37] = -0.0516143; qual_match_simple_bayesian[13][38] = -0.051574; qual_match_simple_bayesian[13][39] = -0.051542; qual_match_simple_bayesian[13][40] = -0.0515165; qual_match_simple_bayesian[13][41] = -0.0514963; qual_match_simple_bayesian[13][42] = -0.0514803; qual_match_simple_bayesian[13][43] = -0.0514675; qual_match_simple_bayesian[13][44] = -0.0514574; qual_match_simple_bayesian[13][45] = -0.0514493; qual_match_simple_bayesian[13][46] = -0.051443;
        qual_match_simple_bayesian[14][0] = -4.32223; qual_match_simple_bayesian[14][1] = -1.5701; qual_match_simple_bayesian[14][2] = -1.01411; qual_match_simple_bayesian[14][3] = -0.722359; qual_match_simple_bayesian[14][4] = -0.539201; qual_match_simple_bayesian[14][5] = -0.414384; qual_match_simple_bayesian[14][6] = -0.325268; qual_match_simple_bayesian[14][7] = -0.259737; qual_match_simple_bayesian[14][8] = -0.210582; qual_match_simple_bayesian[14][9] = -0.173188; qual_match_simple_bayesian[14][10] = -0.144451; qual_match_simple_bayesian[14][11] = -0.122198; qual_match_simple_bayesian[14][12] = -0.104869; qual_match_simple_bayesian[14][13] = -0.0913142; qual_match_simple_bayesian[14][14] = -0.0806768; qual_match_simple_bayesian[14][15] = -0.0723072; qual_match_simple_bayesian[14][16] = -0.0657085; qual_match_simple_bayesian[14][17] = -0.0604979; qual_match_simple_bayesian[14][18] = -0.0563782; qual_match_simple_bayesian[14][19] = -0.0531178; qual_match_simple_bayesian[14][20] = -0.0505356; qual_match_simple_bayesian[14][21] = -0.0484892; qual_match_simple_bayesian[14][22] = -0.0468667; qual_match_simple_bayesian[14][23] = -0.0455797; qual_match_simple_bayesian[14][24] = -0.0445586; qual_match_simple_bayesian[14][25] = -0.0437483; qual_match_simple_bayesian[14][26] = -0.0431051; qual_match_simple_bayesian[14][27] = -0.0425945; qual_match_simple_bayesian[14][28] = -0.0421891; qual_match_simple_bayesian[14][29] = -0.0418671; qual_match_simple_bayesian[14][30] = -0.0416115; qual_match_simple_bayesian[14][31] = -0.0414085; qual_match_simple_bayesian[14][32] = -0.0412473; qual_match_simple_bayesian[14][33] = -0.0411192; qual_match_simple_bayesian[14][34] = -0.0410175; qual_match_simple_bayesian[14][35] = -0.0409368; qual_match_simple_bayesian[14][36] = -0.0408726; qual_match_simple_bayesian[14][37] = -0.0408216; qual_match_simple_bayesian[14][38] = -0.0407812; qual_match_simple_bayesian[14][39] = -0.040749; qual_match_simple_bayesian[14][40] = -0.0407235; qual_match_simple_bayesian[14][41] = -0.0407032; qual_match_simple_bayesian[14][42] = -0.0406871; qual_match_simple_bayesian[14][43] = -0.0406743; qual_match_simple_bayesian[14][44] = -0.0406641; qual_match_simple_bayesian[14][45] = -0.040656; qual_match_simple_bayesian[14][46] = -0.0406496;
        qual_match_simple_bayesian[15][0] = -4.55249; qual_match_simple_bayesian[15][1] = -1.57243; qual_match_simple_bayesian[15][2] = -1.01054; qual_match_simple_bayesian[15][3] = -0.71678; qual_match_simple_bayesian[15][4] = -0.532636; qual_match_simple_bayesian[15][5] = -0.407243; qual_match_simple_bayesian[15][6] = -0.317757; qual_match_simple_bayesian[15][7] = -0.251976; qual_match_simple_bayesian[15][8] = -0.202642; qual_match_simple_bayesian[15][9] = -0.165119; qual_match_simple_bayesian[15][10] = -0.136285; qual_match_simple_bayesian[15][11] = -0.11396; qual_match_simple_bayesian[15][12] = -0.096575; qual_match_simple_bayesian[15][13] = -0.0829777; qual_match_simple_bayesian[15][14] = -0.0723072; qual_match_simple_bayesian[15][15] = -0.0639118; qual_match_simple_bayesian[15][16] = -0.0572929; qual_match_simple_bayesian[15][17] = -0.0520664; qual_match_simple_bayesian[15][18] = -0.0479342; qual_match_simple_bayesian[15][19] = -0.044664; qual_match_simple_bayesian[15][20] = -0.042074; qual_match_simple_bayesian[15][21] = -0.0400214; qual_match_simple_bayesian[15][22] = -0.038394; qual_match_simple_bayesian[15][23] = -0.0371032; qual_match_simple_bayesian[15][24] = -0.0360791; qual_match_simple_bayesian[15][25] = -0.0352663; qual_match_simple_bayesian[15][26] = -0.0346212; qual_match_simple_bayesian[15][27] = -0.0341091; qual_match_simple_bayesian[15][28] = -0.0337024; qual_match_simple_bayesian[15][29] = -0.0333796; qual_match_simple_bayesian[15][30] = -0.0331232; qual_match_simple_bayesian[15][31] = -0.0329196; qual_match_simple_bayesian[15][32] = -0.0327579; qual_match_simple_bayesian[15][33] = -0.0326294; qual_match_simple_bayesian[15][34] = -0.0325274; qual_match_simple_bayesian[15][35] = -0.0324464; qual_match_simple_bayesian[15][36] = -0.0323821; qual_match_simple_bayesian[15][37] = -0.0323309; qual_match_simple_bayesian[15][38] = -0.0322904; qual_match_simple_bayesian[15][39] = -0.0322581; qual_match_simple_bayesian[15][40] = -0.0322325; qual_match_simple_bayesian[15][41] = -0.0322121; qual_match_simple_bayesian[15][42] = -0.032196; qual_match_simple_bayesian[15][43] = -0.0321831; qual_match_simple_bayesian[15][44] = -0.032173; qual_match_simple_bayesian[15][45] = -0.0321649; qual_match_simple_bayesian[15][46] = -0.0321584;
        qual_match_simple_bayesian[16][0] = -4.78275; qual_match_simple_bayesian[16][1] = -1.57428; qual_match_simple_bayesian[16][2] = -1.00771; qual_match_simple_bayesian[16][3] = -0.712372; qual_match_simple_bayesian[16][4] = -0.527451; qual_match_simple_bayesian[16][5] = -0.401606; qual_match_simple_bayesian[16][6] = -0.311831; qual_match_simple_bayesian[16][7] = -0.245853; qual_match_simple_bayesian[16][8] = -0.19638; qual_match_simple_bayesian[16][9] = -0.158755; qual_match_simple_bayesian[16][10] = -0.129846; qual_match_simple_bayesian[16][11] = -0.107464; qual_match_simple_bayesian[16][12] = -0.0900357; qual_match_simple_bayesian[16][13] = -0.0764049; qual_match_simple_bayesian[16][14] = -0.0657085; qual_match_simple_bayesian[16][15] = -0.0572929; qual_match_simple_bayesian[16][16] = -0.0506582; qual_match_simple_bayesian[16][17] = -0.0454193; qual_match_simple_bayesian[16][18] = -0.0412773; qual_match_simple_bayesian[16][19] = -0.0379994; qual_match_simple_bayesian[16][20] = -0.0354033; qual_match_simple_bayesian[16][21] = -0.033346; qual_match_simple_bayesian[16][22] = -0.0317148; qual_match_simple_bayesian[16][23] = -0.0304209; qual_match_simple_bayesian[16][24] = -0.0293944; qual_match_simple_bayesian[16][25] = -0.0285798; qual_match_simple_bayesian[16][26] = -0.0279331; qual_match_simple_bayesian[16][27] = -0.0274198; qual_match_simple_bayesian[16][28] = -0.0270122; qual_match_simple_bayesian[16][29] = -0.0266886; qual_match_simple_bayesian[16][30] = -0.0264316; qual_match_simple_bayesian[16][31] = -0.0262275; qual_match_simple_bayesian[16][32] = -0.0260655; qual_match_simple_bayesian[16][33] = -0.0259367; qual_match_simple_bayesian[16][34] = -0.0258345; qual_match_simple_bayesian[16][35] = -0.0257533; qual_match_simple_bayesian[16][36] = -0.0256888; qual_match_simple_bayesian[16][37] = -0.0256376; qual_match_simple_bayesian[16][38] = -0.0255969; qual_match_simple_bayesian[16][39] = -0.0255645; qual_match_simple_bayesian[16][40] = -0.0255389; qual_match_simple_bayesian[16][41] = -0.0255185; qual_match_simple_bayesian[16][42] = -0.0255023; qual_match_simple_bayesian[16][43] = -0.0254894; qual_match_simple_bayesian[16][44] = -0.0254792; qual_match_simple_bayesian[16][45] = -0.0254711; qual_match_simple_bayesian[16][46] = -0.0254646;
        qual_match_simple_bayesian[17][0] = -5.01301; qual_match_simple_bayesian[17][1] = -1.57576; qual_match_simple_bayesian[17][2] = -1.00546; qual_match_simple_bayesian[17][3] = -0.708883; qual_match_simple_bayesian[17][4] = -0.523352; qual_match_simple_bayesian[17][5] = -0.397151; qual_match_simple_bayesian[17][6] = -0.307149; qual_match_simple_bayesian[17][7] = -0.241016; qual_match_simple_bayesian[17][8] = -0.191434; qual_match_simple_bayesian[17][9] = -0.153729; qual_match_simple_bayesian[17][10] = -0.124761; qual_match_simple_bayesian[17][11] = -0.102334; qual_match_simple_bayesian[17][12] = -0.0848716; qual_match_simple_bayesian[17][13] = -0.0712146; qual_match_simple_bayesian[17][14] = -0.0604979; qual_match_simple_bayesian[17][15] = -0.0520664; qual_match_simple_bayesian[17][16] = -0.0454193; qual_match_simple_bayesian[17][17] = -0.0401706; qual_match_simple_bayesian[17][18] = -0.036021; qual_match_simple_bayesian[17][19] = -0.032737; qual_match_simple_bayesian[17][20] = -0.0301362; qual_match_simple_bayesian[17][21] = -0.028075; qual_match_simple_bayesian[17][22] = -0.0264408; qual_match_simple_bayesian[17][23] = -0.0251447; qual_match_simple_bayesian[17][24] = -0.0241163; qual_match_simple_bayesian[17][25] = -0.0233001; qual_match_simple_bayesian[17][26] = -0.0226523; qual_match_simple_bayesian[17][27] = -0.0221381; qual_match_simple_bayesian[17][28] = -0.0217297; qual_match_simple_bayesian[17][29] = -0.0214055; qual_match_simple_bayesian[17][30] = -0.0211481; qual_match_simple_bayesian[17][31] = -0.0209436; qual_match_simple_bayesian[17][32] = -0.0207812; qual_match_simple_bayesian[17][33] = -0.0206523; qual_match_simple_bayesian[17][34] = -0.0205498; qual_match_simple_bayesian[17][35] = -0.0204685; qual_match_simple_bayesian[17][36] = -0.0204039; qual_match_simple_bayesian[17][37] = -0.0203526; qual_match_simple_bayesian[17][38] = -0.0203118; qual_match_simple_bayesian[17][39] = -0.0202794; qual_match_simple_bayesian[17][40] = -0.0202537; qual_match_simple_bayesian[17][41] = -0.0202333; qual_match_simple_bayesian[17][42] = -0.020217; qual_match_simple_bayesian[17][43] = -0.0202041; qual_match_simple_bayesian[17][44] = -0.0201939; qual_match_simple_bayesian[17][45] = -0.0201858; qual_match_simple_bayesian[17][46] = -0.0201793;
        qual_match_simple_bayesian[18][0] = -5.24327; qual_match_simple_bayesian[18][1] = -1.57693; qual_match_simple_bayesian[18][2] = -1.00368; qual_match_simple_bayesian[18][3] = -0.706121; qual_match_simple_bayesian[18][4] = -0.520107; qual_match_simple_bayesian[18][5] = -0.393627; qual_match_simple_bayesian[18][6] = -0.303445; qual_match_simple_bayesian[18][7] = -0.23719; qual_match_simple_bayesian[18][8] = -0.187522; qual_match_simple_bayesian[18][9] = -0.149755; qual_match_simple_bayesian[18][10] = -0.12074; qual_match_simple_bayesian[18][11] = -0.0982781; qual_match_simple_bayesian[18][12] = -0.0807886; qual_match_simple_bayesian[18][13] = -0.0671109; qual_match_simple_bayesian[18][14] = -0.0563782; qual_match_simple_bayesian[18][15] = -0.0479342; qual_match_simple_bayesian[18][16] = -0.0412773; qual_match_simple_bayesian[18][17] = -0.036021; qual_match_simple_bayesian[18][18] = -0.0318653; qual_match_simple_bayesian[18][19] = -0.0285766; qual_match_simple_bayesian[18][20] = -0.025972; qual_match_simple_bayesian[18][21] = -0.0239079; qual_match_simple_bayesian[18][22] = -0.0222713; qual_match_simple_bayesian[18][23] = -0.0209733; qual_match_simple_bayesian[18][24] = -0.0199434; qual_match_simple_bayesian[18][25] = -0.0191261; qual_match_simple_bayesian[18][26] = -0.0184774; qual_match_simple_bayesian[18][27] = -0.0179624; qual_match_simple_bayesian[18][28] = -0.0175535; qual_match_simple_bayesian[18][29] = -0.0172288; qual_match_simple_bayesian[18][30] = -0.016971; qual_match_simple_bayesian[18][31] = -0.0167662; qual_match_simple_bayesian[18][32] = -0.0166036; qual_match_simple_bayesian[18][33] = -0.0164745; qual_match_simple_bayesian[18][34] = -0.0163719; qual_match_simple_bayesian[18][35] = -0.0162904; qual_match_simple_bayesian[18][36] = -0.0162257; qual_match_simple_bayesian[18][37] = -0.0161743; qual_match_simple_bayesian[18][38] = -0.0161335; qual_match_simple_bayesian[18][39] = -0.0161011; qual_match_simple_bayesian[18][40] = -0.0160753; qual_match_simple_bayesian[18][41] = -0.0160549; qual_match_simple_bayesian[18][42] = -0.0160386; qual_match_simple_bayesian[18][43] = -0.0160257; qual_match_simple_bayesian[18][44] = -0.0160155; qual_match_simple_bayesian[18][45] = -0.0160073; qual_match_simple_bayesian[18][46] = -0.0160009;
        qual_match_simple_bayesian[19][0] = -5.47352; qual_match_simple_bayesian[19][1] = -1.57786; qual_match_simple_bayesian[19][2] = -1.00227; qual_match_simple_bayesian[19][3] = -0.703933; qual_match_simple_bayesian[19][4] = -0.517538; qual_match_simple_bayesian[19][5] = -0.390836; qual_match_simple_bayesian[19][6] = -0.300513; qual_match_simple_bayesian[19][7] = -0.234162; qual_match_simple_bayesian[19][8] = -0.184426; qual_match_simple_bayesian[19][9] = -0.146609; qual_match_simple_bayesian[19][10] = -0.117558; qual_match_simple_bayesian[19][11] = -0.0950678; qual_match_simple_bayesian[19][12] = -0.0775572; qual_match_simple_bayesian[19][13] = -0.0638632; qual_match_simple_bayesian[19][14] = -0.0531178; qual_match_simple_bayesian[19][15] = -0.044664; qual_match_simple_bayesian[19][16] = -0.0379994; qual_match_simple_bayesian[19][17] = -0.032737; qual_match_simple_bayesian[19][18] = -0.0285766; qual_match_simple_bayesian[19][19] = -0.0252842; qual_match_simple_bayesian[19][20] = -0.0226766; qual_match_simple_bayesian[19][21] = -0.0206101; qual_match_simple_bayesian[19][22] = -0.0189717; qual_match_simple_bayesian[19][23] = -0.0176722; qual_match_simple_bayesian[19][24] = -0.0166412; qual_match_simple_bayesian[19][25] = -0.015823; qual_match_simple_bayesian[19][26] = -0.0151735; qual_match_simple_bayesian[19][27] = -0.0146579; qual_match_simple_bayesian[19][28] = -0.0142486; qual_match_simple_bayesian[19][29] = -0.0139235; qual_match_simple_bayesian[19][30] = -0.0136654; qual_match_simple_bayesian[19][31] = -0.0134604; qual_match_simple_bayesian[19][32] = -0.0132976; qual_match_simple_bayesian[19][33] = -0.0131684; qual_match_simple_bayesian[19][34] = -0.0130657; qual_match_simple_bayesian[19][35] = -0.0129841; qual_match_simple_bayesian[19][36] = -0.0129193; qual_match_simple_bayesian[19][37] = -0.0128679; qual_match_simple_bayesian[19][38] = -0.012827; qual_match_simple_bayesian[19][39] = -0.0127945; qual_match_simple_bayesian[19][40] = -0.0127688; qual_match_simple_bayesian[19][41] = -0.0127483; qual_match_simple_bayesian[19][42] = -0.012732; qual_match_simple_bayesian[19][43] = -0.0127191; qual_match_simple_bayesian[19][44] = -0.0127088; qual_match_simple_bayesian[19][45] = -0.0127007; qual_match_simple_bayesian[19][46] = -0.0126942;
        qual_match_simple_bayesian[20][0] = -5.70378; qual_match_simple_bayesian[20][1] = -1.5786; qual_match_simple_bayesian[20][2] = -1.00115; qual_match_simple_bayesian[20][3] = -0.702197; qual_match_simple_bayesian[20][4] = -0.515502; qual_match_simple_bayesian[20][5] = -0.388625; qual_match_simple_bayesian[20][6] = -0.29819; qual_match_simple_bayesian[20][7] = -0.231763; qual_match_simple_bayesian[20][8] = -0.181973; qual_match_simple_bayesian[20][9] = -0.144117; qual_match_simple_bayesian[20][10] = -0.115037; qual_match_simple_bayesian[20][11] = -0.0925252; qual_match_simple_bayesian[20][12] = -0.0749978; qual_match_simple_bayesian[20][13] = -0.061291; qual_match_simple_bayesian[20][14] = -0.0505356; qual_match_simple_bayesian[20][15] = -0.042074; qual_match_simple_bayesian[20][16] = -0.0354033; qual_match_simple_bayesian[20][17] = -0.0301362; qual_match_simple_bayesian[20][18] = -0.025972; qual_match_simple_bayesian[20][19] = -0.0226766; qual_match_simple_bayesian[20][20] = -0.0200667; qual_match_simple_bayesian[20][21] = -0.0179984; qual_match_simple_bayesian[20][22] = -0.0163585; qual_match_simple_bayesian[20][23] = -0.0150578; qual_match_simple_bayesian[20][24] = -0.0140259; qual_match_simple_bayesian[20][25] = -0.0132069; qual_match_simple_bayesian[20][26] = -0.0125569; qual_match_simple_bayesian[20][27] = -0.0120409; qual_match_simple_bayesian[20][28] = -0.0116311; qual_match_simple_bayesian[20][29] = -0.0113058; qual_match_simple_bayesian[20][30] = -0.0110475; qual_match_simple_bayesian[20][31] = -0.0108423; qual_match_simple_bayesian[20][32] = -0.0106794; qual_match_simple_bayesian[20][33] = -0.01055; qual_match_simple_bayesian[20][34] = -0.0104472; qual_match_simple_bayesian[20][35] = -0.0103655; qual_match_simple_bayesian[20][36] = -0.0103007; qual_match_simple_bayesian[20][37] = -0.0102492; qual_match_simple_bayesian[20][38] = -0.0102083; qual_match_simple_bayesian[20][39] = -0.0101758; qual_match_simple_bayesian[20][40] = -0.01015; qual_match_simple_bayesian[20][41] = -0.0101295; qual_match_simple_bayesian[20][42] = -0.0101132; qual_match_simple_bayesian[20][43] = -0.0101003; qual_match_simple_bayesian[20][44] = -0.01009; qual_match_simple_bayesian[20][45] = -0.0100819; qual_match_simple_bayesian[20][46] = -0.0100754;
        qual_match_simple_bayesian[21][0] = -5.93404; qual_match_simple_bayesian[21][1] = -1.57919; qual_match_simple_bayesian[21][2] = -1.00027; qual_match_simple_bayesian[21][3] = -0.700821; qual_match_simple_bayesian[21][4] = -0.513887; qual_match_simple_bayesian[21][5] = -0.386872; qual_match_simple_bayesian[21][6] = -0.296348; qual_match_simple_bayesian[21][7] = -0.229861; qual_match_simple_bayesian[21][8] = -0.180029; qual_match_simple_bayesian[21][9] = -0.142143; qual_match_simple_bayesian[21][10] = -0.113039; qual_match_simple_bayesian[21][11] = -0.09051; qual_match_simple_bayesian[21][12] = -0.0729694; qual_match_simple_bayesian[21][13] = -0.0592525; qual_match_simple_bayesian[21][14] = -0.0484892; qual_match_simple_bayesian[21][15] = -0.0400214; qual_match_simple_bayesian[21][16] = -0.033346; qual_match_simple_bayesian[21][17] = -0.028075; qual_match_simple_bayesian[21][18] = -0.0239079; qual_match_simple_bayesian[21][19] = -0.0206101; qual_match_simple_bayesian[21][20] = -0.0179984; qual_match_simple_bayesian[21][21] = -0.0159286; qual_match_simple_bayesian[21][22] = -0.0142876; qual_match_simple_bayesian[21][23] = -0.012986; qual_match_simple_bayesian[21][24] = -0.0119533; qual_match_simple_bayesian[21][25] = -0.0111338; qual_match_simple_bayesian[21][26] = -0.0104833; qual_match_simple_bayesian[21][27] = -0.00996692; qual_match_simple_bayesian[21][28] = -0.00955691; qual_match_simple_bayesian[21][29] = -0.00923135; qual_match_simple_bayesian[21][30] = -0.00897283; qual_match_simple_bayesian[21][31] = -0.00876752; qual_match_simple_bayesian[21][32] = -0.00860447; qual_match_simple_bayesian[21][33] = -0.00847497; qual_match_simple_bayesian[21][34] = -0.00837212; qual_match_simple_bayesian[21][35] = -0.00829043; qual_match_simple_bayesian[21][36] = -0.00822555; qual_match_simple_bayesian[21][37] = -0.00817401; qual_match_simple_bayesian[21][38] = -0.00813308; qual_match_simple_bayesian[21][39] = -0.00810056; qual_match_simple_bayesian[21][40] = -0.00807474; qual_match_simple_bayesian[21][41] = -0.00805422; qual_match_simple_bayesian[21][42] = -0.00803793; qual_match_simple_bayesian[21][43] = -0.00802498; qual_match_simple_bayesian[21][44] = -0.0080147; qual_match_simple_bayesian[21][45] = -0.00800654; qual_match_simple_bayesian[21][46] = -0.00800005;
        qual_match_simple_bayesian[22][0] = -6.1643; qual_match_simple_bayesian[22][1] = -1.57966; qual_match_simple_bayesian[22][2] = -0.99956; qual_match_simple_bayesian[22][3] = -0.69973; qual_match_simple_bayesian[22][4] = -0.512606; qual_match_simple_bayesian[22][5] = -0.385482; qual_match_simple_bayesian[22][6] = -0.294888; qual_match_simple_bayesian[22][7] = -0.228354; qual_match_simple_bayesian[22][8] = -0.178488; qual_match_simple_bayesian[22][9] = -0.140577; qual_match_simple_bayesian[22][10] = -0.111455; qual_match_simple_bayesian[22][11] = -0.0889123; qual_match_simple_bayesian[22][12] = -0.0713612; qual_match_simple_bayesian[22][13] = -0.0576362; qual_match_simple_bayesian[22][14] = -0.0468667; qual_match_simple_bayesian[22][15] = -0.038394; qual_match_simple_bayesian[22][16] = -0.0317148; qual_match_simple_bayesian[22][17] = -0.0264408; qual_match_simple_bayesian[22][18] = -0.0222713; qual_match_simple_bayesian[22][19] = -0.0189717; qual_match_simple_bayesian[22][20] = -0.0163585; qual_match_simple_bayesian[22][21] = -0.0142876; qual_match_simple_bayesian[22][22] = -0.0126457; qual_match_simple_bayesian[22][23] = -0.0113434; qual_match_simple_bayesian[22][24] = -0.0103101; qual_match_simple_bayesian[22][25] = -0.00949014; qual_match_simple_bayesian[22][26] = -0.00883928; qual_match_simple_bayesian[22][27] = -0.00832259; qual_match_simple_bayesian[22][28] = -0.00791235; qual_match_simple_bayesian[22][29] = -0.00758661; qual_match_simple_bayesian[22][30] = -0.00732794; qual_match_simple_bayesian[22][31] = -0.00712252; qual_match_simple_bayesian[22][32] = -0.00695938; qual_match_simple_bayesian[22][33] = -0.00682981; qual_match_simple_bayesian[22][34] = -0.00672691; qual_match_simple_bayesian[22][35] = -0.00664517; qual_match_simple_bayesian[22][36] = -0.00658025; qual_match_simple_bayesian[22][37] = -0.00652869; qual_match_simple_bayesian[22][38] = -0.00648773; qual_match_simple_bayesian[22][39] = -0.0064552; qual_match_simple_bayesian[22][40] = -0.00642936; qual_match_simple_bayesian[22][41] = -0.00640883; qual_match_simple_bayesian[22][42] = -0.00639253; qual_match_simple_bayesian[22][43] = -0.00637958; qual_match_simple_bayesian[22][44] = -0.00636929; qual_match_simple_bayesian[22][45] = -0.00636112; qual_match_simple_bayesian[22][46] = -0.00635463;
        qual_match_simple_bayesian[23][0] = -6.39456; qual_match_simple_bayesian[23][1] = -1.58003; qual_match_simple_bayesian[23][2] = -0.999001; qual_match_simple_bayesian[23][3] = -0.698863; qual_match_simple_bayesian[23][4] = -0.51159; qual_match_simple_bayesian[23][5] = -0.384379; qual_match_simple_bayesian[23][6] = -0.29373; qual_match_simple_bayesian[23][7] = -0.227158; qual_match_simple_bayesian[23][8] = -0.177265; qual_match_simple_bayesian[23][9] = -0.139335; qual_match_simple_bayesian[23][10] = -0.110198; qual_match_simple_bayesian[23][11] = -0.0876449; qual_match_simple_bayesian[23][12] = -0.0700856; qual_match_simple_bayesian[23][13] = -0.0563542; qual_match_simple_bayesian[23][14] = -0.0455797; qual_match_simple_bayesian[23][15] = -0.0371032; qual_match_simple_bayesian[23][16] = -0.0304209; qual_match_simple_bayesian[23][17] = -0.0251447; qual_match_simple_bayesian[23][18] = -0.0209733; qual_match_simple_bayesian[23][19] = -0.0176722; qual_match_simple_bayesian[23][20] = -0.0150578; qual_match_simple_bayesian[23][21] = -0.012986; qual_match_simple_bayesian[23][22] = -0.0113434; qual_match_simple_bayesian[23][23] = -0.0100405; qual_match_simple_bayesian[23][24] = -0.00900678; qual_match_simple_bayesian[23][25] = -0.00818644; qual_match_simple_bayesian[23][26] = -0.00753529; qual_match_simple_bayesian[23][27] = -0.00701837; qual_match_simple_bayesian[23][28] = -0.00660796; qual_match_simple_bayesian[23][29] = -0.00628208; qual_match_simple_bayesian[23][30] = -0.00602329; qual_match_simple_bayesian[23][31] = -0.00581778; qual_match_simple_bayesian[23][32] = -0.00565457; qual_match_simple_bayesian[23][33] = -0.00552494; qual_match_simple_bayesian[23][34] = -0.00542199; qual_match_simple_bayesian[23][35] = -0.00534022; qual_match_simple_bayesian[23][36] = -0.00527527; qual_match_simple_bayesian[23][37] = -0.00522368; qual_match_simple_bayesian[23][38] = -0.00518271; qual_match_simple_bayesian[23][39] = -0.00515016; qual_match_simple_bayesian[23][40] = -0.00512431; qual_match_simple_bayesian[23][41] = -0.00510378; qual_match_simple_bayesian[23][42] = -0.00508747; qual_match_simple_bayesian[23][43] = -0.00507451; qual_match_simple_bayesian[23][44] = -0.00506422; qual_match_simple_bayesian[23][45] = -0.00505604; qual_match_simple_bayesian[23][46] = -0.00504955;
        qual_match_simple_bayesian[24][0] = -6.62482; qual_match_simple_bayesian[24][1] = -1.58033; qual_match_simple_bayesian[24][2] = -0.998557; qual_match_simple_bayesian[24][3] = -0.698176; qual_match_simple_bayesian[24][4] = -0.510784; qual_match_simple_bayesian[24][5] = -0.383503; qual_match_simple_bayesian[24][6] = -0.29281; qual_match_simple_bayesian[24][7] = -0.226208; qual_match_simple_bayesian[24][8] = -0.176295; qual_match_simple_bayesian[24][9] = -0.138349; qual_match_simple_bayesian[24][10] = -0.109202; qual_match_simple_bayesian[24][11] = -0.0866394; qual_match_simple_bayesian[24][12] = -0.0690735; qual_match_simple_bayesian[24][13] = -0.055337; qual_match_simple_bayesian[24][14] = -0.0445586; qual_match_simple_bayesian[24][15] = -0.0360791; qual_match_simple_bayesian[24][16] = -0.0293944; qual_match_simple_bayesian[24][17] = -0.0241163; qual_match_simple_bayesian[24][18] = -0.0199434; qual_match_simple_bayesian[24][19] = -0.0166412; qual_match_simple_bayesian[24][20] = -0.0140259; qual_match_simple_bayesian[24][21] = -0.0119533; qual_match_simple_bayesian[24][22] = -0.0103101; qual_match_simple_bayesian[24][23] = -0.00900678; qual_match_simple_bayesian[24][24] = -0.00797271; qual_match_simple_bayesian[24][25] = -0.00715208; qual_match_simple_bayesian[24][26] = -0.00650071; qual_match_simple_bayesian[24][27] = -0.00598361; qual_match_simple_bayesian[24][28] = -0.00557305; qual_match_simple_bayesian[24][29] = -0.00524706; qual_match_simple_bayesian[24][30] = -0.00498818; qual_match_simple_bayesian[24][31] = -0.0047826; qual_match_simple_bayesian[24][32] = -0.00461933; qual_match_simple_bayesian[24][33] = -0.00448966; qual_match_simple_bayesian[24][34] = -0.00438667; qual_match_simple_bayesian[24][35] = -0.00430487; qual_match_simple_bayesian[24][36] = -0.0042399; qual_match_simple_bayesian[24][37] = -0.0041883; qual_match_simple_bayesian[24][38] = -0.00414731; qual_match_simple_bayesian[24][39] = -0.00411475; qual_match_simple_bayesian[24][40] = -0.00408889; qual_match_simple_bayesian[24][41] = -0.00406835; qual_match_simple_bayesian[24][42] = -0.00405203; qual_match_simple_bayesian[24][43] = -0.00403907; qual_match_simple_bayesian[24][44] = -0.00402878; qual_match_simple_bayesian[24][45] = -0.0040206; qual_match_simple_bayesian[24][46] = -0.0040141;
        qual_match_simple_bayesian[25][0] = -6.85508; qual_match_simple_bayesian[25][1] = -1.58057; qual_match_simple_bayesian[25][2] = -0.998204; qual_match_simple_bayesian[25][3] = -0.69763; qual_match_simple_bayesian[25][4] = -0.510144; qual_match_simple_bayesian[25][5] = -0.382809; qual_match_simple_bayesian[25][6] = -0.292081; qual_match_simple_bayesian[25][7] = -0.225455; qual_match_simple_bayesian[25][8] = -0.175525; qual_match_simple_bayesian[25][9] = -0.137567; qual_match_simple_bayesian[25][10] = -0.10841; qual_match_simple_bayesian[25][11] = -0.0858414; qual_match_simple_bayesian[25][12] = -0.0682703; qual_match_simple_bayesian[25][13] = -0.0545298; qual_match_simple_bayesian[25][14] = -0.0437483; qual_match_simple_bayesian[25][15] = -0.0352663; qual_match_simple_bayesian[25][16] = -0.0285798; qual_match_simple_bayesian[25][17] = -0.0233001; qual_match_simple_bayesian[25][18] = -0.0191261; qual_match_simple_bayesian[25][19] = -0.015823; qual_match_simple_bayesian[25][20] = -0.0132069; qual_match_simple_bayesian[25][21] = -0.0111338; qual_match_simple_bayesian[25][22] = -0.00949014; qual_match_simple_bayesian[25][23] = -0.00818644; qual_match_simple_bayesian[25][24] = -0.00715208; qual_match_simple_bayesian[25][25] = -0.00633122; qual_match_simple_bayesian[25][26] = -0.00567967; qual_match_simple_bayesian[25][27] = -0.00516243; qual_match_simple_bayesian[25][28] = -0.00475176; qual_match_simple_bayesian[25][29] = -0.00442567; qual_match_simple_bayesian[25][30] = -0.00416673; qual_match_simple_bayesian[25][31] = -0.00396109; qual_match_simple_bayesian[25][32] = -0.00379778; qual_match_simple_bayesian[25][33] = -0.00366807; qual_match_simple_bayesian[25][34] = -0.00356505; qual_match_simple_bayesian[25][35] = -0.00348323; qual_match_simple_bayesian[25][36] = -0.00341824; qual_match_simple_bayesian[25][37] = -0.00336662; qual_match_simple_bayesian[25][38] = -0.00332562; qual_match_simple_bayesian[25][39] = -0.00329306; qual_match_simple_bayesian[25][40] = -0.00326719; qual_match_simple_bayesian[25][41] = -0.00324664; qual_match_simple_bayesian[25][42] = -0.00323032; qual_match_simple_bayesian[25][43] = -0.00321736; qual_match_simple_bayesian[25][44] = -0.00320706; qual_match_simple_bayesian[25][45] = -0.00319888; qual_match_simple_bayesian[25][46] = -0.00319238;
        qual_match_simple_bayesian[26][0] = -7.08533; qual_match_simple_bayesian[26][1] = -1.58075; qual_match_simple_bayesian[26][2] = -0.997924; qual_match_simple_bayesian[26][3] = -0.697196; qual_match_simple_bayesian[26][4] = -0.509636; qual_match_simple_bayesian[26][5] = -0.382257; qual_match_simple_bayesian[26][6] = -0.291502; qual_match_simple_bayesian[26][7] = -0.224857; qual_match_simple_bayesian[26][8] = -0.174914; qual_match_simple_bayesian[26][9] = -0.136946; qual_match_simple_bayesian[26][10] = -0.107782; qual_match_simple_bayesian[26][11] = -0.0852079; qual_match_simple_bayesian[26][12] = -0.0676327; qual_match_simple_bayesian[26][13] = -0.053889; qual_match_simple_bayesian[26][14] = -0.0431051; qual_match_simple_bayesian[26][15] = -0.0346212; qual_match_simple_bayesian[26][16] = -0.0279331; qual_match_simple_bayesian[26][17] = -0.0226523; qual_match_simple_bayesian[26][18] = -0.0184774; qual_match_simple_bayesian[26][19] = -0.0151735; qual_match_simple_bayesian[26][20] = -0.0125569; qual_match_simple_bayesian[26][21] = -0.0104833; qual_match_simple_bayesian[26][22] = -0.00883928; qual_match_simple_bayesian[26][23] = -0.00753529; qual_match_simple_bayesian[26][24] = -0.00650071; qual_match_simple_bayesian[26][25] = -0.00567967; qual_match_simple_bayesian[26][26] = -0.00502798; qual_match_simple_bayesian[26][27] = -0.00451062; qual_match_simple_bayesian[26][28] = -0.00409986; qual_match_simple_bayesian[26][29] = -0.00377371; qual_match_simple_bayesian[26][30] = -0.00351471; qual_match_simple_bayesian[26][31] = -0.00330902; qual_match_simple_bayesian[26][32] = -0.00314567; qual_match_simple_bayesian[26][33] = -0.00301594; qual_match_simple_bayesian[26][34] = -0.0029129; qual_match_simple_bayesian[26][35] = -0.00283106; qual_match_simple_bayesian[26][36] = -0.00276606; qual_match_simple_bayesian[26][37] = -0.00271443; qual_match_simple_bayesian[26][38] = -0.00267342; qual_match_simple_bayesian[26][39] = -0.00264084; qual_match_simple_bayesian[26][40] = -0.00261497; qual_match_simple_bayesian[26][41] = -0.00259442; qual_match_simple_bayesian[26][42] = -0.00257809; qual_match_simple_bayesian[26][43] = -0.00256512; qual_match_simple_bayesian[26][44] = -0.00255482; qual_match_simple_bayesian[26][45] = -0.00254664; qual_match_simple_bayesian[26][46] = -0.00254014;
        qual_match_simple_bayesian[27][0] = -7.31559; qual_match_simple_bayesian[27][1] = -1.5809; qual_match_simple_bayesian[27][2] = -0.997702; qual_match_simple_bayesian[27][3] = -0.696852; qual_match_simple_bayesian[27][4] = -0.509232; qual_match_simple_bayesian[27][5] = -0.38182; qual_match_simple_bayesian[27][6] = -0.291042; qual_match_simple_bayesian[27][7] = -0.224383; qual_match_simple_bayesian[27][8] = -0.174428; qual_match_simple_bayesian[27][9] = -0.136453; qual_match_simple_bayesian[27][10] = -0.107284; qual_match_simple_bayesian[27][11] = -0.0847051; qual_match_simple_bayesian[27][12] = -0.0671265; qual_match_simple_bayesian[27][13] = -0.0533804; qual_match_simple_bayesian[27][14] = -0.0425945; qual_match_simple_bayesian[27][15] = -0.0341091; qual_match_simple_bayesian[27][16] = -0.0274198; qual_match_simple_bayesian[27][17] = -0.0221381; qual_match_simple_bayesian[27][18] = -0.0179624; qual_match_simple_bayesian[27][19] = -0.0146579; qual_match_simple_bayesian[27][20] = -0.0120409; qual_match_simple_bayesian[27][21] = -0.00996692; qual_match_simple_bayesian[27][22] = -0.00832259; qual_match_simple_bayesian[27][23] = -0.00701837; qual_match_simple_bayesian[27][24] = -0.00598361; qual_match_simple_bayesian[27][25] = -0.00516243; qual_match_simple_bayesian[27][26] = -0.00451062; qual_match_simple_bayesian[27][27] = -0.00399318; qual_match_simple_bayesian[27][28] = -0.00358235; qual_match_simple_bayesian[27][29] = -0.00325613; qual_match_simple_bayesian[27][30] = -0.00299709; qual_match_simple_bayesian[27][31] = -0.00279137; qual_match_simple_bayesian[27][32] = -0.00262799; qual_match_simple_bayesian[27][33] = -0.00249823; qual_match_simple_bayesian[27][34] = -0.00239518; qual_match_simple_bayesian[27][35] = -0.00231332; qual_match_simple_bayesian[27][36] = -0.00224831; qual_match_simple_bayesian[27][37] = -0.00219667; qual_match_simple_bayesian[27][38] = -0.00215565; qual_match_simple_bayesian[27][39] = -0.00212307; qual_match_simple_bayesian[27][40] = -0.00209719; qual_match_simple_bayesian[27][41] = -0.00207664; qual_match_simple_bayesian[27][42] = -0.00206031; qual_match_simple_bayesian[27][43] = -0.00204734; qual_match_simple_bayesian[27][44] = -0.00203704; qual_match_simple_bayesian[27][45] = -0.00202886; qual_match_simple_bayesian[27][46] = -0.00202236;
        qual_match_simple_bayesian[28][0] = -7.54585; qual_match_simple_bayesian[28][1] = -1.58102; qual_match_simple_bayesian[28][2] = -0.997525; qual_match_simple_bayesian[28][3] = -0.696579; qual_match_simple_bayesian[28][4] = -0.508912; qual_match_simple_bayesian[28][5] = -0.381472; qual_match_simple_bayesian[28][6] = -0.290677; qual_match_simple_bayesian[28][7] = -0.224006; qual_match_simple_bayesian[28][8] = -0.174043; qual_match_simple_bayesian[28][9] = -0.136062; qual_match_simple_bayesian[28][10] = -0.106888; qual_match_simple_bayesian[28][11] = -0.0843058; qual_match_simple_bayesian[28][12] = -0.0667247; qual_match_simple_bayesian[28][13] = -0.0529765; qual_match_simple_bayesian[28][14] = -0.0421891; qual_match_simple_bayesian[28][15] = -0.0337024; qual_match_simple_bayesian[28][16] = -0.0270122; qual_match_simple_bayesian[28][17] = -0.0217297; qual_match_simple_bayesian[28][18] = -0.0175535; qual_match_simple_bayesian[28][19] = -0.0142486; qual_match_simple_bayesian[28][20] = -0.0116311; qual_match_simple_bayesian[28][21] = -0.00955691; qual_match_simple_bayesian[28][22] = -0.00791235; qual_match_simple_bayesian[28][23] = -0.00660796; qual_match_simple_bayesian[28][24] = -0.00557305; qual_match_simple_bayesian[28][25] = -0.00475176; qual_match_simple_bayesian[28][26] = -0.00409986; qual_match_simple_bayesian[28][27] = -0.00358235; qual_match_simple_bayesian[28][28] = -0.00317146; qual_match_simple_bayesian[28][29] = -0.0028452; qual_match_simple_bayesian[28][30] = -0.00258612; qual_match_simple_bayesian[28][31] = -0.00238037; qual_match_simple_bayesian[28][32] = -0.00221697; qual_match_simple_bayesian[28][33] = -0.0020872; qual_match_simple_bayesian[28][34] = -0.00198413; qual_match_simple_bayesian[28][35] = -0.00190226; qual_match_simple_bayesian[28][36] = -0.00183724; qual_match_simple_bayesian[28][37] = -0.00178559; qual_match_simple_bayesian[28][38] = -0.00174457; qual_match_simple_bayesian[28][39] = -0.00171198; qual_match_simple_bayesian[28][40] = -0.0016861; qual_match_simple_bayesian[28][41] = -0.00166554; qual_match_simple_bayesian[28][42] = -0.00164921; qual_match_simple_bayesian[28][43] = -0.00163624; qual_match_simple_bayesian[28][44] = -0.00162594; qual_match_simple_bayesian[28][45] = -0.00161776; qual_match_simple_bayesian[28][46] = -0.00161126;
        qual_match_simple_bayesian[29][0] = -7.77611; qual_match_simple_bayesian[29][1] = -1.58111; qual_match_simple_bayesian[29][2] = -0.997385; qual_match_simple_bayesian[29][3] = -0.696362; qual_match_simple_bayesian[29][4] = -0.508658; qual_match_simple_bayesian[29][5] = -0.381196; qual_match_simple_bayesian[29][6] = -0.290387; qual_match_simple_bayesian[29][7] = -0.223707; qual_match_simple_bayesian[29][8] = -0.173737; qual_match_simple_bayesian[29][9] = -0.135751; qual_match_simple_bayesian[29][10] = -0.106574; qual_match_simple_bayesian[29][11] = -0.0839888; qual_match_simple_bayesian[29][12] = -0.0664056; qual_match_simple_bayesian[29][13] = -0.0526558; qual_match_simple_bayesian[29][14] = -0.0418671; qual_match_simple_bayesian[29][15] = -0.0333796; qual_match_simple_bayesian[29][16] = -0.0266886; qual_match_simple_bayesian[29][17] = -0.0214055; qual_match_simple_bayesian[29][18] = -0.0172288; qual_match_simple_bayesian[29][19] = -0.0139235; qual_match_simple_bayesian[29][20] = -0.0113058; qual_match_simple_bayesian[29][21] = -0.00923135; qual_match_simple_bayesian[29][22] = -0.00758661; qual_match_simple_bayesian[29][23] = -0.00628208; qual_match_simple_bayesian[29][24] = -0.00524706; qual_match_simple_bayesian[29][25] = -0.00442567; qual_match_simple_bayesian[29][26] = -0.00377371; qual_match_simple_bayesian[29][27] = -0.00325613; qual_match_simple_bayesian[29][28] = -0.0028452; qual_match_simple_bayesian[29][29] = -0.00251891; qual_match_simple_bayesian[29][30] = -0.0022598; qual_match_simple_bayesian[29][31] = -0.00205403; qual_match_simple_bayesian[29][32] = -0.00189061; qual_match_simple_bayesian[29][33] = -0.00176082; qual_match_simple_bayesian[29][34] = -0.00165774; qual_match_simple_bayesian[29][35] = -0.00157586; qual_match_simple_bayesian[29][36] = -0.00151083; qual_match_simple_bayesian[29][37] = -0.00145918; qual_match_simple_bayesian[29][38] = -0.00141815; qual_match_simple_bayesian[29][39] = -0.00138557; qual_match_simple_bayesian[29][40] = -0.00135968; qual_match_simple_bayesian[29][41] = -0.00133912; qual_match_simple_bayesian[29][42] = -0.00132279; qual_match_simple_bayesian[29][43] = -0.00130982; qual_match_simple_bayesian[29][44] = -0.00129951; qual_match_simple_bayesian[29][45] = -0.00129133; qual_match_simple_bayesian[29][46] = -0.00128483;
        qual_match_simple_bayesian[30][0] = -8.00637; qual_match_simple_bayesian[30][1] = -1.58119; qual_match_simple_bayesian[30][2] = -0.997273; qual_match_simple_bayesian[30][3] = -0.69619; qual_match_simple_bayesian[30][4] = -0.508456; qual_match_simple_bayesian[30][5] = -0.380977; qual_match_simple_bayesian[30][6] = -0.290157; qual_match_simple_bayesian[30][7] = -0.223469; qual_match_simple_bayesian[30][8] = -0.173494; qual_match_simple_bayesian[30][9] = -0.135504; qual_match_simple_bayesian[30][10] = -0.106324; qual_match_simple_bayesian[30][11] = -0.083737; qual_match_simple_bayesian[30][12] = -0.0661522; qual_match_simple_bayesian[30][13] = -0.0524012; qual_match_simple_bayesian[30][14] = -0.0416115; qual_match_simple_bayesian[30][15] = -0.0331232; qual_match_simple_bayesian[30][16] = -0.0264316; qual_match_simple_bayesian[30][17] = -0.0211481; qual_match_simple_bayesian[30][18] = -0.016971; qual_match_simple_bayesian[30][19] = -0.0136654; qual_match_simple_bayesian[30][20] = -0.0110475; qual_match_simple_bayesian[30][21] = -0.00897283; qual_match_simple_bayesian[30][22] = -0.00732794; qual_match_simple_bayesian[30][23] = -0.00602329; qual_match_simple_bayesian[30][24] = -0.00498818; qual_match_simple_bayesian[30][25] = -0.00416673; qual_match_simple_bayesian[30][26] = -0.00351471; qual_match_simple_bayesian[30][27] = -0.00299709; qual_match_simple_bayesian[30][28] = -0.00258612; qual_match_simple_bayesian[30][29] = -0.0022598; qual_match_simple_bayesian[30][30] = -0.00200067; qual_match_simple_bayesian[30][31] = -0.00179488; qual_match_simple_bayesian[30][32] = -0.00163145; qual_match_simple_bayesian[30][33] = -0.00150165; qual_match_simple_bayesian[30][34] = -0.00139855; qual_match_simple_bayesian[30][35] = -0.00131667; qual_match_simple_bayesian[30][36] = -0.00125164; qual_match_simple_bayesian[30][37] = -0.00119998; qual_match_simple_bayesian[30][38] = -0.00115895; qual_match_simple_bayesian[30][39] = -0.00112636; qual_match_simple_bayesian[30][40] = -0.00110047; qual_match_simple_bayesian[30][41] = -0.00107991; qual_match_simple_bayesian[30][42] = -0.00106358; qual_match_simple_bayesian[30][43] = -0.0010506; qual_match_simple_bayesian[30][44] = -0.0010403; qual_match_simple_bayesian[30][45] = -0.00103211; qual_match_simple_bayesian[30][46] = -0.00102561;
        qual_match_simple_bayesian[31][0] = -8.23663; qual_match_simple_bayesian[31][1] = -1.58125; qual_match_simple_bayesian[31][2] = -0.997185; qual_match_simple_bayesian[31][3] = -0.696053; qual_match_simple_bayesian[31][4] = -0.508295; qual_match_simple_bayesian[31][5] = -0.380803; qual_match_simple_bayesian[31][6] = -0.289974; qual_match_simple_bayesian[31][7] = -0.22328; qual_match_simple_bayesian[31][8] = -0.173301; qual_match_simple_bayesian[31][9] = -0.135308; qual_match_simple_bayesian[31][10] = -0.106126; qual_match_simple_bayesian[31][11] = -0.0835371; qual_match_simple_bayesian[31][12] = -0.065951; qual_match_simple_bayesian[31][13] = -0.0521989; qual_match_simple_bayesian[31][14] = -0.0414085; qual_match_simple_bayesian[31][15] = -0.0329196; qual_match_simple_bayesian[31][16] = -0.0262275; qual_match_simple_bayesian[31][17] = -0.0209436; qual_match_simple_bayesian[31][18] = -0.0167662; qual_match_simple_bayesian[31][19] = -0.0134604; qual_match_simple_bayesian[31][20] = -0.0108423; qual_match_simple_bayesian[31][21] = -0.00876752; qual_match_simple_bayesian[31][22] = -0.00712252; qual_match_simple_bayesian[31][23] = -0.00581778; qual_match_simple_bayesian[31][24] = -0.0047826; qual_match_simple_bayesian[31][25] = -0.00396109; qual_match_simple_bayesian[31][26] = -0.00330902; qual_match_simple_bayesian[31][27] = -0.00279137; qual_match_simple_bayesian[31][28] = -0.00238037; qual_match_simple_bayesian[31][29] = -0.00205403; qual_match_simple_bayesian[31][30] = -0.00179488; qual_match_simple_bayesian[31][31] = -0.00158908; qual_match_simple_bayesian[31][32] = -0.00142563; qual_match_simple_bayesian[31][33] = -0.00129582; qual_match_simple_bayesian[31][34] = -0.00119272; qual_match_simple_bayesian[31][35] = -0.00111084; qual_match_simple_bayesian[31][36] = -0.0010458; qual_match_simple_bayesian[31][37] = -0.000994137; qual_match_simple_bayesian[31][38] = -0.000953104; qual_match_simple_bayesian[31][39] = -0.000920511; qual_match_simple_bayesian[31][40] = -0.000894622; qual_match_simple_bayesian[31][41] = -0.000874059; qual_match_simple_bayesian[31][42] = -0.000857725; qual_match_simple_bayesian[31][43] = -0.000844751; qual_match_simple_bayesian[31][44] = -0.000834445; qual_match_simple_bayesian[31][45] = -0.000826259; qual_match_simple_bayesian[31][46] = -0.000819756;
        qual_match_simple_bayesian[32][0] = -8.46688; qual_match_simple_bayesian[32][1] = -1.58129; qual_match_simple_bayesian[32][2] = -0.997114; qual_match_simple_bayesian[32][3] = -0.695944; qual_match_simple_bayesian[32][4] = -0.508168; qual_match_simple_bayesian[32][5] = -0.380664; qual_match_simple_bayesian[32][6] = -0.289829; qual_match_simple_bayesian[32][7] = -0.22313; qual_match_simple_bayesian[32][8] = -0.173148; qual_match_simple_bayesian[32][9] = -0.135153; qual_match_simple_bayesian[32][10] = -0.105968; qual_match_simple_bayesian[32][11] = -0.0833783; qual_match_simple_bayesian[32][12] = -0.0657912; qual_match_simple_bayesian[32][13] = -0.0520383; qual_match_simple_bayesian[32][14] = -0.0412473; qual_match_simple_bayesian[32][15] = -0.0327579; qual_match_simple_bayesian[32][16] = -0.0260655; qual_match_simple_bayesian[32][17] = -0.0207812; qual_match_simple_bayesian[32][18] = -0.0166036; qual_match_simple_bayesian[32][19] = -0.0132976; qual_match_simple_bayesian[32][20] = -0.0106794; qual_match_simple_bayesian[32][21] = -0.00860447; qual_match_simple_bayesian[32][22] = -0.00695938; qual_match_simple_bayesian[32][23] = -0.00565457; qual_match_simple_bayesian[32][24] = -0.00461933; qual_match_simple_bayesian[32][25] = -0.00379778; qual_match_simple_bayesian[32][26] = -0.00314567; qual_match_simple_bayesian[32][27] = -0.00262799; qual_match_simple_bayesian[32][28] = -0.00221697; qual_match_simple_bayesian[32][29] = -0.00189061; qual_match_simple_bayesian[32][30] = -0.00163145; qual_match_simple_bayesian[32][31] = -0.00142563; qual_match_simple_bayesian[32][32] = -0.00126218; qual_match_simple_bayesian[32][33] = -0.00113236; qual_match_simple_bayesian[32][34] = -0.00102926; qual_match_simple_bayesian[32][35] = -0.000947368; qual_match_simple_bayesian[32][36] = -0.000882324; qual_match_simple_bayesian[32][37] = -0.000830661; qual_match_simple_bayesian[32][38] = -0.000789625; qual_match_simple_bayesian[32][39] = -0.00075703; qual_match_simple_bayesian[32][40] = -0.00073114; qual_match_simple_bayesian[32][41] = -0.000710576; qual_match_simple_bayesian[32][42] = -0.000694241; qual_match_simple_bayesian[32][43] = -0.000681266; qual_match_simple_bayesian[32][44] = -0.00067096; qual_match_simple_bayesian[32][45] = -0.000662773; qual_match_simple_bayesian[32][46] = -0.00065627;
        qual_match_simple_bayesian[33][0] = -8.69714; qual_match_simple_bayesian[33][1] = -1.58133; qual_match_simple_bayesian[33][2] = -0.997059; qual_match_simple_bayesian[33][3] = -0.695858; qual_match_simple_bayesian[33][4] = -0.508067; qual_match_simple_bayesian[33][5] = -0.380554; qual_match_simple_bayesian[33][6] = -0.289713; qual_match_simple_bayesian[33][7] = -0.223011; qual_match_simple_bayesian[33][8] = -0.173026; qual_match_simple_bayesian[33][9] = -0.135029; qual_match_simple_bayesian[33][10] = -0.105843; qual_match_simple_bayesian[33][11] = -0.0832522; qual_match_simple_bayesian[33][12] = -0.0656642; qual_match_simple_bayesian[33][13] = -0.0519108; qual_match_simple_bayesian[33][14] = -0.0411192; qual_match_simple_bayesian[33][15] = -0.0326294; qual_match_simple_bayesian[33][16] = -0.0259367; qual_match_simple_bayesian[33][17] = -0.0206523; qual_match_simple_bayesian[33][18] = -0.0164745; qual_match_simple_bayesian[33][19] = -0.0131684; qual_match_simple_bayesian[33][20] = -0.01055; qual_match_simple_bayesian[33][21] = -0.00847497; qual_match_simple_bayesian[33][22] = -0.00682981; qual_match_simple_bayesian[33][23] = -0.00552494; qual_match_simple_bayesian[33][24] = -0.00448966; qual_match_simple_bayesian[33][25] = -0.00366807; qual_match_simple_bayesian[33][26] = -0.00301594; qual_match_simple_bayesian[33][27] = -0.00249823; qual_match_simple_bayesian[33][28] = -0.0020872; qual_match_simple_bayesian[33][29] = -0.00176082; qual_match_simple_bayesian[33][30] = -0.00150165; qual_match_simple_bayesian[33][31] = -0.00129582; qual_match_simple_bayesian[33][32] = -0.00113236; qual_match_simple_bayesian[33][33] = -0.00100254; qual_match_simple_bayesian[33][34] = -0.000899433; qual_match_simple_bayesian[33][35] = -0.000817538; qual_match_simple_bayesian[33][36] = -0.000752491; qual_match_simple_bayesian[33][37] = -0.000700826; qual_match_simple_bayesian[33][38] = -0.000659788; qual_match_simple_bayesian[33][39] = -0.000627192; qual_match_simple_bayesian[33][40] = -0.000601301; qual_match_simple_bayesian[33][41] = -0.000580736; qual_match_simple_bayesian[33][42] = -0.0005644; qual_match_simple_bayesian[33][43] = -0.000551424; qual_match_simple_bayesian[33][44] = -0.000541118; qual_match_simple_bayesian[33][45] = -0.000532931; qual_match_simple_bayesian[33][46] = -0.000526428;
        qual_match_simple_bayesian[34][0] = -8.9274; qual_match_simple_bayesian[34][1] = -1.58136; qual_match_simple_bayesian[34][2] = -0.997014; qual_match_simple_bayesian[34][3] = -0.695789; qual_match_simple_bayesian[34][4] = -0.507986; qual_match_simple_bayesian[34][5] = -0.380467; qual_match_simple_bayesian[34][6] = -0.289622; qual_match_simple_bayesian[34][7] = -0.222917; qual_match_simple_bayesian[34][8] = -0.17293; qual_match_simple_bayesian[34][9] = -0.134931; qual_match_simple_bayesian[34][10] = -0.105744; qual_match_simple_bayesian[34][11] = -0.083152; qual_match_simple_bayesian[34][12] = -0.0655634; qual_match_simple_bayesian[34][13] = -0.0518095; qual_match_simple_bayesian[34][14] = -0.0410175; qual_match_simple_bayesian[34][15] = -0.0325274; qual_match_simple_bayesian[34][16] = -0.0258345; qual_match_simple_bayesian[34][17] = -0.0205498; qual_match_simple_bayesian[34][18] = -0.0163719; qual_match_simple_bayesian[34][19] = -0.0130657; qual_match_simple_bayesian[34][20] = -0.0104472; qual_match_simple_bayesian[34][21] = -0.00837212; qual_match_simple_bayesian[34][22] = -0.00672691; qual_match_simple_bayesian[34][23] = -0.00542199; qual_match_simple_bayesian[34][24] = -0.00438667; qual_match_simple_bayesian[34][25] = -0.00356505; qual_match_simple_bayesian[34][26] = -0.0029129; qual_match_simple_bayesian[34][27] = -0.00239518; qual_match_simple_bayesian[34][28] = -0.00198413; qual_match_simple_bayesian[34][29] = -0.00165774; qual_match_simple_bayesian[34][30] = -0.00139855; qual_match_simple_bayesian[34][31] = -0.00119272; qual_match_simple_bayesian[34][32] = -0.00102926; qual_match_simple_bayesian[34][33] = -0.000899433; qual_match_simple_bayesian[34][34] = -0.00079632; qual_match_simple_bayesian[34][35] = -0.000714422; qual_match_simple_bayesian[34][36] = -0.000649373; qual_match_simple_bayesian[34][37] = -0.000597706; qual_match_simple_bayesian[34][38] = -0.000556667; qual_match_simple_bayesian[34][39] = -0.00052407; qual_match_simple_bayesian[34][40] = -0.000498178; qual_match_simple_bayesian[34][41] = -0.000477612; qual_match_simple_bayesian[34][42] = -0.000461276; qual_match_simple_bayesian[34][43] = -0.0004483; qual_match_simple_bayesian[34][44] = -0.000437993; qual_match_simple_bayesian[34][45] = -0.000429806; qual_match_simple_bayesian[34][46] = -0.000423302;
        qual_match_simple_bayesian[35][0] = -9.15766; qual_match_simple_bayesian[35][1] = -1.58138; qual_match_simple_bayesian[35][2] = -0.996979; qual_match_simple_bayesian[35][3] = -0.695735; qual_match_simple_bayesian[35][4] = -0.507922; qual_match_simple_bayesian[35][5] = -0.380398; qual_match_simple_bayesian[35][6] = -0.289549; qual_match_simple_bayesian[35][7] = -0.222842; qual_match_simple_bayesian[35][8] = -0.172853; qual_match_simple_bayesian[35][9] = -0.134853; qual_match_simple_bayesian[35][10] = -0.105665; qual_match_simple_bayesian[35][11] = -0.0830725; qual_match_simple_bayesian[35][12] = -0.0654833; qual_match_simple_bayesian[35][13] = -0.051729; qual_match_simple_bayesian[35][14] = -0.0409368; qual_match_simple_bayesian[35][15] = -0.0324464; qual_match_simple_bayesian[35][16] = -0.0257533; qual_match_simple_bayesian[35][17] = -0.0204685; qual_match_simple_bayesian[35][18] = -0.0162904; qual_match_simple_bayesian[35][19] = -0.0129841; qual_match_simple_bayesian[35][20] = -0.0103655; qual_match_simple_bayesian[35][21] = -0.00829043; qual_match_simple_bayesian[35][22] = -0.00664517; qual_match_simple_bayesian[35][23] = -0.00534022; qual_match_simple_bayesian[35][24] = -0.00430487; qual_match_simple_bayesian[35][25] = -0.00348323; qual_match_simple_bayesian[35][26] = -0.00283106; qual_match_simple_bayesian[35][27] = -0.00231332; qual_match_simple_bayesian[35][28] = -0.00190226; qual_match_simple_bayesian[35][29] = -0.00157586; qual_match_simple_bayesian[35][30] = -0.00131667; qual_match_simple_bayesian[35][31] = -0.00111084; qual_match_simple_bayesian[35][32] = -0.000947368; qual_match_simple_bayesian[35][33] = -0.000817538; qual_match_simple_bayesian[35][34] = -0.000714422; qual_match_simple_bayesian[35][35] = -0.000632522; qual_match_simple_bayesian[35][36] = -0.000567471; qual_match_simple_bayesian[35][37] = -0.000515803; qual_match_simple_bayesian[35][38] = -0.000474763; qual_match_simple_bayesian[35][39] = -0.000442165; qual_match_simple_bayesian[35][40] = -0.000416272; qual_match_simple_bayesian[35][41] = -0.000395705; qual_match_simple_bayesian[35][42] = -0.000379369; qual_match_simple_bayesian[35][43] = -0.000366392; qual_match_simple_bayesian[35][44] = -0.000356085; qual_match_simple_bayesian[35][45] = -0.000347898; qual_match_simple_bayesian[35][46] = -0.000341394;
        qual_match_simple_bayesian[36][0] = -9.38792; qual_match_simple_bayesian[36][1] = -1.5814; qual_match_simple_bayesian[36][2] = -0.996951; qual_match_simple_bayesian[36][3] = -0.695692; qual_match_simple_bayesian[36][4] = -0.507872; qual_match_simple_bayesian[36][5] = -0.380343; qual_match_simple_bayesian[36][6] = -0.289491; qual_match_simple_bayesian[36][7] = -0.222782; qual_match_simple_bayesian[36][8] = -0.172792; qual_match_simple_bayesian[36][9] = -0.134791; qual_match_simple_bayesian[36][10] = -0.105602; qual_match_simple_bayesian[36][11] = -0.0830093; qual_match_simple_bayesian[36][12] = -0.0654198; qual_match_simple_bayesian[36][13] = -0.0516651; qual_match_simple_bayesian[36][14] = -0.0408726; qual_match_simple_bayesian[36][15] = -0.0323821; qual_match_simple_bayesian[36][16] = -0.0256888; qual_match_simple_bayesian[36][17] = -0.0204039; qual_match_simple_bayesian[36][18] = -0.0162257; qual_match_simple_bayesian[36][19] = -0.0129193; qual_match_simple_bayesian[36][20] = -0.0103007; qual_match_simple_bayesian[36][21] = -0.00822555; qual_match_simple_bayesian[36][22] = -0.00658025; qual_match_simple_bayesian[36][23] = -0.00527527; qual_match_simple_bayesian[36][24] = -0.0042399; qual_match_simple_bayesian[36][25] = -0.00341824; qual_match_simple_bayesian[36][26] = -0.00276606; qual_match_simple_bayesian[36][27] = -0.00224831; qual_match_simple_bayesian[36][28] = -0.00183724; qual_match_simple_bayesian[36][29] = -0.00151083; qual_match_simple_bayesian[36][30] = -0.00125164; qual_match_simple_bayesian[36][31] = -0.0010458; qual_match_simple_bayesian[36][32] = -0.000882324; qual_match_simple_bayesian[36][33] = -0.000752491; qual_match_simple_bayesian[36][34] = -0.000649373; qual_match_simple_bayesian[36][35] = -0.000567471; qual_match_simple_bayesian[36][36] = -0.000502419; qual_match_simple_bayesian[36][37] = -0.00045075; qual_match_simple_bayesian[36][38] = -0.000409709; qual_match_simple_bayesian[36][39] = -0.00037711; qual_match_simple_bayesian[36][40] = -0.000351217; qual_match_simple_bayesian[36][41] = -0.00033065; qual_match_simple_bayesian[36][42] = -0.000314313; qual_match_simple_bayesian[36][43] = -0.000301336; qual_match_simple_bayesian[36][44] = -0.000291028; qual_match_simple_bayesian[36][45] = -0.000282841; qual_match_simple_bayesian[36][46] = -0.000276337;
        qual_match_simple_bayesian[37][0] = -9.61818; qual_match_simple_bayesian[37][1] = -1.58142; qual_match_simple_bayesian[37][2] = -0.996929; qual_match_simple_bayesian[37][3] = -0.695657; qual_match_simple_bayesian[37][4] = -0.507831; qual_match_simple_bayesian[37][5] = -0.380299; qual_match_simple_bayesian[37][6] = -0.289445; qual_match_simple_bayesian[37][7] = -0.222734; qual_match_simple_bayesian[37][8] = -0.172744; qual_match_simple_bayesian[37][9] = -0.134742; qual_match_simple_bayesian[37][10] = -0.105553; qual_match_simple_bayesian[37][11] = -0.0829591; qual_match_simple_bayesian[37][12] = -0.0653692; qual_match_simple_bayesian[37][13] = -0.0516143; qual_match_simple_bayesian[37][14] = -0.0408216; qual_match_simple_bayesian[37][15] = -0.0323309; qual_match_simple_bayesian[37][16] = -0.0256376; qual_match_simple_bayesian[37][17] = -0.0203526; qual_match_simple_bayesian[37][18] = -0.0161743; qual_match_simple_bayesian[37][19] = -0.0128679; qual_match_simple_bayesian[37][20] = -0.0102492; qual_match_simple_bayesian[37][21] = -0.00817401; qual_match_simple_bayesian[37][22] = -0.00652869; qual_match_simple_bayesian[37][23] = -0.00522368; qual_match_simple_bayesian[37][24] = -0.0041883; qual_match_simple_bayesian[37][25] = -0.00336662; qual_match_simple_bayesian[37][26] = -0.00271443; qual_match_simple_bayesian[37][27] = -0.00219667; qual_match_simple_bayesian[37][28] = -0.00178559; qual_match_simple_bayesian[37][29] = -0.00145918; qual_match_simple_bayesian[37][30] = -0.00119998; qual_match_simple_bayesian[37][31] = -0.000994137; qual_match_simple_bayesian[37][32] = -0.000830661; qual_match_simple_bayesian[37][33] = -0.000700826; qual_match_simple_bayesian[37][34] = -0.000597706; qual_match_simple_bayesian[37][35] = -0.000515803; qual_match_simple_bayesian[37][36] = -0.00045075; qual_match_simple_bayesian[37][37] = -0.000399079; qual_match_simple_bayesian[37][38] = -0.000358037; qual_match_simple_bayesian[37][39] = -0.000325438; qual_match_simple_bayesian[37][40] = -0.000299544; qual_match_simple_bayesian[37][41] = -0.000278977; qual_match_simple_bayesian[37][42] = -0.00026264; qual_match_simple_bayesian[37][43] = -0.000249663; qual_match_simple_bayesian[37][44] = -0.000239355; qual_match_simple_bayesian[37][45] = -0.000231167; qual_match_simple_bayesian[37][46] = -0.000224664;
        qual_match_simple_bayesian[38][0] = -9.84844; qual_match_simple_bayesian[38][1] = -1.58143; qual_match_simple_bayesian[38][2] = -0.996911; qual_match_simple_bayesian[38][3] = -0.69563; qual_match_simple_bayesian[38][4] = -0.507799; qual_match_simple_bayesian[38][5] = -0.380264; qual_match_simple_bayesian[38][6] = -0.289409; qual_match_simple_bayesian[38][7] = -0.222697; qual_match_simple_bayesian[38][8] = -0.172705; qual_match_simple_bayesian[38][9] = -0.134703; qual_match_simple_bayesian[38][10] = -0.105513; qual_match_simple_bayesian[38][11] = -0.0829192; qual_match_simple_bayesian[38][12] = -0.0653291; qual_match_simple_bayesian[38][13] = -0.051574; qual_match_simple_bayesian[38][14] = -0.0407812; qual_match_simple_bayesian[38][15] = -0.0322904; qual_match_simple_bayesian[38][16] = -0.0255969; qual_match_simple_bayesian[38][17] = -0.0203118; qual_match_simple_bayesian[38][18] = -0.0161335; qual_match_simple_bayesian[38][19] = -0.012827; qual_match_simple_bayesian[38][20] = -0.0102083; qual_match_simple_bayesian[38][21] = -0.00813308; qual_match_simple_bayesian[38][22] = -0.00648773; qual_match_simple_bayesian[38][23] = -0.00518271; qual_match_simple_bayesian[38][24] = -0.00414731; qual_match_simple_bayesian[38][25] = -0.00332562; qual_match_simple_bayesian[38][26] = -0.00267342; qual_match_simple_bayesian[38][27] = -0.00215565; qual_match_simple_bayesian[38][28] = -0.00174457; qual_match_simple_bayesian[38][29] = -0.00141815; qual_match_simple_bayesian[38][30] = -0.00115895; qual_match_simple_bayesian[38][31] = -0.000953104; qual_match_simple_bayesian[38][32] = -0.000789625; qual_match_simple_bayesian[38][33] = -0.000659788; qual_match_simple_bayesian[38][34] = -0.000556667; qual_match_simple_bayesian[38][35] = -0.000474763; qual_match_simple_bayesian[38][36] = -0.000409709; qual_match_simple_bayesian[38][37] = -0.000358037; qual_match_simple_bayesian[38][38] = -0.000316995; qual_match_simple_bayesian[38][39] = -0.000284396; qual_match_simple_bayesian[38][40] = -0.000258502; qual_match_simple_bayesian[38][41] = -0.000237934; qual_match_simple_bayesian[38][42] = -0.000221596; qual_match_simple_bayesian[38][43] = -0.000208619; qual_match_simple_bayesian[38][44] = -0.000198311; qual_match_simple_bayesian[38][45] = -0.000190123; qual_match_simple_bayesian[38][46] = -0.00018362;
        qual_match_simple_bayesian[39][0] = -10.0787; qual_match_simple_bayesian[39][1] = -1.58144; qual_match_simple_bayesian[39][2] = -0.996897; qual_match_simple_bayesian[39][3] = -0.695608; qual_match_simple_bayesian[39][4] = -0.507774; qual_match_simple_bayesian[39][5] = -0.380237; qual_match_simple_bayesian[39][6] = -0.28938; qual_match_simple_bayesian[39][7] = -0.222667; qual_match_simple_bayesian[39][8] = -0.172675; qual_match_simple_bayesian[39][9] = -0.134672; qual_match_simple_bayesian[39][10] = -0.105482; qual_match_simple_bayesian[39][11] = -0.0828876; qual_match_simple_bayesian[39][12] = -0.0652972; qual_match_simple_bayesian[39][13] = -0.051542; qual_match_simple_bayesian[39][14] = -0.040749; qual_match_simple_bayesian[39][15] = -0.0322581; qual_match_simple_bayesian[39][16] = -0.0255645; qual_match_simple_bayesian[39][17] = -0.0202794; qual_match_simple_bayesian[39][18] = -0.0161011; qual_match_simple_bayesian[39][19] = -0.0127945; qual_match_simple_bayesian[39][20] = -0.0101758; qual_match_simple_bayesian[39][21] = -0.00810056; qual_match_simple_bayesian[39][22] = -0.0064552; qual_match_simple_bayesian[39][23] = -0.00515016; qual_match_simple_bayesian[39][24] = -0.00411475; qual_match_simple_bayesian[39][25] = -0.00329306; qual_match_simple_bayesian[39][26] = -0.00264084; qual_match_simple_bayesian[39][27] = -0.00212307; qual_match_simple_bayesian[39][28] = -0.00171198; qual_match_simple_bayesian[39][29] = -0.00138557; qual_match_simple_bayesian[39][30] = -0.00112636; qual_match_simple_bayesian[39][31] = -0.000920511; qual_match_simple_bayesian[39][32] = -0.00075703; qual_match_simple_bayesian[39][33] = -0.000627192; qual_match_simple_bayesian[39][34] = -0.00052407; qual_match_simple_bayesian[39][35] = -0.000442165; qual_match_simple_bayesian[39][36] = -0.00037711; qual_match_simple_bayesian[39][37] = -0.000325438; qual_match_simple_bayesian[39][38] = -0.000284396; qual_match_simple_bayesian[39][39] = -0.000251796; qual_match_simple_bayesian[39][40] = -0.000225901; qual_match_simple_bayesian[39][41] = -0.000205333; qual_match_simple_bayesian[39][42] = -0.000188996; qual_match_simple_bayesian[39][43] = -0.000176018; qual_match_simple_bayesian[39][44] = -0.00016571; qual_match_simple_bayesian[39][45] = -0.000157522; qual_match_simple_bayesian[39][46] = -0.000151019;
        qual_match_simple_bayesian[40][0] = -10.309; qual_match_simple_bayesian[40][1] = -1.58145; qual_match_simple_bayesian[40][2] = -0.996886; qual_match_simple_bayesian[40][3] = -0.695591; qual_match_simple_bayesian[40][4] = -0.507754; qual_match_simple_bayesian[40][5] = -0.380215; qual_match_simple_bayesian[40][6] = -0.289357; qual_match_simple_bayesian[40][7] = -0.222643; qual_match_simple_bayesian[40][8] = -0.17265; qual_match_simple_bayesian[40][9] = -0.134647; qual_match_simple_bayesian[40][10] = -0.105457; qual_match_simple_bayesian[40][11] = -0.0828624; qual_match_simple_bayesian[40][12] = -0.0652719; qual_match_simple_bayesian[40][13] = -0.0515165; qual_match_simple_bayesian[40][14] = -0.0407235; qual_match_simple_bayesian[40][15] = -0.0322325; qual_match_simple_bayesian[40][16] = -0.0255389; qual_match_simple_bayesian[40][17] = -0.0202537; qual_match_simple_bayesian[40][18] = -0.0160753; qual_match_simple_bayesian[40][19] = -0.0127688; qual_match_simple_bayesian[40][20] = -0.01015; qual_match_simple_bayesian[40][21] = -0.00807474; qual_match_simple_bayesian[40][22] = -0.00642936; qual_match_simple_bayesian[40][23] = -0.00512431; qual_match_simple_bayesian[40][24] = -0.00408889; qual_match_simple_bayesian[40][25] = -0.00326719; qual_match_simple_bayesian[40][26] = -0.00261497; qual_match_simple_bayesian[40][27] = -0.00209719; qual_match_simple_bayesian[40][28] = -0.0016861; qual_match_simple_bayesian[40][29] = -0.00135968; qual_match_simple_bayesian[40][30] = -0.00110047; qual_match_simple_bayesian[40][31] = -0.000894622; qual_match_simple_bayesian[40][32] = -0.00073114; qual_match_simple_bayesian[40][33] = -0.000601301; qual_match_simple_bayesian[40][34] = -0.000498178; qual_match_simple_bayesian[40][35] = -0.000416272; qual_match_simple_bayesian[40][36] = -0.000351217; qual_match_simple_bayesian[40][37] = -0.000299544; qual_match_simple_bayesian[40][38] = -0.000258502; qual_match_simple_bayesian[40][39] = -0.000225901; qual_match_simple_bayesian[40][40] = -0.000200007; qual_match_simple_bayesian[40][41] = -0.000179438; qual_match_simple_bayesian[40][42] = -0.000163101; qual_match_simple_bayesian[40][43] = -0.000150123; qual_match_simple_bayesian[40][44] = -0.000139815; qual_match_simple_bayesian[40][45] = -0.000131627; qual_match_simple_bayesian[40][46] = -0.000125123;
        qual_match_simple_bayesian[41][0] = -10.5392; qual_match_simple_bayesian[41][1] = -1.58145; qual_match_simple_bayesian[41][2] = -0.996877; qual_match_simple_bayesian[41][3] = -0.695577; qual_match_simple_bayesian[41][4] = -0.507738; qual_match_simple_bayesian[41][5] = -0.380198; qual_match_simple_bayesian[41][6] = -0.289339; qual_match_simple_bayesian[41][7] = -0.222624; qual_match_simple_bayesian[41][8] = -0.172631; qual_match_simple_bayesian[41][9] = -0.134628; qual_match_simple_bayesian[41][10] = -0.105437; qual_match_simple_bayesian[41][11] = -0.0828425; qual_match_simple_bayesian[41][12] = -0.0652518; qual_match_simple_bayesian[41][13] = -0.0514963; qual_match_simple_bayesian[41][14] = -0.0407032; qual_match_simple_bayesian[41][15] = -0.0322121; qual_match_simple_bayesian[41][16] = -0.0255185; qual_match_simple_bayesian[41][17] = -0.0202333; qual_match_simple_bayesian[41][18] = -0.0160549; qual_match_simple_bayesian[41][19] = -0.0127483; qual_match_simple_bayesian[41][20] = -0.0101295; qual_match_simple_bayesian[41][21] = -0.00805422; qual_match_simple_bayesian[41][22] = -0.00640883; qual_match_simple_bayesian[41][23] = -0.00510378; qual_match_simple_bayesian[41][24] = -0.00406835; qual_match_simple_bayesian[41][25] = -0.00324664; qual_match_simple_bayesian[41][26] = -0.00259442; qual_match_simple_bayesian[41][27] = -0.00207664; qual_match_simple_bayesian[41][28] = -0.00166554; qual_match_simple_bayesian[41][29] = -0.00133912; qual_match_simple_bayesian[41][30] = -0.00107991; qual_match_simple_bayesian[41][31] = -0.000874059; qual_match_simple_bayesian[41][32] = -0.000710576; qual_match_simple_bayesian[41][33] = -0.000580736; qual_match_simple_bayesian[41][34] = -0.000477612; qual_match_simple_bayesian[41][35] = -0.000395705; qual_match_simple_bayesian[41][36] = -0.00033065; qual_match_simple_bayesian[41][37] = -0.000278977; qual_match_simple_bayesian[41][38] = -0.000237934; qual_match_simple_bayesian[41][39] = -0.000205333; qual_match_simple_bayesian[41][40] = -0.000179438; qual_match_simple_bayesian[41][41] = -0.00015887; qual_match_simple_bayesian[41][42] = -0.000142532; qual_match_simple_bayesian[41][43] = -0.000129555; qual_match_simple_bayesian[41][44] = -0.000119246; qual_match_simple_bayesian[41][45] = -0.000111058; qual_match_simple_bayesian[41][46] = -0.000104554;
        qual_match_simple_bayesian[42][0] = -10.7695; qual_match_simple_bayesian[42][1] = -1.58146; qual_match_simple_bayesian[42][2] = -0.99687; qual_match_simple_bayesian[42][3] = -0.695566; qual_match_simple_bayesian[42][4] = -0.507725; qual_match_simple_bayesian[42][5] = -0.380184; qual_match_simple_bayesian[42][6] = -0.289324; qual_match_simple_bayesian[42][7] = -0.222609; qual_match_simple_bayesian[42][8] = -0.172616; qual_match_simple_bayesian[42][9] = -0.134612; qual_match_simple_bayesian[42][10] = -0.105421; qual_match_simple_bayesian[42][11] = -0.0828266; qual_match_simple_bayesian[42][12] = -0.0652359; qual_match_simple_bayesian[42][13] = -0.0514803; qual_match_simple_bayesian[42][14] = -0.0406871; qual_match_simple_bayesian[42][15] = -0.032196; qual_match_simple_bayesian[42][16] = -0.0255023; qual_match_simple_bayesian[42][17] = -0.020217; qual_match_simple_bayesian[42][18] = -0.0160386; qual_match_simple_bayesian[42][19] = -0.012732; qual_match_simple_bayesian[42][20] = -0.0101132; qual_match_simple_bayesian[42][21] = -0.00803793; qual_match_simple_bayesian[42][22] = -0.00639253; qual_match_simple_bayesian[42][23] = -0.00508747; qual_match_simple_bayesian[42][24] = -0.00405203; qual_match_simple_bayesian[42][25] = -0.00323032; qual_match_simple_bayesian[42][26] = -0.00257809; qual_match_simple_bayesian[42][27] = -0.00206031; qual_match_simple_bayesian[42][28] = -0.00164921; qual_match_simple_bayesian[42][29] = -0.00132279; qual_match_simple_bayesian[42][30] = -0.00106358; qual_match_simple_bayesian[42][31] = -0.000857725; qual_match_simple_bayesian[42][32] = -0.000694241; qual_match_simple_bayesian[42][33] = -0.0005644; qual_match_simple_bayesian[42][34] = -0.000461276; qual_match_simple_bayesian[42][35] = -0.000379369; qual_match_simple_bayesian[42][36] = -0.000314313; qual_match_simple_bayesian[42][37] = -0.00026264; qual_match_simple_bayesian[42][38] = -0.000221596; qual_match_simple_bayesian[42][39] = -0.000188996; qual_match_simple_bayesian[42][40] = -0.000163101; qual_match_simple_bayesian[42][41] = -0.000142532; qual_match_simple_bayesian[42][42] = -0.000126194; qual_match_simple_bayesian[42][43] = -0.000113217; qual_match_simple_bayesian[42][44] = -0.000102908; qual_match_simple_bayesian[42][45] = -9.47203e-05; qual_match_simple_bayesian[42][46] = -8.82164e-05;
        qual_match_simple_bayesian[43][0] = -10.9997; qual_match_simple_bayesian[43][1] = -1.58146; qual_match_simple_bayesian[43][2] = -0.996865; qual_match_simple_bayesian[43][3] = -0.695558; qual_match_simple_bayesian[43][4] = -0.507715; qual_match_simple_bayesian[43][5] = -0.380173; qual_match_simple_bayesian[43][6] = -0.289313; qual_match_simple_bayesian[43][7] = -0.222597; qual_match_simple_bayesian[43][8] = -0.172604; qual_match_simple_bayesian[43][9] = -0.1346; qual_match_simple_bayesian[43][10] = -0.105409; qual_match_simple_bayesian[43][11] = -0.082814; qual_match_simple_bayesian[43][12] = -0.0652232; qual_match_simple_bayesian[43][13] = -0.0514675; qual_match_simple_bayesian[43][14] = -0.0406743; qual_match_simple_bayesian[43][15] = -0.0321831; qual_match_simple_bayesian[43][16] = -0.0254894; qual_match_simple_bayesian[43][17] = -0.0202041; qual_match_simple_bayesian[43][18] = -0.0160257; qual_match_simple_bayesian[43][19] = -0.0127191; qual_match_simple_bayesian[43][20] = -0.0101003; qual_match_simple_bayesian[43][21] = -0.00802498; qual_match_simple_bayesian[43][22] = -0.00637958; qual_match_simple_bayesian[43][23] = -0.00507451; qual_match_simple_bayesian[43][24] = -0.00403907; qual_match_simple_bayesian[43][25] = -0.00321736; qual_match_simple_bayesian[43][26] = -0.00256512; qual_match_simple_bayesian[43][27] = -0.00204734; qual_match_simple_bayesian[43][28] = -0.00163624; qual_match_simple_bayesian[43][29] = -0.00130982; qual_match_simple_bayesian[43][30] = -0.0010506; qual_match_simple_bayesian[43][31] = -0.000844751; qual_match_simple_bayesian[43][32] = -0.000681266; qual_match_simple_bayesian[43][33] = -0.000551424; qual_match_simple_bayesian[43][34] = -0.0004483; qual_match_simple_bayesian[43][35] = -0.000366392; qual_match_simple_bayesian[43][36] = -0.000301336; qual_match_simple_bayesian[43][37] = -0.000249663; qual_match_simple_bayesian[43][38] = -0.000208619; qual_match_simple_bayesian[43][39] = -0.000176018; qual_match_simple_bayesian[43][40] = -0.000150123; qual_match_simple_bayesian[43][41] = -0.000129555; qual_match_simple_bayesian[43][42] = -0.000113217; qual_match_simple_bayesian[43][43] = -0.000100239; qual_match_simple_bayesian[43][44] = -8.99308e-05; qual_match_simple_bayesian[43][45] = -8.17427e-05; qual_match_simple_bayesian[43][46] = -7.52387e-05;
        qual_match_simple_bayesian[44][0] = -11.23; qual_match_simple_bayesian[44][1] = -1.58146; qual_match_simple_bayesian[44][2] = -0.99686; qual_match_simple_bayesian[44][3] = -0.695551; qual_match_simple_bayesian[44][4] = -0.507707; qual_match_simple_bayesian[44][5] = -0.380164; qual_match_simple_bayesian[44][6] = -0.289304; qual_match_simple_bayesian[44][7] = -0.222588; qual_match_simple_bayesian[44][8] = -0.172594; qual_match_simple_bayesian[44][9] = -0.13459; qual_match_simple_bayesian[44][10] = -0.105399; qual_match_simple_bayesian[44][11] = -0.082804; qual_match_simple_bayesian[44][12] = -0.0652131; qual_match_simple_bayesian[44][13] = -0.0514574; qual_match_simple_bayesian[44][14] = -0.0406641; qual_match_simple_bayesian[44][15] = -0.032173; qual_match_simple_bayesian[44][16] = -0.0254792; qual_match_simple_bayesian[44][17] = -0.0201939; qual_match_simple_bayesian[44][18] = -0.0160155; qual_match_simple_bayesian[44][19] = -0.0127088; qual_match_simple_bayesian[44][20] = -0.01009; qual_match_simple_bayesian[44][21] = -0.0080147; qual_match_simple_bayesian[44][22] = -0.00636929; qual_match_simple_bayesian[44][23] = -0.00506422; qual_match_simple_bayesian[44][24] = -0.00402878; qual_match_simple_bayesian[44][25] = -0.00320706; qual_match_simple_bayesian[44][26] = -0.00255482; qual_match_simple_bayesian[44][27] = -0.00203704; qual_match_simple_bayesian[44][28] = -0.00162594; qual_match_simple_bayesian[44][29] = -0.00129951; qual_match_simple_bayesian[44][30] = -0.0010403; qual_match_simple_bayesian[44][31] = -0.000834445; qual_match_simple_bayesian[44][32] = -0.00067096; qual_match_simple_bayesian[44][33] = -0.000541118; qual_match_simple_bayesian[44][34] = -0.000437993; qual_match_simple_bayesian[44][35] = -0.000356085; qual_match_simple_bayesian[44][36] = -0.000291028; qual_match_simple_bayesian[44][37] = -0.000239355; qual_match_simple_bayesian[44][38] = -0.000198311; qual_match_simple_bayesian[44][39] = -0.00016571; qual_match_simple_bayesian[44][40] = -0.000139815; qual_match_simple_bayesian[44][41] = -0.000119246; qual_match_simple_bayesian[44][42] = -0.000102908; qual_match_simple_bayesian[44][43] = -8.99308e-05; qual_match_simple_bayesian[44][44] = -7.96225e-05; qual_match_simple_bayesian[44][45] = -7.14344e-05; qual_match_simple_bayesian[44][46] = -6.49304e-05;
        qual_match_simple_bayesian[45][0] = -11.4602; qual_match_simple_bayesian[45][1] = -1.58146; qual_match_simple_bayesian[45][2] = -0.996857; qual_match_simple_bayesian[45][3] = -0.695546; qual_match_simple_bayesian[45][4] = -0.507701; qual_match_simple_bayesian[45][5] = -0.380157; qual_match_simple_bayesian[45][6] = -0.289296; qual_match_simple_bayesian[45][7] = -0.222581; qual_match_simple_bayesian[45][8] = -0.172586; qual_match_simple_bayesian[45][9] = -0.134582; qual_match_simple_bayesian[45][10] = -0.105391; qual_match_simple_bayesian[45][11] = -0.082796; qual_match_simple_bayesian[45][12] = -0.0652051; qual_match_simple_bayesian[45][13] = -0.0514493; qual_match_simple_bayesian[45][14] = -0.040656; qual_match_simple_bayesian[45][15] = -0.0321649; qual_match_simple_bayesian[45][16] = -0.0254711; qual_match_simple_bayesian[45][17] = -0.0201858; qual_match_simple_bayesian[45][18] = -0.0160073; qual_match_simple_bayesian[45][19] = -0.0127007; qual_match_simple_bayesian[45][20] = -0.0100819; qual_match_simple_bayesian[45][21] = -0.00800654; qual_match_simple_bayesian[45][22] = -0.00636112; qual_match_simple_bayesian[45][23] = -0.00505604; qual_match_simple_bayesian[45][24] = -0.0040206; qual_match_simple_bayesian[45][25] = -0.00319888; qual_match_simple_bayesian[45][26] = -0.00254664; qual_match_simple_bayesian[45][27] = -0.00202886; qual_match_simple_bayesian[45][28] = -0.00161776; qual_match_simple_bayesian[45][29] = -0.00129133; qual_match_simple_bayesian[45][30] = -0.00103211; qual_match_simple_bayesian[45][31] = -0.000826259; qual_match_simple_bayesian[45][32] = -0.000662773; qual_match_simple_bayesian[45][33] = -0.000532931; qual_match_simple_bayesian[45][34] = -0.000429806; qual_match_simple_bayesian[45][35] = -0.000347898; qual_match_simple_bayesian[45][36] = -0.000282841; qual_match_simple_bayesian[45][37] = -0.000231167; qual_match_simple_bayesian[45][38] = -0.000190123; qual_match_simple_bayesian[45][39] = -0.000157522; qual_match_simple_bayesian[45][40] = -0.000131627; qual_match_simple_bayesian[45][41] = -0.000111058; qual_match_simple_bayesian[45][42] = -9.47203e-05; qual_match_simple_bayesian[45][43] = -8.17427e-05; qual_match_simple_bayesian[45][44] = -7.14344e-05; qual_match_simple_bayesian[45][45] = -6.32462e-05; qual_match_simple_bayesian[45][46] = -5.67422e-05;
        qual_match_simple_bayesian[46][0] = -11.6905; qual_match_simple_bayesian[46][1] = -1.58147; qual_match_simple_bayesian[46][2] = -0.996854; qual_match_simple_bayesian[46][3] = -0.695541; qual_match_simple_bayesian[46][4] = -0.507695; qual_match_simple_bayesian[46][5] = -0.380152; qual_match_simple_bayesian[46][6] = -0.28929; qual_match_simple_bayesian[46][7] = -0.222575; qual_match_simple_bayesian[46][8] = -0.17258; qual_match_simple_bayesian[46][9] = -0.134576; qual_match_simple_bayesian[46][10] = -0.105385; qual_match_simple_bayesian[46][11] = -0.0827897; qual_match_simple_bayesian[46][12] = -0.0651987; qual_match_simple_bayesian[46][13] = -0.051443; qual_match_simple_bayesian[46][14] = -0.0406496; qual_match_simple_bayesian[46][15] = -0.0321584; qual_match_simple_bayesian[46][16] = -0.0254646; qual_match_simple_bayesian[46][17] = -0.0201793; qual_match_simple_bayesian[46][18] = -0.0160009; qual_match_simple_bayesian[46][19] = -0.0126942; qual_match_simple_bayesian[46][20] = -0.0100754; qual_match_simple_bayesian[46][21] = -0.00800005; qual_match_simple_bayesian[46][22] = -0.00635463; qual_match_simple_bayesian[46][23] = -0.00504955; qual_match_simple_bayesian[46][24] = -0.0040141; qual_match_simple_bayesian[46][25] = -0.00319238; qual_match_simple_bayesian[46][26] = -0.00254014; qual_match_simple_bayesian[46][27] = -0.00202236; qual_match_simple_bayesian[46][28] = -0.00161126; qual_match_simple_bayesian[46][29] = -0.00128483; qual_match_simple_bayesian[46][30] = -0.00102561; qual_match_simple_bayesian[46][31] = -0.000819756; qual_match_simple_bayesian[46][32] = -0.00065627; qual_match_simple_bayesian[46][33] = -0.000526428; qual_match_simple_bayesian[46][34] = -0.000423302; qual_match_simple_bayesian[46][35] = -0.000341394; qual_match_simple_bayesian[46][36] = -0.000276337; qual_match_simple_bayesian[46][37] = -0.000224664; qual_match_simple_bayesian[46][38] = -0.00018362; qual_match_simple_bayesian[46][39] = -0.000151019; qual_match_simple_bayesian[46][40] = -0.000125123; qual_match_simple_bayesian[46][41] = -0.000104554; qual_match_simple_bayesian[46][42] = -8.82164e-05; qual_match_simple_bayesian[46][43] = -7.52387e-05; qual_match_simple_bayesian[46][44] = -6.49304e-05; qual_match_simple_bayesian[46][45] = -5.67422e-05; qual_match_simple_bayesian[46][46] = -5.02381e-05;
        
        
        
        qual_mismatch_simple_bayesian[0][0] = -1.50408; qual_mismatch_simple_bayesian[0][1] = -1.40619; qual_mismatch_simple_bayesian[0][2] = -1.33474; qual_mismatch_simple_bayesian[0][3] = -1.28141; qual_mismatch_simple_bayesian[0][4] = -1.24099; qual_mismatch_simple_bayesian[0][5] = -1.21; qual_mismatch_simple_bayesian[0][6] = -1.18606; qual_mismatch_simple_bayesian[0][7] = -1.16744; qual_mismatch_simple_bayesian[0][8] = -1.15289; qual_mismatch_simple_bayesian[0][9] = -1.14148; qual_mismatch_simple_bayesian[0][10] = -1.13251; qual_mismatch_simple_bayesian[0][11] = -1.12545; qual_mismatch_simple_bayesian[0][12] = -1.11987; qual_mismatch_simple_bayesian[0][13] = -1.11546; qual_mismatch_simple_bayesian[0][14] = -1.11197; qual_mismatch_simple_bayesian[0][15] = -1.10921; qual_mismatch_simple_bayesian[0][16] = -1.10702; qual_mismatch_simple_bayesian[0][17] = -1.10529; qual_mismatch_simple_bayesian[0][18] = -1.10391; qual_mismatch_simple_bayesian[0][19] = -1.10282; qual_mismatch_simple_bayesian[0][20] = -1.10195; qual_mismatch_simple_bayesian[0][21] = -1.10126; qual_mismatch_simple_bayesian[0][22] = -1.10072; qual_mismatch_simple_bayesian[0][23] = -1.10028; qual_mismatch_simple_bayesian[0][24] = -1.09994; qual_mismatch_simple_bayesian[0][25] = -1.09967; qual_mismatch_simple_bayesian[0][26] = -1.09945; qual_mismatch_simple_bayesian[0][27] = -1.09928; qual_mismatch_simple_bayesian[0][28] = -1.09914; qual_mismatch_simple_bayesian[0][29] = -1.09903; qual_mismatch_simple_bayesian[0][30] = -1.09895; qual_mismatch_simple_bayesian[0][31] = -1.09888; qual_mismatch_simple_bayesian[0][32] = -1.09882; qual_mismatch_simple_bayesian[0][33] = -1.09878; qual_mismatch_simple_bayesian[0][34] = -1.09874; qual_mismatch_simple_bayesian[0][35] = -1.09872; qual_mismatch_simple_bayesian[0][36] = -1.0987; qual_mismatch_simple_bayesian[0][37] = -1.09868; qual_mismatch_simple_bayesian[0][38] = -1.09867; qual_mismatch_simple_bayesian[0][39] = -1.09865; qual_mismatch_simple_bayesian[0][40] = -1.09865; qual_mismatch_simple_bayesian[0][41] = -1.09864; qual_mismatch_simple_bayesian[0][42] = -1.09863; qual_mismatch_simple_bayesian[0][43] = -1.09863; qual_mismatch_simple_bayesian[0][44] = -1.09863; qual_mismatch_simple_bayesian[0][45] = -1.09862; qual_mismatch_simple_bayesian[0][46] = -1.09862;
        qual_mismatch_simple_bayesian[1][0] = -1.40619; qual_mismatch_simple_bayesian[1][1] = -1.38979; qual_mismatch_simple_bayesian[1][2] = -1.37696; qual_mismatch_simple_bayesian[1][3] = -1.36688; qual_mismatch_simple_bayesian[1][4] = -1.35894; qual_mismatch_simple_bayesian[1][5] = -1.35268; qual_mismatch_simple_bayesian[1][6] = -1.34774; qual_mismatch_simple_bayesian[1][7] = -1.34383; qual_mismatch_simple_bayesian[1][8] = -1.34073; qual_mismatch_simple_bayesian[1][9] = -1.33828; qual_mismatch_simple_bayesian[1][10] = -1.33634; qual_mismatch_simple_bayesian[1][11] = -1.3348; qual_mismatch_simple_bayesian[1][12] = -1.33358; qual_mismatch_simple_bayesian[1][13] = -1.33261; qual_mismatch_simple_bayesian[1][14] = -1.33184; qual_mismatch_simple_bayesian[1][15] = -1.33123; qual_mismatch_simple_bayesian[1][16] = -1.33074; qual_mismatch_simple_bayesian[1][17] = -1.33036; qual_mismatch_simple_bayesian[1][18] = -1.33005; qual_mismatch_simple_bayesian[1][19] = -1.32981; qual_mismatch_simple_bayesian[1][20] = -1.32962; qual_mismatch_simple_bayesian[1][21] = -1.32946; qual_mismatch_simple_bayesian[1][22] = -1.32934; qual_mismatch_simple_bayesian[1][23] = -1.32924; qual_mismatch_simple_bayesian[1][24] = -1.32917; qual_mismatch_simple_bayesian[1][25] = -1.32911; qual_mismatch_simple_bayesian[1][26] = -1.32906; qual_mismatch_simple_bayesian[1][27] = -1.32902; qual_mismatch_simple_bayesian[1][28] = -1.32899; qual_mismatch_simple_bayesian[1][29] = -1.32896; qual_mismatch_simple_bayesian[1][30] = -1.32895; qual_mismatch_simple_bayesian[1][31] = -1.32893; qual_mismatch_simple_bayesian[1][32] = -1.32892; qual_mismatch_simple_bayesian[1][33] = -1.32891; qual_mismatch_simple_bayesian[1][34] = -1.3289; qual_mismatch_simple_bayesian[1][35] = -1.32889; qual_mismatch_simple_bayesian[1][36] = -1.32889; qual_mismatch_simple_bayesian[1][37] = -1.32889; qual_mismatch_simple_bayesian[1][38] = -1.32888; qual_mismatch_simple_bayesian[1][39] = -1.32888; qual_mismatch_simple_bayesian[1][40] = -1.32888; qual_mismatch_simple_bayesian[1][41] = -1.32888; qual_mismatch_simple_bayesian[1][42] = -1.32888; qual_mismatch_simple_bayesian[1][43] = -1.32887; qual_mismatch_simple_bayesian[1][44] = -1.32887; qual_mismatch_simple_bayesian[1][45] = -1.32887; qual_mismatch_simple_bayesian[1][46] = -1.32887;
        qual_mismatch_simple_bayesian[2][0] = -1.33474; qual_mismatch_simple_bayesian[2][1] = -1.37696; qual_mismatch_simple_bayesian[2][2] = -1.41181; qual_mismatch_simple_bayesian[2][3] = -1.44039; qual_mismatch_simple_bayesian[2][4] = -1.46368; qual_mismatch_simple_bayesian[2][5] = -1.48258; qual_mismatch_simple_bayesian[2][6] = -1.49786; qual_mismatch_simple_bayesian[2][7] = -1.51016; qual_mismatch_simple_bayesian[2][8] = -1.52003; qual_mismatch_simple_bayesian[2][9] = -1.52795; qual_mismatch_simple_bayesian[2][10] = -1.53428; qual_mismatch_simple_bayesian[2][11] = -1.53934; qual_mismatch_simple_bayesian[2][12] = -1.54338; qual_mismatch_simple_bayesian[2][13] = -1.5466; qual_mismatch_simple_bayesian[2][14] = -1.54916; qual_mismatch_simple_bayesian[2][15] = -1.55121; qual_mismatch_simple_bayesian[2][16] = -1.55283; qual_mismatch_simple_bayesian[2][17] = -1.55412; qual_mismatch_simple_bayesian[2][18] = -1.55515; qual_mismatch_simple_bayesian[2][19] = -1.55597; qual_mismatch_simple_bayesian[2][20] = -1.55662; qual_mismatch_simple_bayesian[2][21] = -1.55713; qual_mismatch_simple_bayesian[2][22] = -1.55754; qual_mismatch_simple_bayesian[2][23] = -1.55787; qual_mismatch_simple_bayesian[2][24] = -1.55813; qual_mismatch_simple_bayesian[2][25] = -1.55833; qual_mismatch_simple_bayesian[2][26] = -1.5585; qual_mismatch_simple_bayesian[2][27] = -1.55863; qual_mismatch_simple_bayesian[2][28] = -1.55873; qual_mismatch_simple_bayesian[2][29] = -1.55881; qual_mismatch_simple_bayesian[2][30] = -1.55888; qual_mismatch_simple_bayesian[2][31] = -1.55893; qual_mismatch_simple_bayesian[2][32] = -1.55897; qual_mismatch_simple_bayesian[2][33] = -1.559; qual_mismatch_simple_bayesian[2][34] = -1.55903; qual_mismatch_simple_bayesian[2][35] = -1.55905; qual_mismatch_simple_bayesian[2][36] = -1.55907; qual_mismatch_simple_bayesian[2][37] = -1.55908; qual_mismatch_simple_bayesian[2][38] = -1.55909; qual_mismatch_simple_bayesian[2][39] = -1.5591; qual_mismatch_simple_bayesian[2][40] = -1.5591; qual_mismatch_simple_bayesian[2][41] = -1.55911; qual_mismatch_simple_bayesian[2][42] = -1.55911; qual_mismatch_simple_bayesian[2][43] = -1.55912; qual_mismatch_simple_bayesian[2][44] = -1.55912; qual_mismatch_simple_bayesian[2][45] = -1.55912; qual_mismatch_simple_bayesian[2][46] = -1.55912;
        qual_mismatch_simple_bayesian[3][0] = -1.28141; qual_mismatch_simple_bayesian[3][1] = -1.36688; qual_mismatch_simple_bayesian[3][2] = -1.44039; qual_mismatch_simple_bayesian[3][3] = -1.50289; qual_mismatch_simple_bayesian[3][4] = -1.55549; qual_mismatch_simple_bayesian[3][5] = -1.59933; qual_mismatch_simple_bayesian[3][6] = -1.63558; qual_mismatch_simple_bayesian[3][7] = -1.66534; qual_mismatch_simple_bayesian[3][8] = -1.68963; qual_mismatch_simple_bayesian[3][9] = -1.70935; qual_mismatch_simple_bayesian[3][10] = -1.72529; qual_mismatch_simple_bayesian[3][11] = -1.73814; qual_mismatch_simple_bayesian[3][12] = -1.74847; qual_mismatch_simple_bayesian[3][13] = -1.75675; qual_mismatch_simple_bayesian[3][14] = -1.76338; qual_mismatch_simple_bayesian[3][15] = -1.76867; qual_mismatch_simple_bayesian[3][16] = -1.7729; qual_mismatch_simple_bayesian[3][17] = -1.77627; qual_mismatch_simple_bayesian[3][18] = -1.77895; qual_mismatch_simple_bayesian[3][19] = -1.78109; qual_mismatch_simple_bayesian[3][20] = -1.78279; qual_mismatch_simple_bayesian[3][21] = -1.78414; qual_mismatch_simple_bayesian[3][22] = -1.78522; qual_mismatch_simple_bayesian[3][23] = -1.78608; qual_mismatch_simple_bayesian[3][24] = -1.78676; qual_mismatch_simple_bayesian[3][25] = -1.7873; qual_mismatch_simple_bayesian[3][26] = -1.78773; qual_mismatch_simple_bayesian[3][27] = -1.78807; qual_mismatch_simple_bayesian[3][28] = -1.78834; qual_mismatch_simple_bayesian[3][29] = -1.78855; qual_mismatch_simple_bayesian[3][30] = -1.78873; qual_mismatch_simple_bayesian[3][31] = -1.78886; qual_mismatch_simple_bayesian[3][32] = -1.78897; qual_mismatch_simple_bayesian[3][33] = -1.78906; qual_mismatch_simple_bayesian[3][34] = -1.78912; qual_mismatch_simple_bayesian[3][35] = -1.78918; qual_mismatch_simple_bayesian[3][36] = -1.78922; qual_mismatch_simple_bayesian[3][37] = -1.78926; qual_mismatch_simple_bayesian[3][38] = -1.78928; qual_mismatch_simple_bayesian[3][39] = -1.7893; qual_mismatch_simple_bayesian[3][40] = -1.78932; qual_mismatch_simple_bayesian[3][41] = -1.78934; qual_mismatch_simple_bayesian[3][42] = -1.78935; qual_mismatch_simple_bayesian[3][43] = -1.78935; qual_mismatch_simple_bayesian[3][44] = -1.78936; qual_mismatch_simple_bayesian[3][45] = -1.78937; qual_mismatch_simple_bayesian[3][46] = -1.78937;
        qual_mismatch_simple_bayesian[4][0] = -1.24099; qual_mismatch_simple_bayesian[4][1] = -1.35894; qual_mismatch_simple_bayesian[4][2] = -1.46368; qual_mismatch_simple_bayesian[4][3] = -1.55549; qual_mismatch_simple_bayesian[4][4] = -1.63493; qual_mismatch_simple_bayesian[4][5] = -1.70287; qual_mismatch_simple_bayesian[4][6] = -1.76033; qual_mismatch_simple_bayesian[4][7] = -1.80845; qual_mismatch_simple_bayesian[4][8] = -1.8484; qual_mismatch_simple_bayesian[4][9] = -1.8813; qual_mismatch_simple_bayesian[4][10] = -1.90823; qual_mismatch_simple_bayesian[4][11] = -1.93016; qual_mismatch_simple_bayesian[4][12] = -1.94792; qual_mismatch_simple_bayesian[4][13] = -1.96226; qual_mismatch_simple_bayesian[4][14] = -1.97379; qual_mismatch_simple_bayesian[4][15] = -1.98305; qual_mismatch_simple_bayesian[4][16] = -1.99047; qual_mismatch_simple_bayesian[4][17] = -1.9964; qual_mismatch_simple_bayesian[4][18] = -2.00114; qual_mismatch_simple_bayesian[4][19] = -2.00492; qual_mismatch_simple_bayesian[4][20] = -2.00793; qual_mismatch_simple_bayesian[4][21] = -2.01033; qual_mismatch_simple_bayesian[4][22] = -2.01224; qual_mismatch_simple_bayesian[4][23] = -2.01376; qual_mismatch_simple_bayesian[4][24] = -2.01497; qual_mismatch_simple_bayesian[4][25] = -2.01593; qual_mismatch_simple_bayesian[4][26] = -2.01669; qual_mismatch_simple_bayesian[4][27] = -2.0173; qual_mismatch_simple_bayesian[4][28] = -2.01778; qual_mismatch_simple_bayesian[4][29] = -2.01816; qual_mismatch_simple_bayesian[4][30] = -2.01847; qual_mismatch_simple_bayesian[4][31] = -2.01871; qual_mismatch_simple_bayesian[4][32] = -2.0189; qual_mismatch_simple_bayesian[4][33] = -2.01906; qual_mismatch_simple_bayesian[4][34] = -2.01918; qual_mismatch_simple_bayesian[4][35] = -2.01927; qual_mismatch_simple_bayesian[4][36] = -2.01935; qual_mismatch_simple_bayesian[4][37] = -2.01941; qual_mismatch_simple_bayesian[4][38] = -2.01946; qual_mismatch_simple_bayesian[4][39] = -2.0195; qual_mismatch_simple_bayesian[4][40] = -2.01953; qual_mismatch_simple_bayesian[4][41] = -2.01955; qual_mismatch_simple_bayesian[4][42] = -2.01957; qual_mismatch_simple_bayesian[4][43] = -2.01959; qual_mismatch_simple_bayesian[4][44] = -2.0196; qual_mismatch_simple_bayesian[4][45] = -2.01961; qual_mismatch_simple_bayesian[4][46] = -2.01962;
        qual_mismatch_simple_bayesian[5][0] = -1.21; qual_mismatch_simple_bayesian[5][1] = -1.35268; qual_mismatch_simple_bayesian[5][2] = -1.48258; qual_mismatch_simple_bayesian[5][3] = -1.59933; qual_mismatch_simple_bayesian[5][4] = -1.70287; qual_mismatch_simple_bayesian[5][5] = -1.79352; qual_mismatch_simple_bayesian[5][6] = -1.87187; qual_mismatch_simple_bayesian[5][7] = -1.93881; qual_mismatch_simple_bayesian[5][8] = -1.99536; qual_mismatch_simple_bayesian[5][9] = -2.04269; qual_mismatch_simple_bayesian[5][10] = -2.08194; qual_mismatch_simple_bayesian[5][11] = -2.11426; qual_mismatch_simple_bayesian[5][12] = -2.14069; qual_mismatch_simple_bayesian[5][13] = -2.1622; qual_mismatch_simple_bayesian[5][14] = -2.17962; qual_mismatch_simple_bayesian[5][15] = -2.19368; qual_mismatch_simple_bayesian[5][16] = -2.20499; qual_mismatch_simple_bayesian[5][17] = -2.21406; qual_mismatch_simple_bayesian[5][18] = -2.22133; qual_mismatch_simple_bayesian[5][19] = -2.22714; qual_mismatch_simple_bayesian[5][20] = -2.23178; qual_mismatch_simple_bayesian[5][21] = -2.23548; qual_mismatch_simple_bayesian[5][22] = -2.23843; qual_mismatch_simple_bayesian[5][23] = -2.24078; qual_mismatch_simple_bayesian[5][24] = -2.24265; qual_mismatch_simple_bayesian[5][25] = -2.24414; qual_mismatch_simple_bayesian[5][26] = -2.24532; qual_mismatch_simple_bayesian[5][27] = -2.24626; qual_mismatch_simple_bayesian[5][28] = -2.24701; qual_mismatch_simple_bayesian[5][29] = -2.2476; qual_mismatch_simple_bayesian[5][30] = -2.24808; qual_mismatch_simple_bayesian[5][31] = -2.24845; qual_mismatch_simple_bayesian[5][32] = -2.24875; qual_mismatch_simple_bayesian[5][33] = -2.24899; qual_mismatch_simple_bayesian[5][34] = -2.24918; qual_mismatch_simple_bayesian[5][35] = -2.24933; qual_mismatch_simple_bayesian[5][36] = -2.24945; qual_mismatch_simple_bayesian[5][37] = -2.24954; qual_mismatch_simple_bayesian[5][38] = -2.24962; qual_mismatch_simple_bayesian[5][39] = -2.24967; qual_mismatch_simple_bayesian[5][40] = -2.24972; qual_mismatch_simple_bayesian[5][41] = -2.24976; qual_mismatch_simple_bayesian[5][42] = -2.24979; qual_mismatch_simple_bayesian[5][43] = -2.24981; qual_mismatch_simple_bayesian[5][44] = -2.24983; qual_mismatch_simple_bayesian[5][45] = -2.24985; qual_mismatch_simple_bayesian[5][46] = -2.24986;
        qual_mismatch_simple_bayesian[6][0] = -1.18606; qual_mismatch_simple_bayesian[6][1] = -1.34774; qual_mismatch_simple_bayesian[6][2] = -1.49786; qual_mismatch_simple_bayesian[6][3] = -1.63558; qual_mismatch_simple_bayesian[6][4] = -1.76033; qual_mismatch_simple_bayesian[6][5] = -1.87187; qual_mismatch_simple_bayesian[6][6] = -1.97029; qual_mismatch_simple_bayesian[6][7] = -2.05601; qual_mismatch_simple_bayesian[6][8] = -2.12976; qual_mismatch_simple_bayesian[6][9] = -2.19248; qual_mismatch_simple_bayesian[6][10] = -2.24527; qual_mismatch_simple_bayesian[6][11] = -2.28928; qual_mismatch_simple_bayesian[6][12] = -2.32567; qual_mismatch_simple_bayesian[6][13] = -2.35556; qual_mismatch_simple_bayesian[6][14] = -2.37995; qual_mismatch_simple_bayesian[6][15] = -2.39976; qual_mismatch_simple_bayesian[6][16] = -2.41577; qual_mismatch_simple_bayesian[6][17] = -2.42868; qual_mismatch_simple_bayesian[6][18] = -2.43906; qual_mismatch_simple_bayesian[6][19] = -2.44737; qual_mismatch_simple_bayesian[6][20] = -2.45403; qual_mismatch_simple_bayesian[6][21] = -2.45935; qual_mismatch_simple_bayesian[6][22] = -2.4636; qual_mismatch_simple_bayesian[6][23] = -2.46698; qual_mismatch_simple_bayesian[6][24] = -2.46968; qual_mismatch_simple_bayesian[6][25] = -2.47183; qual_mismatch_simple_bayesian[6][26] = -2.47353; qual_mismatch_simple_bayesian[6][27] = -2.47489; qual_mismatch_simple_bayesian[6][28] = -2.47598; qual_mismatch_simple_bayesian[6][29] = -2.47684; qual_mismatch_simple_bayesian[6][30] = -2.47752; qual_mismatch_simple_bayesian[6][31] = -2.47806; qual_mismatch_simple_bayesian[6][32] = -2.47849; qual_mismatch_simple_bayesian[6][33] = -2.47884; qual_mismatch_simple_bayesian[6][34] = -2.47911; qual_mismatch_simple_bayesian[6][35] = -2.47933; qual_mismatch_simple_bayesian[6][36] = -2.4795; qual_mismatch_simple_bayesian[6][37] = -2.47964; qual_mismatch_simple_bayesian[6][38] = -2.47974; qual_mismatch_simple_bayesian[6][39] = -2.47983; qual_mismatch_simple_bayesian[6][40] = -2.4799; qual_mismatch_simple_bayesian[6][41] = -2.47995; qual_mismatch_simple_bayesian[6][42] = -2.48; qual_mismatch_simple_bayesian[6][43] = -2.48003; qual_mismatch_simple_bayesian[6][44] = -2.48006; qual_mismatch_simple_bayesian[6][45] = -2.48008; qual_mismatch_simple_bayesian[6][46] = -2.4801;
        qual_mismatch_simple_bayesian[7][0] = -1.16744; qual_mismatch_simple_bayesian[7][1] = -1.34383; qual_mismatch_simple_bayesian[7][2] = -1.51016; qual_mismatch_simple_bayesian[7][3] = -1.66534; qual_mismatch_simple_bayesian[7][4] = -1.80845; qual_mismatch_simple_bayesian[7][5] = -1.93881; qual_mismatch_simple_bayesian[7][6] = -2.05601; qual_mismatch_simple_bayesian[7][7] = -2.16001; qual_mismatch_simple_bayesian[7][8] = -2.25109; qual_mismatch_simple_bayesian[7][9] = -2.32986; qual_mismatch_simple_bayesian[7][10] = -2.39718; qual_mismatch_simple_bayesian[7][11] = -2.45408; qual_mismatch_simple_bayesian[7][12] = -2.5017; qual_mismatch_simple_bayesian[7][13] = -2.54122; qual_mismatch_simple_bayesian[7][14] = -2.57376; qual_mismatch_simple_bayesian[7][15] = -2.60038; qual_mismatch_simple_bayesian[7][16] = -2.62204; qual_mismatch_simple_bayesian[7][17] = -2.63959; qual_mismatch_simple_bayesian[7][18] = -2.65376; qual_mismatch_simple_bayesian[7][19] = -2.66515; qual_mismatch_simple_bayesian[7][20] = -2.6743; qual_mismatch_simple_bayesian[7][21] = -2.68162; qual_mismatch_simple_bayesian[7][22] = -2.68748; qual_mismatch_simple_bayesian[7][23] = -2.69215; qual_mismatch_simple_bayesian[7][24] = -2.69588; qual_mismatch_simple_bayesian[7][25] = -2.69886; qual_mismatch_simple_bayesian[7][26] = -2.70122; qual_mismatch_simple_bayesian[7][27] = -2.70311; qual_mismatch_simple_bayesian[7][28] = -2.70461; qual_mismatch_simple_bayesian[7][29] = -2.7058; qual_mismatch_simple_bayesian[7][30] = -2.70675; qual_mismatch_simple_bayesian[7][31] = -2.7075; qual_mismatch_simple_bayesian[7][32] = -2.7081; qual_mismatch_simple_bayesian[7][33] = -2.70858; qual_mismatch_simple_bayesian[7][34] = -2.70896; qual_mismatch_simple_bayesian[7][35] = -2.70926; qual_mismatch_simple_bayesian[7][36] = -2.7095; qual_mismatch_simple_bayesian[7][37] = -2.70969; qual_mismatch_simple_bayesian[7][38] = -2.70984; qual_mismatch_simple_bayesian[7][39] = -2.70996; qual_mismatch_simple_bayesian[7][40] = -2.71005; qual_mismatch_simple_bayesian[7][41] = -2.71013; qual_mismatch_simple_bayesian[7][42] = -2.71019; qual_mismatch_simple_bayesian[7][43] = -2.71024; qual_mismatch_simple_bayesian[7][44] = -2.71028; qual_mismatch_simple_bayesian[7][45] = -2.71031; qual_mismatch_simple_bayesian[7][46] = -2.71033;
        qual_mismatch_simple_bayesian[8][0] = -1.15289; qual_mismatch_simple_bayesian[8][1] = -1.34073; qual_mismatch_simple_bayesian[8][2] = -1.52003; qual_mismatch_simple_bayesian[8][3] = -1.68963; qual_mismatch_simple_bayesian[8][4] = -1.8484; qual_mismatch_simple_bayesian[8][5] = -1.99536; qual_mismatch_simple_bayesian[8][6] = -2.12976; qual_mismatch_simple_bayesian[8][7] = -2.25109; qual_mismatch_simple_bayesian[8][8] = -2.3592; qual_mismatch_simple_bayesian[8][9] = -2.45427; qual_mismatch_simple_bayesian[8][10] = -2.5368; qual_mismatch_simple_bayesian[8][11] = -2.60759; qual_mismatch_simple_bayesian[8][12] = -2.66762; qual_mismatch_simple_bayesian[8][13] = -2.71801; qual_mismatch_simple_bayesian[8][14] = -2.75994; qual_mismatch_simple_bayesian[8][15] = -2.79454; qual_mismatch_simple_bayesian[8][16] = -2.8229; qual_mismatch_simple_bayesian[8][17] = -2.84602; qual_mismatch_simple_bayesian[8][18] = -2.86477; qual_mismatch_simple_bayesian[8][19] = -2.87992; qual_mismatch_simple_bayesian[8][20] = -2.89212; qual_mismatch_simple_bayesian[8][21] = -2.90191; qual_mismatch_simple_bayesian[8][22] = -2.90977; qual_mismatch_simple_bayesian[8][23] = -2.91605; qual_mismatch_simple_bayesian[8][24] = -2.92106; qual_mismatch_simple_bayesian[8][25] = -2.92507; qual_mismatch_simple_bayesian[8][26] = -2.92826; qual_mismatch_simple_bayesian[8][27] = -2.9308; qual_mismatch_simple_bayesian[8][28] = -2.93282; qual_mismatch_simple_bayesian[8][29] = -2.93444; qual_mismatch_simple_bayesian[8][30] = -2.93572; qual_mismatch_simple_bayesian[8][31] = -2.93674; qual_mismatch_simple_bayesian[8][32] = -2.93755; qual_mismatch_simple_bayesian[8][33] = -2.93819; qual_mismatch_simple_bayesian[8][34] = -2.9387; qual_mismatch_simple_bayesian[8][35] = -2.93911; qual_mismatch_simple_bayesian[8][36] = -2.93943; qual_mismatch_simple_bayesian[8][37] = -2.93969; qual_mismatch_simple_bayesian[8][38] = -2.93989; qual_mismatch_simple_bayesian[8][39] = -2.94005; qual_mismatch_simple_bayesian[8][40] = -2.94018; qual_mismatch_simple_bayesian[8][41] = -2.94029; qual_mismatch_simple_bayesian[8][42] = -2.94037; qual_mismatch_simple_bayesian[8][43] = -2.94043; qual_mismatch_simple_bayesian[8][44] = -2.94048; qual_mismatch_simple_bayesian[8][45] = -2.94052; qual_mismatch_simple_bayesian[8][46] = -2.94056;
        qual_mismatch_simple_bayesian[9][0] = -1.14148; qual_mismatch_simple_bayesian[9][1] = -1.33828; qual_mismatch_simple_bayesian[9][2] = -1.52795; qual_mismatch_simple_bayesian[9][3] = -1.70935; qual_mismatch_simple_bayesian[9][4] = -1.8813; qual_mismatch_simple_bayesian[9][5] = -2.04269; qual_mismatch_simple_bayesian[9][6] = -2.19248; qual_mismatch_simple_bayesian[9][7] = -2.32986; qual_mismatch_simple_bayesian[9][8] = -2.45427; qual_mismatch_simple_bayesian[9][9] = -2.56545; qual_mismatch_simple_bayesian[9][10] = -2.66352; qual_mismatch_simple_bayesian[9][11] = -2.74891; qual_mismatch_simple_bayesian[9][12] = -2.82235; qual_mismatch_simple_bayesian[9][13] = -2.8848; qual_mismatch_simple_bayesian[9][14] = -2.93733; qual_mismatch_simple_bayesian[9][15] = -2.98112; qual_mismatch_simple_bayesian[9][16] = -3.01733; qual_mismatch_simple_bayesian[9][17] = -3.04705; qual_mismatch_simple_bayesian[9][18] = -3.07131; qual_mismatch_simple_bayesian[9][19] = -3.09101; qual_mismatch_simple_bayesian[9][20] = -3.10693; qual_mismatch_simple_bayesian[9][21] = -3.11977; qual_mismatch_simple_bayesian[9][22] = -3.13008; qual_mismatch_simple_bayesian[9][23] = -3.13835; qual_mismatch_simple_bayesian[9][24] = -3.14496; qual_mismatch_simple_bayesian[9][25] = -3.15025; qual_mismatch_simple_bayesian[9][26] = -3.15447; qual_mismatch_simple_bayesian[9][27] = -3.15784; qual_mismatch_simple_bayesian[9][28] = -3.16052; qual_mismatch_simple_bayesian[9][29] = -3.16265; qual_mismatch_simple_bayesian[9][30] = -3.16435; qual_mismatch_simple_bayesian[9][31] = -3.1657; qual_mismatch_simple_bayesian[9][32] = -3.16678; qual_mismatch_simple_bayesian[9][33] = -3.16763; qual_mismatch_simple_bayesian[9][34] = -3.16831; qual_mismatch_simple_bayesian[9][35] = -3.16885; qual_mismatch_simple_bayesian[9][36] = -3.16928; qual_mismatch_simple_bayesian[9][37] = -3.16962; qual_mismatch_simple_bayesian[9][38] = -3.16989; qual_mismatch_simple_bayesian[9][39] = -3.17011; qual_mismatch_simple_bayesian[9][40] = -3.17028; qual_mismatch_simple_bayesian[9][41] = -3.17041; qual_mismatch_simple_bayesian[9][42] = -3.17052; qual_mismatch_simple_bayesian[9][43] = -3.17061; qual_mismatch_simple_bayesian[9][44] = -3.17068; qual_mismatch_simple_bayesian[9][45] = -3.17073; qual_mismatch_simple_bayesian[9][46] = -3.17077;
        qual_mismatch_simple_bayesian[10][0] = -1.13251; qual_mismatch_simple_bayesian[10][1] = -1.33634; qual_mismatch_simple_bayesian[10][2] = -1.53428; qual_mismatch_simple_bayesian[10][3] = -1.72529; qual_mismatch_simple_bayesian[10][4] = -1.90823; qual_mismatch_simple_bayesian[10][5] = -2.08194; qual_mismatch_simple_bayesian[10][6] = -2.24527; qual_mismatch_simple_bayesian[10][7] = -2.39718; qual_mismatch_simple_bayesian[10][8] = -2.5368; qual_mismatch_simple_bayesian[10][9] = -2.66352; qual_mismatch_simple_bayesian[10][10] = -2.77704; qual_mismatch_simple_bayesian[10][11] = -2.87741; qual_mismatch_simple_bayesian[10][12] = -2.96499; qual_mismatch_simple_bayesian[10][13] = -3.04048; qual_mismatch_simple_bayesian[10][14] = -3.10478; qual_mismatch_simple_bayesian[10][15] = -3.15899; qual_mismatch_simple_bayesian[10][16] = -3.20424; qual_mismatch_simple_bayesian[10][17] = -3.2417; qual_mismatch_simple_bayesian[10][18] = -3.27249; qual_mismatch_simple_bayesian[10][19] = -3.29764; qual_mismatch_simple_bayesian[10][20] = -3.31808; qual_mismatch_simple_bayesian[10][21] = -3.33462; qual_mismatch_simple_bayesian[10][22] = -3.34796; qual_mismatch_simple_bayesian[10][23] = -3.35868; qual_mismatch_simple_bayesian[10][24] = -3.36728; qual_mismatch_simple_bayesian[10][25] = -3.37416; qual_mismatch_simple_bayesian[10][26] = -3.37966; qual_mismatch_simple_bayesian[10][27] = -3.38405; qual_mismatch_simple_bayesian[10][28] = -3.38756; qual_mismatch_simple_bayesian[10][29] = -3.39035; qual_mismatch_simple_bayesian[10][30] = -3.39257; qual_mismatch_simple_bayesian[10][31] = -3.39434; qual_mismatch_simple_bayesian[10][32] = -3.39574; qual_mismatch_simple_bayesian[10][33] = -3.39686; qual_mismatch_simple_bayesian[10][34] = -3.39775; qual_mismatch_simple_bayesian[10][35] = -3.39846; qual_mismatch_simple_bayesian[10][36] = -3.39902; qual_mismatch_simple_bayesian[10][37] = -3.39947; qual_mismatch_simple_bayesian[10][38] = -3.39982; qual_mismatch_simple_bayesian[10][39] = -3.40011; qual_mismatch_simple_bayesian[10][40] = -3.40033; qual_mismatch_simple_bayesian[10][41] = -3.40051; qual_mismatch_simple_bayesian[10][42] = -3.40065; qual_mismatch_simple_bayesian[10][43] = -3.40076; qual_mismatch_simple_bayesian[10][44] = -3.40085; qual_mismatch_simple_bayesian[10][45] = -3.40092; qual_mismatch_simple_bayesian[10][46] = -3.40098;
        qual_mismatch_simple_bayesian[11][0] = -1.12545; qual_mismatch_simple_bayesian[11][1] = -1.3348; qual_mismatch_simple_bayesian[11][2] = -1.53934; qual_mismatch_simple_bayesian[11][3] = -1.73814; qual_mismatch_simple_bayesian[11][4] = -1.93016; qual_mismatch_simple_bayesian[11][5] = -2.11426; qual_mismatch_simple_bayesian[11][6] = -2.28928; qual_mismatch_simple_bayesian[11][7] = -2.45408; qual_mismatch_simple_bayesian[11][8] = -2.60759; qual_mismatch_simple_bayesian[11][9] = -2.74891; qual_mismatch_simple_bayesian[11][10] = -2.87741; qual_mismatch_simple_bayesian[11][11] = -2.99272; qual_mismatch_simple_bayesian[11][12] = -3.09485; qual_mismatch_simple_bayesian[11][13] = -3.18412; qual_mismatch_simple_bayesian[11][14] = -3.2612; qual_mismatch_simple_bayesian[11][15] = -3.32696; qual_mismatch_simple_bayesian[11][16] = -3.38246; qual_mismatch_simple_bayesian[11][17] = -3.42885; qual_mismatch_simple_bayesian[11][18] = -3.4673; qual_mismatch_simple_bayesian[11][19] = -3.49893; qual_mismatch_simple_bayesian[11][20] = -3.52479; qual_mismatch_simple_bayesian[11][21] = -3.54582; qual_mismatch_simple_bayesian[11][22] = -3.56284; qual_mismatch_simple_bayesian[11][23] = -3.57658; qual_mismatch_simple_bayesian[11][24] = -3.58762; qual_mismatch_simple_bayesian[11][25] = -3.59648; qual_mismatch_simple_bayesian[11][26] = -3.60357; qual_mismatch_simple_bayesian[11][27] = -3.60925; qual_mismatch_simple_bayesian[11][28] = -3.61377; qual_mismatch_simple_bayesian[11][29] = -3.61738; qual_mismatch_simple_bayesian[11][30] = -3.62026; qual_mismatch_simple_bayesian[11][31] = -3.62255; qual_mismatch_simple_bayesian[11][32] = -3.62438; qual_mismatch_simple_bayesian[11][33] = -3.62583; qual_mismatch_simple_bayesian[11][34] = -3.62698; qual_mismatch_simple_bayesian[11][35] = -3.6279; qual_mismatch_simple_bayesian[11][36] = -3.62863; qual_mismatch_simple_bayesian[11][37] = -3.62921; qual_mismatch_simple_bayesian[11][38] = -3.62967; qual_mismatch_simple_bayesian[11][39] = -3.63004; qual_mismatch_simple_bayesian[11][40] = -3.63033; qual_mismatch_simple_bayesian[11][41] = -3.63056; qual_mismatch_simple_bayesian[11][42] = -3.63075; qual_mismatch_simple_bayesian[11][43] = -3.63089; qual_mismatch_simple_bayesian[11][44] = -3.63101; qual_mismatch_simple_bayesian[11][45] = -3.6311; qual_mismatch_simple_bayesian[11][46] = -3.63117;
        qual_mismatch_simple_bayesian[12][0] = -1.11987; qual_mismatch_simple_bayesian[12][1] = -1.33358; qual_mismatch_simple_bayesian[12][2] = -1.54338; qual_mismatch_simple_bayesian[12][3] = -1.74847; qual_mismatch_simple_bayesian[12][4] = -1.94792; qual_mismatch_simple_bayesian[12][5] = -2.14069; qual_mismatch_simple_bayesian[12][6] = -2.32567; qual_mismatch_simple_bayesian[12][7] = -2.5017; qual_mismatch_simple_bayesian[12][8] = -2.66762; qual_mismatch_simple_bayesian[12][9] = -2.82235; qual_mismatch_simple_bayesian[12][10] = -2.96499; qual_mismatch_simple_bayesian[12][11] = -3.09485; qual_mismatch_simple_bayesian[12][12] = -3.21154; qual_mismatch_simple_bayesian[12][13] = -3.31504; qual_mismatch_simple_bayesian[12][14] = -3.40563; qual_mismatch_simple_bayesian[12][15] = -3.48395; qual_mismatch_simple_bayesian[12][16] = -3.55084; qual_mismatch_simple_bayesian[12][17] = -3.60736; qual_mismatch_simple_bayesian[12][18] = -3.65465; qual_mismatch_simple_bayesian[12][19] = -3.69388; qual_mismatch_simple_bayesian[12][20] = -3.72617; qual_mismatch_simple_bayesian[12][21] = -3.75259; qual_mismatch_simple_bayesian[12][22] = -3.77408; qual_mismatch_simple_bayesian[12][23] = -3.79149; qual_mismatch_simple_bayesian[12][24] = -3.80553; qual_mismatch_simple_bayesian[12][25] = -3.81683; qual_mismatch_simple_bayesian[12][26] = -3.8259; qual_mismatch_simple_bayesian[12][27] = -3.83316; qual_mismatch_simple_bayesian[12][28] = -3.83897; qual_mismatch_simple_bayesian[12][29] = -3.84361; qual_mismatch_simple_bayesian[12][30] = -3.8473; qual_mismatch_simple_bayesian[12][31] = -3.85025; qual_mismatch_simple_bayesian[12][32] = -3.8526; qual_mismatch_simple_bayesian[12][33] = -3.85447; qual_mismatch_simple_bayesian[12][34] = -3.85595; qual_mismatch_simple_bayesian[12][35] = -3.85713; qual_mismatch_simple_bayesian[12][36] = -3.85807; qual_mismatch_simple_bayesian[12][37] = -3.85882; qual_mismatch_simple_bayesian[12][38] = -3.85942; qual_mismatch_simple_bayesian[12][39] = -3.85989; qual_mismatch_simple_bayesian[12][40] = -3.86026; qual_mismatch_simple_bayesian[12][41] = -3.86056; qual_mismatch_simple_bayesian[12][42] = -3.8608; qual_mismatch_simple_bayesian[12][43] = -3.86099; qual_mismatch_simple_bayesian[12][44] = -3.86114; qual_mismatch_simple_bayesian[12][45] = -3.86126; qual_mismatch_simple_bayesian[12][46] = -3.86135;
        qual_mismatch_simple_bayesian[13][0] = -1.11546; qual_mismatch_simple_bayesian[13][1] = -1.33261; qual_mismatch_simple_bayesian[13][2] = -1.5466; qual_mismatch_simple_bayesian[13][3] = -1.75675; qual_mismatch_simple_bayesian[13][4] = -1.96226; qual_mismatch_simple_bayesian[13][5] = -2.1622; qual_mismatch_simple_bayesian[13][6] = -2.35556; qual_mismatch_simple_bayesian[13][7] = -2.54122; qual_mismatch_simple_bayesian[13][8] = -2.71801; qual_mismatch_simple_bayesian[13][9] = -2.8848; qual_mismatch_simple_bayesian[13][10] = -3.04048; qual_mismatch_simple_bayesian[13][11] = -3.18412; qual_mismatch_simple_bayesian[13][12] = -3.31504; qual_mismatch_simple_bayesian[13][13] = -3.43281; qual_mismatch_simple_bayesian[13][14] = -3.53737; qual_mismatch_simple_bayesian[13][15] = -3.629; qual_mismatch_simple_bayesian[13][16] = -3.70828; qual_mismatch_simple_bayesian[13][17] = -3.77607; qual_mismatch_simple_bayesian[13][18] = -3.83339; qual_mismatch_simple_bayesian[13][19] = -3.88139; qual_mismatch_simple_bayesian[13][20] = -3.92122; qual_mismatch_simple_bayesian[13][21] = -3.95404; qual_mismatch_simple_bayesian[13][22] = -3.9809; qual_mismatch_simple_bayesian[13][23] = -4.00276; qual_mismatch_simple_bayesian[13][24] = -4.02047; qual_mismatch_simple_bayesian[13][25] = -4.03476; qual_mismatch_simple_bayesian[13][26] = -4.04626; qual_mismatch_simple_bayesian[13][27] = -4.0555; qual_mismatch_simple_bayesian[13][28] = -4.06289; qual_mismatch_simple_bayesian[13][29] = -4.0688; qual_mismatch_simple_bayesian[13][30] = -4.07352; qual_mismatch_simple_bayesian[13][31] = -4.07729; qual_mismatch_simple_bayesian[13][32] = -4.08029; qual_mismatch_simple_bayesian[13][33] = -4.08268; qual_mismatch_simple_bayesian[13][34] = -4.08459; qual_mismatch_simple_bayesian[13][35] = -4.0861; qual_mismatch_simple_bayesian[13][36] = -4.08731; qual_mismatch_simple_bayesian[13][37] = -4.08826; qual_mismatch_simple_bayesian[13][38] = -4.08903; qual_mismatch_simple_bayesian[13][39] = -4.08963; qual_mismatch_simple_bayesian[13][40] = -4.09011; qual_mismatch_simple_bayesian[13][41] = -4.0905; qual_mismatch_simple_bayesian[13][42] = -4.0908; qual_mismatch_simple_bayesian[13][43] = -4.09104; qual_mismatch_simple_bayesian[13][44] = -4.09123; qual_mismatch_simple_bayesian[13][45] = -4.09138; qual_mismatch_simple_bayesian[13][46] = -4.09151;
        qual_mismatch_simple_bayesian[14][0] = -1.11197; qual_mismatch_simple_bayesian[14][1] = -1.33184; qual_mismatch_simple_bayesian[14][2] = -1.54916; qual_mismatch_simple_bayesian[14][3] = -1.76338; qual_mismatch_simple_bayesian[14][4] = -1.97379; qual_mismatch_simple_bayesian[14][5] = -2.17962; qual_mismatch_simple_bayesian[14][6] = -2.37995; qual_mismatch_simple_bayesian[14][7] = -2.57376; qual_mismatch_simple_bayesian[14][8] = -2.75994; qual_mismatch_simple_bayesian[14][9] = -2.93733; qual_mismatch_simple_bayesian[14][10] = -3.10478; qual_mismatch_simple_bayesian[14][11] = -3.2612; qual_mismatch_simple_bayesian[14][12] = -3.40563; qual_mismatch_simple_bayesian[14][13] = -3.53737; qual_mismatch_simple_bayesian[14][14] = -3.65598; qual_mismatch_simple_bayesian[14][15] = -3.76138; qual_mismatch_simple_bayesian[14][16] = -3.85381; qual_mismatch_simple_bayesian[14][17] = -3.93386; qual_mismatch_simple_bayesian[14][18] = -4.00234; qual_mismatch_simple_bayesian[14][19] = -4.0603; qual_mismatch_simple_bayesian[14][20] = -4.10885; qual_mismatch_simple_bayesian[14][21] = -4.14917; qual_mismatch_simple_bayesian[14][22] = -4.1824; qual_mismatch_simple_bayesian[14][23] = -4.20961; qual_mismatch_simple_bayesian[14][24] = -4.23176; qual_mismatch_simple_bayesian[14][25] = -4.24971; qual_mismatch_simple_bayesian[14][26] = -4.2642; qual_mismatch_simple_bayesian[14][27] = -4.27586; qual_mismatch_simple_bayesian[14][28] = -4.28523; qual_mismatch_simple_bayesian[14][29] = -4.29273; qual_mismatch_simple_bayesian[14][30] = -4.29872; qual_mismatch_simple_bayesian[14][31] = -4.30351; qual_mismatch_simple_bayesian[14][32] = -4.30734; qual_mismatch_simple_bayesian[14][33] = -4.31038; qual_mismatch_simple_bayesian[14][34] = -4.31281; qual_mismatch_simple_bayesian[14][35] = -4.31474; qual_mismatch_simple_bayesian[14][36] = -4.31627; qual_mismatch_simple_bayesian[14][37] = -4.3175; qual_mismatch_simple_bayesian[14][38] = -4.31847; qual_mismatch_simple_bayesian[14][39] = -4.31924; qual_mismatch_simple_bayesian[14][40] = -4.31986; qual_mismatch_simple_bayesian[14][41] = -4.32034; qual_mismatch_simple_bayesian[14][42] = -4.32073; qual_mismatch_simple_bayesian[14][43] = -4.32104; qual_mismatch_simple_bayesian[14][44] = -4.32128; qual_mismatch_simple_bayesian[14][45] = -4.32148; qual_mismatch_simple_bayesian[14][46] = -4.32163;
        qual_mismatch_simple_bayesian[15][0] = -1.10921; qual_mismatch_simple_bayesian[15][1] = -1.33123; qual_mismatch_simple_bayesian[15][2] = -1.55121; qual_mismatch_simple_bayesian[15][3] = -1.76867; qual_mismatch_simple_bayesian[15][4] = -1.98305; qual_mismatch_simple_bayesian[15][5] = -2.19368; qual_mismatch_simple_bayesian[15][6] = -2.39976; qual_mismatch_simple_bayesian[15][7] = -2.60038; qual_mismatch_simple_bayesian[15][8] = -2.79454; qual_mismatch_simple_bayesian[15][9] = -2.98112; qual_mismatch_simple_bayesian[15][10] = -3.15899; qual_mismatch_simple_bayesian[15][11] = -3.32696; qual_mismatch_simple_bayesian[15][12] = -3.48395; qual_mismatch_simple_bayesian[15][13] = -3.629; qual_mismatch_simple_bayesian[15][14] = -3.76138; qual_mismatch_simple_bayesian[15][15] = -3.88065; qual_mismatch_simple_bayesian[15][16] = -3.9867; qual_mismatch_simple_bayesian[15][17] = -4.07977; qual_mismatch_simple_bayesian[15][18] = -4.16041; qual_mismatch_simple_bayesian[15][19] = -4.22945; qual_mismatch_simple_bayesian[15][20] = -4.2879; qual_mismatch_simple_bayesian[15][21] = -4.3369; qual_mismatch_simple_bayesian[15][22] = -4.3776; qual_mismatch_simple_bayesian[15][23] = -4.41116; qual_mismatch_simple_bayesian[15][24] = -4.43864; qual_mismatch_simple_bayesian[15][25] = -4.46102; qual_mismatch_simple_bayesian[15][26] = -4.47916; qual_mismatch_simple_bayesian[15][27] = -4.49381; qual_mismatch_simple_bayesian[15][28] = -4.5056; qual_mismatch_simple_bayesian[15][29] = -4.51507; qual_mismatch_simple_bayesian[15][30] = -4.52265; qual_mismatch_simple_bayesian[15][31] = -4.52872; qual_mismatch_simple_bayesian[15][32] = -4.53356; qual_mismatch_simple_bayesian[15][33] = -4.53742; qual_mismatch_simple_bayesian[15][34] = -4.5405; qual_mismatch_simple_bayesian[15][35] = -4.54296; qual_mismatch_simple_bayesian[15][36] = -4.54491; qual_mismatch_simple_bayesian[15][37] = -4.54646; qual_mismatch_simple_bayesian[15][38] = -4.5477; qual_mismatch_simple_bayesian[15][39] = -4.54868; qual_mismatch_simple_bayesian[15][40] = -4.54947; qual_mismatch_simple_bayesian[15][41] = -4.55009; qual_mismatch_simple_bayesian[15][42] = -4.55058; qual_mismatch_simple_bayesian[15][43] = -4.55097; qual_mismatch_simple_bayesian[15][44] = -4.55128; qual_mismatch_simple_bayesian[15][45] = -4.55153; qual_mismatch_simple_bayesian[15][46] = -4.55173;
        qual_mismatch_simple_bayesian[16][0] = -1.10702; qual_mismatch_simple_bayesian[16][1] = -1.33074; qual_mismatch_simple_bayesian[16][2] = -1.55283; qual_mismatch_simple_bayesian[16][3] = -1.7729; qual_mismatch_simple_bayesian[16][4] = -1.99047; qual_mismatch_simple_bayesian[16][5] = -2.20499; qual_mismatch_simple_bayesian[16][6] = -2.41577; qual_mismatch_simple_bayesian[16][7] = -2.62204; qual_mismatch_simple_bayesian[16][8] = -2.8229; qual_mismatch_simple_bayesian[16][9] = -3.01733; qual_mismatch_simple_bayesian[16][10] = -3.20424; qual_mismatch_simple_bayesian[16][11] = -3.38246; qual_mismatch_simple_bayesian[16][12] = -3.55084; qual_mismatch_simple_bayesian[16][13] = -3.70828; qual_mismatch_simple_bayesian[16][14] = -3.85381; qual_mismatch_simple_bayesian[16][15] = -3.9867; qual_mismatch_simple_bayesian[16][16] = -4.10649; qual_mismatch_simple_bayesian[16][17] = -4.21306; qual_mismatch_simple_bayesian[16][18] = -4.30662; qual_mismatch_simple_bayesian[16][19] = -4.38774; qual_mismatch_simple_bayesian[16][20] = -4.45721; qual_mismatch_simple_bayesian[16][21] = -4.51606; qual_mismatch_simple_bayesian[16][22] = -4.5654; qual_mismatch_simple_bayesian[16][23] = -4.60641; qual_mismatch_simple_bayesian[16][24] = -4.64022; qual_mismatch_simple_bayesian[16][25] = -4.66792; qual_mismatch_simple_bayesian[16][26] = -4.69049; qual_mismatch_simple_bayesian[16][27] = -4.70878; qual_mismatch_simple_bayesian[16][28] = -4.72355; qual_mismatch_simple_bayesian[16][29] = -4.73544; qual_mismatch_simple_bayesian[16][30] = -4.74499; qual_mismatch_simple_bayesian[16][31] = -4.75264; qual_mismatch_simple_bayesian[16][32] = -4.75876; qual_mismatch_simple_bayesian[16][33] = -4.76365; qual_mismatch_simple_bayesian[16][34] = -4.76755; qual_mismatch_simple_bayesian[16][35] = -4.77065; qual_mismatch_simple_bayesian[16][36] = -4.77313; qual_mismatch_simple_bayesian[16][37] = -4.7751; qual_mismatch_simple_bayesian[16][38] = -4.77667; qual_mismatch_simple_bayesian[16][39] = -4.77792; qual_mismatch_simple_bayesian[16][40] = -4.77891; qual_mismatch_simple_bayesian[16][41] = -4.7797; qual_mismatch_simple_bayesian[16][42] = -4.78032; qual_mismatch_simple_bayesian[16][43] = -4.78082; qual_mismatch_simple_bayesian[16][44] = -4.78122; qual_mismatch_simple_bayesian[16][45] = -4.78153; qual_mismatch_simple_bayesian[16][46] = -4.78178;
        qual_mismatch_simple_bayesian[17][0] = -1.10529; qual_mismatch_simple_bayesian[17][1] = -1.33036; qual_mismatch_simple_bayesian[17][2] = -1.55412; qual_mismatch_simple_bayesian[17][3] = -1.77627; qual_mismatch_simple_bayesian[17][4] = -1.9964; qual_mismatch_simple_bayesian[17][5] = -2.21406; qual_mismatch_simple_bayesian[17][6] = -2.42868; qual_mismatch_simple_bayesian[17][7] = -2.63959; qual_mismatch_simple_bayesian[17][8] = -2.84602; qual_mismatch_simple_bayesian[17][9] = -3.04705; qual_mismatch_simple_bayesian[17][10] = -3.2417; qual_mismatch_simple_bayesian[17][11] = -3.42885; qual_mismatch_simple_bayesian[17][12] = -3.60736; qual_mismatch_simple_bayesian[17][13] = -3.77607; qual_mismatch_simple_bayesian[17][14] = -3.93386; qual_mismatch_simple_bayesian[17][15] = -4.07977; qual_mismatch_simple_bayesian[17][16] = -4.21306; qual_mismatch_simple_bayesian[17][17] = -4.33325; qual_mismatch_simple_bayesian[17][18] = -4.44022; qual_mismatch_simple_bayesian[17][19] = -4.53419; qual_mismatch_simple_bayesian[17][20] = -4.61567; qual_mismatch_simple_bayesian[17][21] = -4.68549; qual_mismatch_simple_bayesian[17][22] = -4.74465; qual_mismatch_simple_bayesian[17][23] = -4.79427; qual_mismatch_simple_bayesian[17][24] = -4.83552; qual_mismatch_simple_bayesian[17][25] = -4.86954; qual_mismatch_simple_bayesian[17][26] = -4.89741; qual_mismatch_simple_bayesian[17][27] = -4.92012; qual_mismatch_simple_bayesian[17][28] = -4.93853; qual_mismatch_simple_bayesian[17][29] = -4.9534; qual_mismatch_simple_bayesian[17][30] = -4.96537; qual_mismatch_simple_bayesian[17][31] = -4.97499; qual_mismatch_simple_bayesian[17][32] = -4.98269; qual_mismatch_simple_bayesian[17][33] = -4.98885; qual_mismatch_simple_bayesian[17][34] = -4.99377; qual_mismatch_simple_bayesian[17][35] = -4.9977; qual_mismatch_simple_bayesian[17][36] = -5.00083; qual_mismatch_simple_bayesian[17][37] = -5.00332; qual_mismatch_simple_bayesian[17][38] = -5.0053; qual_mismatch_simple_bayesian[17][39] = -5.00688; qual_mismatch_simple_bayesian[17][40] = -5.00814; qual_mismatch_simple_bayesian[17][41] = -5.00914; qual_mismatch_simple_bayesian[17][42] = -5.00993; qual_mismatch_simple_bayesian[17][43] = -5.01056; qual_mismatch_simple_bayesian[17][44] = -5.01107; qual_mismatch_simple_bayesian[17][45] = -5.01147; qual_mismatch_simple_bayesian[17][46] = -5.01178;
        qual_mismatch_simple_bayesian[18][0] = -1.10391; qual_mismatch_simple_bayesian[18][1] = -1.33005; qual_mismatch_simple_bayesian[18][2] = -1.55515; qual_mismatch_simple_bayesian[18][3] = -1.77895; qual_mismatch_simple_bayesian[18][4] = -2.00114; qual_mismatch_simple_bayesian[18][5] = -2.22133; qual_mismatch_simple_bayesian[18][6] = -2.43906; qual_mismatch_simple_bayesian[18][7] = -2.65376; qual_mismatch_simple_bayesian[18][8] = -2.86477; qual_mismatch_simple_bayesian[18][9] = -3.07131; qual_mismatch_simple_bayesian[18][10] = -3.27249; qual_mismatch_simple_bayesian[18][11] = -3.4673; qual_mismatch_simple_bayesian[18][12] = -3.65465; qual_mismatch_simple_bayesian[18][13] = -3.83339; qual_mismatch_simple_bayesian[18][14] = -4.00234; qual_mismatch_simple_bayesian[18][15] = -4.16041; qual_mismatch_simple_bayesian[18][16] = -4.30662; qual_mismatch_simple_bayesian[18][17] = -4.44022; qual_mismatch_simple_bayesian[18][18] = -4.56074; qual_mismatch_simple_bayesian[18][19] = -4.66803; qual_mismatch_simple_bayesian[18][20] = -4.76231; qual_mismatch_simple_bayesian[18][21] = -4.84409; qual_mismatch_simple_bayesian[18][22] = -4.91418; qual_mismatch_simple_bayesian[18][23] = -4.97359; qual_mismatch_simple_bayesian[18][24] = -5.02342; qual_mismatch_simple_bayesian[18][25] = -5.06486; qual_mismatch_simple_bayesian[18][26] = -5.09904; qual_mismatch_simple_bayesian[18][27] = -5.12706; qual_mismatch_simple_bayesian[18][28] = -5.14988; qual_mismatch_simple_bayesian[18][29] = -5.16839; qual_mismatch_simple_bayesian[18][30] = -5.18334; qual_mismatch_simple_bayesian[18][31] = -5.19537; qual_mismatch_simple_bayesian[18][32] = -5.20504; qual_mismatch_simple_bayesian[18][33] = -5.21278; qual_mismatch_simple_bayesian[18][34] = -5.21897; qual_mismatch_simple_bayesian[18][35] = -5.22392; qual_mismatch_simple_bayesian[18][36] = -5.22787; qual_mismatch_simple_bayesian[18][37] = -5.23102; qual_mismatch_simple_bayesian[18][38] = -5.23352; qual_mismatch_simple_bayesian[18][39] = -5.23552; qual_mismatch_simple_bayesian[18][40] = -5.23711; qual_mismatch_simple_bayesian[18][41] = -5.23837; qual_mismatch_simple_bayesian[18][42] = -5.23938; qual_mismatch_simple_bayesian[18][43] = -5.24017; qual_mismatch_simple_bayesian[18][44] = -5.24081; qual_mismatch_simple_bayesian[18][45] = -5.24131; qual_mismatch_simple_bayesian[18][46] = -5.24172;
        qual_mismatch_simple_bayesian[19][0] = -1.10282; qual_mismatch_simple_bayesian[19][1] = -1.32981; qual_mismatch_simple_bayesian[19][2] = -1.55597; qual_mismatch_simple_bayesian[19][3] = -1.78109; qual_mismatch_simple_bayesian[19][4] = -2.00492; qual_mismatch_simple_bayesian[19][5] = -2.22714; qual_mismatch_simple_bayesian[19][6] = -2.44737; qual_mismatch_simple_bayesian[19][7] = -2.66515; qual_mismatch_simple_bayesian[19][8] = -2.87992; qual_mismatch_simple_bayesian[19][9] = -3.09101; qual_mismatch_simple_bayesian[19][10] = -3.29764; qual_mismatch_simple_bayesian[19][11] = -3.49893; qual_mismatch_simple_bayesian[19][12] = -3.69388; qual_mismatch_simple_bayesian[19][13] = -3.88139; qual_mismatch_simple_bayesian[19][14] = -4.0603; qual_mismatch_simple_bayesian[19][15] = -4.22945; qual_mismatch_simple_bayesian[19][16] = -4.38774; qual_mismatch_simple_bayesian[19][17] = -4.53419; qual_mismatch_simple_bayesian[19][18] = -4.66803; qual_mismatch_simple_bayesian[19][19] = -4.78881; qual_mismatch_simple_bayesian[19][20] = -4.89635; qual_mismatch_simple_bayesian[19][21] = -4.99087; qual_mismatch_simple_bayesian[19][22] = -5.07289; qual_mismatch_simple_bayesian[19][23] = -5.1432; qual_mismatch_simple_bayesian[19][24] = -5.2028; qual_mismatch_simple_bayesian[19][25] = -5.25281; qual_mismatch_simple_bayesian[19][26] = -5.29439; qual_mismatch_simple_bayesian[19][27] = -5.32871; qual_mismatch_simple_bayesian[19][28] = -5.35683; qual_mismatch_simple_bayesian[19][29] = -5.37974; qual_mismatch_simple_bayesian[19][30] = -5.39832; qual_mismatch_simple_bayesian[19][31] = -5.41334; qual_mismatch_simple_bayesian[19][32] = -5.42542; qual_mismatch_simple_bayesian[19][33] = -5.43513; qual_mismatch_simple_bayesian[19][34] = -5.44291; qual_mismatch_simple_bayesian[19][35] = -5.44913; qual_mismatch_simple_bayesian[19][36] = -5.4541; qual_mismatch_simple_bayesian[19][37] = -5.45806; qual_mismatch_simple_bayesian[19][38] = -5.46122; qual_mismatch_simple_bayesian[19][39] = -5.46374; qual_mismatch_simple_bayesian[19][40] = -5.46574; qual_mismatch_simple_bayesian[19][41] = -5.46734; qual_mismatch_simple_bayesian[19][42] = -5.46861; qual_mismatch_simple_bayesian[19][43] = -5.46962; qual_mismatch_simple_bayesian[19][44] = -5.47042; qual_mismatch_simple_bayesian[19][45] = -5.47106; qual_mismatch_simple_bayesian[19][46] = -5.47156;
        qual_mismatch_simple_bayesian[20][0] = -1.10195; qual_mismatch_simple_bayesian[20][1] = -1.32962; qual_mismatch_simple_bayesian[20][2] = -1.55662; qual_mismatch_simple_bayesian[20][3] = -1.78279; qual_mismatch_simple_bayesian[20][4] = -2.00793; qual_mismatch_simple_bayesian[20][5] = -2.23178; qual_mismatch_simple_bayesian[20][6] = -2.45403; qual_mismatch_simple_bayesian[20][7] = -2.6743; qual_mismatch_simple_bayesian[20][8] = -2.89212; qual_mismatch_simple_bayesian[20][9] = -3.10693; qual_mismatch_simple_bayesian[20][10] = -3.31808; qual_mismatch_simple_bayesian[20][11] = -3.52479; qual_mismatch_simple_bayesian[20][12] = -3.72617; qual_mismatch_simple_bayesian[20][13] = -3.92122; qual_mismatch_simple_bayesian[20][14] = -4.10885; qual_mismatch_simple_bayesian[20][15] = -4.2879; qual_mismatch_simple_bayesian[20][16] = -4.45721; qual_mismatch_simple_bayesian[20][17] = -4.61567; qual_mismatch_simple_bayesian[20][18] = -4.76231; qual_mismatch_simple_bayesian[20][19] = -4.89635; qual_mismatch_simple_bayesian[20][20] = -5.01732; qual_mismatch_simple_bayesian[20][21] = -5.12507; qual_mismatch_simple_bayesian[20][22] = -5.21979; qual_mismatch_simple_bayesian[20][23] = -5.30199; qual_mismatch_simple_bayesian[20][24] = -5.37247; qual_mismatch_simple_bayesian[20][25] = -5.43222; qual_mismatch_simple_bayesian[20][26] = -5.48237; qual_mismatch_simple_bayesian[20][27] = -5.52408; qual_mismatch_simple_bayesian[20][28] = -5.55849; qual_mismatch_simple_bayesian[20][29] = -5.5867; qual_mismatch_simple_bayesian[20][30] = -5.60969; qual_mismatch_simple_bayesian[20][31] = -5.62833; qual_mismatch_simple_bayesian[20][32] = -5.64339; qual_mismatch_simple_bayesian[20][33] = -5.65552; qual_mismatch_simple_bayesian[20][34] = -5.66525; qual_mismatch_simple_bayesian[20][35] = -5.67306; qual_mismatch_simple_bayesian[20][36] = -5.6793; qual_mismatch_simple_bayesian[20][37] = -5.68429; qual_mismatch_simple_bayesian[20][38] = -5.68827; qual_mismatch_simple_bayesian[20][39] = -5.69144; qual_mismatch_simple_bayesian[20][40] = -5.69396; qual_mismatch_simple_bayesian[20][41] = -5.69598; qual_mismatch_simple_bayesian[20][42] = -5.69758; qual_mismatch_simple_bayesian[20][43] = -5.69885; qual_mismatch_simple_bayesian[20][44] = -5.69986; qual_mismatch_simple_bayesian[20][45] = -5.70067; qual_mismatch_simple_bayesian[20][46] = -5.70131;
        qual_mismatch_simple_bayesian[21][0] = -1.10126; qual_mismatch_simple_bayesian[21][1] = -1.32946; qual_mismatch_simple_bayesian[21][2] = -1.55713; qual_mismatch_simple_bayesian[21][3] = -1.78414; qual_mismatch_simple_bayesian[21][4] = -2.01033; qual_mismatch_simple_bayesian[21][5] = -2.23548; qual_mismatch_simple_bayesian[21][6] = -2.45935; qual_mismatch_simple_bayesian[21][7] = -2.68162; qual_mismatch_simple_bayesian[21][8] = -2.90191; qual_mismatch_simple_bayesian[21][9] = -3.11977; qual_mismatch_simple_bayesian[21][10] = -3.33462; qual_mismatch_simple_bayesian[21][11] = -3.54582; qual_mismatch_simple_bayesian[21][12] = -3.75259; qual_mismatch_simple_bayesian[21][13] = -3.95404; qual_mismatch_simple_bayesian[21][14] = -4.14917; qual_mismatch_simple_bayesian[21][15] = -4.3369; qual_mismatch_simple_bayesian[21][16] = -4.51606; qual_mismatch_simple_bayesian[21][17] = -4.68549; qual_mismatch_simple_bayesian[21][18] = -4.84409; qual_mismatch_simple_bayesian[21][19] = -4.99087; qual_mismatch_simple_bayesian[21][20] = -5.12507; qual_mismatch_simple_bayesian[21][21] = -5.2462; qual_mismatch_simple_bayesian[21][22] = -5.35411; qual_mismatch_simple_bayesian[21][23] = -5.44898; qual_mismatch_simple_bayesian[21][24] = -5.53133; qual_mismatch_simple_bayesian[21][25] = -5.60194; qual_mismatch_simple_bayesian[21][26] = -5.66182; qual_mismatch_simple_bayesian[21][27] = -5.71208; qual_mismatch_simple_bayesian[21][28] = -5.75388; qual_mismatch_simple_bayesian[21][29] = -5.78837; qual_mismatch_simple_bayesian[21][30] = -5.81665; qual_mismatch_simple_bayesian[21][31] = -5.83969; qual_mismatch_simple_bayesian[21][32] = -5.85838; qual_mismatch_simple_bayesian[21][33] = -5.87348; qual_mismatch_simple_bayesian[21][34] = -5.88564; qual_mismatch_simple_bayesian[21][35] = -5.89541; qual_mismatch_simple_bayesian[21][36] = -5.90323; qual_mismatch_simple_bayesian[21][37] = -5.90949; qual_mismatch_simple_bayesian[21][38] = -5.91449; qual_mismatch_simple_bayesian[21][39] = -5.91848; qual_mismatch_simple_bayesian[21][40] = -5.92166; qual_mismatch_simple_bayesian[21][41] = -5.9242; qual_mismatch_simple_bayesian[21][42] = -5.92621; qual_mismatch_simple_bayesian[21][43] = -5.92782; qual_mismatch_simple_bayesian[21][44] = -5.92909; qual_mismatch_simple_bayesian[21][45] = -5.93011; qual_mismatch_simple_bayesian[21][46] = -5.93092;
        qual_mismatch_simple_bayesian[22][0] = -1.10072; qual_mismatch_simple_bayesian[22][1] = -1.32934; qual_mismatch_simple_bayesian[22][2] = -1.55754; qual_mismatch_simple_bayesian[22][3] = -1.78522; qual_mismatch_simple_bayesian[22][4] = -2.01224; qual_mismatch_simple_bayesian[22][5] = -2.23843; qual_mismatch_simple_bayesian[22][6] = -2.4636; qual_mismatch_simple_bayesian[22][7] = -2.68748; qual_mismatch_simple_bayesian[22][8] = -2.90977; qual_mismatch_simple_bayesian[22][9] = -3.13008; qual_mismatch_simple_bayesian[22][10] = -3.34796; qual_mismatch_simple_bayesian[22][11] = -3.56284; qual_mismatch_simple_bayesian[22][12] = -3.77408; qual_mismatch_simple_bayesian[22][13] = -3.9809; qual_mismatch_simple_bayesian[22][14] = -4.1824; qual_mismatch_simple_bayesian[22][15] = -4.3776; qual_mismatch_simple_bayesian[22][16] = -4.5654; qual_mismatch_simple_bayesian[22][17] = -4.74465; qual_mismatch_simple_bayesian[22][18] = -4.91418; qual_mismatch_simple_bayesian[22][19] = -5.07289; qual_mismatch_simple_bayesian[22][20] = -5.21979; qual_mismatch_simple_bayesian[22][21] = -5.35411; qual_mismatch_simple_bayesian[22][22] = -5.47537; qual_mismatch_simple_bayesian[22][23] = -5.5834; qual_mismatch_simple_bayesian[22][24] = -5.67839; qual_mismatch_simple_bayesian[22][25] = -5.76086; qual_mismatch_simple_bayesian[22][26] = -5.83158; qual_mismatch_simple_bayesian[22][27] = -5.89155; qual_mismatch_simple_bayesian[22][28] = -5.9419; qual_mismatch_simple_bayesian[22][29] = -5.98377; qual_mismatch_simple_bayesian[22][30] = -6.01833; qual_mismatch_simple_bayesian[22][31] = -6.04666; qual_mismatch_simple_bayesian[22][32] = -6.06975; qual_mismatch_simple_bayesian[22][33] = -6.08848; qual_mismatch_simple_bayesian[22][34] = -6.10361; qual_mismatch_simple_bayesian[22][35] = -6.1158; qual_mismatch_simple_bayesian[22][36] = -6.12558; qual_mismatch_simple_bayesian[22][37] = -6.13342; qual_mismatch_simple_bayesian[22][38] = -6.1397; qual_mismatch_simple_bayesian[22][39] = -6.14471; qual_mismatch_simple_bayesian[22][40] = -6.14871; qual_mismatch_simple_bayesian[22][41] = -6.15189; qual_mismatch_simple_bayesian[22][42] = -6.15443; qual_mismatch_simple_bayesian[22][43] = -6.15645; qual_mismatch_simple_bayesian[22][44] = -6.15806; qual_mismatch_simple_bayesian[22][45] = -6.15934; qual_mismatch_simple_bayesian[22][46] = -6.16036;
        qual_mismatch_simple_bayesian[23][0] = -1.10028; qual_mismatch_simple_bayesian[23][1] = -1.32924; qual_mismatch_simple_bayesian[23][2] = -1.55787; qual_mismatch_simple_bayesian[23][3] = -1.78608; qual_mismatch_simple_bayesian[23][4] = -2.01376; qual_mismatch_simple_bayesian[23][5] = -2.24078; qual_mismatch_simple_bayesian[23][6] = -2.46698; qual_mismatch_simple_bayesian[23][7] = -2.69215; qual_mismatch_simple_bayesian[23][8] = -2.91605; qual_mismatch_simple_bayesian[23][9] = -3.13835; qual_mismatch_simple_bayesian[23][10] = -3.35868; qual_mismatch_simple_bayesian[23][11] = -3.57658; qual_mismatch_simple_bayesian[23][12] = -3.79149; qual_mismatch_simple_bayesian[23][13] = -4.00276; qual_mismatch_simple_bayesian[23][14] = -4.20961; qual_mismatch_simple_bayesian[23][15] = -4.41116; qual_mismatch_simple_bayesian[23][16] = -4.60641; qual_mismatch_simple_bayesian[23][17] = -4.79427; qual_mismatch_simple_bayesian[23][18] = -4.97359; qual_mismatch_simple_bayesian[23][19] = -5.1432; qual_mismatch_simple_bayesian[23][20] = -5.30199; qual_mismatch_simple_bayesian[23][21] = -5.44898; qual_mismatch_simple_bayesian[23][22] = -5.5834; qual_mismatch_simple_bayesian[23][23] = -5.70476; qual_mismatch_simple_bayesian[23][24] = -5.81289; qual_mismatch_simple_bayesian[23][25] = -5.90798; qual_mismatch_simple_bayesian[23][26] = -5.99054; qual_mismatch_simple_bayesian[23][27] = -6.06134; qual_mismatch_simple_bayesian[23][28] = -6.12139; qual_mismatch_simple_bayesian[23][29] = -6.17181; qual_mismatch_simple_bayesian[23][30] = -6.21374; qual_mismatch_simple_bayesian[23][31] = -6.24836; qual_mismatch_simple_bayesian[23][32] = -6.27673; qual_mismatch_simple_bayesian[23][33] = -6.29986; qual_mismatch_simple_bayesian[23][34] = -6.31861; qual_mismatch_simple_bayesian[23][35] = -6.33377; qual_mismatch_simple_bayesian[23][36] = -6.34597; qual_mismatch_simple_bayesian[23][37] = -6.35578; qual_mismatch_simple_bayesian[23][38] = -6.36363; qual_mismatch_simple_bayesian[23][39] = -6.36991; qual_mismatch_simple_bayesian[23][40] = -6.37493; qual_mismatch_simple_bayesian[23][41] = -6.37894; qual_mismatch_simple_bayesian[23][42] = -6.38213; qual_mismatch_simple_bayesian[23][43] = -6.38467; qual_mismatch_simple_bayesian[23][44] = -6.3867; qual_mismatch_simple_bayesian[23][45] = -6.38831; qual_mismatch_simple_bayesian[23][46] = -6.38959;
        qual_mismatch_simple_bayesian[24][0] = -1.09994; qual_mismatch_simple_bayesian[24][1] = -1.32917; qual_mismatch_simple_bayesian[24][2] = -1.55813; qual_mismatch_simple_bayesian[24][3] = -1.78676; qual_mismatch_simple_bayesian[24][4] = -2.01497; qual_mismatch_simple_bayesian[24][5] = -2.24265; qual_mismatch_simple_bayesian[24][6] = -2.46968; qual_mismatch_simple_bayesian[24][7] = -2.69588; qual_mismatch_simple_bayesian[24][8] = -2.92106; qual_mismatch_simple_bayesian[24][9] = -3.14496; qual_mismatch_simple_bayesian[24][10] = -3.36728; qual_mismatch_simple_bayesian[24][11] = -3.58762; qual_mismatch_simple_bayesian[24][12] = -3.80553; qual_mismatch_simple_bayesian[24][13] = -4.02047; qual_mismatch_simple_bayesian[24][14] = -4.23176; qual_mismatch_simple_bayesian[24][15] = -4.43864; qual_mismatch_simple_bayesian[24][16] = -4.64022; qual_mismatch_simple_bayesian[24][17] = -4.83552; qual_mismatch_simple_bayesian[24][18] = -5.02342; qual_mismatch_simple_bayesian[24][19] = -5.2028; qual_mismatch_simple_bayesian[24][20] = -5.37247; qual_mismatch_simple_bayesian[24][21] = -5.53133; qual_mismatch_simple_bayesian[24][22] = -5.67839; qual_mismatch_simple_bayesian[24][23] = -5.81289; qual_mismatch_simple_bayesian[24][24] = -5.93433; qual_mismatch_simple_bayesian[24][25] = -6.04254; qual_mismatch_simple_bayesian[24][26] = -6.1377; qual_mismatch_simple_bayesian[24][27] = -6.22033; qual_mismatch_simple_bayesian[24][28] = -6.29121; qual_mismatch_simple_bayesian[24][29] = -6.35132; qual_mismatch_simple_bayesian[24][30] = -6.40179; qual_mismatch_simple_bayesian[24][31] = -6.44377; qual_mismatch_simple_bayesian[24][32] = -6.47843; qual_mismatch_simple_bayesian[24][33] = -6.50683; qual_mismatch_simple_bayesian[24][34] = -6.52999; qual_mismatch_simple_bayesian[24][35] = -6.54877; qual_mismatch_simple_bayesian[24][36] = -6.56395; qual_mismatch_simple_bayesian[24][37] = -6.57617; qual_mismatch_simple_bayesian[24][38] = -6.58598; qual_mismatch_simple_bayesian[24][39] = -6.59385; qual_mismatch_simple_bayesian[24][40] = -6.60014; qual_mismatch_simple_bayesian[24][41] = -6.60516; qual_mismatch_simple_bayesian[24][42] = -6.60917; qual_mismatch_simple_bayesian[24][43] = -6.61237; qual_mismatch_simple_bayesian[24][44] = -6.61492; qual_mismatch_simple_bayesian[24][45] = -6.61695; qual_mismatch_simple_bayesian[24][46] = -6.61856;
        qual_mismatch_simple_bayesian[25][0] = -1.09967; qual_mismatch_simple_bayesian[25][1] = -1.32911; qual_mismatch_simple_bayesian[25][2] = -1.55833; qual_mismatch_simple_bayesian[25][3] = -1.7873; qual_mismatch_simple_bayesian[25][4] = -2.01593; qual_mismatch_simple_bayesian[25][5] = -2.24414; qual_mismatch_simple_bayesian[25][6] = -2.47183; qual_mismatch_simple_bayesian[25][7] = -2.69886; qual_mismatch_simple_bayesian[25][8] = -2.92507; qual_mismatch_simple_bayesian[25][9] = -3.15025; qual_mismatch_simple_bayesian[25][10] = -3.37416; qual_mismatch_simple_bayesian[25][11] = -3.59648; qual_mismatch_simple_bayesian[25][12] = -3.81683; qual_mismatch_simple_bayesian[25][13] = -4.03476; qual_mismatch_simple_bayesian[25][14] = -4.24971; qual_mismatch_simple_bayesian[25][15] = -4.46102; qual_mismatch_simple_bayesian[25][16] = -4.66792; qual_mismatch_simple_bayesian[25][17] = -4.86954; qual_mismatch_simple_bayesian[25][18] = -5.06486; qual_mismatch_simple_bayesian[25][19] = -5.25281; qual_mismatch_simple_bayesian[25][20] = -5.43222; qual_mismatch_simple_bayesian[25][21] = -5.60194; qual_mismatch_simple_bayesian[25][22] = -5.76086; qual_mismatch_simple_bayesian[25][23] = -5.90798; qual_mismatch_simple_bayesian[25][24] = -6.04254; qual_mismatch_simple_bayesian[25][25] = -6.16404; qual_mismatch_simple_bayesian[25][26] = -6.27231; qual_mismatch_simple_bayesian[25][27] = -6.36754; qual_mismatch_simple_bayesian[25][28] = -6.45023; qual_mismatch_simple_bayesian[25][29] = -6.52116; qual_mismatch_simple_bayesian[25][30] = -6.58132; qual_mismatch_simple_bayesian[25][31] = -6.63183; qual_mismatch_simple_bayesian[25][32] = -6.67385; qual_mismatch_simple_bayesian[25][33] = -6.70854; qual_mismatch_simple_bayesian[25][34] = -6.73697; qual_mismatch_simple_bayesian[25][35] = -6.76015; qual_mismatch_simple_bayesian[25][36] = -6.77895; qual_mismatch_simple_bayesian[25][37] = -6.79414; qual_mismatch_simple_bayesian[25][38] = -6.80637; qual_mismatch_simple_bayesian[25][39] = -6.8162; qual_mismatch_simple_bayesian[25][40] = -6.82407; qual_mismatch_simple_bayesian[25][41] = -6.83037; qual_mismatch_simple_bayesian[25][42] = -6.8354; qual_mismatch_simple_bayesian[25][43] = -6.83942; qual_mismatch_simple_bayesian[25][44] = -6.84262; qual_mismatch_simple_bayesian[25][45] = -6.84517; qual_mismatch_simple_bayesian[25][46] = -6.8472;
        qual_mismatch_simple_bayesian[26][0] = -1.09945; qual_mismatch_simple_bayesian[26][1] = -1.32906; qual_mismatch_simple_bayesian[26][2] = -1.5585; qual_mismatch_simple_bayesian[26][3] = -1.78773; qual_mismatch_simple_bayesian[26][4] = -2.01669; qual_mismatch_simple_bayesian[26][5] = -2.24532; qual_mismatch_simple_bayesian[26][6] = -2.47353; qual_mismatch_simple_bayesian[26][7] = -2.70122; qual_mismatch_simple_bayesian[26][8] = -2.92826; qual_mismatch_simple_bayesian[26][9] = -3.15447; qual_mismatch_simple_bayesian[26][10] = -3.37966; qual_mismatch_simple_bayesian[26][11] = -3.60357; qual_mismatch_simple_bayesian[26][12] = -3.8259; qual_mismatch_simple_bayesian[26][13] = -4.04626; qual_mismatch_simple_bayesian[26][14] = -4.2642; qual_mismatch_simple_bayesian[26][15] = -4.47916; qual_mismatch_simple_bayesian[26][16] = -4.69049; qual_mismatch_simple_bayesian[26][17] = -4.89741; qual_mismatch_simple_bayesian[26][18] = -5.09904; qual_mismatch_simple_bayesian[26][19] = -5.29439; qual_mismatch_simple_bayesian[26][20] = -5.48237; qual_mismatch_simple_bayesian[26][21] = -5.66182; qual_mismatch_simple_bayesian[26][22] = -5.83158; qual_mismatch_simple_bayesian[26][23] = -5.99054; qual_mismatch_simple_bayesian[26][24] = -6.1377; qual_mismatch_simple_bayesian[26][25] = -6.27231; qual_mismatch_simple_bayesian[26][26] = -6.39386; qual_mismatch_simple_bayesian[26][27] = -6.50219; qual_mismatch_simple_bayesian[26][28] = -6.59746; qual_mismatch_simple_bayesian[26][29] = -6.6802; qual_mismatch_simple_bayesian[26][30] = -6.75117; qual_mismatch_simple_bayesian[26][31] = -6.81137; qual_mismatch_simple_bayesian[26][32] = -6.86191; qual_mismatch_simple_bayesian[26][33] = -6.90396; qual_mismatch_simple_bayesian[26][34] = -6.93867; qual_mismatch_simple_bayesian[26][35] = -6.96713; qual_mismatch_simple_bayesian[26][36] = -6.99033; qual_mismatch_simple_bayesian[26][37] = -7.00914; qual_mismatch_simple_bayesian[26][38] = -7.02435; qual_mismatch_simple_bayesian[26][39] = -7.03659; qual_mismatch_simple_bayesian[26][40] = -7.04642; qual_mismatch_simple_bayesian[26][41] = -7.0543; qual_mismatch_simple_bayesian[26][42] = -7.06061; qual_mismatch_simple_bayesian[26][43] = -7.06564; qual_mismatch_simple_bayesian[26][44] = -7.06966; qual_mismatch_simple_bayesian[26][45] = -7.07286; qual_mismatch_simple_bayesian[26][46] = -7.07542;
        qual_mismatch_simple_bayesian[27][0] = -1.09928; qual_mismatch_simple_bayesian[27][1] = -1.32902; qual_mismatch_simple_bayesian[27][2] = -1.55863; qual_mismatch_simple_bayesian[27][3] = -1.78807; qual_mismatch_simple_bayesian[27][4] = -2.0173; qual_mismatch_simple_bayesian[27][5] = -2.24626; qual_mismatch_simple_bayesian[27][6] = -2.47489; qual_mismatch_simple_bayesian[27][7] = -2.70311; qual_mismatch_simple_bayesian[27][8] = -2.9308; qual_mismatch_simple_bayesian[27][9] = -3.15784; qual_mismatch_simple_bayesian[27][10] = -3.38405; qual_mismatch_simple_bayesian[27][11] = -3.60925; qual_mismatch_simple_bayesian[27][12] = -3.83316; qual_mismatch_simple_bayesian[27][13] = -4.0555; qual_mismatch_simple_bayesian[27][14] = -4.27586; qual_mismatch_simple_bayesian[27][15] = -4.49381; qual_mismatch_simple_bayesian[27][16] = -4.70878; qual_mismatch_simple_bayesian[27][17] = -4.92012; qual_mismatch_simple_bayesian[27][18] = -5.12706; qual_mismatch_simple_bayesian[27][19] = -5.32871; qual_mismatch_simple_bayesian[27][20] = -5.52408; qual_mismatch_simple_bayesian[27][21] = -5.71208; qual_mismatch_simple_bayesian[27][22] = -5.89155; qual_mismatch_simple_bayesian[27][23] = -6.06134; qual_mismatch_simple_bayesian[27][24] = -6.22033; qual_mismatch_simple_bayesian[27][25] = -6.36754; qual_mismatch_simple_bayesian[27][26] = -6.50219; qual_mismatch_simple_bayesian[27][27] = -6.62378; qual_mismatch_simple_bayesian[27][28] = -6.73214; qual_mismatch_simple_bayesian[27][29] = -6.82745; qual_mismatch_simple_bayesian[27][30] = -6.91022; qual_mismatch_simple_bayesian[27][31] = -6.98123; qual_mismatch_simple_bayesian[27][32] = -7.04146; qual_mismatch_simple_bayesian[27][33] = -7.09203; qual_mismatch_simple_bayesian[27][34] = -7.13411; qual_mismatch_simple_bayesian[27][35] = -7.16884; qual_mismatch_simple_bayesian[27][36] = -7.19731; qual_mismatch_simple_bayesian[27][37] = -7.22052; qual_mismatch_simple_bayesian[27][38] = -7.23935; qual_mismatch_simple_bayesian[27][39] = -7.25456; qual_mismatch_simple_bayesian[27][40] = -7.26682; qual_mismatch_simple_bayesian[27][41] = -7.27666; qual_mismatch_simple_bayesian[27][42] = -7.28454; qual_mismatch_simple_bayesian[27][43] = -7.29085; qual_mismatch_simple_bayesian[27][44] = -7.29589; qual_mismatch_simple_bayesian[27][45] = -7.29991; qual_mismatch_simple_bayesian[27][46] = -7.30311;
        qual_mismatch_simple_bayesian[28][0] = -1.09914; qual_mismatch_simple_bayesian[28][1] = -1.32899; qual_mismatch_simple_bayesian[28][2] = -1.55873; qual_mismatch_simple_bayesian[28][3] = -1.78834; qual_mismatch_simple_bayesian[28][4] = -2.01778; qual_mismatch_simple_bayesian[28][5] = -2.24701; qual_mismatch_simple_bayesian[28][6] = -2.47598; qual_mismatch_simple_bayesian[28][7] = -2.70461; qual_mismatch_simple_bayesian[28][8] = -2.93282; qual_mismatch_simple_bayesian[28][9] = -3.16052; qual_mismatch_simple_bayesian[28][10] = -3.38756; qual_mismatch_simple_bayesian[28][11] = -3.61377; qual_mismatch_simple_bayesian[28][12] = -3.83897; qual_mismatch_simple_bayesian[28][13] = -4.06289; qual_mismatch_simple_bayesian[28][14] = -4.28523; qual_mismatch_simple_bayesian[28][15] = -4.5056; qual_mismatch_simple_bayesian[28][16] = -4.72355; qual_mismatch_simple_bayesian[28][17] = -4.93853; qual_mismatch_simple_bayesian[28][18] = -5.14988; qual_mismatch_simple_bayesian[28][19] = -5.35683; qual_mismatch_simple_bayesian[28][20] = -5.55849; qual_mismatch_simple_bayesian[28][21] = -5.75388; qual_mismatch_simple_bayesian[28][22] = -5.9419; qual_mismatch_simple_bayesian[28][23] = -6.12139; qual_mismatch_simple_bayesian[28][24] = -6.29121; qual_mismatch_simple_bayesian[28][25] = -6.45023; qual_mismatch_simple_bayesian[28][26] = -6.59746; qual_mismatch_simple_bayesian[28][27] = -6.73214; qual_mismatch_simple_bayesian[28][28] = -6.85376; qual_mismatch_simple_bayesian[28][29] = -6.96216; qual_mismatch_simple_bayesian[28][30] = -7.0575; qual_mismatch_simple_bayesian[28][31] = -7.1403; qual_mismatch_simple_bayesian[28][32] = -7.21133; qual_mismatch_simple_bayesian[28][33] = -7.27159; qual_mismatch_simple_bayesian[28][34] = -7.32218; qual_mismatch_simple_bayesian[28][35] = -7.36428; qual_mismatch_simple_bayesian[28][36] = -7.39902; qual_mismatch_simple_bayesian[28][37] = -7.42751; qual_mismatch_simple_bayesian[28][38] = -7.45073; qual_mismatch_simple_bayesian[28][39] = -7.46957; qual_mismatch_simple_bayesian[28][40] = -7.48479; qual_mismatch_simple_bayesian[28][41] = -7.49705; qual_mismatch_simple_bayesian[28][42] = -7.50689; qual_mismatch_simple_bayesian[28][43] = -7.51478; qual_mismatch_simple_bayesian[28][44] = -7.52109; qual_mismatch_simple_bayesian[28][45] = -7.52614; qual_mismatch_simple_bayesian[28][46] = -7.53016;
        qual_mismatch_simple_bayesian[29][0] = -1.09903; qual_mismatch_simple_bayesian[29][1] = -1.32896; qual_mismatch_simple_bayesian[29][2] = -1.55881; qual_mismatch_simple_bayesian[29][3] = -1.78855; qual_mismatch_simple_bayesian[29][4] = -2.01816; qual_mismatch_simple_bayesian[29][5] = -2.2476; qual_mismatch_simple_bayesian[29][6] = -2.47684; qual_mismatch_simple_bayesian[29][7] = -2.7058; qual_mismatch_simple_bayesian[29][8] = -2.93444; qual_mismatch_simple_bayesian[29][9] = -3.16265; qual_mismatch_simple_bayesian[29][10] = -3.39035; qual_mismatch_simple_bayesian[29][11] = -3.61738; qual_mismatch_simple_bayesian[29][12] = -3.84361; qual_mismatch_simple_bayesian[29][13] = -4.0688; qual_mismatch_simple_bayesian[29][14] = -4.29273; qual_mismatch_simple_bayesian[29][15] = -4.51507; qual_mismatch_simple_bayesian[29][16] = -4.73544; qual_mismatch_simple_bayesian[29][17] = -4.9534; qual_mismatch_simple_bayesian[29][18] = -5.16839; qual_mismatch_simple_bayesian[29][19] = -5.37974; qual_mismatch_simple_bayesian[29][20] = -5.5867; qual_mismatch_simple_bayesian[29][21] = -5.78837; qual_mismatch_simple_bayesian[29][22] = -5.98377; qual_mismatch_simple_bayesian[29][23] = -6.17181; qual_mismatch_simple_bayesian[29][24] = -6.35132; qual_mismatch_simple_bayesian[29][25] = -6.52116; qual_mismatch_simple_bayesian[29][26] = -6.6802; qual_mismatch_simple_bayesian[29][27] = -6.82745; qual_mismatch_simple_bayesian[29][28] = -6.96216; qual_mismatch_simple_bayesian[29][29] = -7.0838; qual_mismatch_simple_bayesian[29][30] = -7.19222; qual_mismatch_simple_bayesian[29][31] = -7.28759; qual_mismatch_simple_bayesian[29][32] = -7.37041; qual_mismatch_simple_bayesian[29][33] = -7.44147; qual_mismatch_simple_bayesian[29][34] = -7.50174; qual_mismatch_simple_bayesian[29][35] = -7.55235; qual_mismatch_simple_bayesian[29][36] = -7.59446; qual_mismatch_simple_bayesian[29][37] = -7.62922; qual_mismatch_simple_bayesian[29][38] = -7.65772; qual_mismatch_simple_bayesian[29][39] = -7.68095; qual_mismatch_simple_bayesian[29][40] = -7.6998; qual_mismatch_simple_bayesian[29][41] = -7.71502; qual_mismatch_simple_bayesian[29][42] = -7.72729; qual_mismatch_simple_bayesian[29][43] = -7.73713; qual_mismatch_simple_bayesian[29][44] = -7.74503; qual_mismatch_simple_bayesian[29][45] = -7.75134; qual_mismatch_simple_bayesian[29][46] = -7.75639;
        qual_mismatch_simple_bayesian[30][0] = -1.09895; qual_mismatch_simple_bayesian[30][1] = -1.32895; qual_mismatch_simple_bayesian[30][2] = -1.55888; qual_mismatch_simple_bayesian[30][3] = -1.78873; qual_mismatch_simple_bayesian[30][4] = -2.01847; qual_mismatch_simple_bayesian[30][5] = -2.24808; qual_mismatch_simple_bayesian[30][6] = -2.47752; qual_mismatch_simple_bayesian[30][7] = -2.70675; qual_mismatch_simple_bayesian[30][8] = -2.93572; qual_mismatch_simple_bayesian[30][9] = -3.16435; qual_mismatch_simple_bayesian[30][10] = -3.39257; qual_mismatch_simple_bayesian[30][11] = -3.62026; qual_mismatch_simple_bayesian[30][12] = -3.8473; qual_mismatch_simple_bayesian[30][13] = -4.07352; qual_mismatch_simple_bayesian[30][14] = -4.29872; qual_mismatch_simple_bayesian[30][15] = -4.52265; qual_mismatch_simple_bayesian[30][16] = -4.74499; qual_mismatch_simple_bayesian[30][17] = -4.96537; qual_mismatch_simple_bayesian[30][18] = -5.18334; qual_mismatch_simple_bayesian[30][19] = -5.39832; qual_mismatch_simple_bayesian[30][20] = -5.60969; qual_mismatch_simple_bayesian[30][21] = -5.81665; qual_mismatch_simple_bayesian[30][22] = -6.01833; qual_mismatch_simple_bayesian[30][23] = -6.21374; qual_mismatch_simple_bayesian[30][24] = -6.40179; qual_mismatch_simple_bayesian[30][25] = -6.58132; qual_mismatch_simple_bayesian[30][26] = -6.75117; qual_mismatch_simple_bayesian[30][27] = -6.91022; qual_mismatch_simple_bayesian[30][28] = -7.0575; qual_mismatch_simple_bayesian[30][29] = -7.19222; qual_mismatch_simple_bayesian[30][30] = -7.31389; qual_mismatch_simple_bayesian[30][31] = -7.42233; qual_mismatch_simple_bayesian[30][32] = -7.51772; qual_mismatch_simple_bayesian[30][33] = -7.60056; qual_mismatch_simple_bayesian[30][34] = -7.67163; qual_mismatch_simple_bayesian[30][35] = -7.73192; qual_mismatch_simple_bayesian[30][36] = -7.78254; qual_mismatch_simple_bayesian[30][37] = -7.82466; qual_mismatch_simple_bayesian[30][38] = -7.85943; qual_mismatch_simple_bayesian[30][39] = -7.88794; qual_mismatch_simple_bayesian[30][40] = -7.91118; qual_mismatch_simple_bayesian[30][41] = -7.93003; qual_mismatch_simple_bayesian[30][42] = -7.94526; qual_mismatch_simple_bayesian[30][43] = -7.95753; qual_mismatch_simple_bayesian[30][44] = -7.96738; qual_mismatch_simple_bayesian[30][45] = -7.97528; qual_mismatch_simple_bayesian[30][46] = -7.98159;
        qual_mismatch_simple_bayesian[31][0] = -1.09888; qual_mismatch_simple_bayesian[31][1] = -1.32893; qual_mismatch_simple_bayesian[31][2] = -1.55893; qual_mismatch_simple_bayesian[31][3] = -1.78886; qual_mismatch_simple_bayesian[31][4] = -2.01871; qual_mismatch_simple_bayesian[31][5] = -2.24845; qual_mismatch_simple_bayesian[31][6] = -2.47806; qual_mismatch_simple_bayesian[31][7] = -2.7075; qual_mismatch_simple_bayesian[31][8] = -2.93674; qual_mismatch_simple_bayesian[31][9] = -3.1657; qual_mismatch_simple_bayesian[31][10] = -3.39434; qual_mismatch_simple_bayesian[31][11] = -3.62255; qual_mismatch_simple_bayesian[31][12] = -3.85025; qual_mismatch_simple_bayesian[31][13] = -4.07729; qual_mismatch_simple_bayesian[31][14] = -4.30351; qual_mismatch_simple_bayesian[31][15] = -4.52872; qual_mismatch_simple_bayesian[31][16] = -4.75264; qual_mismatch_simple_bayesian[31][17] = -4.97499; qual_mismatch_simple_bayesian[31][18] = -5.19537; qual_mismatch_simple_bayesian[31][19] = -5.41334; qual_mismatch_simple_bayesian[31][20] = -5.62833; qual_mismatch_simple_bayesian[31][21] = -5.83969; qual_mismatch_simple_bayesian[31][22] = -6.04666; qual_mismatch_simple_bayesian[31][23] = -6.24836; qual_mismatch_simple_bayesian[31][24] = -6.44377; qual_mismatch_simple_bayesian[31][25] = -6.63183; qual_mismatch_simple_bayesian[31][26] = -6.81137; qual_mismatch_simple_bayesian[31][27] = -6.98123; qual_mismatch_simple_bayesian[31][28] = -7.1403; qual_mismatch_simple_bayesian[31][29] = -7.28759; qual_mismatch_simple_bayesian[31][30] = -7.42233; qual_mismatch_simple_bayesian[31][31] = -7.54401; qual_mismatch_simple_bayesian[31][32] = -7.65246; qual_mismatch_simple_bayesian[31][33] = -7.74787; qual_mismatch_simple_bayesian[31][34] = -7.83072; qual_mismatch_simple_bayesian[31][35] = -7.90181; qual_mismatch_simple_bayesian[31][36] = -7.96211; qual_mismatch_simple_bayesian[31][37] = -8.01274; qual_mismatch_simple_bayesian[31][38] = -8.05488; qual_mismatch_simple_bayesian[31][39] = -8.08965; qual_mismatch_simple_bayesian[31][40] = -8.11817; qual_mismatch_simple_bayesian[31][41] = -8.14141; qual_mismatch_simple_bayesian[31][42] = -8.16027; qual_mismatch_simple_bayesian[31][43] = -8.1755; qual_mismatch_simple_bayesian[31][44] = -8.18777; qual_mismatch_simple_bayesian[31][45] = -8.19763; qual_mismatch_simple_bayesian[31][46] = -8.20553;
        qual_mismatch_simple_bayesian[32][0] = -1.09882; qual_mismatch_simple_bayesian[32][1] = -1.32892; qual_mismatch_simple_bayesian[32][2] = -1.55897; qual_mismatch_simple_bayesian[32][3] = -1.78897; qual_mismatch_simple_bayesian[32][4] = -2.0189; qual_mismatch_simple_bayesian[32][5] = -2.24875; qual_mismatch_simple_bayesian[32][6] = -2.47849; qual_mismatch_simple_bayesian[32][7] = -2.7081; qual_mismatch_simple_bayesian[32][8] = -2.93755; qual_mismatch_simple_bayesian[32][9] = -3.16678; qual_mismatch_simple_bayesian[32][10] = -3.39574; qual_mismatch_simple_bayesian[32][11] = -3.62438; qual_mismatch_simple_bayesian[32][12] = -3.8526; qual_mismatch_simple_bayesian[32][13] = -4.08029; qual_mismatch_simple_bayesian[32][14] = -4.30734; qual_mismatch_simple_bayesian[32][15] = -4.53356; qual_mismatch_simple_bayesian[32][16] = -4.75876; qual_mismatch_simple_bayesian[32][17] = -4.98269; qual_mismatch_simple_bayesian[32][18] = -5.20504; qual_mismatch_simple_bayesian[32][19] = -5.42542; qual_mismatch_simple_bayesian[32][20] = -5.64339; qual_mismatch_simple_bayesian[32][21] = -5.85838; qual_mismatch_simple_bayesian[32][22] = -6.06975; qual_mismatch_simple_bayesian[32][23] = -6.27673; qual_mismatch_simple_bayesian[32][24] = -6.47843; qual_mismatch_simple_bayesian[32][25] = -6.67385; qual_mismatch_simple_bayesian[32][26] = -6.86191; qual_mismatch_simple_bayesian[32][27] = -7.04146; qual_mismatch_simple_bayesian[32][28] = -7.21133; qual_mismatch_simple_bayesian[32][29] = -7.37041; qual_mismatch_simple_bayesian[32][30] = -7.51772; qual_mismatch_simple_bayesian[32][31] = -7.65246; qual_mismatch_simple_bayesian[32][32] = -7.77416; qual_mismatch_simple_bayesian[32][33] = -7.88263; qual_mismatch_simple_bayesian[32][34] = -7.97804; qual_mismatch_simple_bayesian[32][35] = -8.06091; qual_mismatch_simple_bayesian[32][36] = -8.132; qual_mismatch_simple_bayesian[32][37] = -8.19232; qual_mismatch_simple_bayesian[32][38] = -8.24296; qual_mismatch_simple_bayesian[32][39] = -8.2851; qual_mismatch_simple_bayesian[32][40] = -8.31988; qual_mismatch_simple_bayesian[32][41] = -8.3484; qual_mismatch_simple_bayesian[32][42] = -8.37165; qual_mismatch_simple_bayesian[32][43] = -8.39051; qual_mismatch_simple_bayesian[32][44] = -8.40575; qual_mismatch_simple_bayesian[32][45] = -8.41802; qual_mismatch_simple_bayesian[32][46] = -8.42788;
        qual_mismatch_simple_bayesian[33][0] = -1.09878; qual_mismatch_simple_bayesian[33][1] = -1.32891; qual_mismatch_simple_bayesian[33][2] = -1.559; qual_mismatch_simple_bayesian[33][3] = -1.78906; qual_mismatch_simple_bayesian[33][4] = -2.01906; qual_mismatch_simple_bayesian[33][5] = -2.24899; qual_mismatch_simple_bayesian[33][6] = -2.47884; qual_mismatch_simple_bayesian[33][7] = -2.70858; qual_mismatch_simple_bayesian[33][8] = -2.93819; qual_mismatch_simple_bayesian[33][9] = -3.16763; qual_mismatch_simple_bayesian[33][10] = -3.39686; qual_mismatch_simple_bayesian[33][11] = -3.62583; qual_mismatch_simple_bayesian[33][12] = -3.85447; qual_mismatch_simple_bayesian[33][13] = -4.08268; qual_mismatch_simple_bayesian[33][14] = -4.31038; qual_mismatch_simple_bayesian[33][15] = -4.53742; qual_mismatch_simple_bayesian[33][16] = -4.76365; qual_mismatch_simple_bayesian[33][17] = -4.98885; qual_mismatch_simple_bayesian[33][18] = -5.21278; qual_mismatch_simple_bayesian[33][19] = -5.43513; qual_mismatch_simple_bayesian[33][20] = -5.65552; qual_mismatch_simple_bayesian[33][21] = -5.87348; qual_mismatch_simple_bayesian[33][22] = -6.08848; qual_mismatch_simple_bayesian[33][23] = -6.29986; qual_mismatch_simple_bayesian[33][24] = -6.50683; qual_mismatch_simple_bayesian[33][25] = -6.70854; qual_mismatch_simple_bayesian[33][26] = -6.90396; qual_mismatch_simple_bayesian[33][27] = -7.09203; qual_mismatch_simple_bayesian[33][28] = -7.27159; qual_mismatch_simple_bayesian[33][29] = -7.44147; qual_mismatch_simple_bayesian[33][30] = -7.60056; qual_mismatch_simple_bayesian[33][31] = -7.74787; qual_mismatch_simple_bayesian[33][32] = -7.88263; qual_mismatch_simple_bayesian[33][33] = -8.00433; qual_mismatch_simple_bayesian[33][34] = -8.11281; qual_mismatch_simple_bayesian[33][35] = -8.20823; qual_mismatch_simple_bayesian[33][36] = -8.29111; qual_mismatch_simple_bayesian[33][37] = -8.36221; qual_mismatch_simple_bayesian[33][38] = -8.42253; qual_mismatch_simple_bayesian[33][39] = -8.47318; qual_mismatch_simple_bayesian[33][40] = -8.51533; qual_mismatch_simple_bayesian[33][41] = -8.55012; qual_mismatch_simple_bayesian[33][42] = -8.57864; qual_mismatch_simple_bayesian[33][43] = -8.60189; qual_mismatch_simple_bayesian[33][44] = -8.62076; qual_mismatch_simple_bayesian[33][45] = -8.636; qual_mismatch_simple_bayesian[33][46] = -8.64827;
        qual_mismatch_simple_bayesian[34][0] = -1.09874; qual_mismatch_simple_bayesian[34][1] = -1.3289; qual_mismatch_simple_bayesian[34][2] = -1.55903; qual_mismatch_simple_bayesian[34][3] = -1.78912; qual_mismatch_simple_bayesian[34][4] = -2.01918; qual_mismatch_simple_bayesian[34][5] = -2.24918; qual_mismatch_simple_bayesian[34][6] = -2.47911; qual_mismatch_simple_bayesian[34][7] = -2.70896; qual_mismatch_simple_bayesian[34][8] = -2.9387; qual_mismatch_simple_bayesian[34][9] = -3.16831; qual_mismatch_simple_bayesian[34][10] = -3.39775; qual_mismatch_simple_bayesian[34][11] = -3.62698; qual_mismatch_simple_bayesian[34][12] = -3.85595; qual_mismatch_simple_bayesian[34][13] = -4.08459; qual_mismatch_simple_bayesian[34][14] = -4.31281; qual_mismatch_simple_bayesian[34][15] = -4.5405; qual_mismatch_simple_bayesian[34][16] = -4.76755; qual_mismatch_simple_bayesian[34][17] = -4.99377; qual_mismatch_simple_bayesian[34][18] = -5.21897; qual_mismatch_simple_bayesian[34][19] = -5.44291; qual_mismatch_simple_bayesian[34][20] = -5.66525; qual_mismatch_simple_bayesian[34][21] = -5.88564; qual_mismatch_simple_bayesian[34][22] = -6.10361; qual_mismatch_simple_bayesian[34][23] = -6.31861; qual_mismatch_simple_bayesian[34][24] = -6.52999; qual_mismatch_simple_bayesian[34][25] = -6.73697; qual_mismatch_simple_bayesian[34][26] = -6.93867; qual_mismatch_simple_bayesian[34][27] = -7.13411; qual_mismatch_simple_bayesian[34][28] = -7.32218; qual_mismatch_simple_bayesian[34][29] = -7.50174; qual_mismatch_simple_bayesian[34][30] = -7.67163; qual_mismatch_simple_bayesian[34][31] = -7.83072; qual_mismatch_simple_bayesian[34][32] = -7.97804; qual_mismatch_simple_bayesian[34][33] = -8.11281; qual_mismatch_simple_bayesian[34][34] = -8.23452; qual_mismatch_simple_bayesian[34][35] = -8.34301; qual_mismatch_simple_bayesian[34][36] = -8.43844; qual_mismatch_simple_bayesian[34][37] = -8.52132; qual_mismatch_simple_bayesian[34][38] = -8.59243; qual_mismatch_simple_bayesian[34][39] = -8.65276; qual_mismatch_simple_bayesian[34][40] = -8.70341; qual_mismatch_simple_bayesian[34][41] = -8.74556; qual_mismatch_simple_bayesian[34][42] = -8.78036; qual_mismatch_simple_bayesian[34][43] = -8.80888; qual_mismatch_simple_bayesian[34][44] = -8.83214; qual_mismatch_simple_bayesian[34][45] = -8.851; qual_mismatch_simple_bayesian[34][46] = -8.86625;
        qual_mismatch_simple_bayesian[35][0] = -1.09872; qual_mismatch_simple_bayesian[35][1] = -1.32889; qual_mismatch_simple_bayesian[35][2] = -1.55905; qual_mismatch_simple_bayesian[35][3] = -1.78918; qual_mismatch_simple_bayesian[35][4] = -2.01927; qual_mismatch_simple_bayesian[35][5] = -2.24933; qual_mismatch_simple_bayesian[35][6] = -2.47933; qual_mismatch_simple_bayesian[35][7] = -2.70926; qual_mismatch_simple_bayesian[35][8] = -2.93911; qual_mismatch_simple_bayesian[35][9] = -3.16885; qual_mismatch_simple_bayesian[35][10] = -3.39846; qual_mismatch_simple_bayesian[35][11] = -3.6279; qual_mismatch_simple_bayesian[35][12] = -3.85713; qual_mismatch_simple_bayesian[35][13] = -4.0861; qual_mismatch_simple_bayesian[35][14] = -4.31474; qual_mismatch_simple_bayesian[35][15] = -4.54296; qual_mismatch_simple_bayesian[35][16] = -4.77065; qual_mismatch_simple_bayesian[35][17] = -4.9977; qual_mismatch_simple_bayesian[35][18] = -5.22392; qual_mismatch_simple_bayesian[35][19] = -5.44913; qual_mismatch_simple_bayesian[35][20] = -5.67306; qual_mismatch_simple_bayesian[35][21] = -5.89541; qual_mismatch_simple_bayesian[35][22] = -6.1158; qual_mismatch_simple_bayesian[35][23] = -6.33377; qual_mismatch_simple_bayesian[35][24] = -6.54877; qual_mismatch_simple_bayesian[35][25] = -6.76015; qual_mismatch_simple_bayesian[35][26] = -6.96713; qual_mismatch_simple_bayesian[35][27] = -7.16884; qual_mismatch_simple_bayesian[35][28] = -7.36428; qual_mismatch_simple_bayesian[35][29] = -7.55235; qual_mismatch_simple_bayesian[35][30] = -7.73192; qual_mismatch_simple_bayesian[35][31] = -7.90181; qual_mismatch_simple_bayesian[35][32] = -8.06091; qual_mismatch_simple_bayesian[35][33] = -8.20823; qual_mismatch_simple_bayesian[35][34] = -8.34301; qual_mismatch_simple_bayesian[35][35] = -8.46472; qual_mismatch_simple_bayesian[35][36] = -8.57322; qual_mismatch_simple_bayesian[35][37] = -8.66866; qual_mismatch_simple_bayesian[35][38] = -8.75154; qual_mismatch_simple_bayesian[35][39] = -8.82266; qual_mismatch_simple_bayesian[35][40] = -8.88299; qual_mismatch_simple_bayesian[35][41] = -8.93365; qual_mismatch_simple_bayesian[35][42] = -8.9758; qual_mismatch_simple_bayesian[35][43] = -9.0106; qual_mismatch_simple_bayesian[35][44] = -9.03913; qual_mismatch_simple_bayesian[35][45] = -9.06239; qual_mismatch_simple_bayesian[35][46] = -9.08126;
        qual_mismatch_simple_bayesian[36][0] = -1.0987; qual_mismatch_simple_bayesian[36][1] = -1.32889; qual_mismatch_simple_bayesian[36][2] = -1.55907; qual_mismatch_simple_bayesian[36][3] = -1.78922; qual_mismatch_simple_bayesian[36][4] = -2.01935; qual_mismatch_simple_bayesian[36][5] = -2.24945; qual_mismatch_simple_bayesian[36][6] = -2.4795; qual_mismatch_simple_bayesian[36][7] = -2.7095; qual_mismatch_simple_bayesian[36][8] = -2.93943; qual_mismatch_simple_bayesian[36][9] = -3.16928; qual_mismatch_simple_bayesian[36][10] = -3.39902; qual_mismatch_simple_bayesian[36][11] = -3.62863; qual_mismatch_simple_bayesian[36][12] = -3.85807; qual_mismatch_simple_bayesian[36][13] = -4.08731; qual_mismatch_simple_bayesian[36][14] = -4.31627; qual_mismatch_simple_bayesian[36][15] = -4.54491; qual_mismatch_simple_bayesian[36][16] = -4.77313; qual_mismatch_simple_bayesian[36][17] = -5.00083; qual_mismatch_simple_bayesian[36][18] = -5.22787; qual_mismatch_simple_bayesian[36][19] = -5.4541; qual_mismatch_simple_bayesian[36][20] = -5.6793; qual_mismatch_simple_bayesian[36][21] = -5.90323; qual_mismatch_simple_bayesian[36][22] = -6.12558; qual_mismatch_simple_bayesian[36][23] = -6.34597; qual_mismatch_simple_bayesian[36][24] = -6.56395; qual_mismatch_simple_bayesian[36][25] = -6.77895; qual_mismatch_simple_bayesian[36][26] = -6.99033; qual_mismatch_simple_bayesian[36][27] = -7.19731; qual_mismatch_simple_bayesian[36][28] = -7.39902; qual_mismatch_simple_bayesian[36][29] = -7.59446; qual_mismatch_simple_bayesian[36][30] = -7.78254; qual_mismatch_simple_bayesian[36][31] = -7.96211; qual_mismatch_simple_bayesian[36][32] = -8.132; qual_mismatch_simple_bayesian[36][33] = -8.29111; qual_mismatch_simple_bayesian[36][34] = -8.43844; qual_mismatch_simple_bayesian[36][35] = -8.57322; qual_mismatch_simple_bayesian[36][36] = -8.69494; qual_mismatch_simple_bayesian[36][37] = -8.80344; qual_mismatch_simple_bayesian[36][38] = -8.89888; qual_mismatch_simple_bayesian[36][39] = -8.98177; qual_mismatch_simple_bayesian[36][40] = -9.05289; qual_mismatch_simple_bayesian[36][41] = -9.11323; qual_mismatch_simple_bayesian[36][42] = -9.16389; qual_mismatch_simple_bayesian[36][43] = -9.20605; qual_mismatch_simple_bayesian[36][44] = -9.24085; qual_mismatch_simple_bayesian[36][45] = -9.26938; qual_mismatch_simple_bayesian[36][46] = -9.29264;
        qual_mismatch_simple_bayesian[37][0] = -1.09868; qual_mismatch_simple_bayesian[37][1] = -1.32889; qual_mismatch_simple_bayesian[37][2] = -1.55908; qual_mismatch_simple_bayesian[37][3] = -1.78926; qual_mismatch_simple_bayesian[37][4] = -2.01941; qual_mismatch_simple_bayesian[37][5] = -2.24954; qual_mismatch_simple_bayesian[37][6] = -2.47964; qual_mismatch_simple_bayesian[37][7] = -2.70969; qual_mismatch_simple_bayesian[37][8] = -2.93969; qual_mismatch_simple_bayesian[37][9] = -3.16962; qual_mismatch_simple_bayesian[37][10] = -3.39947; qual_mismatch_simple_bayesian[37][11] = -3.62921; qual_mismatch_simple_bayesian[37][12] = -3.85882; qual_mismatch_simple_bayesian[37][13] = -4.08826; qual_mismatch_simple_bayesian[37][14] = -4.3175; qual_mismatch_simple_bayesian[37][15] = -4.54646; qual_mismatch_simple_bayesian[37][16] = -4.7751; qual_mismatch_simple_bayesian[37][17] = -5.00332; qual_mismatch_simple_bayesian[37][18] = -5.23102; qual_mismatch_simple_bayesian[37][19] = -5.45806; qual_mismatch_simple_bayesian[37][20] = -5.68429; qual_mismatch_simple_bayesian[37][21] = -5.90949; qual_mismatch_simple_bayesian[37][22] = -6.13342; qual_mismatch_simple_bayesian[37][23] = -6.35578; qual_mismatch_simple_bayesian[37][24] = -6.57617; qual_mismatch_simple_bayesian[37][25] = -6.79414; qual_mismatch_simple_bayesian[37][26] = -7.00914; qual_mismatch_simple_bayesian[37][27] = -7.22052; qual_mismatch_simple_bayesian[37][28] = -7.42751; qual_mismatch_simple_bayesian[37][29] = -7.62922; qual_mismatch_simple_bayesian[37][30] = -7.82466; qual_mismatch_simple_bayesian[37][31] = -8.01274; qual_mismatch_simple_bayesian[37][32] = -8.19232; qual_mismatch_simple_bayesian[37][33] = -8.36221; qual_mismatch_simple_bayesian[37][34] = -8.52132; qual_mismatch_simple_bayesian[37][35] = -8.66866; qual_mismatch_simple_bayesian[37][36] = -8.80344; qual_mismatch_simple_bayesian[37][37] = -8.92516; qual_mismatch_simple_bayesian[37][38] = -9.03366; qual_mismatch_simple_bayesian[37][39] = -9.12911; qual_mismatch_simple_bayesian[37][40] = -9.21201; qual_mismatch_simple_bayesian[37][41] = -9.28313; qual_mismatch_simple_bayesian[37][42] = -9.34347; qual_mismatch_simple_bayesian[37][43] = -9.39414; qual_mismatch_simple_bayesian[37][44] = -9.43629; qual_mismatch_simple_bayesian[37][45] = -9.4711; qual_mismatch_simple_bayesian[37][46] = -9.49963;
        qual_mismatch_simple_bayesian[38][0] = -1.09867; qual_mismatch_simple_bayesian[38][1] = -1.32888; qual_mismatch_simple_bayesian[38][2] = -1.55909; qual_mismatch_simple_bayesian[38][3] = -1.78928; qual_mismatch_simple_bayesian[38][4] = -2.01946; qual_mismatch_simple_bayesian[38][5] = -2.24962; qual_mismatch_simple_bayesian[38][6] = -2.47974; qual_mismatch_simple_bayesian[38][7] = -2.70984; qual_mismatch_simple_bayesian[38][8] = -2.93989; qual_mismatch_simple_bayesian[38][9] = -3.16989; qual_mismatch_simple_bayesian[38][10] = -3.39982; qual_mismatch_simple_bayesian[38][11] = -3.62967; qual_mismatch_simple_bayesian[38][12] = -3.85942; qual_mismatch_simple_bayesian[38][13] = -4.08903; qual_mismatch_simple_bayesian[38][14] = -4.31847; qual_mismatch_simple_bayesian[38][15] = -4.5477; qual_mismatch_simple_bayesian[38][16] = -4.77667; qual_mismatch_simple_bayesian[38][17] = -5.0053; qual_mismatch_simple_bayesian[38][18] = -5.23352; qual_mismatch_simple_bayesian[38][19] = -5.46122; qual_mismatch_simple_bayesian[38][20] = -5.68827; qual_mismatch_simple_bayesian[38][21] = -5.91449; qual_mismatch_simple_bayesian[38][22] = -6.1397; qual_mismatch_simple_bayesian[38][23] = -6.36363; qual_mismatch_simple_bayesian[38][24] = -6.58598; qual_mismatch_simple_bayesian[38][25] = -6.80637; qual_mismatch_simple_bayesian[38][26] = -7.02435; qual_mismatch_simple_bayesian[38][27] = -7.23935; qual_mismatch_simple_bayesian[38][28] = -7.45073; qual_mismatch_simple_bayesian[38][29] = -7.65772; qual_mismatch_simple_bayesian[38][30] = -7.85943; qual_mismatch_simple_bayesian[38][31] = -8.05488; qual_mismatch_simple_bayesian[38][32] = -8.24296; qual_mismatch_simple_bayesian[38][33] = -8.42253; qual_mismatch_simple_bayesian[38][34] = -8.59243; qual_mismatch_simple_bayesian[38][35] = -8.75154; qual_mismatch_simple_bayesian[38][36] = -8.89888; qual_mismatch_simple_bayesian[38][37] = -9.03366; qual_mismatch_simple_bayesian[38][38] = -9.15539; qual_mismatch_simple_bayesian[38][39] = -9.2639; qual_mismatch_simple_bayesian[38][40] = -9.35935; qual_mismatch_simple_bayesian[38][41] = -9.44225; qual_mismatch_simple_bayesian[38][42] = -9.51338; qual_mismatch_simple_bayesian[38][43] = -9.57372; qual_mismatch_simple_bayesian[38][44] = -9.62438; qual_mismatch_simple_bayesian[38][45] = -9.66654; qual_mismatch_simple_bayesian[38][46] = -9.70135;
        qual_mismatch_simple_bayesian[39][0] = -1.09865; qual_mismatch_simple_bayesian[39][1] = -1.32888; qual_mismatch_simple_bayesian[39][2] = -1.5591; qual_mismatch_simple_bayesian[39][3] = -1.7893; qual_mismatch_simple_bayesian[39][4] = -2.0195; qual_mismatch_simple_bayesian[39][5] = -2.24967; qual_mismatch_simple_bayesian[39][6] = -2.47983; qual_mismatch_simple_bayesian[39][7] = -2.70996; qual_mismatch_simple_bayesian[39][8] = -2.94005; qual_mismatch_simple_bayesian[39][9] = -3.17011; qual_mismatch_simple_bayesian[39][10] = -3.40011; qual_mismatch_simple_bayesian[39][11] = -3.63004; qual_mismatch_simple_bayesian[39][12] = -3.85989; qual_mismatch_simple_bayesian[39][13] = -4.08963; qual_mismatch_simple_bayesian[39][14] = -4.31924; qual_mismatch_simple_bayesian[39][15] = -4.54868; qual_mismatch_simple_bayesian[39][16] = -4.77792; qual_mismatch_simple_bayesian[39][17] = -5.00688; qual_mismatch_simple_bayesian[39][18] = -5.23552; qual_mismatch_simple_bayesian[39][19] = -5.46374; qual_mismatch_simple_bayesian[39][20] = -5.69144; qual_mismatch_simple_bayesian[39][21] = -5.91848; qual_mismatch_simple_bayesian[39][22] = -6.14471; qual_mismatch_simple_bayesian[39][23] = -6.36991; qual_mismatch_simple_bayesian[39][24] = -6.59385; qual_mismatch_simple_bayesian[39][25] = -6.8162; qual_mismatch_simple_bayesian[39][26] = -7.03659; qual_mismatch_simple_bayesian[39][27] = -7.25456; qual_mismatch_simple_bayesian[39][28] = -7.46957; qual_mismatch_simple_bayesian[39][29] = -7.68095; qual_mismatch_simple_bayesian[39][30] = -7.88794; qual_mismatch_simple_bayesian[39][31] = -8.08965; qual_mismatch_simple_bayesian[39][32] = -8.2851; qual_mismatch_simple_bayesian[39][33] = -8.47318; qual_mismatch_simple_bayesian[39][34] = -8.65276; qual_mismatch_simple_bayesian[39][35] = -8.82266; qual_mismatch_simple_bayesian[39][36] = -8.98177; qual_mismatch_simple_bayesian[39][37] = -9.12911; qual_mismatch_simple_bayesian[39][38] = -9.2639; qual_mismatch_simple_bayesian[39][39] = -9.38563; qual_mismatch_simple_bayesian[39][40] = -9.49414; qual_mismatch_simple_bayesian[39][41] = -9.58959; qual_mismatch_simple_bayesian[39][42] = -9.67249; qual_mismatch_simple_bayesian[39][43] = -9.74362; qual_mismatch_simple_bayesian[39][44] = -9.80396; qual_mismatch_simple_bayesian[39][45] = -9.85463; qual_mismatch_simple_bayesian[39][46] = -9.8968;
        qual_mismatch_simple_bayesian[40][0] = -1.09865; qual_mismatch_simple_bayesian[40][1] = -1.32888; qual_mismatch_simple_bayesian[40][2] = -1.5591; qual_mismatch_simple_bayesian[40][3] = -1.78932; qual_mismatch_simple_bayesian[40][4] = -2.01953; qual_mismatch_simple_bayesian[40][5] = -2.24972; qual_mismatch_simple_bayesian[40][6] = -2.4799; qual_mismatch_simple_bayesian[40][7] = -2.71005; qual_mismatch_simple_bayesian[40][8] = -2.94018; qual_mismatch_simple_bayesian[40][9] = -3.17028; qual_mismatch_simple_bayesian[40][10] = -3.40033; qual_mismatch_simple_bayesian[40][11] = -3.63033; qual_mismatch_simple_bayesian[40][12] = -3.86026; qual_mismatch_simple_bayesian[40][13] = -4.09011; qual_mismatch_simple_bayesian[40][14] = -4.31986; qual_mismatch_simple_bayesian[40][15] = -4.54947; qual_mismatch_simple_bayesian[40][16] = -4.77891; qual_mismatch_simple_bayesian[40][17] = -5.00814; qual_mismatch_simple_bayesian[40][18] = -5.23711; qual_mismatch_simple_bayesian[40][19] = -5.46574; qual_mismatch_simple_bayesian[40][20] = -5.69396; qual_mismatch_simple_bayesian[40][21] = -5.92166; qual_mismatch_simple_bayesian[40][22] = -6.14871; qual_mismatch_simple_bayesian[40][23] = -6.37493; qual_mismatch_simple_bayesian[40][24] = -6.60014; qual_mismatch_simple_bayesian[40][25] = -6.82407; qual_mismatch_simple_bayesian[40][26] = -7.04642; qual_mismatch_simple_bayesian[40][27] = -7.26682; qual_mismatch_simple_bayesian[40][28] = -7.48479; qual_mismatch_simple_bayesian[40][29] = -7.6998; qual_mismatch_simple_bayesian[40][30] = -7.91118; qual_mismatch_simple_bayesian[40][31] = -8.11817; qual_mismatch_simple_bayesian[40][32] = -8.31988; qual_mismatch_simple_bayesian[40][33] = -8.51533; qual_mismatch_simple_bayesian[40][34] = -8.70341; qual_mismatch_simple_bayesian[40][35] = -8.88299; qual_mismatch_simple_bayesian[40][36] = -9.05289; qual_mismatch_simple_bayesian[40][37] = -9.21201; qual_mismatch_simple_bayesian[40][38] = -9.35935; qual_mismatch_simple_bayesian[40][39] = -9.49414; qual_mismatch_simple_bayesian[40][40] = -9.61587; qual_mismatch_simple_bayesian[40][41] = -9.72438; qual_mismatch_simple_bayesian[40][42] = -9.81984; qual_mismatch_simple_bayesian[40][43] = -9.90274; qual_mismatch_simple_bayesian[40][44] = -9.97387; qual_mismatch_simple_bayesian[40][45] = -10.0342; qual_mismatch_simple_bayesian[40][46] = -10.0849;
        qual_mismatch_simple_bayesian[41][0] = -1.09864; qual_mismatch_simple_bayesian[41][1] = -1.32888; qual_mismatch_simple_bayesian[41][2] = -1.55911; qual_mismatch_simple_bayesian[41][3] = -1.78934; qual_mismatch_simple_bayesian[41][4] = -2.01955; qual_mismatch_simple_bayesian[41][5] = -2.24976; qual_mismatch_simple_bayesian[41][6] = -2.47995; qual_mismatch_simple_bayesian[41][7] = -2.71013; qual_mismatch_simple_bayesian[41][8] = -2.94029; qual_mismatch_simple_bayesian[41][9] = -3.17041; qual_mismatch_simple_bayesian[41][10] = -3.40051; qual_mismatch_simple_bayesian[41][11] = -3.63056; qual_mismatch_simple_bayesian[41][12] = -3.86056; qual_mismatch_simple_bayesian[41][13] = -4.0905; qual_mismatch_simple_bayesian[41][14] = -4.32034; qual_mismatch_simple_bayesian[41][15] = -4.55009; qual_mismatch_simple_bayesian[41][16] = -4.7797; qual_mismatch_simple_bayesian[41][17] = -5.00914; qual_mismatch_simple_bayesian[41][18] = -5.23837; qual_mismatch_simple_bayesian[41][19] = -5.46734; qual_mismatch_simple_bayesian[41][20] = -5.69598; qual_mismatch_simple_bayesian[41][21] = -5.9242; qual_mismatch_simple_bayesian[41][22] = -6.15189; qual_mismatch_simple_bayesian[41][23] = -6.37894; qual_mismatch_simple_bayesian[41][24] = -6.60516; qual_mismatch_simple_bayesian[41][25] = -6.83037; qual_mismatch_simple_bayesian[41][26] = -7.0543; qual_mismatch_simple_bayesian[41][27] = -7.27666; qual_mismatch_simple_bayesian[41][28] = -7.49705; qual_mismatch_simple_bayesian[41][29] = -7.71502; qual_mismatch_simple_bayesian[41][30] = -7.93003; qual_mismatch_simple_bayesian[41][31] = -8.14141; qual_mismatch_simple_bayesian[41][32] = -8.3484; qual_mismatch_simple_bayesian[41][33] = -8.55012; qual_mismatch_simple_bayesian[41][34] = -8.74556; qual_mismatch_simple_bayesian[41][35] = -8.93365; qual_mismatch_simple_bayesian[41][36] = -9.11323; qual_mismatch_simple_bayesian[41][37] = -9.28313; qual_mismatch_simple_bayesian[41][38] = -9.44225; qual_mismatch_simple_bayesian[41][39] = -9.58959; qual_mismatch_simple_bayesian[41][40] = -9.72438; qual_mismatch_simple_bayesian[41][41] = -9.84612; qual_mismatch_simple_bayesian[41][42] = -9.95463; qual_mismatch_simple_bayesian[41][43] = -10.0501; qual_mismatch_simple_bayesian[41][44] = -10.133; qual_mismatch_simple_bayesian[41][45] = -10.2041; qual_mismatch_simple_bayesian[41][46] = -10.2645;
        qual_mismatch_simple_bayesian[42][0] = -1.09863; qual_mismatch_simple_bayesian[42][1] = -1.32888; qual_mismatch_simple_bayesian[42][2] = -1.55911; qual_mismatch_simple_bayesian[42][3] = -1.78935; qual_mismatch_simple_bayesian[42][4] = -2.01957; qual_mismatch_simple_bayesian[42][5] = -2.24979; qual_mismatch_simple_bayesian[42][6] = -2.48; qual_mismatch_simple_bayesian[42][7] = -2.71019; qual_mismatch_simple_bayesian[42][8] = -2.94037; qual_mismatch_simple_bayesian[42][9] = -3.17052; qual_mismatch_simple_bayesian[42][10] = -3.40065; qual_mismatch_simple_bayesian[42][11] = -3.63075; qual_mismatch_simple_bayesian[42][12] = -3.8608; qual_mismatch_simple_bayesian[42][13] = -4.0908; qual_mismatch_simple_bayesian[42][14] = -4.32073; qual_mismatch_simple_bayesian[42][15] = -4.55058; qual_mismatch_simple_bayesian[42][16] = -4.78032; qual_mismatch_simple_bayesian[42][17] = -5.00993; qual_mismatch_simple_bayesian[42][18] = -5.23938; qual_mismatch_simple_bayesian[42][19] = -5.46861; qual_mismatch_simple_bayesian[42][20] = -5.69758; qual_mismatch_simple_bayesian[42][21] = -5.92621; qual_mismatch_simple_bayesian[42][22] = -6.15443; qual_mismatch_simple_bayesian[42][23] = -6.38213; qual_mismatch_simple_bayesian[42][24] = -6.60917; qual_mismatch_simple_bayesian[42][25] = -6.8354; qual_mismatch_simple_bayesian[42][26] = -7.06061; qual_mismatch_simple_bayesian[42][27] = -7.28454; qual_mismatch_simple_bayesian[42][28] = -7.50689; qual_mismatch_simple_bayesian[42][29] = -7.72729; qual_mismatch_simple_bayesian[42][30] = -7.94526; qual_mismatch_simple_bayesian[42][31] = -8.16027; qual_mismatch_simple_bayesian[42][32] = -8.37165; qual_mismatch_simple_bayesian[42][33] = -8.57864; qual_mismatch_simple_bayesian[42][34] = -8.78036; qual_mismatch_simple_bayesian[42][35] = -8.9758; qual_mismatch_simple_bayesian[42][36] = -9.16389; qual_mismatch_simple_bayesian[42][37] = -9.34347; qual_mismatch_simple_bayesian[42][38] = -9.51338; qual_mismatch_simple_bayesian[42][39] = -9.67249; qual_mismatch_simple_bayesian[42][40] = -9.81984; qual_mismatch_simple_bayesian[42][41] = -9.95463; qual_mismatch_simple_bayesian[42][42] = -10.0764; qual_mismatch_simple_bayesian[42][43] = -10.1849; qual_mismatch_simple_bayesian[42][44] = -10.2803; qual_mismatch_simple_bayesian[42][45] = -10.3632; qual_mismatch_simple_bayesian[42][46] = -10.4344;
        qual_mismatch_simple_bayesian[43][0] = -1.09863; qual_mismatch_simple_bayesian[43][1] = -1.32887; qual_mismatch_simple_bayesian[43][2] = -1.55912; qual_mismatch_simple_bayesian[43][3] = -1.78935; qual_mismatch_simple_bayesian[43][4] = -2.01959; qual_mismatch_simple_bayesian[43][5] = -2.24981; qual_mismatch_simple_bayesian[43][6] = -2.48003; qual_mismatch_simple_bayesian[43][7] = -2.71024; qual_mismatch_simple_bayesian[43][8] = -2.94043; qual_mismatch_simple_bayesian[43][9] = -3.17061; qual_mismatch_simple_bayesian[43][10] = -3.40076; qual_mismatch_simple_bayesian[43][11] = -3.63089; qual_mismatch_simple_bayesian[43][12] = -3.86099; qual_mismatch_simple_bayesian[43][13] = -4.09104; qual_mismatch_simple_bayesian[43][14] = -4.32104; qual_mismatch_simple_bayesian[43][15] = -4.55097; qual_mismatch_simple_bayesian[43][16] = -4.78082; qual_mismatch_simple_bayesian[43][17] = -5.01056; qual_mismatch_simple_bayesian[43][18] = -5.24017; qual_mismatch_simple_bayesian[43][19] = -5.46962; qual_mismatch_simple_bayesian[43][20] = -5.69885; qual_mismatch_simple_bayesian[43][21] = -5.92782; qual_mismatch_simple_bayesian[43][22] = -6.15645; qual_mismatch_simple_bayesian[43][23] = -6.38467; qual_mismatch_simple_bayesian[43][24] = -6.61237; qual_mismatch_simple_bayesian[43][25] = -6.83942; qual_mismatch_simple_bayesian[43][26] = -7.06564; qual_mismatch_simple_bayesian[43][27] = -7.29085; qual_mismatch_simple_bayesian[43][28] = -7.51478; qual_mismatch_simple_bayesian[43][29] = -7.73713; qual_mismatch_simple_bayesian[43][30] = -7.95753; qual_mismatch_simple_bayesian[43][31] = -8.1755; qual_mismatch_simple_bayesian[43][32] = -8.39051; qual_mismatch_simple_bayesian[43][33] = -8.60189; qual_mismatch_simple_bayesian[43][34] = -8.80888; qual_mismatch_simple_bayesian[43][35] = -9.0106; qual_mismatch_simple_bayesian[43][36] = -9.20605; qual_mismatch_simple_bayesian[43][37] = -9.39414; qual_mismatch_simple_bayesian[43][38] = -9.57372; qual_mismatch_simple_bayesian[43][39] = -9.74362; qual_mismatch_simple_bayesian[43][40] = -9.90274; qual_mismatch_simple_bayesian[43][41] = -10.0501; qual_mismatch_simple_bayesian[43][42] = -10.1849; qual_mismatch_simple_bayesian[43][43] = -10.3066; qual_mismatch_simple_bayesian[43][44] = -10.4151; qual_mismatch_simple_bayesian[43][45] = -10.5106; qual_mismatch_simple_bayesian[43][46] = -10.5935;
        qual_mismatch_simple_bayesian[44][0] = -1.09863; qual_mismatch_simple_bayesian[44][1] = -1.32887; qual_mismatch_simple_bayesian[44][2] = -1.55912; qual_mismatch_simple_bayesian[44][3] = -1.78936; qual_mismatch_simple_bayesian[44][4] = -2.0196; qual_mismatch_simple_bayesian[44][5] = -2.24983; qual_mismatch_simple_bayesian[44][6] = -2.48006; qual_mismatch_simple_bayesian[44][7] = -2.71028; qual_mismatch_simple_bayesian[44][8] = -2.94048; qual_mismatch_simple_bayesian[44][9] = -3.17068; qual_mismatch_simple_bayesian[44][10] = -3.40085; qual_mismatch_simple_bayesian[44][11] = -3.63101; qual_mismatch_simple_bayesian[44][12] = -3.86114; qual_mismatch_simple_bayesian[44][13] = -4.09123; qual_mismatch_simple_bayesian[44][14] = -4.32128; qual_mismatch_simple_bayesian[44][15] = -4.55128; qual_mismatch_simple_bayesian[44][16] = -4.78122; qual_mismatch_simple_bayesian[44][17] = -5.01107; qual_mismatch_simple_bayesian[44][18] = -5.24081; qual_mismatch_simple_bayesian[44][19] = -5.47042; qual_mismatch_simple_bayesian[44][20] = -5.69986; qual_mismatch_simple_bayesian[44][21] = -5.92909; qual_mismatch_simple_bayesian[44][22] = -6.15806; qual_mismatch_simple_bayesian[44][23] = -6.3867; qual_mismatch_simple_bayesian[44][24] = -6.61492; qual_mismatch_simple_bayesian[44][25] = -6.84262; qual_mismatch_simple_bayesian[44][26] = -7.06966; qual_mismatch_simple_bayesian[44][27] = -7.29589; qual_mismatch_simple_bayesian[44][28] = -7.52109; qual_mismatch_simple_bayesian[44][29] = -7.74503; qual_mismatch_simple_bayesian[44][30] = -7.96738; qual_mismatch_simple_bayesian[44][31] = -8.18777; qual_mismatch_simple_bayesian[44][32] = -8.40575; qual_mismatch_simple_bayesian[44][33] = -8.62076; qual_mismatch_simple_bayesian[44][34] = -8.83214; qual_mismatch_simple_bayesian[44][35] = -9.03913; qual_mismatch_simple_bayesian[44][36] = -9.24085; qual_mismatch_simple_bayesian[44][37] = -9.43629; qual_mismatch_simple_bayesian[44][38] = -9.62438; qual_mismatch_simple_bayesian[44][39] = -9.80396; qual_mismatch_simple_bayesian[44][40] = -9.97387; qual_mismatch_simple_bayesian[44][41] = -10.133; qual_mismatch_simple_bayesian[44][42] = -10.2803; qual_mismatch_simple_bayesian[44][43] = -10.4151; qual_mismatch_simple_bayesian[44][44] = -10.5369; qual_mismatch_simple_bayesian[44][45] = -10.6454; qual_mismatch_simple_bayesian[44][46] = -10.7408;
        qual_mismatch_simple_bayesian[45][0] = -1.09862; qual_mismatch_simple_bayesian[45][1] = -1.32887; qual_mismatch_simple_bayesian[45][2] = -1.55912; qual_mismatch_simple_bayesian[45][3] = -1.78937; qual_mismatch_simple_bayesian[45][4] = -2.01961; qual_mismatch_simple_bayesian[45][5] = -2.24985; qual_mismatch_simple_bayesian[45][6] = -2.48008; qual_mismatch_simple_bayesian[45][7] = -2.71031; qual_mismatch_simple_bayesian[45][8] = -2.94052; qual_mismatch_simple_bayesian[45][9] = -3.17073; qual_mismatch_simple_bayesian[45][10] = -3.40092; qual_mismatch_simple_bayesian[45][11] = -3.6311; qual_mismatch_simple_bayesian[45][12] = -3.86126; qual_mismatch_simple_bayesian[45][13] = -4.09138; qual_mismatch_simple_bayesian[45][14] = -4.32148; qual_mismatch_simple_bayesian[45][15] = -4.55153; qual_mismatch_simple_bayesian[45][16] = -4.78153; qual_mismatch_simple_bayesian[45][17] = -5.01147; qual_mismatch_simple_bayesian[45][18] = -5.24131; qual_mismatch_simple_bayesian[45][19] = -5.47106; qual_mismatch_simple_bayesian[45][20] = -5.70067; qual_mismatch_simple_bayesian[45][21] = -5.93011; qual_mismatch_simple_bayesian[45][22] = -6.15934; qual_mismatch_simple_bayesian[45][23] = -6.38831; qual_mismatch_simple_bayesian[45][24] = -6.61695; qual_mismatch_simple_bayesian[45][25] = -6.84517; qual_mismatch_simple_bayesian[45][26] = -7.07286; qual_mismatch_simple_bayesian[45][27] = -7.29991; qual_mismatch_simple_bayesian[45][28] = -7.52614; qual_mismatch_simple_bayesian[45][29] = -7.75134; qual_mismatch_simple_bayesian[45][30] = -7.97528; qual_mismatch_simple_bayesian[45][31] = -8.19763; qual_mismatch_simple_bayesian[45][32] = -8.41802; qual_mismatch_simple_bayesian[45][33] = -8.636; qual_mismatch_simple_bayesian[45][34] = -8.851; qual_mismatch_simple_bayesian[45][35] = -9.06239; qual_mismatch_simple_bayesian[45][36] = -9.26938; qual_mismatch_simple_bayesian[45][37] = -9.4711; qual_mismatch_simple_bayesian[45][38] = -9.66654; qual_mismatch_simple_bayesian[45][39] = -9.85463; qual_mismatch_simple_bayesian[45][40] = -10.0342; qual_mismatch_simple_bayesian[45][41] = -10.2041; qual_mismatch_simple_bayesian[45][42] = -10.3632; qual_mismatch_simple_bayesian[45][43] = -10.5106; qual_mismatch_simple_bayesian[45][44] = -10.6454; qual_mismatch_simple_bayesian[45][45] = -10.7671; qual_mismatch_simple_bayesian[45][46] = -10.8756;
        qual_mismatch_simple_bayesian[46][0] = -1.09862; qual_mismatch_simple_bayesian[46][1] = -1.32887; qual_mismatch_simple_bayesian[46][2] = -1.55912; qual_mismatch_simple_bayesian[46][3] = -1.78937; qual_mismatch_simple_bayesian[46][4] = -2.01962; qual_mismatch_simple_bayesian[46][5] = -2.24986; qual_mismatch_simple_bayesian[46][6] = -2.4801; qual_mismatch_simple_bayesian[46][7] = -2.71033; qual_mismatch_simple_bayesian[46][8] = -2.94056; qual_mismatch_simple_bayesian[46][9] = -3.17077; qual_mismatch_simple_bayesian[46][10] = -3.40098; qual_mismatch_simple_bayesian[46][11] = -3.63117; qual_mismatch_simple_bayesian[46][12] = -3.86135; qual_mismatch_simple_bayesian[46][13] = -4.09151; qual_mismatch_simple_bayesian[46][14] = -4.32163; qual_mismatch_simple_bayesian[46][15] = -4.55173; qual_mismatch_simple_bayesian[46][16] = -4.78178; qual_mismatch_simple_bayesian[46][17] = -5.01178; qual_mismatch_simple_bayesian[46][18] = -5.24172; qual_mismatch_simple_bayesian[46][19] = -5.47156; qual_mismatch_simple_bayesian[46][20] = -5.70131; qual_mismatch_simple_bayesian[46][21] = -5.93092; qual_mismatch_simple_bayesian[46][22] = -6.16036; qual_mismatch_simple_bayesian[46][23] = -6.38959; qual_mismatch_simple_bayesian[46][24] = -6.61856; qual_mismatch_simple_bayesian[46][25] = -6.8472; qual_mismatch_simple_bayesian[46][26] = -7.07542; qual_mismatch_simple_bayesian[46][27] = -7.30311; qual_mismatch_simple_bayesian[46][28] = -7.53016; qual_mismatch_simple_bayesian[46][29] = -7.75639; qual_mismatch_simple_bayesian[46][30] = -7.98159; qual_mismatch_simple_bayesian[46][31] = -8.20553; qual_mismatch_simple_bayesian[46][32] = -8.42788; qual_mismatch_simple_bayesian[46][33] = -8.64827; qual_mismatch_simple_bayesian[46][34] = -8.86625; qual_mismatch_simple_bayesian[46][35] = -9.08126; qual_mismatch_simple_bayesian[46][36] = -9.29264; qual_mismatch_simple_bayesian[46][37] = -9.49963; qual_mismatch_simple_bayesian[46][38] = -9.70135; qual_mismatch_simple_bayesian[46][39] = -9.8968; qual_mismatch_simple_bayesian[46][40] = -10.0849; qual_mismatch_simple_bayesian[46][41] = -10.2645; qual_mismatch_simple_bayesian[46][42] = -10.4344; qual_mismatch_simple_bayesian[46][43] = -10.5935; qual_mismatch_simple_bayesian[46][44] = -10.7408; qual_mismatch_simple_bayesian[46][45] = -10.8756; qual_mismatch_simple_bayesian[46][46] = -10.9974;

        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "loadQmatchValues");
        exit(1);
    }
}
//**********************************************************************************************************************

