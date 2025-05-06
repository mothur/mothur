//
//  makecontigscommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/15/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "makecontigscommand.h"
#include "contigsreport.hpp"
#include "counttable.h"

//**************************************************************************************************

/**
 * Convert the probability to a quality score.
 */

double convertProbToQ(double prob){
		try {
        return round(-10*log10(prob));
    }
    catch(exception& e) {
        MothurOut* m; m = MothurOut::getInstance();
        m->errorOut(e, "MakeContigsCommand", "convertProbToQ");
        exit(1);
    }
}

//**************************************************************************************************

/**
 * Convert the quality score to a probability.
 */

double convertQToProb(double Q){
		try {
            double value = pow(10,(-Q/10));
            if (isnan(value) || isinf(value)) { value = 0.0; }
            return value;
    }
    catch(exception& e) {
        MothurOut* m; m = MothurOut::getInstance();
        m->errorOut(e, "MakeContigsCommand", "convertQToProb");
        exit(1);
    }
}

//**************************************************************************************************

int loadQmatchValues(vector< vector<double> >& qual_match_simple_bayesian, vector< vector<double> >& qual_mismatch_simple_bayesian){
    try {
				vector<double> probs(47);
				for(int i=0;i<probs.size();i++){
					probs[i] = convertQToProb(i);
				}

        //qual_match_simple_bayesian - naturally symmetric
        //qual_mismatch_simple_bayesian - force symmetry
				for(int i=0;i<qual_match_simple_bayesian.size(); i++){
					for(int j=0;j<=i;j++){
						qual_match_simple_bayesian[i][j] = convertProbToQ((probs[i]*probs[j]/3)/(1-probs[i]-probs[j]+4*probs[i]*probs[j]/3));

						qual_mismatch_simple_bayesian[i][j] = convertProbToQ(probs[i]*(1-probs[j]/3)/(probs[i]+probs[j]-4*probs[i]*probs[j]/3));
					}
					qual_mismatch_simple_bayesian[i][i] = 2.0;
				}

				for(int i=0;i<qual_match_simple_bayesian.size(); i++){
					for(int j=0;j<i;j++){
						qual_match_simple_bayesian[j][i] = qual_match_simple_bayesian[i][j];
						qual_mismatch_simple_bayesian[j][i] = qual_mismatch_simple_bayesian[i][j];
					}
				}

        return 0;
    }
    catch(exception& e) {
        MothurOut* m; m = MothurOut::getInstance();
        m->errorOut(e, "MakeContigsCommand", "loadQmatchValues");
        exit(1);
    }
}

//**************************************************************************************************


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
        CommandParameter pqfile("qfile", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pqfile);
				CommandParameter ppdiffs("pdiffs", "Number", "", "0", "", "", "","",false,false,true); parameters.push_back(ppdiffs);
				CommandParameter pbdiffs("bdiffs", "Number", "", "0", "", "", "","",false,false,true); parameters.push_back(pbdiffs);
        CommandParameter ptdiffs("tdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(ptdiffs);
        CommandParameter preorient("checkorient", "Boolean", "", "T", "", "", "","",false,false,true); parameters.push_back(preorient);
        CommandParameter palign("align", "Multiple", "needleman-gotoh-kmer", "needleman", "", "", "","",false,false); parameters.push_back(palign);
        CommandParameter pallfiles("allfiles", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pallfiles);
        CommandParameter ptrimoverlap("trimoverlap", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(ptrimoverlap);
				CommandParameter pmatch("match", "Number", "", "1.0", "", "", "","",false,false); parameters.push_back(pmatch);
				CommandParameter pmismatch("mismatch", "Number", "", "-1.0", "", "", "","",false,false); parameters.push_back(pmismatch);
				CommandParameter pgapopen("gapopen", "Number", "", "-2.0", "", "", "","",false,false); parameters.push_back(pgapopen);
				CommandParameter pgapextend("gapextend", "Number", "", "-1.0", "", "", "","",false,false); parameters.push_back(pgapextend);
        CommandParameter pthreshold("insert", "Number", "", "20", "", "", "","",false,false); parameters.push_back(pthreshold);
        CommandParameter pdeltaq("deltaq", "Number", "", "6", "", "", "","",false,false); parameters.push_back(pdeltaq);
        CommandParameter maxee("maxee", "Number", "", "10000", "", "", "","",false,false); parameters.push_back(maxee);
				CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
        CommandParameter pformat("format", "Multiple", "sanger-illumina-solexa-illumina1.8+", "illumina1.8+", "", "", "","",false,false,true); parameters.push_back(pformat);
        CommandParameter pksize("ksize", "Number", "", "8", "", "", "","",false,false); parameters.push_back(pksize);
        CommandParameter pmaxambig("maxambig", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pmaxambig);
        CommandParameter pmaxhomop("maxhomop", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pmaxhomop);
        CommandParameter pmaxlength("maxlength", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pmaxlength);

        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);

        abort = false; calledHelp = false;
        createFileGroup = false; createOligosGroup = false; gz = false;
        
        //initialize outputTypes
        vector<string> tempOutNames;
        outputTypes["fasta"] = tempOutNames;
        outputTypes["qfile"] = tempOutNames;
        outputTypes["report"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
        
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
		helpString += "The make.contigs command reads a file, forward fastq file and a reverse fastq file or forward fasta and reverse fasta files and outputs a fasta file. \n";
        helpString += "If an oligos file is provided barcodes and primers will be trimmed, and a count file will be created.\n";
        helpString += "If a forward index or reverse index file is provided barcodes be trimmed, and a group file will be created. The oligos parameter is required if an index file is given.\n";
		helpString += "The make.contigs command parameters are file, ffastq, rfastq, ffasta, rfasta, fqfile, rqfile, oligos, findex, rindex, qfile, format, tdiffs, bdiffs, pdiffs, align, match, mismatch, gapopen, gapextend, insert, deltaq, maxee, allfiles and processors.\n";
		helpString += "The ffastq and rfastq, file, or ffasta and rfasta parameters are required.\n";
        helpString += "The file parameter is 2, 3 or 4 column file containing the forward fastq files in the first column and their matching reverse fastq files in the second column, or a groupName then forward fastq file and reverse fastq file, or forward fastq file then reverse fastq then forward index and reverse index file.  If you only have one index file add 'none' for the other one.  Mothur will process each pair and create a combined fasta and report file with all the sequences.\n";
        helpString += "The ffastq and rfastq parameters are used to provide a forward fastq and reverse fastq file to process.  If you provide one, you must provide the other.\n";
        helpString += "The ffasta and rfasta parameters are used to provide a forward fasta and reverse fasta file to process.  If you provide one, you must provide the other.\n";
        helpString += "The fqfile and rqfile parameters are used to provide a forward quality and reverse quality files to process with the ffasta and rfasta parameters.  If you provide one, you must provide the other.\n";
		helpString += "The format parameter is used to indicate whether your sequences are sanger, solexa, illumina1.8+ or illumina, default=illumina1.8+.\n";
        helpString += "The findex and rindex parameters are used to provide a forward index and reverse index files to process.  \n";
        helpString += "The qfile parameter is used to indicate you want a quality file assembled. Default=f. NOTE: The assembled quality scores outputted by mothur cannot be used for downstream quality screening. The score calculations are modeled after pandseq's method. Here's a link to the explanation from their documentation, https://github.com/neufeld/pandaseq#the-scores-of-the-output-bases-seem-really-low-whats-wrong. \n";
        helpString += "The align parameter allows you to specify the alignment method to use.  Your options are: kmer, gotoh and needleman. The default is needleman.\n";
        helpString += "The ksize parameter allows you to set the kmer size if you are doing align=kmer. Default=8.\n";
        helpString += "The tdiffs parameter is used to specify the total number of differences allowed in the sequence. The default is pdiffs + bdiffs + sdiffs + ldiffs.\n";
		helpString += "The bdiffs parameter is used to specify the number of differences allowed in the barcode. The default is 0.\n";
		helpString += "The pdiffs parameter is used to specify the number of differences allowed in the primer. The default is 0.\n";
		helpString += "The match parameter allows you to specify the bonus for having the same base. The default is 1.0.\n";
		helpString += "The mistmatch parameter allows you to specify the penalty for having different bases.  The default is -1.0.\n";
        helpString += "The checkorient parameter will look for the reverse compliment of the barcode or primer in the sequence. If found the sequence is flipped. The default is true.\n";
        helpString += "The deltaq parameter allows you to specify the delta allowed between quality scores of a mismatched base.  For example in the overlap, if deltaq=5 and in the alignment seqA, pos 200 has a quality score of 30 and the same position in seqB has a quality score of 20, you take the base from seqA (30-20 >= 5).  If the quality score in seqB is 28 then the base in the consensus will be an N (30-28<5). The default is 6.\n";
				helpString += "The maxee parameter allows you to specify the maximum number of errors to allow in a sequence. Makes sense to use with deltaq=0. This number is a decimal number. The expected number of errors is based on Edgar's approach used in USEARCH/VSEARCH.";
		helpString += "The gapopen parameter allows you to specify the penalty for opening a gap in an alignment. The default is -2.0.\n";
		helpString += "The gapextend parameter allows you to specify the penalty for extending a gap in an alignment.  The default is -1.0.\n";
        helpString += "The insert parameter allows you to set a quality scores threshold. In the case where we are trying to decide whether to keep a base or remove it because the base is compared to a gap in the other fragment, if the base has a quality score equal to or below the threshold we eliminate it. Default=20.\n";
        helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is all available.\n";
        helpString += "The allfiles parameter will create separate group and fasta file for each grouping. The default is F.\n";

        helpString += "The trimoverlap parameter allows you to trim the sequences to only the overlapping section. The default is F.\n";
        helpString += "The maxambig parameter allows you to set the maximum number of ambiguous bases allowed. The default is -1, meaning ignore.\n";
        helpString += "The maxhomop parameter allows you to set a maximum homopolymer length. The default is -1, meaning ignore.\n";
        helpString += "The maxlength parameter allows you to set a maximum length of your sequences. The default is -1, meaning ignore.\n";
        helpString += "The make.contigs command should be in the following format: \n";
		helpString += "make.contigs(ffastq=yourForwardFastqFile, rfastq=yourReverseFastqFile, align=yourAlignmentMethod) \n";
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
        else if (type == "count") {  pattern = "[filename],[tag],contigs.count_table"; }
        else if (type == "report") {  pattern = "[filename],[tag],contigs_report"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }

        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
MakeContigsCommand::MakeContigsCommand(string option) : Command()  {
	try {

		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }

		else {
			OptionParser parser(option, setParameters());
			map<string, string> parameters = parser.getParameters();

			ValidParameters validParameter;
            ffastqfile = validParameter.validFile(parameters, "ffastq");
			if (ffastqfile == "not open") {  abort = true; }
			else if (ffastqfile == "not found") { ffastqfile = ""; }

			rfastqfile = validParameter.validFile(parameters, "rfastq");
			if (rfastqfile == "not open") {  abort = true; }
			else if (rfastqfile == "not found") { rfastqfile = "";  }

            ffastafile = validParameter.validFile(parameters, "ffasta");
			if (ffastafile == "not open") {  abort = true; }
			else if (ffastafile == "not found") { ffastafile = ""; }

			rfastafile = validParameter.validFile(parameters, "rfasta");
			if (rfastafile == "not open") {  abort = true; }
			else if (rfastafile == "not found") { rfastafile = "";  }

            fqualfile = validParameter.validFile(parameters, "fqfile");
			if (fqualfile == "not open") {  abort = true; }
			else if (fqualfile == "not found") { fqualfile = ""; }

			rqualfile = validParameter.validFile(parameters, "rqfile");
			if (rqualfile == "not open") {  abort = true; }
			else if (rqualfile == "not found") { rqualfile = "";  }

            file = validParameter.validFile(parameters, "file");
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

            oligosfile = validParameter.validFile(parameters, "oligos");
			if (oligosfile == "not found")      {	oligosfile = "";	}
			else if(oligosfile == "not open")   {	abort = true;       }
			else {	 current->setOligosFile(oligosfile);		}

            findexfile = validParameter.validFile(parameters, "findex");
			if (findexfile == "not found")      {	findexfile = "";	}
			else if(findexfile == "not open")   {	abort = true;       }

            rindexfile = validParameter.validFile(parameters, "rindex");
			if (rindexfile == "not found")      {	rindexfile = "";	}
			else if(rindexfile == "not open")   {	abort = true;       }

            if ((rindexfile != "") || (findexfile != "")) {
				if (oligosfile == ""){
					oligosfile = current->getOligosFile();
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

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.valid(parameters, "match");		if (temp == "not found"){	temp = "1.0";			}
			util.mothurConvert(temp, match);

			temp = validParameter.valid(parameters, "mismatch");		if (temp == "not found"){	temp = "-1.0";			}
			util.mothurConvert(temp, misMatch);
            if (misMatch > 0) { m->mothurOut("[ERROR]: mismatch must be negative.\n"); abort=true; }

			temp = validParameter.valid(parameters, "gapopen");		if (temp == "not found"){	temp = "-2.0";			}
			util.mothurConvert(temp, gapOpen);
            if (gapOpen > 0) { m->mothurOut("[ERROR]: gapopen must be negative.\n"); abort=true; }

			temp = validParameter.valid(parameters, "gapextend");	if (temp == "not found"){	temp = "-1.0";			}
			util.mothurConvert(temp, gapExtend);
            if (gapExtend > 0) { m->mothurOut("[ERROR]: gapextend must be negative.\n"); abort=true; }

            temp = validParameter.valid(parameters, "insert");	if (temp == "not found"){	temp = "20";			}
			util.mothurConvert(temp, insert);
            if ((insert < 0) || (insert > 40)) { m->mothurOut("[ERROR]: insert must be between 0 and 40.\n"); abort=true; }

            temp = validParameter.valid(parameters, "deltaq");	if (temp == "not found"){	temp = "6";			}
			util.mothurConvert(temp, deltaq);

            temp = validParameter.valid(parameters, "maxee");	if (temp == "not found"){	temp = "10000";			}
			util.mothurConvert(temp, maxee);

			temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);

            temp = validParameter.valid(parameters, "bdiffs");		if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, bdiffs);

			temp = validParameter.valid(parameters, "pdiffs");		if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, pdiffs);

			temp = validParameter.valid(parameters, "tdiffs");		if (temp == "not found") { int tempTotal = pdiffs + bdiffs;  temp = toString(tempTotal); }
			util.mothurConvert(temp, tdiffs);

			if(tdiffs == 0){	tdiffs = bdiffs + pdiffs;	} 

            temp = validParameter.valid(parameters, "allfiles");		if (temp == "not found") { temp = "F"; }
			allFiles = util.isTrue(temp);

            temp = validParameter.valid(parameters, "ksize");	if (temp == "not found"){	temp = "8";			}
            util.mothurConvert(temp, kmerSize);

            temp = validParameter.valid(parameters, "trimoverlap");		if (temp == "not found") { temp = "F"; }
			trimOverlap = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "qfile");
            if (temp == "not found") { temp = "F"; }
            else { temp = util.getSimpleName(temp); }
            makeQualFile = util.isTrue(temp);

			align = validParameter.valid(parameters, "align");		if (align == "not found"){	align = "needleman";	}
			if ((align != "needleman") && (align != "gotoh") && (align != "kmer")) { m->mothurOut(align + " is not a valid alignment method. Options are kmer, needleman or gotoh. I will use needleman.\n");  align = "needleman"; }

            format = validParameter.valid(parameters, "format");		if (format == "not found"){	format = "illumina1.8+";	}

            if ((format != "sanger") && (format != "illumina") && (format != "illumina1.8+") && (format != "solexa"))  {
				m->mothurOut(format + " is not a valid format. Your format choices are sanger, solexa, illumina1.8+ and illumina, aborting." ); m->mothurOutEndLine();
				abort=true;
			}

            temp = validParameter.valid(parameters, "checkorient");		if (temp == "not found") { temp = "T"; }
			reorient = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "maxambig");        if (temp == "not found") { temp = "-1"; }
            util.mothurConvert(temp, maxAmbig);

            temp = validParameter.valid(parameters, "maxhomop");        if (temp == "not found") { temp = "-1"; }
            util.mothurConvert(temp, maxHomoP);

            temp = validParameter.valid(parameters, "maxlength");    if (temp == "not found") { temp = "-1"; }
            util.mothurConvert(temp, maxLength);
            
            if ((maxLength == -1) && (maxHomoP == -1) && (maxAmbig == -1)) { screenSequences = false; }
            else { screenSequences = true; }
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
        bool debugIndex = false;
        if (debugIndex) { //allows you to run the oligos and index file independantly to check for barcode issues. make.contigs(findex=yourIndexFile, bdiffs=1, oligos=yourOligosFile, checkorient=t). just used for user support
            debugFunction();
        }
        
        if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        unsigned long long numReads = 0;
        long start = time(nullptr);
        string outFastaFile, outScrapFastaFile, outQualFile, outScrapQualFile, outMisMatchFile, inputFile;

        if (file != "")                                     { numReads = processMultipleFileOption(outFastaFile, outMisMatchFile);   inputFile = file; }
        else if ((ffastqfile != "") || (ffastafile != ""))  {
            numReads = processSingleFileOption(outFastaFile, outScrapFastaFile, outQualFile, outScrapQualFile, outMisMatchFile, "");
            inputFile = ffastqfile;
            if (ffastafile != "") { inputFile = ffastafile; }
        } else {  return 0; }

        if (groupMap.size() != 0) {
            string thisOutputDir = outputdir;
            if (outputdir == "") {  thisOutputDir = util.hasPath(inputFile); }
            map<string, string> vars;
            vars["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(inputFile));
            vars["[tag]"] = "";
            string outputCountFileName = getOutputFileName("count",vars);
            outputNames.push_back(outputCountFileName); outputTypes["count"].push_back(outputCountFileName);
            createCountFile(outputCountFileName, outFastaFile);
        }

        //add headers to mismatch file
        ofstream out; util.openOutputFile(outMisMatchFile+".temp", out);
        ContigsReport report; report.printHeaders(out); out.close();//print Headers
        util.appendFilesFront(outMisMatchFile+".temp", outMisMatchFile); //removes temp

        if (m->getControl_pressed()) {	for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0;	}

        string currentFasta = "";
        itTypes = outputTypes.find("fasta");
        if (itTypes != outputTypes.end()) { if ((itTypes->second).size() != 0) { currentFasta = (itTypes->second)[0]; current->setFastaFile(currentFasta); } }

        string currentCount = "";
        itTypes = outputTypes.find("count");
        if (itTypes != outputTypes.end()) { if ((itTypes->second).size() != 0) { currentCount = (itTypes->second)[0]; current->setCountFile(currentCount); } }

        string currentQual = "";
        itTypes = outputTypes.find("qfile");
        if (itTypes != outputTypes.end()) { if ((itTypes->second).size() != 0) { currentQual = (itTypes->second)[0]; current->setQualFile(currentQual); } }

        string currentReport = "";
        itTypes = outputTypes.find("report");
        if (itTypes != outputTypes.end()) { if ((itTypes->second).size() != 0) { currentReport = (itTypes->second)[0]; current->setContigsReportFile(currentReport); } }

        if (m->getControl_pressed()) {	for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0;	}

		//output group counts
		int total = 0;
		if (groupCounts.size() != 0) {  m->mothurOut("\nGroup count: \n");  }
		for (map<string, int>::iterator it = groupCounts.begin(); it != groupCounts.end(); it++) { total += it->second; m->mothurOut(it->first + "\t" + toString(it->second) + "\n"); }
		if (total != 0) { m->mothurOut("\nTotal of all groups is " + toString(total) + "\n"); }

        m->mothurOut("\nIt took " + toString(time(nullptr) - start) + " secs to process " + toString(numReads) + " sequences.\n");

        //output files created by command
		m->mothurOut("\nOutput File Names: \n");
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	} m->mothurOutEndLine();

        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "MakeContigsCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/
int MakeContigsCommand::createCountFile(string outputGroupFile, string resultFastafile) {
    try {
        
        CountTable ct; ct.createTable(groupMap);
        ct.printCompressedTable(outputGroupFile);
        
        if(allFiles){
            //run split.groups command
            //use unique.seqs to create new name and fastafile
            string inputString = "fasta=" + resultFastafile + ", count=" + outputGroupFile;
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
        m->errorOut(e, "MakeContigsCommand", "createCountFile");
        exit(1);
    }
}
//**********************************************************************************************************************
bool testGZReadable(vector<string>& fileInputs, vector<string>& indexInputs, bool& decompressionHelped, string format, MothurOut* m) {
    try {
        
        bool error = false; bool readable = true;
        decompressionHelped = false;
        
#ifdef USE_BOOST
        boost::iostreams::filtering_istream inFF, inRF;
        ifstream inForward, inReverse;
        string forwardFile = fileInputs[0];
        string reverseFile = fileInputs[1];
        
        Utils util;
        util.openInputFileBinary(forwardFile, inForward, inFF);
        util.openInputFileBinary(reverseFile, inReverse, inRF);
        
        FastqRead fread(inFF, error, format);
        FastqRead rread(inRF, error, format);
        inFF.pop(); inRF.pop();
        
        //error=true; to force test of decompression
        if (error) { //error reading fastq files, try unzipping
            
            string forwardOutput = util.getRootName(forwardFile) + "mothurTest_forward.fastq";
            string reverseOutput = util.getRootName(reverseFile) + "mothurTest_reverse.fastq";
            
            string unzipCommand = "gunzip < " + forwardFile + " > " + forwardOutput;
            system(unzipCommand.c_str());
            unzipCommand = "gunzip < " + reverseFile + " > " + reverseOutput;
            system(unzipCommand.c_str());
            
            ifstream inForward1, inReverse1;
            util.openInputFile(forwardOutput, inForward1);
            util.openInputFile(reverseOutput, inReverse1);
            
            FastqRead fread(inForward1, error, format);
            FastqRead rread(inReverse1, error, format);
            
            if (!error) {
                m->mothurOut("[WARNING]: mothur is unable to read your compressed fastq files. Decompressing files and continuing to process.\n\n");
                fileInputs[0] = forwardOutput;
                fileInputs[1] = reverseOutput;
                
                if (indexInputs.size() != 0) {
                    if ((indexInputs[0] != "NONE") && (indexInputs[0] != "")){
                        string forwardIndex = util.getRootName(indexInputs[0]) + "mothurTest_forward_index.fastq";
                        string unzipCommand = "gunzip < " + indexInputs[0] + " > " + forwardIndex;
                        system(unzipCommand.c_str());
                        indexInputs[0] = forwardIndex;
                    }
                    if ((indexInputs[1] != "NONE") && (indexInputs[1] != "")) {
                        string reverseIndex = util.getRootName(indexInputs[1]) + "mothurTest_reverse_index.fastq";
                        unzipCommand = "gunzip < " + indexInputs[1] + " > " + reverseIndex;
                        system(unzipCommand.c_str());
                        indexInputs[1] = reverseIndex;
                    }
                }
                
                decompressionHelped = true;
            }
            else { readable = false; }
            
        }
#endif
        
        return readable;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "testGZReadable");
        exit(1);
    }
}
//**********************************************************************************************************************
unsigned long long MakeContigsCommand::processSingleFileOption(string& outFastaFile, string& outScrapFastaFile, string& outQualFile, string& outScrapQualFile, string& outMisMatchFile, string group) {
    try {
        
        unsigned long long numReads = 0;
        string inputFile = "";
        vector<string> fileInputs;
        vector<string> qualOrIndexInputs;
        delim = '>';
        map<string, string> variables;
        string thisOutputDir = outputdir;

        if (ffastafile != "") {
            inputFile = ffastafile;
            if (outputdir == "") {  thisOutputDir = util.hasPath(inputFile); }
            fileInputs.push_back(ffastafile); fileInputs.push_back(rfastafile);

            if (fqualfile != "") {
                qualOrIndexInputs.push_back(fqualfile); qualOrIndexInputs.push_back(rqualfile);
                variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fqualfile));
                variables["[tag]"] = "trim";
                outQualFile = getOutputFileName("qfile",variables);
                variables["[tag]"] = "scrap";
                outScrapQualFile = getOutputFileName("qfile",variables);
            }else {
                outQualFile = ""; outScrapQualFile = ""; makeQualFile = false;
            }

            variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(inputFile));
            delim = '>';
        }else { //ffastqfile
            inputFile = ffastqfile;
            if (outputdir == "") {  thisOutputDir = util.hasPath(inputFile); }
            variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(inputFile));
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

        bool allGZ = true;
#ifdef USE_BOOST
        bool allPlainTxt = true;
        if (util.isGZ(fileInputs[0])[1]) { allPlainTxt = false;  }
        else {   allGZ = false;  }
        if (util.isGZ(fileInputs[1])[1]) { allPlainTxt = false;  }
        else {   allGZ = false;  }
        if (qualOrIndexInputs.size() != 0) {
            if (qualOrIndexInputs[0] != "NONE") {
                if (util.isGZ(qualOrIndexInputs[0])[1]) { allPlainTxt = false;  }
                else {   allGZ = false;  }
            }
            if (qualOrIndexInputs[1] != "NONE") {
                if (util.isGZ(qualOrIndexInputs[1])[1]) { allPlainTxt = false;  }
                else {   allGZ = false;  }
            }
            if (!allGZ && !allPlainTxt) { //mixed bag of files, uh oh...
                m->mothurOut("[ERROR]: Your files must all be in compressed .gz form or all in plain text form.  Please correct. \n"); m->setControl_pressed(true);
            }
        }
#else
        allGZ = false;
#endif
        bool decompressionHelped = false;
        if (allGZ)      {
            gz = true;
            //test to make sure you can read the gz files
            bool readable = testGZReadable(fileInputs, qualOrIndexInputs, decompressionHelped, format, m);
            
            if (readable) {
                if (decompressionHelped) { gz = false; }
            }else {
                m->mothurOut("[ERROR]: Unable to read compressed .gz files, please decompress and run make.contigs again. \n"); m->setControl_pressed(true); return 0;
            }
        }
        else            { gz = false;   }

        variables["[tag]"] = "trim";
        outFastaFile = getOutputFileName("fasta",variables);
        variables["[tag]"] = "scrap";
        outScrapFastaFile = getOutputFileName("fasta",variables);
        variables["[tag]"] = "";
        outMisMatchFile = getOutputFileName("report",variables);

        vector<vector<string> > fastaFileNames, qualFileNames;
        map<string, string> uniqueFastaNames;// so we don't add the same groupfile multiple times
        createOligosGroup = false;
        map<int, oligosPair> pairedPrimers, rpairedPrimers, revpairedPrimers, pairedBarcodes, rpairedBarcodes, revpairedBarcodes;
        vector<string> barcodeNames, primerNames;

        if(oligosfile != "")                        {       createOligosGroup = getOligos(pairedPrimers, rpairedPrimers, revpairedPrimers, pairedBarcodes, rpairedBarcodes, revpairedBarcodes, barcodeNames, primerNames);    }

        //give group in file file precedence
        if (createFileGroup) {  createOligosGroup = false; }

        m->mothurOut("Making contigs...\n");
        numReads = createProcesses(fileInputs, qualOrIndexInputs, outFastaFile, outScrapFastaFile, outQualFile, outScrapQualFile, outMisMatchFile, fastaFileNames, qualFileNames, group, pairedPrimers, rpairedPrimers, revpairedPrimers, pairedBarcodes, rpairedBarcodes, revpairedBarcodes, barcodeNames, primerNames);
        
        if (decompressionHelped) { util.mothurRemove(fileInputs[0]); util.mothurRemove(fileInputs[1]); }
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); }  return 0; }
        
        
        if (file == "") {
            outputNames.push_back(outFastaFile); outputTypes["fasta"].push_back(outFastaFile); outputNames.push_back(outScrapFastaFile); outputTypes["fasta"].push_back(outScrapFastaFile);
            if (makeQualFile) {
                 outputNames.push_back(outQualFile); outputTypes["qfile"].push_back(outQualFile); outputNames.push_back(outScrapQualFile); outputTypes["qfile"].push_back(outScrapQualFile); }
            else {
                if (outQualFile != "")      { util.mothurRemove(outQualFile);       }
                if (outScrapQualFile != "") { util.mothurRemove(outScrapQualFile);  }
            }
            outputNames.push_back(outMisMatchFile); outputTypes["report"].push_back(outMisMatchFile);
        }
        m->mothurOut("Done.\n");

        return numReads;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "processSingleFileOption");
        exit(1);
    }
}
/**************************************************************************************************/
struct contigsData {
    MothurOut* m;
    Utils util;
    OutputWriter* trimFileName;
    OutputWriter* scrapFileName;
    OutputWriter* trimQFileName;
    OutputWriter* scrapQFileName;
    OutputWriter* misMatchesFile;
    string align, group, format;
    float match, misMatch, gapOpen, gapExtend;
    bool gz, reorient, trimOverlap, createGroupFromOligos, createGroupFromFilePairs, makeQualFile, screenSequences;
    char delim;
    int nameType, offByOneTrimLength, pdiffs, bdiffs, tdiffs, kmerSize, insert, deltaq, maxee, maxAmbig, maxHomoP, maxLength;
    vector<string> inputFiles, qualOrIndexFiles, outputNames;
    set<string> badNames;
    linePair linesInput;
    linePair linesInputReverse;
    linePair qlinesInput;
    linePair qlinesInputReverse;
    long long count;

    vector<string> primerNameVector;
    vector<string> barcodeNameVector;
    map<string, int> groupCounts;
    map<string, string> groupMap;
    map<int, oligosPair> pairedBarcodes, reorientedPairedBarcodes, reversedPairedBarcodes;
    map<int, oligosPair> pairedPrimers, reorientedPairedPrimers, reversedPairedPrimers;


    contigsData(){}
    ~contigsData(){}
    contigsData(OutputWriter* tn, OutputWriter* sn, OutputWriter* tqn, OutputWriter* sqn, OutputWriter* mmf) {
        trimFileName = tn;
        scrapFileName = sn;
        trimQFileName = tqn;
        scrapQFileName = sqn;
        misMatchesFile = mmf;
        m = MothurOut::getInstance();
        count = 0;
        makeQualFile = true;
        if (trimQFileName == nullptr) { makeQualFile = false; }
    }

    contigsData(OutputWriter* tn, OutputWriter* sn, OutputWriter* tqn, OutputWriter* sqn, OutputWriter* mmf, vector<string> ifn, vector<string> qif, linePair li, linePair lir, linePair qli, linePair qlir) {
        trimFileName = tn;
        scrapFileName = sn;
        trimQFileName = tqn;
        scrapQFileName = sqn;
        misMatchesFile = mmf;
        m = MothurOut::getInstance();
        inputFiles = ifn;
        qualOrIndexFiles = qif;
        linesInput = li;
        linesInputReverse = lir;
        qlinesInput = qli;
        qlinesInputReverse = qlir;
        count = 0;
        makeQualFile = true;
        if (trimQFileName == nullptr) { makeQualFile = false; }
    }
    void setVariables(bool isgz, char de, int nt, int offby, map<int, oligosPair> pbr, map<int, oligosPair> ppr, map<int, oligosPair> rpbr, map<int, oligosPair> rppr, map<int, oligosPair> repbr, map<int, oligosPair> reppr, vector<string> priNameVector, vector<string> barNameVector, bool ro, int pdf, int bdf, int tdf, string al, float ma, float misMa, float gapO, float gapE, int thr, int delt, double maxe, int km, string form, bool to, bool cfg, bool cgff, string gp, bool screen, int maxH, int maxL, int maxAm) {
        gz = isgz;
        delim = de;
        nameType = nt;
        offByOneTrimLength = offby;
        pairedPrimers = ppr;
        pairedBarcodes = pbr;
        reorientedPairedPrimers = rppr;
        reorientedPairedBarcodes = rpbr;
        reversedPairedPrimers = reppr;
        reversedPairedBarcodes = repbr;
        primerNameVector = priNameVector;
        barcodeNameVector = barNameVector;
        group = gp;
        createGroupFromOligos = cfg;
        createGroupFromFilePairs = cgff;
        pdiffs = pdf;
        bdiffs = bdf;
        tdiffs = tdf;
        reorient = ro;
        match = ma;
        misMatch = misMa;
        gapOpen = gapO;
        gapExtend = gapE;
        insert = thr;
        kmerSize = km;
        align = al;
        deltaq = delt;
        maxee = maxe;
        format = form;
        trimOverlap = to;
        screenSequences = screen;
        maxHomoP = maxH;
        maxLength = maxL;
        maxAmbig = maxAm;
    }
    void copyVariables(contigsData* copy) {
        gz = copy->gz;
        delim = copy->delim;
        nameType = copy->nameType;
        offByOneTrimLength = copy->offByOneTrimLength;
        pairedPrimers = copy->pairedPrimers;
        pairedBarcodes = copy->pairedBarcodes;
        reorientedPairedPrimers = copy->reorientedPairedPrimers;
        reorientedPairedBarcodes = copy->reorientedPairedBarcodes;
        reversedPairedPrimers = copy->reversedPairedPrimers;
        reversedPairedBarcodes = copy->reversedPairedBarcodes;
        primerNameVector = copy->primerNameVector;
        barcodeNameVector = copy->barcodeNameVector;
        group = copy->group;
        createGroupFromOligos = copy->createGroupFromOligos;
        createGroupFromFilePairs = copy->createGroupFromFilePairs;
        pdiffs = copy->pdiffs;
        bdiffs = copy->bdiffs;
        tdiffs = copy->tdiffs;
        reorient = copy->reorient;
        match = copy->match;
        misMatch = copy->misMatch;
        gapOpen = copy->gapOpen;
        gapExtend = copy->gapExtend;
        insert = copy->insert;
        kmerSize = copy->kmerSize;
        align = copy->align;
        deltaq = copy->deltaq;
        maxee = copy->maxee;
        format = copy->format;
        trimOverlap = copy->trimOverlap;
        screenSequences = copy->screenSequences;
        maxHomoP = copy->maxHomoP;
        maxLength = copy->maxLength;
        maxAmbig = copy->maxAmbig;
    }
};
/**************************************************************************************************/
struct groupContigsData {
    MothurOut* m;
    Utils util;
    int start, end;
    vector< vector<string> > fileInputs;
    set<string> badNames;
    vector<string> file2Groups;
    contigsData* bundle;
    long long count;

    groupContigsData() = default;
    groupContigsData(vector< vector<string> > fi, int s, int e, contigsData* cd, vector<string> f2g) {
        fileInputs = fi;
        start = s;
        end = e;
        bundle = cd;
        file2Groups = f2g;
        count = 0;
        m = MothurOut::getInstance();
    }
    ~groupContigsData() { delete bundle; }
};
/**************************************************************************************************/

int setNameType(string forward, string reverse, int& offByOneTrimLength) {
    MothurOut* m = MothurOut::getInstance();
    try {
        int type = 0;
        Utils util;

        if (forward == reverse) {  type = perfectMatch;  }
        else {
            int pos = forward.find_last_of('#');
            string tempForward = forward;
            if (pos != string::npos) {  tempForward = forward.substr(0, pos);   }

            int pos2 = reverse.find_last_of('#');
            string tempReverse = reverse;
            if (pos2 != string::npos) {  tempReverse = reverse.substr(0, pos2);   }

            if (tempForward == tempReverse) { type = poundMatch;    }
            else {
                char delim = ':';
                if (m->getChangedSeqNames()) { delim = '_'; }
                vector<char> delims; delims.push_back(delim); delims.push_back('/');

                for (int j = 0; j < delims.size(); j++) {
                    delim = delims[j];
                    int pos = forward.find_last_of(delim);
                    string tempForward = forward;
                    string tempForwardEnd = forward;
                    if (pos != string::npos) {
                        tempForwardEnd = forward.substr(pos+1);
                    }

                    int pos2 = reverse.find_last_of(delim);
                    string tempReverse = reverse;
                    string tempReverseEnd = reverse;
                    if (pos2 != string::npos) {
                        tempReverseEnd = reverse.substr(pos2+1);
                    }

                    if (tempForwardEnd != tempReverseEnd) {
                        if ((util.isAllAlphaNumerics(tempForwardEnd)) && (util.isAllAlphaNumerics(tempReverseEnd))) {
                            //check for off by one on rest of name
                            if (tempForward.length() == tempReverse.length()) {
                                int numDiffs = 0;
                                char forwardDiff = ' '; char reverseDiff = ' '; int spot = 0;
                                for (int i = 0; i < tempForward.length(); i++) {
                                    if (tempForward[i] != tempReverse[i]) {
                                        numDiffs++;
                                        forwardDiff = tempForward[i];
                                        reverseDiff = tempReverse[i];
                                        spot = i;
                                    }
                                }
                                if (numDiffs == 1) {
                                    if ((forwardDiff == '1') && (reverseDiff == '2')) { type = offByOne; offByOneTrimLength = tempForward.length()-spot+1; }
                                }
                            }
                        }
                    }
                }
            }
        }

        return type;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "setNameType");
        exit(1);
    }
}
/**************************************************************************************************/

int setNameType(string forwardFile, string reverseFile, char delim, int& offByOneTrimLength, bool gz, string format) {
    MothurOut* m = MothurOut::getInstance();
    try {
        int type = 0; bool error = false;
        string forward = ""; string reverse = "";
        Utils util;

        ifstream inForward, inReverse;
#ifdef USE_BOOST
        boost::iostreams::filtering_istream inFF, inRF;
#endif
        if (!gz) { //plain text files
            util.openInputFile(forwardFile, inForward);
            util.openInputFile(reverseFile, inReverse);

            if (delim == '>') {
                Sequence fread(inForward);
                forward = fread.getName();
                Sequence rread(inReverse);
                reverse = rread.getName();
            }else {
                FastqRead fread(inForward, error, format);
                forward = fread.getName();
                FastqRead rread(inReverse, error, format);
                reverse = rread.getName();
            }
            inForward.close(); inReverse.close();
        }else { //compressed files
#ifdef USE_BOOST
            util.openInputFileBinary(forwardFile, inForward, inFF);
            util.openInputFileBinary(reverseFile, inReverse, inRF);

            if (delim == '>') {
                Sequence fread(inFF);
                forward = fread.getName();
                Sequence rread(inRF);
                reverse = rread.getName();
            }else {
                FastqRead fread(inFF, error, format);
                forward = fread.getName();
                FastqRead rread(inRF, error, format);
                reverse = rread.getName();
            }
            inFF.pop(); inRF.pop();
#endif
        }

        type = setNameType(forward, reverse, offByOneTrimLength);

        return type;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "setNameType");
        exit(1);
    }
}

//**********************************************************************************************************************
unsigned long long MakeContigsCommand::processMultipleFileOption(string& compositeFastaFile, string& compositeMisMatchFile) {
    try {
        //read file
        vector<string> file2Group;
        vector< vector<string> > fileInputs = readFileNames(file, file2Group);  if (m->getControl_pressed()) { return 0; }

        unsigned long long numReads = 0;

        map<string, string> cvars;
        string compOutputDir = outputdir;
        if (outputdir == "") { compOutputDir = util.hasPath(file); }
        cvars["[filename]"] = compOutputDir + util.getRootName(util.getSimpleName(file));
        cvars["[tag]"] = "trim";
        compositeFastaFile = getOutputFileName("fasta",cvars);
        cvars["[tag]"] = "scrap";
        string compositeScrapFastaFile = getOutputFileName("fasta",cvars);
        cvars["[tag]"] = "trim";
        string compositeQualFile = getOutputFileName("qfile",cvars);
        cvars["[tag]"] = "scrap";
        string compositeScrapQualFile = getOutputFileName("qfile",cvars);
        
        cvars["[tag]"] = "";
        compositeMisMatchFile = getOutputFileName("report",cvars);

        ofstream outCTFasta, outCTQual, outCSFasta, outCSQual, outCMisMatch;
        util.openOutputFile(compositeFastaFile, outCTFasta); outCTFasta.close(); outputNames.push_back(compositeFastaFile); outputTypes["fasta"].push_back(compositeFastaFile);
        
        util.openOutputFile(compositeScrapQualFile, outCSQual); outCSQual.close();
        util.openOutputFile(compositeQualFile, outCTQual); outCTQual.close();
        if (makeQualFile) {
            outputNames.push_back(compositeScrapQualFile); outputTypes["qfile"].push_back(compositeScrapQualFile);
            outputNames.push_back(compositeQualFile); outputTypes["qfile"].push_back(compositeQualFile);
        }
        util.openOutputFile(compositeScrapFastaFile, outCSFasta); outCSFasta.close(); outputNames.push_back(compositeScrapFastaFile); outputTypes["fasta"].push_back(compositeScrapFastaFile);
        
        util.openOutputFile(compositeMisMatchFile, outCMisMatch); outCMisMatch.close(); outputNames.push_back(compositeMisMatchFile); outputTypes["report"].push_back(compositeMisMatchFile);
        
        if (gz) {
            numReads = createProcessesGroups(fileInputs, compositeFastaFile, compositeScrapFastaFile, compositeQualFile, compositeScrapQualFile, compositeMisMatchFile, file2Group);
        }else {
            for (int l = 0; l < fileInputs.size(); l++) {
                if (m->getControl_pressed()) { break; }

                int startTime = time(nullptr);

                m->mothurOut("\n>>>>>\tProcessing file pair " + fileInputs[l][0] + " - " + fileInputs[l][1] + " (files " + toString(l+1) + " of " + toString(fileInputs.size()) + ")\t<<<<<\n");

                ffastqfile = fileInputs[l][0]; rfastqfile = fileInputs[l][1]; findexfile = fileInputs[l][2]; rindexfile = fileInputs[l][3];

                //run file as if it was a single
                string outFastaFile, outScrapFastaFile, outQualFile, outScrapQualFile, outMisMatchFile;
                int thisNumReads = processSingleFileOption(outFastaFile, outScrapFastaFile, outQualFile, outScrapQualFile, outMisMatchFile, file2Group[l]);
                numReads += thisNumReads;

                util.appendFiles(outMisMatchFile, compositeMisMatchFile); util.mothurRemove(outMisMatchFile);
                util.appendFiles(outFastaFile, compositeFastaFile);  util.mothurRemove(outFastaFile);
                util.appendFiles(outScrapFastaFile, compositeScrapFastaFile); util.mothurRemove(outScrapFastaFile);
                if (makeQualFile) {
                    util.appendFiles(outQualFile, compositeQualFile);
                    util.appendFiles(outScrapQualFile, compositeScrapQualFile);
                }
                util.mothurRemove(outQualFile); util.mothurRemove(outScrapQualFile);
                m->mothurOut("\nIt took " + toString(time(nullptr) - startTime) + " secs to assemble " + toString(thisNumReads) + " reads.\n\n");
            }
        }
        if (!makeQualFile) { util.mothurRemove(compositeQualFile); util.mothurRemove(compositeScrapQualFile); }
        
        return numReads;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "processMultipleFileOption");
        exit(1);
    }
}
//***************************************************************************************************************
/**
 * checks for minor diffs @MS7_15058:1:1101:11899:1633#8/1 @MS7_15058:1:1101:11899:1633#8/2 should match
 */
bool fixName(string& forward, int nameType, int offByOneTrimLength){
    try {
        bool match = false;

        if (nameType == poundMatch) {
            match = true;
            int pos = forward.find_last_of('#');
            if (pos != string::npos) {  forward = forward.substr(0, pos);   }
        }else if (nameType == perfectMatch) { match = true; }
        else if (nameType == offByOne) {
            match = true;
            forward = forward.substr(0, (forward.length()-offByOneTrimLength));
        }

        return match;
    }
    catch(exception& e) {
        MothurOut* m; m = MothurOut::getInstance();
        m->errorOut(e, "MakeContigsCommand", "fixName");
        exit(1);
    }
}
//***************************************************************************************************************
/**
 * checks for minor diffs @MS7_15058:1:1101:11899:1633#8/1 @MS7_15058:1:1101:11899:1633#8/2 should match
 */
bool checkName(FastqRead& forward, FastqRead& reverse, int nameType, int offByOneTrimLength){
    try {
        bool match = false;

        string forwardName = forward.getName();
        string reverseName = reverse.getName();

        if (nameType == poundMatch) {
            match = true;

            int pos = forwardName.find_last_of('#');
            if (pos != string::npos) {  forwardName = forwardName.substr(0, pos);   }

            int pos2 = reverseName.find_last_of('#');
            if (pos2 != string::npos) {  reverseName = reverseName.substr(0, pos2);   }

            if (forwardName == reverseName) {
                forward.setName(forwardName);
                reverse.setName(reverseName);
            }else{
                match = false;
            }

        }else if (nameType == perfectMatch) {
            if (forwardName == reverseName) { match = true; }
        }
        else if (nameType == offByOne) {

            match = true;

            reverseName = reverseName.substr(0, (reverseName.length()-offByOneTrimLength));
            forwardName = forwardName.substr(0, (forwardName.length()-offByOneTrimLength));

            if (forwardName == reverseName) {
                forward.setName(forwardName);
                reverse.setName(reverseName);
            }else{
                match = false;
            }
        }
        return match;
    }
    catch(exception& e) {
        MothurOut* m; m = MothurOut::getInstance();
        m->errorOut(e, "MakeContigsCommand", "ckeckName");
        exit(1);
    }
}
//***************************************************************************************************************
/**
 * checks for minor diffs @MS7_15058:1:1101:11899:1633#8/1 @MS7_15058:1:1101:11899:1633#8/2 should match
 */
bool checkName(Sequence& forward, Sequence& reverse, int nameType, int offByOneTrimLength){
    try {
        bool match = false;
        string forwardName = forward.getName();
        string reverseName = reverse.getName();

        if (nameType == poundMatch) {
            match = true;

            int pos = forwardName.find_last_of('#');
            if (pos != string::npos) {  forwardName = forwardName.substr(0, pos);   }

            int pos2 = reverseName.find_last_of('#');
            if (pos2 != string::npos) {  reverseName = reverseName.substr(0, pos2);   }

            if (forwardName == reverseName) {
                forward.setName(forwardName);
                reverse.setName(reverseName);
            }else{
                match = false;
            }
        }else if (nameType == perfectMatch) { if (forwardName == reverseName) { match = true; } }
        else if (nameType == offByOne) {

            match = true;

            reverseName = reverseName.substr(0, (reverseName.length()-offByOneTrimLength));
            forwardName = forwardName.substr(0, (forwardName.length()-offByOneTrimLength));

            if (forwardName == reverseName) {
                forward.setName(forwardName);
                reverse.setName(reverseName);
            }else{
                match = false;
            }
        }

        return match;
    }
    catch(exception& e) {
        MothurOut* m; m = MothurOut::getInstance();
        m->errorOut(e, "MakeContigsCommand", "ckeckName");
        exit(1);
    }
}
//***************************************************************************************************************
/**
 * checks for minor diffs @MS7_15058:1:1101:11899:1633#8/1 @MS7_15058:1:1101:11899:1633#8/2 should match
 */
bool checkName(QualityScores& forward, QualityScores& reverse, int nameType, int offByOneTrimLength){
    try {
        bool match = false;
        string forwardName = forward.getName();
        string reverseName = reverse.getName();

        if (nameType == poundMatch) {
            match = true;

            int pos = forwardName.find_last_of('#');
            if (pos != string::npos) {  forwardName = forwardName.substr(0, pos);   }

            int pos2 = reverseName.find_last_of('#');
            if (pos2 != string::npos) {  reverseName = reverseName.substr(0, pos2);   }

            if (forwardName == reverseName) {
                forward.setName(forwardName);
                reverse.setName(reverseName);
            }else{
                match = false;
            }

        }else if (nameType == perfectMatch) { if (forwardName == reverseName) { match = true; } }
        else if (nameType == offByOne) {

            match = true;

            reverseName = reverseName.substr(0, (reverseName.length()-offByOneTrimLength));
            forwardName = forwardName.substr(0, (forwardName.length()-offByOneTrimLength));

            if (forwardName == reverseName) {
                forward.setName(forwardName);
                reverse.setName(reverseName);
            }else{
                match = false;
            }
        }


        return match;

    }
    catch(exception& e) {
        MothurOut* m; m = MothurOut::getInstance();
        m->errorOut(e, "MakeContigsCommand", "ckeckName");
        exit(1);
    }
}

//***************************************************************************************************************
/**
 * checks for minor diffs @MS7_15058:1:1101:11899:1633#8/1 @MS7_15058:1:1101:11899:1633#8/2 should match
 */
bool checkName(Sequence& forward, QualityScores& reverse, int nameType, int offByOneTrimLength){
    try {
        bool match = false;
        string forwardName = forward.getName();
        string reverseName = reverse.getName();

        if (nameType == poundMatch) {
            match = true;

            string forwardName = forward.getName();
            string reverseName = reverse.getName();

            int pos = forwardName.find_last_of('#');
            if (pos != string::npos) {  forwardName = forwardName.substr(0, pos);   }

            int pos2 = reverseName.find_last_of('#');
            if (pos2 != string::npos) {  reverseName = reverseName.substr(0, pos2);   }

            if (forwardName == reverseName) {
                forward.setName(forwardName);
                reverse.setName(reverseName);
            }else{
                match = false;
            }

        }else if (nameType == perfectMatch) { if (forwardName == reverseName) { match = true; } }
        else if (nameType == offByOne) {

            match = true;

            reverseName = reverseName.substr(0, (reverseName.length()-offByOneTrimLength));
            forwardName = forwardName.substr(0, (forwardName.length()-offByOneTrimLength));

            if (forwardName == reverseName) {
                forward.setName(forwardName);
                reverse.setName(reverseName);
            }else{
                match = false;
            }
        }

        return match;
    }
    catch(exception& e) {
        MothurOut* m; m = MothurOut::getInstance();
        m->errorOut(e, "MakeContigsCommand", "ckeckName");
        exit(1);
    }
}

/**************************************************************************************************/

//vector<int> contigScores = assembleFragments(qual_match_simple_bayesian, qual_mismatch_simple_bayesian, fSeq, rSeq, alignment, contig);
vector<int> assembleFragments(vector< vector<double> >&qual_match_simple_bayesian, vector< vector<double> >& qual_mismatch_simple_bayesian, Sequence& fSeq, Sequence& rSeq, vector<int> scores1, vector<int> scores2, bool hasQuality, Alignment*& alignment, string& contig, string& trashCode, int& oend, int& oStart, int& numMismatches, int insert, int deltaq, bool trimOverlap) {
    MothurOut* m; m = MothurOut::getInstance();
    try {
        vector<int> contigScores;

        //flip the reverse reads
        rSeq.reverseComplement();
        
        //pairwise align
        alignment->align(fSeq.getUnaligned(), rSeq.getUnaligned(), true);
        map<int, int> ABaseMap = alignment->getSeqAAlnBaseMap();
        map<int, int> BBaseMap = alignment->getSeqBAlnBaseMap();
        fSeq.setAligned(alignment->getSeqAAln());
        rSeq.setAligned(alignment->getSeqBAln());
        int length = fSeq.getAligned().length();

        //traverse alignments merging into one contiguous seq
        string seq1 = fSeq.getAligned();
        string seq2 = rSeq.getAligned();
        
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

        int firstForward = 0; int seq2FirstForward = 0; int lastReverse = seq1.length(); int seq2lastReverse = seq2.length(); bool firstChooseSeq1 = false; bool lastChooseSeq1 = false;
        if (hasQuality) {
            for (int i = 0; i < seq1.length(); i++) { if ((seq1[i] != '.') && (seq1[i] != '-')) { if (scores1[ABaseMap[i]] == 2) { firstForward++; }else { break; } } }
            for (int i = 0; i < seq2.length(); i++) { if ((seq2[i] != '.') && (seq2[i] != '-')) { if (scores2[BBaseMap[i]] == 2) { seq2FirstForward++; }else { break; } } }
            if (seq2FirstForward > firstForward) { firstForward = seq2FirstForward; firstChooseSeq1 = true; }
            for (int i = seq1.length()-1; i >= 0; i--) { if ((seq1[i] != '.') && (seq1[i] != '-')) { if (scores1[ABaseMap[i]] == 2) { lastReverse--; }else { break; } } }
            for (int i = seq2.length()-1; i >= 0; i--) { if ((seq2[i] != '.') && (seq2[i] != '-')) { if (scores2[BBaseMap[i]] == 2) { seq2lastReverse--; }else { break; } } }
            if (lastReverse > seq2lastReverse) { lastReverse = seq2lastReverse; lastChooseSeq1 = true; }
        }

        for (int i = overlapStart; i < overlapEnd; i++) {
            if (seq1[i] == seq2[i]) {
                contig += seq1[i];
                if (hasQuality) {
                    contigScores.push_back(qual_match_simple_bayesian[PHREDCLAMP(scores1[ABaseMap[i]])][PHREDCLAMP(scores2[BBaseMap[i]])]);
                }
            }else if (((seq1[i] == '.') || (seq1[i] == '-')) && ((seq2[i] != '-') && (seq2[i] != '.'))) { //seq1 is a gap and seq2 is a base, choose seq2, unless quality score for base is below insert. In that case eliminate base
                if (hasQuality) {
                    if (scores2[BBaseMap[i]] <= insert) { } //
                    else {
                        contig += seq2[i];
                        contigScores.push_back(scores2[BBaseMap[i]]);
                    }
                } else {  contig += seq2[i]; } //with no quality info, then we keep it?
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
                            contigScores.push_back(qual_mismatch_simple_bayesian[PHREDCLAMP(scores1[ABaseMap[i]])][PHREDCLAMP(scores2[BBaseMap[i]])]);
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

        if (trimOverlap) {
            contig = contig.substr(overlapStart, oend-oStart);
            if (hasQuality) {
                vector<int> newContigScores;
                for (int i = overlapStart; i < oend; i++)  { newContigScores.push_back(contigScores[i]);  }
                contigScores = newContigScores;
            }
        }

        if (contig == "") { trashCode += "l"; contig = "NNNN"; contigScores.push_back(2); contigScores.push_back(2); contigScores.push_back(2); contigScores.push_back(2); }

        return contigScores;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "assembleFragments");
        exit(1);
    }
}
/**************************************************************************************************/
void screenSequence(string& contig, string& trashCode, contigsData* params) {
    try {
        //if you failed before screening, don't bother screening
        if (trashCode.length() != 0) { return; }
        else {
            bool goodSeq = true;        //    innocent until proven guilty
            
            Sequence currSeq("dummy", contig);
            
            if(params->maxAmbig != -1 && params->maxAmbig <    currSeq.getAmbigBases())        {    goodSeq = false;    trashCode += "ambig|";  }
            if(params->maxHomoP != -1 && params->maxHomoP < currSeq.getLongHomoPolymer())    {    goodSeq = false;    trashCode += "homop|";  }
            if(params->maxLength != -1 && params->maxLength < currSeq.getNumBases())        {    goodSeq = false;    trashCode += ">length|";}
        }
        
    }
    catch(exception& e) {
        params->m->errorOut(e, "MakeContigsCommand", "screenSequence");
        exit(1);
    }
}
/**************************************************************************************************/
#ifdef USE_BOOST
//ignore = read(fSeq, rSeq, fQual, rQual, savedFQual, savedRQual, findexBarcode, rindexBarcode, delim,  inFF, inRF, inFQ, inRQ);
bool read(Sequence& fSeq, Sequence& rSeq, QualityScores*& fQual, QualityScores*& rQual, Sequence& findexBarcode, Sequence& rindexBarcode, char delim, boost::iostreams::filtering_istream& inFF, boost::iostreams::filtering_istream& inRF, boost::iostreams::filtering_istream& inFQ, boost::iostreams::filtering_istream& inRQ, string thisfqualindexfile, string thisrqualindexfile, string format, int nameType, int offByOneTrimLength, MothurOut* m) {
    try {
        bool ignore = false;
        Utils util;
        if (delim == '@') { //fastq files
            bool tignore = false;
            FastqRead fread(inFF, tignore, format);  gobble(inFF);
            FastqRead rread(inRF, ignore, format); gobble(inRF);
            
            if (!checkName(fread, rread, nameType, offByOneTrimLength)) {
                FastqRead f2read(inFF, tignore, format);
                if (!checkName(f2read, rread, nameType, offByOneTrimLength)) {
                    FastqRead r2read(inRF, ignore, format);
                    if (!checkName(fread, r2read, nameType, offByOneTrimLength)) {
                        m->mothurOut("[WARNING]: name mismatch in forward and reverse fastq file. Ignoring, " + fread.getName() + ".\n"); ignore = true;
                    }else { rread = r2read; }
                }else { fread = f2read; }
            }
            if (tignore) { ignore=true; }
            fSeq.setName(fread.getName()); fSeq.setAligned(fread.getSeq());
            rSeq.setName(rread.getName()); rSeq.setAligned(rread.getSeq());
            
            fQual = new QualityScores(fread.getName(), fread.getScores());
            rQual = new QualityScores(rread.getName(), rread.getScores());
            if (thisfqualindexfile != "") { //forward index file
                FastqRead firead(inFQ, tignore, format);
                if (tignore) { ignore=true; }
                findexBarcode.setAligned(firead.getSeq());
                if (!checkName(fread, firead, nameType, offByOneTrimLength)) {
                    FastqRead f2iread(inFQ, tignore, format);
                    if (tignore) { ignore=true; }
                    if (!checkName(fread, f2iread, nameType, offByOneTrimLength)) {
                        m->mothurOut("[WARNING]: name mismatch in forward index file. Ignoring, " + fread.getName() + ".\n"); ignore = true;
                    }else { firead = f2iread; findexBarcode.setAligned(firead.getSeq()); }
                }
            }
            if (thisrqualindexfile != "") { //reverse index file
                FastqRead riread(inRQ, tignore, format);
                if (tignore) { ignore=true; }
                rindexBarcode.setAligned(riread.getSeq());
                if (!checkName(fread, riread, nameType, offByOneTrimLength)) {
                    FastqRead r2iread(inRQ, tignore, format); gobble(inRQ);
                    if (tignore) { ignore=true; }
                    if (!checkName(fread, r2iread, nameType, offByOneTrimLength)) {
                        m->mothurOut("[WARNING]: name mismatch in reverse index file. Ignoring, " + fread.getName() + ".\n"); ignore = true;
                    }else { riread = r2iread; rindexBarcode.setAligned(riread.getSeq()); }
                }
            }
        }else { //reading fasta and maybe qual
            Sequence tfSeq(inFF);
            Sequence trSeq(inRF);
            if (!checkName(tfSeq, trSeq, nameType, offByOneTrimLength)) {
                Sequence t2fSeq(inFF);
                if (!checkName(t2fSeq, trSeq, nameType, offByOneTrimLength)) {
                    Sequence t2rSeq(inRF);
                    if (!checkName(tfSeq, t2rSeq, nameType, offByOneTrimLength)) {
                        m->mothurOut("[WARNING]: name mismatch in forward and reverse fasta file. Ignoring, " + tfSeq.getName() + ".\n"); ignore = true;
                    }else { trSeq = t2fSeq; }
                }else { tfSeq = t2fSeq; }
            }
            fSeq.setName(tfSeq.getName()); fSeq.setAligned(tfSeq.getAligned());
            rSeq.setName(trSeq.getName()); rSeq.setAligned(trSeq.getAligned());
            if (thisfqualindexfile != "") {
                fQual = new QualityScores(inFQ); gobble(inFQ);
                rQual = new QualityScores(inRQ); gobble(inRQ);
                if (!checkName(*fQual, *rQual, nameType, offByOneTrimLength)) {
                    m->mothurOut("[WARNING]: name mismatch in forward and reverse qual file. Ignoring, " + fQual->getName() + ".\n"); ignore = true;
                }
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
#endif
/**************************************************************************************************/

bool read(Sequence& fSeq, Sequence& rSeq, QualityScores*& fQual, QualityScores*& rQual,Sequence& findexBarcode, Sequence& rindexBarcode, char delim, ifstream& inFFasta, ifstream& inRFasta, ifstream& inFQualIndex, ifstream& inRQualIndex, string thisfqualindexfile, string thisrqualindexfile, string format, int nameType, int offByOneTrimLength, MothurOut* m) {
    try {
        bool ignore = false;
        Utils util;
        if (delim == '@') { //fastq files
            bool tignore;
            FastqRead fread(inFFasta, tignore, format); gobble(inFFasta);
            FastqRead rread(inRFasta, ignore, format); gobble(inRFasta);
            
            if (!checkName(fread, rread, nameType, offByOneTrimLength)) {
                FastqRead f2read(inFFasta, tignore, format); gobble(inFFasta);
                if (!checkName(f2read, rread, nameType, offByOneTrimLength)) {
                    FastqRead r2read(inRFasta, ignore, format); gobble(inRFasta);
                    if (!checkName(fread, r2read, nameType, offByOneTrimLength)) {
                        m->mothurOut("[WARNING]: name mismatch in forward and reverse fastq file. Ignoring, " + fread.getName() + ".\n"); ignore = true;
                    }else { rread = r2read; }
                }else { fread = f2read; }
            }
            if (tignore) { ignore=true; }
            fSeq.setName(fread.getName()); fSeq.setAligned(fread.getSeq());
            rSeq.setName(rread.getName()); rSeq.setAligned(rread.getSeq());
            fQual = new QualityScores(fread.getName(), fread.getScores());
            rQual = new QualityScores(rread.getName(), rread.getScores());
            if (thisfqualindexfile != "") { //forward index file
                FastqRead firead(inFQualIndex, tignore, format); gobble(inFQualIndex);
                if (tignore) { ignore=true; }
                findexBarcode.setAligned(firead.getSeq());
                if (!checkName(fread, firead, nameType, offByOneTrimLength)) {
                    FastqRead f2iread(inFQualIndex, tignore, format); gobble(inFQualIndex);
                    if (tignore) { ignore=true; }
                    if (!checkName(fread, f2iread, nameType, offByOneTrimLength)) {
                        m->mothurOut("[WARNING]: name mismatch in forward index file. Ignoring, " + fread.getName() + ".\n"); ignore = true;
                    }else { firead = f2iread; findexBarcode.setAligned(firead.getSeq()); }
                }
            }
            if (thisrqualindexfile != "") { //reverse index file
                FastqRead riread(inRQualIndex, tignore, format); gobble(inRQualIndex);
                if (tignore) { ignore=true; }
                rindexBarcode.setAligned(riread.getSeq());
                if (!checkName(fread, riread, nameType, offByOneTrimLength)) {
                    FastqRead r2iread(inRQualIndex, tignore, format); gobble(inRQualIndex);
                    if (tignore) { ignore=true; }
                    if (!checkName(fread, r2iread, nameType, offByOneTrimLength)) {
                        m->mothurOut("[WARNING]: name mismatch in reverse index file. Ignoring, " + fread.getName() + ".\n"); ignore = true;
                    }else { riread = r2iread; rindexBarcode.setAligned(riread.getSeq()); }
                }
            }
        }else { //reading fasta and maybe qual
            Sequence tfSeq(inFFasta); gobble(inFFasta);
            Sequence trSeq(inRFasta); gobble(inRFasta);
            if (!checkName(tfSeq, trSeq, nameType, offByOneTrimLength)) {
                Sequence t2fSeq(inFFasta); gobble(inFFasta);
                if (!checkName(t2fSeq, trSeq, nameType, offByOneTrimLength)) {
                    Sequence t2rSeq(inRFasta); gobble(inRFasta);
                    if (!checkName(tfSeq, t2rSeq, nameType, offByOneTrimLength)) {
                        m->mothurOut("[WARNING]: name mismatch in forward and reverse fasta file. Ignoring, " + tfSeq.getName() + ".\n"); ignore = true;
                    }else { trSeq = t2fSeq; }
                }else { tfSeq = t2fSeq; }
            }
            fSeq.setName(tfSeq.getName()); fSeq.setAligned(tfSeq.getAligned());
            rSeq.setName(trSeq.getName()); rSeq.setAligned(trSeq.getAligned());
            if (thisfqualindexfile != "") {
                fQual = new QualityScores(inFQualIndex); gobble(inFQualIndex);
                rQual = new QualityScores(inRQualIndex); gobble(inRQualIndex);
                if (!checkName(*fQual, *rQual, nameType, offByOneTrimLength)) {
                    m->mothurOut("[WARNING]: name mismatch in forward and reverse qual file. Ignoring, " + fQual->getName() + ".\n"); ignore = true;
                }
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

vector<int> trimBarCodesAndPrimers(Sequence& fSeq, Sequence& rSeq, QualityScores*& fQual, QualityScores*& rQual, Sequence& findexBarcode, Sequence& rindexBarcode, vector<TrimOligos*> trims, vector<string>& codes, int numBarcodes, int numPrimers, bool hasQuality, bool hasIndex, int pdiffs, int bdiffs, int tdiffs, MothurOut* m) {
    try {
        
        vector<int> oligosResults; oligosResults.resize(2, 0); //barcodeIndex, primerIndex
        codes.resize(2, "");
        
        for (int i = 0; i < trims.size(); i++) {
            Sequence savedFSeq(fSeq.getName(), fSeq.getAligned());
            Sequence savedRSeq(rSeq.getName(), rSeq.getAligned());
            Sequence savedFindex(findexBarcode.getName(), findexBarcode.getAligned());
            Sequence savedRIndex(rindexBarcode.getName(), rindexBarcode.getAligned());
            
            QualityScores* savedFQual = nullptr;
            QualityScores* savedRQual = nullptr;
            if (hasQuality) {
                savedFQual = new QualityScores(fQual->getName(), fQual->getScores());
                savedRQual = new QualityScores(rQual->getName(), rQual->getScores());
            }
            
            string trashCode = "";
            string commentString = "";
            int currentSeqsDiffs = 0;
            int barcodeIndex = 0;
            int primerIndex = 0;
            
            if(numBarcodes != 0){
                vector<int> results;
                if (hasQuality) {
                    if (hasIndex)   {  results = trims[i]->stripBarcode(savedFindex, savedRIndex, *savedFQual, *savedRQual, barcodeIndex);    }
                    else            {  results = trims[i]->stripBarcode(savedFSeq, savedRSeq, *savedFQual, *savedRQual, barcodeIndex);                            }
                }else {
                    results = trims[i]->stripBarcode(savedFSeq, savedRSeq, barcodeIndex);
                }
                int success = results[0] + results[2];
                commentString += "fbdiffs=" + toString(results[0]) + "(" + trims[i]->getCodeValue(results[1], bdiffs) + "), rbdiffs=" + toString(results[2]) + "(" + trims[i]->getCodeValue(results[3], bdiffs) + ") ";
                if(success > bdiffs)        {    trashCode += 'b';    }
                else{ currentSeqsDiffs += success;  }
            }
            
            if(numPrimers != 0){
                vector<int> results;
                if (hasQuality)     { results = trims[i]->stripForward(savedFSeq, savedRSeq, *savedFQual, *savedRQual, primerIndex);   }
                else                { results = trims[i]->stripForward(savedFSeq, savedRSeq, primerIndex);                   }
                int success = results[0] + results[2];
                commentString += "fpdiffs=" + toString(results[0]) + "(" + trims[i]->getCodeValue(results[1], pdiffs) + "), rpdiffs=" + toString(results[2]) + "(" + trims[i]->getCodeValue(results[3], pdiffs) + ") ";
                if(success > pdiffs)        {    trashCode += 'f';    }
                else{ currentSeqsDiffs += success;  }
            }
            
            if (currentSeqsDiffs > tdiffs)    {    trashCode += 't';   }
            
            if (trashCode == "") {
                oligosResults[0] = barcodeIndex;
                oligosResults[1] = primerIndex;
                
                codes[0] = "";  
                codes[1] = commentString;
               
                if (i > 0) { //checkOrient trimOligos - reoriented and reversed
                    savedFSeq.reverseComplement();
                    savedRSeq.reverseComplement();
                }
                fSeq.setAligned(savedFSeq.getUnaligned());
                rSeq.setAligned(savedRSeq.getUnaligned());
                if(hasQuality){
                    if (i > 0) { //checkOrient trimOligos - reoriented and reversed
                        savedFQual->flipQScores();
                        savedRQual->flipQScores();
                    }
                    fQual->setScores(savedFQual->getScores());
                    rQual->setScores(savedRQual->getScores());
                    delete savedRQual; delete savedFQual;
                }
                break;
            }else {
                if (codes[0] == "") { codes[0] = trashCode;                 }
                else                { codes[0] += "(" + trashCode + ")";    }
                
                codes[1] = commentString;
                
                if(hasQuality){  delete savedRQual; delete savedFQual;  }
            }
        }
        
        if (hasQuality) { rQual->flipQScores(); }
        
        return oligosResults;

    }catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "trimBarCodesAndPrimers");
        exit(1);
    }
}
//**********************************************************************************************************************
//vector<vector<string> > fastaFileNames, vector<vector<string> > qualFileNames, , string group
void driverContigs(contigsData* params){
    try {
        vector< vector<double> > qual_match_simple_bayesian;
        qual_match_simple_bayesian.resize(47);
        for (int i = 0; i < qual_match_simple_bayesian.size(); i++) { qual_match_simple_bayesian[i].resize(47);  }

        vector< vector<double> > qual_mismatch_simple_bayesian;
        qual_mismatch_simple_bayesian.resize(47);
        for (int i = 0; i < qual_mismatch_simple_bayesian.size(); i++) { qual_mismatch_simple_bayesian[i].resize(47);  }

        loadQmatchValues(qual_match_simple_bayesian, qual_mismatch_simple_bayesian);

        params->count = 0;
        string thisfqualindexfile, thisrqualindexfile, thisffastafile, thisrfastafile;
        thisfqualindexfile = ""; thisrqualindexfile = "";
        thisffastafile = params->inputFiles[0]; thisrfastafile = params->inputFiles[1];
        if (params->qualOrIndexFiles.size() != 0) {
            thisfqualindexfile = params->qualOrIndexFiles[0];
            thisrqualindexfile = params->qualOrIndexFiles[1];
        }

        if (params->m->getDebug()) {  params->m->mothurOut("[DEBUG]: ffasta = " + thisffastafile + ".\n[DEBUG]: rfasta = " + thisrfastafile + ".\n[DEBUG]: fqualindex = " + thisfqualindexfile + ".\n[DEBUG]: rqualindex = " + thisfqualindexfile + ".\n"); }

        ifstream inFFasta, inRFasta, inFQualIndex, inRQualIndex;
#ifdef USE_BOOST
        boost::iostreams::filtering_istream inFF, inRF, inFQ, inRQ;
#endif
        if (!params->gz) { //plain text files
            params->util.openInputFile(thisffastafile, inFFasta);
            params->util.openInputFile(thisrfastafile, inRFasta);

            inFFasta.seekg(params->linesInput.start);
            inRFasta.seekg(params->linesInputReverse.start);
        }else { //compressed files - no need to seekg because compressed files divide workload differently
#ifdef USE_BOOST
            params->util.openInputFileBinary(thisffastafile, inFFasta, inFF);
            params->util.openInputFileBinary(thisrfastafile, inRFasta, inRF);
#endif
        }

        ofstream outFasta, outMisMatch, outScrapFasta, outQual, outScrapQual;
        if (thisfqualindexfile != "") {
            if (thisfqualindexfile != "NONE") {
                if (!params->gz) { //plain text files
                    params->util.openInputFile(thisfqualindexfile, inFQualIndex);
                    inFQualIndex.seekg(params->qlinesInput.start);
                }else {
#ifdef USE_BOOST
                    params->util.openInputFileBinary(thisfqualindexfile, inFQualIndex, inFQ);
#endif
                } //compressed files - no need to seekg because compressed files divide workload differently
            }
            else {  thisfqualindexfile = ""; }
            if (thisrqualindexfile != "NONE") {
                if (!params->gz) { //plain text files
                    params->util.openInputFile(thisrqualindexfile, inRQualIndex);
                    inRQualIndex.seekg(params->qlinesInputReverse.start);
                }else {
#ifdef USE_BOOST
                    params->util.openInputFileBinary(thisrqualindexfile, inRQualIndex, inRQ);
#endif
                } //compressed files - no need to seekg because compressed files divide workload differently
            }
            else { thisrqualindexfile = ""; }
        }

        bool hasQuality = false;
        bool hasIndex = false;
        if (params->delim == '@') { //fastq files so make an output quality
            hasQuality = true;
            if (thisfqualindexfile != "")           { if (thisfqualindexfile != "NONE") {  hasIndex = true; } }
            if (thisrqualindexfile != "")           { if (thisrqualindexfile != "NONE") {  hasIndex = true; } }
        }else if ((params->delim == '>') && (params->qualOrIndexFiles.size() != 0)) { hasQuality = true; }

        if (params->m->getDebug()) { if (hasQuality) { params->m->mothurOut("[DEBUG]: hasQuality = true\n");  } else { params->m->mothurOut("[DEBUG]: hasQuality = false\n"); } }

        int numPrimers = params->pairedPrimers.size();
        int numBarcodes = params->pairedBarcodes.size();
        
        vector<TrimOligos*> trims;
        if ((numPrimers != 0) || (numBarcodes != 0)) {
            //standard
            trims.push_back(new TrimOligos(params->pdiffs, params->bdiffs, 0, 0, params->pairedPrimers, params->pairedBarcodes, hasIndex));
    
            if (params->reorient) {
                //reoriented
                trims.push_back(new TrimOligos(params->pdiffs, params->bdiffs, 0, 0, params->reorientedPairedPrimers, params->reorientedPairedBarcodes, hasIndex));
                //reversed
                trims.push_back(new TrimOligos(params->pdiffs, params->bdiffs, 0, 0, params->reversedPairedPrimers, params->reversedPairedBarcodes, hasIndex));
            }
        }
        Alignment* alignment; int longestBase = 1000;
        if(params->align == "gotoh")			{	alignment = new GotohOverlap(params->gapOpen, params->gapExtend, params->match, params->misMatch, longestBase);			}
        else if(params->align == "needleman")	{	alignment = new NeedlemanOverlap(params->gapOpen, params->match, params->misMatch, longestBase);                        }
        else if(params->align == "kmer")        {   alignment = new KmerAlign(params->kmerSize);                                                                            }

        bool good = true;
        while (good) {

            if (params->m->getControl_pressed()) { break; }

            bool ignore = false;
            Sequence fSeq, rSeq;
            QualityScores* fQual = nullptr; QualityScores* rQual = nullptr;
            Sequence findexBarcode("findex", "NONE");  Sequence rindexBarcode("rindex", "NONE");

            //read from input files
            if (params->gz) {
#ifdef USE_BOOST
                ignore = read(fSeq, rSeq, fQual, rQual, findexBarcode, rindexBarcode, params->delim, inFF, inRF, inFQ, inRQ, thisfqualindexfile, thisrqualindexfile, params->format, params->nameType, params->offByOneTrimLength, params->m);
#endif
            }else    {
                ignore = read(fSeq, rSeq, fQual, rQual, findexBarcode, rindexBarcode, params->delim, inFFasta, inRFasta, inFQualIndex, inRQualIndex, thisfqualindexfile, thisrqualindexfile, params->format, params->nameType, params->offByOneTrimLength, params->m);
            }

            if (!ignore) {
                //remove primers and barcodes if neccessary
                vector<string> codes;
                vector<int> oligosResults = trimBarCodesAndPrimers(fSeq, rSeq, fQual, rQual, findexBarcode, rindexBarcode, trims, codes, numBarcodes, numPrimers, hasQuality, hasIndex, params->pdiffs, params->bdiffs, params->tdiffs, params->m);
                
                string trashCode = codes[0];
                string commentString = codes[1];
                int barcodeIndex = oligosResults[0];
                int primerIndex = oligosResults[1];
                
                //assemble reads
                string contig = "";
                int oend, oStart;
                int numMismatches = 0;
                vector<int> scores1, scores2;
                if(hasQuality){
                    scores1 = fQual->getScores(); scores2 = rQual->getScores();
                    delete fQual; delete rQual;
                }
                
                vector<int> contigScores = assembleFragments(qual_match_simple_bayesian, qual_mismatch_simple_bayesian, fSeq, rSeq, scores1, scores2, hasQuality, alignment, contig, trashCode, oend, oStart, numMismatches, params->insert, params->deltaq, params->trimOverlap);

								//Note that usearch/vsearch cap the maximum Q value at 41 - perhaps due to ascii
								//limits? we leave this value unbounded. if two sequences have a 40 then the
								//assembled quality score will be 85. If two 250 nt reads are all 40 and they
 								//perfectly match each other, then the difference in the number of expected errors
								//between using 85 and 41 all the way across will be 0.01986 - this is a "worst"
								//case scenario

								double expected_errors = 0;
								for(int i=0;i<contigScores.size();i++){
									expected_errors += convertQToProb(contigScores[i]);
								}

								if(expected_errors > params->maxee) { trashCode += 'e' ;}
                
                if (params->screenSequences) { screenSequence(contig, trashCode, params); }
                
                if(trashCode.length() == 0){
                    string thisGroup = params->group; //group from file file
                    if (params->createGroupFromOligos) { //overwrite file file group for oligos group
                        if(numBarcodes != 0){
                            thisGroup = params->barcodeNameVector[barcodeIndex];
                            if (numPrimers != 0) {
                                if (params->primerNameVector[primerIndex] != "") {
                                    if(thisGroup != "") { thisGroup += "." + params->primerNameVector[primerIndex]; }
                                    else                { thisGroup = params->primerNameVector[primerIndex];        }
                                }
                            }
                        }
                    }

                    int pos = thisGroup.find("ignore");
                    if (pos == string::npos) {
                        if (thisGroup != "") {
                            params->groupMap[fSeq.getName()] = thisGroup;

                            map<string, int>::iterator it = params->groupCounts.find(thisGroup);
                            if (it == params->groupCounts.end()) {	params->groupCounts[thisGroup] = 1; }
                            else { params->groupCounts[it->first] ++; }
                        }
                    }else { ignore = true; }

                    //print good stuff
                    if(!ignore){
                        //output
                        string output = ">" + fSeq.getName() + '\t' + "ee=" + toString(expected_errors) + '\t' + commentString + "\n" + contig + "\n";
                        params->trimFileName->write(output);
                        if (hasQuality && params->makeQualFile) {
                            output = ">" + fSeq.getName() + '\t' + "ee=" + toString(expected_errors) + '\t' + commentString +"\n";
                            for (int i = 0; i < contigScores.size(); i++) { output += toString(contigScores[i]) + " "; }  output += "\n";
                            params->trimQFileName->write(output);
                        }
                        int numNs = 0;
                        for (int i = 0; i < contig.length(); i++) { if (contig[i] == 'N') { numNs++; }  }
                        output = fSeq.getName() + '\t' + toString(contig.length()) + '\t' + toString(oend-oStart) + '\t' + toString(oStart) + '\t' + toString(oend) + '\t' + toString(numMismatches) + '\t' + toString(numNs) + '\t' + toString(expected_errors) + "\n";
                        params->misMatchesFile->write(output);
                    }
                }else{
                    params->badNames.insert(fSeq.getName());

                    string output = ">" + fSeq.getName() + " | " + trashCode + '\t' + "ee=" +  toString(expected_errors) + '\t' + commentString + "\n" + contig + "\n";
                    params->scrapFileName->write(output);

                    if (hasQuality && params->makeQualFile) {
                        output = ">" + fSeq.getName() + " | " + trashCode + '\t' + "ee=" + toString(expected_errors) + '\t' + commentString + "\n";
                        for (int i = 0; i < contigScores.size(); i++) { output += toString(contigScores[i]) + " "; }  output += "\n";
                        params->scrapQFileName->write(output);
                    }
                }
                if (params->m->getDebug()) { params->m->mothurOut("\n"); }
            }
            params->count++;

#if defined NON_WINDOWS
            if (!params->gz) {
                double pos = inFFasta.tellg();
                if (params->util.isEqual(pos,-1) || (pos >= params->linesInput.end)) { good = false; break; }
            }else {
    #ifdef USE_BOOST
                if (inFF.eof() || inRF.eof()) { good = false; break; }
    #endif
            }
#else
            if (!params->gz) {
                if (params->count >= params->linesInput.end) { good = false; break; }
            }else {
    #ifdef USE_BOOST
                if (inFF.eof() || inRF.eof()) { good = false; break; }
    #endif
            }
#endif

            //report progress
            if((params->count) % 1000 == 0){	params->m->mothurOutJustToScreen(toString(params->count)+"\n"); }
        }

        //report progress
        if((params->count) % 1000 != 0){	params->m->mothurOutJustToScreen(toString(params->count)+"\n"); }

        //close files
        inFFasta.close();
        inRFasta.close();
        if (params->gz) {
#ifdef USE_BOOST
            inFF.pop(); inRF.pop();
#endif
        }

        if (params->delim == '@') {
            if (thisfqualindexfile != "") { inFQualIndex.close();
                if (params->gz) {
#ifdef USE_BOOST
                    inFQ.pop();
#endif
                }
            }
            if (thisrqualindexfile != "") { inRQualIndex.close();
                if (params->gz) {
#ifdef USE_BOOST
                    inRQ.pop();
#endif
                }
            }
        }else{
            if (hasQuality) {
                inFQualIndex.close();
                inRQualIndex.close();
                if (params->gz) {
#ifdef USE_BOOST
                    inFQ.pop(); inRQ.pop();
#endif
                }
            }
        }

        //cleanup memory
        for (int i = 0; i < trims.size(); i++) {  delete trims[i]; }
        delete alignment;
        
    }
    catch(exception& e) {
        params->m->errorOut(e, "MakeContigsCommand", "driverContigs");
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
unsigned long long MakeContigsCommand::createProcesses(vector<string> fileInputs, vector<string> qualOrIndexFiles, string outputFasta, string outputScrapFasta, string outputQual, string outputScrapQual, string outputMisMatches, vector<vector<string> > fastaFileNames, vector<vector<string> > qualFileNames, string group, map<int, oligosPair>& pairedPrimers, map<int, oligosPair>& rpairedPrimers, map<int, oligosPair>& revpairedPrimers, map<int, oligosPair>& pairedBarcodes, map<int, oligosPair>& rpairedBarcodes, map<int, oligosPair>& revpairedBarcodes, vector<string>& barcodeNames, vector<string>& primerNames) {
    try {
        vector<linePair> lines;
        vector<linePair> qLines;

        if (gz)  {
            nameType = setNameType(fileInputs[0], fileInputs[1], delim, offByOneTrimLength,  gz, format);
            for (int i = 0; i < fileInputs.size(); i++) {
                //fake out lines - we are just going to check for end of file. Work is divided by number of files per processor.
                lines.push_back(linePair(0, 1000));
                qLines.push_back(linePair(0, 1000));
            }
            processors = fileInputs.size() / 2;
        }else        {
            //divides the files so that the processors can share the workload.
            setLines(fileInputs, qualOrIndexFiles, lines, qLines, delim);
        }

        bool hasQuality = false;
        if (delim == '@')                                           { hasQuality = true; }
        else if ((delim == '>') && (qualOrIndexFiles.size() != 0))  { hasQuality = true; }

        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<contigsData*> data;

        auto synchronizedOutputFastaTrimFile = std::make_shared<SynchronizedOutputFile>(outputFasta);
        auto synchronizedOutputFastaScrapFile = std::make_shared<SynchronizedOutputFile>(outputScrapFasta);
        auto synchronizedMisMatchFile = std::make_shared<SynchronizedOutputFile>(outputMisMatches);
        auto synchronizedOutputQTrimFile = std::make_shared<SynchronizedOutputFile>(outputQual);
        auto synchronizedOutputQScrapFile = std::make_shared<SynchronizedOutputFile>(outputScrapQual);
       

        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            OutputWriter* threadFastaTrimWriter = new OutputWriter(synchronizedOutputFastaTrimFile);
            OutputWriter* threadFastaScrapWriter = new OutputWriter(synchronizedOutputFastaScrapFile);
            OutputWriter* threadMismatchWriter = new OutputWriter(synchronizedMisMatchFile);
            OutputWriter* threadQTrimWriter = nullptr;
            OutputWriter* threadQScrapWriter = nullptr;
            if (makeQualFile) {
                threadQTrimWriter = new OutputWriter(synchronizedOutputQTrimFile);
                threadQScrapWriter = new OutputWriter(synchronizedOutputQScrapFile);
            }

            int spot = (i+1)*2;
            contigsData* dataBundle = new contigsData(threadFastaTrimWriter, threadFastaScrapWriter, threadQTrimWriter, threadQScrapWriter, threadMismatchWriter, fileInputs, qualOrIndexFiles, lines[spot], lines[spot+1], qLines[spot], qLines[spot+1]);
            dataBundle->setVariables(gz, delim, nameType, offByOneTrimLength, pairedBarcodes, pairedPrimers, rpairedBarcodes, rpairedPrimers, revpairedBarcodes, revpairedPrimers, primerNames, barcodeNames, reorient, pdiffs, bdiffs, tdiffs, align, match, misMatch, gapOpen, gapExtend, insert, deltaq, maxee, kmerSize, format, trimOverlap, createOligosGroup, createFileGroup, group, screenSequences, maxHomoP, maxLength, maxAmbig);
            data.push_back(dataBundle);

            workerThreads.push_back(new std::thread(driverContigs, dataBundle));
        }

        OutputWriter* threadMisMatchWriter = new OutputWriter(synchronizedMisMatchFile);
        OutputWriter* threadFastaTrimWriter = new OutputWriter(synchronizedOutputFastaTrimFile);
        OutputWriter* threadFastaScrapWriter = new OutputWriter(synchronizedOutputFastaScrapFile);
        OutputWriter* threadQTrimWriter = nullptr;
        OutputWriter* threadQScrapWriter = nullptr;
        if (makeQualFile) {
            threadQTrimWriter = new OutputWriter(synchronizedOutputQTrimFile);
            threadQScrapWriter = new OutputWriter(synchronizedOutputQScrapFile);
        }
        contigsData* dataBundle = new contigsData(threadFastaTrimWriter, threadFastaScrapWriter, threadQTrimWriter, threadQScrapWriter, threadMisMatchWriter, fileInputs, qualOrIndexFiles, lines[0], lines[1], qLines[0], qLines[1]);
        dataBundle->setVariables(gz, delim, nameType, offByOneTrimLength, pairedBarcodes, pairedPrimers, rpairedBarcodes, rpairedPrimers, revpairedBarcodes, revpairedPrimers, primerNames, barcodeNames, reorient, pdiffs, bdiffs, tdiffs, align, match, misMatch, gapOpen, gapExtend, insert, deltaq, maxee, kmerSize, format, trimOverlap, createOligosGroup, createFileGroup, group, screenSequences, maxHomoP, maxLength, maxAmbig);

        driverContigs(dataBundle);

        long long num = dataBundle->count;
        badNames.insert(dataBundle->badNames.begin(), dataBundle->badNames.end());
        groupMap.insert(dataBundle->groupMap.begin(), dataBundle->groupMap.end());
        for (map<string, int>::iterator it = dataBundle->groupCounts.begin(); it != dataBundle->groupCounts.end(); it++) {
            map<string, int>::iterator itMine = groupCounts.find(it->first);
            if (itMine != groupCounts.end()) { itMine->second += it->second; }
            else { groupCounts[it->first] = it->second; }
        }

        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            num += data[i]->count;

            delete data[i]->trimFileName;
            delete data[i]->scrapFileName;
            delete data[i]->misMatchesFile;
            if (makeQualFile) {
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

        delete threadFastaTrimWriter;
        delete threadFastaScrapWriter;
        delete threadMisMatchWriter;
        if (makeQualFile) {
            delete threadQTrimWriter;
            delete threadQScrapWriter;
        }
        delete dataBundle;

        return num;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "createProcesses");
        exit(1);
    }
}
//**********************************************************************************************************************
//process one file at a time, only get here with gz=true
void driverContigsGroups(groupContigsData* gparams) {
    try {
        gparams->count = 0;
        gparams->bundle->delim = '@';

        for (int l = gparams->start; l < gparams->end; l++) {
            int startTime = time(nullptr);

            if (gparams->bundle->m->getControl_pressed()) { break; }

            gparams->bundle->m->mothurOut("\n>>>>>\tProcessing file pair " + gparams->fileInputs[l][0] + " - " + gparams->fileInputs[l][1] + " (files " + toString(l+1) + " of " + toString(gparams->fileInputs.size()) + ")\t<<<<<\n");

            vector<string>  theseFileInputs;
            vector<string>  theseQIInputs;
            string ffastqfile = gparams->fileInputs[l][0]; theseFileInputs.push_back(ffastqfile);
            string rfastqfile = gparams->fileInputs[l][1]; theseFileInputs.push_back(rfastqfile);
            string findexfile = gparams->fileInputs[l][2]; theseQIInputs.push_back(findexfile); //could be blank, "NONE" or filename
            string rindexfile = gparams->fileInputs[l][3]; theseQIInputs.push_back(rindexfile); //could be blank, "NONE" or filename
            gparams->bundle->group = gparams->file2Groups[l]; //blank if no group assigned to file pair
            
            bool decompressionHelped = false;
            
            //test to make sure you can read the gz files
            bool readable = testGZReadable(theseFileInputs, theseQIInputs, decompressionHelped, gparams->bundle->format, gparams->bundle->m);
                
            if (readable) {
                if (decompressionHelped) { gparams->bundle->gz = false; }
            }else {
                gparams->bundle->m->mothurOut("[ERROR]: Unable to read compressed .gz files, please decompress and run make.contigs again. \n"); gparams->bundle->m->setControl_pressed(true); break;
            }

            //find read name type to speed read matching later
            gparams->bundle->nameType = setNameType(theseFileInputs[0], theseFileInputs[1], gparams->bundle->delim, gparams->bundle->offByOneTrimLength, gparams->bundle->gz, gparams->bundle->format);

            //fake out lines - we are just going to check for end of file. Work is divided by number of files per processor.
            vector<linePair> thisLines; vector<linePair> thisQLines;
            if (decompressionHelped) {
                //set file positions for file
                int processors = 1;
                vector<double> fastaFilePos = gparams->bundle->util.divideFile(theseFileInputs[0], processors, gparams->bundle->delim);
                thisLines.push_back(linePair(fastaFilePos[0], fastaFilePos[1])); //forward fastq
                fastaFilePos = gparams->bundle->util.divideFile(theseFileInputs[1], processors, gparams->bundle->delim);
                thisLines.push_back(linePair(fastaFilePos[0], fastaFilePos[1])); //reverse fastq
                
                if ((theseQIInputs[0] != "") && (theseQIInputs[0] != "NONE")){
                    fastaFilePos = gparams->bundle->util.divideFile(theseQIInputs[0], processors, gparams->bundle->delim);
                    thisQLines.push_back(linePair(fastaFilePos[0], fastaFilePos[1])); //forward index
                }
                if ((theseQIInputs[1] != "") && (theseQIInputs[1] != "NONE")){
                    fastaFilePos = gparams->bundle->util.divideFile(theseQIInputs[1], processors, gparams->bundle->delim);
                    thisQLines.push_back(linePair(fastaFilePos[0], fastaFilePos[1])); //forward index
                }
                if (thisQLines.size() == 0) { thisQLines = thisLines; }
            }
            else {
                thisLines.push_back(linePair(0, 1000)); thisLines.push_back(linePair(0, 1000)); //fasta[0], fasta[1] - forward and reverse
                thisQLines.push_back(linePair(0, 1000)); thisQLines.push_back(linePair(0, 1000));  //qual[0], qual[1] - forward and reverse
            }
            gparams->bundle->m->mothurOut("Making contigs...\n");

            contigsData* dataBundle = new contigsData(gparams->bundle->trimFileName, gparams->bundle->scrapFileName, gparams->bundle->trimQFileName, gparams->bundle->scrapQFileName, gparams->bundle->misMatchesFile, theseFileInputs, theseQIInputs, thisLines[0], thisLines[1], thisQLines[0], thisQLines[1]);
            dataBundle->copyVariables(gparams->bundle);
            driverContigs(dataBundle);

            if (decompressionHelped) {
                gparams->bundle->util.mothurRemove(theseFileInputs[0]); gparams->bundle->util.mothurRemove(theseFileInputs[1]);
                if (theseQIInputs[0] != "NONE") { gparams->bundle->util.mothurRemove(theseQIInputs[0]); }
                if (theseQIInputs[1] != "NONE") { gparams->bundle->util.mothurRemove(theseQIInputs[1]); }
            }
            
            gparams->count += dataBundle->count;
            gparams->badNames.insert(dataBundle->badNames.begin(), dataBundle->badNames.end());
            gparams->bundle->groupMap.insert(dataBundle->groupMap.begin(), dataBundle->groupMap.end());
            for (map<string, int>::iterator it = dataBundle->groupCounts.begin(); it != dataBundle->groupCounts.end(); it++) {
                map<string, int>::iterator itMine = gparams->bundle->groupCounts.find(it->first);
                if (itMine != gparams->bundle->groupCounts.end()) { itMine->second += it->second; }
                else { gparams->bundle->groupCounts[it->first] = it->second; }
            }
            gparams->bundle->m->mothurOut("Done.\n\nIt took " + toString(time(nullptr) - startTime) + " secs to assemble " + toString(dataBundle->count) + " reads.\n\n");
            delete dataBundle;
        }
    }
    catch(exception& e) {
        gparams->bundle->m->errorOut(e, "MakeContigsCommand", "driverContigsGroups");
        exit(1);
    }
}
//**********************************************************************************************************************
//only getting here is gz=true
unsigned long long MakeContigsCommand::createProcessesGroups(vector< vector<string> > fileInputs, string compositeFastaFile, string compositeScrapFastaFile, string compositeQualFile, string compositeScrapQualFile, string compositeMisMatchFile, vector<string>& file2Groups) {
    try {
        map<int, oligosPair> pairedPrimers, rpairedPrimers, revpairedPrimers, pairedBarcodes, rpairedBarcodes, revpairedBarcodes;
        vector<string> barcodeNames, primerNames;

        if(oligosfile != "")  {   createOligosGroup = getOligos(pairedPrimers, rpairedPrimers, revpairedPrimers, pairedBarcodes, rpairedBarcodes, revpairedBarcodes, barcodeNames, primerNames);    }
        
        //give group in file file precedence
        if (createFileGroup) {  createOligosGroup = false; }

        vector<std::thread*> workerThreads;
        vector<groupContigsData*> data;

        //divide files between processors
        vector<linePair> startEndIndexes;
        int remainingPairs = fileInputs.size();
        if (remainingPairs < processors) { processors = remainingPairs; }
        int startIndex = 0;
        for (int remainingProcessors = processors; remainingProcessors > 0; remainingProcessors--) {
            int numPairs = remainingPairs; //case for last processor
            if (remainingProcessors != 1) { numPairs = ceil(remainingPairs / remainingProcessors); }
            startEndIndexes.push_back(linePair(startIndex, (startIndex+numPairs))); //startIndex, endIndex
            startIndex = startIndex + numPairs;
            remainingPairs = remainingPairs - numPairs;
        }

        auto synchronizedOutputFastaTrimFile = std::make_shared<SynchronizedOutputFile>(compositeFastaFile);
        auto synchronizedOutputFastaScrapFile = std::make_shared<SynchronizedOutputFile>(compositeScrapFastaFile);
        auto synchronizedOutputQTrimFile = std::make_shared<SynchronizedOutputFile>(compositeQualFile);
        auto synchronizedOutputQScrapFile = std::make_shared<SynchronizedOutputFile>(compositeScrapQualFile);
        auto synchronizedMisMatchFile = std::make_shared<SynchronizedOutputFile>(compositeMisMatchFile);

        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            OutputWriter* threadFastaTrimWriter = new OutputWriter(synchronizedOutputFastaTrimFile);
            OutputWriter* threadFastaScrapWriter = new OutputWriter(synchronizedOutputFastaScrapFile);
            OutputWriter* threadMismatchWriter = new OutputWriter(synchronizedMisMatchFile);
            OutputWriter* threadQTrimWriter = nullptr;
            OutputWriter* threadQScrapWriter = nullptr;
            if (makeQualFile) {
                threadQTrimWriter = new OutputWriter(synchronizedOutputQTrimFile);
                threadQScrapWriter = new OutputWriter(synchronizedOutputQScrapFile);
            }

            contigsData* dataBundle = new contigsData(threadFastaTrimWriter, threadFastaScrapWriter, threadQTrimWriter, threadQScrapWriter, threadMismatchWriter);
            dataBundle->setVariables(gz, delim, nameType, offByOneTrimLength, pairedBarcodes, pairedPrimers, rpairedBarcodes, rpairedPrimers, revpairedBarcodes, revpairedPrimers, primerNames, barcodeNames, reorient, pdiffs, bdiffs, tdiffs, align, match, misMatch, gapOpen, gapExtend, insert, deltaq, maxee, kmerSize, format, trimOverlap, createOligosGroup, createFileGroup, "", screenSequences, maxHomoP, maxLength, maxAmbig);
            groupContigsData* groupDataBundle = new groupContigsData(fileInputs, startEndIndexes[i+1].start, startEndIndexes[i+1].end, dataBundle, file2Groups);
            data.push_back(groupDataBundle);

            workerThreads.push_back(new std::thread(driverContigsGroups, groupDataBundle));
        }

        OutputWriter* threadMisMatchWriter = new OutputWriter(synchronizedMisMatchFile);
        OutputWriter* threadFastaTrimWriter = new OutputWriter(synchronizedOutputFastaTrimFile);
        OutputWriter* threadFastaScrapWriter = new OutputWriter(synchronizedOutputFastaScrapFile);
        OutputWriter* threadQTrimWriter = nullptr;
        OutputWriter* threadQScrapWriter = nullptr;
        if (makeQualFile) {
            threadQTrimWriter = new OutputWriter(synchronizedOutputQTrimFile);
            threadQScrapWriter = new OutputWriter(synchronizedOutputQScrapFile);
        }
        contigsData* dataBundle = new contigsData(threadFastaTrimWriter, threadFastaScrapWriter, threadQTrimWriter, threadQScrapWriter, threadMisMatchWriter);
        dataBundle->setVariables(gz, delim, nameType, offByOneTrimLength, pairedBarcodes, pairedPrimers, rpairedBarcodes, rpairedPrimers, revpairedBarcodes, revpairedPrimers, primerNames, barcodeNames, reorient, pdiffs, bdiffs, tdiffs, align, match, misMatch, gapOpen, gapExtend, insert, deltaq, maxee, kmerSize, format, trimOverlap, createOligosGroup, createFileGroup, "", screenSequences, maxHomoP, maxLength, maxAmbig);
        groupContigsData* groupDataBundle = new groupContigsData(fileInputs, startEndIndexes[0].start, startEndIndexes[0].end, dataBundle, file2Groups);
        driverContigsGroups(groupDataBundle);

        delete threadFastaTrimWriter;
        delete threadFastaScrapWriter;
        delete threadMisMatchWriter;
        if (makeQualFile) {
            delete threadQTrimWriter;
            delete threadQScrapWriter;
        }
        long long num = groupDataBundle->count;
        badNames.insert(dataBundle->badNames.begin(), dataBundle->badNames.end());
        groupMap.insert(groupDataBundle->bundle->groupMap.begin(), groupDataBundle->bundle->groupMap.end());
        for (map<string, int>::iterator it = dataBundle->groupCounts.begin(); it != dataBundle->groupCounts.end(); it++) {
            map<string, int>::iterator itMine = groupCounts.find(it->first);
            if (itMine != groupCounts.end()) { itMine->second += it->second; }
            else { groupCounts[it->first] = it->second; }
        }
        delete groupDataBundle;

        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            num += data[i]->count;

            delete data[i]->bundle->trimFileName;
            delete data[i]->bundle->scrapFileName;
            delete data[i]->bundle->misMatchesFile;
            if (makeQualFile) {
                delete data[i]->bundle->trimQFileName;
                delete data[i]->bundle->scrapQFileName;
            }
            badNames.insert(data[i]->badNames.begin(), data[i]->badNames.end());
            groupMap.insert(data[i]->bundle->groupMap.begin(), data[i]->bundle->groupMap.end());
            //merge counts
            for (map<string, int>::iterator it = data[i]->bundle->groupCounts.begin(); it != data[i]->bundle->groupCounts.end(); it++) {
                map<string, int>::iterator itMine = groupCounts.find(it->first);
                if (itMine != groupCounts.end()) { itMine->second += it->second; }
                else { groupCounts[it->first] = it->second; }
            }

            delete data[i];
            delete workerThreads[i];
        }

        return num;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeContigsCommand", "createProcessesGroups");
        exit(1);
    }
}
/**************************************************************************************************/

int MakeContigsCommand::setLines(vector<string> fasta, vector<string> qual, vector<linePair>& lines, vector<linePair>& qLines, char delim) {
    try {
        lines.clear();
        qLines.clear();
        vector<double> fastaFilePos;
        vector<double> qfileFilePos;
        vector<double> temp;

        nameType = setNameType(fasta[0], fasta[1], delim, offByOneTrimLength, gz, format);

#if defined NON_WINDOWS
        //set file positions for fasta file
        fastaFilePos = util.divideFile(fasta[0], processors, delim);

        //get name of first sequence in each chunk
        map<string, int> firstSeqNames; map<string, int> trimmedNames;
        for (int i = 0; i < (fastaFilePos.size()-1); i++) {
            ifstream in; util.openInputFile(fasta[0], in);
            in.seekg(fastaFilePos[i]);
            string name = "";

            if (delim == '>') {
                Sequence temp(in);
                name = temp.getName();
            }else {
                string line = util.getline(in); gobble(in);
                vector<string> pieces = util.splitWhiteSpace(line);
                name = pieces[0];
                name = name.substr(1);
                util.checkName(name);
            }
            fixName(name, nameType, offByOneTrimLength);
            firstSeqNames[name] = i;
            in.close();
        }

        map<string, int> copy;
        if (qual.size() != 0) { copy = firstSeqNames; }

        //look for match in reverse file
        ifstream in2;
        util.openInputFile(fasta[1], in2);

        string input;
        while(!in2.eof()){
            input = util.getline(in2); gobble(in2);

            if (input.length() != 0) {
                if(input[0] == delim){ //this is a name line
                    vector<string> pieces = util.splitWhiteSpace(input);
                    string name = pieces[0];
                    name = name.substr(1);
                    util.checkName(name);
                    fixName(name, nameType, offByOneTrimLength);

                    map<string, int>::iterator it = firstSeqNames.find(name);

                    if (it != firstSeqNames.end())  { //this is the start of a new chunk
                        double pos = in2.tellg();
                        qfileFilePos.push_back(pos - input.length() - 1);
                        firstSeqNames.erase(it);
                    }
                }
            }

            if ((firstSeqNames.size() == 0)) { break; }
        }
        in2.close();

        //get last file position of reverse fasta[1]
        FILE * pFile;
        double size;

        //get num bytes in file
        fasta[1] = util.getFullPathName(fasta[1]);
        pFile = fopen (fasta[1].c_str(),"rb");
        if (pFile==nullptr) perror ("Error opening file");
        else{
            fseek (pFile, 0, SEEK_END);
            size=ftell (pFile);
            fclose (pFile);
        }
        qfileFilePos.push_back(size);


        if ((firstSeqNames.size() != 0)){
            for (map<string, int>::iterator it = firstSeqNames.begin(); it != firstSeqNames.end(); it++) {
                if (delim == '>') {
                    m->mothurOut(it->first + " is in your forward fasta file and not in your reverse file, please remove it using the remove.seqs command before proceeding.\n"); 
                }else {
                    m->mothurOut(it->first + " is in your forward fastq file and not in your reverse file, please remove it using the remove.seqs command before proceeding.\n"); 
                }
            }
            m->setControl_pressed(true);
            return processors;
        }

        //fill lines with paired forward and reverse fasta lines
        for (int i = 0; i < (fastaFilePos.size()-1); i++) {
            if (m->getDebug()) { m->mothurOut("[DEBUG]: forward " + toString(i) +'\t' + toString(fastaFilePos[i]) + '\t' + toString(fastaFilePos[i+1]) + '\n'); }
            lines.push_back(linePair(fastaFilePos[i], fastaFilePos[(i+1)]));
            if (m->getDebug()) { m->mothurOut("[DEBUG]: reverse " + toString(i) +'\t' + toString(qfileFilePos[i]) + '\t' + toString(qfileFilePos[i+1]) + '\n'); }
            lines.push_back(linePair(qfileFilePos[i], qfileFilePos[(i+1)]));
        }

        qfileFilePos.clear();


        if (qual.size() != 0) {
            firstSeqNames = copy;

            if (qual[0] != "NONE") {
                //seach for filePos of each first name in the qfile and save in qfileFilePos
                ifstream inQual;
                util.openInputFile(qual[0], inQual);

                string input;
                while(!inQual.eof()){
                    input = util.getline(inQual); gobble(inQual);

                    if (input.length() != 0) {
                        if(input[0] == delim){ //this is a sequence name line
                            vector<string> pieces = util.splitWhiteSpace(input);
                            string name = pieces[0];
                            name = name.substr(1);
                            util.checkName(name);
                            fixName(name, nameType, offByOneTrimLength);

                            map<string, int>::iterator it = firstSeqNames.find(name);

                            if(it != firstSeqNames.end()) { //this is the start of a new chunk
                                double pos = inQual.tellg();
                                qfileFilePos.push_back(pos - input.length() - 1);
                                firstSeqNames.erase(it);
                            }                        }
                    }

                    if ((firstSeqNames.size() == 0)) { break; }
                }
                inQual.close();

                //get last file position of reverse qual[0]
                FILE * pFile;
                double size;

                //get num bytes in file
                qual[0] = util.getFullPathName(qual[0]);
                pFile = fopen (qual[0].c_str(),"rb");
                if (pFile==nullptr) perror ("Error opening file");
                else{
                    fseek (pFile, 0, SEEK_END);
                    size=ftell (pFile);
                    fclose (pFile);
                }
                qfileFilePos.push_back(size);


                if ((firstSeqNames.size() != 0)){
                    for (map<string, int>::iterator it = firstSeqNames.begin(); it != firstSeqNames.end(); it++) {
                        if (delim == '>') {
                            m->mothurOut(it->first + " is in your forward fasta file and reverse fasta file, but not your forward qfile, please remove it using the remove.seqs command before proceeding.\n"); 
                        }else {
                            m->mothurOut(it->first + " is in your forward fastq file and reverse fastq file, but not your forward index, please remove it using the remove.seqs command before proceeding.\n"); 
                        }
                    }
                    m->setControl_pressed(true);
                    return processors;
                }
            }
            firstSeqNames = copy;

            if (qual[1] != "NONE") {
                ifstream inQual2;
                util.openInputFile(qual[1], inQual2);

                while(!inQual2.eof()){
                    input = util.getline(inQual2); gobble(inQual2);

                    if (input.length() != 0) {
                        if(input[0] == delim){ //this is a sequence name line
                            vector<string> pieces = util.splitWhiteSpace(input);
                            string name = pieces[0];
                            name = name.substr(1);

                            util.checkName(name);
                            fixName(name, nameType, offByOneTrimLength);

                            map<string, int>::iterator it = firstSeqNames.find(name);

                            if(it != firstSeqNames.end()) { //this is the start of a new chunk
                                double pos = inQual2.tellg();
                                temp.push_back(pos - input.length() - 1);
                                firstSeqNames.erase(it);
                            }
                        }
                    }

                    if ((firstSeqNames.size() == 0)) { break; }
                }
                inQual2.close();

                //get last file position of reverse qual[1]
                FILE * pFile2;

                //get num bytes in file
                qual[1] = util.getFullPathName(qual[1]);
                pFile2 = fopen (qual[1].c_str(),"rb");
                if (pFile2==nullptr) perror ("Error opening file");
                else{
                    fseek (pFile2, 0, SEEK_END);
                    size=ftell (pFile2);
                    fclose (pFile2);
                }
                temp.push_back(size);


                if ((firstSeqNames.size() != 0)){
                    for (map<string, int>::iterator it = firstSeqNames.begin(); it != firstSeqNames.end(); it++) {
                        if (delim == '>') {
                            m->mothurOut(it->first + " is in your forward fasta file, reverse fasta file, and forward qfile but not your reverse qfile, please remove it using the remove.seqs command before proceeding.\n"); 
                        }else {
                            if (qual[0] != "NONE") {
                                m->mothurOut(it->first + " is in your forward fastq file, reverse fastq file, and forward index but not your reverse index, please remove it using the remove.seqs command before proceeding.\n"); 
                            }else {
                                m->mothurOut(it->first + " is in your forward fastq file, reverse fastq file, but not your reverse index, please remove it using the remove.seqs command before proceeding.\n"); 
                            }
                        }
                    }
                    m->setControl_pressed(true);
                    return processors;
                }
            }

            if (qual[0] == "NONE") { qfileFilePos = temp; } //fill with duds, if both were NONE then qual.size() == 0
            if (qual[1] == "NONE") { temp = qfileFilePos; } //fill with duds, if both were NONE then qual.size() == 0


            //fill lines with paired forward and reverse fasta lines
            for (int i = 0; i < (fastaFilePos.size()-1); i++) {
                if (m->getDebug()) { m->mothurOut("[DEBUG]: forward " + toString(i) +'\t' + toString(qfileFilePos[i]) + '\t' + toString(qfileFilePos[i+1]) + '\n'); }
                qLines.push_back(linePair(qfileFilePos[i], qfileFilePos[(i+1)]));
                if (m->getDebug()) { m->mothurOut("[DEBUG]: reverse " + toString(i) +'\t' + toString(temp[i]) + '\t' + toString(temp[i+1]) + '\n'); }
                qLines.push_back(linePair(temp[i], temp[(i+1)]));
            }
        }else {  qLines = lines;	} //files with duds


        return processors;

#else

        
            long long numFastaSeqs = 0;
            fastaFilePos = util.setFilePosFasta(fasta[0], numFastaSeqs, delim); //forward
            if (numFastaSeqs < processors) { processors = numFastaSeqs; }

            long long numRFastaSeqs = 0;
            qfileFilePos = util.setFilePosFasta(fasta[1], numRFastaSeqs, delim); //reverse

            if (numFastaSeqs != numRFastaSeqs) {
                if (delim == '>') {
                    m->mothurOut("[ERROR]: You have " + toString(numFastaSeqs) + " sequences in your forward fasta file, but " + toString(numRFastaSeqs) + " sequences in your reverse fasta file. Please use the list.seqs and get.seqs commands to make the files match before proceeding. list.seqs(fasta=yourForward.fasta-yourReverse.fasta);get.seqs(fasta=yourForward.fasta, accnos=current);get.seqs(fasta=yourReverse.fasta, accnos=current);\n");  m->setControl_pressed(true); return processors;
                }else {
                    m->mothurOut("[ERROR]: You have " + toString(numFastaSeqs) + " sequences in your forward fastq file, but " + toString(numRFastaSeqs) + " sequences in your reverse fastq file. Please use the list.seqs and get.seqs commands to make the files match before proceeding. list.seqs(fastq=yourForward.fastq-yourReverse.fastq);get.seqs(fastq=yourForward.fastq, accnos=current);get.seqs(fastq=yourReverse.fastq, accnos=current);\n"); m->setControl_pressed(true); return processors;
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

                if (qual[0] != "NONE") {  fastaFilePos = util.setFilePosFasta(qual[0], numFQualSeqs, delim);  } //forward index or qual file
                if (qual[1] != "NONE") {  qfileFilePos = util.setFilePosFasta(qual[1], numRQualSeqs, delim);  }//reverse index or qual file

                if (qual[0] == "NONE") { fastaFilePos = qfileFilePos; numFQualSeqs = numRQualSeqs; } //fill with duds, if both were NONE then qual.size() == 0
                if (qual[1] == "NONE") { qfileFilePos = fastaFilePos; numRQualSeqs = numFQualSeqs; } //fill with duds, if both were NONE then qual.size() == 0

                if ((numFQualSeqs != numRQualSeqs) || (numFQualSeqs != numFastaSeqs)){
                    if (delim == '>') {
                        m->mothurOut("[ERROR]: You have " + toString(numFastaSeqs) + " sequences in your forward fasta file, " + toString(numRFastaSeqs) + " sequences in your reverse fasta file, " + toString(numFQualSeqs) + " sequences in your forward qual file, " + toString(numRQualSeqs) + " sequences in your reverse qual file. Please use the list.seqs and get.seqs commands to make the files match before proceeding.\n");  m->setControl_pressed(true); return processors;
                    }else {
                        if (qual[0] != "NONE") {
                            m->mothurOut("[ERROR]: You have " + toString(numFastaSeqs) + " sequences in your forward fastq file, " + toString(numRFastaSeqs) + " sequences in your reverse fastq file and " + toString(numRQualSeqs) + " sequences in your reverse index file. Please use the list.seqs and get.seqs commands to make the files match before proceeding.\n");  m->setControl_pressed(true); return processors;
                        }else if (qual[1] != "NONE") {
                            m->mothurOut("[ERROR]: You have " + toString(numFastaSeqs) + " sequences in your forward fastq file, " + toString(numRFastaSeqs) + " sequences in your reverse fastq file and " + toString(numFQualSeqs) + " sequences in your forward index file. Please use the list.seqs and get.seqs commands to make the files match before proceeding.\n");  m->setControl_pressed(true); return processors;
                        }else {
                            m->mothurOut("[ERROR]: You have " + toString(numFastaSeqs) + " sequences in your forward fastq file, " + toString(numRFastaSeqs) + " sequences in your reverse fastq file, " + toString(numFQualSeqs) + " sequences in your forward index file, " + toString(numRQualSeqs) + " sequences in your reverse index file. Please use the list.seqs and get.seqs commands to make the files match before proceeding.\n");  m->setControl_pressed(true); return processors;
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
    
        if(qual.size() == 0)	{	qLines = lines;	} //files with duds
        
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
vector< vector<string> > MakeContigsCommand::readFileNames(string filename, vector<string>& file2Group){
	try {
        
        FileFile dataFile(filename, "contigs");
        vector< vector<string> > files = dataFile.getFiles();
        gz = dataFile.isGZ();
        file2Group = dataFile.getGroupNames();
        createFileGroup = dataFile.isColumnWithGroupNames();
        if (dataFile.containsIndexFiles() && (oligosfile == "")) { m->mothurOut("[ERROR]: You need to provide an oligos file if you are going to use an index file.\n"); m->setControl_pressed(true);  }
        
        if (files.size() == 0) { m->setControl_pressed(true); }

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
bool MakeContigsCommand::getOligos(map<int, oligosPair>& pairedPrimers, map<int, oligosPair>& rpairedPrimers, map<int, oligosPair>& revpairedPrimers, map<int, oligosPair>& pairedBarcodes, map<int, oligosPair>& rpairedBarcodes, map<int, oligosPair>& revpairedBarcodes, vector<string>& barcodeNames, vector<string>& primerNames){
	try {
        if (m->getDebug()) { m->mothurOut("[DEBUG]: oligosfile = " + oligosfile + "\n"); }

        bool allBlank = false;
        Oligos oligos; oligos.read(oligosfile, false);

        if (m->getControl_pressed()) { return false; } //error in reading oligos

        if (oligos.hasPairedBarcodes() || oligos.hasPairedPrimers()) {
            pairedPrimers = oligos.getPairedPrimers();
            rpairedPrimers = oligos.getReorientedPairedPrimers();
            revpairedPrimers = oligos.getReversedPairedPrimers();
            primerNames = oligos.getPrimerNames();
            pairedBarcodes = oligos.getPairedBarcodes();
            rpairedBarcodes = oligos.getReorientedPairedBarcodes();
            revpairedBarcodes = oligos.getReversedPairedBarcodes();
            barcodeNames = oligos.getBarcodeNames();
            
            if (m->getDebug()) {
                map<int, oligosPair>::iterator it;
                m->mothurOut("\n[DEBUG]: paired primers - \n");
                for (it = pairedPrimers.begin(); it != pairedPrimers.end(); it++) {
                    m->mothurOut("[DEBUG]: " + primerNames[it->first] + "\t" + it->second.forward + "\t" + it->second.reverse + "\n");
                }
                m->mothurOut("\n[DEBUG]: paired reoriented primers - \n");
                for (it = rpairedPrimers.begin(); it != rpairedPrimers.end(); it++) {
                    m->mothurOut("[DEBUG]: " + primerNames[it->first] + "\t" + it->second.forward + "\t" + it->second.reverse + "\n");
                }
                m->mothurOut("\n[DEBUG]: paired reversed primers - \n");
                for (it = revpairedPrimers.begin(); it != revpairedPrimers.end(); it++) {
                    m->mothurOut("[DEBUG]: " + primerNames[it->first] + "\t" + it->second.forward + "\t" + it->second.reverse + "\n");
                }
                
                m->mothurOut("\n[DEBUG]: paired barcodes - \n");
                for (it = pairedBarcodes.begin(); it != pairedBarcodes.end(); it++) {
                    m->mothurOut("[DEBUG]: " + barcodeNames[it->first] + "\t" + it->second.forward + "\t" + it->second.reverse + "\n");
                }
                m->mothurOut("\n[DEBUG]: paired reoriented barcodes - \n");
                for (it = rpairedBarcodes.begin(); it != rpairedBarcodes.end(); it++) {
                    m->mothurOut("[DEBUG]: " + barcodeNames[it->first] + "\t" + it->second.forward + "\t" + it->second.reverse + "\n");
                }
                m->mothurOut("\n[DEBUG]: paired reversed barcodes - \n");
                for (it = revpairedBarcodes.begin(); it != revpairedBarcodes.end(); it++) {
                    m->mothurOut("[DEBUG]: " + barcodeNames[it->first] + "\t" + it->second.forward + "\t" + it->second.reverse + "\n");
                }
            }
        }else {
            m->mothurOut("[ERROR]: make.contigs requires paired barcodes and primers. You can set one end to NONE if you are using an index file.\n"); m->setControl_pressed(true);
        }

        if (m->getControl_pressed()) { return false; }

        int numLinkers = oligos.getLinkers().size();
        int numSpacers = oligos.getSpacers().size();
        if (numLinkers != 0) { m->mothurOut("[WARNING]: make.contigs is not setup to remove linkers, ignoring.\n"); }
        if (numSpacers != 0) { m->mothurOut("[WARNING]: make.contigs is not setup to remove spacers, ignoring.\n"); }

        vector<string> groupNames = oligos.getGroupNames();
        if (groupNames.size() == 0) { allFiles = false; allBlank = true;  }

        if (allBlank) {
            m->mothurOut("[WARNING]: your oligos file does not contain any group names.  mothur will not create a groupfile.\n"); 
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
void MakeContigsCommand::debugFunction() {
    try{
//allows you to run the oligos and index file independantly to check for barcode issues. make.contigs(findex=yourIndexFile, bdiffs=1, oligos=yourOligosFile, checkorient=t). just used for user support
    
    map<int, oligosPair> pairedPrimers, pairedBarcodes, reorientedPairedBarcodes, reorientedPairedPrimers, reversedPairedBarcodes, reverseedPairedPrimers;
    vector<string> barcodeNames, primerNames;
    
    if(oligosfile != "")                        {       createOligosGroup = getOligos(pairedPrimers, reorientedPairedPrimers, reverseedPairedPrimers, pairedBarcodes, reorientedPairedBarcodes, reversedPairedBarcodes,  barcodeNames, primerNames);    }
    
    int numPrimers = pairedPrimers.size();
    TrimOligos trimOligos(pdiffs, bdiffs, 0, 0, pairedPrimers, pairedBarcodes, true);
    int numBarcodes = pairedBarcodes.size();
    TrimOligos* rtrimOligos = nullptr;
    if (reorient) {  rtrimOligos = new TrimOligos(pdiffs, bdiffs, 0, 0, reorientedPairedPrimers, reorientedPairedBarcodes, true); numBarcodes = reorientedPairedBarcodes.size();   numPrimers = reorientedPairedPrimers.size();  }
    
    ifstream in; util.openInputFile(findexfile, in);
    
    while (!in.eof()) {
        if (m->getControl_pressed()) { break; }
        
        bool ignore = false;
        FastqRead index(in, ignore, format); gobble(in);
        
        int success = 1;
        string trashCode = "";
        string commentString = "";
        int currentSeqsDiffs = 0;
        int barcodeIndex = 0;
        
        Sequence fSeq, rSeq;
        QualityScores* fQual = nullptr; QualityScores* rQual = nullptr;
        QualityScores* savedFQual = nullptr; QualityScores* savedRQual = nullptr;
        Sequence findexBarcode("findex", index.getSeq());  Sequence rindexBarcode("rindex", "NONE");
        Sequence savedFindex("findex", index.getSeq());  Sequence savedRIndex("rindex", "NONE");
        
        if(numBarcodes != 0){
            vector<int> results;
            
            results = trimOligos.stripBarcode(findexBarcode, rindexBarcode, *fQual, *rQual, barcodeIndex);
            
            success = results[0] + results[2];
            commentString += "fbdiffs=" + toString(results[0]) + "(" + trimOligos.getCodeValue(results[1], bdiffs) + "), rbdiffs=" + toString(results[2]) + "(" + trimOligos.getCodeValue(results[3], bdiffs) + ") ";
            if(success > bdiffs)        {    trashCode += 'b';    }
            else{ currentSeqsDiffs += success;  }
        }
        
        if (reorient && (trashCode != "")) { //if you failed and want to check the reverse
            int thisSuccess = 0;
            string thisTrashCode = "";
            string thiscommentString = "";
            int thisCurrentSeqsDiffs = 0;
            
            int thisBarcodeIndex = 0;
            
            if(numBarcodes != 0){
                vector<int> results;
                
                results = rtrimOligos->stripBarcode(savedFindex, savedRIndex, *savedFQual, *savedRQual, thisBarcodeIndex);
                
                thisSuccess = results[0] + results[2];
                thiscommentString += "fbdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], bdiffs) + "), rbdiffs=" + toString(results[2]) + "(" + rtrimOligos->getCodeValue(results[3], bdiffs) + ") ";
                if(thisSuccess > bdiffs)        {    thisTrashCode += 'b';    }
                else{ thisCurrentSeqsDiffs += thisSuccess;  }
            }
            
            if (thisTrashCode == "") {
                trashCode = thisTrashCode;
                success = thisSuccess;
                currentSeqsDiffs = thisCurrentSeqsDiffs;
                commentString = thiscommentString;
                barcodeIndex = thisBarcodeIndex;
            }else { trashCode += "(" + thisTrashCode + ")";  }
        }
        
        
        if (trashCode == "") {
            string thisGroup = "";
            if(numBarcodes != 0){ thisGroup = barcodeNames[barcodeIndex]; }
            
            
            int pos = thisGroup.find("ignore");
            if (pos == string::npos) {
                if (thisGroup != "") {
                    groupMap[index.getName()] = thisGroup;
                    
                    map<string, int>::iterator it = groupCounts.find(thisGroup);
                    if (it == groupCounts.end()) {    groupCounts[thisGroup] = 1; }
                    else { groupCounts[it->first] ++; }
                }
            }
        }
        }
        in.close();
        
        int total = 0;
        if (groupCounts.size() != 0) {  m->mothurOut("\nGroup count: \n");  }
        for (map<string, int>::iterator it = groupCounts.begin(); it != groupCounts.end(); it++) { total += it->second; m->mothurOut(it->first + "\t" + toString(it->second) + "\n"); }
        if (total != 0) { m->mothurOut("\nTotal of all groups is " + toString(total) + "\n"); }
        
        exit(1);
        
}
catch(exception& e) {
    m->errorOut(e, "MakeContigsCommand", "debugFunction");
    exit(1);
}
}
//**********************************************************************************************************************
