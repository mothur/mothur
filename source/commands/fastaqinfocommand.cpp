/*
 *  parsefastaqcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 9/30/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "fastaqinfocommand.h"
#include "sequence.hpp"
#include "counttable.h"

//**********************************************************************************************************************
vector<string> ParseFastaQCommand::setParameters(){	
	try {
        CommandParameter pfile("file", "InputTypes", "", "", "fastqFile", "fastqFile", "none","",false,false,true); parameters.push_back(pfile);
		CommandParameter pfastq("fastq", "InputTypes", "", "", "fastqFile", "fastqFile", "none","",false,false,true); parameters.push_back(pfastq);
        CommandParameter poligos("oligos", "InputTypes", "", "", "oligosGroup", "none", "none","",false,false); parameters.push_back(poligos);
        CommandParameter pgroup("group", "InputTypes", "", "", "oligosGroup", "none", "none","",false,false); parameters.push_back(pgroup);
        CommandParameter preorient("checkorient", "Boolean", "", "F", "", "", "","",false,false,true); parameters.push_back(preorient);
        CommandParameter ppdiffs("pdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(ppdiffs);
		CommandParameter pbdiffs("bdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pbdiffs);
        CommandParameter pldiffs("ldiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pldiffs);
		CommandParameter psdiffs("sdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(psdiffs);
        CommandParameter ptdiffs("tdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(ptdiffs);
		CommandParameter pfasta("fasta", "Boolean", "", "T", "", "", "","fasta",false,false); parameters.push_back(pfasta);
		CommandParameter pqual("qfile", "Boolean", "", "T", "", "", "","qfile",false,false); parameters.push_back(pqual);
        CommandParameter ppacbio("pacbio", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(ppacbio);
 		CommandParameter pformat("format", "Multiple", "sanger-illumina-solexa-illumina1.8+", "illumina1.8+", "", "", "","",false,false,true); parameters.push_back(pformat);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false; fileOption = 0; createFileGroup = false; hasIndex = false;
        split = 1;
        
        vector<string> tempOutNames;
        outputTypes["fasta"] = tempOutNames;
        outputTypes["qfile"] = tempOutNames;
        outputTypes["fastq"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ParseFastaQCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The fastq.info command reads a fastq file and creates a fasta and quality file or can be used to parse fastq files by sample.\n";
		helpString += "The fastq.info command parameters are file, fastq, fasta, qfile, oligos, group and format; file or fastq is required.\n";
        helpString += "The fastq.info command should be in the following format: fastq.info(fastaq=yourFastaQFile).\n";
        helpString += "The oligos parameter allows you to provide an oligos file to split your fastq file into separate fastq files by barcode and primers. \n";
        helpString += "The group parameter allows you to provide a group file to split your fastq file into separate fastq files by group. \n";
        helpString += "The tdiffs parameter is used to specify the total number of differences allowed in the reads. The default is pdiffs + bdiffs + sdiffs + ldiffs.\n";
		helpString += "The bdiffs parameter is used to specify the number of differences allowed in the barcode. The default is 0.\n";
		helpString += "The pdiffs parameter is used to specify the number of differences allowed in the primer. The default is 0.\n";
        helpString += "The ldiffs parameter is used to specify the number of differences allowed in the linker. The default is 0.\n";
		helpString += "The sdiffs parameter is used to specify the number of differences allowed in the spacer. The default is 0.\n";
        helpString += "The checkorient parameter will check look for the reverse compliment of the barcode or primer in the sequence. If found the sequence is flipped. The default is false.\n";
		helpString += "The format parameter is used to indicate whether your sequences are sanger, solexa, illumina1.8+ or illumina, default=illumina1.8+.\n";
        helpString += "The fasta parameter allows you to indicate whether you want a fasta file generated. Default=T.\n";
        helpString += "The qfile parameter allows you to indicate whether you want a quality file generated. Default=T.\n";
        helpString += "The pacbio parameter allows you to indicate .... When set to true, quality scores of 0 will results in a corresponding base of N. Default=F.\n";
		helpString += "Example fastq.info(fastaq=test.fastaq).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ParseFastaQCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],fasta-[filename],[sample],[tag],fasta-[filename],[sample],fasta"; }
        else if (type == "count") {  pattern = "[filename],count_table"; }
        else if (type == "qfile") {  pattern = "[filename],qual-[filename],[sample],[tag],qual-[filename],[sample],qual"; }
        else if (type == "fastq") {  pattern = "[filename],[sample],fastq-[filename],[sample],[tag],fastq"; } //make.sra assumes the [filename],[sample],[tag],fastq format for the 4 column file option. If this changes, may have to modify fixMap function.
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ParseFastaQCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ParseFastaQCommand::ParseFastaQCommand(string option) : Command(){
	try {

		if(option == "help") {	help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			fastaQFile = validParameter.validFile(parameters, "fastq");
			if (fastaQFile == "not found") {	fastaQFile= "";	}
			else if (fastaQFile == "not open")	{	fastaQFile = ""; abort = true;	}
            else { inputfile = fastaQFile; }
            
            file = validParameter.validFile(parameters, "file");
			if (file == "not found") {	file = "";	}
			else if (file == "not open")	{	file = ""; abort = true;	}
            else { inputfile = file; fileOption = true; }
            
            if ((file == "") && (fastaQFile == "")) {  m->mothurOut("You must provide a file or fastq option.\n"); abort = true;  }

            oligosfile = validParameter.validFile(parameters, "oligos");
			if (oligosfile == "not found") {	oligosfile = "";	}
			else if (oligosfile == "not open")	{	oligosfile = ""; abort = true;	}
            else { current->setOligosFile(oligosfile); split = 2; }
            
            groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not found") {	groupfile = "";	}
			else if (groupfile == "not open")	{	groupfile = ""; abort = true;	}
            else { current->setGroupFile(groupfile); split = 2; }
            
            if ((groupfile != "") && (oligosfile != "")) { m->mothurOut("You must enter ONLY ONE of the following: oligos or group.\n");  abort = true;  }
			
			 
				if (outputdir == ""){    outputdir = util.hasPath(inputfile); 	}
			
			string temp;
			temp = validParameter.valid(parameters, "fasta");	if(temp == "not found"){	temp = "T";	}
			fasta = util.isTrue(temp); 

			temp = validParameter.valid(parameters, "qfile");	if(temp == "not found"){	temp = "T";	}
			qual = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "pacbio");	if(temp == "not found"){	temp = "F";	}
			pacbio = util.isTrue(temp);

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

			
            format = validParameter.valid(parameters, "format");		if (format == "not found"){	format = "illumina1.8+";	}
            
            if ((format != "sanger") && (format != "illumina") && (format != "illumina1.8+") && (format != "solexa"))  { 
				m->mothurOut(format + " is not a valid format. Your format choices are sanger, solexa, illumina1.8+ and illumina, aborting." ); m->mothurOutEndLine();
				abort=true;
			}

            if ((!fasta) && (!qual) && (file == "") && (fastaQFile == "") && (oligosfile == "")) { m->mothurOut("[ERROR]: no outputs selected. Aborting.\n");  abort=true; }
            temp = validParameter.valid(parameters, "checkorient");		if (temp == "not found") { temp = "F"; }
			reorient = util.isTrue(temp);

		}		
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "ParseFastaQCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int ParseFastaQCommand::execute(){
	try {
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
        vector< vector<string> > files;
        if (file != "") { files = readFile(); }
        
        if (m->getControl_pressed()) { return 0; }
        
        TrimOligos* trimOligos = nullptr; TrimOligos* rtrimOligos = nullptr;
        pairedOligos = false; numBarcodes = 0; numPrimers= 0; numLinkers= 0; numSpacers = 0; numRPrimers = 0;
        if (oligosfile != "")       {
            readOligos(oligosfile);
            //find group read belongs to
            if (pairedOligos)   {   trimOligos = new TrimOligos(pdiffs, bdiffs, 0, 0, oligos.getPairedPrimers(), oligos.getPairedBarcodes(), hasIndex); numBarcodes = oligos.getPairedBarcodes().size(); numPrimers = oligos.getPairedPrimers().size(); }
            else                {   trimOligos = new TrimOligos(pdiffs, bdiffs, ldiffs, sdiffs, oligos.getPrimers(), oligos.getBarcodes(), oligos.getReversePrimers(), oligos.getLinkers(), oligos.getSpacers());  numPrimers = oligos.getPrimers().size(); numBarcodes = oligos.getBarcodes().size();  }
            
            if (reorient) {
                rtrimOligos = new TrimOligos(pdiffs, bdiffs, 0, 0, oligos.getReorientedPairedPrimers(), oligos.getReorientedPairedBarcodes(), hasIndex); numBarcodes = oligos.getReorientedPairedBarcodes().size();
            }
            
        }else if (groupfile != "")   { readGroup(groupfile);     }
        
        string inputFile = "";
        if (file != "") {
            
            inputFile = file;
            
            if (m->getControl_pressed()) {
                if (groupfile != "")        { delete groupMap;      }
                if (oligosfile != "")  { delete trimOligos; if (reorient) { delete rtrimOligos; }   }
                return 0;
            }
            
            //groupfile name for pacbio with option 2
            map<string, string> variables;
            variables["[filename]"] = util.getRootName(file);
            string pacbioFastaFileName = getOutputFileName("fasta", variables);
            string pacbioQualFileName = getOutputFileName("qfile", variables);
            if ((fileOption == 2) && pacbio) {
                seqGroups.clear();
                
                if (fasta) {
                    ofstream temppbf; util.openOutputFile(pacbioFastaFileName, temppbf);
                    temppbf.close();
                    outputNames.push_back(pacbioFastaFileName); outputTypes["fasta"].push_back(pacbioFastaFileName);
                }
                
                if (qual) {
                    ofstream temppbq; util.openOutputFile(pacbioQualFileName, temppbq);
                    temppbq.close();
                    outputNames.push_back(pacbioQualFileName); outputTypes["qfile"].push_back(pacbioQualFileName);
                }
                
            } //clear old file for append
            
            for (int i = 0; i < files.size(); i++) { //process each pair
                
                if (m->getControl_pressed()) { break; }
                
                if (((fileOption == 2) || (fileOption == 4)) && !pacbio)  { //2 column and 4 column format file file
                    processFile(files[i], trimOligos, rtrimOligos);
                }else if ((fileOption == 2) && pacbio)  { //pacbio with group filename option
                    split = 1;
                    
                    if (current->getMothurCalling()) {
                        //add group names to fastq files and make copies - for sra command parse
                        ofstream temp;
                        map<string, string> variables;
                        variables["[filename]"] = util.getRootName(files[i][0]);
                        variables["[sample]"] = file2Group[i];
                        variables["[tag]"] = "";
                        string newfqFile = getOutputFileName("fastq", variables);
                        util.openOutputFile(newfqFile, temp);        temp.close();
                        util.appendFiles(files[i][0], newfqFile);
                        outputNames.push_back(newfqFile); outputTypes["fastq"].push_back(newfqFile);
                    }
                    
                    inputFile = files[i][0];
                    
                    //process each file to create fasta and qual files
                    set<string> seqNames;
                    if (fasta || qual) {  seqNames = processFile(inputFile, trimOligos, rtrimOligos); }  //split = 1, so no parsing by group will be done.
                    
                    if (seqNames.size() != 0) {
                        string pacbioGroup = file2Group[i];
                        for (set<string>::iterator it = seqNames.begin(); it != seqNames.end(); it++) {
                            seqGroups[*it] = pacbioGroup;
                        }
                        groupCounts[pacbioGroup] = seqNames.size();
                        
                        map<string, string> variables;
                        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(inputFile));
                        string fastaFile = getOutputFileName("fasta",variables);
                        string qualFile = getOutputFileName("qfile",variables);
                        if (fasta) { util.appendFiles(fastaFile, pacbioFastaFileName); }
                        if (qual)  { util.appendFiles(qualFile, pacbioQualFileName);   }
                    }
                }else if (fileOption == 3) { //3 column file option with sample names
                    if (current->getMothurCalling()) {
                        //add group names to fastq files and make copies - for sra command parse
                        ofstream temp, temp2;
                        map<string, string> variables;
                        variables["[filename]"] = util.getRootName(files[i][0]);
                        variables["[sample]"] = file2Group[i];
                        variables["[tag]"] = "forward";
                        string newffqFile = getOutputFileName("fastq", variables);
                        util.openOutputFile(newffqFile, temp);		temp.close();
                        util.appendFiles(files[i][0], newffqFile);
                        outputNames.push_back(newffqFile); outputTypes["fastq"].push_back(newffqFile);
                        
                        variables["[filename]"] = util.getRootName(files[i][1]);
                        variables["[sample]"] = file2Group[i];
                        variables["[tag]"] = "reverse";
                        string newfrqFile = getOutputFileName("fastq", variables);
                        util.openOutputFile(newfrqFile, temp2);		temp2.close();
                        util.appendFiles(files[i][1], newfrqFile);
                        outputNames.push_back(newfrqFile); outputTypes["fastq"].push_back(newfrqFile);
                    }
                    
                    //if requested, make fasta and qual
                    if (fasta || qual) {  processFile(files[i], trimOligos, rtrimOligos); }  //split = 1, so no parsing by group will be done.
                }
            }
        }else {
            inputFile = fastaQFile;
            processFile(fastaQFile, trimOligos, rtrimOligos);
            vector<string> filesFakeOut; filesFakeOut.push_back(fastaQFile);
            files.push_back(filesFakeOut);
        }
        
        if ((fileOption == 2) && pacbio) {
            map<string, string> variables;
            variables["[filename]"] = util.getRootName(file);
            string pacbioGroupFileName = getOutputFileName("count", variables);
            outputNames.push_back(pacbioGroupFileName); outputTypes["count"].push_back(pacbioGroupFileName);

            CountTable ct; ct.createTable(seqGroups);
            ct.printCompressedTable(pacbioGroupFileName);
        }
		
        if (split > 1) {
            
            string thisOutputDir = outputdir;
            if (outputdir == "") {  thisOutputDir = util.hasPath(inputFile); }
            map<string, string> vars;
            vars["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(inputFile));
            string outputGroupFileName = getOutputFileName("count",vars);

            if (seqGroups.size() != 0) { //Create count file
                outputNames.push_back(outputGroupFileName); outputTypes["count"].push_back(outputGroupFileName);
                
                CountTable ct; ct.createTable(seqGroups);
                ct.printCompressedTable(outputGroupFileName);
            }
            
            for (int i = 0; i < files.size(); i++) { //process each pair
                for (int j = 0; j < files[i].size(); j++) {
                    if (files[i][j] != "") {
                        map<string, vector<string> > filenames = splitFastqFile(outputGroupFileName, files[i][j]);
                        map<string, vector<string> >::iterator it = filenames.find("fastq");
                        if (it != filenames.end()) { //fastq files were produced
                            for (int k = 0; k < it->second.size(); k++) {  outputNames.push_back(it->second[k]); outputTypes["fastq"].push_back(it->second[k]);   }
                        }
                    }
                }
            }
            
            //ffqnoMatchFile, rfqnoMatchFile, ffnoMatchFile, rfnoMatchFile, fqnoMatchFile, rqnoMatchFile
            if(util.isBlank(ffqnoMatchFile)){  util.mothurRemove(ffqnoMatchFile); }
            else { outputNames.push_back(ffqnoMatchFile); outputTypes["fastq"].push_back(ffqnoMatchFile); }
            
            if(fasta){
                if(util.isBlank(ffnoMatchFile)){  util.mothurRemove(ffnoMatchFile); }
                else { outputNames.push_back(ffnoMatchFile); outputTypes["fasta"].push_back(ffnoMatchFile); }
            }
            
            if(qual){
                if(util.isBlank(fqnoMatchFile)){  util.mothurRemove(fqnoMatchFile); }
                else { outputNames.push_back(fqnoMatchFile); outputTypes["qfile"].push_back(fqnoMatchFile); }
            }
            
            if (pairedOligos) {
                if (fileOption > 0) {
                    if(util.isBlank(rfqnoMatchFile)){  util.mothurRemove(rfqnoMatchFile); }
                    else { outputNames.push_back(rfqnoMatchFile); outputTypes["fastq"].push_back(rfqnoMatchFile); }
                    
                    if(fasta){
                        if(util.isBlank(rfnoMatchFile)){  util.mothurRemove(rfnoMatchFile); }
                        else { outputNames.push_back(rfnoMatchFile); outputTypes["fasta"].push_back(rfnoMatchFile); }
                    }
                    
                    if(qual){
                        if(util.isBlank(rqnoMatchFile)){  util.mothurRemove(rqnoMatchFile); }
                        else { outputNames.push_back(rqnoMatchFile); outputTypes["qfile"].push_back(rqnoMatchFile); }
                    }
                }
            }
        }

        //output group counts
        int total = 0;
        if (groupCounts.size() != 0) {  m->mothurOut("\nGroup count: \n");  }
        for (map<string, long long>::iterator it = groupCounts.begin(); it != groupCounts.end(); it++) { total += it->second; m->mothurOut(it->first + "\t" + toString(it->second) + "\n"); }
        if (total != 0) { m->mothurOut("\nTotal of all groups is " + toString(total) + "\n"); }
        
        if (groupfile != "")        { delete groupMap;      }
        if (oligosfile != "")  { delete trimOligos; if (reorient) { delete rtrimOligos; }   }
        
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); }  outputTypes.clear(); outputNames.clear();  return 0; }
		
		//set fasta file as new current fastafile
		string currentName = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
		}
        
        itTypes = outputTypes.find("count");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
        }
		
		itTypes = outputTypes.find("qfile");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setQualFile(currentName); }
		}		
		
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/
map<string, vector<string> > ParseFastaQCommand::splitFastqFile(string outputGroupFile, string resultFastqfile) {
    try {
        
        //run split.groups command
        //use unique.seqs to create new name and fastafile
        string inputString = "fastq=" + resultFastqfile + ", count=" + outputGroupFile;
        m->mothurOut("/******************************************/\n");
        m->mothurOut("Generating parsed files... Running command: split.groups(" + inputString + ")\n");
        current->setMothurCalling(true);
        
        Command* splitCommand = new SplitGroupCommand(inputString);
        splitCommand->execute();
        
        map<string, vector<string> > filenames = splitCommand->getOutputFiles();
        
        delete splitCommand;
        current->setMothurCalling(false);
        m->mothurOut("/******************************************/\n");
        
        if (fasta || qual) { //do we need to create a fasta and qual file for these split fastq files
            
            string fastaBool = "false"; string qualBool = "false";
            if (fasta)  { fastaBool = "true"; }
            if (qual)   { qualBool = "true";  }
            
            map<string, vector<string> >::iterator it = filenames.find("fastq");
            if (it != filenames.end()) { //fastq files were produced
                for (int k = 0; k < it->second.size(); k++) {
                    
                    string inputString = "fastq=" + it->second[k] + ", fasta=" + fastaBool + ", qfile=" + qualBool;
                    m->mothurOut("/******************************************/\n");
                    m->mothurOut("Generating parsed fasta and qual files... Running command: fastq.info(" + inputString + ")\n");
                    current->setMothurCalling(true);
                    
                    Command* fastqCommand = new ParseFastaQCommand(inputString);
                    fastqCommand->execute();
                    
                    map<string, vector<string> > fnames = fastqCommand->getOutputFiles();
                    
                    delete fastqCommand;
                    current->setMothurCalling(false);
                    m->mothurOut("/******************************************/\n");

                    
                    if (fasta) {
                        map<string, vector<string> >::iterator itFastaName = fnames.find("fasta");
                        if (itFastaName != fnames.end()) {
                            string fName = itFastaName->second[0];
                            outputNames.push_back(fName); outputTypes["fasta"].push_back(fName);
                        }
                    }
                    
                    if (qual) {
                        map<string, vector<string> >::iterator itQualName = fnames.find("qfile");
                        if (itQualName != fnames.end()) {
                            string qName = itQualName->second[0];
                            outputNames.push_back(qName); outputTypes["qfile"].push_back(qName);
                        }
                    }
                }
            }
        }
        
        return filenames;
    }
    catch(exception& e) {
        m->errorOut(e, "ParseFastaQCommand", "splitFastqFile");
        exit(1);
    }
}
//**********************************************************************************************************************
//assumes file option was used.
//Adds reads to seqGroup and groupCounts for use with split.groups command later.
//Outputs fasta and qual files for file pair if desired
//Appends scrap files if needed
int ParseFastaQCommand::processFile(vector<string> files, TrimOligos*& trimOligos, TrimOligos*& rtrimOligos){
    try {
        string inputfile = files[0]; string inputReverse = files[1];
        
        //open Output Files
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(inputfile));
        string ffastaFile = getOutputFileName("fasta",variables);
        string fqualFile = getOutputFileName("qfile",variables);
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(inputReverse));
        string rfastaFile = getOutputFileName("fasta",variables);
        string rqualFile = getOutputFileName("qfile",variables);
        ofstream outfFasta, outfQual, outrFasta, outrQual;
        
        if (fasta) { util.openOutputFile(ffastaFile, outfFasta);  outputNames.push_back(ffastaFile); outputTypes["fasta"].push_back(ffastaFile);	util.openOutputFile(rfastaFile, outrFasta);  outputNames.push_back(rfastaFile); outputTypes["fasta"].push_back(rfastaFile);}
        if (qual) { util.openOutputFile(fqualFile, outfQual);	outputNames.push_back(fqualFile);  outputTypes["qfile"].push_back(fqualFile);	util.openOutputFile(rqualFile, outrQual);	outputNames.push_back(rqualFile);  outputTypes["qfile"].push_back(rqualFile);	}
        
        ifstream inf; util.openInputFile(inputfile, inf);
        
        ifstream inr; util.openInputFile(inputReverse, inr);
        
        ifstream inFIndex, inRIndex;
        if (files[2] != "") { util.openInputFile(files[2], inFIndex);  }
        if (files[3] != "") { util.openInputFile(files[3], inRIndex);  }
        
        int count = 0;
        while (!inf.eof() && !inr.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            bool ignoref, ignorer;
            FastqRead thisfRead(inf, ignoref, format);
            FastqRead thisrRead(inr, ignorer, format);
            
            if (!ignoref && ! ignorer) {
                if (qual) {
                    thisfRead.getQuality().printQScores(outfQual);
                    thisrRead.getQuality().printQScores(outrQual);
                }
                
                if (pacbio) { //change sequence bases with 0 quality scores to N
                    vector<int> fqual = thisfRead.getScores();
                    vector<int> rqual = thisrRead.getScores();
                    string fseq = thisfRead.getSeq();
                    string rseq = thisrRead.getSeq();
                    
                    for (int i = 0; i < fqual.size(); i++) { if (fqual[i] == 0){ fseq[i] = 'N'; } }
                    thisfRead.setSeq(fseq);
                
                    for (int i = 0; i < rqual.size(); i++) { if (rqual[i] == 0){ rseq[i] = 'N'; } }
                    thisrRead.setSeq(rseq);
                }
                FastqRead copyForward = thisfRead;
                FastqRead copyReverse = thisrRead;
                
                //print sequence info to files
                if (fasta) {
                    thisfRead.getSequence().printSequence(outfFasta);
                    thisrRead.getSequence().printSequence(outrFasta);
                }
                
                if (m->getControl_pressed()) { break; }
                
                if (split > 1) {
                    
                    Sequence findexBarcode("findex", "NONE");  Sequence rindexBarcode("rindex", "NONE");
                    if (fileOption == 4) {
                        bool ignorefi, ignoreri;
                    
                        if (files[2] != "") {
                            FastqRead thisfiRead(inFIndex, ignorefi, format);
                            if (!ignorefi) {  findexBarcode.setAligned(thisfiRead.getSequence().getAligned());  }
                        }
                        
                        if (files[3] != "") {
                            FastqRead thisriRead(inRIndex, ignoreri, format);
                            if (!ignoreri) {  rindexBarcode.setAligned(thisriRead.getSequence().getAligned());  }
                        }
                    }
                    
                    int trashCodeLength; string thisGroup = "ignore";
                    if (oligosfile != "") {
                        QualityScores tempF = thisfRead.getQuality();
                        QualityScores tempR = thisrRead.getQuality();
                        if ((files[2] != "") || (files[3] != "")) { //has index files
                            //barcode already removed so no need to reset sequence to trimmed version
                            trashCodeLength = findGroup(findexBarcode, tempF, rindexBarcode, tempR, thisGroup, trimOligos, rtrimOligos, numBarcodes, numPrimers);
                        }else {
                            Sequence tempSeqF = thisfRead.getSequence();
                            Sequence tempSeqR = thisrRead.getSequence();
                            trashCodeLength = findGroup(tempSeqF, tempF, tempSeqR, tempR, thisGroup, trimOligos, rtrimOligos, numBarcodes, numPrimers);
                            thisfRead.setSeq(tempSeqF.getUnaligned());
                            thisrRead.setSeq(tempSeqR.getUnaligned());
                        }
                        thisfRead.setScores(tempF.getScores()); //set to trimmed scores
                        thisrRead.setScores(tempR.getScores());
                    }else if (groupfile != "")  {  trashCodeLength = findGroup(thisfRead.getSequence(), thisGroup, "groupMode");   }
                    else {  m->mothurOut("[ERROR]: uh oh, we shouldn't be here...\n"); }
                    
                    bool addToScrap = false;
                    if(trashCodeLength == 0){
                        int pos = thisGroup.find("ignore");
                        if (pos == string::npos) {
                            if (thisGroup != "") {
                                seqGroups[copyForward.getName()] = thisGroup;
                                
                                map<string, long long>::iterator it = groupCounts.find(thisGroup);
                                if (it == groupCounts.end()) {	groupCounts[thisGroup] = 1; }
                                else { groupCounts[it->first]++; }
                            }
                        }else { addToScrap = true; }

                    }else{ addToScrap = true; }
                    
                    if (addToScrap) {
                        //print no match fastq
                        ofstream out, out2;
                        util.openOutputFileAppend(ffqnoMatchFile, out);
                        copyForward.printFastq(out);
                        out.close();
                        
                        util.openOutputFileAppend(rfqnoMatchFile, out2);
                        copyReverse.printFastq(out2);
                        out2.close();
                        
                        //print no match fasta, if wanted
                        if (fasta) {
                            ofstream outf, outr;
                            util.openOutputFileAppend(ffnoMatchFile, outf);
                            thisfRead.getSequence().printSequence(outf);
                            outf.close();
                            
                            util.openOutputFileAppend(rfnoMatchFile, outr);
                            thisrRead.getSequence().printSequence(outr);
                            outr.close();
                        }
                        
                        //print no match quality parse, if wanted
                        if (qual) {
                            ofstream outq, outq2;
                            util.openOutputFileAppend(fqnoMatchFile, outq);
                            thisfRead.getQuality().printQScores(outq);
                            outq.close();
                            
                            util.openOutputFileAppend(rqnoMatchFile, outq2);
                            thisrRead.getQuality().printQScores(outq2);
                            outq2.close();
                        }
                    }
                }
                //report progress
                if((count+1) % 10000 == 0){	m->mothurOut(toString(count+1)); m->mothurOutEndLine();		}
                count++;
            }
        }
        
        inf.close(); inr.close();
        if (files[2] != "") { inFIndex.close();  }
        if (files[3] != "") { inRIndex.close();  }
        
        if (fasta)	{ outfFasta.close(); outrFasta.close();	}
        if (qual)	{ outfQual.close();	outrQual.close();   }
        
        //report progress
        if (!m->getControl_pressed()) {   if((count) % 10000 != 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}  }
        
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "processFile");
		exit(1);
	}
}
//**********************************************************************************************************************
set<string> ParseFastaQCommand::processFile(string inputfile, TrimOligos*& trimOligos, TrimOligos*& rtrimOligos){
    try {
        //fill convert table - goes from solexa to sanger. Used fq_all2std.pl as a reference.
        for (int i = -64; i < 65; i++) {
            char temp = (char) ((int)(33 + 10*log(1+pow(10,(i/10.0)))/log(10)+0.499));
            convertTable.push_back(temp);
        }
        
        //open Output Files
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(inputfile));
        string fastaFile = getOutputFileName("fasta",variables);
        string qualFile = getOutputFileName("qfile",variables);
        ofstream outFasta, outQual;
        
        //fasta and quality files for whole input file
        if (fasta)  { util.openOutputFile(fastaFile, outFasta);  outputNames.push_back(fastaFile); outputTypes["fasta"].push_back(fastaFile);       }
        if (qual)   { util.openOutputFile(qualFile, outQual);	outputNames.push_back(qualFile);  outputTypes["qfile"].push_back(qualFile);         }
        
        ifstream in; util.openInputFile(inputfile, in);
        
        int count = 0;
        set<string> names;
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            bool ignore;
            FastqRead thisRead(in, ignore, format);
            
            if (!ignore) {
                if (qual) {  thisRead.getQuality().printQScores(outQual); }
                
                if (pacbio) {
                    vector<int> qual = thisRead.getScores();
                    string seq = thisRead.getSeq();
                    
                    for (int i = 0; i < qual.size(); i++) { if (qual[i] == 0){ seq[i] = 'N'; } }
                    thisRead.setSeq(seq);
                    names.insert(thisRead.getName());
                }
                
                FastqRead copy = thisRead;
                
                //print sequence info to files
                if (fasta) { thisRead.getSequence().printSequence(outFasta); }
                
                if (m->getControl_pressed()) { break; }
                
                if (split > 1) {
                    int trashCodeLength = 0; string thisGroup = "ignore";
                    if (oligosfile != "")      {
                        Sequence tempSeq = thisRead.getSequence();
                        QualityScores tempQual = thisRead.getQuality();
                        trashCodeLength = findGroup(tempSeq, tempQual, thisGroup, trimOligos, rtrimOligos, numBarcodes, numPrimers);
                        thisRead.setSeq(tempSeq.getUnaligned());
                        thisRead.setScores(tempQual.getScores());
                    }
                    else if (groupfile != "")  {  trashCodeLength = findGroup(thisRead.getSequence(), thisGroup, "groupMode");   }
                    else {  m->mothurOut("[ERROR]: uh oh, we shouldn't be here...\n"); }
                    
                    bool addToScrap = false;
                    if(trashCodeLength == 0){
                        int pos = thisGroup.find("ignore");
                        if (pos == string::npos) {
                            if (thisGroup != "") {
                                seqGroups[copy.getName()] = thisGroup;
                                
                                map<string, long long>::iterator it = groupCounts.find(thisGroup);
                                if (it == groupCounts.end()) {	groupCounts[thisGroup] = 1; }
                                else { groupCounts[it->first]++; }
                            }
                        }else { addToScrap = true; }
                        
                    }else{ addToScrap = true; }
                    
                    if (addToScrap) {
                        //print no match fastq
                        ofstream out;
                        util.openOutputFileAppend(ffqnoMatchFile, out);
                        copy.printFastq(out);
                        out.close();
                        
                        //print no match fasta, if wanted
                        if (fasta) {
                            ofstream outf;
                            util.openOutputFileAppend(ffnoMatchFile, outf);
                            thisRead.getSequence().printSequence(outf);
                            outf.close();
                        }
                        
                        //print no match quality parse, if wanted
                        if (qual) {
                            ofstream outq;
                            util.openOutputFileAppend(fqnoMatchFile, outq);
                            thisRead.getQuality().printQScores(outq);
                            outq.close();
                        }
                    }
                }
                //report progress
                if((count+1) % 10000 == 0){	m->mothurOut(toString(count+1)); m->mothurOutEndLine();		}
                count++;
            }
        }
        
        in.close();
        if (fasta)	{ outFasta.close();	}
        if (qual)	{ outQual.close();	}
        
        //report progress
        if (!m->getControl_pressed()){   if((count) % 10000 != 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}  }
        
        return names;
    }
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "processFile");
		exit(1);
	}
}
//**********************************************************************************************************************
int ParseFastaQCommand::findGroup(Sequence& currSeq, QualityScores& currQual, string& thisGroup, TrimOligos*& trimOligos, TrimOligos*& rtrimOligos, int numBarcodes, int numPrimers) {
	try {
        int success = 1; int barcode, primer;
        string trashCode = "";
        int currentSeqsDiffs = 0;
        
        
        //for reorient
        Sequence savedSeq(currSeq.getName(), currSeq.getAligned());
        QualityScores savedQual(currQual.getName(), currQual.getScores());
        
        if(numLinkers != 0){
            success = trimOligos->stripLinker(currSeq, currQual);
            if(success > ldiffs)		{	trashCode += 'k';	}
            else{ currentSeqsDiffs += success;  }
            
        }
        
        if(numBarcodes != 0){
            vector<int> results = trimOligos->stripBarcode(currSeq, currQual, barcode);
            if (pairedOligos)   {  success = results[0] + results[2];   }
            else                { success = results[0];                 }
            if(success > bdiffs)		{	trashCode += 'b';	}
            else{ currentSeqsDiffs += success;  }
        }
        
        if(numSpacers != 0){
            success = trimOligos->stripSpacer(currSeq, currQual);
            if(success > sdiffs)		{	trashCode += 's';	}
            else{ currentSeqsDiffs += success;  }
            
        }
        
        if(numPrimers != 0){
            vector<int> results = trimOligos->stripForward(currSeq, currQual, primer, true);
            if (pairedOligos)   {  success = results[0] + results[2];   }
            else                {  success = results[0];                }
            if(success > pdiffs)		{	trashCode += 'f';	}
            else{ currentSeqsDiffs += success;  }
        }
        
        if(numRPrimers != 0){
            vector<int> results = trimOligos->stripReverse(currSeq, currQual);
            success = results[0];
            if(success > pdiffs)		{	trashCode += 'r';	}
            else{ currentSeqsDiffs += success;  }
        }
        
        if (currentSeqsDiffs > tdiffs)	{	trashCode += 't';   }
        
        if (reorient && (trashCode != "")) { //if you failed and want to check the reverse
            int thisSuccess = 0;
            string thisTrashCode = "";
            int thisCurrentSeqsDiffs = 0;
            
            int thisBarcodeIndex = 0;
            int thisPrimerIndex = 0;
            
            if(numBarcodes != 0){
                vector<int> results = rtrimOligos->stripBarcode(savedSeq, savedQual, thisBarcodeIndex);
                if (pairedOligos)   {  thisSuccess = results[0] + results[2];   }
                else                {  thisSuccess = results[0];                }
                if(thisSuccess > bdiffs)		{ thisTrashCode += "b"; }
                else{ thisCurrentSeqsDiffs += thisSuccess;  }
            }
            
            if(numPrimers != 0){
                vector<int> results = rtrimOligos->stripForward(savedSeq, savedQual, thisPrimerIndex, true);
                if (pairedOligos)   {  thisSuccess = results[0] + results[2];   }
                else                {  thisSuccess = results[0];                }
                if(thisSuccess > pdiffs)		{ thisTrashCode += "f"; }
                else{ thisCurrentSeqsDiffs += thisSuccess;  }
            }
            
            if (thisCurrentSeqsDiffs > tdiffs)	{	thisTrashCode += 't';   }
            
            if (thisTrashCode == "") {
                trashCode = thisTrashCode;
                success = thisSuccess;
                currentSeqsDiffs = thisCurrentSeqsDiffs;
                barcode = thisBarcodeIndex;
                primer = thisPrimerIndex;
            }else { trashCode += "(" + thisTrashCode + ")";  }
        }
        
        if (trashCode.length() == 0) { //is this sequence in the ignore group
            thisGroup = oligos.getGroupName(barcode, primer);
            
            int pos = thisGroup.find("ignore");
            if (pos != string::npos) {  trashCode += "i"; }
        }

        
        return trashCode.length();
    }
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "findGroup");
		exit(1);
	}
}
//**********************************************************************************************************************
int ParseFastaQCommand::findGroup(Sequence seq, string& group, string groupMode) {
	try {
        string trashCode = "";
        
        group = groupMap->getGroup(seq.getName());
        if (group == "not found") {     trashCode += "g";   } //scrap for group
    
        return trashCode.length();
    }
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "findGroup");
		exit(1);
	}
}
//**********************************************************************************************************************
int ParseFastaQCommand::findGroup(Sequence& fcurrSeq, QualityScores& fcurrQual, Sequence& rcurrSeq, QualityScores& rcurrQual, string& thisGroup, TrimOligos*& trimOligos, TrimOligos*& rtrimOligos, int numBarcodes, int numPrimers) {
	try {
        int success = 1; int barcode, primer;
        string trashCode = "";
        int currentSeqsDiffs = 0;
        
        //for reorient
        Sequence fsavedSeq(fcurrSeq.getName(), fcurrSeq.getAligned());
        QualityScores fsavedQual(fcurrQual.getName(), fcurrQual.getScores());
        Sequence rsavedSeq(rcurrSeq.getName(), rcurrSeq.getAligned());
        QualityScores rsavedQual(rcurrQual.getName(), rcurrQual.getScores());
        
        if(numBarcodes != 0){
            vector<int> results = trimOligos->stripBarcode(fcurrSeq, rcurrSeq, fcurrQual, rcurrQual, barcode);
            if (pairedOligos)   {  success = results[0] + results[2];   }
            else                { success = results[0];                 }
            if(success > bdiffs)		{	trashCode += 'b';	}
            else{ currentSeqsDiffs += success;  }
        }
        
        if(numPrimers != 0){
            vector<int> results = trimOligos->stripForward(fcurrSeq, rcurrSeq, fcurrQual, rcurrQual, primer);
            if (pairedOligos)   {  success = results[0] + results[2];   }
            else                { success = results[0];                 }
            if(success > pdiffs)		{	trashCode += 'f';	}
            else{ currentSeqsDiffs += success;  }
        }
        
        if (currentSeqsDiffs > tdiffs)	{	trashCode += 't';   }
        
        if (reorient && (trashCode != "")) { //if you failed and want to check the reverse
            int thisSuccess = 0;
            string thisTrashCode = "";
            int thisCurrentSeqsDiffs = 0;
            
            int thisBarcodeIndex = 0;
            int thisPrimerIndex = 0;
            
            if(numBarcodes != 0){
                
                vector<int> results = rtrimOligos->stripBarcode(fsavedSeq, rsavedSeq, fsavedQual, rsavedQual, thisBarcodeIndex);
                if (pairedOligos)   {  thisSuccess = results[0] + results[2];   }
                else                {  thisSuccess = results[0];                 }
                
                if(thisSuccess > bdiffs)		{	thisTrashCode += 'b';	}
                else{ thisCurrentSeqsDiffs += thisSuccess;  }
            }
            
            if(numPrimers != 0){
                
                vector<int> results = rtrimOligos->stripForward(fsavedSeq, rsavedSeq, fsavedQual, rsavedQual, thisPrimerIndex);
                if (pairedOligos)   {  thisSuccess = results[0] + results[2];   }
                else                {  thisSuccess = results[0];                 }
                
                if(thisSuccess > pdiffs)		{	thisTrashCode += 'f';	}
                else{ thisCurrentSeqsDiffs += thisSuccess;  }
            }
            
            if (thisCurrentSeqsDiffs > tdiffs)	{	thisTrashCode += 't';   }
            
            if (thisTrashCode == "") {
                trashCode = thisTrashCode;
                success = thisSuccess;
                currentSeqsDiffs = thisCurrentSeqsDiffs;
                barcode = thisBarcodeIndex;
                primer = thisPrimerIndex;
            }else { trashCode += "(" + thisTrashCode + ")";  }
        }
        
        if (trashCode.length() == 0) { //is this sequence in the ignore group
            thisGroup = oligos.getGroupName(barcode, primer);
            
            int pos = thisGroup.find("ignore");
            if (pos != string::npos) {  trashCode += "i"; }
        }
        
        return trashCode.length();
    }
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "findGroup");
		exit(1);
	}
}
//**********************************************************************************************************************
/*
 
 file option 1
 
 ffastqfile1 rfastqfile1
 ffastqfile2 rfastqfile2
 ...
 
 file option 2
 
 group ffastqfile  rfastqfile
 group ffastqfile  rfastqfile
 group ffastqfile  rfastqfile
 ...
 
 file option 3 
 
 My.forward.fastq My.reverse.fastq none My.rindex.fastq //none is an option is no forward or reverse index file
 
 */

//lines can be 2, 3, or 4 columns
// forward.fastq reverse.fastq -> 2 column
// groupName forward.fastq reverse.fastq -> 3 column
// forward.fastq reverse.fastq forward.index.fastq  reverse.index.fastq  -> 4 column
// forward.fastq reverse.fastq none  reverse.index.fastq  -> 4 column
// forward.fastq reverse.fastq forward.index.fastq  none  -> 4 column
vector< vector<string> > ParseFastaQCommand::readFile(){
	try {
        string mode = "parseFastq";
        if (pacbio) { mode = "parsefastqpacbio"; } //reads 2 column option as group filename
        
        FileFile dataFile(inputfile, mode);
        vector< vector<string> > files = dataFile.getFiles(); //if pacbio 2 columns, files[x][0] = filename, files[x][1] = "", files[x][2] = "", files[x][3] = "",
        file2Group = dataFile.getFile2Group();
        createFileGroup = dataFile.isColumnWithGroupNames();
        hasIndex = dataFile.containsIndexFiles();
        int dataFileFormat = dataFile.getFileFormat();
        if (hasIndex && (oligosfile == "")) { m->mothurOut("[ERROR]: You need to provide an oligos file if you are going to use an index file.\n"); m->setControl_pressed(true);  }
        if ((oligosfile != "") && (dataFileFormat == 2)) { m->mothurOut("[ERROR]: You cannot have an oligosfile and 3 column file option at the same time. Aborting. \n"); m->setControl_pressed(true); }
        if ((oligosfile != "") && (dataFileFormat == 1) && pacbio) { m->mothurOut("[ERROR]: You cannot have an oligosfile and 2 column pacbio file option at the same time. Aborting. \n"); m->setControl_pressed(true); }
        if ((groupfile != "")  && (dataFileFormat == 2)){ m->mothurOut("[ERROR]: You cannot have an groupfile and 3 column file option at the same time. Aborting. \n"); m->setControl_pressed(true); }

        for (int i = 0; i < files.size(); i++) {
            string group = "";
            string forward, reverse, findex, rindex;
            forward = files[i][0]; reverse = files[i][1]; findex = files[i][2]; rindex = files[i][3];
            
            if (dataFileFormat == 1) { //2 column
                fileOption = 2;
            }else if (dataFileFormat == 2) { //3 column
                fileOption = 3;
            }else if (dataFileFormat == 3) { //4 column
                fileOption = 4;
                if ((findex == "none") || (findex == "NONE")){ files[i][2] = ""; }
                if ((rindex == "none") || (rindex == "NONE")){ files[i][3] = ""; }
            }
        }
        
        if (files.size() == 0) { m->setControl_pressed(true); }
        
        return files;
    }
    catch(exception& e) {
        m->errorOut(e, "ParseFastaQCommand", "readFileNames");
        exit(1);
    }
}
//***************************************************************************************************************

bool ParseFastaQCommand::readOligos(string oligoFile){
	try {
        bool allBlank = false;
        
        if (fileOption > 0) { oligos.read(oligosfile, false);  } // like make.contigs
        else {  oligos.read(oligosfile);  }
        
        if (m->getControl_pressed()) { return false; } //error in reading oligos
        
        if (oligos.hasPairedPrimers() || oligos.hasPairedBarcodes()) {
            pairedOligos = true;
            numPrimers = oligos.getPairedPrimers().size();
            numBarcodes = oligos.getPairedBarcodes().size();
        }else {
            pairedOligos = false;
            numPrimers = oligos.getPrimers().size();
            numBarcodes = oligos.getBarcodes().size();
        }
        
        numLinkers = oligos.getLinkers().size();
        numSpacers = oligos.getSpacers().size();
        numRPrimers = oligos.getReversePrimers().size();
        
        vector<string> groupNames = oligos.getSRAGroupNames();
        if (groupNames.size() == 0) { allBlank = true;  }
        
        if (allBlank) { m->mothurOut("[WARNING]: your oligos file does not contain any group names.  mothur will not create a groupfile.\n"); return false; }
        
        //make blank files for scrap matches
        ofstream temp, tempff, tempfq, rtemp, temprf, temprq;
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(inputfile));
        variables["[sample]"] = "scrap";
        if (fileOption > 0) {  variables["[tag]"] = "forward"; }
        ffqnoMatchFile = getOutputFileName("fastq", variables);
        util.openOutputFile(ffqnoMatchFile, temp);		temp.close();
        
        if (fileOption > 0) {
            variables["[tag]"] = "reverse";
            rfqnoMatchFile = getOutputFileName("fastq", variables);
            util.openOutputFile(rfqnoMatchFile, rtemp);		rtemp.close();
        }
        
        if (fasta) {
            if (fileOption > 0) {  variables["[tag]"] = "forward"; }
            ffnoMatchFile = getOutputFileName("fasta", variables);
            util.openOutputFile(ffnoMatchFile, tempff);		tempff.close();
            
            if (fileOption > 0) {
                variables["[tag]"] = "reverse";
                rfnoMatchFile = getOutputFileName("fasta", variables);
                util.openOutputFile(rfnoMatchFile, temprf);		temprf.close();
            }
        }
        
        if (qual) {
            if (fileOption > 0) {  variables["[tag]"] = "forward"; }
            fqnoMatchFile = getOutputFileName("qfile", variables);
            util.openOutputFile(fqnoMatchFile, tempfq);		tempfq.close();
            
            if (fileOption > 0) {
                variables["[tag]"] = "reverse";
                rqnoMatchFile = getOutputFileName("qfile", variables);
                util.openOutputFile(rqnoMatchFile, temprq);		temprq.close();
            }
        }

       
		return true;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "getOligos");
		exit(1);
	}
}
//***************************************************************************************************************
bool ParseFastaQCommand::readGroup(string groupfile){
	try {
        groupMap = new GroupMap();
        groupMap->readMap(groupfile);
        
        vector<string> groups = groupMap->getNamesOfGroups();
		
		return true;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "readGroup");
		exit(1);
	}
}
//**********************************************************************************************************************



