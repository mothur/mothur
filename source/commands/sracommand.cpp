//
//  sracommand.cpp
//  Mothur
//
//  Created by SarahsWork on 10/28/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "sracommand.h"
#include "sffinfocommand.h"
#include "parsefastaqcommand.h"

//**********************************************************************************************************************
vector<string> SRACommand::setParameters(){
	try {
        CommandParameter psff("sff", "InputTypes", "", "", "sffFastQFile", "sffFastQFile", "none","xml",false,false); parameters.push_back(psff);
        CommandParameter poligos("oligos", "InputTypes", "", "", "oligos", "none", "none","",false,false,true); parameters.push_back(poligos);
        CommandParameter pfile("file", "InputTypes", "", "", "sffFastQFile-oligos", "sffFastQFile", "none","xml",false,false); parameters.push_back(pfile);
		CommandParameter pfastq("fastq", "InputTypes", "", "", "sffFastQFile", "sffFastQFile", "none","xml",false,false); parameters.push_back(pfastq);
        CommandParameter pcontact("project", "InputTypes", "", "", "none", "none", "none","xml",false,true,true); parameters.push_back(pcontact);
        CommandParameter preorient("checkorient", "Boolean", "", "F", "", "", "","",false,false,true); parameters.push_back(preorient);
        CommandParameter pmimark("mimark", "InputTypes", "", "", "none", "none", "none","xml",false,true,true); parameters.push_back(pmimark);
        //choose only one multiple options
        CommandParameter pplatform("platform", "Multiple", "_LS454-ILLUMINA-ION_TORRENT-PACBIO_SMRT", "_LS454", "", "", "","",false,false); parameters.push_back(pplatform);
        CommandParameter pinstrument("instrument", "Multiple", "454_GS-454_GS_20-454_GS_FLX-454_GS_FLX_Titanium-454_GS_Junior-Illumina_Genome_Analyzer-Illumina_Genome_Analyzer_II-Illumina_Genome_Analyzer_IIx-Illumina_HiSeq_2000-Illumina_HiSeq_1000-Illumina_MiSeq-PacBio_RS-Ion_Torrent_PGM-unspecified", "454_GS", "", "", "","",false,false); parameters.push_back(pinstrument);
        CommandParameter plibstrategy("libstrategy", "String", "AMPLICON", "", "", "", "","",false,false); parameters.push_back(plibstrategy);
        CommandParameter pdatatype("datatype", "String", "METAGENOME", "", "", "", "","",false,false); parameters.push_back(pdatatype);
        CommandParameter plibsource("libsource", "String", "METAGENOMIC", "", "", "", "","",false,false); parameters.push_back(plibsource);
        CommandParameter plibselection("libselection", "String", "PCR", "", "", "", "","",false,false); parameters.push_back(plibselection);
        CommandParameter porientation("orientation", "Multiple", "forward-reverse", "forward", "", "", "","",false,false); parameters.push_back(porientation);
        CommandParameter ppdiffs("pdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(ppdiffs);
		CommandParameter pbdiffs("bdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pbdiffs);
        CommandParameter pldiffs("ldiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pldiffs);
		CommandParameter psdiffs("sdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(psdiffs);
        CommandParameter ptdiffs("tdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(ptdiffs);
        
         //every command must have inputdir and outputdir.  This allows mothur users to redirect input and output files.
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SRACommand::getHelpString(){
	try {
		string helpString = "";
		helpString += "The make.sra command creates the necessary files for a NCBI submission. The xml file and individual sff or fastq files parsed from the original sff or fastq file.\n";
		helpString += "The make.sra command parameters are: sff, fastq, file, oligos, project, mimarksfile, pdiffs, bdiffs, ldiffs, sdiffs, tdiffs, checkorient, platform, orientation, libstrategy, datatype, libsource, libselection and instrument.\n";
        helpString += "The sff parameter is used to provide the original sff file.\n";
		helpString += "The fastq parameter is used to provide the original fastq file.\n";
        helpString += "The project parameter is used to provide your project file.\n";
        helpString += "The oligos parameter is used to provide an oligos file to parse your sff or fastq file by. It is required and must contain barcodes and primers, or you must provide a file option. \n";
        helpString += "The mimark parameter is used to provide your mimarks file.  You can create the template for this file using the get.mimarkspackage command.\n";
		helpString += "The file parameter is used to provide a file containing a list of individual fastq or sff files or paired fastq files with a group assignment. File lines can be 2, 3 or 4 columns. The 2 column files are sff file then oligos or fastqfile then oligos or ffastq and rfastq. You may have multiple lines in the file.  The 3 column files are for paired read libraries. The format is groupName, forwardFastqFile reverseFastqFile. Four column files are for inputting file pairs with index files. Example: My.forward.fastq My.reverse.fastq NONE My.rindex.fastq. The keyword NONE can be used when there is not a index file for either the forward or reverse file. \n";
        helpString += "The tdiffs parameter is used to specify the total number of differences allowed in the sequence. The default is pdiffs + bdiffs + sdiffs + ldiffs.\n";
		helpString += "The bdiffs parameter is used to specify the number of differences allowed in the barcode. The default is 0.\n";
		helpString += "The pdiffs parameter is used to specify the number of differences allowed in the primer. The default is 0.\n";
        helpString += "The ldiffs parameter is used to specify the number of differences allowed in the linker. The default is 0.\n";
		helpString += "The sdiffs parameter is used to specify the number of differences allowed in the spacer. The default is 0.\n";
        helpString += "The checkorient parameter will check look for the reverse compliment of the barcode or primer in the sequence. The default is false.\n";
        helpString += "The platform parameter is used to specify platform you are using choices are: _LS454,ILLUMINA,ION_TORRENT,PACBIO_SMRT. Default=_LS454. This is a controlled vocabulary section in the XML file that will be generated.\n";
        helpString += "The orientation parameter is used to specify sequence orientation. Choices are: forward and reverse. Default=forward. This is a controlled vocabulary section in the XML file that will be generated.\n";
        helpString += "The instrument parameter is used to specify instrument. Choices are 454_GS-454_GS_20-454_GS_FLX-454_GS_FLX_Titanium-454_GS_Junior-Illumina_Genome_Analyzer-Illumina_Genome_Analyzer_II-Illumina_Genome_Analyzer_IIx-Illumina_HiSeq_2000-Illumina_HiSeq_1000-Illumina_MiSeq-PacBio_RS-Ion_Torrent_PGM-unspecified. Default=454_GS. This is a controlled vocabulary section in the XML file that will be generated. \n";
        helpString += "The libstrategy parameter is used to specify library strategy. Default=AMPLICON. Choices are AMPLICON,WGA,WGS,WGX,RNA-Seq,miRNA-Seq,WCS,CLONE,POOLCLONE,CLONEEND,FINISHING,ChIP-Seq,MNase-Seq,DNase-Hypersensitivity,Bisulfite-Seq,Tn-Seq,EST,FL-cDNA,CTS,MRE-Seq,MeDIP-Seq,MBD-Seq,OTHER. This is a controlled vocabulary section in the XML file that will be generated.  \n";
        helpString += "The libsource parameter is used to specify library source. Default=METAGENOMIC. Choices are METAGENOMIC,GENOMIC,TRANSCRIPTOMIC,METATRANSCRIPTOMIC,SYNTHETIC,VIRAL_RNA,OTHER. This is a controlled vocabulary section in the XML file that will be generated. \n";
        helpString += "The libselection parameter is used to specify library selection. Default=PCR. Choices are PCR,RANDOM,RANDOM_PCR,RT-PCR,HMPR,MF,CF-S,CF-H,CF-T,CF-M,MDA,MSLL,cDNA,ChIP,MNase,DNAse,Hybrid_Selection,Reduced_Representation,Restriction_Digest,5-methylcytidine_antibody,MBD2_protein_methyl-CpG_binding_domain,CAGE,RACE,size_fractionation,Padlock_probes_capture_method,other,unspecified. This is a controlled vocabulary section in the XML file that will be generated. \n";
        helpString += "The datatype parameter is used to specify datatype. Default=METAGENOME. Choices are METAGENOME,GENOME_SEQUENCING,METAGENOMIC_ASSEMBLY,ASSEMBLY,TRANSCRIPTOME,PROTEOMIC,MAP,CLONE_ENDS,TARGETED_LOCI,RANDOM_SURVEY,EXOME,VARIATION,EPIGENOMICS,PHENOTYPE,GENOTYPE,OTHER. This is a controlled vocabulary section in the XML file that will be generated. \n";
		helpString += "make.sra(sff=sff=GHL4YHV01.sff, GHL4YHV01.oligos, project=test.project, mimark=MIMarksData.txt)\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string SRACommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        if (type == "xml") {  pattern = "[filename],xml"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SRACommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
SRACommand::SRACommand(){
	try {
		abort = true; calledHelp = true;
		setParameters();
        vector<string> tempOutNames;
		outputTypes["xml"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "SRACommand");
		exit(1);
	}
}
//**********************************************************************************************************************
SRACommand::SRACommand(string option)  {
	try {
        abort = false; calledHelp = false; fileOption = 0;
        libLayout = "single"; //controlled vocab
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			//valid paramters for this command
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) {
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
            vector<string> tempOutNames;
            outputTypes["xml"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter
			inputDir = validParameter.validFile(parameters, "inputdir", false);
			if (inputDir == "not found"){	inputDir = "";		}
			else {
            
                string path;
				it = parameters.find("sff");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["sff"] = inputDir + it->second;		}
				}
				
				it = parameters.find("fastq");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fastq"] = inputDir + it->second;		}
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
                
                it = parameters.find("project");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["project"] = inputDir + it->second;		}
				}
                
                it = parameters.find("mimark");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["mimark"] = inputDir + it->second;		}
				}
            }
            
			//check for parameters
            fastqfile = validParameter.validFile(parameters, "fastq", true);
			if (fastqfile == "not open") { fastqfile = "";  abort = true; }
			else if (fastqfile == "not found") { fastqfile = ""; }
			
			sfffile = validParameter.validFile(parameters, "sff", true);
			if (sfffile == "not open") {  sfffile = "";  abort = true; }
			else if (sfffile == "not found") { sfffile = ""; }
            
            setOligosParameter = false;
            oligosfile = validParameter.validFile(parameters, "oligos", true);
			if (oligosfile == "not found")      {  oligosfile = "";     }
            else if(oligosfile == "not open")	{	abort = true;		}
			else {	m->setOligosFile(oligosfile); setOligosParameter = true; }
            
            contactfile = validParameter.validFile(parameters, "project", true);
			if (contactfile == "not found")      {	contactfile = ""; m->mothurOut("[ERROR]: You must provide a project file before you can use the sra command."); m->mothurOutEndLine(); abort = true;	}
			else if(contactfile == "not open")	{	abort = true;		}
            
            mimarksfile = validParameter.validFile(parameters, "mimark", true);
			if (mimarksfile == "not found")      {	mimarksfile = ""; m->mothurOut("[ERROR]: You must provide a mimark file before you can use the sra command. You can create a template for this file using the get.mimarkspackage command."); m->mothurOutEndLine(); abort = true;	}
			else if(mimarksfile == "not open")	{	abort = true;		}
            
            file = validParameter.validFile(parameters, "file", true);
			if (file == "not open") {  file = "";  abort = true; }
			else if (file == "not found") { file = ""; }
            else {  fileOption = findFileOption();  }
			
            if ((file == "") && (oligosfile == "")) {
                m->mothurOut("[ERROR]: You must provide an oligos file or file with oligos files in them before you can use the sra command."); m->mothurOutEndLine(); abort = true;
            }
            
			if ((fastqfile == "") && (file == "") && (sfffile == "")) {
                m->mothurOut("[ERROR]: You must provide a file, sff file or fastq file before you can use the sra command."); m->mothurOutEndLine(); abort = true;
            }
            
            //use only one Mutliple type _LS454-ILLUMINA-ION_TORRENT-PACBIO_SMRT
			platform = validParameter.validFile(parameters, "platform", false);         if (platform == "not found") { platform = "_LS454"; }
			if (!checkCasesPlatforms(platform)) { abort = true; } //error message in function
			         
            if (!abort) { //don't check instrument model is platform is bad
                //454_GS-454_GS_20-454_GS_FLX-454_GS_FLX_Titanium-454_GS_Junior-Illumina_Genome_Analyzer-Illumina_Genome_Analyzer_II-Illumina_Genome_Analyzer_IIx-Illumina_HiSeq_2000-Illumina_HiSeq_1000-Illumina_MiSeq-PacBio_RS-Ion_Torrent_PGM-unspecified
                instrumentModel = validParameter.validFile(parameters, "instrument", false);         if (instrumentModel == "not found") { instrumentModel = "454_GS"; }
                if (!checkCasesInstrumentModels(instrumentModel)) { abort = true; } //error message in function
            }
            //turn _ to spaces mothur's work around
            for (int i = 0; i < instrumentModel.length(); i++) { if (instrumentModel[i] == '_') { instrumentModel[i] = ' '; } }
            
            libStrategy = validParameter.validFile(parameters, "libstrategy", false);         if (libStrategy == "not found") { libStrategy = "AMPLICON"; }
            if (!checkCasesLibStrategy(libStrategy)) { abort = true; } //error message in function

            //turn _ to spaces mothur's work around
            for (int i = 0; i < libStrategy.length(); i++) { if (libStrategy[i] == '_') { libStrategy[i] = ' '; }  }
            
            libSource = validParameter.validFile(parameters, "libsource", false);         if (libSource == "not found") { libSource = "METAGENOMIC"; }
            if (!checkCasesLibSource(libSource)) { abort = true; } //error message in function
            
            //turn _ to spaces mothur's work around
            for (int i = 0; i < libSource.length(); i++) { if (libSource[i] == '_') { libSource[i] = ' '; }  }
            
            libSelection = validParameter.validFile(parameters, "libselection", false);         if (libSelection == "not found") { libSelection = "PCR"; }
            if (!checkCasesLibSelection(libSelection)) { abort = true; } //error message in function
            
            //turn _ to spaces mothur's work around
            for (int i = 0; i < libSelection.length(); i++) { if (libSelection[i] == '_') { libSelection[i] = ' '; }  }
            
            dataType = validParameter.validFile(parameters, "datatype", false);         if (dataType == "not found") { dataType = "metagenome"; }
            if (!checkCasesDataType(dataType)) { abort = true; } //error message in function
            
            //turn _ to spaces mothur's work around
            for (int i = 0; i < dataType.length(); i++) { if (dataType[i] == '_') { dataType[i] = ' '; }  }
            
            orientation = validParameter.validFile(parameters, "orientation", false);         if (orientation == "not found") { orientation = "forward"; }
            
            if ((orientation == "forward") || (orientation == "reverse")) {  }
            else {  m->mothurOut("[ERROR]: " + orientation + " is not a valid orientation option. Choices are: forward and reverse.\n"); m->mothurOutEndLine(); abort = true; }

            
            string temp = validParameter.validFile(parameters, "bdiffs", false);		if (temp == "not found"){	temp = "0";		}
			m->mothurConvert(temp, bdiffs);
			
			temp = validParameter.validFile(parameters, "pdiffs", false);		if (temp == "not found"){	temp = "0";		}
			m->mothurConvert(temp, pdiffs);
			
            temp = validParameter.validFile(parameters, "ldiffs", false);		if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, ldiffs);
            
            temp = validParameter.validFile(parameters, "sdiffs", false);		if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, sdiffs);
			
			temp = validParameter.validFile(parameters, "tdiffs", false);		if (temp == "not found") { int tempTotal = pdiffs + bdiffs + ldiffs + sdiffs;  temp = toString(tempTotal); }
			m->mothurConvert(temp, tdiffs);
			
			if(tdiffs == 0){	tdiffs = bdiffs + pdiffs + ldiffs + sdiffs;	}
            
            checkorient = validParameter.validFile(parameters, "checkorient", false);		if (temp == "not found") { temp = "F"; }
            			
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "SRACommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int SRACommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        uniqueNames.insert("scrap");
        
        readContactFile();
        if (m->debug) { m->mothurOut("[DEBUG]: read contact file.\n"); }
        readMIMarksFile();
        if (m->debug) { m->mothurOut("[DEBUG]: read mimarks file.\n"); }
        if (oligosfile != "") { readOligos(); }
        if (m->debug) { m->mothurOut("[DEBUG]: read oligos file.\n"); }
        if (m->control_pressed) { return 0; }
        
        //parse files
        map<string, vector<string> > filesBySample;
        isSFF = false;
        
        if (file != "")             {       readFile(filesBySample);        }
        else if (sfffile != "")     {       parseSffFile(filesBySample);    }
        else if (fastqfile != "")   {       parseFastqFile(filesBySample);  }
        
        //cout << "files by sample size = " << filesBySample.size() << endl;
        //checks groups and files returned from parse - removes any groups that did not get reads assigned to them, orders files.
        checkGroups(filesBySample);
        
        sanityCheckMiMarksGroups();
        
        if (m->debug) { m->mothurOut("[DEBUG]: finished sanity check.\n"); }
        
        //create xml file
        string thisOutputDir = outputDir;
        if (outputDir == "") {  thisOutputDir += m->hasPath(inputfile);  }
		map<string, string> variables;
        variables["[filename]"] = thisOutputDir + "submission.";
        string outputFileName = getOutputFileName("xml", variables);
        outputNames.push_back(outputFileName); outputTypes["xml"].push_back(outputFileName);
        ofstream out;
        m->openOutputFile(outputFileName, out);
        
        string blankFile = thisOutputDir + "submission.ready";
        ofstream outT;
        m->openOutputFile(blankFile, outT); outT.close();
        
        //contacts portion
        ////////////////////////////////////////////////////////
        out << "<Submission>\n";
        out << "\t<Description>\n";
        out << "\t\t<Comment> New Submission. Generated by mothur version " + m->getVersion() + "</Comment> \n";
        out << "\t\t<Submitter user_name=\"" + submissionName + "\"/>\n";
        out << "\t\t<Organization role=\"" + ownership + "\" type=\"" + centerType + "\">\n";
        out << "\t\t<Name>" + centerName + "</Name>\n";
        out << "\t\t<Contact email=\"" + email + "\">\n";
        out << "\t\t\t<Name>\n";
        out << "\t\t\t\t<First>" + firstName + "</First>\n";
        out << "\t\t\t\t<Last>" + lastName + "</Last>\n";
        out << "\t\t\t</Name>\n";
        out << "\t\t</Contact>\n";
        out << "\t\t</Organization>\n";
        out << "\t</Description>\n";
        ////////////////////////////////////////////////////////
        
        
        //bioproject
        ////////////////////////////////////////////////////////
        out << "\t<Action>\n";
        out << "\t\t<AddData target_db=\"BioProject\">\n";
        out << "\t\t\t<Data content_type=\"XML\">\n";
        out << "\t\t\t\t<XmlContent>\n";
        out << "\t\t\t\t\t<Project schema_version=\"2.0\">\n";
        out << "\t\t\t\t\t\t<ProjectID>\n";
        out << "\t\t\t\t\t\t<SPUID spuid_namespace=\"" + centerName + "\">" + projectName + "</SPUID> \n";
        out << "\t\t\t\t\t\t</ProjectID>\n";
        out << "\t\t\t\t\t\t<Descriptor>\n";
        out << "\t\t\t\t\t\t\t<Title>" + projectTitle + " </Title> \n";
        out << "\t\t\t\t\t\t\t<Description><p>" + description + "</p></Description> \n";
        if (website != "") {
            out << "\t\t\t\t\t\t\t<ExternalLink label=\"Website name\">\n";
            out << "\t\t\t\t\t\t\t\t<URL>" + website + "</URL>\n";
            out << "\t\t\t\t\t\t\t</ExternalLink>\n";
        }
        if (Grants.size() != 0) {
            for (int i = 0; i < Grants.size(); i++) {
                out << "\t\t\t\t\t\t\t<Grant GrantId=\"" + Grants[i].grantId + "\">\n";
                out << "\t\t\t\t\t\t\t\t<Agency>" + Grants[i].grantAgency + "</Agency>\n";
                if (Grants[i].grantTitle != "") {
                    out << "\t\t\t\t\t\t\t\t<Title>" + Grants[i].grantTitle + "</Title>\n";
                }
                out << "\t\t\t\t\t\t\t</Grant>\n";
            }
        }
        out << "\t\t\t\t\t\t</Descriptor>\n";
        out << "\t\t\t\t\t\t<ProjectType>\n";
        out << "\t\t\t\t\t\t\t<ProjectTypeSubmission sample_scope=\"eEnvironment\">\n";
        out << "\t\t\t\t\t\t\t\t<IntendedDataTypeSet>\n";
        out << "\t\t\t\t\t\t\t\t\t<DataType>" + dataType + " </DataType> \n";
        out << "\t\t\t\t\t\t\t\t</IntendedDataTypeSet>\n";
        out << "\t\t\t\t\t\t\t</ProjectTypeSubmission>\n";
        out << "\t\t\t\t\t\t</ProjectType>\n";
        out << "\t\t\t\t\t</Project>\n";
        out << "\t\t\t\t</XmlContent>\n";
        out << "\t\t\t</Data>\n";
        out << "\t\t\t<Identifier>\n";
        out << "\t\t\t\t\t\t<SPUID spuid_namespace=\"" + centerName + "\">" + projectName + "</SPUID> \n";
        out << "\t\t\t</Identifier>\n";
        out << "\t\t</AddData>\n";
        out << "\t</Action>\n";
        ////////////////////////////////////////////////////////
        
        //bioSample
        ////////////////////////////////////////////////////////
        for (int i = 0; i < Groups.size(); i++) {
            
            if (m->control_pressed) { break; }
            out << "\t<Action>\n";
            out << "\t\t<AddData target_db=\"BioSample\">\n";
            out << "\t\t\t<Data content_type=\"XML\">\n";
            out << "\t\t\t\t<XmlContent>\n";
            out << "\t\t\t\t\t<BioSample schema_version=\"2.0\">\n";
            out << "\t\t\t\t\t\t<SampleId>\n";
            out << "\t\t\t\t\t\t<SPUID spuid_namespace=\"" + centerName + "\">" + Groups[i] + "</SPUID> \n";
            out << "\t\t\t\t\t\t</SampleId>\n";
            out << "\t\t\t\t\t\t<Organism>\n";
            string organismName = "metagenome";
            map<string, string>::iterator itOrganism = Group2Organism.find(Groups[i]);
            if (itOrganism != Group2Organism.end()) { organismName = itOrganism->second; } //user supplied acceptable organism, so use it.
            out << "\t\t\t\t\t\t\t<OrganismName>" + organismName + "</OrganismName> \n";
            out << "\t\t\t\t\t\t</Organism>\n";
            out << "\t\t\t\t\t\t<Descriptor>\n";
            out << "\t\t\t\t\t\t\t<Title>" + mimarks[Groups[i]]["sample_title"] + "</Title> \n";
            out << "\t\t\t\t\t\t\t<Description><p>" + mimarks[Groups[i]]["description"] + "</p></Description> \n";
            out << "\t\t\t\t\t\t</Descriptor>\n";
            out << "\t\t\t\t\t\t<Package>" + packageType + "</Package>\n";
            out << "\t\t\t\t\t\t<Attributes>\n";
            //add biosample required attributes
            map<string, map<string, string> >:: iterator it = mimarks.find(Groups[i]);
            if (it != mimarks.end()) {
                map<string, string> categories = it->second;
                for (map<string, string>:: iterator it2 = categories.begin(); it2 != categories.end(); it2++) {
                    if (m->control_pressed) { break; }
                    out << "\t\t\t\t\t\t\t<Attribute attribute_name=\"" + it2->first + "\">" + it2->second + "</Attribute>\n";
                }
            }
            out << "\t\t\t\t\t\t</Attributes>\n";
            out << "\t\t\t\t\t</BioSample>\n";
            out << "\t\t\t\t</XmlContent>\n";
            out << "\t\t\t</Data>\n";
            out << "\t\t\t<Identifier>\n";
            out << "\t\t\t\t<SPUID spuid_namespace=\"" + centerName + "\">" + Groups[i] + "</SPUID>\n";
            out << "\t\t\t</Identifier>\n";
            out << "\t\t</AddData>\n";
            out << "\t</Action>\n";
        }
        
        map<string, string>::iterator itGroup;
        //File objects
        ////////////////////////////////////////////////////////
        for (int i = 0; i < Groups.size(); i++) {
            
            vector<string> thisGroupsFiles = filesBySample[Groups[i]];
            string thisGroupsBarcode, thisGroupsPrimer;
            if (libLayout == "paired") {  thisGroupsBarcode = "."; thisGroupsPrimer = "."; }
            else { thisGroupsBarcode = ""; thisGroupsPrimer = ""; }
            
            itGroup = Group2Barcode.find(Groups[i]);
            if (itGroup != Group2Barcode.end()) {
                if (fileOption != 5) { thisGroupsBarcode = itGroup->second;  } //don't include barcodes if using index files.
            }
            
            itGroup = Group2Primer.find(Groups[i]);
            if (itGroup != Group2Primer.end()) { thisGroupsPrimer = itGroup->second;  }
            
            //cout << Groups[i] << '\t' << thisGroupsFiles.size() << endl;
            
            for (int j = 0; j < thisGroupsFiles.size(); j++) {
                string libId = m->getSimpleName(thisGroupsFiles[j]) + "." + Groups[i];
                
                if (m->control_pressed) { break; }
                out << "\t<Action>\n";
                out << "\t\t<AddFiles target_db=\"SRA\">\n";
                if (libLayout == "paired") { //adjust the libID because the thisGroupsFiles[j] contains two filenames
                    vector<string> pieces = m->splitWhiteSpace(thisGroupsFiles[j]);
                    libId = m->getSimpleName(pieces[0]) + "." + Groups[i];
                    out << "\t\t\t<File file_path=\"" + m->getSimpleName(pieces[0]) + "\">\n";
                    out << "\t\t\t\t<DataType>generic-data</DataType> \n";
                    out << "\t\t\t</File>\n";
                    out << "\t\t\t<File file_path=\"" + m->getSimpleName(pieces[1]) + "\">\n";
                    out << "\t\t\t\t<DataType>generic-data</DataType> \n";
                    out << "\t\t\t</File>\n";
                    
                    //attributes
                    if (linkers.size() != 0) {
                        string linkerString = "";
                        //linker size forced to 1
                        for (int k = 0; k < linkers.size(); k++) {  linkerString += linkers[k] + ";"; }  linkerString = linkerString.substr(0, linkerString.length()-1);
                        out << "\t\t\t<Attribute name=\"Linker\">" + linkerString + "</Attribute>\n";
                        out << "\t\t\t<Attribute name=\"Linker_max_mismatch\">" + toString(ldiffs) + "</Attribute>\n";
                    }
                    
                    if (thisGroupsBarcode != ".") {
                        string barcodeString = "";
                        vector<string> thisBarcodes; m->splitAtChar(thisGroupsBarcode, thisBarcodes, '.');
                        if (thisBarcodes[0] != "NONE") { barcodeString += thisBarcodes[0] + ";"; }
                        if (thisBarcodes[1] != "NONE") { barcodeString += thisBarcodes[1] + ";"; }//forward barcode + reverse barcode
                        
                        barcodeString = barcodeString.substr(0, barcodeString.length()-1);
                        out << "\t\t\t<Attribute name=\"BarCode\">" + barcodeString + "</Attribute>\n";
                        out << "\t\t\t<Attribute name=\"BarCode_max_mismatch\">" + toString(bdiffs) + "</Attribute>\n";
                    }
                    if (spacers.size() != 0) {
                        string spacerString = "";
                        //spacer size forced to 1
                        for (int k = 0; k < spacers.size(); k++) {  spacerString += spacers[k] + ";"; }  spacerString = spacerString.substr(0, spacerString.length()-1);
                        out << "\t\t\t<Attribute name=\"Adapter\">" + spacerString + "</Attribute>\n";
                        out << "\t\t\t<Attribute name=\"Adapter_max_mismatch\">" + toString(sdiffs) + "</Attribute>\n";
                    }
                    
                    if (thisGroupsPrimer != ".") {
                        string primerString = "";
                        
                        vector<string> thisPrimers; m->splitAtChar(thisGroupsPrimer, thisPrimers, '.');
                        if (thisPrimers[0] != "") { primerString += thisPrimers[0] + ";"; }
                        if (thisPrimers[1] != "") { primerString += thisPrimers[1] + ";"; }
                        
                        if (primerString != "") {
                            primerString = primerString.substr(0, primerString.length()-1);
                            out << "\t\t\t<Attribute name=\"Primer\">" + primerString + "</Attribute>\n";
                            out << "\t\t\t<Attribute name=\"Primer_max_mismatch\">" + toString(pdiffs) + "</Attribute>\n";
                        }
                    }
                    out << "\t\t\t<Attribute name=\"library_name\">" + libId + "</Attribute>\n";
                    out << "\t\t\t<Attribute name=\"library_strategy\">" + libStrategy + "</Attribute>\n";
                    out << "\t\t\t<Attribute name=\"library_source\">" + libSource + "</Attribute>\n";
                    out << "\t\t\t<Attribute name=\"library_selection\">" + libSelection + "</Attribute>\n";
                    out << "\t\t\t<Attribute name=\"library_layout\">" + libLayout + "</Attribute>\n";
                    out << "\t\t\t<Attribute name=\"instrument_model\">" + instrumentModel + "</Attribute>\n";
                    out << "\t\t\t<Attribute name=\"library_construction_protocol\">" + mimarks[Groups[i]]["seq_methods"] + "</Attribute>\n";

                }else { //single
                    out << "\t\t\t<File file_path=\"" + m->getSimpleName(thisGroupsFiles[j]) + "\">\n";
                    out << "\t\t\t\t<DataType>generic-data</DataType> \n";
                    out << "\t\t\t</File>\n";
                    //attributes
                    //linkers -> barcodes -> spacers -> primers
                   
                    if (linkers.size() != 0) {
                        string linkerString = "";
                        for (int k = 0; k < linkers.size(); k++) {  linkerString += linkers[k] + ";"; }  linkerString = linkerString.substr(0, linkerString.length()-1);
                        out << "\t\t\t<Attribute name=\"Linker\">" + linkerString + "</Attribute>\n";
                        out << "\t\t\t<Attribute name=\"Linker_max_mismatch\">" + toString(ldiffs) + "</Attribute>\n";
                    }
                    
                    if (thisGroupsBarcode != "") {
                        out << "\t\t\t<Attribute name=\"BarCode\">" + thisGroupsBarcode + "</Attribute>\n";
                        out << "\t\t\t<Attribute name=\"BarCode_max_mismatch\">" + toString(bdiffs) + "</Attribute>\n";
                    }
                    if (spacers.size() != 0) {
                        string spacerString = "";
                        for (int k = 0; k < spacers.size(); k++) {  spacerString += spacers[k] + ";"; }  spacerString = spacerString.substr(0, spacerString.length()-1);
                        out << "\t\t\t<Attribute name=\"Adapter\">" + spacerString + "</Attribute>\n";
                        out << "\t\t\t<Attribute name=\"Adapter_max_mismatch\">" + toString(sdiffs) + "</Attribute>\n";
                    }
                    
                    if (thisGroupsPrimer != "") {
                        out << "\t\t\t<Attribute name=\"Primer\">" + thisGroupsPrimer + "</Attribute>\n";
                        out << "\t\t\t<Attribute name=\"Primer_max_mismatch\">" + toString(pdiffs) + "</Attribute>\n";
                        
                    }
                    //out << "\t\t\t<Attribute name=\"read_type\">" + orientation + "</Attribute>\n";
                    out << "\t\t\t<Attribute name=\"library_name\">" + libId + "</Attribute>\n";
                    out << "\t\t\t<Attribute name=\"library_strategy\">" + libStrategy + "</Attribute>\n";
                    out << "\t\t\t<Attribute name=\"library_source\">" + libSource + "</Attribute>\n";
                    out << "\t\t\t<Attribute name=\"library_selection\">" + libSelection + "</Attribute>\n";
                    out << "\t\t\t<Attribute name=\"library_layout\">" + libLayout + "</Attribute>\n";
                    out << "\t\t\t<Attribute name=\"instrument_model\">" + instrumentModel + "</Attribute>\n";
                    out << "\t\t\t<Attribute name=\"library_construction_protocol\">" + mimarks[Groups[i]]["seq_methods"] + "</Attribute>\n";

                }
                ///////////////////bioProject info
                out << "\t\t\t<AttributeRefId name=\"BioProject\">\n";
                out << "\t\t\t\t<RefId>\n";
                out << "\t\t\t\t\t<SPUID spuid_namespace=\"" + centerName + "\">" + projectName + "</SPUID> \n";
                out << "\t\t\t\t</RefId>\n";
                out << "\t\t\t</AttributeRefId>\n";
                //////////////////bioSample info
                out << "\t\t\t<AttributeRefId name=\"BioSample\">\n";
                out << "\t\t\t\t<RefId>\n";
                out << "\t\t\t\t\t<SPUID spuid_namespace=\"" + centerName + "\">" + Groups[i] + "</SPUID>\n";
                out << "\t\t\t\t</RefId>\n";
                out << "\t\t\t</AttributeRefId>\n";
                //libID
                out << "\t\t\t<Identifier>\n";
                if (libLayout == "paired") { //adjust the libID because the thisGroupsFiles[j] contains two filenames
                    vector<string> pieces = m->splitWhiteSpace(thisGroupsFiles[j]);
                    libId = m->getSimpleName(pieces[0]) + "." + Groups[i];
                }
                out << "\t\t\t\t<LocalId>" + libId + "</LocalId>\n";
                out << "\t\t\t</Identifier>\n";
                out << "\t\t</AddFiles>\n";
                out << "\t</Action>\n";
            }
        }
        out << "</Submission>\n";
        out.close();
        
        if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } return 0; }
		
        //output files created by command
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
        return 0;
		
    }
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int SRACommand::readContactFile(){
	try {
        lastName = ""; firstName = ""; submissionName = ""; email = ""; centerName = ""; centerType = ""; description = ""; website = ""; projectName = "";
        projectTitle = "";  ownership = "owner";
        
        ifstream in;
        m->openInputFile(contactfile, in);
        
        while(!in.eof()) {
            
            if (m->control_pressed) { break; }
            
            string key, value;
            in >> key; m->gobble(in);
            value = m->getline(in); m->gobble(in);
            
            for (int i = 0; i < key.length(); i++) { key[i] = toupper(key[i]); }
            
            if (key == "USERNAME")          {   submissionName = value; }
            else if (key == "LAST")         {   lastName = value;       }
            else if (key == "FIRST")        {   firstName = value;      }
            else if (key == "EMAIL")        {   email = value;          }
            else if (key == "CENTER")       {   centerName = value;     }
            else if (key == "TYPE")         {
                centerType = value;
                for (int i = 0; i < centerType.length(); i++) { centerType[i] = tolower(centerType[i]); }
                if ((centerType == "consortium") || (centerType == "center") ||  (centerType == "institute") ||  (centerType == "lab")) {}
                else { m->mothurOut("[ERROR]: " + centerType + " is not a center type option.  Valid center type options are consortium, center, institute and lab. This is a controlled vocabulary section in the XML file that will be generated."); m->mothurOutEndLine(); m->control_pressed = true; }
            }else if (key == "OWNERSHIP")         {
                    ownership = value;
                    for (int i = 0; i < centerType.length(); i++) { centerType[i] = tolower(ownership[i]); }
                    if ((ownership == "owner") || (ownership == "participant")) {}
                    else { m->mothurOut("[ERROR]: " + ownership + " is not a ownership option.  Valid ownership options are owner or participant. This is a controlled vocabulary section in the XML file that will be generated."); m->mothurOutEndLine(); m->control_pressed = true; }
            }else if (key == "DESCRIPTION")     {   description = value;    }
            else if (key == "WEBSITE")          {   website = value;        }
            else if (key == "PROJECTNAME")      {   projectName = value;    }
            else if (key == "PROJECTTITLE")     {   projectTitle = value;   }
            else if (key == "GRANT")            {
                string temp = value;
                vector<string> values;
                m->splitAtComma(temp, values);
                Grant thisGrant;
                for (int i = 0; i < values.size(); i++) {
                    vector<string> items;
                    m->splitAtChar(values[i], items, '=');
                    if (items.size() != 2) { m->mothurOut("[ERROR]: error parsing grant info for line \"" + value + "\", skipping it.\n"); break; }
                    else {
                        //remove any leading spaces in tag
                        int i;
                        for (i = 0; i < items[0].length(); i++) {  if (isspace(items[0][i])) {}else {break;} }
                        items[0] = items[0].substr(i);
                        
                        if (items[0] == "id") { thisGrant.grantId = items[1];  }
                        else if (items[0] == "title") { thisGrant.grantTitle = items[1];  }
                        else if (items[0] == "agency") { thisGrant.grantAgency = items[1];  }
                        else { m->mothurOut("[ERROR]: unknown identifier '" + items[0] + "', skipping it.\n");  }
                    }
                }
                if ((thisGrant.grantId == "") || (thisGrant.grantAgency == "")) { m->mothurOut("[ERROR]: Missing info for line \"" + value + "\", skipping it. Note: the id and agency fields are required. Example: Grant id=yourID, agency=yourAgency.\n"); }
                else { Grants.push_back(thisGrant); }
            }
            
        }
        in.close();
        
        if (lastName == "") { m->mothurOut("[ERROR]: missing last name from project file, quitting."); m->mothurOutEndLine(); m->control_pressed = true; }
        if (firstName == "") { m->mothurOut("[ERROR]: missing first name from project file, quitting."); m->mothurOutEndLine(); m->control_pressed = true; }
        if (submissionName == "") { m->mothurOut("[ERROR]: missing submission name from project file, quitting."); m->mothurOutEndLine(); m->control_pressed = true; }
        if (email == "") { m->mothurOut("[ERROR]: missing email from project file, quitting."); m->mothurOutEndLine(); m->control_pressed = true; }
        if (centerName == "") { m->mothurOut("[ERROR]: missing center name from project file, quitting."); m->mothurOutEndLine(); m->control_pressed = true; }
        if (centerType == "") { m->mothurOut("[ERROR]: missing center type from project file, quitting."); m->mothurOutEndLine(); m->control_pressed = true; }
        if (description == "") { m->mothurOut("[ERROR]: missing description from project file, quitting."); m->mothurOutEndLine(); m->control_pressed = true; }
        if (projectTitle == "") { m->mothurOut("[ERROR]: missing project title from project file, quitting."); m->mothurOutEndLine(); m->control_pressed = true; }
        if (projectName == "") { m->mothurOut("[ERROR]: missing project name from project file, quitting."); m->mothurOutEndLine(); m->control_pressed = true; }

        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "readContactFile");
		exit(1);
	}
}
//**********************************************************************************************************************
//air, host_associated, human_associated, human_gut, human_oral, human_skin, human_vaginal, microbial, miscellaneous, plant_associated, sediment, soil, wastewater or water
//all packages require: *sample_name	*organism	*collection_date	*biome	*feature	*material	*geo_loc_name	*lat_lon
//air: *altitude
//host_associated, human_associated, human_gut, human_oral, human_skin, human_vaginal, plant_associated: *host
//microbial, sediment, soil: *depth	*elev
//water: *depth
int SRACommand::readMIMarksFile(){
	try {
        //acceptable organisms
        vector<string> acceptableOrganisms;
        bool organismError = false;
        //ecological
        acceptableOrganisms.push_back("activated carbon metagenome"); acceptableOrganisms.push_back("activated sludge metagenome"); acceptableOrganisms.push_back("air metagenome"); acceptableOrganisms.push_back("anaerobic digester metagenome"); acceptableOrganisms.push_back("ant fungus garden metagenome"); acceptableOrganisms.push_back("aquatic metagenome"); acceptableOrganisms.push_back("activated carbon metagenome"); acceptableOrganisms.push_back("activated sludge metagenome"); acceptableOrganisms.push_back("beach sand metagenome"); acceptableOrganisms.push_back("biofilm metagenome"); acceptableOrganisms.push_back("biofilter metagenome"); acceptableOrganisms.push_back("biogas fermenter metagenome"); acceptableOrganisms.push_back("bioreactor metagenome"); acceptableOrganisms.push_back("bioreactor sludge metagenome"); acceptableOrganisms.push_back("clinical metagenome"); acceptableOrganisms.push_back("coal metagenome"); acceptableOrganisms.push_back("compost metagenome"); acceptableOrganisms.push_back("dust metagenome"); acceptableOrganisms.push_back("fermentation metagenome"); acceptableOrganisms.push_back("food fermentation metagenome"); acceptableOrganisms.push_back("food metagenome"); acceptableOrganisms.push_back("freshwater metagenome"); acceptableOrganisms.push_back("freshwater sediment metagenome"); acceptableOrganisms.push_back("groundwater metagenome"); acceptableOrganisms.push_back("halite metagenome"); acceptableOrganisms.push_back("hot springs metagenome"); acceptableOrganisms.push_back("hydrocarbon metagenome"); acceptableOrganisms.push_back("hydrothermal vent metagenome"); acceptableOrganisms.push_back("hypersaline lake metagenome"); acceptableOrganisms.push_back("ice metagenome"); acceptableOrganisms.push_back("indoor metagenome"); acceptableOrganisms.push_back("industrial waste metagenome"); acceptableOrganisms.push_back("mangrove metagenome"); acceptableOrganisms.push_back("marine metagenome"); acceptableOrganisms.push_back("marine sediment metagenome"); acceptableOrganisms.push_back("microbial mat metagenome"); acceptableOrganisms.push_back("mine drainage metagenome"); acceptableOrganisms.push_back("mixed culture metagenome"); acceptableOrganisms.push_back("oil production facility metagenome"); acceptableOrganisms.push_back("paper pulp metagenome"); acceptableOrganisms.push_back("permafrost metagenome"); acceptableOrganisms.push_back("plastisphere metagenome"); acceptableOrganisms.push_back("power plant metagenome"); acceptableOrganisms.push_back("retting rhizosphere metagenome"); acceptableOrganisms.push_back("rock metagenome"); acceptableOrganisms.push_back("salt lake metagenome"); acceptableOrganisms.push_back("saltern metagenome"); acceptableOrganisms.push_back("sediment metagenome"); acceptableOrganisms.push_back("snow metagenome"); acceptableOrganisms.push_back("soil metagenome"); acceptableOrganisms.push_back("stromatolite metagenome"); acceptableOrganisms.push_back("terrestrial metagenome"); acceptableOrganisms.push_back("tomb wall metagenome"); acceptableOrganisms.push_back("wastewater metagenome"); acceptableOrganisms.push_back("wetland metagenome"); acceptableOrganisms.push_back("whale fall metagenome");
        //oganismal
        acceptableOrganisms.push_back("algae metagenome"); acceptableOrganisms.push_back("ant metagenome"); acceptableOrganisms.push_back("bat metagenome"); acceptableOrganisms.push_back("beetle metagenome"); acceptableOrganisms.push_back("bovine gut metagenome"); acceptableOrganisms.push_back("bovine metagenome"); acceptableOrganisms.push_back("chicken gut metagenome"); acceptableOrganisms.push_back("coral metagenome"); acceptableOrganisms.push_back("echinoderm metagenome"); acceptableOrganisms.push_back("endophyte metagenome"); acceptableOrganisms.push_back("epibiont metagenome"); acceptableOrganisms.push_back("fish metagenome"); acceptableOrganisms.push_back("fossil metagenome"); acceptableOrganisms.push_back("gill metagenome"); acceptableOrganisms.push_back("gut metagenome"); acceptableOrganisms.push_back("honeybee metagenome"); acceptableOrganisms.push_back("human gut metagenome"); acceptableOrganisms.push_back("human lung metagenome"); acceptableOrganisms.push_back("human metagenome"); acceptableOrganisms.push_back("human nasal/pharyngeal metagenome"); acceptableOrganisms.push_back("human oral metagenome"); acceptableOrganisms.push_back("human skin metagenome"); acceptableOrganisms.push_back("insect gut metagenome"); acceptableOrganisms.push_back("insect metagenome"); acceptableOrganisms.push_back("mollusc metagenome"); acceptableOrganisms.push_back("mosquito metagenome"); acceptableOrganisms.push_back("mouse gut metagenome"); acceptableOrganisms.push_back("mouse metagenome"); acceptableOrganisms.push_back("mouse skin metagenome"); acceptableOrganisms.push_back("nematode metagenome"); acceptableOrganisms.push_back("oral metagenome"); acceptableOrganisms.push_back("phyllosphere metagenome"); acceptableOrganisms.push_back("pig metagenome"); acceptableOrganisms.push_back("plant metagenome"); acceptableOrganisms.push_back("primate metagenome"); acceptableOrganisms.push_back("rat metagenome"); acceptableOrganisms.push_back("root metagenome"); acceptableOrganisms.push_back("sea squirt metagenome"); acceptableOrganisms.push_back("seed metagenome"); acceptableOrganisms.push_back("shoot metagenome"); acceptableOrganisms.push_back("skin metagenome"); acceptableOrganisms.push_back("snake metagenome"); acceptableOrganisms.push_back("sponge metagenome"); acceptableOrganisms.push_back("stomach metagenome"); acceptableOrganisms.push_back("symbiont metagenome"); acceptableOrganisms.push_back("termite gut metagenome"); acceptableOrganisms.push_back("termite metagenome"); acceptableOrganisms.push_back("upper respiratory tract metagenome"); acceptableOrganisms.push_back("urine metagenome"); acceptableOrganisms.push_back("viral metagenome"); acceptableOrganisms.push_back("wallaby gut metagenome"); acceptableOrganisms.push_back("wasp metagenome"); acceptableOrganisms.push_back("sythetic metagenome"); acceptableOrganisms.push_back("metagenome");
        
        vector<string> requiredFieldsForPackage;
        requiredFieldsForPackage.push_back("sample_name");
        requiredFieldsForPackage.push_back("description");
        requiredFieldsForPackage.push_back("sample_title");
        requiredFieldsForPackage.push_back("collection_date");
        requiredFieldsForPackage.push_back("env_biome");
        requiredFieldsForPackage.push_back("env_feature");
        requiredFieldsForPackage.push_back("env_material");
        requiredFieldsForPackage.push_back("geo_loc_name");
        requiredFieldsForPackage.push_back("lat_lon");
        requiredFieldsForPackage.push_back("seq_methods");
        requiredFieldsForPackage.push_back("organism");
        
        ifstream in;
        m->openInputFile(mimarksfile, in);
        
        //read comments
        string temp; packageType = "";
        while(!in.eof()) {
            
            if (m->control_pressed) { break; }
            temp = m->getline(in); m->gobble(in);
            
            if (m->debug) { m->mothurOut("[DEBUG]: " + temp + "\n"); }
            
            if (temp[0] == '#') {
                int pos = temp.find("MIMARKS.survey");
                if (pos != string::npos) {
                    packageType = temp.substr(1);
                }
            }
            else{ break; } //hit headers line
         }
        
        //in future may want to add parsing of format header....
        
        vector<string> headers; m->splitAtChar(temp, headers, '\t');
        m->removeBlanks(headers);
        //remove * from required's
        for (int i = 0; i < headers.size(); i++) {
            if (headers[i][0] == '*') { headers[i] = headers[i].substr(1); }
            if (m->debug) { m->mothurOut("[DEBUG]: " + headers[i] + "\n"); }
        }
        
        if (m->debug) {  m->mothurOut("[DEBUG]: packageType = '" + packageType + "'\n");   }
        
        if (packageType == "MIMARKS.survey.air.4.0") {	requiredFieldsForPackage.push_back("altitude");	}
        if (packageType == "MIMARKS.survey.host-associated.4.0") {		requiredFieldsForPackage.push_back("host");	}
        if (packageType == "MIMARKS.survey.human-associated.4.0") {		requiredFieldsForPackage.push_back("host");	}
        if (packageType == "MIMARKS.survey.human-gut.4.0") {		requiredFieldsForPackage.push_back("host");	}
        if (packageType == "MIMARKS.survey.human-oral.4.0") {		requiredFieldsForPackage.push_back("host");	}
        if (packageType == "MIMARKS.survey.human-skin.4.0") {	requiredFieldsForPackage.push_back("host");	}
        if (packageType == "MIMARKS.survey.human-vaginal.4.0") {	requiredFieldsForPackage.push_back("host");	}
        if (packageType == "MIMARKS.survey.microbial.4.0") {	requiredFieldsForPackage.push_back("depth");  requiredFieldsForPackage.push_back("elev");	}
        if (packageType == "MIMARKS.survey.miscellaneous.4.0") {}
        if (packageType == "MIMARKS.survey.plant-associated.4.0") {	requiredFieldsForPackage.push_back("host");	}
        if (packageType == "MIMARKS.survey.sediment.4.0") {	requiredFieldsForPackage.push_back("depth");	requiredFieldsForPackage.push_back("elev");	}
        if (packageType == "MIMARKS.survey.soil.4.0") {	requiredFieldsForPackage.push_back("depth");	requiredFieldsForPackage.push_back("elev");	}
        if (packageType == "MIMARKS.survey.wastewater.4.0") {}
        if (packageType == "MIMARKS.survey.water.4.0") {	requiredFieldsForPackage.push_back("depth"); }

        
        if (!m->isSubset(headers, requiredFieldsForPackage)){
            string requiredFields = "";
            for (int i = 0; i < requiredFieldsForPackage.size()-1; i++) { requiredFields += requiredFieldsForPackage[i] + ", "; } requiredFields += requiredFieldsForPackage[requiredFieldsForPackage.size()-1];
            m->mothurOut("[ERROR]: missing required fields for package, please correct. Required fields are " + requiredFields + ".\n"); m->control_pressed = true; in.close(); return 0;
            
            
        }
        
        //if (m->debug) {  m->mothurOut("[DEBUG]: chooseAtLeastOneForPackage.size() = " + toString(chooseAtLeastOneForPackage.size()) + "\n");   }
        
        //if (!m->inUsersGroups(chooseAtLeastOneForPackage, headers)){ //returns true if any of the choose at least ones are in headers
            //string requiredFields = "";
            //for (int i = 0; i < chooseAtLeastOneForPackage.size()-1; i++) { requiredFields += chooseAtLeastOneForPackage[i] + ", "; cout << chooseAtLeastOneForPackage[i] << endl; }
            //if (chooseAtLeastOneForPackage.size() < 1) { requiredFields += chooseAtLeastOneForPackage[chooseAtLeastOneForPackage.size()-1]; }
            //m->mothurOut("[ERROR]: missing a choose at least one fields for the package, please correct. These are marked with '**'. Required fields are " + requiredFields + ".\n"); m->control_pressed = true; in.close(); return 0;
       // }
        
        map<string, bool> allNA;  for (int i = 1; i < headers.size(); i++) {  allNA[headers[i]] = true; }
        while(!in.eof()) {
            
            if (m->control_pressed) { break; }
            
            temp = m->getline(in);  m->gobble(in);
            //cout << temp << endl;
            if (m->debug) { m->mothurOut("[DEBUG]: " + temp + "\n"); }
            
            string original = temp;
            vector<string> linePieces; m->splitAtChar(temp, linePieces, '\t');
            m->removeBlanks(linePieces);
            
            if (linePieces.size() != headers.size()) { m->mothurOut("[ERROR]: line: " + original + " contains " + toString(linePieces.size()) + " columns, but you have " + toString(headers.size()) + " column headers, please correct.\n"); m->control_pressed = true; }
            else {
                map<string, map<string, string> >:: iterator it = mimarks.find(linePieces[0]);
                
                if (it == mimarks.end()) {
                    map<string, string> categories;
                    //start after *sample_name
                    for (int i = 1; i < headers.size(); i++) {
                        categories[headers[i]] = linePieces[i];
                        //check the users inputs for appropriate organisms
                        if (headers[i] == "organism") {
                            if (!m->inUsersGroups(linePieces[i], acceptableOrganisms)) { //not an acceptable organism
                                organismError = true;
                                m->mothurOut("[WARNING]: " + linePieces[i]+ " is not an acceptable organism, changing to acceptable 'metagenome'. NCBI will allow you to modify the organism after submission.\n"); linePieces[i] = "metagenome"; categories[headers[i]] = linePieces[i];
                            }else {
                                if (linePieces[i] == "metagenome") {
                                    m->mothurOut("[WARNING]: metagenome is an acceptable organism, but NCBI would prefer a more specific choice if possible. Here is a link to the organism choices and descriptions, http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=408169&lvl=3&keep=1&srchmode=1&unlock. To request the addition of a taxonomy to the list, please contact Anjanette Johnston at johnston@ncbi.nlm.nih.gov.\n");
                                }
                            }
                            Group2Organism[linePieces[0]] = linePieces[i];
                        }
                        if (linePieces[i] != "missing") {  allNA[headers[i]] = false;     }
                    }
                    
                    //does this sample already match an existing sample?
                    bool isOkaySample = true;
                    for (map<string, map<string, string> >:: iterator it2 = mimarks.begin(); it2 != mimarks.end(); it2++) {
                        if (m->control_pressed) { break; }
                        bool allSame = true;
                        for (int i = 1; i < headers.size(); i++) {
                            if ((it2->second)[headers[i]] != categories[headers[i]]) { allSame = false; }
                        }
                        if (allSame) { m->mothurOut("[ERROR]: " + linePieces[0]+ " is a duplicate sample to " + it2->first + ". It has all the same attributes in the MIMarks file. Samples must have distinguishing features to be uploaded to the NCBI library, please correct.\n"); m->control_pressed = true; isOkaySample = false; }
                    }
                    if (isOkaySample) { mimarks[linePieces[0]] = categories; }
                }else {
                    m->mothurOut("[ERROR]: " + linePieces[0]+ " is a duplicate sampleName. Sample names must be unique, please correct.\n"); m->control_pressed = true;
                }
            }
        }
        in.close();
        
        //add in values for "scrap" group
        map<string, string> categories;
        //start after *sample_name
        for (int i = 1; i < headers.size(); i++) {
            categories[headers[i]] = "missing";
            if (headers[i] == "organism")       { categories[headers[i]] = "metagenome"; }
            if (headers[i] == "description")    { categories[headers[i]] = "these sequences were scrapped"; }
            if (headers[i] == "sample_title")          { categories[headers[i]] = "these sequences were scrapped"; }
        }
        mimarks["scrap"] = categories;
        Group2Organism["scrap"] = "metagenome";
        
        if (organismError) {
            string organismTypes = "";
            for (int i = 0; i < acceptableOrganisms.size()-1; i++) { organismTypes += acceptableOrganisms[i] + ", "; }
            organismTypes += acceptableOrganisms[acceptableOrganisms.size()-1];
            m->mothurOut("\n[WARNING]: The acceptable organism choices are: " + organismTypes + ".\n\n\n");
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "readMIMarksFile");
		exit(1);
	}
}

//**********************************************************************************************************************

/*
 file option 1
 
 sfffile1   oligosfile1
 sfffile2   oligosfile2
 ...
 
 file option 2
 
 fastqfile1 oligosfile1
 fastqfile2 oligosfile2
 ...
 
 file option 3
 
 ffastqfile1 rfastqfile1
 ffastqfile2 rfastqfile2
 ...
 
 file option 4
 
 group fastqfile  fastqfile
 group fastqfile  fastqfile
 group fastqfile  fastqfile
 ...
 
 file option 5
  
 My.forward.fastq My.reverse.fastq none My.rindex.fastq //none is an option is no forward or reverse index file
 ...
 
 
*/

int SRACommand::readFile(map<string, vector<string> >& files){
	try {
        bool runParseFastqFile = false;
        bool using3NONE = false;
        inputfile = file;
        files.clear();
        
        ifstream in;
        m->openInputFile(file, in);
        
        fileOption = 0;
        
        while(!in.eof()) {
            
            if (m->control_pressed) { return 0; }
            
            string line = m->getline(in);  m->gobble(in);
            vector<string> pieces = m->splitWhiteSpace(line);
            
            string group = "";
            string thisFileName1, thisFileName2, findex, rindex; thisFileName1 = ""; thisFileName2 = ""; findex = ""; rindex = "";
            if (pieces.size() == 2) {
                thisFileName1 = pieces[0];
                thisFileName2 = pieces[1];
            }else if (pieces.size() == 3) {
                thisFileName1 = pieces[1];
                thisFileName2 = pieces[2];
                group = pieces[0];
                if (setOligosParameter) { m->mothurOut("[ERROR]: You cannot have an oligosfile and 3 column file option at the same time. Aborting. \n"); m->control_pressed = true; }
                if ((thisFileName2 != "none") && (thisFileName2 != "NONE" )) {  if (!using3NONE) { libLayout = "paired"; } else { m->mothurOut("[ERROR]: You cannot have a 3 column file with paired and unpaired files at the same time. Aborting. \n"); m->control_pressed = true; } }
                else {  thisFileName2 = ""; libLayout = "single"; using3NONE = true; }
            }else if (pieces.size() == 4) {
                if (!setOligosParameter) { m->mothurOut("[ERROR]: You must have an oligosfile with the index file option. Aborting. \n"); m->control_pressed = true; }
                thisFileName1 = pieces[0];
                thisFileName2 = pieces[1];
                findex = pieces[2];
                rindex = pieces[3];
                if ((findex == "none") || (findex == "NONE")){ findex = ""; }
                if ((rindex == "none") || (rindex == "NONE")){ rindex = ""; }
            }else {
                m->mothurOut("[ERROR]: file lines can be 2, 3 or 4 columns. The 2 column files are sff file then oligos or fastqfile then oligos or ffastq and rfastq. You may have multiple lines in the file.  The 3 column files are for paired read libraries. The format is groupName, forwardFastqFile reverseFastqFile. Four column files are for inputting file pairs with index files. Example: My.forward.fastq My.reverse.fastq NONE My.rindex.fastq. The keyword NONE can be used when there is not a index file for either the forward or reverse file.\n"); m->control_pressed = true;
            }
            
            if (m->debug) { m->mothurOut("[DEBUG]: group = " + group + ", thisFileName1 = " + thisFileName1 + ", thisFileName2 = " + thisFileName2  + ".\n"); }
            
            if (inputDir != "") {
                string path = m->hasPath(thisFileName1);
                if (path == "") {  thisFileName1 = inputDir + thisFileName1;  }
                
                if (thisFileName2 != "") {
                    path = m->hasPath(thisFileName2);
                    if (path == "") {  thisFileName2 = inputDir + thisFileName2;  }
                }
                
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
            int openForward = m->openInputFile(thisFileName1, in2, "noerror");
            
            //if you can't open it, try default location
            if (openForward == 1) {
                
                if (m->getDefaultPath() != "") { //default path is set
                    string tryPath = m->getDefaultPath() + m->getSimpleName(thisFileName1);
                    m->mothurOut("Unable to open " + thisFileName1 + ". Trying default " + tryPath); m->mothurOutEndLine();
                    ifstream in3;
                    openForward = m->openInputFile(tryPath, in3, "noerror");
                    in3.close();
                    thisFileName1 = tryPath;
                }
            }
            
            //if you can't open it, try output location
            if (openForward == 1) {
                if (m->getOutputDir() != "") { //default path is set
                    string tryPath = m->getOutputDir() + m->getSimpleName(thisFileName1);
                    m->mothurOut("Unable to open " + thisFileName1 + ". Trying output directory " + tryPath); m->mothurOutEndLine();
                    ifstream in4;
                    openForward = m->openInputFile(tryPath, in4, "noerror");
                    thisFileName1 = tryPath;
                    in4.close();
                }
            }
            
            if (openForward == 1) { //can't find it
                m->mothurOut("[WARNING]: can't find " + thisFileName1 + ", ignoring.\n");
            }else{  in2.close();  }
            
            int openReverse = 0;
            
            if (thisFileName2 != "") {
                ifstream in3;
                openReverse = m->openInputFile(thisFileName2, in3, "noerror");
                
                //if you can't open it, try default location
                if (openReverse == 1) {
                    if (m->getDefaultPath() != "") { //default path is set
                        string tryPath = m->getDefaultPath() + m->getSimpleName(thisFileName2);
                        m->mothurOut("Unable to open " + thisFileName2 + ". Trying default " + tryPath); m->mothurOutEndLine();
                        ifstream in3;
                        openReverse = m->openInputFile(tryPath, in3, "noerror");
                        in3.close();
                        thisFileName2 = tryPath;
                    }
                }
                
                //if you can't open it, try output location
                if (openReverse == 1) {
                    if (m->getOutputDir() != "") { //default path is set
                        string tryPath = m->getOutputDir() + m->getSimpleName(thisFileName2);
                        m->mothurOut("Unable to open " + thisFileName2 + ". Trying output directory " + tryPath); m->mothurOutEndLine();
                        ifstream in4;
                        openReverse = m->openInputFile(tryPath, in4, "noerror");
                        thisFileName2 = tryPath;
                        in4.close();
                    }
                }
                
                if (openReverse == 1) { //can't find it
                    m->mothurOut("[WARNING]: can't find " + thisFileName2 + ", ignoring pair.\n");
                }else{  in3.close();  }
            }
            
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

            
            if ((pieces.size() == 2) && (openForward != 1) && (openReverse != 1)) { //good pair and sff or fastq and oligos
                libLayout = "single";
                if (!setOligosParameter) {
                    //process pair
                    int pos = thisFileName1.find(".sff");
                    if (pos != string::npos) {//these files are sff files
                        fileOption = 1;
                        isSFF = true;
                        sfffile = thisFileName1; oligosfile = thisFileName2;
                        if (m->debug) { m->mothurOut("[DEBUG]: about to read oligos\n"); }
                        readOligos();
                        if (m->debug) { m->mothurOut("[DEBUG]: about to parse\n"); }
                        parseSffFile(files);
                        if (m->debug) { m->mothurOut("[DEBUG]: done parsing " + sfffile + "\n"); }
                    }else{
                        fileOption = 2;
                        isSFF = false;
                        fastqfile = thisFileName1; oligosfile = thisFileName2;
                        if (m->debug) { m->mothurOut("[DEBUG]: about to read oligos\n"); }
                        readOligos();
                        if (m->debug) { m->mothurOut("[DEBUG]: about to parse\n"); }
                        parseFastqFile(files);
                        if (m->debug) { m->mothurOut("[DEBUG]: done parsing " + fastqfile + "\n"); }
                    }
                }else {  runParseFastqFile = true;  libLayout = "paired"; fileOption = 3; }
            }else if((pieces.size() == 3) && (openForward != 1) && (openReverse != 1)) { //good pair and paired read
                Groups.push_back(group);
                string thisname = thisFileName1 + " " + thisFileName2;
                if (using3NONE) { thisname = thisFileName1;  }
                map<string, vector<string> >::iterator it = files.find(group);
                if (it == files.end()) {
                    vector<string> temp; temp.push_back(thisname); files[group] = temp;
                }else {
                    files[group].push_back(thisname);
                }
                fileOption = 4;
            }else if ((pieces.size() == 4) && (openForward != 1) && (openReverse != 1) && (openFindex != 1) && (openRindex != 1)) {
                libLayout = "paired"; runParseFastqFile = true; fileOption = 5;
            }
        }
        in.close();
        
        if (runParseFastqFile) {
            
            vector<string> theseFiles;
            string commandString = "fasta=f, qfile=f, file=" + file;
            
            commandString += ", oligos=" + oligosfile;
            //add in pdiffs, bdiffs, ldiffs, sdiffs, tdiffs
            if (pdiffs != 0) { commandString += ", pdiffs=" + toString(pdiffs); }
            if (bdiffs != 0) { commandString += ", bdiffs=" + toString(bdiffs); }
            if (ldiffs != 0) { commandString += ", ldiffs=" + toString(ldiffs); }
            if (sdiffs != 0) { commandString += ", sdiffs=" + toString(sdiffs); }
            if (tdiffs != 0) { commandString += ", tdiffs=" + toString(tdiffs); }
            if (m->isTrue(checkorient)) { commandString += ", checkorient=" + checkorient; }
            
            m->mothurOutEndLine();
            m->mothurOut("/******************************************/"); m->mothurOutEndLine();
            m->mothurOut("Running command: fastq.info(" + commandString + ")"); m->mothurOutEndLine();
            m->mothurCalling = true;
            
            Command* fastqinfoCommand = new ParseFastaQCommand(commandString);
            fastqinfoCommand->execute();
            
            map<string, vector<string> > filenames = fastqinfoCommand->getOutputFiles();
            map<string, vector<string> >::iterator it = filenames.find("fastq");
            if (it != filenames.end()) { theseFiles = it->second; }
            else { m->control_pressed = true; } // error in sffinfo
            
            delete fastqinfoCommand;
            m->mothurCalling = false;
            m->mothurOut("/******************************************/"); m->mothurOutEndLine();
            
            for (int i = 0; i < theseFiles.size(); i++) { outputNames.push_back(theseFiles[i]); }
            
            mapGroupToFile(files, theseFiles);
            fixMap(files);
        }
        
        inputfile = file;
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "readFile");
		exit(1);
	}
}
//**********************************************************************************************************************
int SRACommand::parseSffFile(map<string, vector<string> >& files){
	try {
        vector<string> theseFiles;
        inputfile = sfffile;
        libLayout = "single"; //controlled vocab
        
        isSFF = true;
        //run sffinfo to parse sff file into individual sampled sff files
        string commandString = "sff=" + sfffile;
        
        commandString += ", oligos=" + oligosfile;
        //add in pdiffs, bdiffs, ldiffs, sdiffs, tdiffs
        if (pdiffs != 0) { commandString += ", pdiffs=" + toString(pdiffs); }
        if (bdiffs != 0) { commandString += ", bdiffs=" + toString(bdiffs); }
        if (ldiffs != 0) { commandString += ", ldiffs=" + toString(ldiffs); }
        if (sdiffs != 0) { commandString += ", sdiffs=" + toString(sdiffs); }
        if (tdiffs != 0) { commandString += ", tdiffs=" + toString(tdiffs); }
        if (m->isTrue(checkorient)) { commandString += ", checkorient=" + checkorient; }
        
        m->mothurOutEndLine();
        m->mothurOut("/******************************************/"); m->mothurOutEndLine();
        m->mothurOut("Running command: sffinfo(" + commandString + ")"); m->mothurOutEndLine();
        m->mothurCalling = true;
        
        Command* sffinfoCommand = new SffInfoCommand(commandString);
        sffinfoCommand->execute();
        
        map<string, vector<string> > filenames = sffinfoCommand->getOutputFiles();
        map<string, vector<string> >::iterator it = filenames.find("sff");
        if (it != filenames.end()) { theseFiles = it->second; }
        else { m->control_pressed = true; } // error in sffinfo
        
        delete sffinfoCommand;
        m->mothurCalling = false;
        m->mothurOut("/******************************************/"); m->mothurOutEndLine();
        
        for (int i = 0; i < theseFiles.size(); i++) { outputNames.push_back(theseFiles[i]); }
        
        mapGroupToFile(files, theseFiles);
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "readFile");
		exit(1);
	}
}

//**********************************************************************************************************************
int SRACommand::parseFastqFile(map<string, vector<string> >& files){
	try {
        vector<string> theseFiles;
        inputfile = fastqfile;
        
        //run sffinfo to parse sff file into individual sampled sff files
        string commandString = "fastq=" + fastqfile;
        
        commandString += ", oligos=" + oligosfile;
        //add in pdiffs, bdiffs, ldiffs, sdiffs, tdiffs
        if (pdiffs != 0) { commandString += ", pdiffs=" + toString(pdiffs); }
        if (bdiffs != 0) { commandString += ", bdiffs=" + toString(bdiffs); }
        if (ldiffs != 0) { commandString += ", ldiffs=" + toString(ldiffs); }
        if (sdiffs != 0) { commandString += ", sdiffs=" + toString(sdiffs); }
        if (tdiffs != 0) { commandString += ", tdiffs=" + toString(tdiffs); }
        if (m->isTrue(checkorient)) { commandString += ", checkorient=" + checkorient; }
       
        m->mothurOutEndLine();
        m->mothurOut("/******************************************/"); m->mothurOutEndLine();
        m->mothurOut("Running command: fastq.info(" + commandString + ")"); m->mothurOutEndLine();
        m->mothurCalling = true;
        
        Command* fastqinfoCommand = new ParseFastaQCommand(commandString);
        fastqinfoCommand->execute();
        
        map<string, vector<string> > filenames = fastqinfoCommand->getOutputFiles();
        map<string, vector<string> >::iterator it = filenames.find("fastq");
        if (it != filenames.end()) { theseFiles = it->second; }
        else { m->control_pressed = true; } // error in sffinfo
        
        delete fastqinfoCommand;
        m->mothurCalling = false;
        m->mothurOut("/******************************************/"); m->mothurOutEndLine();
        
        for (int i = 0; i < theseFiles.size(); i++) { outputNames.push_back(theseFiles[i]); }
        
        mapGroupToFile(files, theseFiles);
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "readFile");
		exit(1);
	}
}
//***************************************************************************************************************
//maps group to file
int SRACommand::mapGroupToFile(map<string, vector<string> >& files, vector<string> theseFiles){
	try {
        for (int i = 0; i < Groups.size(); i++) {
            
            for (int j = 0; j < theseFiles.size(); j++) {
                
                string tempName = m->getSimpleName(theseFiles[j]);
                //if ((tempName == "GZGO5KL01.F006D146.sff") || (tempName == "G3BMWHG01.F008D021.sff") || (tempName == "GO2JXTW01.M002D125.sff") || (tempName == "GO5715J01.M003D125.sff")) { cout << Groups[i] << '\t' << theseFiles[j] << endl; }
                //cout << i << '\t' << j << '\t' << Groups[i] << '\t' << theseFiles[j] << endl;
                int pos = theseFiles[j].find(Groups[i]);
                if (pos != string::npos) { //you have a potential match, make sure you dont have a case of partial name
                    if (theseFiles[j][pos+Groups[i].length()] == '.') { //final.soil.sff vs final.soil2.sff both would match soil.
                        map<string, vector<string> >::iterator it = files.find(Groups[i]);
                        if (it == files.end()) {
                            vector<string> temp; temp.push_back(theseFiles[j]);
                            files[Groups[i]] = temp;
                        }else {
                            files[Groups[i]].push_back(theseFiles[j]);
                        }
                    }
                }
            }
        }
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "mapGroupToFile");
		exit(1);
	}
}
//***************************************************************************************************************
//fixes map to files for index files parse
int SRACommand::fixMap(map<string, vector<string> >& files){
    try {
        for (map<string, vector<string> >::iterator it = files.begin(); it != files.end(); it++) {
            
            vector<string> theseFiles = it->second;
            
            if (theseFiles.size() != 2) { m->mothurOut("[ERROR]: unexpected number of files, quitting. \n."); m->control_pressed = true; }
            
            if (m->control_pressed) { return 0; }
            
            vector<string> temp; temp.resize(1, "");
            
            for (int j = 0; j < theseFiles.size(); j++) {
                
                string tempName = m->getSimpleName(theseFiles[j]);
                int pos = theseFiles[j].find("forward.fastq");
                if (pos != string::npos) { //you have a potential match for the forward file
                    if (temp[0] == "") {
                        temp[0] = theseFiles[j];
                    }else {
                        string reverse = temp[0];
                        temp[0] = theseFiles[j] + " " + reverse;
                    }
                }else {
                    pos = theseFiles[j].find("reverse.fastq");
                    if (pos != string::npos) { //you have a potential match for the reverse file
                        if (temp[0] == "") {
                            temp[0] = theseFiles[j];
                        }else {
                            temp[0] += " " + theseFiles[j];
                        }
                    }else {
                        m->mothurOut("[ERROR]: unexpected parsing results, quitting. \n."); m->control_pressed = true;  //shouldn't get here unless the fastq.info changes the format of the output filenames???
                    }
                }
            }
            
            //cout << it->first << '\t' << temp[0] << endl;
            
            it->second = temp;
        }
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "SRACommand", "mapGroupToFile");
        exit(1);
    }
}
//***************************************************************************************************************
//checks groups and files returned from parse - removes any groups that did not get reads assigned to them, orders files.
int SRACommand::checkGroups(map<string, vector<string> >& files){
	try {
        vector<string> newGroups;
        for (int i = 0; i < Groups.size(); i++) {
            if (m->debug) { m->mothurOut("[DEBUG]: group " + toString(i) + " = " + Groups[i] + "\n"); }
            
            map<string, vector<string> >::iterator it = files.find(Groups[i]);
             //no files for this group, remove it
            if (it == files.end()) {  }
            else { newGroups.push_back(Groups[i]); }
        }
        
        Groups = newGroups;

        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "checkGroups");
		exit(1);
	}
}
//***************************************************************************************************************
int SRACommand::readOligos(){
	try {
        Oligos oligos;
        if ((fileOption == 3) || (fileOption == 5)) { oligos.read(oligosfile, false);  } //like make.contigs
        else {  oligos.read(oligosfile);  }
        
        if (m->control_pressed) { return false; } //error in reading oligos
        
        if (oligos.hasPairedPrimers() || oligos.hasPairedBarcodes())    {  pairedOligos = true;    libLayout = "paired"; }
        else                                                            {  pairedOligos = false;    libLayout = "single"; }
        
        vector<string> thisFilesLinkers = oligos.getLinkers();
        for (int i = 0; i < thisFilesLinkers.size(); i++) { linkers.push_back(thisFilesLinkers[i]); break; }
        if (thisFilesLinkers.size() > 1) { m->mothurOut("[WARNING]: the make.sra command only allows for the use of one linker at a time, disregarding all but first one.\n");  }
        
        vector<string> thisFilesSpacers = oligos.getSpacers();
        for (int i = 0; i < thisFilesSpacers.size(); i++) { spacers.push_back(thisFilesSpacers[i]); break; }
        if (thisFilesSpacers.size() > 1) { m->mothurOut("[WARNING]: the make.sra command only allows for the use of one spacer at a time, disregarding all but first one.\n");  }
        
        if (pairedOligos) {
            map<int, oligosPair> barcodes = oligos.getPairedBarcodes();
            map<int, oligosPair> primers = oligos.getPairedPrimers();
            for(map<int, oligosPair>::iterator itBar = barcodes.begin();itBar != barcodes.end();itBar++){
                for(map<int, oligosPair>::iterator itPrimer = primers.begin();itPrimer != primers.end(); itPrimer++){
                    
                    string primerName = oligos.getPrimerName(itPrimer->first);
                    string barcodeName = oligos.getBarcodeName(itBar->first);
                    
                    if ((primerName == "ignore") || (barcodeName == "ignore")) { } //do nothing
                    else if ((primerName == "") && (barcodeName == "")) { } //do nothing
                    else {
                        string comboGroupName = "";
                        string comboName = "";
                        
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
                        
                        if(((itPrimer->second).forward+(itPrimer->second).reverse) == ""){
                            if ((itBar->second).forward != "NONE") { comboName += (itBar->second).forward; }
                            if ((itBar->second).reverse != "NONE") {
                                if (comboName == "") {  comboName += (itBar->second).reverse; }
                                else {  comboName += ("."+(itBar->second).reverse);  }
                            }
                        }else{
                            if(((itBar->second).forward+(itBar->second).reverse) == ""){
                                if ((itPrimer->second).forward != "NONE") { comboName += (itPrimer->second).forward; }
                                if ((itPrimer->second).reverse != "NONE") {
                                    if (comboName == "") {  comboName += (itPrimer->second).reverse; }
                                    else {  comboName += ("."+(itPrimer->second).reverse);  }
                                }
                            }
                            else{
                                if ((itBar->second).forward != "NONE") { comboName += (itBar->second).forward; }
                                if ((itBar->second).reverse != "NONE") {
                                    if (comboName == "") {  comboName += (itBar->second).reverse; }
                                    else {  comboName += ("."+(itBar->second).reverse);  }
                                }
                                if ((itPrimer->second).forward != "NONE") {
                                    if (comboName == "") {  comboName += (itPrimer->second).forward; }
                                    else {  comboName += ("."+(itPrimer->second).forward);  }
                                }
                                if ((itPrimer->second).reverse != "NONE") {
                                    if (comboName == "") {  comboName += (itPrimer->second).reverse; }
                                    else {  comboName += ("."+(itPrimer->second).reverse);  }
                                }
                            }
                        }
                        
                        if (comboName != "") {  comboGroupName +=  "_" + comboName;  }
                        uniqueNames.insert(comboGroupName);
                        
                        map<string, string>::iterator itGroup2Barcode = Group2Barcode.find(comboGroupName);
                        if (itGroup2Barcode == Group2Barcode.end()) {
                            string temp = (itBar->second).forward+"."+(itBar->second).reverse;
                            Group2Barcode[comboGroupName] = temp;
                        }else {
                            string temp = (itBar->second).forward+"."+(itBar->second).reverse;
                            if ((temp != ".") && (temp != itGroup2Barcode->second)) {
                                m->mothurOut("[ERROR]: group and barcodes/primers not unique. Should never get here.\n");
                            }
                        }
                        
                        itGroup2Barcode = Group2Primer.find(comboGroupName);
                        if (itGroup2Barcode == Group2Primer.end()) {
                            string temp = ((itPrimer->second).forward+"."+(itPrimer->second).reverse);
                            Group2Primer[comboGroupName] = temp;
                        }else {
                            string temp = ((itPrimer->second).forward+"."+(itPrimer->second).reverse);
                            if ((temp != ".") && (temp != itGroup2Barcode->second)) {
                                m->mothurOut("[ERROR]: group and barcodes/primers not unique. Should never get here.\n");
                            }
                        }
                    }
                }
            }
        }else {
            map<string, int> barcodes = oligos.getBarcodes() ;
            map<string, int> primers = oligos.getPrimers();
            for(map<string, int>::iterator itBar = barcodes.begin();itBar != barcodes.end();itBar++){
                for(map<string, int>::iterator itPrimer = primers.begin();itPrimer != primers.end(); itPrimer++){
                    
                    string primerName = oligos.getPrimerName(itPrimer->second);
                    string barcodeName = oligos.getBarcodeName(itBar->second);
                    
                    if ((primerName == "ignore") || (barcodeName == "ignore")) { } //do nothing
                    else if ((primerName == "") && (barcodeName == "")) { } //do nothing
                    else {
                        string comboGroupName = "";
                        string comboName = "";
                        
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
                        
                        if(itPrimer->first == ""){
                            comboName = itBar->first;
                        }else{
                            if(itBar->first == ""){
                                comboName = itPrimer->first;
                            }
                            else{
                                comboName = itBar->first + "." + itPrimer->first;
                            }
                        }
                        //cout << comboGroupName << '\t' << comboName << endl;
                        if (comboName != "") {  comboGroupName +=  "_" + comboName;  }
                        uniqueNames.insert(comboGroupName);
                        
                        map<string, string >::iterator itGroup2Barcode = Group2Barcode.find(comboGroupName);
                        if (itGroup2Barcode == Group2Barcode.end()) {
                            string temp = (itBar->first);
                            Group2Barcode[comboGroupName] = temp;
                        }else {
                            string temp = (itBar->first);
                           if ((temp != ".") && (temp != itGroup2Barcode->second)) {
                                m->mothurOut("[ERROR]: group and barcodes/primers not unique. Should never get here.\n");
                            }
                        }
                        
                        itGroup2Barcode = Group2Primer.find(comboGroupName);
                        if (itGroup2Barcode == Group2Primer.end()) {
                            string temp = (itPrimer->first);
                            Group2Primer[comboGroupName] = temp;
                        }else {
                            string temp = (itPrimer->first);
                            if ((temp != ".") && (temp != itGroup2Barcode->second)) {
                                m->mothurOut("[ERROR]: group and barcodes/primers not unique. Should never get here.\n");
                            }
                        }
                    }
                }
            }
        }
        
        if (m->debug) { int count = 0; for (set<string>::iterator it = uniqueNames.begin(); it != uniqueNames.end(); it++) { m->mothurOut("[DEBUG]: " + toString(count) + " groupName = " + *it + "\n"); count++; } }
        
        Groups.clear();
        for (set<string>::iterator it = uniqueNames.begin(); it != uniqueNames.end(); it++) {  Groups.push_back(*it);  }
        
        
		return true;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "readOligos");
		exit(1);
	}
}
//********************************************************************/
//_LS454-ILLUMINA-ION_TORRENT-PACBIO_SMRT
bool SRACommand::checkCasesPlatforms(string& platform){
	try {
        string original = platform;
        bool isOkay = true;
        
        //remove users possible case errors
        for (int i = 0; i < platform.size(); i++) { platform[i] = toupper(platform[i]); }
        
        //_LS454-ILLUMINA-ION_TORRENT-PACBIO_SMRT
        
            if ((platform == "_LS454") || (platform == "ILLUMINA") || (platform == "ION_TORRENT") || (platform == "PACBIO_SMRT") || (platform == "454")) { }
            else { isOkay = false; }
        
            if (isOkay) {
                if (platform == "454")   {  platform = "_LS454"; }
            }else {
                m->mothurOut("[ERROR]: " + original + " is not a valid platform option.  Valid platform options are _LS454, ILLUMINA-ION, TORRENT or PACBIO_SMRT."); m->mothurOutEndLine(); abort = true;
            }
            
            return isOkay;
    }
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "checkCasesPlatforms");
		exit(1);
	}
}
//********************************************************************/
//454_GS-454_GS_20-454_GS_FLX-454_GS_FLX_Titanium-454_GS_Junior-Illumina_Genome_Analyzer-Illumina_Genome_Analyzer_II-Illumina_Genome_Analyzer_IIx-Illumina_HiSeq_2000-Illumina_HiSeq_1000-Illumina_MiSeq-PacBio_RS-Ion_Torrent_PGM-unspecified
bool SRACommand::checkCasesInstrumentModels(string& instrumentModel){
	try {
        string original = instrumentModel;
        bool isOkay = true;
        
        //remove users possible case errors
        for (int i = 0; i < instrumentModel.size(); i++) { instrumentModel[i] = toupper(instrumentModel[i]); }
        
        //_LS454-ILLUMINA-ION_TORRENT-PACBIO_SMRT
        if (platform == "_LS454") { //instrument model options are 454_GS-454_GS_20-454_GS_FLX-454_GS_FLX_Titanium-454_GS_Junior-unspecified
            if ((instrumentModel == "454_GS") || (instrumentModel == "454_GS_20") || (instrumentModel == "454_GS_FLX") || (instrumentModel == "454_GS_FLX_TITANIUM") || (instrumentModel == "454_GS_JUNIOR") || (instrumentModel == "UNSPECIFIED")) { }
            else { isOkay = false; }
            if (isOkay) {
                if (instrumentModel == "454_GS_FLX_TITANIUM")   {  instrumentModel = "454_GS_FLX_Titanium"; }
                if (instrumentModel == "454_GS_JUNIOR")         {  instrumentModel = "454_GS_Junior";       }
                if (instrumentModel == "UNSPECIFIED")           {  instrumentModel = "unspecified";         }
            }else {
                m->mothurOut("[ERROR]: " + original + " is not a valid instrument option for the " + platform + " platform.  Valid instrument options are 454_GS, 454_GS_20, 454_GS_FLX, 454_GS_FLX_Titanium, 454_GS_Junior or unspecified."); m->mothurOutEndLine(); abort = true;
            }
            
        }else if (platform == "ILLUMINA") { //instrument model options are Illumina_Genome_Analyzer-Illumina_Genome_Analyzer_II-Illumina_Genome_Analyzer_IIx-Illumina_HiSeq_2000-Illumina_HiSeq_1000-Illumina_MiSeq-unspecified
            if ((instrumentModel == "ILLUMINA_GENOME_ANALYZER") || (instrumentModel == "ILLUMINA_GENOME_ANALYZER_II") || (instrumentModel == "ILLUMINA_GENOME_ANALYZER_IIX") || (instrumentModel == "ILLUMINA_HISEQ_2000") || (instrumentModel == "ILLUMINA_HISEQ_1000") || (instrumentModel == "ILLUMINA_MISEQ") || (instrumentModel == "UNSPECIFIED")) { }
            else { isOkay = false; }
            
            if (isOkay) {
                if (instrumentModel == "ILLUMINA_GENOME_ANALYZER")          {  instrumentModel = "Illumina_Genome_Analyzer";        }
                if (instrumentModel == "ILLUMINA_GENOME_ANALYZER_II")       {  instrumentModel = "Illumina_Genome_Analyzer_II";     }
                if (instrumentModel == "ILLUMINA_GENOME_ANALYZER_IIX")      {  instrumentModel = "Illumina_Genome_Analyzer_IIx";    }
                if (instrumentModel == "ILLUMINA_HISEQ_2000")               {  instrumentModel = "Illumina_HiSeq_2000";             }
                if (instrumentModel == "ILLUMINA_HISEQ_1000")               {  instrumentModel = "Illumina_HiSeq_1000";             }
                if (instrumentModel == "ILLUMINA_MISEQ")                    {  instrumentModel = "Illumina_MiSeq";                  }
                if (instrumentModel == "UNSPECIFIED")                       {  instrumentModel = "unspecified";                     }
            }else {
                m->mothurOut("[ERROR]: " + original + " is not a valid instrument option for the " + platform + " platform.  Valid instrument options are Illumina_Genome_Analyzer, Illumina_Genome_Analyzer_II, Illumina_Genome_Analyzer_IIx, Illumina_HiSeq_2000, Illumina_HiSeq_1000, Illumina_MiSeq or unspecified."); m->mothurOutEndLine(); abort = true;
            }
            
        }else if (platform == "ION_TORRENT") { //instrument model options are Ion_Torrent_PGM-unspecified
            if ((instrumentModel == "ION_TORRENT_PGM")  || (instrumentModel == "UNSPECIFIED")) { }
            else { isOkay = false; }
            
            if (isOkay) {
                if (instrumentModel == "ION_TORRENT_PGM")          {  instrumentModel = "Ion_Torrent_PGM";        }
                if (instrumentModel == "UNSPECIFIED")              {  instrumentModel = "unspecified";            }
            }else {
                m->mothurOut("[ERROR]: " + original + " is not a valid instrument option for the " + platform + " platform.  Valid instrument options are Ion_Torrent_PGM or unspecified."); m->mothurOutEndLine(); abort = true;
            }
        }else if (platform == "PACBIO_SMRT") { //instrument model options are PacBio_RS-unspecified
            if ((instrumentModel == "PACBIO_RS")  || (instrumentModel == "UNSPECIFIED")) { }
            else { isOkay = false; }
            
            if (isOkay) {
                if (instrumentModel == "PACBIO_RS")          {  instrumentModel = "PacBio_RS";        }
                if (instrumentModel == "UNSPECIFIED")        {  instrumentModel = "unspecified";      }
            }else {
                m->mothurOut("[ERROR]: " + original + " is not a valid instrument option for the " + platform + " platform.  Valid instrument options are PacBio_RS or unspecified."); m->mothurOutEndLine(); abort = true;
            }
        }
        return isOkay;
    }
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "checkCasesInstrumentModels");
		exit(1);
	}
}
//**********************************************************************************************************************
//AMPLICON,WGA,WGS,WGX,RNA-Seq,miRNA-Seq,WCS,CLONE,POOLCLONE,CLONEEND,FINISHING,ChIP-Seq,MNase-Seq,DNase-Hypersensitivity,Bisulfite-Seq,Tn-Seq,EST,FL-cDNA,CTS,MRE-Seq,MeDIP-Seq,MBD-Seq,OTHER
bool SRACommand::checkCasesLibStrategy(string& libStrategy){
	try {
        string original = libStrategy;
        bool isOkay = true;
        
        //remove users possible case errors
        for (int i = 0; i < libStrategy.size(); i++) { libStrategy[i] = toupper(libStrategy[i]); }
        
        if ((libStrategy == "AMPLICON") || (libStrategy == "WGA") || (libStrategy == "WGS") || (libStrategy == "WGX") || (libStrategy == "RNA-SEQ") || (libStrategy == "MIRNA-SEQ") || (libStrategy == "WCS") || (libStrategy == "CLONE") || (libStrategy == "POOLCLONE") || (libStrategy == "CLONEEND") || (libStrategy == "FINISHING") || (libStrategy == "CHIP-SEQ") || (libStrategy == "MNASE-SEQ") || (libStrategy == "DNASE-HYPERSENSITIVITY") || (libStrategy == "BISULFITE-SEQ") || (libStrategy == "TN-SEQ") || (libStrategy == "EST") || (libStrategy == "FL-CDNA") || (libStrategy == "CTS") || (libStrategy == "MRE-SEQ")|| (libStrategy == "MEDIP-SEQ") || (libStrategy == "MBD-SEQ") || (libStrategy == "OTHER")) { }
        else { isOkay = false; }
        
        if (isOkay) {
            if (libStrategy == "RNA-SEQ")                   {  libStrategy = "RNA-Seq";                 }
            if (libStrategy == "MIRNA-SEQ")                 {  libStrategy = "miRNA-Seq";               }
            if (libStrategy == "CHIP-SEQ")                  {  libStrategy = "ChIP-Seq";                }
            if (libStrategy == "MNASE-SEQ")                 {  libStrategy = "MNase-Seq";               }
            if (libStrategy == "DNASE-HYPERSENSITIVITY")    {  libStrategy = "DNase-Hypersensitivity";  }
            if (libStrategy == "BISULFITE-SEQ")             {  libStrategy = "Bisulfite-Seq";           }
            if (libStrategy == "TN-SEQ")                    {  libStrategy = "Tn-Seq";                  }
            if (libStrategy == "FL-CDNA")                   {  libStrategy = "FL-cDNA";                 }
            if (libStrategy == "MRE-SEQ")                   {  libStrategy = "MRE-Seq";                 }
            if (libStrategy == "MEDIP-SEQ")                 {  libStrategy = "MeDIP-Seq";               }
            }else {
            m->mothurOut("[ERROR]: " + original + " is not a valid libstrategy option.  Valid libstrategy options are AMPLICON,WGA,WGS,WGX,RNA-Seq,miRNA-Seq,WCS,CLONE,POOLCLONE,CLONEEND,FINISHING,ChIP-Seq,MNase-Seq,DNase-Hypersensitivity,Bisulfite-Seq,Tn-Seq,EST,FL-cDNA,CTS,MRE-Seq,MeDIP-Seq,MBD-Seq or OTHER."); m->mothurOutEndLine(); abort = true;
        }
        
        return isOkay;
    }
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "checkCasesLibStrategy");
		exit(1);
	}
}

//**********************************************************************************************************************
//METAGENOMIC,GENOMIC,TRANSCRIPTOMIC,METATRANSCRIPTOMIC,SYNTHETIC,VIRAL_RNA,OTHER
bool SRACommand::checkCasesLibSource(string& libSource){
	try {
        string original = libSource;
        bool isOkay = true;
        
        //remove users possible case errors
        for (int i = 0; i < libSource.size(); i++) { libSource[i] = toupper(libSource[i]); }
        
        if ((libSource == "METAGENOMIC") || (libSource == "GENOMIC") || (libSource == "TRANSCRIPTOMIC") || (libSource == "METATRANSCRIPTOMIC") || (libSource == "SYNTHETIC") || (libSource == "VIRAL_RNA") || (libSource == "OTHER")) { }
        else { isOkay = false; }
        
        if (isOkay) {
            
        }else {
            m->mothurOut("[ERROR]: " + original + " is not a valid libsource option.  Valid libsource options are METAGENOMIC,GENOMIC,TRANSCRIPTOMIC,METATRANSCRIPTOMIC,SYNTHETIC,VIRAL_RNA or OTHER."); m->mothurOutEndLine(); abort = true;
        }
        
        return isOkay;
    }
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "checkCasesLibStrategy");
		exit(1);
	}
}

//**********************************************************************************************************************
//PCR,RANDOM,RANDOM_PCR,RT-PCR,HMPR,MF,CF-S,CF-H,CF-T,CF-M,MDA,MSLL,cDNA,ChIP,MNase,DNAse,Hybrid_Selection,Reduced_Representation,Restriction_Digest,5-methylcytidine_antibody,MBD2_protein_methyl-CpG_binding_domain,CAGE,RACE,size_fractionation,Padlock_probes_capture_method,other,unspecified
bool SRACommand::checkCasesLibSelection(string& libSelection){
	try {
        string original = libSelection;
        bool isOkay = true;
        
        //remove users possible case errors
        for (int i = 0; i < libSelection.size(); i++) { libSelection[i] = toupper(libSelection[i]); }
        
        if ((libSelection == "PCR") || (libSelection == "RANDOM") || (libSelection == "RANDOM_PCR") || (libSelection == "RT-PCR") || (libSelection == "HMPR") || (libSelection == "MF") || (libSelection == "CF-S") || (libSelection == "CF-H") || (libSelection == "CF-T") || (libSelection == "CF-M") || (libSelection == "MDA") || (libSelection == "MSLL") || (libSelection == "CDNA") || (libSelection == "CHIP") || (libSelection == "MNASE") || (libSelection == "DNASE") || (libSelection == "HYBRID_SELECTION") || (libSelection == "REDUCED_REPRESENTATION") || (libSelection == "RESTRICTION_DIGEST") || (libSelection == "5-METHYLCYTIDINE_ANTIBODY") || (libSelection == "MBD2_PROTEIN_METHYL-CPG_BINDING_DOMAIN") || (libSelection == "CAGE") || (libSelection == "RACE") || (libSelection == "SIZE_FRACTIONATION") || (libSelection == "PADLOCK_PROBES_CAPTURE_METHOD") || (libSelection == "OTHER") || (libSelection == "UNSPECIFIED")) { }
        else { isOkay = false; }
        
        if (isOkay) {
            if (libSelection == "CDNA")                                         {  libSelection = "cDNA";                                       }
            if (libSelection == "CHIP")                                         {  libSelection = "ChIP";                                       }
            if (libSelection == "MNASE")                                        {  libSelection = "MNase";                                      }
            if (libSelection == "DNASE")                                        {  libSelection = "DNAse";                                      }
            if (libSelection == "HYBRID_SELECTION")                             {  libSelection = "Hybrid_Selection";                           }
            if (libSelection == "REDUCED_REPRESENTATION")                       {  libSelection = "Reduced_Representation";                     }
            if (libSelection == "RESTRICTION_DIGEST")                           {  libSelection = "Restriction_Digest";                         }
            if (libSelection == "5-METHYLCYTIDINE_ANTIBODY")                    {  libSelection = "5-methylcytidine_antibody";                  }
            if (libSelection == "MBD2_PROTEIN_METHYL-CPG_BINDING_DOMAIN")       {  libSelection = "MBD2_protein_methyl-CpG_binding_domain";     }
            if (libSelection == "SIZE_FRACTIONATION")                           {  libSelection = "size_fractionation";                         }
            if (libSelection == "PADLOCK_PROBES_CAPTURE_METHOD")                {  libSelection = "Padlock_probes_capture_method";              }
            if (libSelection == "OTHER")                                        {  libSelection = "other";                                      }
            if (libSelection == "UNSPECIFIED")                                  {  libSelection = "unspecified";                                }
            
        }else {
            m->mothurOut("[ERROR]: " + original + " is not a valid libselection option.  Valid libselection options are PCR,RANDOM,RANDOM_PCR,RT-PCR,HMPR,MF,CF-S,CF-H,CF-T,CF-M,MDA,MSLL,cDNA,ChIP,MNase,DNAse,Hybrid_Selection,Reduced_Representation,Restriction_Digest,5-methylcytidine_antibody,MBD2_protein_methyl-CpG_binding_domain,CAGE,RACE,size_fractionation,Padlock_probes_capture_method,other or unspecified."); m->mothurOutEndLine(); abort = true;
        }
        
        return isOkay;
    }
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "checkCasesLibSelection");
		exit(1);
	}
}
//**********************************************************************************************************************
//METAGENOME,GENOME_SEQUENCING,METAGENOMIC_ASSEMBLY,ASSEMBLY,TRANSCRIPTOME,PROTEOMIC,MAP,CLONE_ENDS,TARGETED_LOCI,RANDOM_SURVEY,EXOME,VARIATION,EPIGENOMICS,PHENOTYPE,GENOTYPE,OTHER
bool SRACommand::checkCasesDataType(string& dataType){
	try {
        string original = dataType;
        bool isOkay = true;
        
        //remove users possible case errors
        for (int i = 0; i < dataType.size(); i++) { dataType[i] = toupper(dataType[i]); }
        
        if ((dataType == "METAGENOME") || (dataType == "GENOME_SEQUENCING") || (dataType == "METAGENOMIC_ASSEMBLY") || (dataType == "ASSEMBLY") || (dataType == "TRANSCRIPTOME") || (dataType == "PROTEOMIC") || (dataType == "MAP") || (dataType == "CLONE_ENDS") || (dataType == "TARGETED_LOCI") || (dataType == "RANDOM_SURVEY") || (dataType == "EXOME") || (dataType == "VARIATION") || (dataType == "EPIGENOMICS") || (dataType == "PHENOTYPE") || (dataType == "GENOTYPE") || (dataType == "OTHER"))
        { dataType = original; }
        else { isOkay = false; }
        
        if (isOkay) {
            
        }else {
            m->mothurOut("[ERROR]: " + original + " is not a valid datatype option.  Valid datatype options are METAGENOME,GENOME_SEQUENCING,METAGENOMIC_ASSEMBLY,ASSEMBLY,TRANSCRIPTOME,PROTEOMIC,MAP,CLONE_ENDS,TARGETED_LOCI,RANDOM_SURVEY,EXOME,VARIATION,EPIGENOMICS,PHENOTYPE,GENOTYPE,OTHER."); m->mothurOutEndLine(); abort = true;
        }
        
        return isOkay;
    }
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "checkCasesDataType");
		exit(1);
	}
}
//**********************************************************************************************************************
bool SRACommand::sanityCheckMiMarksGroups(){
	try {
        bool isOkay = true;
        
        for (int i = 0; i < Groups.size(); i++) {
            if (m->control_pressed) { break; }
            
            map<string, map<string, string> >::iterator it = mimarks.find(Groups[i]);
            if (it == mimarks.end()) {
                isOkay = false;
                m->mothurOut("[ERROR]: MIMarks file is missing group " + Groups[i] + ", please correct.\n");
            }
        }
        
        if (!isOkay) { m->control_pressed = true; }
        
        return isOkay;
    }
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "sanityCheckMiMarksGroups");
		exit(1);
	}
}
//**********************************************************************************************************************
/*
 file option 1
 
 sfffile1   oligosfile1
 sfffile2   oligosfile2
 ...
 
 file option 2
 
 fastqfile1 oligosfile1
 fastqfile2 oligosfile2
 ...
 
 file option 3
 
 ffastqfile1 rfastqfile1
 ffastqfile2 rfastqfile2
 ...
 
 file option 4
 
 group fastqfile  fastqfile
 group fastqfile  fastqfile
 group fastqfile  fastqfile
 ...
 
 file option 5
 
 My.forward.fastq My.reverse.fastq none My.rindex.fastq //none is an option is no forward or reverse index file
 ...
*/

int SRACommand::findFileOption(){
    try {
        ifstream in;
        m->openInputFile(file, in);
        
        fileOption = 0;
        
        while(!in.eof()) {
            
            if (m->control_pressed) { return 0; }
            
            string line = m->getline(in);  m->gobble(in);
            vector<string> pieces = m->splitWhiteSpace(line);
            
            if (pieces.size() == 2) { //good pair and sff or fastq and oligos
                if (!setOligosParameter) {
                    fileOption = 12; //1 or 2
                }else {  fileOption = 3; }
            }else if(pieces.size() == 3) { //good pair and paired read
                fileOption = 4;
            }else if (pieces.size() == 4) {
                fileOption = 5;
            }
            break;
        }
        in.close();
        
        return fileOption;
    }
    catch(exception& e) {
        m->errorOut(e, "SRACommand", "findFileOption");
        exit(1);
    }
}

//**********************************************************************************************************************
