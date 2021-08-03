//
//  prcseqscommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/14/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "pcrseqscommand.h"

//**********************************************************************************************************************
vector<string> PcrSeqsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta",false,true,true); parameters.push_back(pfasta);
		CommandParameter poligos("oligos", "InputTypes", "", "", "ecolioligos", "none", "none","",false,false,true); parameters.push_back(poligos);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","name",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","count",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","group",false,false,true); parameters.push_back(pgroup);
        CommandParameter ptax("taxonomy", "InputTypes", "", "", "none", "none", "none","taxonomy",false,false,true); parameters.push_back(ptax);
        CommandParameter preorient("checkorient", "Boolean", "", "T", "", "", "","",false,false,true); parameters.push_back(preorient);
        CommandParameter pecoli("ecoli", "InputTypes", "", "", "ecolioligos", "none", "none","",false,false); parameters.push_back(pecoli);
		CommandParameter pstart("start", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pstart);
		CommandParameter pend("end", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pend);
 		CommandParameter pnomatch("nomatch", "Multiple", "reject-keep", "reject", "", "", "","",false,false); parameters.push_back(pnomatch);
        CommandParameter ppdiffs("pdiffs", "Number", "", "0", "", "", "","",false,false,true); parameters.push_back(ppdiffs);
        CommandParameter prdiffs("rdiffs", "Number", "", "0", "", "", "","",false,false,true); parameters.push_back(prdiffs);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pkeepprimer("keepprimer", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pkeepprimer);
        CommandParameter pkeepdots("keepdots", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pkeepdots);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["fasta"] = tempOutNames;
        outputTypes["taxonomy"] = tempOutNames;
        outputTypes["group"] = tempOutNames;
        outputTypes["name"] = tempOutNames;
        outputTypes["accnos"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
        
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string PcrSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The pcr.seqs command reads a fasta file.\n";
        helpString += "The pcr.seqs command parameters are fasta, oligos, name, group, count, taxonomy, ecoli, start, end, nomatch, pdiffs, rdiffs, processors, keepprimer and keepdots.\n";
		helpString += "The ecoli parameter is used to provide a fasta file containing a single reference sequence (e.g. for e. coli) this must be aligned. Mothur will trim to the start and end positions of the reference sequence.\n";
        helpString += "The start parameter allows you to provide a starting position to trim to.\n";
        helpString += "The end parameter allows you to provide a ending position to trim from.\n";
        helpString += "The nomatch parameter allows you to decide what to do with sequences where the primer is not found. Default=reject, meaning remove from fasta file.  if nomatch=true, then do nothing to sequence.\n";
        helpString += "The checkorient parameter will look for the reverse compliment of the barcode or primer in the sequence. If found the sequence is flipped. The default is true.\n";
        helpString += "The processors parameter allows you to use multiple processors.\n";
        helpString += "The keepprimer parameter allows you to keep the primer, default=false.\n";
        helpString += "The keepdots parameter allows you to keep the leading and trailing .'s, default=true.\n";
        helpString += "The pdiffs parameter is used to specify the number of differences allowed in the forward primer. The default is 0.\n";
        helpString += "The rdiffs parameter is used to specify the number of differences allowed in the reverse primer. The default is 0.\n";
		;
		helpString += "For more details please check out the wiki http://www.mothur.org/wiki/Pcr.seqs .\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string PcrSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta")            {   pattern = "[filename],pcr,[extension]-[filename],[tag],pcr,[extension]";    }
        else if (type == "taxonomy")    {   pattern = "[filename],pcr,[extension]";    }
        else if (type == "name")        {   pattern = "[filename],pcr,[extension]";    }
        else if (type == "group")       {   pattern = "[filename],pcr,[extension]";    }
        else if (type == "count")       {   pattern = "[filename],pcr,[extension]";    }
        else if (type == "accnos")      {   pattern = "[filename],bad.accnos";    }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "PcrSeqsCommand", "getOutputPattern");
        exit(1);
    }
}
//***************************************************************************************************************
PcrSeqsCommand::PcrSeqsCommand(string option) : Command()  {
	try {
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not found") { 				
				fastafile = current->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter.\n");  }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required.\n"); abort = true; }
			}else if (fastafile == "not open") { fastafile = ""; abort = true; }	
			else { current->setFastaFile(fastafile); }
			
             
					if (outputdir == ""){    outputdir = util.hasPath(fastafile);	}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.valid(parameters, "keepprimer");  if (temp == "not found")    {	temp = "f";	}
			keepprimer = util.isTrue(temp);	
            
            temp = validParameter.valid(parameters, "keepdots");  if (temp == "not found")    {	temp = "t";	}
			keepdots = util.isTrue(temp);	
            
			temp = validParameter.validFile(parameters, "oligos");
			if (temp == "not found"){	oligosfile = "";		}
			else if(temp == "not open"){	oligosfile = ""; abort = true;	} 
			else					{	oligosfile = temp; current->setOligosFile(oligosfile);		}
			
            ecolifile = validParameter.validFile(parameters, "ecoli");
			if (ecolifile == "not found"){	ecolifile = "";		}
			else if(ecolifile == "not open"){	ecolifile = ""; abort = true;	} 
			
            namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not found"){	namefile = "";		}
			else if(namefile == "not open"){	namefile = ""; abort = true;	} 
            else { current->setNameFile(namefile); }
            
            groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not found"){	groupfile = "";		}
			else if(groupfile == "not open"){	groupfile = ""; abort = true;	} 
            else { current->setGroupFile(groupfile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { countfile = ""; abort = true; }
			else if (countfile == "not found") { countfile = "";  }	
			else { current->setCountFile(countfile); }
            
            if ((namefile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: name or count.\n"); abort = true;
            }
			
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count.\n"); abort = true;
            }
            
            taxfile = validParameter.validFile(parameters, "taxonomy");
			if (taxfile == "not found"){	taxfile = "";		}
			else if(taxfile == "not open"){	taxfile = ""; abort = true;	} 
            else { current->setTaxonomyFile(taxfile); }
			 			
			temp = validParameter.valid(parameters, "start");	if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, start);
            
            temp = validParameter.valid(parameters, "end");	if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, end);
			
			temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
                    
            temp = validParameter.valid(parameters, "pdiffs");		if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, pdiffs);
            
            temp = validParameter.valid(parameters, "rdiffs");		if (temp == "not found") { temp = "0"; }
            util.mothurConvert(temp, rdiffs);
			
            temp = validParameter.valid(parameters, "checkorient");        if (temp == "not found") { temp = "T"; }
            reorient = util.isTrue(temp);
            
            nomatch = validParameter.valid(parameters, "nomatch");	if (nomatch == "not found") { nomatch = "reject"; }
			
            if ((nomatch != "reject") && (nomatch != "keep")) { m->mothurOut("[ERROR]: " + nomatch + " is not a valid entry for nomatch. Choices are reject and keep.\n");  abort = true; }
            
            //didnt set anything
			if ((oligosfile == "") && (ecolifile == "") && (start == -1) && (end == -1)) {
                m->mothurOut("[ERROR]: You did not set any options. Please provide an oligos or ecoli file, or set start or end.\n"); abort = true;
            }
            
            if ((oligosfile == "") && (ecolifile == "") && (start < 0) && (end == -1)) { m->mothurOut("[ERROR]: Invalid start value.\n"); abort = true; }
            
            if ((ecolifile != "") && (start != -1) && (end != -1)) {
                m->mothurOut("[ERROR]: You provided an ecoli file , but set the start or end parameters. Unsure what you intend.  When you provide the ecoli file, mothur thinks you want to use the start and end of the sequence in the ecoli file.\n"); abort = true;
            }

            
            if ((oligosfile != "") && (ecolifile != "")) {
                 m->mothurOut("[ERROR]: You can not use an ecoli file at the same time as an oligos file.\n"); abort = true;
            }
            
			//check to make sure you didn't forget the name file by mistake			
			if (countfile == "") { 
                if (namefile == "") {
                    vector<string> files; files.push_back(fastafile);
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
                }
            }
		}
        
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "PcrSeqsCommand");
		exit(1);
	}
}
//***************************************************************************************************************

int PcrSeqsCommand::execute(){
	try{
        
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
        long start = time(NULL);
        fileAligned = true; pairedOligos = false;
        
        string thisOutputDir = outputdir;
		if (outputdir == "") {  thisOutputDir += util.hasPath(fastafile);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastafile));
        variables["[extension]"] = util.getExtension(fastafile);
		string trimSeqFile = getOutputFileName("fasta",variables);
		outputNames.push_back(trimSeqFile); outputTypes["fasta"].push_back(trimSeqFile);
        variables["[tag]"] = "scrap";
        string badSeqFile = getOutputFileName("fasta",variables);
        length = 0;
        
        if (m->getControl_pressed()) {  return 0; }

        set<string> badNames;
        long long numFastaSeqs = createProcesses(fastafile, trimSeqFile, badSeqFile, badNames);
		
		if (m->getControl_pressed()) {  return 0; }
        
        thisOutputDir = outputdir;
        if (outputdir == "") {  thisOutputDir += util.hasPath(fastafile);  }
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastafile));
        string outputFileName = getOutputFileName("accnos",variables);

        //don't write or keep if blank
        bool wroteAccnos = false;
        if (badNames.size() != 0)   { writeAccnos(badNames, outputFileName);    wroteAccnos = true;   outputNames.push_back(outputFileName); outputTypes["accnos"].push_back(outputFileName);  }
        else { m->mothurOut("[NOTE]: no sequences were bad, removing " + outputFileName + "\n\n"); }
        
        if (util.isBlank(badSeqFile)) { util.mothurRemove(badSeqFile);  }
        else { outputNames.push_back(badSeqFile); outputTypes["fasta"].push_back(badSeqFile); }
        
        if (wroteAccnos) {
            string inputStringTemp = "";
            if (countfile != "")            {   inputStringTemp += ", count=" + countfile;  }
            else{
                if (namefile != "")         {   inputStringTemp += ", name=" + namefile;    }
                if (groupfile != "")        {   inputStringTemp += ", group=" + groupfile;  }
            }
            if (taxfile != "")              {   inputStringTemp += ", taxonomy=" + taxfile;  }
            string inputString = "accnos=" + outputFileName + inputStringTemp;
            
            if (inputStringTemp != "") {
                m->mothurOut("/******************************************/\n");
                m->mothurOut("Running command: remove.seqs(" + inputString + ")\n");
                current->setMothurCalling(true);
                
                Command* removeCommand = new RemoveSeqsCommand(inputString);
                removeCommand->execute();
                
                map<string, vector<string> > filenames = removeCommand->getOutputFiles();
                
                delete removeCommand;
                current->setMothurCalling(false);
                
                if (groupfile != "") {
                    thisOutputDir = outputdir;
                    if (outputdir == "") {  thisOutputDir += util.hasPath(groupfile);  }
                    variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(groupfile));
                    variables["[extension]"] = util.getExtension(groupfile);
                    string outGroup = getOutputFileName("group", variables);
                    util.renameFile(filenames["group"][0], outGroup);
                    outputNames.push_back(outGroup); outputTypes["group"].push_back(outGroup);
                }
                
                if (namefile != "") {
                    thisOutputDir = outputdir;
                    if (outputdir == "") {  thisOutputDir += util.hasPath(namefile);  }
                    variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(namefile));
                    variables["[extension]"] = util.getExtension(namefile);
                    string outName = getOutputFileName("name", variables);
                    util.renameFile(filenames["name"][0], outName);
                    outputNames.push_back(outName); outputTypes["name"].push_back(outName);
                }
                
                if (countfile != "") {
                    thisOutputDir = outputdir;
                    if (outputdir == "") {  thisOutputDir += util.hasPath(countfile);  }
                    variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(countfile));
                    variables["[extension]"] = util.getExtension(countfile);
                    string outCount = getOutputFileName("count", variables);
                    util.renameFile(filenames["count"][0], outCount);
                    outputNames.push_back(outCount); outputTypes["count"].push_back(outCount);
                }
                
                if (taxfile != "")              {
                    thisOutputDir = outputdir;
                    if (outputdir == "") {  thisOutputDir += util.hasPath(taxfile);  }
                    variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(taxfile));
                    variables["[extension]"] = util.getExtension(taxfile);
                    string outputFileName = getOutputFileName("taxonomy", variables);
                    util.renameFile(filenames["taxonomy"][0], outputFileName);
                    outputNames.push_back(outputFileName); outputTypes["taxonomy"].push_back(outputFileName);
                }
                m->mothurOut("/******************************************/\n"); 
            }
        }
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
        
        m->mothurOut("It took " + toString(time(NULL) - start) + " secs to screen " + toString(numFastaSeqs) + " sequences.\n");
        
		m->mothurOut("\nOutput File Names: \n");
		for (int i = 0; i < outputNames.size(); i++) { m->mothurOut(outputNames[i]); m->mothurOutEndLine(); }
		m->mothurOutEndLine();
		m->mothurOutEndLine();
		
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
		
		itTypes = outputTypes.find("group");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setGroupFile(currentName); }
		}
		
		itTypes = outputTypes.find("accnos");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setAccnosFile(currentName); }
		}
		
		itTypes = outputTypes.find("taxonomy");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setTaxonomyFile(currentName); }
		}
        
        itTypes = outputTypes.find("count");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
		}
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "execute");
		exit(1);
	}
}
//***************************************************************************************************************
bool readOligos(Oligos& oligos, string oligosfile, bool& pairedOligos, int& numFPrimers, int& numRPrimers, MothurOut* m){
    try {
        oligos.read(oligosfile);
        
        if (m->getControl_pressed()) { return false; } //error in reading oligos
        
        if (oligos.hasPairedPrimers()) {
            pairedOligos = true;
            numFPrimers = oligos.getPairedPrimers().size();
            numRPrimers = numFPrimers;
        }else {
            pairedOligos = false;
            numFPrimers = oligos.getPrimers().size();
            numRPrimers = oligos.getReversePrimers().size();
        }
        
        if (oligos.getLinkers().size() != 0) { m->mothurOut("[WARNING]: pcr.seqs is not setup to remove linkers, ignoring.\n"); }
        if (oligos.getSpacers().size() != 0) { m->mothurOut("[WARNING]: pcr.seqs is not setup to remove spacers, ignoring.\n"); }
        
        return true;
    }
    catch(exception& e) {
        m->errorOut(e, "PcrSeqsCommand", "readOligos");
        exit(1);
    }
}
//********************************************************************/
bool isAligned(string seq, map<int, int>& aligned){
    aligned.clear();
    bool isAligned = false;
    
    int countBases = 0;
    for (int i = 0; i < seq.length(); i++) {
        if (!isalpha(seq[i])) { isAligned = true; }
        else { aligned[countBases] = i; countBases++; } //maps location in unaligned -> location in aligned.
    }                                                   //ie. the 3rd base may be at spot 10 in the alignment
    //later when we trim we want to trim from spot 10.
    return isAligned;
}
/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct pcrData {
    string filename;
    string oligosfile, nomatch;
    OutputWriter* goodFasta;
    OutputWriter* badFasta;
    unsigned long long fstart;
    unsigned long long fend;
    int count, start, end, length, pdiffs, pstart, pend, rdiffs;
    int numFPrimers, numRPrimers;
    MothurOut* m;
    set<string> badSeqNames;
    bool keepprimer, keepdots, fileAligned, adjustNeeded, reorient;
    Utils util;
    map<string, vector<int> > locations;
    Sequence ecoli;
    
    pcrData(){}
    pcrData(string f, string ol, OutputWriter* gf, OutputWriter* bfn, Sequence ec, string nm, bool kp, bool kd, int pd, int rd, bool re, unsigned long long fst, unsigned long long fen, int st, int en) {
        filename = f;
        goodFasta = gf;
        badFasta = bfn;
        m = MothurOut::getInstance();
        oligosfile = ol;
        start = st;
        end = en;
        ecoli = ec;
        reorient = re;
        if (ecoli.getName() != "filler") {
            length = ecoli.getAligned().length();
            start = ecoli.getStartPos();
            end = ecoli.getEndPos();
        }
        nomatch = nm;
        keepprimer = kp;
        keepdots = kd;
        fstart = fst;
        fend = fen;
        pdiffs = pd;
        rdiffs = rd;
        count = 0;
        fileAligned = true;
        adjustNeeded = false;
        pstart = -1;
        pend = -1;
    }
};
/**************************************************************************************************/
vector<TrimOligos*> fillTrims(pcrData* params, bool& pairedOligos) {
    try {
        
        vector<TrimOligos*> trims;
        
        if (params->oligosfile != "") {
            
            Oligos oligos;
            params->numFPrimers = 0; params->numRPrimers = 0;
            map<string, int> barcodes; //not used
            
            readOligos(oligos, params->oligosfile, pairedOligos, params->numFPrimers, params->numRPrimers, params->m);
            
            if (pairedOligos) {
                map<string, int> primers; vector<string> revPrimer;
                map<int, oligosPair> primerPairs = oligos.getPairedPrimers();
                for (map<int, oligosPair>::iterator it = primerPairs.begin(); it != primerPairs.end(); it++) {
                    primers[(it->second).forward] = it->first;
                    revPrimer.push_back((it->second).reverse);
                }
                
                //standard
                trims.push_back(new TrimOligos(params->pdiffs, params->rdiffs, 0, primers, barcodes, revPrimer));

            }else{
                map<string, int> primers; vector<string> revPrimer;
                primers = oligos.getPrimers();
                revPrimer = oligos.getReversePrimers();
                
                //standard
                trims.push_back(new TrimOligos(params->pdiffs, params->rdiffs, 0, primers, barcodes, revPrimer));
            }
            
                
            if (params->reorient) {
                
                //reoriented
                if (pairedOligos) {
                    map<string, int> primers; vector<string> revPrimer;
                    map<int, oligosPair> rprimerPairs = oligos.getReorientedPairedPrimers();
                    for (map<int, oligosPair>::iterator it = rprimerPairs.begin(); it != rprimerPairs.end(); it++) {
                        primers[(it->second).forward] = it->first;
                        revPrimer.push_back((it->second).reverse);
                    }
                    
                    //reoriented
                    trims.push_back(new TrimOligos(params->pdiffs, params->rdiffs, 0, primers, barcodes, revPrimer));
                    
                    primers.clear(); revPrimer.clear();
                    
                    map<int, oligosPair> revprimerPairs = oligos.getReversedPairedPrimers();
                    for (map<int, oligosPair>::iterator it = revprimerPairs.begin(); it != revprimerPairs.end(); it++) {
                        primers[(it->second).forward] = it->first;
                        revPrimer.push_back((it->second).reverse);
                    }
                    
                    //reversed
                    trims.push_back(new TrimOligos(params->pdiffs, params->rdiffs, 0, primers, barcodes, revPrimer));
                    
                }else{
                    map<string, int> primers; vector<string> revPrimer;
                    
                    primers = oligos.getReorientedPrimers();
                    revPrimer = oligos.getReorientedReversePrimers();
                    
                    //reoriented
                    trims.push_back(new TrimOligos(params->pdiffs, params->rdiffs, 0, primers, barcodes, revPrimer));
                    
                    primers.clear(); revPrimer.clear();
                    
                    primers = oligos.getReversedPrimers();
                    revPrimer = oligos.getReversedReversePrimers();
                    
                    //reversed
                    trims.push_back(new TrimOligos(params->pdiffs, params->rdiffs, 0, primers, barcodes, revPrimer));
                }
            }
        }
        
        return trims;
        
    }catch(exception& e) {
        params->m->errorOut(e, "PcrSeqsCommand", "fillTrims");
        exit(1);
    }
}
/**************************************************************************************************/
bool trimStartEnd(Sequence& seq, pcrData* params) {
    try {
        bool good = true;
        
        //make sure the seqs are aligned
        if (!params->fileAligned) { params->m->mothurOut("[ERROR]: seqs are not aligned. When using start and end your sequences must be aligned.\n"); params->m->setControl_pressed(true); good = false; }
        else {
            string alignedString = seq.getAligned();
            
            if ((seq.getStartPos() > params->start) || (seq.getEndPos() < params->end)) {
                good = false;
                if (params->m->getDebug()) {
                    params->m->mothurOut("[DEBUG]: " + seq.getName()+ " values at locations (" + toString(params->start) + "," + toString(params->end) + ") = (" + alignedString[params->start] + "," + alignedString[params->end] + ")\n");
                    
                }
            }
            else {
                if (params->end != -1) {
                    if (params->end > seq.getAligned().length()) {  params->m->mothurOut("[ERROR]: end of " +toString(params->end)  + " is longer than " + seq.getName() + " length of " +toString(seq.getAligned().length()) + ", aborting.\n"); params->m->setControl_pressed(true); good = false; }
                    else {
                        if (params->keepdots)   { seq.filterFromPos(params->end); }
                        else {
                            string seqString = seq.getAligned().substr(0, (params->end));
                            seq.setAligned(seqString);
                        }
                    }
                }
                
                if (params->start != -1) {
                    if (params->keepdots)   {  seq.filterToPos(params->start-1);  }
                    else {
                        string seqString = seq.getAligned().substr(params->start);
                        seq.setAligned(seqString);
                    }
                }
            }
        }
        
        return good;
        
    }catch(exception& e) {
        params->m->errorOut(e, "PcrSeqsCommand", "trimStartEnd");
        exit(1);
    }
}
/**************************************************************************************************/
vector<string> trimPrimers(Sequence& seq, vector<TrimOligos*> trims, vector<int>& thisSeqsLocations, int& thisPStart, int& thisPEnd, pcrData* params) {
    try {
        
        vector<string> codes; codes.resize(2, "");
    
        for (int i = 0; i < trims.size(); i++) {
            Sequence savedSeq(seq.getName(), seq.getAligned());
            
            map<int, int> mapAligned;
            bool aligned = isAligned(savedSeq.getAligned(), mapAligned);
            
            string trashCode = ""; string commentString = "";
            int currentSeqsDiffs = 0; int reverseIndex = 0; int primerIndex = 0;
            bool goodSeq = true;
            
            if(params->numFPrimers != 0){
                int primerStart = 0; int primerEnd = 0;
                vector<int> results = trims[i]->findForward(savedSeq, primerStart, primerEnd);
                bool good = true;
                if (results[0] > params->pdiffs) { good = false; }
                currentSeqsDiffs += results[0];
                commentString += "fpdiffs=" + toString(results[0]) + "(" + trims[i]->getCodeValue(results[1], params->pdiffs) + ") ";
                
                if(!good){    if (params->nomatch == "reject") { goodSeq = false; } trashCode += "f";    }
                else{
                    //are you aligned
                    if (aligned) {
                        if (!params->keepprimer)    {
                            if (params->keepdots)   { savedSeq.filterToPos(mapAligned[primerEnd-1]+1);   } //mapAligned[primerEnd-1] is the location of the last base in the primer. we want to trim to the space just after that.  The -1 & +1 ensures if the primer is followed by gaps they are not trimmed causing an aligned sequence dataset to become unaligned.
                            else{
                                savedSeq.setAligned(savedSeq.getAligned().substr(mapAligned[primerEnd-1]+1));
                                if (params->fileAligned) {
                                    thisPStart = mapAligned[primerEnd-1]+1; //locations[0].insert(mapAligned[primerEnd-1]+1);
                                    thisSeqsLocations.push_back(thisPStart);
                                }
                            }
                        }else                {
                            if (params->keepdots)   { savedSeq.filterToPos(mapAligned[primerStart]);  }
                            else            {
                                savedSeq.setAligned(savedSeq.getAligned().substr(mapAligned[primerStart]));
                                if (params->fileAligned) {
                                    thisPStart = mapAligned[primerStart]; //locations[0].insert(mapAligned[primerStart]);
                                    thisSeqsLocations.push_back(thisPStart);
                                }
                            }
                        }
                        isAligned(savedSeq.getAligned(), mapAligned);
                    }else {
                        if (!params->keepprimer)    { savedSeq.setAligned(savedSeq.getUnaligned().substr(primerEnd));     }
                        else                        { savedSeq.setAligned(savedSeq.getUnaligned().substr(primerStart));   }
                    }
                }
            }
            
            if(params->numRPrimers != 0){
                int primerStart = 0; int primerEnd = 0;
                vector<int> results = trims[i]->findReverse(savedSeq, primerStart, primerEnd);
                bool good = true;
                if (results[0] > params->rdiffs) { good = false; }
                currentSeqsDiffs += results[0];
                commentString += "rpdiffs=" + toString(results[0]) + "(" + trims[i]->getCodeValue(results[1], params->rdiffs) + ") ";
                
                if(!good){    if (params->nomatch == "reject") { goodSeq = false; } trashCode += "r";    }
                else{
                    //are you aligned
                    if (aligned) {
                        if (!params->keepprimer)    {
                            if (params->keepdots)   { savedSeq.filterFromPos(mapAligned[primerStart]); }
                            else            {
                                savedSeq.setAligned(savedSeq.getAligned().substr(0, mapAligned[primerStart]));
                                if (params->fileAligned) {
                                    thisPEnd = mapAligned[primerStart]; //locations[1].insert(mapAligned[primerStart]);
                                    thisSeqsLocations.push_back(thisPEnd);
                                }
                            }
                        }
                        else                {
                            if (params->keepdots)   { savedSeq.filterFromPos(mapAligned[primerEnd-1]+1); }
                            else            {
                                savedSeq.setAligned(savedSeq.getAligned().substr(0, mapAligned[primerEnd-1]+1));
                                if (params->fileAligned) {
                                    thisPEnd = mapAligned[primerEnd-1]+1; //locations[1].insert(mapAligned[primerEnd-1]+1);
                                    thisSeqsLocations.push_back(thisPEnd);
                                }
                            }
                        }
                    }
                    else {
                        if (!params->keepprimer)    { savedSeq.setAligned(savedSeq.getUnaligned().substr(0, primerStart));   }
                        else                { savedSeq.setAligned(savedSeq.getUnaligned().substr(0, primerEnd));     }
                    }
                }
            }
            
            if (currentSeqsDiffs > (params->pdiffs + params->rdiffs))    {    trashCode += 't';   }
            
            if (trashCode == "") {
                codes[0] = "";
                codes[1] = commentString;
               
                if (i > 0) { //checkOrient trimOligos - reoriented and reversed
                    savedSeq.reverseComplement();
                }
                seq.setAligned(savedSeq.getAligned());
                break;
            }else {
                if (codes[0] == "") { codes[0] = trashCode;                 }
                else                { codes[0] += "(" + trashCode + ")";    }
                
                codes[1] = commentString;
            }
        }
        
        return codes;

    }catch(exception& e) {
        params->m->errorOut(e, "PcrSeqsCommand", "trimPrimers");
        exit(1);
    }
}
//**********************************************************************************************************************
int driverPcr(pcrData* params){
    try {
        ifstream inFASTA; params->util.openInputFile(params->filename, inFASTA);
        inFASTA.seekg(params->fstart);
        
        bool done = false; params->count = 0;
        set<int> lengths; set<int> startLocations;  set<int> endLocations;
        
        bool pairedOligos = false;
        vector<TrimOligos*> trims = fillTrims(params, pairedOligos); //standard, if reorient parameter then (reorient & reverse) as well
        
        while (!done) {
            
            if (params->m->getControl_pressed()) {  break; }
            
            Sequence currSeq(inFASTA); params->util.gobble(inFASTA);
            
            if (params->fileAligned) { //assume aligned until proven otherwise
                lengths.insert(currSeq.getAligned().length());
                if (lengths.size() > 1) { params->fileAligned = false; }
            }
        
            if (params->m->getControl_pressed()) {  break; }
            
            string trashCode = ""; string commentString = "";
            int thisPStart = -1; int thisPEnd = -1;
            
            if (currSeq.getName() != "") {
                
                if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: seq name = " + currSeq.getName() + ".\n"); }
                
                bool goodSeq = true;
                vector<int> thisSeqsLocations;
                
                if (params->oligosfile != "") { //removing primers
                    
                    vector<string> results = trimPrimers(currSeq, trims,  thisSeqsLocations, thisPStart, thisPEnd, params);
                    
                    trashCode = results[0]; commentString = results[1];
                    
                    if (commentString != "") {
                        string seqComment = currSeq.getComment();
                        currSeq.setComment("\t" + commentString + "\t" + seqComment);
                    }
                    
                    if (trashCode != "") { goodSeq = false; }
     
                }else if (params->ecoli.getName() != "filler") {
                    //make sure the seqs are aligned
                    if (!params->fileAligned) { params->m->mothurOut("[ERROR]: seqs are not aligned. When using start and end your sequences must be aligned.\n"); params->m->setControl_pressed(true); break; }
                    else if (currSeq.getAligned().length() != params->length) {
                        params->m->mothurOut("[ERROR]: seqs are not the same length as ecoli seq. When using ecoli option your sequences must be aligned and the same length as the ecoli sequence.\n"); params->m->setControl_pressed(true); break;
                    }else {
                        if (params->keepdots)   {
                            currSeq.filterFromPos(params->end);
                            currSeq.filterToPos(params->start-1);
                        }else {
                            string seqString = currSeq.getAligned().substr(0, params->end);
                            seqString = seqString.substr(params->start);
                            currSeq.setAligned(seqString);
                        }
                    }
                }else{ //using start and end to trim
                    goodSeq = trimStartEnd(currSeq, params); //error message if seqs unaligned
                }
                
                //remove super short reads
                if (currSeq.getUnaligned() == "") { goodSeq = false;  currSeq.setAligned("NNNNNNN"); }
                
                if(goodSeq)    {
                    currSeq.printSequence(params->goodFasta);
                    if (thisPStart != -1)   { startLocations.insert(thisPStart);    }
                    if (thisPEnd != -1)     { endLocations.insert(thisPEnd);        }
                    if (thisSeqsLocations.size() != 0) { params->locations[currSeq.getName()] = thisSeqsLocations; }
                }
                else {
                    params->badSeqNames.insert(currSeq.getName());
                    currSeq.setName(currSeq.getName() + '|' + trashCode);
                    currSeq.printSequence(params->badFasta);
                }
                params->count++;
            }
            
#if defined NON_WINDOWS
            unsigned long long pos = inFASTA.tellg();
            if ((pos == -1) || (pos >= params->fend)) { break; }
#else
            if ((params->count == params->fend) || (inFASTA.eof())) { break; }
#endif
            
            //report progress
            if((params->count) % 1000 == 0){	params->m->mothurOutJustToScreen(toString(params->count)+"\n");		}
        }
        //report progress
        if((params->count) % 1000 != 0){	params->m->mothurOutJustToScreen(toString(params->count)+"\n"); 	}
        
        inFASTA.close();
        for (int i = 0; i < trims.size(); i++) { delete trims[i]; }
        
        if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: fileAligned = " + toString(params->fileAligned) +'\n'); }
        
        if (params->fileAligned && !params->keepdots) { //print out smallest start value and largest end value
            if (startLocations.size() > 1)  { params->adjustNeeded = true; }
            if (endLocations.size() > 1)    { params->adjustNeeded = true; }
            if (params->numFPrimers != 0)        {
                set<int>::iterator it = startLocations.begin();  params->pstart = *it;
                if (params->m->getDebug()) {  params->m->mothurOut("[DEBUG]: " + params->util.getStringFromSet(startLocations, " ")+"\n"); }
            }
            if (params->numRPrimers != 0)      {
                set<int>::reverse_iterator it = endLocations.rbegin();  params->pend = *it;
                if (params->m->getDebug()) {  params->m->mothurOut("[DEBUG]: " + params->util.getStringFromSet(endLocations, " ")+"\n"); }
            }
        }
        
        return params->count;
    }
    catch(exception& e) {
        params->m->errorOut(e, "PcrSeqsCommand", "driverPcr");
        exit(1);
    }
}

/**************************************************************************************************/
long long PcrSeqsCommand::createProcesses(string filename, string goodFileName, string badFileName, set<string>& badSeqNames) {
	try {
        Sequence ecoliSeq("filler","NNNN");
        if(ecolifile != "") {    ecoliSeq = readEcoli();      }  if (m->getControl_pressed()) {  return 0; }
        
        vector<double> positions;
        vector<linePair> lines;
        long long numFastaSeqs = 0;
#if defined NON_WINDOWS
        positions = util.divideFile(filename, processors);
        for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(linePair(positions[i], positions[(i+1)]));	}
#else
        
            positions = util.setFilePosFasta(fastafile, numFastaSeqs);
            if (numFastaSeqs < processors) { processors = numFastaSeqs; }
            
            //figure out how many sequences you have to process
            int numSeqsPerProcessor = numFastaSeqs / processors;
            for (int i = 0; i < processors; i++) {
                int startIndex =  i * numSeqsPerProcessor;
                if(i == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor; 	}
                lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
            }
        
#endif
        
        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<pcrData*> data;
        
        auto synchronizedGoodFastaFile = make_shared<SynchronizedOutputFile>(goodFileName);
        auto synchronizedBadFastaFile = make_shared<SynchronizedOutputFile>(badFileName);
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            OutputWriter* threadFastaWriter = new OutputWriter(synchronizedGoodFastaFile);
            OutputWriter* threadFastaScrapWriter = new OutputWriter(synchronizedBadFastaFile);
            
            pcrData* dataBundle = new pcrData(filename, oligosfile, threadFastaWriter, threadFastaScrapWriter, ecoliSeq, nomatch, keepprimer, keepdots, pdiffs, rdiffs, reorient, lines[i+1].start, lines[i+1].end, start, end);
            data.push_back(dataBundle);
            
            workerThreads.push_back(new std::thread(driverPcr, dataBundle));
        }
        
        OutputWriter* threadFastaWriter = new OutputWriter(synchronizedGoodFastaFile);
        OutputWriter* threadFastaScrapWriter = new OutputWriter(synchronizedBadFastaFile);
        
        pcrData* dataBundle = new pcrData(filename, oligosfile, threadFastaWriter, threadFastaScrapWriter, ecoliSeq, nomatch, keepprimer, keepdots, pdiffs, rdiffs, reorient, lines[0].start, lines[0].end, start, end);
        
        driverPcr(dataBundle);
        numFastaSeqs = dataBundle->count;
        
        badSeqNames = dataBundle->badSeqNames;
        map<string, vector<int> > locations = dataBundle->locations;
        bool adjustNeeded = dataBundle->adjustNeeded;
        int pstart = -1; int pend = -1;
        pstart = dataBundle->pstart; pend = dataBundle->pend;
        bool hasFPrimers = false; if (dataBundle->numFPrimers != 0) { hasFPrimers = true;  }
        bool hasRPrimers = false; if (dataBundle->numRPrimers != 0) { hasRPrimers = true;  }
        
        if (m->getDebug()) {  m->mothurOut("[DEBUG]: pstart = " + toString(pstart) + ", pend = " + toString(pend) + "\n");  }
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            numFastaSeqs += data[i]->count;
            
            delete data[i]->goodFasta;
            delete data[i]->badFasta;
            if (data[i]->adjustNeeded) { adjustNeeded = true; }
            if (data[i]->pstart != -1)   {
                if (data[i]->pstart != pstart)  { adjustNeeded = true;          }
                if (data[i]->pstart < pstart)   { pstart = data[i]->pstart;     }
            } //smallest start
            
            if (data[i]->pend != -1)   {
                if (data[i]->pend != pend)      { adjustNeeded = true;          }
                if (data[i]->pend > pend)       { pend = data[i]->pend;         }
            }//largest end
            
            if (m->getDebug()) {  m->mothurOut("[DEBUG]: process " + toString(i) + " pstart = " + toString(data[i]->pstart) + ", pend = " + toString(data[i]->pend) + "\n");  }
            
            badSeqNames.insert(data[i]->badSeqNames.begin(), data[i]->badSeqNames.end());
            locations.insert(data[i]->locations.begin(), data[i]->locations.end());
            
            delete data[i];
            delete workerThreads[i];
        }
        synchronizedGoodFastaFile->close(); //must explicitly close or file may still be open when reading and writing in adjustDots
        synchronizedBadFastaFile->close(); //must explicitly close or file may still be open when reading and writing in adjustDots
        delete threadFastaWriter;
        delete threadFastaScrapWriter;
        delete dataBundle;
        
        if (m->getDebug()) {  m->mothurOut("[DEBUG]: pstart = " + toString(pstart) + ", pend = " + toString(pend) + "\n");  }

        if (fileAligned && adjustNeeded) {
            //find pend - pend is the biggest ending value, but we must account for when we adjust the start.  That adjustment may make the "new" end larger then the largest end. So lets find out what that "new" end will be.
            for (map<string, vector<int> >::iterator it = locations.begin(); it != locations.end(); it++) {
                if (m->getControl_pressed()) { break; }
                
                string name = it->first;
                int thisStart = -1; int thisEnd = -1;
                if (hasFPrimers)       { thisStart = it->second[0];    }
                if (hasRPrimers)       { thisEnd = it->second[1];      }
                else { pend = -1; break; }
                
                int myDiff = 0;
                if (thisStart != -1) { //my start
                    if (thisStart != pstart) { //my start is after where the first start occurs, so I need to pad in the front
                        myDiff += (thisStart - pstart); //size of my pad
                    }
                }
                
                int myEnd = thisEnd + myDiff;
                if (thisEnd != -1) {
                    if (myEnd > pend) { pend = myEnd; }
                }
            }
            
            if (m->getDebug()) {  m->mothurOut("[DEBUG]: pstart = " + toString(pstart) + ", pend = " + toString(pend) + "\n");  }
            
            adjustDots(goodFileName, locations, pstart, pend, hasFPrimers, hasRPrimers); 
        }
        
        return numFastaSeqs;
	}
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "createProcesses");
		exit(1);
	}
}
//**********************************************************************************************************************
int PcrSeqsCommand::adjustDots(string goodFasta, map<string, vector<int> > locations, int pstart, int pend, bool hasFPrimers, bool hasRPrimers){
    try {
        ifstream inFasta;
        util.openInputFile(goodFasta, inFasta);
        
        ofstream out;
        util.openOutputFile(goodFasta+".temp", out);
        
        set<int> lengths;
        while(!inFasta.eof()) {
            if(m->getControl_pressed()) { break; }
            
            Sequence seq(inFasta); util.gobble(inFasta);
            
            string name = seq.getName();
            int thisStart = -1; int thisEnd = -1;
            map<string, vector<int> >::iterator it = locations.find(name);
            if (it != locations.end()) {
                if (hasFPrimers)       { thisStart = it->second[0];    }
                if (hasRPrimers)       { thisEnd = it->second[1];      }
            }else { m->mothurOut("[ERROR]: should never get here in pcr.seqs.\n"); }
            
            if (name != seq.getName()) { m->mothurOut("[ERROR]: name mismatch in pcr.seqs.\n"); }
            else {
                
                string forwardPad = "";  string reversePad = "";
                
                if ((pstart != -1) && (thisStart != -1) && (thisStart != pstart)) {
                    for (int i = pstart; i < thisStart; i++) { forwardPad += "."; }
                    thisEnd += forwardPad.length();
                }
                
                if ((pend != -1) && (thisEnd != -1) && (thisEnd != pend)) {
                    for (int i = thisEnd; i < pend; i++) { reversePad += "."; }
                }
                
                string aligned = forwardPad + seq.getAligned() + reversePad;
                seq.setAligned(aligned);
                lengths.insert(seq.getAligned().length());
            }
            
            seq.printSequence(out);
        }
        
        inFasta.close();
        out.close();
        util.mothurRemove(goodFasta);
        util.renameFile(goodFasta+".temp", goodFasta);
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "PcrSeqsCommand", "adjustDots");
        exit(1);
    }
}
//***************************************************************************************************************
Sequence PcrSeqsCommand::readEcoli(){
	try {
		ifstream in; util.openInputFile(ecolifile, in);
		
        //read seq
        Sequence result;
        if (!in.eof()){ 
            Sequence ecoli(in);
            result.setName(ecoli.getName()); result.setAligned(ecoli.getAligned());
        }else {  m->setControl_pressed(true);  }
        in.close();    
			
        return result;
    }
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "readEcoli");
		exit(1);
	}
    
}
//***************************************************************************************************************
int PcrSeqsCommand::writeAccnos(set<string> badNames, string outputFileName){
	try {
        ofstream out;
        util.openOutputFile(outputFileName, out);
        
        for (set<string>::iterator it = badNames.begin(); it != badNames.end(); it++) {
            if (m->getControl_pressed()) { break; }
            out << (*it) << endl;
        }
        
        out.close();
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "PcrSeqsCommand", "writeAccnos");
		exit(1);
	}
    
}
/**************************************************************************************/


