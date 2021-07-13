/*
 *  distancecommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "distancecommand.h"

//**********************************************************************************************************************
vector<string> DistanceCommand::setParameters(){	
	try {
		CommandParameter pcolumn("column", "InputTypes", "", "", "none", "none", "OldFastaColumn","column",false,false); parameters.push_back(pcolumn);
		CommandParameter poldfasta("oldfasta", "InputTypes", "", "", "none", "none", "OldFastaColumn","",false,false); parameters.push_back(poldfasta);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","phylip-column",false,true, true); parameters.push_back(pfasta);
		CommandParameter poutput("output", "Multiple", "column-lt-square", "column", "", "", "","phylip-column",false,false, true); parameters.push_back(poutput);
		CommandParameter pcalc("calc", "Multiple", "nogaps-eachgap-onegap-jtt-pmb-pam-kimura", "onegap", "", "", "","",false,false); parameters.push_back(pcalc);
		CommandParameter pcountends("countends", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pcountends);
        CommandParameter pfitcalc("fitcalc", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pfitcalc);
		CommandParameter pcompress("compress", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pcompress);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false, true); parameters.push_back(pprocessors);
		CommandParameter pcutoff("cutoff", "Number", "", "1.0", "", "", "","",false,false, true); parameters.push_back(pcutoff);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
        abort = false; calledHelp = false;
       
        vector<string> tempOutNames;
        outputTypes["phylip"] = tempOutNames;
        outputTypes["column"] = tempOutNames;
        
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string DistanceCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The dist.seqs command reads a file containing sequences and creates a distance file.\n";
		helpString += "The dist.seqs command parameters are fasta, oldfasta, column, calc, countends, output, compress, cutoff and processors.  \n";
		helpString += "The fasta parameter is required, unless you have a valid current fasta file.\n";
		helpString += "The oldfasta and column parameters allow you to append the distances calculated to the column file.\n";
		helpString += "The calc parameter allows you to specify the method of calculating the distances.  Your options are: nogaps, onegap or eachgap for dna/rna sequences. If using protein sequences, your calc options are jtt, pmb, pam and kimura. The default is onegap.\n";
		helpString += "The countends parameter allows you to specify whether to include terminal gaps in distance.  Your options are: T or F. The default is T.\n";
		helpString += "The cutoff parameter allows you to specify maximum distance to keep. The default is 1.0.\n";
		helpString += "The output parameter allows you to specify format of your distance matrix. Options are column, lt, and square. The default is column.\n";
		helpString += "The processors parameter allows you to specify number of processors to use.  The default is 1.\n";
		helpString += "The compress parameter allows you to indicate that you want the resulting distance file compressed.  The default is false.\n";
		helpString += "The dist.seqs command should be in the following format: \n";
		helpString += "dist.seqs(fasta=yourFastaFile, calc=yourCalc, countends=yourEnds, cutoff= yourCutOff, processors=yourProcessors) \n";
		helpString += "Example dist.seqs(fasta=amazon.fasta, calc=eachgap, countends=F, cutoff= 2.0, processors=3).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string DistanceCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "phylip") {  pattern = "[filename],[outputtag],dist"; } 
        else if (type == "column") { pattern = "[filename],dist-[filename],[outputtag],dist"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "DistanceCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
DistanceCommand::DistanceCommand(string option) {
	try {
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not found") { 				
				fastafile = current->getFastaFile(); 
                if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter.\n");  }
                else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required.\n"); abort = true; }
			}else if (fastafile == "not open") { abort = true; }	
			else{ current->setFastaFile(fastafile); }
			
			oldfastafile = validParameter.validFile(parameters, "oldfasta");
			if (oldfastafile == "not found") { oldfastafile = ""; }
			else if (oldfastafile == "not open") { abort = true; }	
			
			column = validParameter.validFile(parameters, "column");
			if (column == "not found") { column = ""; }
			else if (column == "not open") { abort = true; }	
			else { current->setColumnFile(column); }
			
            if (outputdir == ""){ outputdir += util.hasPath(fastafile);  }

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			calc = validParameter.valid(parameters, "calc");
			if (calc == "not found") { calc = "onegap";  }
			else { 
				 if (calc == "default")  {  calc = "onegap";  }
            }

			string temp;
			temp = validParameter.valid(parameters, "countends");	if(temp == "not found"){	temp = "T";	}
            countends = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "fitcalc");	if(temp == "not found"){	temp = "F";	}
            fitCalc = util.isTrue(temp);
			
			temp = validParameter.valid(parameters, "cutoff");		if(temp == "not found"){	temp = "1.0"; }
			util.mothurConvert(temp, cutoff); 
			
			temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
			
			temp = validParameter.valid(parameters, "compress");		if(temp == "not found"){  temp = "F"; }
			convert(temp, compress);

			output = validParameter.valid(parameters, "output");		if(output == "not found"){	output = "column"; }
            if (output == "phylip") { output = "lt";  }
            
			if (((column != "") && (oldfastafile == "")) || ((column == "") && (oldfastafile != ""))) { m->mothurOut("If you provide column or oldfasta, you must provide both.\n");  abort=true; }
			
			if ((column != "") && (oldfastafile != "") && (output != "column")) { m->mothurOut("You have provided column and oldfasta, indicating you want to append distances to your column file. Your output must be in column format to do so.\n"); abort=true; }
			
			if ((output != "column") && (output != "lt") && (output != "square")) { m->mothurOut(output + " is not a valid output form. Options are column, lt and square. I will use column.\n");  output = "column"; }
            
            if ((calc != "onegap") && (calc != "eachgap") && (calc != "nogaps") && (calc != "jtt") && (calc != "pmb") && (calc != "pam") && (calc != "kimura")) { m->mothurOut(calc + " is not a valid calc. Options are eachgap, onegap, nogaps, jtt, pmb, pam and kimura. I'll use onegap.\n");  calc = "onegap";  }
            
            prot = false; //not using protein sequences
            if ((calc == "jtt") || (calc == "pmb") || (calc == "pam") || (calc == "kimura")) { prot = true; }

		}
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "DistanceCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
DistanceCommand::DistanceCommand(StorageDatabase*& storageDB, string outputFileRoot, double cut, string outputformat, int proc) {
    try {
        abort = false; calledHelp = false;
        vector<string> tempOutNames;
        outputTypes["phylip"] = tempOutNames;
        outputTypes["column"] = tempOutNames;
        
        calc = "onegap";
        countends = true;
        fitCalc = false;
        cutoff = cut;
        processors = proc;
        compress = false;
        output = outputformat;
        prot = false; //not using protein sequences
            
        numDistsBelowCutoff = 0;
        db = storageDB;
        numNewFasta = db->getNumSeqs();
        numSeqs = db->getNumSeqs();
        
        if (!db->sameLength()) {  m->mothurOut("[ERROR]: your sequences are not the same length, aborting.\n");  return; }
        if (numSeqs < 2) {  m->mothurOut("[ERROR]: you must have at least 2 sequences to calculate the distances, aborting.\n");  return; }
        
        string outputFile;
        map<string, string> variables;
        variables["[filename]"] = outputFileRoot;
        
        if (output == "lt") { //does the user want lower triangle phylip formatted file
            variables["[outputtag]"] = "phylip";
            outputFile = getOutputFileName("phylip", variables);
            util.mothurRemove(outputFile); outputTypes["phylip"].push_back(outputFile);
            
            //output numSeqs to phylip formatted dist file
        }else if (output == "column") { //user wants column format
            
            outputFile = getOutputFileName("column", variables);
            outputTypes["column"].push_back(outputFile);
            util.mothurRemove(outputFile);
        }
       
        
        m->mothurOut("\nSequence\tTime\tNum_Dists_Below_Cutoff\n");
                     
        createProcesses(outputFile);
        
        m->mothurOut("\nOutput File Names:\n"); m->mothurOut(outputFile+"\n\n");
        
    }
    catch(exception& e) {
        m->errorOut(e, "DistanceCommand", "DistanceCommand");
        exit(1);
    }
}
//**********************************************************************************************************************

int DistanceCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
        numDistsBelowCutoff = 0;
        
        ifstream inFASTA; util.openInputFile(fastafile, inFASTA);
        if (prot) { db = new ProteinDB(inFASTA);  }
        else      { db = new SequenceDB(inFASTA); }
        inFASTA.close();
		
		//save number of new sequence
		numNewFasta = db->getNumSeqs();
        
		//sanity check the oldfasta and column file as well as add oldfasta sequences to db
        if ((oldfastafile != "") && (column != ""))  {	if (!(sanityCheck())) { return 0; }  }
		
        if (m->getControl_pressed()) { delete db; return 0; }
		
		numSeqs = db->getNumSeqs();
		
		if (!db->sameLength()) {  m->mothurOut("[ERROR]: your sequences are not the same length, aborting.\n");  return 0; }
        if (numSeqs < 2) {  m->mothurOut("[ERROR]: you must have at least 2 sequences to calculate the distances, aborting.\n");  return 0; }

        
		string outputFile;
        map<string, string> variables; 
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(fastafile));
        if ((oldfastafile != "") && (column != ""))  { variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(oldfastafile)); }
            
		if (output == "lt") { //does the user want lower triangle phylip formatted file 
            variables["[outputtag]"] = "phylip";
			outputFile = getOutputFileName("phylip", variables);
			util.mothurRemove(outputFile); outputTypes["phylip"].push_back(outputFile);
			
			//output numSeqs to phylip formatted dist file
		}else if (output == "column") { //user wants column format
            if (fitCalc) {  variables["[outputtag]"] = "fit";  }
			outputFile = getOutputFileName("column", variables);
			outputTypes["column"].push_back(outputFile);
            
			//so we don't accidentally overwrite
			if (outputFile == column) { 
				string tempcolumn = column + ".old"; 
				rename(column.c_str(), tempcolumn.c_str());
			}
			
			util.mothurRemove(outputFile);
		}else { //assume square
			variables["[outputtag]"] = "square";
			outputFile = getOutputFileName("phylip", variables);
			util.mothurRemove(outputFile);
			outputTypes["phylip"].push_back(outputFile);
		}
        m->mothurOut("\nSequence\tTime\tNum_Dists_Below_Cutoff\n");
                     
        createProcesses(outputFile);
        
		if (m->getControl_pressed()) { outputTypes.clear();  util.mothurRemove(outputFile); return 0; }
		
		ifstream fileHandle; fileHandle.open(outputFile.c_str());
		if(fileHandle) {
			util.gobble(fileHandle);
			if (fileHandle.eof()) { m->mothurOut(outputFile + " is blank. This can result if there are no distances below your cutoff.\n"); }
		}
		
		//append the old column file to the new one
		if ((oldfastafile != "") && (column != ""))  {
			//we had to rename the column file so we didnt overwrite above, but we want to keep old name
			if (outputFile == column) { 
				string tempcolumn = column + ".old";
				util.appendFiles(tempcolumn, outputFile);
				util.mothurRemove(tempcolumn);
			}else{
                if (!fitCalc) {
                    util.appendFiles(outputFile, column);
                    util.mothurRemove(outputFile);
                    outputFile = column;
                }
			}
            outputTypes["column"].clear(); outputTypes["column"].push_back(outputFile);
        }
		
		if (m->getControl_pressed()) { outputTypes.clear();  util.mothurRemove(outputFile); return 0; }
		
		//set phylip file as new current phylipfile
		string currentName = "";
		itTypes = outputTypes.find("phylip");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setPhylipFile(currentName); }
		}
		
		//set column file as new current columnfile
		itTypes = outputTypes.find("column");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setColumnFile(currentName); }
		}
		
		m->mothurOut("\nOutput File Names: \n");
		m->mothurOut(outputFile+"\n\n");
		
		if (util.isTrue(compress)) {
			m->mothurOut("Compressing...\n"); 
			m->mothurOut("(Replacing " + outputFile + " with " + outputFile + ".gz)\n"); 
			system(("gzip -v " + outputFile).c_str());
			outputNames.push_back(outputFile + ".gz");
		}else { outputNames.push_back(outputFile); }

		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/
void driverColumn(distanceData* params){
    try {
        ValidCalculators validCalculator;
        DistCalc* distCalculator;

        if (!params->prot) {
            if (params->countends) {
                if (validCalculator.isValidCalculator("distance", params->calc) ) {
                    if (params->calc == "nogaps")			{	distCalculator = new ignoreGaps(params->cutoff);	}
                    else if (params->calc == "eachgap")	{	distCalculator = new eachGapDist(params->cutoff);	}
                    else if (params->calc == "onegap")		{	distCalculator = new oneGapDist(params->cutoff);	}
                }
            }else {
                if (validCalculator.isValidCalculator("distance", params->calc) ) {
                    if (params->calc == "nogaps")		{	distCalculator = new ignoreGaps(params->cutoff);					}
                    else if (params->calc == "eachgap"){	distCalculator = new eachGapIgnoreTermGapDist(params->cutoff);	}
                    else if (params->calc == "onegap")	{	distCalculator = new oneGapIgnoreTermGapDist(params->cutoff);		}
                }
            }
        }else {
            if (validCalculator.isValidCalculator("protdist", params->calc) ) {
                if (params->calc == "jtt")        {    distCalculator = new JTT(params->cutoff);                    }
                else if (params->calc == "pmb")        {    distCalculator = new PMB(params->cutoff);               }
                else if (params->calc == "pam")        {    distCalculator = new PAM(params->cutoff);               }
                else if (params->calc == "kimura")        {    distCalculator = new Kimura(params->cutoff);               }
            }
        }
            
        
        int startTime = time(NULL);
       
        params->count = 0;
        string buffer = "";
        
        for(int i=params->startLine;i<params->endLine;i++){
            
            Sequence seqI; Protein seqIP; string nameI = "";
            if (params->prot)   { seqIP = params->db->getProt(i);   nameI = seqIP.getName();    }
            else                { seqI = params->db->getSeq(i);     nameI = seqI.getName();     }
            
            for(int j=0;j<i;j++){
                
                if (params->m->getControl_pressed()) { break;  }
                
                if ((i >= params->numNewFasta) && (j >= params->numNewFasta)) { break; }
                
                double dist = 1.0; string nameJ = "";
                if (params->prot)   { Protein seqJP = params->db->getProt(j); nameJ = seqJP.getName(); dist = distCalculator->calcDist(seqIP, seqJP);   }
                else                { Sequence seqJ = params->db->getSeq(j);  nameJ = seqJ.getName(); dist = distCalculator->calcDist(seqI, seqJ);       }
                
                
                if(dist <= params->cutoff){
                    buffer += (nameI + " " + nameJ + " " + toString(dist) + "\n");
                    params->count++;
                }
            }
            
            if(i % 100 == 0){ params->threadWriter->write(buffer);  buffer = ""; params->m->mothurOutJustToScreen(toString(i) + "\t" + toString(time(NULL) - startTime) + "\t" + toString(params->count) +"\n"); }
            
        }
        params->threadWriter->write(buffer);
        
        if((params->endLine-1) % 100 != 0){ params->m->mothurOutJustToScreen(toString(params->endLine-1) + "\t" + toString(time(NULL) - startTime) + "\t" + toString(params->count) +"\n"); }
        
        delete distCalculator;
    }
    catch(exception& e) {
        params->m->errorOut(e, "DistanceCommand", "driverColumn");
        exit(1);
    }
}
/**************************************************************************************************/
void driverLt(distanceData* params){
    try {
        ValidCalculators validCalculator;
        DistCalc* distCalculator;
        double cutoff = 1.0;
        
        if (!params->prot) {
            if (params->countends) {
                if (validCalculator.isValidCalculator("distance", params->calc) ) {
                    if (params->calc == "nogaps")			{	distCalculator = new ignoreGaps(cutoff);	}
                    else if (params->calc == "eachgap")	{	distCalculator = new eachGapDist(cutoff);	}
                    else if (params->calc == "onegap")		{	distCalculator = new oneGapDist(cutoff);	}
                }
            }else {
                if (validCalculator.isValidCalculator("distance", params->calc) ) {
                    if (params->calc == "nogaps")		{	distCalculator = new ignoreGaps(cutoff);					}
                    else if (params->calc == "eachgap"){	distCalculator = new eachGapIgnoreTermGapDist(cutoff);	}
                    else if (params->calc == "onegap")	{	distCalculator = new oneGapIgnoreTermGapDist(cutoff);		}
                }
            }
        }else {
            if (validCalculator.isValidCalculator("protdist", params->calc) ) {
                if (params->calc == "jtt")        {    distCalculator = new JTT(params->cutoff);                    }
                else if (params->calc == "pmb")        {    distCalculator = new PMB(params->cutoff);               }
                else if (params->calc == "pam")        {    distCalculator = new PAM(params->cutoff);               }
                else if (params->calc == "kimura")        {    distCalculator = new Kimura(params->cutoff);         }
            }
        }
        
        int startTime = time(NULL);
        long long numSeqs = params->db->getNumSeqs();
        
        //column file
        ofstream outFile; params->util.openOutputFile(params->outputFileName, outFile);
        outFile.setf(ios::fixed, ios::showpoint); outFile << setprecision(4);
        
        if(params->startLine == 0){	outFile << numSeqs << endl;	}
        
        params->count = 0;
        for(int i=params->startLine;i<params->endLine;i++){
            
            Sequence seqI; Protein seqIP; string nameI = "";
            if (params->prot)   { seqIP = params->db->getProt(i);   nameI = seqIP.getName();    }
            else                { seqI = params->db->getSeq(i);     nameI = seqI.getName();     }
            
            if (nameI.length() < 10) {  while (nameI.length() < 10) {  nameI += " ";  } }
            outFile << nameI;
            
            for(int j=0;j<i;j++){
                
                if (params->m->getControl_pressed()) { break;  }
                
                if ((i >= params->numNewFasta) && (j >= params->numNewFasta)) { break; }
                
                double dist = 1.0;
                if (params->prot)   { Protein seqJP = params->db->getProt(j);  dist = distCalculator->calcDist(seqIP, seqJP);   }
                else                { Sequence seqJ = params->db->getSeq(j);  dist = distCalculator->calcDist(seqI, seqJ);       }
                
                if(dist <= params->cutoff){ params->count++; }
                outFile  << '\t' << dist;
            }
            
            outFile << endl;
            
            if(i % 100 == 0){ params->m->mothurOutJustToScreen(toString(i) + "\t" + toString(time(NULL) - startTime) + "\t" + toString(params->count) +"\n"); }
            
        }
        
        if((params->endLine-1) % 100 != 0){ params->m->mothurOutJustToScreen(toString(params->endLine-1) + "\t" + toString(time(NULL) - startTime) + "\t" + toString(params->count) +"\n"); }
        
        outFile.close();
        delete distCalculator;
    }
    catch(exception& e) {
        params->m->errorOut(e, "DistanceCommand", "driverLt");
        exit(1);
    }
}
/**************************************************************************************************/
void driverSquare(distanceData* params){
    try {
        ValidCalculators validCalculator;
        DistCalc* distCalculator;
        double cutoff = 1.0;
        
        if (!params->prot) {
            if (params->countends) {
                if (validCalculator.isValidCalculator("distance", params->calc) ) {
                    if (params->calc == "nogaps")			{	distCalculator = new ignoreGaps(cutoff);	}
                    else if (params->calc == "eachgap")	{	distCalculator = new eachGapDist(cutoff);	    }
                    else if (params->calc == "onegap")		{	distCalculator = new oneGapDist(cutoff);	}
                }
            }else {
                if (validCalculator.isValidCalculator("distance", params->calc) ) {
                    if (params->calc == "nogaps")		{	distCalculator = new ignoreGaps(cutoff);					}
                    else if (params->calc == "eachgap"){	distCalculator = new eachGapIgnoreTermGapDist(cutoff);	}
                    else if (params->calc == "onegap")	{	distCalculator = new oneGapIgnoreTermGapDist(cutoff);		}
                }
            }
        }else {
            if (validCalculator.isValidCalculator("protdist", params->calc) ) {
                if (params->calc == "jtt")        {    distCalculator = new JTT(params->cutoff);                    }
                else if (params->calc == "pmb")        {    distCalculator = new PMB(params->cutoff);               }
                else if (params->calc == "pam")        {    distCalculator = new PAM(params->cutoff);               }
                else if (params->calc == "kimura")        {    distCalculator = new Kimura(params->cutoff);         }
            }
        }
        
        int startTime = time(NULL);
        
        //column file
        ofstream outFile;
        params->util.openOutputFile(params->outputFileName, outFile);
        outFile.setf(ios::fixed, ios::showpoint); outFile << setprecision(4);
        
        long long numSeqs = params->db->getNumSeqs();
        if(params->startLine == 0){	outFile << numSeqs << endl;	}
        
        params->count = 0;
        for(int i=params->startLine;i<params->endLine;i++){
            
            Sequence seqI; Protein seqIP; string nameI = "";
            if (params->prot)   { seqIP = params->db->getProt(i);   nameI = seqIP.getName();    }
            else                { seqI = params->db->getSeq(i);     nameI = seqI.getName();     }
            
            if (nameI.length() < 10) {  while (nameI.length() < 10) {  nameI += " ";  } }
            outFile << nameI << '\t';
            
            for(int j=0;j<numSeqs;j++){
                
                if (params->m->getControl_pressed()) { break; }
                
                double dist = 1.0;
                if (i == j) { dist = 0.0000; }
                else {
                    if (params->prot)   { Protein seqJP = params->db->getProt(j);  dist = distCalculator->calcDist(seqIP, seqJP);   }
                    else                { Sequence seqJ = params->db->getSeq(j);  dist = distCalculator->calcDist(seqI, seqJ);       }
                }
                
                if(dist <= params->cutoff){ params->count++; }
                
                outFile << dist << '\t'; 
            }
            
            outFile << endl; 
            
            if(i % 100 == 0){ params->m->mothurOutJustToScreen(toString(i) + "\t" + toString(time(NULL) - startTime) + "\t" + toString(params->count) +"\n");  }
        }
        
        if((params->endLine-1) % 100 != 0){ params->m->mothurOutJustToScreen(toString(params->endLine-1) + "\t" + toString(time(NULL) - startTime) + "\t" + toString(params->count) +"\n"); }
        
        outFile.close();
        delete distCalculator;
    }
    catch(exception& e) {
        params->m->errorOut(e, "DistanceCommand", "driverSquare");
        exit(1);
    }
}
/**************************************************************************************************/
void driverFitCalc(distanceData* params){
    try {
        ValidCalculators validCalculator;
        DistCalc* distCalculator;
        
        if (!params->prot) {
            if (params->countends) {
                if (validCalculator.isValidCalculator("distance", params->calc) ) {
                    if (params->calc == "nogaps")			{	distCalculator = new ignoreGaps(params->cutoff);	}
                    else if (params->calc == "eachgap")	{	distCalculator = new eachGapDist(params->cutoff);	}
                    else if (params->calc == "onegap")		{	distCalculator = new oneGapDist(params->cutoff);	}
                }
            }else {
                if (validCalculator.isValidCalculator("distance", params->calc) ) {
                    if (params->calc == "nogaps")		{	distCalculator = new ignoreGaps(params->cutoff);					}
                    else if (params->calc == "eachgap"){	distCalculator = new eachGapIgnoreTermGapDist(params->cutoff);	}
                    else if (params->calc == "onegap")	{	distCalculator = new oneGapIgnoreTermGapDist(params->cutoff);		}
                }
            }
        }else {
            if (validCalculator.isValidCalculator("protdist", params->calc) ) {
                if (params->calc == "jtt")        {    distCalculator = new JTT(params->cutoff);                    }
                else if (params->calc == "pmb")        {    distCalculator = new PMB(params->cutoff);               }
                else if (params->calc == "pam")        {    distCalculator = new PAM(params->cutoff);               }
                else if (params->calc == "kimura")        {    distCalculator = new Kimura(params->cutoff);         }
            }
        }
        
        int startTime = time(NULL);
        params->count = 0;
        string buffer = "";
        for(int i=params->startLine;i<params->endLine;i++){
            
            Sequence seqI; Protein seqIP; string nameI = "";
            if (params->prot)   { seqIP = params->oldFastaDB->getProt(i);   nameI = seqIP.getName();    }
            else                { seqI = params->oldFastaDB->getSeq(i);     nameI = seqI.getName();     }
            
            
            for(int j = 0; j < params->db->getNumSeqs(); j++){
                
                if (params->m->getControl_pressed()) { break;  }
                
                double dist = 1.0; string nameJ = "";
                
                if (params->prot)   { Protein seqJP = params->db->getProt(j); nameJ = seqJP.getName(); dist = distCalculator->calcDist(seqIP, seqJP);   }
                else                { Sequence seqJ = params->db->getSeq(j);  nameJ = seqJ.getName(); dist = distCalculator->calcDist(seqI, seqJ);       }
                
                if(dist <= params->cutoff){
                    buffer += nameI + " " + nameJ + " " + toString(dist) + "\n";
                    params->count++;
                }
            }
            
            if(i % 100 == 0){ params->threadWriter->write(buffer);  buffer = ""; params->m->mothurOutJustToScreen(toString(i) + "\t" + toString(time(NULL) - startTime) + "\t" + toString(params->count) +"\n"); }
            
        }
        params->threadWriter->write(buffer);
        
        if((params->endLine-1) % 100 != 0){ params->m->mothurOutJustToScreen(toString(params->endLine-1) + "\t" + toString(time(NULL) - startTime) + "\t" + toString(params->count) +"\n"); }
        
        delete distCalculator;

    }
    catch(exception& e) {
        params->m->errorOut(e, "DistanceCommand", "driverFitCalc");
        exit(1);
    }
}
/**************************************************************************************************/
void DistanceCommand::createProcesses(string filename) {
    try {
        long long num = db->getNumSeqs();
        long long distsBelowCutoff = 0;
        time_t start, end;
        time(&start);
        
        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<distanceData*> data;
        
        double numDists = 0;
        
        if (output == "square") { numDists = numSeqs; }
        else { for(int i=0;i<numSeqs;i++){ for(int j=0;j<i;j++){ numDists++; if (numDists > processors) { break; } } } }
        if (numDists < processors) { processors = numDists; }
        
        vector<linePair> lines;
        for (int i = 0; i < processors; i++) {
            linePair tempLine;
            lines.push_back(tempLine);
            if (output != "square") {
                lines[i].start = int (sqrt(float(i)/float(processors)) * numSeqs);
                lines[i].end = int (sqrt(float(i+1)/float(processors)) * numSeqs);
            }else{
                lines[i].start = int ((float(i)/float(processors)) * numSeqs);
                lines[i].end = int ((float(i+1)/float(processors)) * numSeqs);
            }
        }
        
        auto synchronizedOutputFile = std::make_shared<SynchronizedOutputFile>(filename);
        synchronizedOutputFile->setFixedShowPoint(); synchronizedOutputFile->setPrecision(4);
        
        StorageDatabase* oldFastaDB;
        if (fitCalc) {
            ifstream inFASTA; util.openInputFile(oldfastafile, inFASTA);
            if (!prot) { oldFastaDB = new SequenceDB(inFASTA); }
            else                    { oldFastaDB = new ProteinDB(inFASTA); }
            inFASTA.close();
            
            lines.clear();
            if (processors > oldFastaDB->getNumSeqs()) { processors = oldFastaDB->getNumSeqs(); }
            int remainingSeqs = oldFastaDB->getNumSeqs();
            int startIndex = 0;
            for (int remainingProcessors = processors; remainingProcessors > 0; remainingProcessors--) {
                int numSeqsToFit = remainingSeqs; //case for last processor
                if (remainingProcessors != 1) { numSeqsToFit = ceil(remainingSeqs / remainingProcessors); }
                lines.push_back(linePair(startIndex, (startIndex+numSeqsToFit))); //startIndex, endIndex
                startIndex = startIndex + numSeqsToFit;
                remainingSeqs -= numSeqsToFit;
            }
        }
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            OutputWriter* threadWriter = NULL;
            distanceData* dataBundle = NULL;
            string extension = toString(i+1) + ".temp";
            if (output == "column") {
                threadWriter = new OutputWriter(synchronizedOutputFile);
                dataBundle = new distanceData(threadWriter);
            }else { dataBundle = new distanceData(filename+extension); }
            dataBundle->setVariables(lines[i+1].start, lines[i+1].end, cutoff, db, oldFastaDB, calc, prot, numNewFasta, countends);
            data.push_back(dataBundle);
            
            std::thread* thisThread = NULL;
            if (output == "column")     {
                if (fitCalc)    { thisThread = new std::thread(driverFitCalc, dataBundle);   }
                else            {  thisThread = new std::thread(driverColumn, dataBundle);   }
            }
            else if (output == "lt")    { thisThread = new std::thread(driverLt, dataBundle);            }
            else                        { thisThread = new std::thread(driverSquare, dataBundle);        }
            workerThreads.push_back(thisThread);
        }
        
        OutputWriter* threadWriter = NULL;
        distanceData* dataBundle = NULL;
        if (output == "column") {
            threadWriter = new OutputWriter(synchronizedOutputFile);
            dataBundle = new distanceData(threadWriter);
        }else { dataBundle = new distanceData(filename); }
        dataBundle->setVariables(lines[0].start, lines[0].end, cutoff, db, oldFastaDB, calc, prot, numNewFasta, countends);
        
        if (output == "column")     {
            if (fitCalc)    { driverFitCalc(dataBundle);    }
            else            { driverColumn(dataBundle);     }
        }
        else if (output == "lt")    { driverLt(dataBundle);            }
        else                        { driverSquare(dataBundle);        }
        distsBelowCutoff = dataBundle->count;
        
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            
            distsBelowCutoff += data[i]->count;
            if (output == "column") {  delete data[i]->threadWriter; }
            else {
                string extension = toString(i+1) + ".temp";
                util.appendFiles((filename+extension), filename);
                util.mothurRemove(filename+extension);
            }
            delete data[i];
            delete workerThreads[i];
        }
        if (output == "column")     { synchronizedOutputFile->close(); delete threadWriter; }
        delete dataBundle;
        
        time(&end);
        m->mothurOut("\nIt took " + toString(difftime(end, start)) + " secs to find distances for " + toString(num) + " sequences. " + toString(distsBelowCutoff+numDistsBelowCutoff) + " distances below cutoff " + toString(cutoff) + ".\n\n");
        
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "createProcesses");
		exit(1);
	}
}

/**************************************************************************************************/
//its okay if the column file does not contain all the names in the fasta file, since some distance may have been above a cutoff,
//but no sequences can be in the column file that are not in oldfasta. also, if a distance is above the cutoff given then remove it.
//also check to make sure the 2 files have the same alignment length.
bool DistanceCommand::sanityCheck() {
	try{
		bool good = true;
		
		//make sure the 2 fasta files have the same alignment length
		ifstream in; util.openInputFile(fastafile, in);
		int fastaAlignLength = 0;
        
		if (in) {
            if (!prot) {
                Sequence tempIn(in);
                fastaAlignLength = tempIn.getAligned().length();
            }else {
                Protein tempIn(in);
                fastaAlignLength = tempIn.getAligned().size();
            }
		}
		in.close();
		
		ifstream in2; util.openInputFile(oldfastafile, in2);
		int oldfastaAlignLength = 0;
		if (in2) {
            if (!prot) {
                Sequence tempIn(in2);
                oldfastaAlignLength = tempIn.getAligned().length();
            }else {
                Protein tempIn(in2);
                oldfastaAlignLength = tempIn.getAligned().size();
            }
		}
		in2.close();
		
		if (fastaAlignLength != oldfastaAlignLength) { m->mothurOut("fasta files do not have the same alignment length.\n");  return false;  }
		
        //read fasta file and save names as well as adding them to the alignDB
        set<string> namesOldFasta;
        
        ifstream inFasta; util.openInputFile(oldfastafile, inFasta);
        
        while (!inFasta.eof()) {
            if (m->getControl_pressed()) {  inFasta.close(); return good;  }
            
            if (!prot) {
                Sequence temp(inFasta);  util.gobble(inFasta);
                if (temp.getName() != "") {
                    namesOldFasta.insert(temp.getName());  //save name
                    if (!fitCalc) { db->push_back(temp);  }//add to DB
                }
            }else {
                Protein temp(inFasta);  util.gobble(inFasta);
                if (temp.getName() != "") {
                    namesOldFasta.insert(temp.getName());  //save name
                    if (!fitCalc) { db->push_back(temp);  }//add to DB
                }
            }
        }
        inFasta.close();
        
		//read through the column file checking names and removing distances above the cutoff
		ifstream inDist;
		util.openInputFile(column, inDist);
		
		ofstream outDist;
		string outputFile = column + ".temp";
		util.openOutputFile(outputFile, outDist);
		
		string name1, name2;
		float dist;
		while (!inDist.eof()) {
			if (m->getControl_pressed()) {  inDist.close(); outDist.close(); util.mothurRemove(outputFile); return good;  }
		
            inDist >> name1; util.gobble(inDist); inDist >> name2; util.gobble(inDist); inDist >> dist; util.gobble(inDist);
			
			//both names are in fasta file and distance is below cutoff
			if ((namesOldFasta.count(name1) == 0) || (namesOldFasta.count(name2) == 0)) {  good = false; break;  }
            else{ if (dist <= cutoff) { numDistsBelowCutoff++; outDist << name1 << '\t' << name2 << '\t' << dist << endl; } }
		}
		
		inDist.close();
		outDist.close();
		
		if (good) {
			util.mothurRemove(column);
			rename(outputFile.c_str(), column.c_str());
		}else{
			util.mothurRemove(outputFile); //temp file is bad because file mismatch above
		}
		
		return good;
		
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceCommand", "sanityCheck");
		exit(1);
	}
}
/**************************************************************************************************/




