/*
 *  chimeraperseuscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 10/26/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "chimeraperseuscommand.h"
#include "deconvolutecommand.h"
#include "sequence.hpp"
#include "counttable.h"
#include "sequencecountparser.h"
#include "removeseqscommand.h"

//**********************************************************************************************************************
vector<string> ChimeraPerseusCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","chimera-accnos",false,true,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "NameCount", "NameCount", "none","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "NameCount", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
        CommandParameter pdups("dereplicate", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pdups);
        CommandParameter premovechimeras("removechimeras", "Boolean", "", "t", "", "", "","fasta",false,false); parameters.push_back(premovechimeras);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		CommandParameter pcutoff("cutoff", "Number", "", "0.5", "", "", "","",false,false); parameters.push_back(pcutoff);
		CommandParameter palpha("alpha", "Number", "", "-5.54", "", "", "","",false,false); parameters.push_back(palpha);
		CommandParameter pbeta("beta", "Number", "", "0.33", "", "", "","",false,false); parameters.push_back(pbeta);
        
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["chimera"] = tempOutNames;
        outputTypes["accnos"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
        outputTypes["fasta"] = tempOutNames;
			
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraPerseusCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The chimera.perseus command reads a fastafile and namefile or countfile and outputs potentially chimeric sequences.\n";
		helpString += "The chimera.perseus command parameters are fasta, name, group, cutoff, processors, dereplicate, alpha and beta.\n";
		helpString += "The fasta parameter allows you to enter the fasta file containing your potentially chimeric sequences, and is required, unless you have a valid current fasta file. \n";
		helpString += "The name parameter allows you to provide a name file associated with your fasta file.\n";
        helpString += "The count parameter allows you to provide a count file associated with your fasta file. A count or name file is required. When you use a count file with group info and dereplicate=T, mothur will create a *.pick.count_table file containing seqeunces after chimeras are removed.\n";
		helpString += "The group parameter allows you to provide a group file.  When checking sequences, only sequences from the same group as the query sequence will be used as the reference. \n";
		helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n";
        helpString += "If the dereplicate parameter is false, then if one group finds the seqeunce to be chimeric, then all groups find it to be chimeric, default=f.\n";
        helpString += "The removechimeras parameter allows you to indicate you would like to automatically remove the sequences that are flagged as chimeric. Default=t.\n";
		helpString += "The alpha parameter ....  The default is -5.54. \n";
		helpString += "The beta parameter ....  The default is 0.33. \n";
		helpString += "The cutoff parameter ....  The default is 0.50. \n";
		helpString += "The chimera.perseus command should be in the following format: \n";
		helpString += "chimera.perseus(fasta=yourFastaFile, name=yourNameFile) \n";
		helpString += "Example: chimera.perseus(fasta=AD.align, name=AD.names) \n";
			
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraPerseusCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "chimera") {  pattern = "[filename],perseus.chimeras"; } 
        else if (type == "accnos") {  pattern = "[filename],perseus.accnos"; }
        else if (type == "fasta") {  pattern = "[filename],perseus.fasta"; }
        else if (type == "count") {  pattern = "[filename],perseus.count_table"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraPerseusCommand", "getOutputPattern");
        exit(1);
    }
}
//***************************************************************************************************************
ChimeraPerseusCommand::ChimeraPerseusCommand(string option) : Command()  {
	try {
        hasCount = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();

			//check for required parameters
            ValidParameters validParameter;
            fastafile = validParameter.validFile(parameters, "fasta");
            if (fastafile == "not found") {
                fastafile = current->getFastaFile();
                if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter.\n"); }
                else { 	m->mothurOut("[ERROR]: You have no current fasta file and the fasta parameter is required.\n");  abort = true; }
            }
            else if (fastafile == "not open") { abort = true; }
            else { current->setFastaFile(fastafile); }
			
            bool hasName = false;
            string namefile = validParameter.validFile(parameters, "name");
            if (namefile == "not open") { namefile = ""; abort = true; }
            else if (namefile == "not found") {  namefile = "";  }
            else { current->setNameFile(namefile); }
            if (namefile != "") { hasName = true; }
            
            countfile = validParameter.validFile(parameters, "count");
            if (countfile == "not open") { countfile = ""; abort = true; }
            else if (countfile == "not found") { countfile = "";  }
            else { current->setCountFile(countfile); }
            if (countfile != "") { hasCount = true; }
            
			//make sure there is at least one valid file left
            if (hasName && hasCount) { m->mothurOut("[ERROR]: You must enter ONLY ONE of the following: count or name.\n");  abort = true; }
            
            if (!hasName && !hasCount) { 
                //if there is a current name file, use it, else look for current count file
				string filename = current->getNameFile();
				if (filename != "") { hasName = true; namefile = filename; m->mothurOut("Using " + filename + " as input file for the name parameter.\n"); }
				else { 
                    filename = current->getCountFile();
                    if (filename != "") { hasCount = true; countfile = filename; m->mothurOut("Using " + filename + " as input file for the count parameter.\n"); }
                    else { m->mothurOut("[ERROR]: You must provide a count or name file.\n");  abort = true;  }
                }
            }
            
			bool hasGroup = false;
            string groupfile = validParameter.validFile(parameters, "group");
            if (groupfile == "not open") { abort = true; }
            else if (groupfile == "not found") {  groupfile = "";  }
            else { current->setGroupFile(groupfile); hasGroup = true; }
			
            if (hasGroup && hasCount) { m->mothurOut("[ERROR]: You must enter ONLY ONE of the following: count or group.\n");  abort = true; }
            
			string temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
			
			temp = validParameter.valid(parameters, "cutoff");	if (temp == "not found"){	temp = "0.50";	}
			util.mothurConvert(temp, cutoff);
			
			temp = validParameter.valid(parameters, "alpha");	if (temp == "not found"){	temp = "-5.54";	}
			util.mothurConvert(temp, alpha);
			
			temp = validParameter.valid(parameters, "beta");	if (temp == "not found"){	temp = "0.33";	}
			util.mothurConvert(temp, beta);
            
			temp = validParameter.valid(parameters, "dereplicate");	
			if (temp == "not found") { temp = "false";			}
			dups = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "removechimeras");            if (temp == "not found") { temp = "t"; }
            removeChimeras = util.isTrue(temp);
            
            if (!abort) {
                if ((namefile != "") || (groupfile != "")) { //convert to count
                    
                    string rootFileName = namefile;
                    if (rootFileName == "") { rootFileName = groupfile; }
                    
                    if (outputdir == "") { outputdir = util.hasPath(rootFileName); }
                    string outputFileName = outputdir + util.getRootName(util.getSimpleName(rootFileName)) + "count_table";
                    
                    CountTable ct; ct.createTable(namefile, groupfile, nullVector); ct.printCompressedTable(outputFileName);
                    outputNames.push_back(outputFileName); 
                    
                    current->setCountFile(outputFileName);
                    countfile = outputFileName;
                    hasCount = true;
                }
            }
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "ChimeraPerseusCommand");
		exit(1);
	}
}
/**************************************************************************************************/
struct perseusData {
    Utils util;
    MothurOut* m;
    vector<seqData> sequences;
    string group;
    int count, numChimeras;
    string chimeraFileName;
    string accnosFileName;
    double alpha, beta, cutoff;
    
    perseusData(string cf, string ac, double a, double b, double c){
        m = MothurOut::getInstance();
        count = 0;
        numChimeras = 0;
        accnosFileName = ac;
        chimeraFileName = cf;
        alpha = a;
        beta = b;
        cutoff = c;
    }
};
//**********************************************************************************************************************
//void driver(string chimeraFileName, vector<seqData>& sequences, string accnosFileName, int& numChimeras){
void driver(perseusData* params){
    try {
        vector<vector<double> > correctModel(4);	//could be an option in the future to input own model matrix
        for(int i=0;i<4;i++){	correctModel[i].resize(4);	}
        
        correctModel[0][0] = 0.000000;	//AA
        correctModel[1][0] = 11.619259;	//CA
        correctModel[2][0] = 11.694004;	//TA
        correctModel[3][0] = 7.748623;	//GA
        
        correctModel[1][1] = 0.000000;	//CC
        correctModel[2][1] = 7.619657;	//TC
        correctModel[3][1] = 12.852562;	//GC
        
        correctModel[2][2] = 0.000000;	//TT
        correctModel[3][2] = 10.964048;	//TG
        
        correctModel[3][3] = 0.000000;	//GG
        
        for(int i=0;i<4;i++){ for(int j=0;j<i;j++){ correctModel[j][i] = correctModel[i][j]; } }
        
        int numSeqs = params->sequences.size();
        int alignLength = params->sequences[0].sequence.size();
        
        ofstream chimeraFile;
        ofstream accnosFile;
        params->util.openOutputFile(params->chimeraFileName, chimeraFile);
        params->util.openOutputFile(params->accnosFileName, accnosFile);
        
        Perseus myPerseus;
        vector<vector<double> > binMatrix = myPerseus.binomial(alignLength);
        
        chimeraFile << "SequenceIndex\tName\tDiffsToBestMatch\tBestMatchIndex\tBestMatchName\tDiffstToChimera\tIndexofLeftParent\tIndexOfRightParent\tNameOfLeftParent\tNameOfRightParent\tDistanceToBestMatch\tcIndex\t(cIndex - singleDist)\tloonIndex\tMismatchesToChimera\tMismatchToTrimera\tChimeraBreakPoint\tLogisticProbability\tTypeOfSequence\n";
        
        vector<bool> chimeras(numSeqs, 0);
        
        for(int i=0;i<numSeqs;i++){
            if (params->m->getControl_pressed()) { chimeraFile.close(); accnosFile.close(); break; }
            
            vector<bool> restricted = chimeras;
            
            vector<vector<int> > leftDiffs(numSeqs);
            vector<vector<int> > leftMaps(numSeqs);
            vector<vector<int> > rightDiffs(numSeqs);
            vector<vector<int> > rightMaps(numSeqs);
            
            vector<int> singleLeft, bestLeft;
            vector<int> singleRight, bestRight;
            
            int bestSingleIndex, bestSingleDiff;
            vector<pwAlign> alignments(numSeqs);
            
            int comparisons = myPerseus.getAlignments(i, params->sequences, alignments, leftDiffs, leftMaps, rightDiffs, rightMaps, bestSingleIndex, bestSingleDiff, restricted);
            if (params->m->getControl_pressed()) { chimeraFile.close(); accnosFile.close(); break; }
            
            int minMismatchToChimera, leftParentBi, rightParentBi, breakPointBi;
            
            string dummyA, dummyB;
            
            if (params->sequences[i].sequence.size() < 3) {
                chimeraFile << i << '\t' << params->sequences[i].seqName << "\t0\t0\tNull\t0\t0\t0\tNull\tNull\t0.0\t0.0\t0.0\t0\t0\t0\t0.0\t0.0\tgood" << endl;
            }else if(comparisons >= 2){
                minMismatchToChimera = myPerseus.getChimera(params->sequences, leftDiffs, rightDiffs, leftParentBi, rightParentBi, breakPointBi, singleLeft, bestLeft, singleRight, bestRight, restricted);
                if (params->m->getControl_pressed()) { chimeraFile.close(); accnosFile.close(); break; }
                
                int minMismatchToTrimera = numeric_limits<int>::max();
                int leftParentTri, middleParentTri, rightParentTri, breakPointTriA, breakPointTriB;
                
                if(minMismatchToChimera >= 3 && comparisons >= 3){
                    minMismatchToTrimera = myPerseus.getTrimera(params->sequences, leftDiffs, leftParentTri, middleParentTri, rightParentTri, breakPointTriA, breakPointTriB, singleLeft, bestLeft, singleRight, bestRight, restricted);
                    if (params->m->getControl_pressed()) { chimeraFile.close(); accnosFile.close(); break; }
                }
                
                double singleDist = myPerseus.modeledPairwiseAlignSeqs(params->sequences[i].sequence, params->sequences[bestSingleIndex].sequence, dummyA, dummyB, correctModel);
                
                if (params->m->getControl_pressed()) { chimeraFile.close(); accnosFile.close(); break; }
                
                string type;
                string chimeraRefSeq;
                
                if(minMismatchToChimera - minMismatchToTrimera >= 3){
                    type = "trimera";
                    chimeraRefSeq = myPerseus.stitchTrimera(alignments, leftParentTri, middleParentTri, rightParentTri, breakPointTriA, breakPointTriB, leftMaps, rightMaps);
                }
                else{
                    type = "chimera";
                    chimeraRefSeq = myPerseus.stitchBimera(alignments, leftParentBi, rightParentBi, breakPointBi, leftMaps, rightMaps);
                }
                
                if (params->m->getControl_pressed()) { chimeraFile.close(); accnosFile.close(); break; }
                
                double chimeraDist = myPerseus.modeledPairwiseAlignSeqs(params->sequences[i].sequence, chimeraRefSeq, dummyA, dummyB, correctModel);
                
                if (params->m->getControl_pressed()) { chimeraFile.close(); accnosFile.close(); break; }
                
                double cIndex = chimeraDist;//modeledPairwiseAlignSeqs(sequences[i].sequence, chimeraRefSeq);
                double loonIndex = myPerseus.calcLoonIndex(params->sequences[i].sequence, params->sequences[leftParentBi].sequence, params->sequences[rightParentBi].sequence, breakPointBi, binMatrix);
                
                if (params->m->getControl_pressed()) { chimeraFile.close(); accnosFile.close(); break; }
                
                chimeraFile << i << '\t' << params->sequences[i].seqName << '\t' << bestSingleDiff << '\t' << bestSingleIndex << '\t' << params->sequences[bestSingleIndex].seqName << '\t';
                chimeraFile << minMismatchToChimera << '\t' << leftParentBi << '\t' << rightParentBi << '\t' << params->sequences[leftParentBi].seqName << '\t' << params->sequences[rightParentBi].seqName << '\t';
                chimeraFile << singleDist << '\t' << cIndex << '\t' << (cIndex - singleDist) << '\t' << loonIndex << '\t';
                chimeraFile << minMismatchToChimera << '\t' << minMismatchToTrimera << '\t' << breakPointBi << '\t';
                
                double probability = myPerseus.classifyChimera(singleDist, cIndex, loonIndex, params->alpha, params->beta);
                
                chimeraFile << probability << '\t';
                
                if(probability > params->cutoff){
                    chimeraFile << type << endl;
                    accnosFile << params->sequences[i].seqName << endl;
                    chimeras[i] = 1;
                    params->numChimeras++;
                }
                else{ chimeraFile << "good" << endl; }
            }
            else{
                chimeraFile << i << '\t' << params->sequences[i].seqName << "\t0\t0\tNull\t0\t0\t0\tNull\tNull\t0.0\t0.0\t0.0\t0\t0\t0\t0.0\t0.0\tgood" << endl;
            }
            
            //report progress
            if((i+1) % 100 == 0){ 	params->m->mothurOutJustToScreen("Processing sequence: " + toString(i+1) + "\n");		}
            params->count++; //# of sequences completed. Used by calling function to check for failure
        }
        
        if((numSeqs) % 100 != 0){ 	params->m->mothurOutJustToScreen("Processing sequence: " + toString(numSeqs) + "\n");		}
        
        if (!params->m->getControl_pressed()) { chimeraFile.close(); accnosFile.close(); }
    }
    catch(exception& e) {
        params->m->errorOut(e, "ChimeraPerseusCommand", "driver");
        exit(1);
    }
}
//***************************************************************************************************************

int ChimeraPerseusCommand::execute(){
	try{
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        m->mothurOut("Checking sequences from " + fastafile + " ...\n" );
        
        long start = time(NULL);
        if (outputdir == "") { outputdir = util.hasPath(fastafile);  }
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(fastafile));
        string outputFileName = getOutputFileName("chimera", variables);
        string accnosFileName = getOutputFileName("accnos", variables);
        string newCountFile = "";
        
        if (countfile == "") { countfile = getCountFile(fastafile); hasCount=true; }
        
        if (m->getControl_pressed()) {  return 0;	}
        
        int numSeqs = 0; int numChimeras = 0;
        
        if (hasCount) {
            CountTable ct;
            vector<string> groups;
            if (ct.testGroups(countfile, groups)) { //fills groups if count file has them
                
                variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(countfile));
                newCountFile = getOutputFileName("count", variables);
                
                vector<string> groups;
                map<string, vector<string> > group2Files;
                if (hasCount) {
                    current->setMothurCalling(true);
                    SequenceCountParser cparser(countfile, fastafile, nullVector);
                    current->setMothurCalling(false);
                    groups = cparser.getNamesOfGroups();
                    group2Files = cparser.getFiles();
                }
                
                if (m->getControl_pressed()) { return 0; }
                
                //clears files
                ofstream out, out1, out2;
                util.openOutputFile(outputFileName, out); out.close();
                util.openOutputFile(accnosFileName, out1); out1.close();
                
                string countlist = accnosFileName+".byCount";
                numSeqs = createProcessesGroups(group2Files, outputFileName, countlist, accnosFileName, newCountFile, groups, fastafile, countfile, numChimeras);
                
                if (m->getControl_pressed()) {  for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0;	}
                
                if (!dups) {
                    numChimeras = deconvoluteResults(outputFileName, accnosFileName);
                }else {
                    CountTable newCount; newCount.readTable(countfile, true, false);
                    
                    if (!util.isBlank(countlist)) {
                        ifstream in2; util.openInputFile(countlist, in2);
                        
                        string name, group;
                        while (!in2.eof()) {
                            in2 >> name; util.gobble(in2); in2 >> group; util.gobble(in2);
                            newCount.setAbund(name, group, 0);
                        }
                        in2.close();
                    }
                    util.mothurRemove(countlist);
                    
                    //print new *.pick.count_table
                    vector<string> namesInTable = newCount.printTable(newCountFile);  //returns non zeroed names
                    outputNames.push_back(newCountFile); outputTypes["count"].push_back(newCountFile);
                    
                    set<string> doNotRemove = util.mothurConvert(namesInTable);
                   
                    //remove names we want to keep from accnos file.
                    set<string> accnosNames = util.readAccnos(accnosFileName);
                    ofstream out2; util.openOutputFile(accnosFileName, out2);
                    for (set<string>::iterator it = accnosNames.begin(); it != accnosNames.end(); it++) { if (doNotRemove.count(*it) == 0) {  out2 << (*it) << endl; } }
                    out2.close();
                }
                
                util.mothurRemove(countlist);
                m->mothurOut("The number of sequences checked may be larger than the number of unique sequences because some sequences are found in several samples.\n");
                
                if (m->getControl_pressed()) {   for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0;  }
                
            }else {
                if (processors != 1) { m->mothurOut("Your count file does not contain group information, mothur can only use 1 processor, continuing.\n");  processors = 1; }
                
                //read sequences and store sorted by frequency
                ct.readTable(countfile, false, false);
                vector<seqData> sequences = readFiles(fastafile, ct.getNameMap());
                
                if (m->getControl_pressed()) {  for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	} return 0; }
                
                perseusData* dataBundle = new perseusData(outputFileName, accnosFileName, alpha, beta, cutoff);
                dataBundle->sequences = sequences;
                driver(dataBundle);
                numSeqs = dataBundle->count; numChimeras = dataBundle->numChimeras;
                delete dataBundle;
            }
        }
        
        if (m->getControl_pressed()) { for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	} return 0; }
        
        m->mothurOut("\nIt took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences. " + toString(numChimeras) + " chimeras were found.\n");
        outputNames.push_back(outputFileName); outputTypes["chimera"].push_back(outputFileName);
        outputNames.push_back(accnosFileName); outputTypes["accnos"].push_back(accnosFileName);
        
        if (removeChimeras) {
            if (!util.isBlank(accnosFileName)) {
                m->mothurOut("\nRemoving chimeras from your input files:\n");
                
                string inputString = "fasta=" + fastafile + ", accnos=" + accnosFileName;
                if ((countfile != "") && (!dups))   {   inputString += ", count=" + countfile;  }
                
                m->mothurOut("/******************************************/\n");
                m->mothurOut("Running command: remove.seqs(" + inputString + ")\n");
                current->setMothurCalling(true);
                
                Command* removeCommand = new RemoveSeqsCommand(inputString);
                removeCommand->execute();
                
                map<string, vector<string> > filenames = removeCommand->getOutputFiles();
                
                delete removeCommand;
                current->setMothurCalling(false);
                m->mothurOut("/******************************************/\n");

                if (countfile != "") {
                    if (!dups) { //dereplicate=f, so remove sequences where any sample found the reads to be chimeric
                        map<string, string> variables;
                        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(countfile));
                        string currentName = getOutputFileName("count", variables);

                        util.renameFile(filenames["count"][0], currentName);
                        util.mothurRemove(filenames["count"][0]);
                        outputNames.push_back(currentName); outputTypes["count"].push_back(currentName);
                    }//else, mothur created a modified count file removing chimeras by sample. No need to include count file on remove.seqs command. Deconvolute function created modified count table already
                }
                
                map<string, string> variables;
                variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(fastafile));
                string currentName = getOutputFileName("fasta", variables);

                util.renameFile(filenames["fasta"][0], currentName);
                util.mothurRemove(filenames["fasta"][0]);
                
                outputNames.push_back(currentName); outputTypes["fasta"].push_back(currentName);
            }else { m->mothurOut("\nNo chimeras found, skipping remove.seqs.\n"); }
        }
        
		//set accnos file as new current accnosfile
		string currentName = "";
		itTypes = outputTypes.find("accnos");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setAccnosFile(currentName); }
		}
        
        itTypes = outputTypes.find("count");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
		}
		
        itTypes = outputTypes.find("fasta");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
        }
        
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}	
		m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraPerseusCommand::getCountFile(string& inputFile){
	try {
		string countFile = "";
		
		m->mothurOut("\nNo count file given, running unique.seqs command to generate one.\n\n");
		
		//use unique.seqs to create new name and fastafile
		string inputString = "format=count, fasta=" + inputFile;
		m->mothurOut("/******************************************/\n");
		m->mothurOut("Running command: unique.seqs(" + inputString + ")\n");
		current->setMothurCalling(true);
        
		Command* uniqueCommand = new DeconvoluteCommand(inputString);
		uniqueCommand->execute();
		
		map<string, vector<string> > filenames = uniqueCommand->getOutputFiles();
		
		delete uniqueCommand;
		current->setMothurCalling(false);
		m->mothurOut("/******************************************/\n");
		
		countFile = filenames["count"][0];
		inputFile = filenames["fasta"][0];
		
		return countFile;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "getNamesFile");
		exit(1);
	}
}
/**************************************************************************************************/
struct perseusGroupsData {
    string fastafile;
    string dupsfile;
    string chimeraFileName;
    string accnosFileName;
    string countlist;
    map<string, vector<string> > parsedFiles;
    
    bool hasCount, dups;
    int threadID, count, numChimeras;
    double alpha, beta, cutoff;
    vector<string> groups;
    Utils util;
    MothurOut* m;
    
    perseusGroupsData(){}
    perseusGroupsData(map<string, vector<string> >& g2f,bool dps, bool hc, double a, double b, double c, string o,  string f, string n, string ac, string ctlist, vector<string> gr, int tid) {
        alpha = a;
        beta = b;
        cutoff = c;
        fastafile = f;
        dupsfile = n;
        chimeraFileName = o;
        countlist = ctlist;
        accnosFileName = ac;
        m = MothurOut::getInstance();
        threadID = tid;
        groups = gr;
        hasCount = hc;
        dups = dps;
        count = 0;
        numChimeras = 0;
        parsedFiles = g2f;
    }
};
//**********************************************************************************************************************
vector<seqData> loadSequences(map<string, int>& nameMap, string thisGroupsFastaFile, perseusGroupsData* params){
    try {
        bool error = false;
        vector<seqData> sequences;
        
        ifstream in;
        params->util.openInputFile(thisGroupsFastaFile, in);
        
        vector<seqPriorityNode> nameVector;
        map<string, int>::iterator itNameMap;
        while (!in.eof()) {
            if (params->m->getControl_pressed()) { break; }
            
            Sequence seq(in); params->util.gobble(in);
            
            itNameMap = nameMap.find(seq.getName());
            
            if (itNameMap == nameMap.end()){
                error = true;
                params->m->mothurOut("[ERROR]: " + seq.getName() + " is in your fastafile, but is not in your name or count file, please correct.\n");
            }else {
                int num = itNameMap->second;
                
                seq.setAligned(params->util.removeNs(seq.getUnaligned()));
                sequences.push_back(seqData(seq.getName(), seq.getUnaligned(), num));
            }
        }
        in.close();

        if (error) { params->m->setControl_pressed(true); }
        
        //sort by frequency
        sort(sequences.rbegin(), sequences.rend());
        
        return sequences;
    }
    catch(exception& e) {
        params->m->errorOut(e, "ChimeraPerseusCommand", "loadSequences");
        exit(1);
    }
}

//**********************************************************************************************************************
//string outputFName, string accnos, string countlist, int start, int end, vector<string> groups
void driverGroups(perseusGroupsData* params){
	try {
        //clears files
        ofstream out, out1, out2;
        params->util.openOutputFile(params->chimeraFileName, out); out.close();
        params->util.openOutputFile(params->accnosFileName, out1); out1.close();
        
		int totalSeqs = 0;
        ofstream outCountList;
        if (params->hasCount && params->dups) { params->util.openOutputFile(params->countlist, outCountList); }
		
        for (map<string, vector<string> >::iterator it = params->parsedFiles.begin(); it != params->parsedFiles.end(); it++) {
            long start = time(NULL);	 if (params->m->getControl_pressed()) {  break; }
            
            
            string thisGroup = it->first;
            
            map<string, int> nameMap;
            if (params->hasCount) {
                CountTable ct; ct.readTable(it->second[1], false, true);
                nameMap = ct.getNameMap();
            }

			params->m->mothurOut("\nChecking sequences from group " + thisGroup + "...\n");
			
            perseusData* driverParams = new perseusData((params->chimeraFileName+thisGroup), (params->accnosFileName+thisGroup), params->alpha, params->beta, params->cutoff);
			driverParams->sequences = loadSequences(nameMap, it->second[0], params);
			
            if (params->m->getControl_pressed()) { break; }
			
			driver(driverParams);
			totalSeqs += driverParams->count;
            params->numChimeras += driverParams->numChimeras;
			
			if (params->m->getControl_pressed()) { break; }
            
            if (params->dups) {
                if (!params->util.isBlank(driverParams->accnosFileName)) {
                    ifstream in;
                    params->util.openInputFile(driverParams->accnosFileName, in);
                    string name;
                    if (params->hasCount) {
                        while (!in.eof()) {
                            in >> name; params->util.gobble(in);
                            outCountList << name << '\t' << thisGroup << endl;
                        }
                        in.close();
                    }
                }
            }
			
			//append files
			params->util.appendFiles(driverParams->chimeraFileName, params->chimeraFileName); params->util.mothurRemove(driverParams->chimeraFileName);
			params->util.appendFiles(driverParams->accnosFileName, params->accnosFileName); params->util.mothurRemove(driverParams->accnosFileName);
			
			params->m->mothurOut("\nIt took " + toString(time(NULL) - start) + " secs to check " + toString(driverParams->count) + " sequences from group " + thisGroup + ".\n");
            delete driverParams;
		}	
		
        if (params->hasCount && params->dups) { outCountList.close(); }
        
		params->count = totalSeqs;
		
	}
	catch(exception& e) {
		params->m->errorOut(e, "ChimeraPerseusCommand", "driverGroups");
		exit(1);
	}
}	
//**********************************************************************************************************************
vector<seqData> ChimeraPerseusCommand::readFiles(string inputFile, map<string, int> nameMap){
	try {
		map<string, int>::iterator it;
		
		//read fasta file and create sequenceData structure - checking for file mismatches
		vector<seqData> sequences;
		bool error = false;
		ifstream in;
		util.openInputFile(inputFile, in);
		alignLength = 0;
        
		while (!in.eof()) {
			
			if (m->getControl_pressed()) { in.close(); return sequences; }
			
			Sequence temp(in); util.gobble(in);
			
			it = nameMap.find(temp.getName());
			if (it == nameMap.end()) { error = true; m->mothurOut("[ERROR]: " + temp.getName() + " is in your fasta file and not in your namefile, please correct.\n");  }
			else {
                temp.setAligned(util.removeNs(temp.getUnaligned()));
				sequences.push_back(seqData(temp.getName(), temp.getUnaligned(), it->second));
                if (temp.getUnaligned().length() > alignLength) { alignLength = temp.getUnaligned().length(); }
			}
		}
		in.close();
		
		if (error) { m->setControl_pressed(true); }
		
		//sort by frequency
		sort(sequences.rbegin(), sequences.rend());
		
		return sequences;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "readFiles");
		exit(1);
	}
}
/**************************************************************************************************/
//perseusData(vector<seqData>& s, double a, double b, double c, string o, string ac, MothurOut* mout)
//numSeqs = createProcessesGroups(outputFileName, countlist, accnosFileName, newCountFile, groups, fastafile, countfile, numChimeras);
int ChimeraPerseusCommand::createProcessesGroups(map<string, vector<string> >& parsedFiles, string outputFName, string countlisttemp, string accnos, string newCountFile, vector<string> groups, string fasta, string dupsFile, int& numChimeras) {
	try {
        numChimeras = 0;
        
		//sanity check
		if (groups.size() < processors) { processors = groups.size(); m->mothurOut("Reducing processors to " + toString(groups.size()) + ".\n"); }
		
		//divide the groups between the processors
		vector<linePair> lines;
		int remainingPairs = groups.size();
        int startIndex = 0;
        for (int remainingProcessors = processors; remainingProcessors > 0; remainingProcessors--) {
            int numPairs = remainingPairs; //case for last processor
            if (remainingProcessors != 1) { numPairs = ceil(remainingPairs / remainingProcessors); }
            lines.push_back(linePair(startIndex, (startIndex+numPairs))); //startIndex, endIndex
            startIndex = startIndex + numPairs;
            remainingPairs = remainingPairs - numPairs;
        }

        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<perseusGroupsData*> data;
        
        long long num = 0;
        time_t start, end;
        time(&start);
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            string extension = toString(i+1) + ".temp";
            vector<string> thisGroups;
            map<string, vector<string> > thisGroupsParsedFiles;
            for (int j = lines[i+1].start; j < lines[i+1].end; j++) {
                
                map<string, vector<string> >::iterator it = parsedFiles.find(groups[j]);
                if (it != parsedFiles.end()) {
                    thisGroupsParsedFiles[groups[j]] = (it->second);
                    thisGroups.push_back(groups[j]);
                }
                else { m->mothurOut("[ERROR]: missing files for group " + groups[j] + ", skipping\n"); }
            }
            perseusGroupsData* dataBundle = new perseusGroupsData(thisGroupsParsedFiles, dups, hasCount, alpha, beta, cutoff, (outputFName+extension), fasta, dupsFile,  (accnos+extension), (countlisttemp+extension), thisGroups, (i+1));
            data.push_back(dataBundle);
            
            workerThreads.push_back(new std::thread(driverGroups, dataBundle));
        }
        
        vector<string> thisGroups;
        map<string, vector<string> > thisGroupsParsedFiles;
        for (int j = lines[0].start; j < lines[0].end; j++) {
            
            map<string, vector<string> >::iterator it = parsedFiles.find(groups[j]);
            if (it != parsedFiles.end()) {
                thisGroupsParsedFiles[groups[j]] = (it->second);
                thisGroups.push_back(groups[j]);
            }
            else { m->mothurOut("[ERROR]: missing files for group " + groups[j] + ", skipping\n"); }
        }
        perseusGroupsData* dataBundle = new perseusGroupsData(thisGroupsParsedFiles, dups, hasCount, alpha, beta, cutoff, outputFName, fasta, dupsFile,  accnos, countlisttemp, thisGroups, 0);
        driverGroups(dataBundle);
        num = dataBundle->count;
        numChimeras = dataBundle->numChimeras;

        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            num += data[i]->count;
            numChimeras += data[i]->numChimeras;
            
            string extension = toString(i+1) + ".temp";
            util.appendFiles((outputFName+extension), outputFName);
            util.mothurRemove((outputFName+extension));
            
            util.appendFiles((accnos+extension), accnos);
            util.mothurRemove((accnos+extension));
            
            util.appendFiles((countlisttemp+extension), countlisttemp);
            util.mothurRemove((countlisttemp+extension));
            
            delete data[i];
            delete workerThreads[i];
        }
        delete dataBundle;
        
        time(&end);
        m->mothurOut("It took " + toString(difftime(end, start)) + " secs to check " + toString(num) + " sequences.\n\n");

		return num;	
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "createProcessesGroups");
		exit(1);
	}
}
//**********************************************************************************************************************
int ChimeraPerseusCommand::deconvoluteResults(string outputFileName, string accnosFileName){
	try {
        int total = 0;
		
        set<string> chimerasInFile = util.readAccnos(accnosFileName);//this is so if a sequence is found to be chimera in several samples we dont write it to the results file more than once
        util.printAccnos(accnosFileName, chimerasInFile);
    
		//edit chimera file
		ifstream in; 
		util.openInputFile(outputFileName, in);
		
		ofstream out;
		util.openOutputFile(outputFileName+".temp", out); out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		int DiffsToBestMatch, BestMatchIndex, DiffstToChimera, IndexofLeftParent, IndexOfRightParent;
		float temp1,temp2, temp3, temp4, temp5, temp6, temp7, temp8;
		string index, BestMatchName, parent1, parent2, flag;
		string name = "";
        set<string> namesInFile;
			
		//assumptions - in file each read will always look like 
		/*										
		 SequenceIndex	Name	DiffsToBestMatch	BestMatchIndex	BestMatchName	DiffstToChimera	IndexofLeftParent	IndexOfRightParent	NameOfLeftParent	NameOfRightParent	DistanceToBestMatch	cIndex	(cIndex - singleDist)	loonIndex	MismatchesToChimera	MismatchToTrimera	ChimeraBreakPoint	LogisticProbability	TypeOfSequence
		 0	F01QG4L02JVBQY	0	0	Null	0	0	0	Null	Null	0.0	0.0	0.0	0.0	0	0	0	0.0	0.0	good
		 1	F01QG4L02ICTC6	0	0	Null	0	0	0	Null	Null	0.0	0.0	0.0	0.0	0	0	0	0.0	0.0	good
		 2	F01QG4L02JZOEC	48	0	F01QG4L02JVBQY	47	0	0	F01QG4L02JVBQY	F01QG4L02JVBQY	2.0449	2.03545	-0.00944493	0	47	2147483647	138	0	good
		 3	F01QG4L02G7JEC	42	0	F01QG4L02JVBQY	40	1	0	F01QG4L02ICTC6	F01QG4L02JVBQY	1.87477	1.81113	-0.0636404	5.80145	40	2147483647	25	0	good
		 */
		
		//get and print headers
		BestMatchName = util.getline(in); util.gobble(in);
		out << BestMatchName << endl;
		
		while (!in.eof()) {
			
			if (m->getControl_pressed()) { in.close(); out.close(); util.mothurRemove((outputFileName+".temp")); return 0; }
			
			bool print = false;
			in >> index;	util.gobble(in);
			
			if (index != "SequenceIndex") { //if you are not a header line, there will be a header line for each group if group file is given
				in >> name;		util.gobble(in);
				in >> DiffsToBestMatch; util.gobble(in);
				in >> BestMatchIndex; util.gobble(in);
				in >> BestMatchName; util.gobble(in);
				in >> DiffstToChimera; util.gobble(in);
				in >> IndexofLeftParent; util.gobble(in);
				in >> IndexOfRightParent; util.gobble(in);
				in >> parent1;	util.gobble(in);
				in >> parent2;	util.gobble(in);
				in >> temp1 >> temp2 >> temp3 >> temp4 >> temp5 >> temp6 >> temp7 >> temp8 >> flag; util.gobble(in);
				
				
                //is this name already in the file
                set<string>::iterator itNames = namesInFile.find((name));
                
                if (itNames == namesInFile.end()) { //no not in file
                    if (flag == "good") { //are you really a no??
                        //is this sequence really not chimeric??
                        set<string>::iterator itChimeras = chimerasInFile.find(name);
                        
                        //then you really are a no so print, otherwise skip
                        if (itChimeras == chimerasInFile.end()) { print = true; }
                    }else{ print = true; }
                }
                
				if (print) {
                    namesInFile.insert(name);
					out << index << '\t' << name  << '\t' << DiffsToBestMatch << '\t' << BestMatchIndex << '\t';
					out << BestMatchName << '\t' << DiffstToChimera << '\t' << IndexofLeftParent << '\t' << IndexOfRightParent << '\t' << parent1 << '\t' << parent2 << '\t';
					out << temp1 << '\t' << temp2 << '\t' << temp3 << '\t' << temp4 << '\t' << temp5 << '\t' << temp6 << '\t' << temp7 << '\t' << temp8 << '\t' << flag << endl;	
				}
			}else { index = util.getline(in); util.gobble(in); }
		}
		in.close();
		out.close();
		
		util.mothurRemove(outputFileName);
		rename((outputFileName+".temp").c_str(), outputFileName.c_str());
		
		return total;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "deconvoluteResults");
		exit(1);
	}
}	
//**********************************************************************************************************************


