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
//**********************************************************************************************************************
vector<string> ChimeraPerseusCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","chimera-accnos",false,true,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "NameCount", "NameCount", "none","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "NameCount", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
        CommandParameter pdups("dereplicate", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pdups);

		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		CommandParameter pcutoff("cutoff", "Number", "", "0.5", "", "", "","",false,false); parameters.push_back(pcutoff);
		CommandParameter palpha("alpha", "Number", "", "-5.54", "", "", "","",false,false); parameters.push_back(palpha);
		CommandParameter pbeta("beta", "Number", "", "0.33", "", "", "","",false,false); parameters.push_back(pbeta);
			
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
		helpString += "You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amazon.fasta \n";
		helpString += "The group parameter allows you to provide a group file.  When checking sequences, only sequences from the same group as the query sequence will be used as the reference. \n";
		helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n";
        helpString += "If the dereplicate parameter is false, then if one group finds the seqeunce to be chimeric, then all groups find it to be chimeric, default=f.\n";
		helpString += "The alpha parameter ....  The default is -5.54. \n";
		helpString += "The beta parameter ....  The default is 0.33. \n";
		helpString += "The cutoff parameter ....  The default is 0.50. \n";
		helpString += "The chimera.perseus command should be in the following format: \n";
		helpString += "chimera.perseus(fasta=yourFastaFile, name=yourNameFile) \n";
		helpString += "Example: chimera.perseus(fasta=AD.align, name=AD.names) \n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n";	
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
        else if (type == "count") {  pattern = "[filename],perseus.pick.count_table"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraPerseusCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ChimeraPerseusCommand::ChimeraPerseusCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
		vector<string> tempOutNames;
		outputTypes["chimera"] = tempOutNames;
		outputTypes["accnos"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "ChimeraPerseusCommand");
		exit(1);
	}
}
//***************************************************************************************************************
ChimeraPerseusCommand::ChimeraPerseusCommand(string option)  {
	try {
		abort = false; calledHelp = false; 
        hasCount = false;
        hasName = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter("chimera.perseus");
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			vector<string> tempOutNames;
			outputTypes["chimera"] = tempOutNames;
			outputTypes["accnos"] = tempOutNames;
            outputTypes["count"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			
			//check for required parameters
			fastafile = validParameter.valid(parameters, "fasta");
			if (fastafile == "not found") { 				
				//if there is a current fasta file, use it
				string filename = current->getFastaFile(); 
				if (filename != "") { fastaFileNames.push_back(filename); m->mothurOut("Using " + filename + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { 
				util.splitAtDash(fastafile, fastaFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastaFileNames.size(); i++) {
					
					bool ignore = false;
					if (fastaFileNames[i] == "current") { 
						fastaFileNames[i] = current->getFastaFile(); 
						if (fastaFileNames[i] != "") {  m->mothurOut("Using " + fastaFileNames[i] + " as input file for the fasta parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current fastafile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							fastaFileNames.erase(fastaFileNames.begin()+i);
							i--;
						}
					}
					
                    if (!ignore) {
                        if (util.checkLocations(fastaFileNames[i], current->getLocations())) { current->setFastaFile(fastaFileNames[i]); }
                        else { fastaFileNames.erase(fastaFileNames.begin()+i); i--; } //erase from file list
                    }
                }
				
				//make sure there is at least one valid file left
				if (fastaFileNames.size() == 0) { m->mothurOut("[ERROR]: no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			
			//check for required parameters
			namefile = validParameter.valid(parameters, "name");
			if (namefile == "not found") { namefile = "";  	}
			else { 
				util.splitAtDash(namefile, nameFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < nameFileNames.size(); i++) {
					
					bool ignore = false;
					if (nameFileNames[i] == "current") { 
						nameFileNames[i] = current->getNameFile();
						if (nameFileNames[i] != "") {  m->mothurOut("Using " + nameFileNames[i] + " as input file for the name parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current namefile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							nameFileNames.erase(nameFileNames.begin()+i);
							i--;
						}
					}
					
                    if (!ignore) {
                        if (util.checkLocations(nameFileNames[i], current->getLocations())) { current->setNameFile(nameFileNames[i]); }
                        else { nameFileNames.erase(nameFileNames.begin()+i); i--; } //erase from file list
                    }
                }
			}
            
            if (nameFileNames.size() != 0) { hasName = true; }
            
            //check for required parameters
            vector<string> countfileNames;
			countfile = validParameter.valid(parameters, "count");
			if (countfile == "not found") { 
                countfile = "";  
			}else { 
				util.splitAtDash(countfile, countfileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < countfileNames.size(); i++) {
					
					bool ignore = false;
					if (countfileNames[i] == "current") { 
						countfileNames[i] = current->getCountFile();
						if (countfileNames[i] != "") {  m->mothurOut("Using " + countfileNames[i] + " as input file for the count parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current count file, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							countfileNames.erase(countfileNames.begin()+i);
							i--;
						}
					}
					
                    if (!ignore) {
                        if (util.checkLocations(countfileNames[i], current->getLocations())) { current->setCountFile(countfileNames[i]); }
                        else { countfileNames.erase(countfileNames.begin()+i); i--; } //erase from file list
                    }
				}
			}
            
            if (countfileNames.size() != 0) { hasCount = true; }
            
			//make sure there is at least one valid file left
            if (hasName && hasCount) { m->mothurOut("[ERROR]: You must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
            
            if (!hasName && !hasCount) { 
                //if there is a current name file, use it, else look for current count file
				string filename = current->getNameFile();
				if (filename != "") { hasName = true; nameFileNames.push_back(filename); m->mothurOut("Using " + filename + " as input file for the name parameter."); m->mothurOutEndLine(); }
				else { 
                    filename = current->getCountFile();
                    if (filename != "") { hasCount = true; countfileNames.push_back(filename); m->mothurOut("Using " + filename + " as input file for the count parameter."); m->mothurOutEndLine(); }
                    else { m->mothurOut("[ERROR]: You must provide a count or name file."); m->mothurOutEndLine(); abort = true;  }
                }
            }
            if (!hasName && hasCount) { nameFileNames = countfileNames; }
            
			if (nameFileNames.size() != fastaFileNames.size()) { m->mothurOut("[ERROR]: The number of name or count files does not match the number of fastafiles, please correct."); m->mothurOutEndLine(); abort=true; }
			
			bool hasGroup = true;
			groupfile = validParameter.valid(parameters, "group");
			if (groupfile == "not found") { groupfile = "";  hasGroup = false; }
			else { 
				util.splitAtDash(groupfile, groupFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < groupFileNames.size(); i++) {
					
					bool ignore = false;
					if (groupFileNames[i] == "current") { 
						groupFileNames[i] = current->getGroupFile();
						if (groupFileNames[i] != "") {  m->mothurOut("Using " + groupFileNames[i] + " as input file for the group parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current namefile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							groupFileNames.erase(groupFileNames.begin()+i);
							i--;
						}
					}
					
                    if (!ignore) {
                        if (util.checkLocations(groupFileNames[i], current->getLocations())) { current->setGroupFile(groupFileNames[i]); }
                        else { groupFileNames.erase(groupFileNames.begin()+i); i--; } //erase from file list
                    }
                }
				
				//make sure there is at least one valid file left
				if (groupFileNames.size() == 0) { m->mothurOut("[ERROR]: no valid group files."); m->mothurOutEndLine(); abort = true; }
			}
			
			if (hasGroup && (groupFileNames.size() != fastaFileNames.size())) { m->mothurOut("[ERROR]: The number of groupfiles does not match the number of fastafiles, please correct."); m->mothurOutEndLine(); abort=true; }
			
            if (hasGroup && hasCount) { m->mothurOut("[ERROR]: You must enter ONLY ONE of the following: count or group."); m->mothurOutEndLine(); abort = true; }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";	}
			
			string temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
			
			temp = validParameter.valid(parameters, "cutoff");	if (temp == "not found"){	temp = "0.50";	}
			util.mothurConvert(temp, cutoff);
			
			temp = validParameter.valid(parameters, "alpha");	if (temp == "not found"){	temp = "-5.54";	}
			util.mothurConvert(temp, alpha);
			
			temp = validParameter.valid(parameters, "cutoff");	if (temp == "not found"){	temp = "0.33";	}
			util.mothurConvert(temp, beta);
            
			temp = validParameter.validFile(parameters, "dereplicate");	
			if (temp == "not found") { temp = "false";			}
			dups = util.isTrue(temp);
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
		
				
		//process each file
		for (int s = 0; s < fastaFileNames.size(); s++) {
			
			m->mothurOut("Checking sequences from " + fastaFileNames[s] + " ..." ); m->mothurOutEndLine();
			
			long start = time(NULL);	
			if (outputDir == "") { outputDir = util.hasPath(fastaFileNames[s]);  }//if user entered a file with a path then preserve it	
			map<string, string> variables;
			variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastaFileNames[s]));
			string outputFileName = getOutputFileName("chimera", variables);
			string accnosFileName = getOutputFileName("accnos", variables);
            string newCountFile = "";

			//string newFasta = util.getRootName(fastaFileNames[s]) + "temp";
			
			//you provided a groupfile
			string groupFile = "";
			if (groupFileNames.size() != 0) { groupFile = groupFileNames[s]; }
            
			string nameFile = "";
			if (nameFileNames.size() != 0) { //you provided a namefile and we don't need to create one
				nameFile = nameFileNames[s];
			}else { nameFile = getNamesFile(fastaFileNames[s]); }
			
			if (m->getControl_pressed()) {  for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	} return 0;	}			
			
			int numSeqs = 0;
			int numChimeras = 0;
            
            if (hasCount) {
                CountTable* ct = new CountTable();
                ct->readTable(nameFile, true, false);
                
                if (ct->hasGroupInfo()) {
                    vector<string> temp;
                    SequenceCountParser* cparser = new SequenceCountParser(fastaFileNames[s], *ct, temp);
                    variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(nameFile));
                    newCountFile = getOutputFileName("count", variables);
                    
                    vector<string> groups = cparser->getNamesOfGroups();
                    
                    if (m->getControl_pressed()) { delete ct; delete cparser; for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0; }
                    
                    //clears files
                    ofstream out, out1, out2;
                    util.openOutputFile(outputFileName, out); out.close(); 
                    util.openOutputFile(accnosFileName, out1); out1.close();
                    
                    numSeqs = createProcessesGroups(outputFileName, accnosFileName, newCountFile, groups, groupFile, fastaFileNames[s], nameFile, numChimeras);
                    
                    if (m->getControl_pressed()) {  delete ct; delete cparser; for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0;	}				
                    map<string, string> uniqueNames = cparser->getAllSeqsMap();
                    if (!dups) { 
                        numChimeras = deconvoluteResults(uniqueNames, outputFileName, accnosFileName);
                    }else {
                        set<string> doNotRemove;
                        CountTable c; c.readTable(newCountFile, true, true);
                        vector<string> namesInTable = c.getNamesOfSeqs();
                        for (int i = 0; i < namesInTable.size(); i++) {
                            int temp = c.getNumSeqs(namesInTable[i]);
                            if (temp == 0) {  c.remove(namesInTable[i]);  }
                            else { doNotRemove.insert((namesInTable[i])); }
                        }
                        //remove names we want to keep from accnos file.
                        set<string> accnosNames = util.readAccnos(accnosFileName);
                        ofstream out2;
                        util.openOutputFile(accnosFileName, out2);
                        for (set<string>::iterator it = accnosNames.begin(); it != accnosNames.end(); it++) {
                            if (doNotRemove.count(*it) == 0) {  out2 << (*it) << endl; }
                        }
                        out2.close();
                        c.printTable(newCountFile);
                        outputNames.push_back(newCountFile); outputTypes["count"].push_back(newCountFile);

                    }
                    delete cparser;

                    m->mothurOut("The number of sequences checked may be larger than the number of unique sequences because some sequences are found in several samples."); m->mothurOutEndLine(); 
                    
                    if (m->getControl_pressed()) {  delete ct; for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0;  }	
                    
                }else {
                    if (processors != 1) { m->mothurOut("Your count file does not contain group information, mothur can only use 1 processor, continuing."); m->mothurOutEndLine(); processors = 1; }
                    
                    //read sequences and store sorted by frequency
                    vector<seqData> sequences = readFiles(fastaFileNames[s], ct);
                    
                    if (m->getControl_pressed()) { delete ct; for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	} return 0; }
                    
                    perseusData* dataBundle = new perseusData(outputFileName, accnosFileName, alpha, beta, cutoff);
                    dataBundle->sequences = sequences;
                    driver(dataBundle);
                    numSeqs = dataBundle->count; numChimeras = dataBundle->numChimeras;
                    delete dataBundle;
                }
                delete ct;
            }else {
                if (groupFile != "") {
                    //Parse sequences by group
                    vector<string> temp;
                    SequenceParser* parser = new SequenceParser(groupFile, fastaFileNames[s], nameFile, temp);
                    vector<string> groups = parser->getNamesOfGroups();
                    
                    if (m->getControl_pressed()) { delete parser; for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0; }
                    
                    //clears files
                    ofstream out, out1, out2;
                    util.openOutputFile(outputFileName, out); out.close(); 
                    util.openOutputFile(accnosFileName, out1); out1.close();
                    
                    numSeqs = createProcessesGroups(outputFileName, accnosFileName, "", groups, groupFile, fastaFileNames[s], nameFile, numChimeras);
                    
                    if (m->getControl_pressed()) {  delete parser; for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0;	}				
                    map<string, string> uniqueNames = parser->getAllSeqsMap();
                    if (!dups) { 
                        numChimeras = deconvoluteResults(uniqueNames, outputFileName, accnosFileName);
                    }
                    delete parser;
                    
                    m->mothurOut("The number of sequences checked may be larger than the number of unique sequences because some sequences are found in several samples."); m->mothurOutEndLine(); 
                    
                    if (m->getControl_pressed()) {  for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0;  }		
                }else{
                    if (processors != 1) { m->mothurOut("Without a groupfile, mothur can only use 1 processor, continuing."); m->mothurOutEndLine(); processors = 1; }
                    
                    //read sequences and store sorted by frequency
                    vector<seqData> sequences = readFiles(fastaFileNames[s], nameFile);
                    
                    if (m->getControl_pressed()) { for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	} return 0; }
                    
                    perseusData* dataBundle = new perseusData(outputFileName, accnosFileName, alpha, beta, cutoff);
                    dataBundle->sequences = sequences;
                    driver(dataBundle);
                    numSeqs = dataBundle->count; numChimeras = dataBundle->numChimeras;
                    delete dataBundle;
                }
			}
            
			if (m->getControl_pressed()) { for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	} return 0; }
			
			m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences. " + toString(numChimeras) + " chimeras were found.");	m->mothurOutEndLine();
			outputNames.push_back(outputFileName); outputTypes["chimera"].push_back(outputFileName);
			outputNames.push_back(accnosFileName); outputTypes["accnos"].push_back(accnosFileName);
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
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
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
string ChimeraPerseusCommand::getNamesFile(string& inputFile){
	try {
		string nameFile = "";
		
		m->mothurOutEndLine(); m->mothurOut("No namesfile given, running unique.seqs command to generate one."); m->mothurOutEndLine(); m->mothurOutEndLine();
		
		//use unique.seqs to create new name and fastafile
		string inputString = "fasta=" + inputFile;
		m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
		m->mothurOut("Running command: unique.seqs(" + inputString + ")"); m->mothurOutEndLine(); 
		current->setMothurCalling(true);
        
		Command* uniqueCommand = new DeconvoluteCommand(inputString);
		uniqueCommand->execute();
		
		map<string, vector<string> > filenames = uniqueCommand->getOutputFiles();
		
		delete uniqueCommand;
		current->setMothurCalling(false);
		m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
		
		nameFile = filenames["name"][0];
		inputFile = filenames["fasta"][0];
		
		return nameFile;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "getNamesFile");
		exit(1);
	}
}
/**************************************************************************************************/
struct parserData {
    string group;
    SequenceParser* parser;
    SequenceCountParser* cparser;
    bool hasCount;
    Utils util;
    MothurOut* m;
    int alignLength;
    
    
    parserData(){}
    parserData(SequenceParser* p, SequenceCountParser* cp, bool hc) {
        parser = p;
        cparser = cp;
        hasCount = hc;
        group = "";
        
        m = MothurOut::getInstance();
    }
};
//**********************************************************************************************************************
vector<seqData> loadSequences(parserData* params){
    try {
        bool error = false;
        params->alignLength = 0;
        vector<seqData> sequences;
        if (params->hasCount) {
            vector<Sequence> thisGroupsSeqs = params->cparser->getSeqs(params->group);
            map<string, int> counts = params->cparser->getCountTable(params->group);
            map<string, int>::iterator it;
            
            for (int i = 0; i < thisGroupsSeqs.size(); i++) {
                
                if (params->m->getControl_pressed()) {  return sequences; }
                
                it = counts.find(thisGroupsSeqs[i].getName());
                if (it == counts.end()) { error = true; params->m->mothurOut("[ERROR]: " + thisGroupsSeqs[i].getName() + " is in your fasta file and not in your count file, please correct.\n");  }
                else {
                    thisGroupsSeqs[i].setAligned(params->util.removeNs(thisGroupsSeqs[i].getUnaligned()));
                    sequences.push_back(seqData(thisGroupsSeqs[i].getName(), thisGroupsSeqs[i].getUnaligned(), it->second));
                    if (thisGroupsSeqs[i].getUnaligned().length() > params->alignLength) { params->alignLength = thisGroupsSeqs[i].getUnaligned().length(); }
                }
            }
        }else{
            vector<Sequence> thisGroupsSeqs = params->parser->getSeqs(params->group);
            map<string, string> nameMap = params->parser->getNameMap(params->group);
            map<string, string>::iterator it;
            
            for (int i = 0; i < thisGroupsSeqs.size(); i++) {
                
                if (params->m->getControl_pressed()) {  return sequences; }
                
                it = nameMap.find(thisGroupsSeqs[i].getName());
                if (it == nameMap.end()) { error = true; params->m->mothurOut("[ERROR]: " + thisGroupsSeqs[i].getName() + " is in your fasta file and not in your namefile, please correct.\n");  }
                else {
                    int num = params->util.getNumNames(it->second);
                    thisGroupsSeqs[i].setAligned(params->util.removeNs(thisGroupsSeqs[i].getUnaligned()));
                    sequences.push_back(seqData(thisGroupsSeqs[i].getName(), thisGroupsSeqs[i].getUnaligned(), num));
                    if (thisGroupsSeqs[i].getUnaligned().length() > params->alignLength) { params->alignLength = thisGroupsSeqs[i].getUnaligned().length(); }
                }
            }
            
        }
        
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
/**************************************************************************************************/
struct perseusGroupsData {
    string fastafile;
    string namefile;
    string groupfile;
    string chimeraFileName;
    string accnosFileName;
    string countlist;
    
    bool hasName, hasCount, dups;
    int threadID, count, numChimeras;
    double alpha, beta, cutoff;
    vector<string> groups;
    Utils util;
    MothurOut* m;
    
    perseusGroupsData(){}
    perseusGroupsData(bool dps, bool hn, bool hc, double a, double b, double c, string o,  string f, string n, string g, string ac, string ctlist, vector<string> gr, int tid) {
        alpha = a;
        beta = b;
        cutoff = c;
        fastafile = f;
        namefile = n;
        groupfile = g;
        chimeraFileName = o;
        countlist = ctlist;
        accnosFileName = ac;
        m = MothurOut::getInstance();
        threadID = tid;
        groups = gr;
        hasName = hn;
        hasCount = hc;
        dups = dps;
        count = 0;
        numChimeras = 0;
    }
};
//**********************************************************************************************************************
//string outputFName, string accnos, string countlist, int start, int end, vector<string> groups
void driverGroups(perseusGroupsData* params){
	try {
        //clears files
        ofstream out, out1, out2;
        params->util.openOutputFile(params->chimeraFileName, out); out.close();
        params->util.openOutputFile(params->accnosFileName, out1); out1.close();
        
        //parse fasta and name file by group
        SequenceParser* parser;
        SequenceCountParser* cparser;
        if (params->hasCount) {
            CountTable* ct = new CountTable();
            ct->readTable(params->namefile, true, false);
            cparser = new SequenceCountParser(params->fastafile, *ct, params->groups);
            delete ct;
        }else {
            if (params->namefile != "") { parser = new SequenceParser(params->groupfile, params->fastafile, params->namefile, params->groups); }
            else						{ parser = new SequenceParser(params->groupfile, params->fastafile, params->groups);	 }
        }
        parserData* sequenceLoader = new parserData(parser, cparser, params->hasCount);
        
		int totalSeqs = 0;
        ofstream outCountList;
        if (params->hasCount && params->dups) { params->util.openOutputFile(params->countlist, outCountList); }
		
		for (int i = 0; i < params->groups.size(); i++) {
			
			params->m->mothurOut("\nChecking sequences from group " + params->groups[i] + "...\n");
			
            long start = time(NULL);	 if (params->m->getControl_pressed()) {  break; }
			
            sequenceLoader->group = params->groups[i];
            perseusData* driverParams = new perseusData((params->chimeraFileName+params->groups[i]), (params->accnosFileName+params->groups[i]), params->alpha, params->beta, params->cutoff);
			driverParams->sequences = loadSequences(sequenceLoader);
			
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
                            outCountList << name << '\t' << params->groups[i] << endl;
                        }
                        in.close();
                    }else {
                        map<string, string> thisnamemap = parser->getNameMap(params->groups[i]);
                        map<string, string>::iterator itN;
                        ofstream out;
                        params->util.openOutputFile(driverParams->accnosFileName+".temp", out);
                        while (!in.eof()) {
                            in >> name; params->util.gobble(in);
                            itN = thisnamemap.find(name);
                            if (itN != thisnamemap.end()) {
                                vector<string> tempNames; params->util.splitAtComma(itN->second, tempNames);
                                for (int j = 0; j < tempNames.size(); j++) { out << tempNames[j] << endl; }
                                
                            }else { params->m->mothurOut("[ERROR]: parsing cannot find " + name + ".\n"); params->m->setControl_pressed(true); }
                        }
                        out.close();
                        in.close();
                        params->util.renameFile(driverParams->accnosFileName+".temp", driverParams->accnosFileName);
                    }
                    
                }
            }
			
			//append files
			params->util.appendFiles(driverParams->chimeraFileName, params->chimeraFileName); params->util.mothurRemove(driverParams->chimeraFileName);
			params->util.appendFiles(driverParams->accnosFileName, params->accnosFileName); params->util.mothurRemove(driverParams->accnosFileName);
			
			params->m->mothurOut("\nIt took " + toString(time(NULL) - start) + " secs to check " + toString(driverParams->count) + " sequences from group " + params->groups[i] + ".\n");
            delete driverParams;
		}	
		
        delete sequenceLoader;
        if (params->hasCount && params->dups) { outCountList.close(); }
        
		params->count = totalSeqs;
		
	}
	catch(exception& e) {
		params->m->errorOut(e, "ChimeraPerseusCommand", "driverGroups");
		exit(1);
	}
}	
//**********************************************************************************************************************
vector<seqData> ChimeraPerseusCommand::readFiles(string inputFile, string name){
	try {
		map<string, int>::iterator it;
		map<string, int> nameMap = util.readNames(name);
		
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
			if (it == nameMap.end()) { error = true; m->mothurOut("[ERROR]: " + temp.getName() + " is in your fasta file and not in your namefile, please correct."); m->mothurOutEndLine(); }
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
//**********************************************************************************************************************
vector<seqData> ChimeraPerseusCommand::readFiles(string inputFile, CountTable* ct){
	try {		
		//read fasta file and create sequenceData structure - checking for file mismatches
		vector<seqData> sequences;
		ifstream in;
		util.openInputFile(inputFile, in);
		alignLength = 0;
        
		while (!in.eof()) {
            Sequence temp(in); util.gobble(in);
			
			int count = ct->getNumSeqs(temp.getName());
			if (m->getControl_pressed()) { break; }
			else {
                temp.setAligned(util.removeNs(temp.getUnaligned()));
				sequences.push_back(seqData(temp.getName(), temp.getUnaligned(), count));
                if (temp.getUnaligned().length() > alignLength) { alignLength = temp.getUnaligned().length(); }
			}
		}
		in.close();
		
		//sort by frequency
		sort(sequences.rbegin(), sequences.rend());
		
		return sequences;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "getNamesFile");
		exit(1);
	}
}
/**************************************************************************************************/
//perseusData(vector<seqData>& s, double a, double b, double c, string o, string ac, MothurOut* mout)
int ChimeraPerseusCommand::createProcessesGroups(string outputFName, string accnos, string newCountFile, vector<string> groups, string group, string fasta, string name, int& numChimeras) {
	try {
        numChimeras = 0;
        CountTable newCount;
        if (hasCount && dups) { newCount.readTable(name, true, false); }
        
		//sanity check
		if (groups.size() < processors) { processors = groups.size(); }
		
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
        vector<thread*> workerThreads;
        vector<perseusGroupsData*> data;
        
        long long num = 0;
        time_t start, end;
        time(&start);
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            string extension = toString(i+1) + ".temp";
            vector<string> thisGroups;
            for (int j = lines[i+1].start; j < lines[i+1].end; j++) { thisGroups.push_back(groups[j]); }
            perseusGroupsData* dataBundle = new perseusGroupsData(dups, hasName, hasCount, alpha, beta, cutoff, (outputFName+extension), fasta, name, group, (accnos+extension), (accnos+".byCount."+extension), thisGroups, (i+1));
            data.push_back(dataBundle);
            
            workerThreads.push_back(new thread(driverGroups, dataBundle));
        }
        
        vector<string> thisGroups;
        for (int j = lines[0].start; j < lines[0].end; j++) { thisGroups.push_back(groups[j]); }
        perseusGroupsData* dataBundle = new perseusGroupsData(dups, hasName, hasCount, alpha, beta, cutoff, outputFName, fasta, name, group, accnos, (accnos+".byCount.temp"), thisGroups, 0);
        driverGroups(dataBundle);
        num = dataBundle->count;
        numChimeras = dataBundle->numChimeras;
        
        delete dataBundle;
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            num += data[i]->count;
            numChimeras += data[i]->numChimeras;
            
            delete data[i];
            delete workerThreads[i];
        }
        
        time(&end);
        m->mothurOut("It took " + toString(difftime(end, start)) + " secs to check " + toString(num) + " sequences.\n\n");

		//read my own
        if (hasCount && dups) {
            string countlisttemp = accnos+".byCount.temp";
            if (!util.isBlank(countlisttemp)) {
                ifstream in2;
                util.openInputFile(countlisttemp, in2);
                
                string name, group;
                while (!in2.eof()) {
                    in2 >> name >> group; util.gobble(in2);
                    newCount.setAbund(name, group, 0);
                }
                in2.close();
            }
            util.mothurRemove(countlisttemp);
        }
		
		//append output files
		 for (int i = 0; i < processors-1; i++) {
            string extension = toString(i+1) + ".temp";
			util.appendFiles((outputFName+extension), outputFName);
			util.mothurRemove((outputFName+extension));
			
			util.appendFiles((accnos+extension), accnos);
			util.mothurRemove((accnos+extension));
            
            if (hasCount && dups) {
                if (!util.isBlank(accnos+".byCount."+extension)) {
                    ifstream in2;
                    util.openInputFile(accnos+".byCount."+extension, in2);
                    
                    string name, group;
                    while (!in2.eof()) {
                        in2 >> name >> group; util.gobble(in2);
                        newCount.setAbund(name, group, 0);
                    }
                    in2.close();
                }
                util.mothurRemove(accnos+".byCount."+extension);
            }

		}
		
        //print new *.pick.count_table
        if (hasCount && dups) {  newCount.printTable(newCountFile);   }

		return num;	
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPerseusCommand", "createProcessesGroups");
		exit(1);
	}
}
//**********************************************************************************************************************
int ChimeraPerseusCommand::deconvoluteResults(map<string, string>& uniqueNames, string outputFileName, string accnosFileName){
	try {
		map<string, string>::iterator itUnique;
		int total = 0;
		
		//edit accnos file
		ifstream in2; 
		util.openInputFile(accnosFileName, in2);
		
		ofstream out2;
		util.openOutputFile(accnosFileName+".temp", out2);
		
		string name;
		set<string> namesInFile; //this is so if a sequence is found to be chimera in several samples we dont write it to the results file more than once
		set<string>::iterator itNames;
		set<string> chimerasInFile;
		set<string>::iterator itChimeras;
		
		
		while (!in2.eof()) {
			if (m->getControl_pressed()) { in2.close(); out2.close(); util.mothurRemove(outputFileName); util.mothurRemove((accnosFileName+".temp")); return 0; }
			
			in2 >> name; util.gobble(in2);
			
			//find unique name
			itUnique = uniqueNames.find(name);
			
			if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing accnos results. Cannot find "+ name + "."); m->mothurOutEndLine(); m->setControl_pressed(true); }
			else {
				itChimeras = chimerasInFile.find((itUnique->second));
				
				if (itChimeras == chimerasInFile.end()) {
					out2 << itUnique->second << endl;
					chimerasInFile.insert((itUnique->second));
					total++;
				}
			}
		}
		in2.close();
		out2.close();
		
		util.mothurRemove(accnosFileName);
		rename((accnosFileName+".temp").c_str(), accnosFileName.c_str());
		
		//edit chimera file
		ifstream in; 
		util.openInputFile(outputFileName, in);
		
		ofstream out;
		util.openOutputFile(outputFileName+".temp", out); out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		int DiffsToBestMatch, BestMatchIndex, DiffstToChimera, IndexofLeftParent, IndexOfRightParent;
		float temp1,temp2, temp3, temp4, temp5, temp6, temp7, temp8;
		string index, BestMatchName, parent1, parent2, flag;
		name = "";
		namesInFile.clear();	
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
				
				//find unique name
				itUnique = uniqueNames.find(name);
				
				if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find "+ name + "."); m->mothurOutEndLine(); m->setControl_pressed(true); }
				else {
					name = itUnique->second;
					//is this name already in the file
					itNames = namesInFile.find((name));
					
					if (itNames == namesInFile.end()) { //no not in file
						if (flag == "good") { //are you really a no??
							//is this sequence really not chimeric??
							itChimeras = chimerasInFile.find(name);
							
							//then you really are a no so print, otherwise skip
							if (itChimeras == chimerasInFile.end()) { print = true; }
						}else{ print = true; }
					}
				}
				
				if (print) {
					out << index << '\t' << name  << '\t' << DiffsToBestMatch << '\t' << BestMatchIndex << '\t';
					namesInFile.insert(name);
					
					if (BestMatchName != "Null") {
						itUnique = uniqueNames.find(BestMatchName);
						if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find BestMatchName "+ BestMatchName + "."); m->mothurOutEndLine(); m->setControl_pressed(true); }
						else {	out << itUnique->second << '\t';	}					
					}else { out << "Null" << '\t'; }
					
					out << DiffstToChimera << '\t' << IndexofLeftParent << '\t' << IndexOfRightParent << '\t';
					
					if (parent1 != "Null") {
						itUnique = uniqueNames.find(parent1);
						if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find parent1 "+ parent1 + "."); m->mothurOutEndLine(); m->setControl_pressed(true); }
						else {	out << itUnique->second << '\t';	}
					}else { out << "Null" << '\t'; }
					
					if (parent1 != "Null") {
						itUnique = uniqueNames.find(parent2);
						if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find parent2 "+ parent2 + "."); m->mothurOutEndLine(); m->setControl_pressed(true); }
						else {	out << itUnique->second << '\t';	}
					}else { out << "Null" << '\t'; }
					
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


