//
//  chimeravsearchcommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 6/16/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#include "chimeravsearchcommand.h"
#include "deconvolutecommand.h"
#include "sequence.hpp"
#include "systemcommand.h"
#include "degapseqscommand.h"

//**********************************************************************************************************************
vector<string> ChimeraVsearchCommand::setParameters(){
    try {
        CommandParameter ptemplate("reference", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(ptemplate);
        CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","chimera-accnos",false,true,true); parameters.push_back(pfasta);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","",false,false,true); parameters.push_back(pcount);
        CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
        CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        CommandParameter pabskew("abskew", "Number", "", "1.9", "", "", "","",false,false); parameters.push_back(pabskew);
        CommandParameter pchimealns("uchimealns", "Boolean", "", "F", "", "", "","alns",false,false); parameters.push_back(pchimealns);
        CommandParameter pminh("minh", "Number", "", "0.28", "", "", "","",false,false); parameters.push_back(pminh);
        CommandParameter pmindiv("mindiv", "Number", "", "0.8", "", "", "","",false,false); parameters.push_back(pmindiv);
        CommandParameter pxn("xn", "Number", "", "8.0", "", "", "","",false,false); parameters.push_back(pxn);
        CommandParameter pdn("dn", "Number", "", "1.4", "", "", "","",false,false); parameters.push_back(pdn);
        CommandParameter pmindiffs("mindiffs", "Number", "", "3", "", "", "","",false,false); parameters.push_back(pmindiffs);
        CommandParameter pdups("dereplicate", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pdups);
        
        vector<string> myArray;
        for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
        return myArray;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraVsearchCommand", "setParameters");
        exit(1);
    }
}
//**********************************************************************************************************************
string ChimeraVsearchCommand::getHelpString(){
    try {
        string helpString = "";
        helpString += "The chimera.vsearch command reads a fastafile and referencefile and outputs potentially chimeric sequences.\n";
        helpString += "This command is a wrapper for vsearch https://github.com/torognes/vsearch.\n";
        helpString += "The chimera.vsearch command parameters are fasta, name, count, reference, processors, dereplicate, abskew, uchimealns, minh, mindiv, xn, dn, mindiffs.\n";
        helpString += "The fasta parameter allows you to enter the fasta file containing your potentially chimeric sequences, and is required, unless you have a valid current fasta file. \n";
        helpString += "The name parameter allows you to provide a name file, if you are using template=self. \n";
        helpString += "The count parameter allows you to provide a count file, if you are using template=self. When you use a count file with group info and dereplicate=T, mothur will create a *.pick.count_table file containing seqeunces after chimeras are removed. \n";
        helpString += "You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amazon.fasta \n";
        helpString += "The group parameter allows you to provide a group file. The group file can be used with a namesfile and reference=self. When checking sequences, only sequences from the same group as the query sequence will be used as the reference. \n";
        helpString += "If the dereplicate parameter is false, then if one group finds the sequence to be chimeric, then all groups find it to be chimeric, default=f.\n";
        helpString += "The reference parameter allows you to enter a reference file containing known non-chimeric sequences, and is required. You may also set template=self, in this case the abundant sequences will be used as potential parents. \n";
        helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n";
        helpString += "The abskew parameter can only be used with template=self. Minimum abundance skew. Default 1.9. Abundance skew is: min [ abund(parent1), abund(parent2) ] / abund(query).\n";
        helpString += "The uchimealns parameter allows you to indicate you would like a file containing multiple alignments of query sequences to parents in human readable format. Alignments show columns with differences that support or contradict a chimeric model.\n";
        helpString += "The minh parameter - mininum score to report chimera. Default 0.3. Values from 0.1 to 5 might be reasonable. Lower values increase sensitivity but may report more false positives. If you decrease xn you may need to increase minh, and vice versa.\n";
        helpString += "The mindiv parameter - minimum divergence ratio, default 0.5. Div ratio is 100%% - %%identity between query sequence and the closest candidate for being a parent. If you don't care about very close chimeras, then you could increase mindiv to, say, 1.0 or 2.0, and also decrease minh, say to 0.1, to increase sensitivity. How well this works will depend on your data. Best is to tune parameters on a good benchmark.\n";
        helpString += "The xn parameter - weight of a no vote. Default 8.0. Decreasing this weight to around 3 or 4 may give better performance on denoised data.\n";
        helpString += "The dn parameter - pseudo-count prior on number of no votes. Default 1.4. Probably no good reason to change this unless you can retune to a good benchmark for your data. Reasonable values are probably in the range from 0.2 to 2.\n";
        helpString += "The mindiffs parameter - minimum number of differences in segment Default = (3).\n";
        helpString += "The chimera.vsearch command should be in the following format: \n";
        helpString += "chimera.vsearch(fasta=yourFastaFile, reference=yourTemplate) \n";
        helpString += "Example: chimera.vsearch(fasta=AD.align, reference=silva.gold.align) \n";
        
        return helpString;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraVsearchCommand", "getHelpString");
        exit(1);
    }
}
//**********************************************************************************************************************
string ChimeraVsearchCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "chimera") {  pattern = "[filename],[tag],vsearch.chimeras"; }
        else if (type == "accnos") {  pattern = "[filename],[tag],vsearch.accnos"; }
        else if (type == "alns") {  pattern = "[filename],[tag],vsearch.alns"; }
        else if (type == "count") {  pattern = "[filename],[tag],vsearch.pick.count_table"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraVsearchCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ChimeraVsearchCommand::ChimeraVsearchCommand() : Command(){
    try {
        abort = true; calledHelp = true;
        setParameters();
        vector<string> tempOutNames;
        outputTypes["chimera"] = tempOutNames;
        outputTypes["accnos"] = tempOutNames;
        outputTypes["alns"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraVsearchCommand", "ChimeraVsearchCommand");
        exit(1);
    }
}
//***************************************************************************************************************
ChimeraVsearchCommand::ChimeraVsearchCommand(string option) : Command() {
    try {
        abort = false; calledHelp = false; hasName=false; hasCount=false;
        
        //allow user to run help
        if(option == "help") { help(); abort = true; calledHelp = true; }
        else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        
        else {
            vector<string> myArray = setParameters();
            
            OptionParser parser(option);
            map<string,string> parameters = parser.getParameters();
            
            ValidParameters validParameter("chimera.vsearch");
            map<string,string>::iterator it;
            
            //check to make sure all parameters are valid for command
            for (it = parameters.begin(); it != parameters.end(); it++) {
                if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
            }
            
            vector<string> tempOutNames;
            outputTypes["chimera"] = tempOutNames;
            outputTypes["accnos"] = tempOutNames;
            outputTypes["alns"] = tempOutNames;
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
                        if (nameFileNames[i] != "") {  m->mothurOut("Using " + countfileNames[i] + " as input file for the count parameter where you had given current."); m->mothurOutEndLine(); }
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
            
            if (!hasName && hasCount) { nameFileNames = countfileNames; }
            
            if ((hasCount || hasName) && (nameFileNames.size() != fastaFileNames.size())) { m->mothurOut("[ERROR]: The number of name or count files does not match the number of fastafiles, please correct."); m->mothurOutEndLine(); abort=true; }
            
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
            
            
            //if the user changes the output directory command factory will send this info to us in the output parameter
            outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";	}
            
            string path;
            it = parameters.find("reference");
            //user has given a template file
            if(it != parameters.end()){
                if (it->second == "self") { templatefile = "self"; }
                else {
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["reference"] = inputDir + it->second;		}
                    
                    templatefile = validParameter.validFile(parameters, "reference");
                    if (templatefile == "not open") { abort = true; }
                    else if (templatefile == "not found") { //check for saved reference sequences
                        m->mothurOut("[ERROR]: The reference parameter is a required.\n"); abort = true;
                    }
                }
            }else if (hasName) {  templatefile = "self"; }
            else if (hasCount) {  templatefile = "self"; }
            else {
                m->mothurOut("[ERROR]: The reference parameter is a required.");
                
                templatefile = ""; abort = true;
            }
            
            string temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
            processors = current->setProcessors(temp);
            
            abskew = validParameter.valid(parameters, "abskew");	if (abskew == "not found"){	useAbskew = false;  abskew = "1.9";	}else{  useAbskew = true;  }
            if (useAbskew && templatefile != "self") { m->mothurOut("The abskew parameter is only valid with template=self, ignoring."); m->mothurOutEndLine(); useAbskew = false; }
            
            temp = validParameter.valid(parameters, "chimealns");			if (temp == "not found") { temp = "f"; }
            chimealns = util.isTrue(temp);
            
            minh = validParameter.valid(parameters, "minh");						if (minh == "not found")			{ useMinH = false; minh = "0.28";					}	else{ useMinH = true;			}
            mindiv = validParameter.valid(parameters, "mindiv");					if (mindiv == "not found")			{ useMindiv = false; mindiv = "0.8";				}	else{ useMindiv = true;			}
            xn = validParameter.valid(parameters, "xn");							if (xn == "not found")				{ useXn = false; xn = "8.0";						}	else{ useXn = true;				}
            dn = validParameter.valid(parameters, "dn");							if (dn == "not found")				{ useDn = false; dn = "1.4";						}	else{ useDn = true;				}
            mindiffs = validParameter.valid(parameters, "mindiffs");				if (mindiffs == "not found")				{ useMindiffs = false; mindiffs = "3";							}	else{ useMindiffs = true;				}
            
            temp = validParameter.valid(parameters, "dereplicate");
            if (temp == "not found") { temp = "false";			}
            dups = util.isTrue(temp);
            
            
            if (hasName && (templatefile != "self")) { m->mothurOut("You have provided a namefile and the reference parameter is not set to self. I am not sure what reference you are trying to use, aborting."); m->mothurOutEndLine(); abort=true; }
            if (hasCount && (templatefile != "self")) { m->mothurOut("You have provided a countfile and the reference parameter is not set to self. I am not sure what reference you are trying to use, aborting."); m->mothurOutEndLine(); abort=true; }
            if (hasGroup && (templatefile != "self")) { m->mothurOut("You have provided a group file and the reference parameter is not set to self. I am not sure what reference you are trying to use, aborting."); m->mothurOutEndLine(); abort=true; }

            //look for uchime exe
            path = current->getProgramPath();

            string vsearchCommand;
#if defined NON_WINDOWS
            vsearchCommand = path + "vsearch";	//	format the database, -o option gives us the ability
            if (m->getDebug()) {
                m->mothurOut("[DEBUG]: vsearch location using \"which vsearch\" = ");
                Command* newCommand = new SystemCommand("which vsearch\n");
                newCommand->execute();
                delete newCommand;
                m->mothurOut("[DEBUG]: Mothur's location using \"which mothur\" = ");
                newCommand = new SystemCommand("which mothur\n");
                newCommand->execute();
                delete newCommand;
            }
#else
            vsearchCommand = path + "\\vsearch.exe";
#endif
            
            //test to make sure uchime exists
            ifstream in;
            vsearchCommand = util.getFullPathName(vsearchCommand);
            bool ableToOpen = util.openInputFile(vsearchCommand, in, "no error"); in.close();
            if(!ableToOpen) {
                m->mothurOut(vsearchCommand + " file does not exist. Checking path... \n");
                //check to see if uchime is in the path??
                
                ifstream in2;
                string programName = "vsearch"; programName += EXECUTABLE_EXT;
                string uLocation = util.findProgramPath(programName);
                uLocation += programName;
                ableToOpen = util.openInputFile(uLocation, in2, "no error"); in2.close();
                
                if(!ableToOpen) { m->mothurOut("[ERROR]: " + uLocation + " file does not exist. mothur requires the vsearch executable.\n");  abort = true; }
                else {  m->mothurOut("Found vsearch in your path, using " + uLocation + "\n");vsearchLocation = uLocation; }
            }else {  vsearchLocation = vsearchCommand; }
            
            vsearchLocation = util.getFullPathName(vsearchLocation);
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: vsearch location using " + vsearchLocation + "\n"); }
        }
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraVsearchCommand", "ChimeraVsearchCommand");
        exit(1);
    }
}
//***************************************************************************************************************

int ChimeraVsearchCommand::execute(){
    try{
        
        if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        for (int s = 0; s < fastaFileNames.size(); s++) {
            
            m->mothurOut("Checking sequences from " + fastaFileNames[s] + " ..." ); m->mothurOutEndLine();
            
            long start = time(NULL);
            string nameFile = "";
            if (outputDir == "") { outputDir = util.hasPath(fastaFileNames[s]);  }//if user entered a file with a path then preserve it
            map<string, string> variables;
            variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastaFileNames[s]));
            variables["[tag]"] = "denovo";
            if (templatefile != "self") { variables["[tag]"] = "ref"; }
            string outputFileName = getOutputFileName("chimera", variables);
            string accnosFileName = getOutputFileName("accnos", variables);
            string alnsFileName = getOutputFileName("alns", variables);
            string newFasta = util.getRootName(fastaFileNames[s]) + "temp";
            string newCountFile = "";
            
            //you provided a groupfile
            string groupFile = "";
            bool hasGroup = false;
            int numSeqs = 0;
            if (groupFileNames.size() != 0) { groupFile = groupFileNames[s]; hasGroup = true; }
            else if (hasCount) {
                CountTable ct;
                if (ct.testGroups(nameFileNames[s])) { hasGroup = true; }
                variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(nameFileNames[s]));
                newCountFile = getOutputFileName("count", variables);
            }
            
            //setup fasta file if denovo and no groups
            if ((templatefile == "self") && (!hasGroup)) { //you want to run uchime with a template=self and no groups
                
                if (processors != 1) { m->mothurOut("When using template=self, mothur can only use 1 processor, continuing."); m->mothurOutEndLine(); processors = 1; }
                if (nameFileNames.size() != 0) { //you provided a namefile and we don't need to create one
                    nameFile = nameFileNames[s];
                }else { nameFile = getNamesFile(fastaFileNames[s]); }
                
                map<string, string> seqs;
                numSeqs = readFasta(fastaFileNames[s], seqs);  if (m->getControl_pressed()) { for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0; }
                
                //read namefile
                vector<seqPriorityNode> nameMapCount;
                int error;
                if (hasCount) {
                    CountTable ct; ct.readTable(nameFile, true, false);
                    for(map<string, string>::iterator it = seqs.begin(); it != seqs.end(); it++) {
                        int num = ct.getNumSeqs(it->first);
                        if (num == 0) { error = 1; }
                        else { seqPriorityNode temp(num, it->second, it->first); nameMapCount.push_back(temp); }
                    }
                }else { error = util.readNames(nameFile, nameMapCount, seqs); if (m->getControl_pressed() || (error == 1)) { for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}return 0; }
                }
                if (seqs.size() != nameMapCount.size()) { m->mothurOut( "The number of sequences in your fastafile does not match the number of sequences in your namefile, aborting."); m->mothurOutEndLine(); for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0; }
                
                util.printVsearchFile(nameMapCount, newFasta, ";size=", ";");
                fastaFileNames[s] = newFasta;
            }
            
            if (m->getControl_pressed()) {  for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0;	}
            
            if (hasGroup) {
                if (nameFileNames.size() != 0) { //you provided a namefile and we don't need to create one
                    nameFile = nameFileNames[s];
                }else { nameFile = getNamesFile(fastaFileNames[s]); }
                
                //Parse sequences by group
                vector<string> groups;
                map<string, string> uniqueNames;
                vector<string> temp;
                if (hasCount) {
                    cparser = new SequenceCountParser(nameFile, fastaFileNames[s], temp);
                    groups = cparser->getNamesOfGroups();
                    uniqueNames = cparser->getAllSeqsMap();
                }else{
                    sparser = new SequenceParser(groupFile, fastaFileNames[s], nameFile, temp);
                    groups = sparser->getNamesOfGroups();
                    uniqueNames = sparser->getAllSeqsMap();
                }
                
                if (m->getControl_pressed()) { for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0; }
                
                //clears files
                ofstream out, out1, out2;
                util.openOutputFile(outputFileName, out); out.close();
                util.openOutputFile(accnosFileName, out1); out1.close();
                if (chimealns) { util.openOutputFile(alnsFileName, out2); out2.close(); }
                
                //paralellized in vsearch
                driverGroups(outputFileName, newFasta, accnosFileName, alnsFileName, newCountFile, 0, groups.size(), groups);
                if (hasCount && dups) {
                    CountTable c; c.readTable(nameFile, true, false);
                    if (!util.isBlank(newCountFile)) {
                        ifstream in2;
                        util.openInputFile(newCountFile, in2);
                        
                        string name, group;
                        while (!in2.eof()) {
                            in2 >> name >> group; util.gobble(in2);
                            c.setAbund(name, group, 0);
                        }
                        in2.close();
                    }
                    util.mothurRemove(newCountFile);
                    c.printTable(newCountFile);
                }
                
                if (m->getControl_pressed()) {  for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0;	}
                
                if (!dups) {
                    int totalChimeras = deconvoluteResults(uniqueNames, outputFileName, accnosFileName, alnsFileName);
                    
                    m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(uniqueNames.size()) + " sequences. " + toString(totalChimeras) + " chimeras were found.");	m->mothurOutEndLine();
                    m->mothurOut("The number of sequences checked may be larger than the number of unique sequences because some sequences are found in several samples."); m->mothurOutEndLine();
                }else {
                    
                    if (hasCount) {
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
                }
                
                if (hasCount) { delete cparser; }
                else { delete sparser; }
                
                if (m->getControl_pressed()) {  for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0;	}
                
            }else{
                if (m->getControl_pressed()) {  for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0;	}
                
                int numChimeras = 0;
                
                //paralellized in vsearch
                driver(outputFileName, fastaFileNames[s], accnosFileName, alnsFileName, numChimeras);
                
                //add headings
                ofstream out;
                util.openOutputFile(outputFileName+".temp", out);
                out << "Score\tQuery\tParentA\tParentB\tIdQM\tIdQA\tIdQB\tIdAB\tIdQT\tLY\tLN\tLA\tRY\tRN\tRA\tDiv\tYN\n";
                out.close();
                
                util.appendFiles(outputFileName, outputFileName+".temp");
                util.mothurRemove(outputFileName); rename((outputFileName+".temp").c_str(), outputFileName.c_str());
                
                if (m->getControl_pressed()) { for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	} return 0; }
                
                //remove file made for uchime
                if (templatefile == "self") {  util.mothurRemove(fastaFileNames[s]); }
                
                m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check your sequences. " + toString(numChimeras) + " chimeras were found.");	m->mothurOutEndLine();
            }
            
            outputNames.push_back(outputFileName); outputTypes["chimera"].push_back(outputFileName);
            outputNames.push_back(accnosFileName); outputTypes["accnos"].push_back(accnosFileName);
            if (chimealns) { outputNames.push_back(alnsFileName); outputTypes["alns"].push_back(alnsFileName); }
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
        m->errorOut(e, "ChimeraVsearchCommand", "execute");
        exit(1);
    }
}
//**********************************************************************************************************************
int ChimeraVsearchCommand::deconvoluteResults(map<string, string>& uniqueNames, string outputFileName, string accnosFileName, string alnsFileName){
    try {
        map<string, string>::iterator itUnique;
        int total = 0;
        
        ofstream out2;
        util.openOutputFile(accnosFileName+".temp", out2);
        
        string name;
        set<string> namesInFile; //this is so if a sequence is found to be chimera in several samples we dont write it to the results file more than once
        set<string>::iterator itNames;
        set<string> chimerasInFile;
        set<string>::iterator itChimeras;
        
        if (!util.isBlank(accnosFileName)) {
            //edit accnos file
            ifstream in2;
            util.openInputFile(accnosFileName, in2);
            
            while (!in2.eof()) {
                if (m->getControl_pressed()) { in2.close(); out2.close(); util.mothurRemove(outputFileName); util.mothurRemove((accnosFileName+".temp")); return 0; }
                
                in2 >> name; util.gobble(in2);
                
                //find unique name
                itUnique = uniqueNames.find(name);
                
                if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing accnos results. Cannot find " + name + "."); m->mothurOutEndLine(); m->setControl_pressed(true); }
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
        }
        out2.close();
        
        util.mothurRemove(accnosFileName);
        rename((accnosFileName+".temp").c_str(), accnosFileName.c_str());
        
        
        
        //edit chimera file
        ifstream in;
        util.openInputFile(outputFileName, in);
        
        ofstream out;
        util.openOutputFile(outputFileName+".temp", out); out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
        //out << "Score\tQuery\tParentA\tParentB\tIdQM\tIdQA\tIdQB\tIdAB\tIdQT\tLY\tLN\tLA\tRY\tRN\tRA\tDiv\tYN\n";
        
        float temp1;
        string parent1, parent2, parent3, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, flag;
        name = "";
        namesInFile.clear();
        //assumptions - in file each read will always look like - if uchime source is updated, revisit this code.
        /*										1	2	3	4	5	6	7	8	9	10	11	12	13	14	15
         0.000000	F11Fcsw_33372/ab=18/		*	*	*	*	*	*	*	*	*	*	*	*	*	*	N
         0.0000	GQY1XT001C296C;size=356;	*	*	*	*	*	*	*	*	0	0	0	0	0	0	*	N
         0.0469	GQY1XT001CPCVN;size=154;	GQY1XT001C296C;size=356;	GQY1XT001C44N8;size=323;	GQY1XT001C44N8;size=323;	93.8	91.5	92.3	92.6	92.3	4	0	7	9	3	7	1.5	N
         0.018300	F11Fcsw_14980/ab=16/		F11Fcsw_1915/ab=35/	F11Fcsw_6032/ab=42/	79.9	78.7	78.2	78.7	79.2	3	0	5	11	10	20	1.46	N
         */
        
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { in.close(); out.close(); util.mothurRemove((outputFileName+".temp")); return 0; }
            
            bool print = false;
            in >> temp1;	util.gobble(in);
            in >> name;		util.gobble(in);
            in >> parent1;	util.gobble(in);
            in >> parent2;	util.gobble(in);
            in >> parent3;	util.gobble(in);
            in >> temp2 >> temp3 >> temp4 >> temp5 >> temp6 >> temp7 >> temp8 >> temp9 >> temp10 >> temp11 >> temp12 >> temp13 >> flag;
            util.gobble(in);
            
            //parse name - name will look like U68590/ab=1/
            string restOfName = "";
            int pos = name.find_first_of(';');
            if (pos != string::npos) {
                restOfName = name.substr(pos);
                name = name.substr(0, pos);
            }
            
            //find unique name
            itUnique = uniqueNames.find(name);
            
            if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find "+ name + "."); m->mothurOutEndLine(); m->setControl_pressed(true); }
            else {
                name = itUnique->second;
                //is this name already in the file
                itNames = namesInFile.find((name));
                
                if (itNames == namesInFile.end()) { //no not in file
                    if (flag == "N") { //are you really a no??
                        //is this sequence really not chimeric??
                        itChimeras = chimerasInFile.find(name);
                        
                        //then you really are a no so print, otherwise skip
                        if (itChimeras == chimerasInFile.end()) { print = true; }
                    }else{ print = true; }
                }
            }
            
            if (print) {
                out << temp1 << '\t' << name << restOfName << '\t';
                namesInFile.insert(name);
                
                //parse parent1 names
                if (parent1 != "*") {
                    restOfName = "";
                    pos = parent1.find_first_of(';');
                    if (pos != string::npos) {
                        restOfName = parent1.substr(pos);
                        parent1 = parent1.substr(0, pos);
                    }
                    
                    itUnique = uniqueNames.find(parent1);
                    if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find parentA "+ parent1 + "."); m->mothurOutEndLine(); m->setControl_pressed(true); }
                    else {	out << itUnique->second << restOfName << '\t';	}
                }else { out << parent1 << '\t'; }
                
                //parse parent2 names
                if (parent2 != "*") {
                    restOfName = "";
                    pos = parent2.find_first_of(';');
                    if (pos != string::npos) {
                        restOfName = parent2.substr(pos);
                        parent2 = parent2.substr(0, pos);
                    }
                    
                    itUnique = uniqueNames.find(parent2);
                    if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find parentB "+ parent2 + "."); m->mothurOutEndLine(); m->setControl_pressed(true); }
                    else {	out << itUnique->second << restOfName << '\t';	}
                }else { out << parent2 << '\t'; }
                
                //parse parent3 names
                if (parent3 != "*") {
                    restOfName = "";
                    pos = parent3.find_first_of(';');
                    if (pos != string::npos) {
                        restOfName = parent3.substr(pos);
                        parent3 = parent3.substr(0, pos);
                    }
                    
                    itUnique = uniqueNames.find(parent3);
                    if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing chimera results. Cannot find parentC "+ parent3 + "."); m->mothurOutEndLine(); m->setControl_pressed(true); }
                    else {	out << itUnique->second << restOfName << '\t';	}
                }else { out << parent3 << '\t'; }
                
                out << temp2 << '\t' << temp3 << '\t' << temp4 << '\t' << temp5 << '\t' << temp6 << '\t' << temp7 << '\t' << temp8 << '\t' << temp9 << '\t' << temp10 << '\t' << temp11 << '\t' << temp12 << '\t' << temp13 << '\t' << flag << endl;
            }
        }
        in.close();
        out.close();
        
        util.mothurRemove(outputFileName);
        rename((outputFileName+".temp").c_str(), outputFileName.c_str());
        
        
        //edit anls file
        //assumptions - in file each read will always look like - if uchime source is updated, revisit this code.
        /*
         ------------------------------------------------------------------------
         Query   (  179 nt) F21Fcsw_11639/ab=591/
         ParentA (  179 nt) F11Fcsw_6529/ab=1625/
         ParentB (  181 nt) F21Fcsw_12128/ab=1827/
         
         A     1 AAGgAAGAtTAATACaagATGgCaTCatgAGtccgCATgTtcAcatGATTAAAG--gTaTtcCGGTagacGATGGGGATG 78
         Q     1 AAGTAAGACTAATACCCAATGACGTCTCTAGAAGACATCTGAAAGAGATTAAAG--ATTTATCGGTGATGGATGGGGATG 78
         B     1 AAGgAAGAtTAATcCaggATGggaTCatgAGttcACATgTccgcatGATTAAAGgtATTTtcCGGTagacGATGGGGATG 80
         Diffs      N    N    A N?N   N N  NNN  N?NB   N ?NaNNN          B B NN    NNNN
         Votes      0    0    + 000   0 0  000  000+   0 00!000            + 00    0000
         Model   AAAAAAAAAAAAAAAAAAAAAAxBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
         
         A    79 CGTtccATTAGaTaGTaGGCGGGGTAACGGCCCACCtAGtCttCGATggaTAGGGGTTCTGAGAGGAAGGTCCCCCACAT 158
         Q    79 CGTCTGATTAGCTTGTTGGCGGGGTAACGGCCCACCAAGGCAACGATCAGTAGGGGTTCTGAGAGGAAGGTCCCCCACAT 158
         B    81 CGTtccATTAGaTaGTaGGCGGGGTAACGGCCCACCtAGtCAACGATggaTAGGGGTTCTGAGAGGAAGGTCCCCCACAT 160
         Diffs      NNN     N N  N                   N  N BB    NNN
         Votes      000     0 0  0                   0  0 ++    000
         Model   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
         
         A   159 TGGAACTGAGACACGGTCCAA 179
         Q   159 TGGAACTGAGACACGGTCCAA 179
         B   161 TGGAACTGAGACACGGTCCAA 181
         Diffs
         Votes
         Model   BBBBBBBBBBBBBBBBBBBBB
         
         Ids.  QA 76.6%, QB 77.7%, AB 93.7%, QModel 78.9%, Div. +1.5%
         Diffs Left 7: N 0, A 6, Y 1 (14.3%); Right 35: N 1, A 30, Y 4 (11.4%), Score 0.0047
         */
        if (chimealns) {
            ifstream in3;
            util.openInputFile(alnsFileName, in3);
            
            ofstream out3;
            util.openOutputFile(alnsFileName+".temp", out3); out3.setf(ios::fixed, ios::floatfield); out3.setf(ios::showpoint);
            
            name = "";
            namesInFile.clear();
            string line = "";
            
            while (!in3.eof()) {
                if (m->getControl_pressed()) { in3.close(); out3.close(); util.mothurRemove(outputFileName); util.mothurRemove((accnosFileName)); util.mothurRemove((alnsFileName+".temp")); return 0; }
                
                line = "";
                line = util.getline(in3);
                string temp = "";
                
                if (line != "") {
                    istringstream iss(line);
                    iss >> temp;
                    
                    //are you a name line
                    if ((temp == "Query") || (temp == "ParentA") || (temp == "ParentB")) {
                        int spot = 0;
                        for (int i = 0; i < line.length(); i++) {
                            spot = i;
                            if (line[i] == ')') { break; }
                            else { out3 << line[i]; }
                        }
                        
                        if (spot == (line.length() - 1)) { m->mothurOut("[ERROR]: could not line sequence name in line " + line + "."); m->mothurOutEndLine(); m->setControl_pressed(true); }
                        else if ((spot+2) > (line.length() - 1)) { m->mothurOut("[ERROR]: could not line sequence name in line " + line + "."); m->mothurOutEndLine(); m->setControl_pressed(true); }
                        else {
                            out << line[spot] << line[spot+1];
                            
                            name = line.substr(spot+2);
                            
                            //parse name - name will either look like U68590/ab=1/ or U68590
                            string restOfName = "";
                            int pos = name.find_first_of(';');
                            if (pos != string::npos) {
                                restOfName = name.substr(pos);
                                name = name.substr(0, pos);
                            }
                            
                            //find unique name
                            itUnique = uniqueNames.find(name);
                            
                            if (itUnique == uniqueNames.end()) { m->mothurOut("[ERROR]: trouble parsing alns results. Cannot find "+ name + "."); m->mothurOutEndLine();m->setControl_pressed(true);  }
                            else {
                                //only limit repeats on query names
                                if (temp == "Query") {
                                    itNames = namesInFile.find((itUnique->second));
                                    
                                    if (itNames == namesInFile.end()) {
                                        out << itUnique->second << restOfName << endl;
                                        namesInFile.insert((itUnique->second));
                                    }
                                }else { out << itUnique->second << restOfName << endl;  }
                            }
                            
                        }
                        
                    }else { //not need to alter line
                        out3 << line << endl;
                    }
                }else { out3 << endl; }
            }
            in3.close();
            out3.close();
            
            util.mothurRemove(alnsFileName);
            rename((alnsFileName+".temp").c_str(), alnsFileName.c_str());
        }
        
        return total;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraVsearchCommand", "deconvoluteResults");
        exit(1);
    }
}		
//**********************************************************************************************************************
int ChimeraVsearchCommand::readFasta(string filename, map<string, string>& seqs){
    try {
        //create input file for uchime
        //read through fastafile and store info
        ifstream in;
        util.openInputFile(filename, in);
        
        int num = 0;
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { in.close(); return 0; }
            
            Sequence seq(in); util.gobble(in);
            seqs[seq.getName()] = seq.getUnaligned();
            num++;
        }
        in.close();
        
        return num;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraVsearchCommand", "readFasta");
        exit(1);
    }
}
//**********************************************************************************************************************

string ChimeraVsearchCommand::getNamesFile(string& inputFile){
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
        m->errorOut(e, "ChimeraVsearchCommand", "getNamesFile");
        exit(1);
    }
}
//**********************************************************************************************************************
int ChimeraVsearchCommand::driverGroups(string outputFName, string filename, string accnos, string alns, string countlist, int start, int end, vector<string> groups){
    try {
        
        int totalSeqs = 0;
        int numChimeras = 0;
        
        
        ofstream outCountList;
        if (hasCount && dups) { util.openOutputFile(countlist, outCountList); }
        
        for (int i = start; i < end; i++) {
            long start = time(NULL);	 if (m->getControl_pressed()) {  outCountList.close(); util.mothurRemove(countlist); return 0; }
            
            int error;
            long long thisGroupsSeqs = 0;
            if (hasCount) { error = cparser->getSeqs(groups[i], filename, ";size=", ";", thisGroupsSeqs, true); if ((error == 1) || m->getControl_pressed()) {  return 0; } }
            else { error = sparser->getSeqs(groups[i], filename, ";size=", ";", thisGroupsSeqs, true); if ((error == 1) || m->getControl_pressed()) {  return 0; } }
            
            totalSeqs += thisGroupsSeqs;
            driver((outputFName + groups[i]), filename, (accnos+groups[i]), (alns+ groups[i]), numChimeras);
            
            if (m->getControl_pressed()) { return 0; }
            
            //remove file made for vsearch
            if (!m->getDebug()) {  util.mothurRemove(filename);  }
            else { m->mothurOut("[DEBUG]: saving file: " + filename + ".\n"); }
            
            //if we provided a count file with group info and set dereplicate=t, then we want to create a *.pick.count_table
            //This table will zero out group counts for seqs determined to be chimeric by that group.
            if (dups) {
                if (!util.isBlank(accnos+groups[i])) {
                    ifstream in;
                    util.openInputFile(accnos+groups[i], in);
                    string name;
                    if (hasCount) {
                        while (!in.eof()) {
                            in >> name; util.gobble(in);
                            outCountList << name << '\t' << groups[i] << endl;
                        }
                        in.close();
                    }else {
                        map<string, string> thisnamemap = sparser->getNameMap(groups[i]);
                        map<string, string>::iterator itN;
                        ofstream out;
                        util.openOutputFile(accnos+groups[i]+".temp", out);
                        while (!in.eof()) {
                            in >> name; util.gobble(in);
                            itN = thisnamemap.find(name);
                            if (itN != thisnamemap.end()) {
                                vector<string> tempNames; util.splitAtComma(itN->second, tempNames);
                                for (int j = 0; j < tempNames.size(); j++) { out << tempNames[j] << endl; }
                                
                            }else { m->mothurOut("[ERROR]: parsing cannot find " + name + ".\n"); m->setControl_pressed(true); }
                        }
                        out.close();
                        in.close();
                        util.renameFile(accnos+groups[i]+".temp", accnos+groups[i]);
                    }
                    
                }
            }
            
            //append files
            util.appendFiles((outputFName+groups[i]), outputFName); util.mothurRemove((outputFName+groups[i]));
            util.appendFiles((accnos+groups[i]), accnos); util.mothurRemove((accnos+groups[i]));
            if (chimealns) { util.appendFiles((alns+groups[i]), alns); util.mothurRemove((alns+groups[i])); }
            
            m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(thisGroupsSeqs) + " sequences from group " + groups[i] + ".");	m->mothurOutEndLine();
        }
        
        if (hasCount && dups) { outCountList.close(); }
        
        return totalSeqs;
        
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraVsearchCommand", "driverGroups");
        exit(1);
    }
}
//**********************************************************************************************************************

int ChimeraVsearchCommand::driver(string outputFName, string filename, string accnos, string alns, int& numChimeras){
    try {
        
        outputFName = util.getFullPathName(outputFName);
        string outputFNamec = util.getFullPathName(outputFName+"vsearch_out");
        filename = util.getFullPathName(filename);
        alns = util.getFullPathName(alns);
        
        //to allow for spaces in the path
        outputFName = "\"" + outputFName + "\"";
        outputFNamec = "\"" + outputFNamec + "\"";
        filename = "\"" + filename + "\"";
        alns = "\"" + alns + "\"";
        
        vector<char*> cPara;
        string vsearchCommand = "";

        vsearchCommand += vsearchLocation;
        vsearchCommand = "\"" + vsearchCommand + "\" ";
        
        char* tempVsearch;
        tempVsearch= new char[vsearchCommand.length()+1];
        *tempVsearch = '\0';
        strncat(tempVsearch, vsearchCommand.c_str(), vsearchCommand.length());
        cPara.push_back(tempVsearch);
        
        //are you using a reference file
        if (templatefile != "self") {
            string outputFileName = filename.substr(1, filename.length()-2) + ".vsearch_formatted";
            prepFile(filename.substr(1, filename.length()-2), outputFileName);
            filename = outputFileName;
            filename = "\"" + filename + "\"";
            //add reference file
            char* tempRef = new char[5];
            //strcpy(tempRef, "--db");
            *tempRef = '\0'; strncat(tempRef, "--db", 4);
            cPara.push_back(tempRef);
            char* tempR = new char[templatefile.length()+1];
            //strcpy(tempR, templatefile.c_str());
            *tempR = '\0'; strncat(tempR, templatefile.c_str(), templatefile.length());
            cPara.push_back(tempR);
            
            char* tempIn = new char[13];
            *tempIn = '\0'; strncat(tempIn, "--uchime_ref", 12);
            cPara.push_back(tempIn);
            char* temp = new char[filename.length()+1];
            *temp = '\0'; strncat(temp, filename.c_str(), filename.length());
            cPara.push_back(temp);

        }else { //denovo
            char* tempIn = new char[16];
            *tempIn = '\0'; strncat(tempIn, "--uchime_denovo", 15);
            cPara.push_back(tempIn);
            char* temp = new char[filename.length()+1];
            *temp = '\0'; strncat(temp, filename.c_str(), filename.length());
            cPara.push_back(temp);
        }
        
        
        char* tempO = new char[11];
        *tempO = '\0'; strncat(tempO, "--chimeras", 10);
        cPara.push_back(tempO);
        char* tempout = new char[outputFNamec.length()+1];
        *tempout = '\0'; strncat(tempout, outputFNamec.c_str(), outputFNamec.length());
        cPara.push_back(tempout);
        
        char* tempchimeraout = new char[12];
        *tempchimeraout = '\0'; strncat(tempchimeraout, "--uchimeout", 11);
        cPara.push_back(tempchimeraout);
        char* tempoutc = new char[outputFName.length()+1];
        *tempoutc = '\0'; strncat(tempoutc, outputFName.c_str(), outputFName.length());
        cPara.push_back(tempoutc);
        
        char* tempxsize = new char[8];
        *tempxsize = '\0'; strncat(tempxsize, "--xsize", 7);
        cPara.push_back(tempxsize);
        
        if (chimealns) {
            char* tempA = new char[13];
            *tempA = '\0'; strncat(tempA, "--uchimealns", 12);
            //strcpy(tempA, "--uchimealns");
            cPara.push_back(tempA);
            char* tempa = new char[alns.length()+1];
            //strcpy(tempa, alns.c_str());
            *tempa = '\0'; strncat(tempa, alns.c_str(), alns.length());
            cPara.push_back(tempa);
        }
        
        
        if (useAbskew) {
            char* tempskew = new char[9];
            *tempskew = '\0'; strncat(tempskew, "--abskew", 8);
            //strcpy(tempskew, "--abskew");
            cPara.push_back(tempskew);
            char* tempSkew = new char[abskew.length()+1];
            //strcpy(tempSkew, abskew.c_str());
            *tempSkew = '\0'; strncat(tempSkew, abskew.c_str(), abskew.length());
            cPara.push_back(tempSkew);
        }
        
        if (useMinH) {
            char* tempminh = new char[7];
            *tempminh = '\0'; strncat(tempminh, "--minh", 6);
            //strcpy(tempminh, "--minh");
            cPara.push_back(tempminh);
            char* tempMinH = new char[minh.length()+1];
            *tempMinH = '\0'; strncat(tempMinH, minh.c_str(), minh.length());
            //strcpy(tempMinH, minh.c_str());
            cPara.push_back(tempMinH);
        }
        
        if (useMindiv) {
            char* tempmindiv = new char[9];
            *tempmindiv = '\0'; strncat(tempmindiv, "--mindiv", 8);
            //strcpy(tempmindiv, "--mindiv");
            cPara.push_back(tempmindiv);
            char* tempMindiv = new char[mindiv.length()+1];
            *tempMindiv = '\0'; strncat(tempMindiv, mindiv.c_str(), mindiv.length());
            //strcpy(tempMindiv, mindiv.c_str());
            cPara.push_back(tempMindiv);
        }
        
        if (useMindiffs) {
            char* tempmindiv = new char[9];
            *tempmindiv = '\0'; strncat(tempmindiv, "--mindiffs", 10);
            cPara.push_back(tempmindiv);
            char* tempMindiv = new char[mindiffs.length()+1];
            *tempMindiv = '\0'; strncat(tempMindiv, mindiffs.c_str(), mindiffs.length());
            cPara.push_back(tempMindiv);
        }
        
        if (useXn) {
            char* tempxn = new char[5];
            //strcpy(tempxn, "--xn");
            *tempxn = '\0'; strncat(tempxn, "--xn", 4);
            cPara.push_back(tempxn);
            char* tempXn = new char[xn.length()+1];
            //strcpy(tempXn, xn.c_str());
            *tempXn = '\0'; strncat(tempXn, xn.c_str(), xn.length());
            cPara.push_back(tempXn);
        }
        
        if (useDn) {
            char* tempdn = new char[5];
            //strcpy(tempdn, "--dn");
            *tempdn = '\0'; strncat(tempdn, "--dn", 4);
            cPara.push_back(tempdn);
            char* tempDn = new char[dn.length()+1];
            *tempDn = '\0'; strncat(tempDn, dn.c_str(), dn.length());
            //strcpy(tempDn, dn.c_str());
            cPara.push_back(tempDn);
        }
        
        //--threads=1
        char* threads = new char[10];  threads[0] = '\0'; strncat(threads, "--threads", 9);
        cPara.push_back(threads);
        string numProcessors = toString(processors);
        char* tempThreads = new char[numProcessors.length()+1];
        *tempThreads = '\0'; strncat(tempThreads, numProcessors.c_str(), numProcessors.length());
        cPara.push_back(tempThreads);
        
        char** vsearchParameters;
        vsearchParameters = new char*[cPara.size()];
        string commandString = "";
        for (int i = 0; i < cPara.size(); i++) {  vsearchParameters[i] = cPara[i];  commandString += toString(cPara[i]) + " "; }
        //int numArgs = cPara.size();
        
#if defined NON_WINDOWS
#else
        commandString = "\"" + commandString + "\"";
#endif
        
        if (m->getDebug()) { m->mothurOut("[DEBUG]: vsearch command = " + commandString + ".\n"); }
        
        system(commandString.c_str());
        
        //free memory
        for(int i = 0; i < cPara.size(); i++)  {  delete cPara[i];  }
        delete[] vsearchParameters;
        
        //remove "" from filenames
        outputFName = outputFName.substr(1, outputFName.length()-2);
        outputFNamec = outputFNamec.substr(1, outputFNamec.length()-2);
        filename = filename.substr(1, filename.length()-2);
        alns = alns.substr(1, alns.length()-2);
        
        if (m->getControl_pressed()) { return 0; }
        
        //create accnos file from uchime results
        ifstream in;
        util.openInputFile(outputFNamec, in, "no error");
        
        ofstream out;
        util.openOutputFile(accnos, out);
        
        numChimeras = 0;
        while(!in.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            Sequence seq(in); util.gobble(in);
            
            string name = seq.getName();
            
            if (templatefile == "self") {
                name = name.substr(0, name.length()-1); //rip off last ;
                name = name.substr(0, name.find_last_of(';'));
            }
            
            out << name << endl; numChimeras++;
        }
        in.close();
        out.close();
        
        util.mothurRemove(outputFNamec);
        
        //if (templatefile != "self") {  util.mothurRemove(filename); }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraVsearchCommand", "driver");
        exit(1);
    }
}
/**************************************************************************************************/
//uchime can't handle some of the things allowed in mothurs fasta files. This functions "cleans up" the file.
int ChimeraVsearchCommand::prepFile(string filename, string output) {
    try {
        
        ifstream in;
        util.openInputFile(filename, in);
        
        ofstream out;
        util.openOutputFile(output, out);
        
        while (!in.eof()) {
            if (m->getControl_pressed()) { break;  }
            
            Sequence seq(in); util.gobble(in);
            
            if (seq.getName() != "") { seq.printUnAlignedSequence(out); }
        }
        in.close();
        out.close();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraVsearchCommand", "prepFile");
        exit(1);
    }
}
/**************************************************************************************************/

