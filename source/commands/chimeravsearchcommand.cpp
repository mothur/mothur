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
        CommandParameter pvsearchlocation("vsearch", "String", "", "", "", "", "","",false,false); parameters.push_back(pvsearchlocation);
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
        helpString += "The vsearch parameter allows you to specify the name and location of your vsearch executable. By default mothur will look in your path and mothur's executable and mothur tools locations.  You can set the vsearch location as follows, vsearch=/usr/bin/vsearch.\n";
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
string ChimeraVsearchCommand::getCommonQuestions(){
    try {
        vector<string> questions, issues, qanswers, ianswers, howtos, hanswers;
        
        string issue = "... vsearch file does not exist. mothur requires the vsearch executable."; issues.push_back(issue);
        string ianswer = "\tThe chimera.vsearch command is a wrapper for the vsearch program, https://github.com/torognes/vsearch. We distribute the vsearch executable with the executable versions of mothur. By default, mothur will look for vsearch in the same location mothur's executable is as well as looking in your $PATH variable.\n"; ianswers.push_back(ianswer);
        
        string commonQuestions = util.getFormattedHelp(questions, qanswers, issues, ianswers, howtos, hanswers);
        
        return commonQuestions;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraVsearchCommand", "getCommonQuestions");
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
        else if (type == "count") {  pattern = "[filename],[tag],vsearch.pick.count_table-[filename],count_table"; }
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
        abort = false; calledHelp = false; hasCount=false;
        
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
                if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
            }
            
            vector<string> tempOutNames;
            outputTypes["chimera"] = tempOutNames;
            outputTypes["accnos"] = tempOutNames;
            outputTypes["alns"] = tempOutNames;
            outputTypes["count"] = tempOutNames;
            
            //if the user changes the input directory command factory will send this info to us in the output parameter
            string inputDir = validParameter.valid(parameters, "inputdir");
            if (inputDir == "not found"){	inputDir = "";		}
            else {
                string path;
                it = parameters.find("count");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["count"] = inputDir + it->second;		}
                }
                
                it = parameters.find("fasta");
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
                }
                
                it = parameters.find("name");
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["name"] = inputDir + it->second;		}
                }
                
                it = parameters.find("group");
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["group"] = inputDir + it->second;		}
                }
            }
            
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
            
            //check for required parameters
            countfile = validParameter.validFile(parameters, "count");
            if (countfile == "not open") { countfile = ""; abort = true; }
            else if (countfile == "not found") { countfile = "";  }
            else { current->setCountFile(countfile); }
            if (countfile != "") { hasCount = true; }
            
            //make sure there is at least one valid file left
            if (hasName && hasCount) { m->mothurOut("[ERROR]: You must enter ONLY ONE of the following: count or name.\n");  abort = true; }
            
            bool hasGroup = false;
            string groupfile = validParameter.validFile(parameters, "group");
            if (groupfile == "not open") { abort = true; }
            else if (groupfile == "not found") {  groupfile = "";  }
            else { current->setGroupFile(groupfile); hasGroup = true; }
            
            if (hasGroup && hasCount) { m->mothurOut("[ERROR]: You must enter ONLY ONE of the following: count or group.\n");  abort = true; }
            
            //if the user changes the output directory command factory will send this info to us in the output parameter
            outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";	}
            
            string path;
            it = parameters.find("reference");
            //user has given a template file
            if(it != parameters.end()){
                if (it->second == "self") {  templatefile = "self";  }
                else {
                    string path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["reference"] = inputDir + it->second;		}
                    
                    templatefile = validParameter.validFile(parameters, "reference");
                    if (templatefile == "not open") { abort = true; }
                    else if (templatefile == "not found") { //check for saved reference sequences
                        m->mothurOut("[ERROR]: The reference parameter is a required, aborting.\n"); abort = true;
                    }
                }
            }else if ((hasName) || (hasCount) || (hasGroup)) {  templatefile = "self"; }
            else {  m->mothurOut("[ERROR]: The reference parameter is a required, aborting.\n"); templatefile = ""; abort = true; }
            
            string temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
            processors = current->setProcessors(temp);
            
            abskew = validParameter.valid(parameters, "abskew");	if (abskew == "not found"){	useAbskew = false;  abskew = "1.9";	}else{  useAbskew = true;  }
            if (useAbskew && templatefile != "self") { m->mothurOut("The abskew parameter is only valid with template=self, ignoring.\n");  useAbskew = false; }
            
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
            
            vector<string> versionOutputs;
            bool foundTool = false;
            path = current->getProgramPath();
            string programName = "vsearch"; programName += EXECUTABLE_EXT;
            
            vsearchLocation = validParameter.valid(parameters, "vsearch");
            if (vsearchLocation == "not found") {
                vsearchLocation = "";
                foundTool = util.findTool(programName, vsearchLocation, path, versionOutputs, current->getLocations());
            }
            else {
                //test to make sure vsearch exists
                ifstream in;
                vsearchLocation = util.getFullPathName(vsearchLocation);
                bool ableToOpen = util.openInputFile(vsearchLocation, in, "no error"); in.close();
                if(!ableToOpen) {
                    m->mothurOut(vsearchLocation + " file does not exist or cannot be opened, ignoring.\n"); vsearchLocation = "";
                    programName = util.getSimpleName(vsearchLocation); vsearchLocation = "";
                    foundTool = util.findTool(programName, vsearchLocation, path, versionOutputs, current->getLocations());
                }
            }

            
            if (hasName && (templatefile != "self")) { m->mothurOut("You have provided a namefile and the reference parameter is not set to self. I am not sure what reference you are trying to use, aborting.\n");  abort=true; }
            if (hasCount && (templatefile != "self")) { m->mothurOut("You have provided a countfile and the reference parameter is not set to self. I am not sure what reference you are trying to use, aborting.\n");  abort=true; }
            if (hasGroup && (templatefile != "self")) { m->mothurOut("You have provided a group file and the reference parameter is not set to self. I am not sure what reference you are trying to use, aborting.\n");  abort=true; }

            //look for vsearch exe
            path = current->getProgramPath();
            
            if (foundTool && !abort) {
                        
                if (versionOutputs.size() != 0) {
                            
                    if (versionOutputs[0] == "vsearch") {
                        if (versionOutputs.size() >= 2) {
                            string version = versionOutputs[1];
                                    
                            int pos = version.find_first_of('_');
                            if (pos != string::npos) { version = version.substr(0, pos); }
                                    
                            if (version != "v2.13.3") {
                                m->mothurOut("[ERROR]: vsearch version found = " + version + ". Mothur requires version v2.13.3 which is distributed with mothur's executable or available on github https://github.com/torognes/vsearch/releases/tag/v2.13.3, please correct. \n");  abort = true;
                            }else { m->mothurOut("Using vsearch version " + version + ".\n"); }
                        }
                    }
                }
            }
                   
            if (!abort) {
                if ((namefile != "") || (groupfile != "")) { //convert to count
                    
                    string rootFileName = namefile;
                    if (rootFileName == "") { rootFileName = groupfile; }
                    
                    if (outputDir == "") { outputDir = util.hasPath(rootFileName); }
                    map<string, string> variables; variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(rootFileName));
                    string outputFileName = getOutputFileName("count", variables);
                    
                    CountTable ct; ct.createTable(namefile, groupfile, nullVector); ct.printCompressedTable(outputFileName);
                    outputNames.push_back(outputFileName);
                    
                    current->setCountFile(outputFileName);
                    countfile = outputFileName;
                    hasCount = true;
                }
            }
            
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

        m->mothurOut("Checking sequences from " + fastafile + " ...\n" );
        
        long start = time(NULL);
        if (outputDir == "") { outputDir = util.hasPath(fastafile);  }//if user entered a file with a path then preserve it
        map<string, string> variables;
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastafile));
        variables["[tag]"] = "denovo";
        if (templatefile != "self") { variables["[tag]"] = "ref"; }
        string outputFileName = getOutputFileName("chimera", variables);
        string accnosFileName = getOutputFileName("accnos", variables);
        string alnsFileName = getOutputFileName("alns", variables);
        string newFasta = util.getRootName(fastafile) + "temp";
        string newCountFile = "";
        
        //you provided a groupfile
        bool hasGroups = false;
        int numSeqs = 0;
        if (hasCount) {
            CountTable ct;
            if (ct.testGroups(countfile)) { hasGroups = true; }
            variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(countfile));
            newCountFile = getOutputFileName("count", variables);
        }
        
        //setup fasta file if denovo and no groups
        if ((templatefile == "self") && (!hasGroups)) { //you want to run vsearch with a template=self and no groups
            
            if (processors != 1) { m->mothurOut("When using template=self, mothur can only use 1 processor, continuing.\n"); processors = 1; }
            
            if (hasCount) { }
            else { countfile = getCountFile(fastafile); hasCount = true; }
            
            map<string, string> seqs;
            numSeqs = readFasta(fastafile, seqs);  if (m->getControl_pressed()) { for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0; }
            
            //read namefile
            vector<seqPriorityNode> nameMapCount;
            //int error;
            if (hasCount) {
                CountTable ct; ct.readTable(countfile, true, false);
                for(map<string, string>::iterator it = seqs.begin(); it != seqs.end(); it++) {
                    int num = ct.getNumSeqs(it->first);
                    if (num != 0) { seqPriorityNode temp(num, it->second, it->first); nameMapCount.push_back(temp); }
                }
            }
            
            if (seqs.size() != nameMapCount.size()) { m->mothurOut( "The number of sequences in your fastafile does not match the number of sequences in your namefile, aborting.\n");  for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0; }
            
            util.printVsearchFile(nameMapCount, newFasta, ";size=", ";");
            fastafile = newFasta;
        }
        
        if (m->getControl_pressed()) {  for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0;	}
        
        if (hasGroups) {
    
            //Parse sequences by group
            vector<string> groups;
            map<string, vector<string> > group2Files;
            if (hasCount) {
                current->setMothurCalling(true);
                SequenceCountParser cparser(countfile, fastafile, nullVector);
                current->setMothurCalling(false);
                groups = cparser.getNamesOfGroups();
                group2Files = cparser.getFiles();
                
            }
            
            if (m->getControl_pressed()) { for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0; }
            
            //clears files
            ofstream out, out1, out2;
            util.openOutputFile(outputFileName, out); out.close();
            util.openOutputFile(accnosFileName, out1); out1.close();
            if (chimealns) { util.openOutputFile(alnsFileName, out2); out2.close(); }
            
            //paralellized in vsearch
            driverGroups(group2Files, outputFileName, newFasta, accnosFileName, alnsFileName, newCountFile);
            if (hasCount && dups) {
                CountTable c; c.readTable(countfile, true, false);
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
                long long numRedund = 0;
                int totalChimeras = deconvoluteResults(outputFileName, accnosFileName, alnsFileName, numRedund);
                
                m->mothurOut("\nIt took " + toString(time(NULL) - start) + " secs to check your sequences. " + toString(totalChimeras) + " chimeras were found.\n");
                m->mothurOut("The number of sequences checked may be larger than the number of unique sequences because some sequences are found in several samples.\n");
            }else {
                
                if (hasCount) {
                    set<string> doNotRemove;
                    CountTable c; c.readTable(newCountFile, true, true);
                    //returns non zeroed names
                    vector<string> namesInTable = c.printTable(newCountFile);
                    outputNames.push_back(newCountFile); outputTypes["count"].push_back(newCountFile);
                    
                    for (int i = 0; i < namesInTable.size(); i++) { doNotRemove.insert(namesInTable[i]); }
                    
                    //remove names we want to keep from accnos file.
                    set<string> accnosNames = util.readAccnos(accnosFileName);
                    ofstream out2;
                    util.openOutputFile(accnosFileName, out2);
                    for (set<string>::iterator it = accnosNames.begin(); it != accnosNames.end(); it++) {
                        if (doNotRemove.count(*it) == 0) {  out2 << (*it) << endl; }
                    }
                    out2.close();
                }
            }
            
            if (m->getControl_pressed()) {  for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0;	}
            
        }else{
            if (m->getControl_pressed()) {  for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0;	}
            
            int numChimeras = 0;
            
            //paralellized in vsearch
            driver(outputFileName, fastafile, accnosFileName, alnsFileName, numChimeras);
            
            //add headings
            ofstream out;
            util.openOutputFile(outputFileName+".temp", out);
            out << "Score\tQuery\tParentA\tParentB\tIdQM\tIdQA\tIdQB\tIdAB\tIdQT\tLY\tLN\tLA\tRY\tRN\tRA\tDiv\tYN\n";
            out.close();
            
            util.appendFiles(outputFileName, outputFileName+".temp");
            util.mothurRemove(outputFileName); rename((outputFileName+".temp").c_str(), outputFileName.c_str());
            
            if (m->getControl_pressed()) { for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	} return 0; }
            
            //remove file made for vsearch
            if (templatefile == "self") {  util.mothurRemove(fastafile); }
            
            m->mothurOut("\nIt took " + toString(time(NULL) - start) + " secs to check your sequences. " + toString(numChimeras) + " chimeras were found.\n");
        }
        
        outputNames.push_back(outputFileName); outputTypes["chimera"].push_back(outputFileName);
        outputNames.push_back(accnosFileName); outputTypes["accnos"].push_back(accnosFileName);
        if (chimealns) { outputNames.push_back(alnsFileName); outputTypes["alns"].push_back(alnsFileName); }
        
        
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
        
        m->mothurOut("\nOutput File Names:\n");
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
int ChimeraVsearchCommand::deconvoluteResults(string outputFileName, string accnosFileName, string alnsFileName, long long& numRedund){
    try {
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
                
                
                itChimeras = chimerasInFile.find(name);
                
                if (itChimeras == chimerasInFile.end()) {
                    out2 << name << endl;
                    chimerasInFile.insert(name);
                    total++;
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
        //assumptions - in file each read will always look like - if vsearch source is updated, revisit this code.
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
            
            
            if (print) {
                namesInFile.insert(name);
                out << temp1 << '\t' << name << '\t' << parent1 << '\t' << parent2 << '\t' << parent3 << '\t' << temp2 << '\t' << temp3 << '\t' << temp4 << '\t' << temp5 << '\t' << temp6 << '\t' << temp7 << '\t' << temp8 << '\t' << temp9 << '\t' << temp10 << '\t' << temp11 << '\t' << temp12 << '\t' << temp13 << '\t' << flag << endl;
            }
        }
        in.close();
        out.close();
        
        util.mothurRemove(outputFileName);
        rename((outputFileName+".temp").c_str(), outputFileName.c_str());
        
        
        //edit anls file
        //assumptions - in file each read will always look like - if vsearch source is updated, revisit this code.
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
                        
                        if (spot == (line.length() - 1)) { m->mothurOut("[ERROR]: could not line sequence name in line " + line + ".\n");  m->setControl_pressed(true); }
                        else if ((spot+2) > (line.length() - 1)) { m->mothurOut("[ERROR]: could not line sequence name in line " + line + ".\n");  m->setControl_pressed(true); }
                        else {
                            out << line[spot] << line[spot+1];
                            
                            name = line.substr(spot+2);
                                                            //only limit repeats on query names
                            if (temp == "Query") {
                                itNames = namesInFile.find(name);
                                
                                if (itNames == namesInFile.end()) {
                                    out << name << endl;
                                    namesInFile.insert(name);
                                }
                            }else { out << name << endl;  }
                            
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
        //create input file for vsearch
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

string ChimeraVsearchCommand::getCountFile(string& inputFile){
    try {
        string countFile = "";
        
        m->mothurOut("\nNo namesfile given, running unique.seqs command to generate one.\n\n");
        
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
        m->errorOut(e, "ChimeraVsearchCommand", "getCountFile");
        exit(1);
    }
}


//**********************************************************************************************************************

int ChimeraVsearchCommand::getSeqs(map<string, int>& nameMap, string thisGroupsFormattedOutputFilename, string tag, string tag2, long long& numSeqs, string thisGroupsFastaFile){
    try {
        int error = 0;
        ifstream in;
        util.openInputFile(thisGroupsFastaFile, in);
        
        vector<seqPriorityNode> nameVector;
        map<string, int>::iterator itNameMap;
        while (!in.eof()) {
            if (m->getControl_pressed()) { break; }
            
            Sequence seq(in); util.gobble(in);
            
            itNameMap = nameMap.find(seq.getName());
            
            if (itNameMap == nameMap.end()){
                error = 1;
                m->mothurOut("[ERROR]: " + seq.getName() + " is in your fastafile, but is not in your name or count file, please correct.\n");
            }else {
                int num = itNameMap->second;
                
                seqPriorityNode temp(num, seq.getUnaligned(), seq.getName());
                nameVector.push_back(temp);
            }
        }
        in.close();
        
        if (error == 1) {  return 1; }
        
        numSeqs = nameVector.size();
        
        util.printVsearchFile(nameVector, thisGroupsFormattedOutputFilename, tag, tag2);
        
        return error;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraVsearchCommand", "getSeqs");
        exit(1);
    }
}
//**********************************************************************************************************************
//out << ">" << nameMapCount[i].name  << tag << nameMapCount[i].numIdentical << tag2 << endl << nameMapCount[i].seq << endl;
int ChimeraVsearchCommand::driverGroups(map<string, vector<string> > parsedFiles, string outputFName, string filename, string accnos, string alns, string countlist){
    try {
        
        int totalSeqs = 0;
        int numChimeras = 0;
        
        ofstream outCountList;
        if (hasCount && dups) { util.openOutputFile(countlist, outCountList); }
        
        for (map<string, vector<string> >::iterator it = parsedFiles.begin(); it != parsedFiles.end(); it++) {
            long start = time(NULL);	 if (m->getControl_pressed()) {  outCountList.close(); util.mothurRemove(countlist); return 0; }
            
            int error;
            long long thisGroupsSeqs = 0;
            string thisGroup = it->first;
            
            map<string, int> nameMap;
            if (hasCount) {
                CountTable ct; ct.readTable(it->second[1], false, true);
                nameMap = ct.getNameMap();
            }
            else { nameMap = util.readNames(it->second[1]); }
            
            error = getSeqs(nameMap, filename, ";size=", ";", thisGroupsSeqs, it->second[0]);
            if ((error == 1) || m->getControl_pressed()) {  return 0; }
            
            totalSeqs += thisGroupsSeqs;
            driver((outputFName + thisGroup), filename, (accnos+thisGroup), (alns+thisGroup), numChimeras);
            
            if (m->getControl_pressed()) { return 0; }
            
            //remove file made for vsearch
            if (!m->getDebug()) {  util.mothurRemove(filename);  }
            else { m->mothurOut("[DEBUG]: saving file: " + filename + ".\n"); }
            
            //if we provided a count file with group info and set dereplicate=t, then we want to create a *.pick.count_table
            //This table will zero out group counts for seqs determined to be chimeric by that group.
            if (dups) {
                if (!util.isBlank(accnos+thisGroup)) {
                    ifstream in;
                    util.openInputFile(accnos+thisGroup, in);
                    string name;
                    if (hasCount) {
                        while (!in.eof()) {
                            in >> name; util.gobble(in);
                            outCountList << name << '\t' << thisGroup << endl;
                        }
                        in.close();
                    }else {
                        map<string, string> thisnamemap; util.readNames(it->second[1], thisnamemap);
                        map<string, string>::iterator itN;
                        ofstream out;
                        util.openOutputFile(accnos+thisGroup+".temp", out);
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
                        util.renameFile(accnos+thisGroup+".temp", accnos+thisGroup);
                    }
                    
                }
            }
            
            //append files
            util.appendFiles((outputFName+thisGroup), outputFName); util.mothurRemove((outputFName+thisGroup));
            util.appendFiles((accnos+thisGroup), accnos); util.mothurRemove((accnos+thisGroup));
            if (chimealns) { util.appendFiles((alns+thisGroup), alns); util.mothurRemove((alns+thisGroup)); }
            
            m->mothurOut("\nIt took " + toString(time(NULL) - start) + " secs to check " + toString(thisGroupsSeqs) + " sequences from group " + thisGroup + ".\n");
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
        
        string fileToRemove = "";
        //are you using a reference file
        if (templatefile != "self") {
            string outputFileName = filename.substr(1, filename.length()-2) + ".vsearch_formatted";
            fileToRemove = outputFileName;
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
        
        /*char* tempRandSeed = new char[11];
        *tempRandSeed = '\0'; strncat(tempRandSeed, "--randseed", 10);
        cPara.push_back(tempRandSeed);
        string mothurSeed = toString(m->getRandomSeed());
        char* tempRS = new char[mothurSeed.length()+1];
        *tempRS = '\0'; strncat(tempRS, mothurSeed.c_str(), mothurSeed.length());
        cPara.push_back(tempRS);*/
        
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
        
        if (fileToRemove != "") { util.mothurRemove(fileToRemove); }
        
        //remove "" from filenames
        outputFName = outputFName.substr(1, outputFName.length()-2);
        outputFNamec = outputFNamec.substr(1, outputFNamec.length()-2);
        filename = filename.substr(1, filename.length()-2);
        alns = alns.substr(1, alns.length()-2);
        
        if (m->getControl_pressed()) { return 0; }
        
        //create accnos file from vsearch results
        ifstream in;
        util.openInputFile(outputFNamec, in, "no error");
        
        ofstream out;
        util.openOutputFile(accnos, out);
        
        numChimeras = 0;
        while(!in.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            Sequence seq(in); util.gobble(in);
            
            out << seq.getName() << endl; numChimeras++;
        }
        in.close();
        out.close();
        
        util.mothurRemove(outputFNamec);
                
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraVsearchCommand", "driver");
        exit(1);
    }
}
/**************************************************************************************************/
//vsearch can't handle some of the things allowed in mothurs fasta files. This functions "cleans up" the file.
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

