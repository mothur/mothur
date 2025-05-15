//
//  chimeravsearchcommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 6/16/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#include "chimeravsearchcommand.h"
#include "uniqueseqscommand.h"
#include "sequence.hpp"
#include "systemcommand.h"
#include "degapseqscommand.h"
#include "removeseqscommand.h"

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
        CommandParameter premovechimeras("removechimeras", "Boolean", "", "t", "", "", "","alns",false,false); parameters.push_back(premovechimeras);
        CommandParameter pminh("minh", "Number", "", "0.28", "", "", "","",false,false); parameters.push_back(pminh);
        CommandParameter pmindiv("mindiv", "Number", "", "0.8", "", "", "","",false,false); parameters.push_back(pmindiv);
        CommandParameter pxn("xn", "Number", "", "8.0", "", "", "","",false,false); parameters.push_back(pxn);
        CommandParameter pdn("dn", "Number", "", "1.4", "", "", "","",false,false); parameters.push_back(pdn);
        CommandParameter pmindiffs("mindiffs", "Number", "", "3", "", "", "","",false,false); parameters.push_back(pmindiffs);
        CommandParameter pvsearchlocation("vsearch", "String", "", "", "", "", "","",false,false); parameters.push_back(pvsearchlocation);
        CommandParameter pdups("dereplicate", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pdups);
        
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["chimera"] = tempOutNames;
        outputTypes["accnos"] = tempOutNames;
        outputTypes["fasta"] = tempOutNames;
        outputTypes["alns"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
        
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
        helpString += "The chimera.vsearch command parameters are fasta, name, count, reference, processors, dereplicate, removechimeras, abskew, uchimealns, minh, mindiv, xn, dn, mindiffs.\n";
        helpString += "The fasta parameter allows you to enter the fasta file containing your potentially chimeric sequences, and is required, unless you have a valid current fasta file. \n";
        helpString += "The name parameter allows you to provide a name file, if you are using template=self. \n";
        helpString += "The count parameter allows you to provide a count file, if you are using template=self. When you use a count file with group info and dereplicate=T, mothur will create a *.pick.count_table file containing seqeunces after chimeras are removed. \n";
         helpString += "The group parameter allows you to provide a group file. The group file can be used with a namesfile and reference=self. When checking sequences, only sequences from the same group as the query sequence will be used as the reference. \n";
        helpString += "If the dereplicate parameter is false, then if one group finds the sequence to be chimeric, then all groups find it to be chimeric, default=f.\n";
        helpString += "The reference parameter allows you to enter a reference file containing known non-chimeric sequences, and is required. You may also set template=self, in this case the abundant sequences will be used as potential parents. \n";
        helpString += "The processors parameter allows you to specify how many processors you would like to use.  The default is 1. \n";
        helpString += "The abskew parameter can only be used with template=self. Minimum abundance skew. Default 1.9. Abundance skew is: min [ abund(parent1), abund(parent2) ] / abund(query).\n";
        helpString += "The uchimealns parameter allows you to indicate you would like a file containing multiple alignments of query sequences to parents in human readable format. Alignments show columns with differences that support or contradict a chimeric model.\n";
        helpString += "The removechimeras parameter allows you to indicate you would like to automatically remove the sequences that are flagged as chimeric. Default=t.\n";
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
        
        string howto = "How do I use the dereplicate parameter?"; howtos.push_back(howto);
        string hanswer = "\tThe dereplicate parameter can be used when checking for chimeras by group. If the dereplicate parameter is false, then if one group finds the sequence to be chimeric, then all groups find it to be chimeric, default=f. If you set dereplicate=t, and then when a sequence is found to be chimeric it is removed from it’s group, not the entire dataset.\n\nNote: When you set dereplicate=t, mothur generates a new count table with the chimeras removed and counts adjusted by sample. It is important to note if you set dereplicate=true, do NOT include the count file with the remove.seqs command. For a detailed example, please reference https://mothur.org/wiki/chimera_dereplicate_example/\n"; hanswers.push_back(hanswer);
        
    
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
        else if (type == "fasta") {  pattern = "[filename],[tag],vsearch.fasta"; }
        else if (type == "alns") {  pattern = "[filename],[tag],vsearch.alns"; }
        else if (type == "count") {  pattern = "[filename],[tag],vsearch.count_table-[filename],count_table"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraVsearchCommand", "getOutputPattern");
        exit(1);
    }
}
//***************************************************************************************************************
ChimeraVsearchCommand::ChimeraVsearchCommand(string option) : Command()  {
    try {
        hasCount=false;
        
        //allow user to run help
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

            string path;
            map<string,string>::iterator it = parameters.find("reference");
            //user has given a template file
            if(it != parameters.end()){
                if (it->second == "self") {  templatefile = "self";  }
                else {
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
            
            temp = validParameter.valid(parameters, "removechimeras");            if (temp == "not found") { temp = "t"; }
            removeChimeras = util.isTrue(temp);
            
            minh = validParameter.valid(parameters, "minh");						if (minh == "not found")			{ useMinH = false; minh = "0.28";					}	else{ useMinH = true;			}
            mindiv = validParameter.valid(parameters, "mindiv");					if (mindiv == "not found")			{ useMindiv = false; mindiv = "0.8";				}	else{ useMindiv = true;			}
            xn = validParameter.valid(parameters, "xn");							if (xn == "not found")				{ useXn = false; xn = "8.0";						}	else{ useXn = true;				}
            dn = validParameter.valid(parameters, "dn");							if (dn == "not found")				{ useDn = false; dn = "1.4";						}	else{ useDn = true;				}
            mindiffs = validParameter.valid(parameters, "mindiffs");				if (mindiffs == "not found")				{ useMindiffs = false; mindiffs = "3";							}	else{ useMindiffs = true;				}
            
            temp = validParameter.valid(parameters, "dereplicate"); if (temp == "not found") { temp = "false";			}
            dups = util.isTrue(temp);
            
            vector<string> versionOutputs;
            bool foundTool = false;
            string programName = "vsearch"; programName += EXECUTABLE_EXT;
            
            vsearchLocation = validParameter.validFile(parameters, "vsearch");
            if (vsearchLocation == "not found") {
                vsearchLocation = "";
                foundTool = util.findTool(programName, vsearchLocation, versionOutputs, current->getLocations());
            }
            else {
                //test to make sure vsearch exists
                ifstream in;
                vsearchLocation = util.getFullPathName(vsearchLocation);
                bool ableToOpen = util.openInputFile(vsearchLocation, in, "no error"); in.close();
                if(!ableToOpen) {
                    m->mothurOut(vsearchLocation + " file does not exist or cannot be opened, ignoring.\n"); vsearchLocation = "";
                    foundTool = util.findTool(programName, vsearchLocation, versionOutputs, current->getLocations());
                }else { foundTool = true; }
            }

            
            if (hasName && (templatefile != "self")) { m->mothurOut("You have provided a namefile and the reference parameter is not set to self. I am not sure what reference you are trying to use, aborting.\n");  abort=true; }
            if (hasCount && (templatefile != "self")) { m->mothurOut("You have provided a countfile and the reference parameter is not set to self. I am not sure what reference you are trying to use, aborting.\n");  abort=true; }
            if (hasGroup && (templatefile != "self")) { m->mothurOut("You have provided a group file and the reference parameter is not set to self. I am not sure what reference you are trying to use, aborting.\n");  abort=true; }

            //look for vsearch exe
            path = current->getProgramPath();
            
            if (!foundTool) { abort = true; }
            if (foundTool && !abort) {
                        
                if (versionOutputs.size() != 0) {
                            
                    if (versionOutputs[0] == "vsearch") {
                        if (versionOutputs.size() >= 2) {
                            string version = versionOutputs[1];
                                    
                            int pos = version.find_first_of('_');
                            if (pos != string::npos) { version = version.substr(0, pos); }
                                    
                            if (!util.isVsearchVersionValid(version, "v2.13.5")) {
                                m->mothurOut("[ERROR]: vsearch version found = " + version + ". Mothur requires minimum version v2.13.5. Vsearch is distributed with mothur's executable or available on github https://github.com/torognes/vsearch/releases/, please correct.\n");  abort = true;
                            }else { m->mothurOut("Using vsearch version " + version + ".\n"); }
                        }
                    }
                }
            }
                   
            if (!abort) {
                if ((namefile != "") || (groupfile != "")) { //convert to count
                    
                    string rootFileName = namefile;
                    if (rootFileName == "") { rootFileName = groupfile; }
                    
                    if (outputdir == "") { outputdir = util.hasPath(rootFileName); }
                    map<string, string> variables; variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(rootFileName));
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
//**********************************************************************************************************************
//string outputFName, string filename, string accnos, string alns, int& numChimeras
void driver(vsearchData* params){
    try {
        params->driverOutputFName = params->util.getFullPathName(params->driverOutputFName);
        params->formattedFastaFilename = params->util.getFullPathName(params->formattedFastaFilename);
        params->driverAlns = params->util.getFullPathName(params->driverAlns);
        string outputFNamec = params->driverOutputFName+"vsearch_out";
        
        //to allow for spaces in the path
        params->driverOutputFName = "\"" + params->driverOutputFName + "\"";
        params->formattedFastaFilename = "\"" + params->formattedFastaFilename + "\"";
        params->driverAlns = "\"" + params->driverAlns + "\"";
        outputFNamec = "\"" + outputFNamec + "\"";
        
        vector<char*> cPara;
        string vsearchCommand = params->vsearchLocation;
        vsearchCommand = "\"" + vsearchCommand + "\" ";
        cPara.push_back(params->util.mothurConvert(vsearchCommand));
       
        string fileToRemove = "";
        string numProcessors = toString(params->processors);
        //are you using a reference file
        if (params->templatefile != "self") {
            
            // degap fasta and reference files
            string rootFileName = params->formattedFastaFilename.substr(1, params->formattedFastaFilename.length()-2);
            string outputFileName = rootFileName + ".vsearch_formatted";
            fileToRemove = outputFileName;
            
            //vsearch cant handle some of the things allowed in mothurs fasta files so we remove them
            ifstream in; params->util.openInputFile(rootFileName, in);
            ofstream out; params->util.openOutputFile(outputFileName, out);
            
            while (!in.eof()) {
                if (params->m->getControl_pressed()) { break;  }
                
                Sequence seq(in); gobble(in);
                
                if (seq.getName() != "") { seq.printUnAlignedSequence(out); }
            }
            in.close(); out.close();
            
            params->formattedFastaFilename = outputFileName;
            params->formattedFastaFilename = "\"" + params->formattedFastaFilename + "\"";
            
            //add reference file
            cPara.push_back(params->util.mothurConvert("--db"));
            cPara.push_back(params->util.mothurConvert(params->templatefile));
            
            //add reference file
            cPara.push_back(params->util.mothurConvert("--uchime_ref"));
            cPara.push_back(params->util.mothurConvert(params->formattedFastaFilename));
            
        }else { //denovo
            numProcessors = toString(1);
            
            cPara.push_back(params->util.mothurConvert("--uchime_denovo"));
            cPara.push_back(params->util.mothurConvert(params->formattedFastaFilename));
        }
        
        //output filename
        cPara.push_back(params->util.mothurConvert("--chimeras"));
        cPara.push_back(params->util.mothurConvert(outputFNamec));
          
        cPara.push_back(params->util.mothurConvert("--uchimeout"));
        cPara.push_back(params->util.mothurConvert(params->driverOutputFName));
        
        cPara.push_back(params->util.mothurConvert("--xsize"));
        
        if (params->vars->chimealns) {
            cPara.push_back(params->util.mothurConvert("--uchimealns"));
            cPara.push_back(params->util.mothurConvert(params->driverAlns));
        }
        
        
        if (params->vars->useAbskew) {
            cPara.push_back(params->util.mothurConvert("--abskew"));
            cPara.push_back(params->util.mothurConvert(params->vars->abskew));
        }
        
        if (params->vars->useMinH) {
            cPara.push_back(params->util.mothurConvert("--minh"));
            cPara.push_back(params->util.mothurConvert(params->vars->minh));
        }
        
        if (params->vars->useMindiv) {
            cPara.push_back(params->util.mothurConvert("--mindiv"));
            cPara.push_back(params->util.mothurConvert(params->vars->mindiv));
        }
        
        if (params->vars->useMindiffs) {
            cPara.push_back(params->util.mothurConvert("--mindiffs"));
            cPara.push_back(params->util.mothurConvert(params->vars->mindiffs));
        }
        
        if (params->vars->useXn) {
            cPara.push_back(params->util.mothurConvert("--xn"));
            cPara.push_back(params->util.mothurConvert(params->vars->xn));
        }
        
        if (params->vars->useDn) {
            cPara.push_back(params->util.mothurConvert("--dn"));
            cPara.push_back(params->util.mothurConvert(params->vars->dn));
        }
        
        //--threads
        cPara.push_back(params->util.mothurConvert("--threads"));
        cPara.push_back(params->util.mothurConvert(numProcessors));
        
        char** vsearchParameters;
        vsearchParameters = new char*[cPara.size()];
        string commandString = "";
        for (int i = 0; i < cPara.size(); i++) {  vsearchParameters[i] = cPara[i];  commandString += toString(cPara[i]) + " "; }
        
#if defined NON_WINDOWS
#else
        commandString = "\"" + commandString + "\"";
#endif
        
        if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: vsearch command = " + commandString + ".\n"); }
       
        system(commandString.c_str());
        
        //free memory
        for(int i = 0; i < cPara.size(); i++)  {  delete cPara[i];  }
        delete[] vsearchParameters;
        
        if (fileToRemove != "") { params->util.mothurRemove(fileToRemove); }
        
        //remove "" from filenames
        params->driverOutputFName = params->driverOutputFName.substr(1, params->driverOutputFName.length()-2);
        outputFNamec = outputFNamec.substr(1, outputFNamec.length()-2);
        params->formattedFastaFilename = params->formattedFastaFilename.substr(1, params->formattedFastaFilename.length()-2);
        params->driverAlns = params->driverAlns.substr(1, params->driverAlns.length()-2);
        
        if (params->m->getControl_pressed()) { return; }
        
        //create accnos file from vsearch results
        ifstream in; params->util.openInputFile(outputFNamec, in, "no error");
        ofstream out; params->util.openOutputFile(params->driverAccnos, out);
        
        params->numChimeras = 0;
        while(!in.eof()) {
            
            if (params->m->getControl_pressed()) { break; }
            
            Sequence seq(in); gobble(in);
            
            out << seq.getName() << endl; params->numChimeras++;
        }
        in.close();
        out.close();
        
        params->util.mothurRemove(outputFNamec);
                
        return;
    }
    catch(exception& e) {
        params->m->errorOut(e, "ChimeraVsearchCommand", "driver");
        exit(1);
    }
}
//***************************************************************************************************************

int ChimeraVsearchCommand::execute(){
    try{
        
        if (abort) { if (calledHelp) { return 0; }  return 2;	}

        m->mothurOut("Checking sequences from " + fastafile + " ...\n" );
        
        long start = time(nullptr);
        if (outputdir == "") { outputdir = util.hasPath(fastafile);  }
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(fastafile));
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
            variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(countfile));
            newCountFile = getOutputFileName("count", variables);
        }
        
        vars = new vsearchVariables();
        vars->setBooleans(dups, useAbskew, chimealns, useMinH, useMindiv, useXn, useDn, useMindiffs, hasCount);
        vars->setVariables(abskew, minh, mindiv, xn, dn, mindiffs);
        
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
            map<string, vector<string> > seqs2RemoveByGroup;
            createProcessesGroups(group2Files, outputFileName, newFasta, accnosFileName, alnsFileName, newCountFile, groups, seqs2RemoveByGroup);
            
            
            if (hasCount && dups) {
                CountTable newCount; newCount.readTable(countfile, true, false);
                
                for (map<string, vector<string> >::iterator it = seqs2RemoveByGroup.begin(); it != seqs2RemoveByGroup.end(); it++) {
                    
                    string group = it->first;
                    for (int k = 0; k < it->second.size(); k++) { newCount.setAbund(it->second[k], group, 0); }
                }
                
                newCount.printTable(newCountFile);
            }
            
            if (m->getControl_pressed()) {  for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0;	}
            
            if (!dups) {
                long long numRedund = 0;
                int totalChimeras = deconvoluteResults(outputFileName, accnosFileName, alnsFileName, numRedund);
                
                m->mothurOut("\nIt took " + toString(time(nullptr) - start) + " secs to check your sequences. " + toString(totalChimeras) + " chimeras were found.\n");
                m->mothurOut("The number of sequences checked may be larger than the number of unique sequences because some sequences are found in several samples.\n");
            }else {
                
                if (hasCount) {
                    unordered_set<string> doNotRemove;
                    CountTable c; c.readTable(newCountFile, true, true);
                    //returns non zeroed names
                    vector<string> namesInTable = c.printTable(newCountFile);
                    outputNames.push_back(newCountFile); outputTypes["count"].push_back(newCountFile);
                    
                    for (int i = 0; i < namesInTable.size(); i++) { doNotRemove.insert(namesInTable[i]); }
                    
                    //remove names we want to keep from accnos file.
                    unordered_set<string> accnosNames = util.readAccnos(accnosFileName);
                    ofstream out2;
                    util.openOutputFile(accnosFileName, out2);
                    for (auto it = accnosNames.begin(); it != accnosNames.end(); it++) {
                        if (doNotRemove.count(*it) == 0) {  out2 << (*it) << endl; }
                    }
                    out2.close();
                }
            }
            
            if (m->getControl_pressed()) {  for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0;	}
            
        }else{
            if (m->getControl_pressed()) {  for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0;	}
            
            map<string, vector<string> > dummay;
            vector<string> dummyGroups;
            
            vsearchData* dataBundle = new vsearchData(processors, dummay, outputFileName, vsearchLocation, templatefile, newFasta, countfile, accnosFileName, alnsFileName, "", dummyGroups, vars);

            dataBundle->setDriverNames(outputFileName, alnsFileName, accnosFileName);
            driver(dataBundle);
            
            int numChimeras = dataBundle->numChimeras;
            
            //add headings
            ofstream out; util.openOutputFile(outputFileName+".temp", out);
            out << "Score\tQuery\tParentA\tParentB\tIdQM\tIdQA\tIdQB\tIdAB\tIdQT\tLY\tLN\tLA\tRY\tRN\tRA\tDiv\tYN\n";
            out.close();
            
            util.appendFiles(outputFileName, outputFileName+".temp");
            util.mothurRemove(outputFileName); rename((outputFileName+".temp").c_str(), outputFileName.c_str());
            
            if (m->getControl_pressed()) { for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	} return 0; }
            
            //remove file made for vsearch
            if (templatefile == "self") {  util.mothurRemove(newFasta); }
            
            m->mothurOut("\nIt took " + toString(time(nullptr) - start) + " secs to check your sequences. " + toString(numChimeras) + " chimeras were found.\n");
        }
        
        outputNames.push_back(outputFileName); outputTypes["chimera"].push_back(outputFileName);
        outputNames.push_back(accnosFileName); outputTypes["accnos"].push_back(accnosFileName);
        if (chimealns) { outputNames.push_back(alnsFileName); outputTypes["alns"].push_back(alnsFileName); }
        
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
                        variables["[tag]"] = "denovo";
                        if (templatefile != "self") { variables["[tag]"] = "ref"; }
                        string currentName = getOutputFileName("count", variables);

                        util.renameFile(filenames["count"][0], currentName);
                        util.mothurRemove(filenames["count"][0]);
                        outputNames.push_back(currentName); outputTypes["count"].push_back(currentName);
                    }//else, mothur created a modified count file removing chimeras by sample. No need to include count file on remove.seqs command. Deconvolute function created modified count table already
                }
                
                map<string, string> variables;
                variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(fastafile));
                variables["[tag]"] = "denovo";
                if (templatefile != "self") { variables["[tag]"] = "ref"; }
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

        ofstream out2; util.openOutputFile(accnosFileName+".temp", out2);
        
        int total = 0;
        string name;
        set<string> namesInFile; //this is so if a sequence is found to be chimera in several samples we dont write it to the results file more than once
        set<string>::iterator itNames;
        set<string> chimerasInFile;
        set<string>::iterator itChimeras;
        
        if (!util.isBlank(accnosFileName)) {
            //edit accnos file
            ifstream in2; util.openInputFile(accnosFileName, in2);
            
            while (!in2.eof()) {
                if (m->getControl_pressed()) { in2.close(); out2.close(); util.mothurRemove(outputFileName); util.mothurRemove((accnosFileName+".temp")); return 0; }
                
                in2 >> name; gobble(in2);
            
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
        ifstream in; util.openInputFile(outputFileName, in);
        ofstream out; util.openOutputFile(outputFileName+".temp", out); out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
        //out << "Score\tQuery\tParentA\tParentB\tIdQM\tIdQA\tIdQB\tIdAB\tIdQT\tLY\tLN\tLA\tRY\tRN\tRA\tDiv\tYN\n";
        
        float temp1;
        string parent1, parent2, parent3, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, flag;
        name = ""; namesInFile.clear();
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
            in >> temp1;	gobble(in);
            in >> name;		gobble(in);
            in >> parent1;	gobble(in);
            in >> parent2;	gobble(in);
            in >> parent3;	gobble(in);
            in >> temp2 >> temp3 >> temp4 >> temp5 >> temp6 >> temp7 >> temp8 >> temp9 >> temp10 >> temp11 >> temp12 >> temp13 >> flag;
            gobble(in);
            
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
        in.close(); out.close();
        
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
            ifstream in3; util.openInputFile(alnsFileName, in3);
            ofstream out3; util.openOutputFile(alnsFileName+".temp", out3); out3.setf(ios::fixed, ios::floatfield); out3.setf(ios::showpoint);
            
            name = "";
            namesInFile.clear();
            string line = "";
            
            while (!in3.eof()) {
                if (m->getControl_pressed()) { in3.close(); out3.close(); util.mothurRemove(outputFileName); util.mothurRemove((accnosFileName)); util.mothurRemove((alnsFileName+".temp")); return 0; }
                
                line = ""; line = util.getline(in3);
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
                    }else { out3 << line << endl; } //not need to alter line
                }else { out3 << endl; }
            }
            in3.close(); out3.close();
            
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
        ifstream in; util.openInputFile(filename, in);
        
        int num = 0;
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { in.close(); return 0; }
            
            Sequence seq(in); gobble(in);
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
        
        Command* uniqueCommand = new UniqueSeqsCommand(inputString);
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

int getSeqsVsearch(map<string, int>& nameMap, string thisGroupsFormattedOutputFilename, string tag, string tag2, long long& numSeqs, string thisGroupsFastaFile, MothurOut* m){
    try {
        int error = 0;
        ifstream in; Utils util; util.openInputFile(thisGroupsFastaFile, in);
        
        vector<seqPriorityNode> nameVector;
        map<string, int>::iterator itNameMap;
        while (!in.eof()) {
            if (m->getControl_pressed()) { break; }
            
            Sequence seq(in); gobble(in);
            
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
        m->errorOut(e, "ChimeraVsearchCommand", "getSeqsVsearch");
        exit(1);
    }
}
//**********************************************************************************************************************
//out << ">" << nameMapCount[i].name  << tag << nameMapCount[i].numIdentical << tag2 << endl << nameMapCount[i].seq << endl;
//map<string, vector<string> > parsedFiles, string outputFName, string filename, string accnos, string alns, string countlist
void driverGroups(vsearchData* params){
    try {
        int totalSeqs = 0;
        
        Utils util;
        for (map<string, vector<string> >::iterator it = params->parsedFiles.begin(); it != params->parsedFiles.end(); it++) {
            long start = time(nullptr);
            
            if (params->m->getControl_pressed()) {  return; }
            
            long long thisGroupsSeqs = 0;
            string thisGroup = it->first;
            
            map<string, int> nameMap;
            if (params->vars->hasCount) {
                CountTable ct; ct.readTable(it->second[1], false, true);
                nameMap = ct.getNameMap();
            }
            else { nameMap = util.readNames(it->second[1]); }
            
            int error = getSeqsVsearch(nameMap, params->formattedFastaFilename, ";size=", ";", thisGroupsSeqs, it->second[0], params->m);
            if ((error == 1) || params->m->getControl_pressed()) {  return; }
            
            totalSeqs += thisGroupsSeqs;
            //driver((outputFName + thisGroup), filename, (accnos+thisGroup), (alns+thisGroup), numChimeras);
            params->setDriverNames((params->outputFName + thisGroup), (params->alns+thisGroup), (params->accnos+thisGroup));
            driver(params);
            
            if (params->m->getControl_pressed()) { return; }
            
            //remove file made for vsearch
            if (!params->m->getDebug()) {  util.mothurRemove(params->formattedFastaFilename);  }
            else { params->m->mothurOut("[DEBUG]: saving file: " + params->formattedFastaFilename + ".\n"); }
            
            //if we provided a count file with group info and set dereplicate=t, then we want to create a *.pick.count_table
            //This table will zero out group counts for seqs determined to be chimeric by that group.
            if (params->vars->dups) {
                if (!util.isBlank(params->accnos+thisGroup)) {
                    ifstream in; util.openInputFile(params->accnos+thisGroup, in);
                    string name;
                    if (params->vars->hasCount) {
                        //add group to seqs2
                        vector<string> namesOfChimeras;
                        while (!in.eof()) {
                            in >> name; gobble(in);
                            namesOfChimeras.push_back(name);
                        }
                        in.close();
                        params->seqs2RemoveByGroup[thisGroup] = namesOfChimeras;
                    }else {
                        map<string, string> thisnamemap; util.readNames(it->second[1], thisnamemap);
                        map<string, string>::iterator itN;
                        ofstream out;
                        util.openOutputFile(params->accnos+thisGroup+".temp", out);
                        while (!in.eof()) {
                            in >> name; gobble(in);
                            itN = thisnamemap.find(name);
                            if (itN != thisnamemap.end()) {
                                vector<string> tempNames; util.splitAtComma(itN->second, tempNames);
                                for (int j = 0; j < tempNames.size(); j++) { out << tempNames[j] << endl; }
                                
                            }else { params->m->mothurOut("[ERROR]: parsing cannot find " + name + ".\n"); params->m->setControl_pressed(true); }
                        }
                        out.close();
                        in.close();
                        util.renameFile(params->accnos+thisGroup+".temp", params->accnos+thisGroup);
                    }
                    
                }
            }
            
            //append files
            util.appendFiles((params->outputFName+thisGroup), params->outputFName); util.mothurRemove((params->outputFName+thisGroup));
            util.appendFiles((params->accnos+thisGroup), params->accnos); util.mothurRemove((params->accnos+thisGroup));
            if (params->vars->chimealns) { util.appendFiles((params->alns+thisGroup), params->alns); util.mothurRemove((params->alns+thisGroup)); }
            
            params->m->mothurOut("\nIt took " + toString(time(nullptr) - start) + " secs to check " + toString(thisGroupsSeqs) + " sequences from group " + thisGroup + ".\n");
        }
        params->count = totalSeqs;
    }
    catch(exception& e) {
        params->m->errorOut(e, "ChimeraVsearchCommand", "driverGroups");
        exit(1);
    }
}
/**************************************************************************************************/
//driverGroups(group2Files, outputFileName, newFasta, accnosFileName, alnsFileName, newCountFile);
int ChimeraVsearchCommand::createProcessesGroups(map<string, vector<string> >& groups2Files, string outputFName, string filename, string accnos, string alns, string newCountFile, vector<string> groups, map<string, vector<string> >& seqs2RemoveByGroup) {
    try {
              
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
        vector<vsearchData*> data;
        
        long long num = 0;
        time_t start, end;
        time(&start);
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            string extension = toString(i+1) + ".temp";
            vector<string> thisGroups;
            map<string, vector<string> > thisGroupsParsedFiles;
            for (int j = lines[i+1].start; j < lines[i+1].end; j++) {
                
                map<string, vector<string> >::iterator it = groups2Files.find(groups[j]);
                if (it != groups2Files.end()) {
                    thisGroupsParsedFiles[groups[j]] = (it->second);
                    thisGroups.push_back(groups[j]);
                }
                else { m->mothurOut("[ERROR]: missing files for group " + groups[j] + ", skipping\n"); }
            }
            vsearchData* dataBundle = new vsearchData(processors, thisGroupsParsedFiles, outputFName+extension, vsearchLocation, templatefile, filename+extension, countfile,  accnos+extension, alns+extension, accnos+".byCount."+extension, thisGroups, vars);
            data.push_back(dataBundle);
            
            workerThreads.push_back(new std::thread(driverGroups, dataBundle));
        }
        
        vector<string> thisGroups;
        map<string, vector<string> > thisGroupsParsedFiles;
        for (int j = lines[0].start; j < lines[0].end; j++) {
            map<string, vector<string> >::iterator it = groups2Files.find(groups[j]);
            if (it != groups2Files.end()) {
                thisGroupsParsedFiles[groups[j]] = (it->second);
                thisGroups.push_back(groups[j]);
            }
            else { m->mothurOut("[ERROR]: missing files for group " + groups[j] + ", skipping\n"); }
        }
        vsearchData* dataBundle = new vsearchData(processors, thisGroupsParsedFiles, outputFName, vsearchLocation, templatefile, filename, countfile, accnos, alns, accnos+".byCount.temp", thisGroups, vars);
        driverGroups(dataBundle);
        num = dataBundle->count;
        int numChimeras = dataBundle->numChimeras;
        seqs2RemoveByGroup = dataBundle->seqs2RemoveByGroup;
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            num += data[i]->count;
            numChimeras += data[i]->numChimeras;
            
            for (map<string, vector<string> >::iterator it = data[i]->seqs2RemoveByGroup.begin(); it != data[i]->seqs2RemoveByGroup.end(); it++) {
                
                map<string, vector<string> >::iterator itSanity = seqs2RemoveByGroup.find(it->first);
                if (itSanity == seqs2RemoveByGroup.end()) { //we haven't seen this group, should always be true
                    seqs2RemoveByGroup[it->first] = it->second;
                }
            }
            
            string extension = toString(i+1) + ".temp";
            util.appendFiles((outputFName+extension), outputFName);
            util.mothurRemove((outputFName+extension));
                       
            util.appendFiles((accnos+extension), accnos);
            util.mothurRemove((accnos+extension));
            
            delete data[i];
            delete workerThreads[i];
        }
        delete dataBundle;
        time(&end);
        m->mothurOut("It took " + toString(difftime(end, start)) + " secs to check " + toString(num) + " sequences.\n\n");

        return num;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraUchimeCommand", "createProcessesGroups");
        exit(1);
    }
}
/**************************************************************************************************/

