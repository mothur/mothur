//
//  renamefilecommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/18/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#include "renamefilecommand.h"
#include "systemcommand.h"

//**********************************************************************************************************************
vector<string> RenameFileCommand::setParameters(){
    try {
        CommandParameter pflow("flow", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pflow);
        CommandParameter pfile("file", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pfile);
        CommandParameter pbiom("biom", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pbiom);
        CommandParameter pphylip("phylip", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pphylip);
        CommandParameter pcolumn("column", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pcolumn);
        CommandParameter psummary("summary", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(psummary);
        CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pfasta);
        CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pname);
        CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pgroup);
        CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(plist);
        CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(ptaxonomy);
        CommandParameter pqfile("qfile", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pqfile);
        CommandParameter paccnos("accnos", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(paccnos);
        CommandParameter prabund("rabund", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(prabund);
        CommandParameter psabund("sabund", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(psabund);
        CommandParameter pdesign("design", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pdesign);
        CommandParameter porder("order", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(porder);
        CommandParameter ptree("tree", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(ptree);
        CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pshared);
        CommandParameter pcount("count", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pcount);
        CommandParameter poutputname("new", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(poutputname);
        CommandParameter pinputname("input", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pinputname);
        CommandParameter prelabund("relabund", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(prelabund);
        CommandParameter psff("sff", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(psff);
        CommandParameter pconstaxonomy("constaxonomy", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(pconstaxonomy);
        CommandParameter poligos("oligos", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(poligos);
        CommandParameter pmothurgenerated("shorten", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pmothurgenerated);
        CommandParameter pdeleteold("deleteold", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pdeleteold);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        CommandParameter pprefix("prefix", "String", "", "", "", "", "","",false,false); parameters.push_back(pprefix);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        vector<string> myArray;
        for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
        return myArray;
    }
    catch(exception& e) {
        m->errorOut(e, "RenameFileCommand", "setParameters");
        exit(1);
    }
}
//**********************************************************************************************************************
string RenameFileCommand::getHelpString(){
    try {
        string helpString = "";
        helpString += "The rename.file command allows you to rename files and updates the current files saved by mothur.\n";
        helpString += "The rename.file command parameters are: phylip, column, list, rabund, sabund, name, group, design, tree, shared, relabund, fasta, qfile, sff, oligos, accnos, biom, count, summary, file, taxonomy, constaxonomy, input, new, prefix, deletedold and shorten.\n";
        helpString += "The new parameter allows you to provide an output file name for the input file you provide.\n";
        helpString += "The shorten parameter is used to inicate you want mothur to generate output file names for you. For example: stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.shared would become stability.an.shared. Default=true.";
        helpString += "The prefix parameter allows you to enter your own prefix for shortened names.";
        helpString += "The deleteold parameter indicates whether you want to delete the old file.  Default=true.";
        helpString += "The rename.file command should be in the following format: \n";
        helpString += "rename.file(fasta=current, name=current, group=current, taxonomy=current, shorten=t)\n";
        return helpString;
    }
    catch(exception& e) {
        m->errorOut(e, "RenameFileCommand", "getHelpString");
        exit(1);
    }
}


//**********************************************************************************************************************
RenameFileCommand::RenameFileCommand(){
    try {
        abort = true; calledHelp = true;
        setParameters();
        vector<string> tempOutNames;
    }
    catch(exception& e) {
        m->errorOut(e, "RenameFileCommand", "RenameFileCommand");
        exit(1);
    }
}
//**********************************************************************************************************************
string RenameFileCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "RenameFileCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
RenameFileCommand::RenameFileCommand(string option)  {
    try {
        abort = false; calledHelp = false;
        
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
            outputTypes["summary"] = tempOutNames;
            
            //if the user changes the input directory command factory will send this info to us in the output parameter
            string inputDir = validParameter.validFile(parameters, "inputdir", false);
            if (inputDir == "not found"){	inputDir = "";		}
            else {
                string path;
                it = parameters.find("phylip");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["phylip"] = inputDir + it->second;		}
                }
                
                it = parameters.find("column");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["column"] = inputDir + it->second;		}
                }
                
                it = parameters.find("fasta");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
                }
                
                it = parameters.find("list");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["list"] = inputDir + it->second;		}
                }
                
                it = parameters.find("rabund");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["rabund"] = inputDir + it->second;		}
                }
                
                it = parameters.find("sabund");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["sabund"] = inputDir + it->second;		}
                }
                
                it = parameters.find("name");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["name"] = inputDir + it->second;		}
                }
                
                it = parameters.find("group");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["group"] = inputDir + it->second;		}
                }
                
                it = parameters.find("design");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["design"] = inputDir + it->second;		}
                }
                
                
                it = parameters.find("tree");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["tree"] = inputDir + it->second;		}
                }
                
                it = parameters.find("shared");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["shared"] = inputDir + it->second;		}
                }
                
                it = parameters.find("input");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["input"] = inputDir + it->second;		}
                }
                
                it = parameters.find("count");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["count"] = inputDir + it->second;		}
                }
                
                it = parameters.find("relabund");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["relabund"] = inputDir + it->second;		}
                }
                
                it = parameters.find("fasta");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
                }
                
                it = parameters.find("qfile");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["qfile"] = inputDir + it->second;		}
                }
                
                it = parameters.find("sff");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["sff"] = inputDir + it->second;		}
                }
                
                it = parameters.find("oligos");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["oligos"] = inputDir + it->second;		}
                }
                
                it = parameters.find("accnos");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["accnos"] = inputDir + it->second;		}
                }
                
                it = parameters.find("taxonomy");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
                }
                
                it = parameters.find("flow");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["flow"] = inputDir + it->second;		}
                }
                
                it = parameters.find("biom");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["biom"] = inputDir + it->second;		}
                }
                
                it = parameters.find("summary");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["summary"] = inputDir + it->second;		}
                }
                
                it = parameters.find("file");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["file"] = inputDir + it->second;		}
                }
                
                it = parameters.find("constaxonomy");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["constaxonomy"] = inputDir + it->second;		}
                }
            }
            
            //if the user changes the output directory command factory will send this info to us in the output parameter
            outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){  outputDir = ""; }
            
            int numFiles = 0;
            //check for parameters
            phylipfile = validParameter.validFile(parameters, "phylip", true);
            if (phylipfile == "not open") { m->mothurOut("Ignoring: " + parameters["phylip"]); m->mothurOutEndLine(); phylipfile = ""; }
            else if (phylipfile == "not found") {  phylipfile = "";  }
            if (phylipfile != "") { numFiles++; }
            
            columnfile = validParameter.validFile(parameters, "column", true);
            if (columnfile == "not open") { m->mothurOut("Ignoring: " + parameters["column"]); m->mothurOutEndLine(); columnfile = ""; }
            else if (columnfile == "not found") {  columnfile = "";  }
            if (columnfile != "") { numFiles++; }
            
            listfile = validParameter.validFile(parameters, "list", true);
            if (listfile == "not open") { m->mothurOut("Ignoring: " + parameters["list"]); m->mothurOutEndLine(); listfile = ""; }
            else if (listfile == "not found") {  listfile = "";  }
            if (listfile != "") { numFiles++; }
            
            rabundfile = validParameter.validFile(parameters, "rabund", true);
            if (rabundfile == "not open") { m->mothurOut("Ignoring: " + parameters["rabund"]); m->mothurOutEndLine(); rabundfile = ""; }
            else if (rabundfile == "not found") {  rabundfile = "";  }
            if (rabundfile != "") { numFiles++; }
            
            sabundfile = validParameter.validFile(parameters, "sabund", true);
            if (sabundfile == "not open") { m->mothurOut("Ignoring: " + parameters["sabund"]); m->mothurOutEndLine(); sabundfile = ""; }
            else if (sabundfile == "not found") {  sabundfile = "";  }
            if (sabundfile != "") { numFiles++; }
            
            namefile = validParameter.validFile(parameters, "name", true);
            if (namefile == "not open") { m->mothurOut("Ignoring: " + parameters["name"]); m->mothurOutEndLine(); namefile = ""; }
            else if (namefile == "not found") {  namefile = "";  }
            if (namefile != "") { numFiles++; }
            
            groupfile = validParameter.validFile(parameters, "group", true);
            if (groupfile == "not open") { m->mothurOut("Ignoring: " + parameters["group"]); m->mothurOutEndLine(); groupfile = ""; }
            else if (groupfile == "not found") {  groupfile = "";  }
            if (groupfile != "") { numFiles++; }
            
            countfile = validParameter.validFile(parameters, "count", true);
            if (countfile == "not open") { m->mothurOut("Ignoring: " + parameters["count"]); m->mothurOutEndLine(); countfile = ""; }
            else if (countfile == "not found") {  countfile = "";  }
            if (countfile != "") { numFiles++; }
            
            designfile = validParameter.validFile(parameters, "design", true);
            if (designfile == "not open") { m->mothurOut("Ignoring: " + parameters["design"]); m->mothurOutEndLine(); designfile = ""; }
            else if (designfile == "not found") {  designfile = "";  }
            if (designfile != "") { numFiles++; }
            
            inputfile = validParameter.validFile(parameters, "input", true);
            if (inputfile == "not open") { m->mothurOut("Ignoring: " + parameters["input"]); m->mothurOutEndLine(); inputfile = ""; }
            else if (inputfile == "not found") {  inputfile = "";  }
            if (inputfile != "") { numFiles++; }
            
            treefile = validParameter.validFile(parameters, "tree", true);
            if (treefile == "not open") { m->mothurOut("Ignoring: " + parameters["tree"]); m->mothurOutEndLine(); treefile = ""; }
            else if (treefile == "not found") {  treefile = "";  }
            if (treefile != "") { numFiles++; }
            
            sharedfile = validParameter.validFile(parameters, "shared", true);
            if (sharedfile == "not open") { m->mothurOut("Ignoring: " + parameters["shared"]); m->mothurOutEndLine(); sharedfile = ""; }
            else if (sharedfile == "not found") {  sharedfile = "";  }
            if (sharedfile != "") { numFiles++; }
            
            relabundfile = validParameter.validFile(parameters, "relabund", true);
            if (relabundfile == "not open") { m->mothurOut("Ignoring: " + parameters["relabund"]); m->mothurOutEndLine(); relabundfile = ""; }
            else if (relabundfile == "not found") {  relabundfile = "";  }
            if (relabundfile != "") { numFiles++; }
            
            fastafile = validParameter.validFile(parameters, "fasta", true);
            if (fastafile == "not open") { m->mothurOut("Ignoring: " + parameters["fasta"]); m->mothurOutEndLine(); fastafile = ""; }
            else if (fastafile == "not found") {  fastafile = "";  }
            if (fastafile != "") { numFiles++; }
            
            qualfile = validParameter.validFile(parameters, "qfile", true);
            if (qualfile == "not open") { m->mothurOut("Ignoring: " + parameters["qfile"]); m->mothurOutEndLine(); qualfile = ""; }
            else if (qualfile == "not found") {  qualfile = "";  }
            if (qualfile != "") { numFiles++; }
            
            sfffile = validParameter.validFile(parameters, "sff", true);
            if (sfffile == "not open") { m->mothurOut("Ignoring: " + parameters["sff"]); m->mothurOutEndLine(); sfffile = ""; }
            else if (sfffile == "not found") {  sfffile = "";  }
            if (sfffile != "") { numFiles++; }
            
            oligosfile = validParameter.validFile(parameters, "oligos", true);
            if (oligosfile == "not open") { m->mothurOut("Ignoring: " + parameters["oligos"]); m->mothurOutEndLine(); oligosfile = ""; }
            else if (oligosfile == "not found") {  oligosfile = "";  }
            if (oligosfile != "") { numFiles++; }
            
            accnosfile = validParameter.validFile(parameters, "accnos", true);
            if (accnosfile == "not open") { m->mothurOut("Ignoring: " + parameters["accnos"]); m->mothurOutEndLine(); accnosfile = ""; }
            else if (accnosfile == "not found") {  accnosfile = "";  }
            if (accnosfile != "") { numFiles++; }
            
            taxonomyfile = validParameter.validFile(parameters, "taxonomy", true);
            if (taxonomyfile == "not open") { m->mothurOut("Ignoring: " + parameters["taxonomy"]); m->mothurOutEndLine(); taxonomyfile = ""; }
            else if (taxonomyfile == "not found") {  taxonomyfile = "";  }
            if (taxonomyfile != "") { numFiles++; }
            
            constaxonomyfile = validParameter.validFile(parameters, "constaxonomy", true);
            if (constaxonomyfile == "not open") { m->mothurOut("Ignoring: " + parameters["constaxonomy"]); m->mothurOutEndLine(); constaxonomyfile = ""; }
            else if (constaxonomyfile == "not found") {  constaxonomyfile = "";  }
            if (constaxonomyfile != "") { numFiles++; }
            
            flowfile = validParameter.validFile(parameters, "flow", true);
            if (flowfile == "not open") { m->mothurOut("Ignoring: " + parameters["flow"]); m->mothurOutEndLine(); flowfile = ""; }
            else if (flowfile == "not found") {  flowfile = "";  }
            if (flowfile != "") { numFiles++; }
            
            biomfile = validParameter.validFile(parameters, "biom", true);
            if (biomfile == "not open") { m->mothurOut("Ignoring: " + parameters["biom"]); m->mothurOutEndLine(); biomfile = ""; }
            else if (biomfile == "not found") {  biomfile = "";  }
            if (biomfile != "") { numFiles++; }
            
            summaryfile = validParameter.validFile(parameters, "summary", true);
            if (summaryfile == "not open") { m->mothurOut("Ignoring: " + parameters["summary"]); m->mothurOutEndLine(); summaryfile = ""; }
            else if (summaryfile == "not found") {  summaryfile = "";  }
            if (summaryfile != "") { numFiles++; }
            
            filefile = validParameter.validFile(parameters, "file", true);
            if (filefile == "not open") { m->mothurOut("Ignoring: " + parameters["file"]); m->mothurOutEndLine(); filefile = ""; }
            else if (filefile == "not found") {  filefile = "";  }
            if (filefile != "") { numFiles++; }
            
            string temp = validParameter.validFile(parameters, "shorten", false);		if (temp == "not found") { temp = "T"; }
            mothurGenerated = m->isTrue(temp);
            
            temp = validParameter.validFile(parameters, "deleteold", false);		if (temp == "not found") { temp = "T"; }
            deleteOld = m->isTrue(temp);
            
            prefix = validParameter.validFile(parameters, "prefix", false);		if (prefix == "not found") { prefix = ""; }
            
            outputfile = validParameter.validFile(parameters, "new", false);
            if (outputfile == "not found") {
                if (!mothurGenerated) { m->mothurOut("[ERROR]: you must enter an output file name"); m->mothurOutEndLine();  abort=true; }
                outputfile = "";
            }else { mothurGenerated=false; if (outputDir != "") { outputfile = outputDir + m->getSimpleName(outputfile);  } }
            
            
            if ((!mothurGenerated) && (numFiles > 1)) {
                m->mothurOut("[ERROR]: You cannot use more than one file parameter unless mothur is generating the output filenames for you.\n"); abort= true;
            }
            
            if ((mothurGenerated) && (outputfile != "") && (numFiles != 1)) {
                m->mothurOut("[ERROR]: You must allow mothur to generate the filenames or input one file at a time with a new name, not both.\n"); abort= true;
            }
            
            if (outputDir != "") { outputfile = outputDir + m->getSimpleName(outputfile);  }
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "RenameFileCommand", "RenameFileCommand");
        exit(1);
    }
}
//**********************************************************************************************************************

int RenameFileCommand::execute(){
    try {
        
        if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        string newName = outputfile;
        
        //look for file types
        if (fastafile != "") {
            newName = getNewName(fastafile, "fasta");
            renameOrCopy(fastafile, newName);
            m->setFastaFile(newName);
        }
        if (qualfile != "") {
            newName = getNewName(qualfile, "qfile");
            renameOrCopy(qualfile, newName);
            m->setQualFile(newName);
        }
        if (phylipfile != "") {
            newName = getNewName(phylipfile, "phylip");
            renameOrCopy(phylipfile, newName);
            m->setPhylipFile(newName);
        }
        if (columnfile != "") {
            newName = getNewName(columnfile, "column");
            renameOrCopy(columnfile, newName);
            m->setColumnFile(newName);
        }
        if (listfile != "") {
            newName = getNewName(listfile, "list");
            renameOrCopy(listfile, newName);
            m->setListFile(newName);
        }
        if (rabundfile != "") {
            newName = getNewName(rabundfile, "rabund");
            renameOrCopy(rabundfile, newName);
            m->setRabundFile(newName);
        }
        if (sabundfile != "") {
            newName = getNewName(sabundfile, "sabund");
            renameOrCopy(sabundfile, newName);
            m->setSabundFile(newName);
        }
        if (namefile != "") {
            newName = getNewName(namefile, "name");
            renameOrCopy(namefile, newName);
            m->setNameFile(newName);
        }
        if (groupfile != "") {
            newName = getNewName(groupfile, "group");
            renameOrCopy(groupfile, newName);
            m->setGroupFile(newName);
        }
        if (treefile != "") {
            newName = getNewName(treefile, "tree");
            renameOrCopy(treefile, newName);
            m->setTreeFile(newName);
        }
        if (sharedfile != "") {
            newName = getNewName(sharedfile, "shared");
            renameOrCopy(sharedfile, newName);
            m->setSharedFile(newName);
        }
        if (relabundfile != "") {
            newName = getNewName(relabundfile, "relabund");
            renameOrCopy(relabundfile, newName);
            m->setRelAbundFile(newName);
        }
        if (designfile != "") {
            newName = getNewName(designfile, "design");
            renameOrCopy(designfile, newName);
            m->setDesignFile(newName);
        }
        if (sfffile != "") {
            newName = getNewName(sfffile, "sff");
            renameOrCopy(sfffile, newName);
            m->setSFFFile(newName);
        }
        if (oligosfile != "") {
            newName = getNewName(oligosfile, "oligos");
            renameOrCopy(oligosfile, newName);
            m->setOligosFile(newName);
        }
        if (accnosfile != "") {
            newName = getNewName(accnosfile, "accnos");
            renameOrCopy(accnosfile, newName);
            m->setAccnosFile(newName);
        }
        if (taxonomyfile != "") {
            newName = getNewName(taxonomyfile, "taxonomy");
            renameOrCopy(taxonomyfile, newName);
            m->setTaxonomyFile(newName);
        }
        if (constaxonomyfile != "") {
            newName = getNewName(constaxonomyfile, "constaxonomy");
            renameOrCopy(constaxonomyfile, newName);
            m->setConsTaxonomyFile(newName);
        }
        if (flowfile != "") {
            newName = getNewName(flowfile, "flow");
            renameOrCopy(flowfile, newName);
            m->setFlowFile(newName);
        }
        if (biomfile != "") {
            newName = getNewName(biomfile, "biom");
            renameOrCopy(biomfile, newName);
            m->setBiomFile(newName);
        }
        if (countfile != "") {
            newName = getNewName(countfile, "count");
            renameOrCopy(countfile, newName);
            m->setCountTableFile(newName);
        }
        if (summaryfile != "") {
            newName = getNewName(summaryfile, "summary");
            renameOrCopy(summaryfile, newName);
            m->setSummaryFile(newName);
        }
        if (filefile != "") {
            newName = getNewName(filefile, "file");
            renameOrCopy(filefile, newName);
            m->setFileFile(newName);
        }
        if (inputfile != "") {
            newName = getNewName(inputfile, "input");
            renameOrCopy(inputfile, newName);
        }
        
        m->mothurOutEndLine(); m->mothurOut("Current files saved by mothur:"); m->mothurOutEndLine();
        if (m->hasCurrentFiles()) {  m->printCurrentFiles(""); }
        
        return 0;	
    }
    
    catch(exception& e) {
        m->errorOut(e, "RenameFileCommand", "execute");
        exit(1);
    }
}
//**********************************************************************************************************************
                         
string RenameFileCommand::getNewName(string name, string type){
    try {
        string newName = outputfile;
        name = m->getFullPathName(name);
        
        if (mothurGenerated) {
            string extension = m->getExtension(name);
            string basicName = "final";
            string tag = "";
            
            if (prefix == "") {
                int pos = name.find_first_of(".");
                if (pos != string::npos) { basicName = name.substr(0, pos); }
            }else { basicName = prefix; }
            
            if ((type == "shared") || (type == "list") || (type == "relabund") || (type == "rabund") || (type == "sabund")) {
                vector<string> tags; tags.push_back(".an."); tags.push_back(".tx.");  tags.push_back(".agc."); tags.push_back(".dgc."); tags.push_back(".nn."); tags.push_back(".fn."); tags.push_back(".wn."); tags.push_back(".opti_");
                
                for (int i = 0; i < tags.size(); i++) {
                    int pos2 = name.find(tags[i]);
                    if (pos2 != string::npos) {
                        int pos3 = name.substr(pos2+1).find_first_of('.');
                        tag = name.substr(pos2+1, pos3);
                        break;
                    }
                }
            }else if (type == "constaxonomy") {
                extension = ".cons.taxonomy";
            }
            
            newName = basicName;
            if (tag != "") { newName += "." + tag; }
            newName += extension;
        }
        
        return newName;
    }
    
    catch(exception& e) {
        m->errorOut(e, "RenameFileCommand", "getNewFileName");
        exit(1);
    }
}
//**********************************************************************************************************************

string RenameFileCommand::renameOrCopy(string oldName, string newName){
    try {
        if (deleteOld) {  m->renameFile(oldName, newName); }
        else {
            string command = "copy ";
            
            #if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
                command = "cp ";
            #endif
            
            string inputString = command + oldName + " " + newName;
            m->mothurOut("/******************************************/"); m->mothurOutEndLine();
            m->mothurOut("Running command: system(" + inputString + ")"); m->mothurOutEndLine();
            m->setMothurCalling(true);
            
            Command* systemCommand = new SystemCommand(inputString);
            systemCommand->execute();
            delete systemCommand;
            m->setMothurCalling(false);
            m->mothurOut("/******************************************/"); m->mothurOutEndLine();
        }
        
        return newName;
    }
    
    catch(exception& e) {
        m->errorOut(e, "RenameFileCommand", "renameOrCopy");
        exit(1);
    }
}
//**********************************************************************************************************************



