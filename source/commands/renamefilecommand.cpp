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
        
        abort = false; calledHelp = false;
        
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
        helpString += "The shorten parameter is used to inicate you want mothur to generate output file names for you. For example: stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.shared would become stability.an.shared. Default=true.";
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
        if(option == "help") { help(); abort = true; calledHelp = true; }
        else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
        
        else {
            OptionParser parser(option, setParameters());
            map<string,string> parameters = parser.getParameters();
            
            ValidParameters validParameter;
            
            int numFiles = 0;
            //check for parameters
            phylipfile = validParameter.validFile(parameters, "phylip");
            if (phylipfile == "not open") { m->mothurOut("Ignoring: " + parameters["phylip"]); m->mothurOutEndLine(); phylipfile = ""; }
            else if (phylipfile == "not found") {  phylipfile = "";  }
            if (phylipfile != "") { numFiles++; }
            
            columnfile = validParameter.validFile(parameters, "column");
            if (columnfile == "not open") { m->mothurOut("Ignoring: " + parameters["column"]); m->mothurOutEndLine(); columnfile = ""; }
            else if (columnfile == "not found") {  columnfile = "";  }
            if (columnfile != "") { numFiles++; }
            
            listfile = validParameter.validFile(parameters, "list");
            if (listfile == "not open") { m->mothurOut("Ignoring: " + parameters["list"]); m->mothurOutEndLine(); listfile = ""; }
            else if (listfile == "not found") {  listfile = "";  }
            if (listfile != "") { numFiles++; }
            
            rabundfile = validParameter.validFile(parameters, "rabund");
            if (rabundfile == "not open") { m->mothurOut("Ignoring: " + parameters["rabund"]); m->mothurOutEndLine(); rabundfile = ""; }
            else if (rabundfile == "not found") {  rabundfile = "";  }
            if (rabundfile != "") { numFiles++; }
            
            sabundfile = validParameter.validFile(parameters, "sabund");
            if (sabundfile == "not open") { m->mothurOut("Ignoring: " + parameters["sabund"]); m->mothurOutEndLine(); sabundfile = ""; }
            else if (sabundfile == "not found") {  sabundfile = "";  }
            if (sabundfile != "") { numFiles++; }
            
            namefile = validParameter.validFile(parameters, "name");
            if (namefile == "not open") { m->mothurOut("Ignoring: " + parameters["name"]); m->mothurOutEndLine(); namefile = ""; }
            else if (namefile == "not found") {  namefile = "";  }
            if (namefile != "") { numFiles++; }
            
            groupfile = validParameter.validFile(parameters, "group");
            if (groupfile == "not open") { m->mothurOut("Ignoring: " + parameters["group"]); m->mothurOutEndLine(); groupfile = ""; }
            else if (groupfile == "not found") {  groupfile = "";  }
            if (groupfile != "") { numFiles++; }
            
            countfile = validParameter.validFile(parameters, "count");
            if (countfile == "not open") { m->mothurOut("Ignoring: " + parameters["count"]); m->mothurOutEndLine(); countfile = ""; }
            else if (countfile == "not found") {  countfile = "";  }
            if (countfile != "") { numFiles++; }
            
            designfile = validParameter.validFile(parameters, "design");
            if (designfile == "not open") { m->mothurOut("Ignoring: " + parameters["design"]); m->mothurOutEndLine(); designfile = ""; }
            else if (designfile == "not found") {  designfile = "";  }
            if (designfile != "") { numFiles++; }
            
            inputfile = validParameter.validFile(parameters, "input");
            if (inputfile == "not open") { m->mothurOut("Ignoring: " + parameters["input"]); m->mothurOutEndLine(); inputfile = ""; }
            else if (inputfile == "not found") {  inputfile = "";  }
            if (inputfile != "") { numFiles++; }
            
            treefile = validParameter.validFile(parameters, "tree");
            if (treefile == "not open") { m->mothurOut("Ignoring: " + parameters["tree"]); m->mothurOutEndLine(); treefile = ""; }
            else if (treefile == "not found") {  treefile = "";  }
            if (treefile != "") { numFiles++; }
            
            sharedfile = validParameter.validFile(parameters, "shared");
            if (sharedfile == "not open") { m->mothurOut("Ignoring: " + parameters["shared"]); m->mothurOutEndLine(); sharedfile = ""; }
            else if (sharedfile == "not found") {  sharedfile = "";  }
            if (sharedfile != "") { numFiles++; }
            
            relabundfile = validParameter.validFile(parameters, "relabund");
            if (relabundfile == "not open") { m->mothurOut("Ignoring: " + parameters["relabund"]); m->mothurOutEndLine(); relabundfile = ""; }
            else if (relabundfile == "not found") {  relabundfile = "";  }
            if (relabundfile != "") { numFiles++; }
            
            fastafile = validParameter.validFile(parameters, "fasta");
            if (fastafile == "not open") { m->mothurOut("Ignoring: " + parameters["fasta"]); m->mothurOutEndLine(); fastafile = ""; }
            else if (fastafile == "not found") {  fastafile = "";  }
            if (fastafile != "") { numFiles++; }
            
            qualfile = validParameter.validFile(parameters, "qfile");
            if (qualfile == "not open") { m->mothurOut("Ignoring: " + parameters["qfile"]); m->mothurOutEndLine(); qualfile = ""; }
            else if (qualfile == "not found") {  qualfile = "";  }
            if (qualfile != "") { numFiles++; }
            
            sfffile = validParameter.validFile(parameters, "sff");
            if (sfffile == "not open") { m->mothurOut("Ignoring: " + parameters["sff"]); m->mothurOutEndLine(); sfffile = ""; }
            else if (sfffile == "not found") {  sfffile = "";  }
            if (sfffile != "") { numFiles++; }
            
            oligosfile = validParameter.validFile(parameters, "oligos");
            if (oligosfile == "not open") { m->mothurOut("Ignoring: " + parameters["oligos"]); m->mothurOutEndLine(); oligosfile = ""; }
            else if (oligosfile == "not found") {  oligosfile = "";  }
            if (oligosfile != "") { numFiles++; }
            
            accnosfile = validParameter.validFile(parameters, "accnos");
            if (accnosfile == "not open") { m->mothurOut("Ignoring: " + parameters["accnos"]); m->mothurOutEndLine(); accnosfile = ""; }
            else if (accnosfile == "not found") {  accnosfile = "";  }
            if (accnosfile != "") { numFiles++; }
            
            taxonomyfile = validParameter.validFile(parameters, "taxonomy");
            if (taxonomyfile == "not open") { m->mothurOut("Ignoring: " + parameters["taxonomy"]); m->mothurOutEndLine(); taxonomyfile = ""; }
            else if (taxonomyfile == "not found") {  taxonomyfile = "";  }
            if (taxonomyfile != "") { numFiles++; }
            
            constaxonomyfile = validParameter.validFile(parameters, "constaxonomy");
            if (constaxonomyfile == "not open") { m->mothurOut("Ignoring: " + parameters["constaxonomy"]); m->mothurOutEndLine(); constaxonomyfile = ""; }
            else if (constaxonomyfile == "not found") {  constaxonomyfile = "";  }
            if (constaxonomyfile != "") { numFiles++; }
            
            flowfile = validParameter.validFile(parameters, "flow");
            if (flowfile == "not open") { m->mothurOut("Ignoring: " + parameters["flow"]); m->mothurOutEndLine(); flowfile = ""; }
            else if (flowfile == "not found") {  flowfile = "";  }
            if (flowfile != "") { numFiles++; }
            
            biomfile = validParameter.validFile(parameters, "biom");
            if (biomfile == "not open") { m->mothurOut("Ignoring: " + parameters["biom"]); m->mothurOutEndLine(); biomfile = ""; }
            else if (biomfile == "not found") {  biomfile = "";  }
            if (biomfile != "") { numFiles++; }
            
            summaryfile = validParameter.validFile(parameters, "summary");
            if (summaryfile == "not open") { m->mothurOut("Ignoring: " + parameters["summary"]); m->mothurOutEndLine(); summaryfile = ""; }
            else if (summaryfile == "not found") {  summaryfile = "";  }
            if (summaryfile != "") { numFiles++; }
            
            filefile = validParameter.validFile(parameters, "file");
            if (filefile == "not open") { m->mothurOut("Ignoring: " + parameters["file"]); m->mothurOutEndLine(); filefile = ""; }
            else if (filefile == "not found") {  filefile = "";  }
            if (filefile != "") { numFiles++; }
            
            string temp = validParameter.valid(parameters, "shorten");		if (temp == "not found") { temp = "T"; }
            mothurGenerated = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "deleteold");		if (temp == "not found") { temp = "T"; }
            deleteOld = util.isTrue(temp);
            
            prefix = validParameter.valid(parameters, "prefix");		if (prefix == "not found") { prefix = ""; }
            
            outputfile = validParameter.validPath(parameters, "new");
            if (outputfile == "not found") {
                if (!mothurGenerated) { m->mothurOut("[ERROR]: you must enter an output file name\n");   abort=true; }
                outputfile = "";
            }else { mothurGenerated=false; if (outputdir != "") { outputfile = outputdir + util.getSimpleName(outputfile);  } }
            
            
            if ((!mothurGenerated) && (numFiles > 1)) {
                m->mothurOut("[ERROR]: You cannot use more than one file parameter unless mothur is generating the output filenames for you.\n"); abort= true;
            }
            
            if ((mothurGenerated) && (outputfile != "") && (numFiles != 1)) {
                m->mothurOut("[ERROR]: You must allow mothur to generate the filenames or input one file at a time with a new name, not both.\n"); abort= true;
            }
            
            if (outputdir != "") { outputfile = outputdir + util.getSimpleName(outputfile);  }
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
        
        if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        string newName = outputfile;
        
        //look for file types
        if (fastafile != "") {
            newName = getNewName(fastafile, "fasta");
            renameOrCopy(fastafile, newName);
            current->setFastaFile(newName);
        }
        if (qualfile != "") {
            newName = getNewName(qualfile, "qfile");
            renameOrCopy(qualfile, newName);
            current->setQualFile(newName);
        }
        if (phylipfile != "") {
            newName = getNewName(phylipfile, "phylip");
            renameOrCopy(phylipfile, newName);
            current->setPhylipFile(newName);
        }
        if (columnfile != "") {
            newName = getNewName(columnfile, "column");
            renameOrCopy(columnfile, newName);
            current->setColumnFile(newName);
        }
        if (listfile != "") {
            newName = getNewName(listfile, "list");
            renameOrCopy(listfile, newName);
            current->setListFile(newName);
        }
        if (rabundfile != "") {
            newName = getNewName(rabundfile, "rabund");
            renameOrCopy(rabundfile, newName);
            current->setRabundFile(newName);
        }
        if (sabundfile != "") {
            newName = getNewName(sabundfile, "sabund");
            renameOrCopy(sabundfile, newName);
            current->setSabundFile(newName);
        }
        if (namefile != "") {
            newName = getNewName(namefile, "name");
            renameOrCopy(namefile, newName);
            current->setNameFile(newName);
        }
        if (groupfile != "") {
            newName = getNewName(groupfile, "group");
            renameOrCopy(groupfile, newName);
            current->setGroupFile(newName);
        }
        if (treefile != "") {
            newName = getNewName(treefile, "tree");
            renameOrCopy(treefile, newName);
            current->setTreeFile(newName);
        }
        if (sharedfile != "") {
            newName = getNewName(sharedfile, "shared");
            renameOrCopy(sharedfile, newName);
            current->setSharedFile(newName);
        }
        if (relabundfile != "") {
            newName = getNewName(relabundfile, "relabund");
            renameOrCopy(relabundfile, newName);
            current->setRelAbundFile(newName);
        }
        if (designfile != "") {
            newName = getNewName(designfile, "design");
            renameOrCopy(designfile, newName);
            current->setDesignFile(newName);
        }
        if (sfffile != "") {
            newName = getNewName(sfffile, "sff");
            renameOrCopy(sfffile, newName);
            current->setSFFFile(newName);
        }
        if (oligosfile != "") {
            newName = getNewName(oligosfile, "oligos");
            renameOrCopy(oligosfile, newName);
            current->setOligosFile(newName);
        }
        if (accnosfile != "") {
            newName = getNewName(accnosfile, "accnos");
            renameOrCopy(accnosfile, newName);
            current->setAccnosFile(newName);
        }
        if (taxonomyfile != "") {
            newName = getNewName(taxonomyfile, "taxonomy");
            renameOrCopy(taxonomyfile, newName);
            current->setTaxonomyFile(newName);
        }
        if (constaxonomyfile != "") {
            newName = getNewName(constaxonomyfile, "constaxonomy");
            renameOrCopy(constaxonomyfile, newName);
            current->setConsTaxonomyFile(newName);
        }
        if (flowfile != "") {
            newName = getNewName(flowfile, "flow");
            renameOrCopy(flowfile, newName);
            current->setFlowFile(newName);
        }
        if (biomfile != "") {
            newName = getNewName(biomfile, "biom");
            renameOrCopy(biomfile, newName);
            current->setBiomFile(newName);
        }
        if (countfile != "") {
            newName = getNewName(countfile, "count");
            renameOrCopy(countfile, newName);
            current->setCountFile(newName);
        }
        if (summaryfile != "") {
            newName = getNewName(summaryfile, "summary");
            renameOrCopy(summaryfile, newName);
            current->setSummaryFile(newName);
        }
        if (filefile != "") {
            newName = getNewName(filefile, "file");
            renameOrCopy(filefile, newName);
            current->setFileFile(newName);
        }
        if (inputfile != "") {
            newName = getNewName(inputfile, "input");
            renameOrCopy(inputfile, newName);
        }
        
        m->mothurOutEndLine(); m->mothurOut("Current files saved by mothur:\n"); 
        if (current->hasCurrentFiles()) {  current->printCurrentFiles(""); }
        
        return 0;	
    }
    
    catch(exception& e) {
        m->errorOut(e, "RenameFileCommand", "execute");
        exit(1);
    }
}
//**********************************************************************************************************************
                         
string RenameFileCommand::getNewName(string inputFileName, string type){
    try {
        string newName = outputfile;
        inputFileName = util.getFullPathName(inputFileName);
        
        if (mothurGenerated) {
            string extension = util.getExtension(inputFileName);
            string basicName = "final";
            string tag = "";
            
            if (prefix == "") {
                int pos = inputFileName.find_first_of(".");
                if (pos != string::npos) { basicName = util.getSimpleName(inputFileName.substr(0, pos)); }
            }else { basicName = prefix; }
            
            if ((type == "shared") || (type == "list") || (type == "relabund") || (type == "rabund") || (type == "sabund")) {
                vector<string> tags; tags.push_back(".an."); tags.push_back(".tx.");  tags.push_back(".agc."); tags.push_back(".dgc."); tags.push_back(".nn."); tags.push_back(".fn."); tags.push_back(".wn."); tags.push_back(".opti_");
                
                for (int i = 0; i < tags.size(); i++) {
                    int pos2 = inputFileName.find(tags[i]);
                    if (pos2 != string::npos) {
                        int pos3 = inputFileName.substr(pos2+1).find_first_of('.');
                        tag = inputFileName.substr(pos2+1, pos3);
                        break;
                    }
                }
            }else if (type == "constaxonomy") {
                extension = ".cons.taxonomy";
            }
            
            string thisOutputDir = outputdir;
            if (outputdir == "") {  thisOutputDir += util.hasPath(inputFileName);  }
            
            newName = thisOutputDir + basicName;
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
        if (deleteOld) {  util.renameFile(oldName, newName); }
        else {
            string command = "copy ";
            
            #if defined NON_WINDOWS
                command = "cp ";
            #endif
            
            string inputString = command + oldName + " " + newName;
            m->mothurOut("/******************************************/\n"); 
            m->mothurOut("Running command: system(" + inputString + ")\n"); 
            current->setMothurCalling(true);
            
            Command* systemCommand = new SystemCommand(inputString);
            systemCommand->execute();
            delete systemCommand;
            current->setMothurCalling(false);
            m->mothurOut("/******************************************/\n"); 
        }
        
        return newName;
    }
    
    catch(exception& e) {
        m->errorOut(e, "RenameFileCommand", "renameOrCopy");
        exit(1);
    }
}
//**********************************************************************************************************************



