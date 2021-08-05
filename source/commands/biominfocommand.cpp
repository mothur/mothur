//
//  biominfocommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 8/5/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#include "biominfocommand.h"

//**********************************************************************************************************************
vector<string> BiomInfoCommand::setParameters(){
    try {
        CommandParameter pbiom("biom", "InputTypes", "", "", "", "", "","",false,true, true); parameters.push_back(pbiom);
        CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
        CommandParameter prelabund("relabund", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(prelabund);
        CommandParameter pbasis("basis", "Multiple", "otu-sequence", "otu", "", "", "","",false,false); parameters.push_back(pbasis);
        CommandParameter pformat("format", "Multiple", "hdf5-simple", "hdf5", "", "", "","",false,false, true); parameters.push_back(pformat);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter poutput("output", "Multiple", "simple-detail", "detail", "", "", "","",false,false, true); parameters.push_back(poutput);
        CommandParameter pprintlevel("printlevel", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pprintlevel);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false; maxLevel = 0;
        
        vector<string> tempOutNames;
        outputTypes["taxonomy"] = tempOutNames;
        outputTypes["shared"] = tempOutNames;
        outputTypes["constaxonomy"] = tempOutNames;
        outputTypes["taxsummary"] = tempOutNames;
        
        vector<string> myArray;
        for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
        return myArray;
    }
    catch(exception& e) {
        m->errorOut(e, "BiomInfoCommand", "setParameters");
        exit(1);
    }
}
//**********************************************************************************************************************
string BiomInfoCommand::getHelpString(){
    try {
        string helpString = "";
        helpString += "The biom.info command reads a biom file creates a shared file. If your biom file contains metadata mothur will also create taxonomy or constaxonomy along with tax.summary files.\n";
        helpString += "The biom.info command parameters are " + getCommandParameters() + ". The biom parameter is required.\n";
        helpString += "The format parameter allows you indicate type of biom file you have. Options hdf5 or simple. Default is hdf5, unless you are running a version without HDF5 libraries.\n";
        helpString += "The label parameter allows you to enter a distance label to be used in the shared file created from your biom file.\n";
        helpString += "The relabund parameter allows you to indicate you want the tax.summary file values to be relative abundances rather than raw abundances. Default=F. \n";
        helpString += "The basis parameter allows you indicate what you want the summary file to represent, options are otu and sequence. Default is otu.\n";
        helpString += "The output parameter allows you to specify format of your summary file. Options are simple and detail. The default is detail.\n";
        helpString += "The printlevel parameter allows you to specify taxlevel of your summary file to print to. Options are 1 to the maz level in the file.  The default is -1, meaning max level.  If you select a level greater than the level your sequences classify to, mothur will print to the level your max level. \n";
        helpString += "For example consider the following basis=sequence could give Clostridiales 3 105, where 105 is the total number of sequences whose OTU classified to Clostridiales. ";
        helpString += "Now for basis=otu could give Clostridiales 3 7, where 7 is the number of OTUs that classified to Clostridiales.\n";
        helpString += "The biom.info command should be in the following format: biom.info(biom=test.biom, label=0.03).\n";
        
        getCommonQuestions();
        
        return helpString;
    }
    catch(exception& e) {
        m->errorOut(e, "BiomInfoCommand", "getHelpString");
        exit(1);
    }
}
//**********************************************************************************************************************
string BiomInfoCommand::getCommonQuestions(){
    try {
        vector<string> questions, issues, qanswers, ianswers, howtos, hanswers;
        
        string issue = "Cannot convert error. What do I do?"; issues.push_back(issue);
        string ianswer = "\tThis issue is caused by a matrix_element_type mismatch. The biom file contains a field called 'matrix_element_type'. This field tells mothur what form your observation data is in: int or float. Mothur expects 'int' (an interger value) because the shared file contains interger value abundance counts. If your file contains float values mothur will round down to the nearest integer value. But if your matrix_element_type=int and yet the file contains integer counts in float form, (ie. 31.0 instead of 31) you will get this error. You can resolve this issue by setting matrix_element_type=float in the biom file.\n"; ianswers.push_back(ianswer);
        
        issue = "Mothur can't read my biom file. What does this mean?"; issues.push_back(issue);
        ianswer = "\tMothur allows for 2 formats: classic (http://biom-format.org/documentation/format_versions/biom-1.0.html) and hdf5 (http://biom-format.org/documentation/format_versions/biom-2.0.html). NOTE: you can only process hdf5 files if you are using our pre-built version or have built your version of mothur with USEHDF5=yes.\n"; ianswers.push_back(ianswer);
        
        string commonQuestions = util.getFormattedHelp(questions, qanswers, issues, ianswers, howtos, hanswers);
        
        return commonQuestions;
    }
    catch(exception& e) {
        m->errorOut(e, "BiomInfoCommand", "getCommonQuestions");
        exit(1);
    }
}

//**********************************************************************************************************************
string BiomInfoCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "shared") {  pattern = "[filename],[tag],shared"; }
        else if (type == "constaxonomy") {  pattern = "[filename],[tag],cons.taxonomy"; }
        else if (type == "taxonomy") {  pattern = "[filename],[tag],taxonomy"; }
        else if (type == "taxsummary") {  pattern = "[filename],[tag],[tag2],tax.summary"; } //tag2 = "" for taxonomy tag2 = cons for constaxonomy
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "BiomInfoCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
BiomInfoCommand::BiomInfoCommand(string option) : Command()  {
    try {
        maxLevel = 0;
        
        //allow user to run help
        if(option == "help") { help(); abort = true; calledHelp = true; }
        else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
        
        else {
            OptionParser parser(option, setParameters());
            map<string, string> parameters = parser.getParameters();
            
            ValidParameters validParameter;
            //check for required parameters
            biomfile = validParameter.validFile(parameters, "biom");
            if (biomfile == "not open") { biomfile = ""; abort = true; }
            else if (biomfile == "not found") { biomfile = ""; m->mothurOut("[ERROR]: You must provide a biom file, please correct.\n");  abort = true;}
            else { current->setBiomFile(biomfile); }
            
            label = validParameter.valid(parameters, "label");
            if (label == "not found") { label = "userLabel"; }
            
            output = validParameter.valid(parameters, "output");		if(output == "not found"){	output = "detail"; }
            if ((output != "simple") && (output != "detail")) { m->mothurOut(output + " is not a valid output form. Options are simple and detail. I will use detail.\n");  output = "detail"; }
            
            string temp = validParameter.valid(parameters, "relabund");
            if (temp == "not found"){	temp = "false";			}
            else { temp = util.getSimpleName(temp); }
            relabund = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "printlevel");		if (temp == "not found"){	temp = "-1";		}
            util.mothurConvert(temp, printlevel);
            
            basis = validParameter.valid(parameters, "basis");
            if (basis == "not found") { basis = "otu"; }
            
            if ((basis != "otu") && (basis != "sequence")) { m->mothurOut("Invalid option for basis. basis options are otu and sequence, using otu.\n"); }
            
            format = validParameter.valid(parameters, "format");
            if (format == "not found") {
                #ifdef USE_HDF5
                if (!abort) {
                    if (util.isHDF5(biomfile)) { format = "hdf5"; }
                    else { format = "simple"; }
                }
                #else
                    format = "simple";
                #endif
            }
            
            if ((format != "hdf5") && (format != "simple")) {
                 m->mothurOut("Invalid option for format. format options are hdf5 and simple, quitting.\n"); abort = true;
            }
            
            if (format == "hdf5") {
                #ifdef USE_HDF5
                //do nothing we have the api
                #else
                    m->mothurOut("[ERROR]: To read HDF5 biom files, you must have the API installed, quitting.\n"); abort=true;
                #endif
            }
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "BiomInfoCommand", "BiomInfoCommand");
        exit(1);
    }
}
//**********************************************************************************************************************

int BiomInfoCommand::execute(){
    try {
        
        if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        long start = time(NULL);
        
        Biom* biom;
        if (format == "hdf5")   { biom = new BiomHDF5(biomfile, label);   }
        else                   { biom = new BiomSimple(biomfile, label);   }
        
        //getting output filename
        string filename = biomfile;
        if (outputdir == "") { outputdir += util.hasPath(filename); }
        fileroot = outputdir + util.getRootName(util.getSimpleName(biomfile));
        
        SharedRAbundVectors* shared = biom->getSharedRAbundVectors();
        if (format == "hdf5")   {  label = shared->getLabel(); }
        if (label == "") { label = "userLabel"; shared->setLabels(label); }
        
        CountTable ct;
        vector< map<string, bool> > otuContainsGroups;
        
        if (shared != NULL) {
            map<string, string> variables;
            variables["[filename]"] = fileroot;
            variables["[tag]"] = label;
            string sharedFilename = getOutputFileName("shared",variables);
            outputNames.push_back(sharedFilename); outputTypes["shared"].push_back(sharedFilename);
            bool printHeaders = true;
            
            ofstream out; util.openOutputFile(sharedFilename, out);
            shared->print(out, printHeaders);
            out.close();
            
            vector<string> groupNames = shared->getNamesGroups();
            for (int j = 0; j < groupNames.size(); j++) {  ct.addGroup(groupNames[j]); }
                
            int numBins = shared->getNumBins();
            
            for (int i = 0; i < numBins; i++) {
                int total = 0;
                map<string, bool> containsGroup;
                vector<int> abunds;
                for (int j = 0; j < shared->size(); j++) {
                    if (m->getControl_pressed()) { break; }
                    int abund = shared->get(i, groupNames[j]);
                    
                    total += abund;
                    containsGroup[groupNames[j]] = abund;
                    
                    if (basis == "otu") { if (abund > 0) { abund = 1;  } } //count presence in otu
                    abunds.push_back(abund);
                }
                ct.push_back(shared->getOTUName(i), abunds);
                otuContainsGroups.push_back(containsGroup);
            }
        }
        
        //print group taxonomies if given
        map<string, string> groupTaxonomies = biom->getGroupTaxonomies();
        
        if (groupTaxonomies.size() != 0) {
            //write taxonomy file
            map<string, string> variables;
            variables["[filename]"] = fileroot;
            variables["[tag]"] = label;
            string taxFilename = getOutputFileName("taxonomy",variables);
            outputNames.push_back(taxFilename); outputTypes["taxonomy"].push_back(taxFilename);
            ofstream outTax; util.openOutputFile(taxFilename, outTax);
            
            GroupMap* g = NULL; PhyloSummary taxSum(g, relabund, printlevel);
            
            //print group taxonomy if given
            for (map<string, string>::iterator it = groupTaxonomies.begin(); it!= groupTaxonomies.end(); it++) {
                outTax << it->first << '\t' << it->second << endl;
                taxSum.addSeqToTree(it->first, it->second);
            }
            outTax.close();
            
            //write taxonomy file
            variables["[tag2]"] = "";
            string taxSumFilename = getOutputFileName("taxsummary",variables);
            outputNames.push_back(taxSumFilename); outputTypes["taxsummary"].push_back(taxSumFilename);
            ofstream outTaxSum; util.openOutputFile(taxSumFilename, outTaxSum);
                
            //write tax.summary
            if (relabund)   {   taxSum.print(outTaxSum, relabund);     }
            else            {   taxSum.print(outTaxSum, output);       }
                
            outTaxSum.close();
        }
        
        
        //print consTaxonomy if given
        vector<Taxonomy> consTax = biom->getConsTaxonomies();
        if (consTax.size() != 0) {

            //write taxonomy file
            map<string, string> variables;
            variables["[filename]"] = fileroot;
            variables["[tag]"] = label;
            string taxFilename = getOutputFileName("constaxonomy",variables);
            outputNames.push_back(taxFilename); outputTypes["constaxonomy"].push_back(taxFilename);
            ofstream outTax; util.openOutputFile(taxFilename, outTax);
            
            outTax << "OTU\tSize\tTaxonomy\n";
            
            PhyloSummary consTaxSum(&ct, relabund, printlevel);
            
            for (int i = 0; i < consTax.size(); i++) {
                
                consTax[i].printConsTax(outTax);
                
                if (basis == "sequence")    {  consTaxSum.addSeqToTree(consTax[i].getName(), consTax[i].getConsTaxString());  }
                else                       { consTaxSum.addSeqToTree(consTax[i].getConsTaxString(), otuContainsGroups[i]); } //add otu
            }
            outTax.close();
            
            variables["[tag2]"] = "cons";
            string taxSumFilename = getOutputFileName("taxsummary",variables);
            outputNames.push_back(taxSumFilename); outputTypes["taxsummary"].push_back(taxSumFilename);
            ofstream outTaxSum; util.openOutputFile(taxSumFilename, outTaxSum);
        
            //write tax.summary
            if (relabund)   {   consTaxSum.print(outTaxSum, relabund);     }
            else            {   consTaxSum.print(outTaxSum, output);       }
            
            outTaxSum.close();
        }
        
        delete biom;
        
        m->mothurOut("\nIt took " + toString(time(NULL) - start) + " create mothur files from your biom file.\n\n");
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); } }
        
        string currentName = "";
        itTypes = outputTypes.find("shared");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setSharedFile(currentName); }
        }
        
        //set taxonomy file as new current taxonomyfile
        itTypes = outputTypes.find("taxonomy");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setTaxonomyFile(currentName); }
        }
        
        //set constaxonomy file as new current constaxonomyfile
        itTypes = outputTypes.find("constaxonomy");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setConsTaxonomyFile(currentName); }
        }
        
        m->mothurOutEndLine();
        m->mothurOut("Output File Names: \n"); 
        for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
        m->mothurOutEndLine();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "BiomInfoCommand", "execute");
        exit(1);
    }
}
//**********************************************************************************************************************

