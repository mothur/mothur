//
//  createdatabasecommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/28/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "createdatabasecommand.h"
#include "inputdata.h"

//**********************************************************************************************************************
vector<string> CreateDatabaseCommand::setParameters(){	
	try {
		CommandParameter pfasta("repfasta", "InputTypes", "", "", "none", "none", "none","database",false,false,true); parameters.push_back(pfasta);
		CommandParameter pname("repname", "InputTypes", "", "", "NameCount", "none", "none","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
		CommandParameter pconstaxonomy("constaxonomy", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(pconstaxonomy);
		CommandParameter plist("list", "InputTypes", "", "", "ListShared", "ListShared", "none","",false,false,true); parameters.push_back(plist);
        CommandParameter pshared("shared", "InputTypes", "", "", "ListShared", "ListShared", "none","",false,false,true); parameters.push_back(pshared);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "CreateDatabaseCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string CreateDatabaseCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The create.database command reads a list file or a shared file, *.cons.taxonomy, and optional *.rep.fasta, *.rep.names, groupfile, or count file and creates a database file.\n";
		helpString += "The create.database command parameters are repfasta, list, shared, repname, constaxonomy, group, count and label. List or shared and constaxonomy are required.\n";
        helpString += "The repfasta file is fasta file outputted by get.oturep(fasta=yourFastaFile, list=yourListfile, column=yourDistFile, name=yourNameFile).\n";
        helpString += "The repname file is the name file outputted by get.oturep(fasta=yourFastaFile, list=yourListfile, column=yourDistFile, name=yourNameFile).\n";
        helpString += "The count file is the count file outputted by get.oturep(fasta=yourFastaFile, list=yourListfile, column=yourDistFile, count=yourCountFile). If it includes group info, mothur will give you the abundance breakdown by group. \n";
        helpString += "The constaxonomy file is the taxonomy file outputted by classify.otu(list=yourListfile, taxonomy=yourTaxonomyFile, name=yourNameFile).\n";
        helpString += "The group file is optional and will just give you the abundance breakdown by group.\n";
        helpString += "The label parameter allows you to specify a label to be used from your listfile.\n";
        helpString += "NOTE: Make SURE the repfasta, repnames and contaxonomy are for the same label as the listfile.\n";
        helpString += "The create.database command should be in the following format: \n";
		helpString += "create.database(repfasta=yourFastaFileFromGetOTURep, repname=yourNameFileFromGetOTURep, contaxonomy=yourConTaxFileFromClassifyOTU, list=yourListFile) \n";	
		helpString += "Example: create.database(repfasta=final.an.0.03.rep.fasta, repname=final.an.0.03.rep.names, list=final.an.list, label=0.03, contaxonomy=final.an.0.03.cons.taxonomy) \n";
		helpString += "Note: No spaces between parameter labels (i.e. repfasta), '=' and parameters (i.e.yourFastaFileFromGetOTURep).\n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "CreateDatabaseCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string CreateDatabaseCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "database") {  pattern = "[filename],database"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "CreateDatabaseCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
CreateDatabaseCommand::CreateDatabaseCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["database"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "CreateDatabaseCommand", "CreateDatabaseCommand");
		exit(1);
	}
}

//**********************************************************************************************************************
CreateDatabaseCommand::CreateDatabaseCommand(string option)  {
	try{
		abort = false; calledHelp = false;   
        
		//allow user to run help
		if (option == "help") { 
			help(); abort = true; calledHelp = true;
		}else if(option == "citation") { citation(); abort = true; calledHelp = true;} 
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
            
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["database"] = tempOutNames;
            
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
				
				it = parameters.find("repname");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["repname"] = inputDir + it->second;		}
				}
				
				it = parameters.find("constaxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["constaxonomy"] = inputDir + it->second;		}
				}
				
				it = parameters.find("repfasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["repfasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
                
                it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
			}
            
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";		}
			
			//check for required parameters
			listfile = validParameter.validFile(parameters, "list");
			if (listfile == "not found") {	listfile = "";			}
			else if (listfile == "not open") { listfile = ""; abort = true; }	
			else { current->setListFile(listfile); }
            
            sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not found") {	sharedfile = "";			}
			else if (sharedfile == "not open") { sharedfile = ""; abort = true; }	
			else { current->setSharedFile(sharedfile); }
            
            if ((sharedfile == "") && (listfile == "")) { 
				//is there are current file available for either of these?
				//give priority to list, then shared
				listfile = current->getListFile(); 
				if (listfile != "") {  m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
				else { 
					sharedfile = current->getSharedFile(); 
					if (sharedfile != "") {  m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("[ERROR]: No valid current files. You must provide a shared or list file before you can use the create.database command."); m->mothurOutEndLine();
						abort = true;
					}
				}
			}
			else if ((sharedfile != "") && (listfile != "")) { m->mothurOut("When executing a create.database command you must enter ONLY ONE of the following: shared or list."); m->mothurOutEndLine(); abort = true; }
            
            if (sharedfile != "") { if (outputDir == "") { outputDir = util.hasPath(sharedfile); } }
            else { if (outputDir == "") { outputDir = util.hasPath(listfile); } }
			
			contaxonomyfile = validParameter.validFile(parameters, "constaxonomy");
			if (contaxonomyfile == "not found") {  //if there is a current list file, use it
               contaxonomyfile = "";  m->mothurOut("The constaxonomy parameter is required, aborting."); m->mothurOutEndLine(); abort = true;
			}
			else if (contaxonomyfile == "not open") { contaxonomyfile = ""; abort = true; }

            repfastafile = validParameter.validFile(parameters, "repfasta");
			if (repfastafile == "not found") { repfastafile = ""; }
			else if (repfastafile == "not open") { repfastafile = ""; abort = true; }

            repnamesfile = validParameter.validFile(parameters, "repname");
			if (repnamesfile == "not found") {  repnamesfile = "";  			}
			else if (repnamesfile == "not open") { repnamesfile = ""; abort = true; }
            
            if ((repnamesfile != "") && (repfastafile == "")) { m->mothurOut("[ERROR]: You must provide a repfasta file if you are using a repnames file."); m->mothurOutEndLine();
                abort = true; }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not found") {  countfile = "";  			}
			else if (countfile == "not open") { countfile = ""; abort = true; }
            
			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not open") { groupfile = ""; abort = true; }	
			else if (groupfile == "not found") { groupfile = ""; }
			else { current->setGroupFile(groupfile); }
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
            label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; m->mothurOut("You did not provide a label, I will use the first label in your listfile.\n");}
        }
	}
	catch(exception& e) {
		m->errorOut(e, "CreateDatabaseCommand", "CreateDatabaseCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int CreateDatabaseCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        //taxonomies holds the taxonomy info for each Otu
        //classifyOtuSizes holds the size info of each Otu to help with error checking
        vector<string> taxonomies;
        vector<string> otuLabels;
        vector<int> classifyOtuSizes = readTax(taxonomies, otuLabels);
        
        if (m->getControl_pressed()) { return 0; }
        
        vector<Sequence> seqs;
        vector<int> repOtusSizes;
        if (repfastafile != "") { repOtusSizes = readFasta(seqs); }
        
        if (m->getControl_pressed()) { return 0; }
        
        //names redundants to uniques. backwards to how we normally do it, but each bin is the list file will be a key entry in the map.
        map<string, string> repNames;
        map<string, int> nameMap;
        int numUniqueNamesFile = 0;
        CountTable ct;
        if (repnamesfile != "") {
            numUniqueNamesFile = util.readNames(repnamesfile, repNames, 1);
            //the repnames file does not have the same order as the list file bins so we need to sort and reassemble for the search below
            map<string, string> tempRepNames;
            for (map<string, string>::iterator it = repNames.begin(); it != repNames.end();) {
                string bin = it->first;
                vector<string> temp;
                util.splitAtChar(bin, temp, ',');
                sort(temp.begin(), temp.end());
                bin = "";
                for (int i = 0; i < temp.size()-1; i++) {
                    bin += temp[i] + ',';
                }
                bin += temp[temp.size()-1];
                tempRepNames[bin] = it->second;
                repNames.erase(it++);
            }
            repNames = tempRepNames;
        }else if (countfile != ""){
            ct.readTable(countfile, true, false);
            numUniqueNamesFile = ct.getNumUniqueSeqs();
            nameMap = ct.getNameMap();
        }
        
        if (m->getControl_pressed()) { return 0; }
        
        if (repfastafile != "") {
            
            //are there the same number of otus in the fasta and name files
            if (repOtusSizes.size() != numUniqueNamesFile) { m->mothurOut("[ERROR]: you have " + toString(numUniqueNamesFile) + " unique seqs in your repname file, but " + toString(repOtusSizes.size()) + " seqs in your repfasta file.  These should match.\n"); m->setControl_pressed(true); }
            
            //are there the same number of OTUs in the tax and fasta file
            if (classifyOtuSizes.size() != repOtusSizes.size()) { m->mothurOut("[ERROR]: you have " + toString(classifyOtuSizes.size()) + " taxonomies in your contaxonomy file, but " + toString(repOtusSizes.size()) + " seqs in your repfasta file.  These should match.\n"); m->setControl_pressed(true); }
            
            if (m->getControl_pressed()) { return 0; }
            
            //at this point we have the same number of OTUs. Are the sizes we have found so far accurate?
            for (int i = 0; i < classifyOtuSizes.size(); i++) {
                if (classifyOtuSizes[i] != repOtusSizes[i]) {
                    m->mothurOut("[ERROR]: OTU size info does not match for bin " + toString(i+1) + ". The contaxonomy file indicated the OTU represented " + toString(classifyOtuSizes[i]) + " sequences, but the repfasta file had " + toString(repOtusSizes[i]) + ".  These should match. Make sure you are using files for the same distance.\n"); m->setControl_pressed(true);
                }
            }
        }
        
        if (m->getControl_pressed()) { return 0; }
        
        
        map<string, string> variables; 
        if (listfile != "") {  variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(listfile)); }
        else { variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(sharedfile)); }
        string outputFileName = getOutputFileName("database", variables); 
        outputNames.push_back(outputFileName); outputTypes["database"].push_back(outputFileName);
        
        ofstream out;
        util.openOutputFile(outputFileName, out);
        
        string header = "OTUNumber\tAbundance";

        
        if (listfile != "") {
            //at this point we are fairly sure the repfasta, repnames and contaxonomy files match so lets proceed with the listfile
            ListVector* list = getList();
            
            if (otuLabels.size() != list->getNumBins()) { 
                m->mothurOut("[ERROR]: you have " + toString(otuLabels.size()) + " otus in your contaxonomy file, but your list file has " + toString(list->getNumBins()) + " otus. These should match. Make sure you are using files for the same distance.\n"); m->setControl_pressed(true);  }
            
            if (m->getControl_pressed()) { delete list; return 0; }
            
            GroupMap* groupmap = NULL;
            if (groupfile != "") {
                groupmap = new GroupMap(groupfile);
                groupmap->readMap();
            }
            
            if (m->getControl_pressed()) { delete list; if (groupfile != "") { delete groupmap; } return 0; }
            
            if (groupfile != "") { 
                header = "OTUNumber";
                for (int i = 0; i < groupmap->getNamesOfGroups().size(); i++) { header += '\t' + (groupmap->getNamesOfGroups())[i]; }
            }else if (countfile != "") {
                if (ct.hasGroupInfo()) {
                    header = "OTUNumber";
                    for (int i = 0; i < ct.getNamesOfGroups().size(); i++) { header += '\t' + (ct.getNamesOfGroups())[i]; }
                }
            }
            
            if (repfastafile != "") {  header += "\trepSeqName\trepSeq";  }
            header += "\tOTUConTaxonomy";
            out << header << endl;
            
            vector<string> binLabels = list->getLabels();
            for (int i = 0; i < list->getNumBins(); i++) {
                
                int index = findIndex(otuLabels, binLabels[i]);
                if (index == -1) {  m->mothurOut("[ERROR]: " + binLabels[i] + " is not in your constaxonomy file, aborting.\n"); m->setControl_pressed(true); }
                
                if (m->getControl_pressed()) { break; }
                
                out << otuLabels[index];
                
                vector<string> binNames;
                string bin = list->get(i);
                util.splitAtComma(bin, binNames);
                
                string seqRepName = "";
                int numSeqsRep = binNames.size();
                
                if (repnamesfile != "") {
                    sort(binNames.begin(), binNames.end());
                    bin = "";
                    for (int j = 0; j < binNames.size()-1; j++) {
                        bin += binNames[j] + ',';
                    }
                    bin += binNames[binNames.size()-1];
                    map<string, string>::iterator it = repNames.find(bin);
                    
                    if (it == repNames.end()) {
                        m->mothurOut("[ERROR: OTU " + otuLabels[index] + " is not in the repnames file. Make sure you are using files for the same distance.\n"); m->setControl_pressed(true);   break;
                    }else { seqRepName = it->second;  numSeqsRep = binNames.size(); }
                    
                    //sanity check
                    if (binNames.size() != classifyOtuSizes[index]) {
                        m->mothurOut("[ERROR: OTU " + otuLabels[index] + " contains " + toString(binNames.size()) + " sequence, but the rep and taxonomy files indicated this OTU should have " + toString(classifyOtuSizes[index]) + ". Make sure you are using files for the same distance.\n"); m->setControl_pressed(true);   break;
                    }
                }else if ((countfile != "") && (repfastafile != "")) {
                    //find rep sequence in bin
                    for (int j = 0; j < binNames.size(); j++) {
                        map<string, int>::iterator itNameMap = nameMap.find(binNames[j]); //if you are in the counttable you must be the rep. because get.oturep with a countfile only includes the rep sequences in the rep.count file.
                        if (itNameMap != nameMap.end()) {
                            seqRepName = itNameMap->first;
                            numSeqsRep = itNameMap->second;
                            j += binNames.size(); //exit loop
                        }
                    }
                    
                    if (seqRepName == "") {
                        m->mothurOut("[ERROR: OTU " + otuLabels[index] + " is not in the count file. Make sure you are using files for the same distance.\n"); m->setControl_pressed(true);   break;
                    }
                    
                    if (numSeqsRep != classifyOtuSizes[i]) {
                        m->mothurOut("[ERROR: OTU " + otuLabels[index] + " contains " + toString(numSeqsRep) + " sequence, but the rep and taxonomy files indicated this OTU should have " + toString(classifyOtuSizes[index]) + ". Make sure you are using files for the same distance.\n"); m->setControl_pressed(true);   break;
                    }
                }
                
                //output abundances
                if (groupfile != "") {
                    string groupAbunds = "";
                    map<string, int> counts;
                    //initialize counts to 0
                    for (int j = 0; j < groupmap->getNamesOfGroups().size(); j++) { counts[(groupmap->getNamesOfGroups())[j]] = 0; }
                    
                    //find abundances by group
                    bool error = false;
                    for (int j = 0; j < binNames.size(); j++) {
                        string group = groupmap->getGroup(binNames[j]);
                        if (group == "not found") {
                            m->mothurOut("[ERROR]: " + binNames[j] + " is not in your groupfile, please correct.\n");
                            error = true;
                        }else { counts[group]++; }
                    }
                    
                    //output counts
                    for (int j = 0; j < groupmap->getNamesOfGroups().size(); j++) { out << '\t' << counts[(groupmap->getNamesOfGroups())[j]];  }
                    
                    if (error) { m->setControl_pressed(true); }
                }else if ((countfile != "") && (repfastafile != "")) {
                    if (ct.hasGroupInfo()) {
                        vector<int> groupCounts = ct.getGroupCounts(seqRepName);
                        for (int j = 0; j < groupCounts.size(); j++) { out << '\t' << groupCounts[j];  }
                    }else { out << '\t' << numSeqsRep; }
                }else if ((countfile != "") && (repfastafile == "")) {
                    if (ct.hasGroupInfo()) {
                        vector<int> groupTotals; groupTotals.resize(ct.getNumGroups(), 0);
                        for (int j = 0; j < binNames.size(); j++) {
                            vector<int> groupCounts = ct.getGroupCounts(binNames[j]);
                            for (int k = 0; k < groupCounts.size(); k++) { groupTotals[k] += groupCounts[k]; }
                        }
                        for (int j = 0; j < groupTotals.size(); j++) { out << '\t' << groupTotals[j];  }
                    }else { out << '\t' << numSeqsRep; }
                }else { out << '\t' << numSeqsRep; }
                
                //output repSeq
                if (repfastafile != "") { out << '\t' << seqRepName << '\t' << seqs[index].getAligned() << '\t' << taxonomies[index] << endl; }
                else { out << '\t' << taxonomies[index] << endl; }
            }
            
            
            delete list;
            if (groupfile != "") { delete groupmap; }
           
        }else {
            SharedRAbundVectors* lookup = getShared();
            vector<string> namesOfGroups = lookup->getNamesGroups();
            
            header = "OTUNumber";
            for (int i = 0; i < namesOfGroups.size(); i++) { header += '\t' + namesOfGroups[i]; }
            if (repfastafile != "") {  header += "\trepSeqName\trepSeq";  }
            header += "\tOTUConTaxonomy";
            out << header << endl;
            
            vector<string> currentLabels = lookup->getOTUNames();
            for (int h = 0; h < lookup->getNumBins(); h++) {
                
                if (m->getControl_pressed()) { break; }
                
                int index = findIndex(otuLabels, currentLabels[h]);
                if (index == -1) {  m->mothurOut("[ERROR]: " + currentLabels[h] + " is not in your constaxonomy file, aborting.\n"); m->setControl_pressed(true); }
                
                if (m->getControl_pressed()) { break; }
                
                out << otuLabels[index];
                
                int totalAbund = 0;
                for (int i = 0; i < lookup->size(); i++) {
                    int abund = lookup->get(h, namesOfGroups[i]);
                    totalAbund += abund; 
                    out  << '\t' << abund;
                }
                
                //sanity check
                if (totalAbund != classifyOtuSizes[index]) {
                    m->mothurOut("[WARNING]: OTU " + currentLabels[h] + " contains " + toString(totalAbund) + " sequence, but the rep and taxonomy files indicated this OTU should have " + toString(classifyOtuSizes[index]) + ". Make sure you are using files for the same distance.\n"); //m->setControl_pressed(true);   break;
                }
                
                //output repSeq
                if (repfastafile != "") { out << '\t' << seqs[index].getName() << '\t' << seqs[index].getAligned() << '\t' << taxonomies[index] << endl; }
                else { out << '\t' << taxonomies[index] << endl; }
            }
        }
        out.close();
        if (m->getControl_pressed()) { util.mothurRemove(outputFileName); return 0; }
        
        m->mothurOut("\nOutput File Names: \n"); 
		m->mothurOut(outputFileName); m->mothurOutEndLine();	
		m->mothurOutEndLine();
        
        return 0;
        
    }
	catch(exception& e) {
		m->errorOut(e, "CreateDatabaseCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int CreateDatabaseCommand::findIndex(vector<string>& otuLabels, string label){
	try {
        int index = -1;
        for (int i = 0; i < otuLabels.size(); i++) {
            if (util.isLabelEquivalent(otuLabels[i],label)) { index = i; break; }
        }
		return index;
    }
	catch(exception& e) {
		m->errorOut(e, "CreateDatabaseCommand", "findIndex");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<int> CreateDatabaseCommand::readTax(vector<string>& taxonomies, vector<string>& otuLabels){
	try {
		
        vector<int> sizes; 
        
        ifstream in;
        util.openInputFile(contaxonomyfile, in);
        
        //read headers
        util.getline(in);
        
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            string otu = ""; string tax = "unknown";
            int size = 0;
            
            in >> otu >> size; util.gobble(in);
            tax = util.getline(in); util.gobble(in);
            
            sizes.push_back(size);
            taxonomies.push_back(tax);
            otuLabels.push_back(otu);
        }
        in.close();
        
        return sizes;
    }
	catch(exception& e) {
		m->errorOut(e, "CreateDatabaseCommand", "readTax");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<int> CreateDatabaseCommand::readFasta(vector<Sequence>& seqs){
	try {
		
        vector<int> sizes; 
        
        ifstream in;
        util.openInputFile(repfastafile, in);
        
        set<int> sanity;
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            string binInfo;
            Sequence seq(in, binInfo, true);  util.gobble(in);
            
            //the binInfo should look like - binNumber|size ie. 1|200 if it is binNumber|size|group then the user gave us the wrong repfasta file
            vector<string> info;
            util.splitAtChar(binInfo, info, '|');
            //if (info.size() != 2) { m->mothurOut("[ERROR]: your repfasta file is not the right format.  The create database command is designed to be used with the output from get.oturep.  When running get.oturep you can not use a group file, because mothur is only expecting one representative sequence per OTU and when you use a group file with get.oturep a representative is found for each group.\n");  m->setControl_pressed(true); break;}
            
            int size = 0;
            util.mothurConvert(info[1], size);
            
            int binNumber = 0;
            string temp = "";
            for (int i = 0; i < info[0].size(); i++) { if (isspace(info[0][i])) {;}else{temp +=info[0][i]; } }
            util.mothurConvert(util.getSimpleLabel(temp), binNumber);
            set<int>::iterator it = sanity.find(binNumber);
            if (it != sanity.end()) {
                m->mothurOut("[ERROR]: your repfasta file is not the right format.  The create database command is designed to be used with the output from get.oturep.  When running get.oturep you can not use a group file, because mothur is only expecting one representative sequence per OTU and when you use a group file with get.oturep a representative is found for each group.\n");  m->setControl_pressed(true); break;
            }else { sanity.insert(binNumber); }
            
            sizes.push_back(size);
            seqs.push_back(seq);
        }
        in.close();
        
        return sizes;
    }
	catch(exception& e) {
		m->errorOut(e, "CreateDatabaseCommand", "readFasta");
		exit(1);
	}
}
//**********************************************************************************************************************
ListVector* CreateDatabaseCommand::getList(){
	try {
		InputData* input = new InputData(listfile, "list", nullVector);
		ListVector* list = input->getListVector();
		string lastLabel = list->getLabel();
		
		if (label == "") { label = lastLabel; delete input; return list; }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> labels; labels.insert(label);
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((list != NULL) && (userLabels.size() != 0)) {
			if (m->getControl_pressed()) {  delete input; return list;  }
			
			if(labels.count(list->getLabel()) == 1){
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				break;
			}
			
			if ((util.anyLabelsToProcess(list->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = list->getLabel();
				
				delete list;
				list = input->getListVector(lastLabel);
				
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				
				//restore real lastlabel to save below
				list->setLabel(saveLabel);
				break;
			}
			
			lastLabel = list->getLabel();			
			
			//get next line to process
			//prevent memory leak
			delete list; 
			list = input->getListVector();
		}
		
		
		if (m->getControl_pressed()) { delete input; return list;  }
		
		//output error messages about any remaining user labels
		bool needToRun = false;
		for (set<string>::iterator it = userLabels.begin(); it != userLabels.end(); it++) {
			m->mothurOut("Your file does not include the label " + *it); 
            if (processedLabels.count(lastLabel) != 1)  { m->mothurOut(". I will use " + lastLabel + ".\n"); needToRun = true;  }
			else                                        { m->mothurOut(". Please refer to " + lastLabel + ".\n");               }
		}
		
		//run last label if you need to
		if (needToRun )  {
			delete list;
			list = input->getListVector(lastLabel);
		}	
		
		delete input;

        return list;
    }
	catch(exception& e) {
		m->errorOut(e, "CreateDatabaseCommand", "getList");
		exit(1);
	}
}
//**********************************************************************************************************************
SharedRAbundVectors* CreateDatabaseCommand::getShared(){
	try {
		InputData input(sharedfile, "sharedfile", nullVector);
		SharedRAbundVectors* lookup = input.getSharedRAbundVectors();
		string lastLabel = lookup->getLabel();
		
		if (label == "") { label = lastLabel; return lookup; }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> labels; labels.insert(label);
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup != NULL) && (userLabels.size() != 0)) {
			if (m->getControl_pressed()) {  return lookup;  }
			
			if(labels.count(lookup->getLabel()) == 1){
				processedLabels.insert(lookup->getLabel()); userLabels.erase(lookup->getLabel());
				break;
			}
			
			if ((util.anyLabelsToProcess(lookup->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup->getLabel();
				
                delete lookup;
				lookup = input.getSharedRAbundVectors(lastLabel);
				
				processedLabels.insert(lookup->getLabel()); userLabels.erase(lookup->getLabel());
				
				//restore real lastlabel to save below
				lookup->setLabels(saveLabel);
				break;
			}
			
			lastLabel = lookup->getLabel();
			
			//get next line to process
			//prevent memory leak
			delete lookup;
			lookup = input.getSharedRAbundVectors();
		}
		
		
		if (m->getControl_pressed()) { return lookup;  }
		
		//output error messages about any remaining user labels
		bool needToRun = false;
		for (set<string>::iterator it = userLabels.begin(); it != userLabels.end(); it++) {
			m->mothurOut("Your file does not include the label " + *it); 
            if (processedLabels.count(lastLabel) != 1)  { m->mothurOut(". I will use " + lastLabel + ".\n"); needToRun = true;  }
			else                                        { m->mothurOut(". Please refer to " + lastLabel + ".\n");               }
		}
		
		//run last label if you need to
		if (needToRun )  {
			delete lookup;
			lookup = input.getSharedRAbundVectors(lastLabel);
		}	
        
        return lookup;
    }
	catch(exception& e) {
		m->errorOut(e, "CreateDatabaseCommand", "getShared");
		exit(1);
	}
}

//**********************************************************************************************************************


