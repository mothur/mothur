/*
 *  clustersplitcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 5/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "clustersplitcommand.h"
#include "systemcommand.h"
#include "sensspeccommand.h"
#include "mcc.hpp"
#include "sensitivity.hpp"
#include "specificity.hpp"
#include "fdr.hpp"
#include "npv.hpp"
#include "ppv.hpp"
#include "f1score.hpp"
#include "tp.hpp"
#include "fp.hpp"
#include "fpfn.hpp"
#include "tptn.hpp"
#include "tn.hpp"
#include "fn.hpp"
#include "accuracy.hpp"

//**********************************************************************************************************************
vector<string> ClusterSplitCommand::setParameters(){	
	try {
        CommandParameter pfile("file", "InputTypes", "", "", "PhylipColumnFasta", "PhylipColumnFasta", "none","",false,false,true); parameters.push_back(pfile);
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "none", "FastaTaxName","",false,false,true); parameters.push_back(ptaxonomy);
		CommandParameter pphylip("phylip", "InputTypes", "", "", "PhylipColumnFasta", "PhylipColumnFasta", "none","list",false,false,true); parameters.push_back(pphylip);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "PhylipColumnFasta", "PhylipColumnFasta", "FastaTaxName","list",false,false,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "ColumnName-FastaTaxName","rabund-sabund",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount", "none", "","",false,false,true); parameters.push_back(pcount);
		CommandParameter pcolumn("column", "InputTypes", "", "", "PhylipColumnFasta", "PhylipColumnFasta", "ColumnName","list",false,false,true); parameters.push_back(pcolumn);
		CommandParameter ptaxlevel("taxlevel", "Number", "", "3", "", "", "","",false,false,true); parameters.push_back(ptaxlevel);
		CommandParameter psplitmethod("splitmethod", "Multiple", "classify-fasta-distance", "distance", "", "", "","",false,false,true); parameters.push_back(psplitmethod);
		CommandParameter plarge("large", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(plarge);
		CommandParameter pshowabund("showabund", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pshowabund);
        CommandParameter prunspenspec("runsensspec", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(prunspenspec);
        CommandParameter pcluster("cluster", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pcluster);
		CommandParameter ptiming("timing", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(ptiming);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pcutoff("cutoff", "Number", "", "0.03", "", "", "","",false,false,true); parameters.push_back(pcutoff);
        CommandParameter pmetriccutoff("delta", "Number", "", "0.0001", "", "", "","",false,false,true); parameters.push_back(pmetriccutoff);
        CommandParameter piters("iters", "Number", "", "100", "", "", "","",false,false,true); parameters.push_back(piters);
        CommandParameter pinitialize("initialize", "Multiple", "oneotu-singleton", "singleton", "", "", "","",false,false,true); parameters.push_back(pinitialize);
        CommandParameter pprecision("precision", "Number", "", "100", "", "", "","",false,false); parameters.push_back(pprecision);
        CommandParameter pmethod("method", "Multiple", "furthest-nearest-average-weighted-agc-dgc-opti", "opti", "", "", "","",false,false,true); parameters.push_back(pmethod);
        CommandParameter pmetric("metric", "Multiple", "mcc-sens-spec-tptn-fpfn-tp-tn-fp-fn-f1score-accuracy-ppv-npv-fdr", "mcc", "", "", "","",false,false,true); parameters.push_back(pmetric);
       CommandParameter pdist("dist", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pdist);
        CommandParameter pislist("islist", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pislist);
        CommandParameter pclassic("classic", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pclassic);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
			
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ClusterSplitCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The cluster.split command parameter options are file, fasta, phylip, column, name, count, cutoff, precision, method, splitmethod, taxonomy, taxlevel, showabund, timing, large, cluster, iters, delta, initialize, dist, processors, runsensspec. Fasta or Phylip or column and name are required.\n";
		helpString += "The cluster.split command can split your files in 3 ways. Splitting by distance file, by classification, or by classification also using a fasta file. \n";
		helpString += "For the distance file method, you need only provide your distance file and mothur will split the file into distinct groups. \n";
		helpString += "For the classification method, you need to provide your distance file and taxonomy file, and set the splitmethod to classify.  \n";
		helpString += "You will also need to set the taxlevel you want to split by. mothur will split the sequences into distinct taxonomy groups, and split the distance file based on those groups. \n";
		helpString += "For the classification method using a fasta file, you need to provide your fasta file, names file and taxonomy file.  \n";
		helpString += "You will also need to set the taxlevel you want to split by. mothur will split the sequence into distinct taxonomy groups, and create distance files for each grouping. \n";
        helpString += "The file option allows you to enter your file containing your list of column and names/count files as well as the singleton file.  This file is mothur generated, when you run cluster.split() with the cluster=f parameter.  This can be helpful when you have a large dataset that you may be able to use all your processors for the splitting step, but have to reduce them for the cluster step due to RAM constraints. For example: cluster.split(fasta=yourFasta, taxonomy=yourTax, count=yourCount, taxlevel=3, cluster=f, processors=8) then cluster.split(file=yourFile, processors=4).  This allows your to maximize your processors during the splitting step.  Also, if you are unsure if the cluster step will have RAM issue with multiple processors, you can avoid running the first part of the command multiple times.\n";
		helpString += "The phylip and column parameter allow you to enter your distance file. \n";
		helpString += "The fasta parameter allows you to enter your aligned fasta file. \n";
		helpString += "The name parameter allows you to enter your name file. \n";
        helpString += "The count parameter allows you to enter your count file. \n A count or name file is required if your distance file is in column format";
        helpString += "The cluster parameter allows you to indicate whether you want to run the clustering or just split the distance matrix, default=t";
        helpString += "The dist parameter allows you to indicate whether you want a column formatted distance matrix outputted along with the list file. Default=F.";
		helpString += "The cutoff parameter allow you to set the distance you want to cluster to, default is 0.03. \n";
		helpString += "The precision parameter allows you specify the precision of the precision of the distances outputted, default=100, meaning 2 decimal places. \n";
        helpString += "The iters parameter allow you to set the maxiters for the opticluster method. \n";
        helpString += "The metric parameter allows to select the metric in the opticluster method. Options are Matthews correlation coefficient (mcc), sensitivity (sens), specificity (spec), true positives + true negatives (tptn), false positives + false negatives (fpfn), true positives (tp), true negative (tn), false positive (fp), false negative (fn), f1score (f1score), accuracy (accuracy), positive predictive value (ppv), negative predictive value (npv), false discovery rate (fdr). Default=mcc.\n";
        helpString += "The delta parameter allows to set the stable value for the metric in the opticluster method. Default=0.0001\n";
        helpString += "The initialize parameter allows to select the initial randomization for the opticluster method. Options are singleton, meaning each sequence is randomly assigned to its own OTU, or oneotu meaning all sequences are assigned to one otu. Default=singleton.\n";
        helpString += "The runsensspec parameter allows to run the sens.spec command on the completed list file. Default=true.\n";
		helpString += "The method parameter allows you to enter your clustering mothod. Options are furthest, nearest, average, weighted, agc, dgc and opti. Default=opti.  The agc and dgc methods require a fasta file.";
		helpString += "The splitmethod parameter allows you to specify how you want to split your distance file before you cluster, default=distance, options distance, classify or fasta. \n";
		helpString += "The taxonomy parameter allows you to enter the taxonomy file for your sequences, this is only valid if you are using splitmethod=classify. Be sure your taxonomy file does not include the probability scores. \n";
		helpString += "The taxlevel parameter allows you to specify the taxonomy level you want to use to split the distance file, default=3, meaning use the first taxon in each list. \n";
		helpString += "The large parameter allows you to indicate that your distance matrix is too large to fit in RAM.  The default value is false.\n";
        helpString += "The classic parameter allows you to indicate that you want to run your files with cluster.classic.  It is only valid with splitmethod=fasta. Default=f.\n";
        helpString += "The processors parameter allows you to specify the number of processors to use. The default is 1.\n";
		helpString += "The cluster.split command should be in the following format: \n";
		helpString += "cluster.split(column=youDistanceFile, name=yourNameFile, method=yourMethod, cutoff=yourCutoff, precision=yourPrecision, splitmethod=yourSplitmethod, taxonomy=yourTaxonomyfile, taxlevel=yourtaxlevel) \n";
		helpString += "Example: cluster.split(column=abrecovery.dist, name=abrecovery.names, method=opti, cutoff=0.10, precision=1000, splitmethod=classify, taxonomy=abrecovery.silva.slv.taxonomy, taxlevel=5) \n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ClusterSplitCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "list") {  pattern = "[filename],[clustertag],list-[filename],[clustertag],[tag2],list"; } 
        else if (type == "rabund") {  pattern = "[filename],[clustertag],rabund"; } 
        else if (type == "sabund") {  pattern = "[filename],[clustertag],sabund"; }
        else if (type == "sensspec") {  pattern = "[filename],[clustertag],sensspec"; }
        else if (type == "column") {  pattern = "[filename],dist"; }
        else if (type == "file")   {  pattern = "[filename],file"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterSplitCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
string ClusterSplitCommand::getCommonQuestions(){
    try {
        vector<string> questions, issues, qanswers, ianswers, howtos, hanswers;
        
        string issue = "Cluster.split crashes after merging individual list files. What do I do?"; issues.push_back(issue);
        string ianswer = "\tAfter merging the split list files, mothur runs the sens.spec command on the entire dataset. The entire dataset's distance matrix may be too large to fit in memory, which causes the crash. You can skip this step by setting the runsensspec parameter to false. Skipping the sens.spec analysis does not effect the OTU assignment, and you can run the sens.spec analysis separately using the sens.spec command. \n"; ianswers.push_back(ianswer);
        
        issue = "Cluster.split crashes while reading the split distance matrices. What should I do?"; issues.push_back(issue);
        ianswer = "\tThe command is crashing because the distance matrices are too large to fit into memory. Why do I have such a large distance matrix? This is most often caused by poor overlap of your reads. When reads have poor overlap, it greatly increases your error rate. Also, sequences that should cluster together don't because the errors appear to be genetic differences when in fact they are not. The quality of the data you are processing can not be overstressed. Error filled reads produce error filled results. To take a step back, if you look through our MiSeq SOP, you’ll see that we go to great pains to only work with the unique sequences to limit the number of sequences we have to align, screen for chimeras, classify, etc. We all know that 20 million reads will never make it through the pipeline without setting your computer on fire. Returning to the question at hand, you can imagine that if the reads do not fully overlap then any error in the 5’ end of the first read will be uncorrected by the 3’ end of the second read. If we assume for now that the errors are random, then every error will generate a new unique sequence. Granted, this happens less than 1% of the time, but multiply that by 20 million reads at whatever length you choose and you’ve got a big number. Viola, a bunch of unique reads and a ginormous distance matrix. \n"; ianswers.push_back(ianswer);
        
        string howto = "How do I cluster my sequences into OTUs at distance 0.03?"; howtos.push_back(howto);
        string hanswer = "\tBy default the cluster.split command will use the opti method to cluster to 0.03. To find OTUs at a different distance set the cutoff parameter. ie. cutoff=0.01 will assemble OTUs for distance 0.01.\n"; hanswers.push_back(hanswer);
        
        string commonQuestions = util.getFormattedHelp(questions, qanswers, issues, ianswers, howtos, hanswers);
        
        return commonQuestions;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterSplitCommand", "getCommonQuestions");
        exit(1);
    }
}
//**********************************************************************************************************************
ClusterSplitCommand::ClusterSplitCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["list"] = tempOutNames;
		outputTypes["rabund"] = tempOutNames;
		outputTypes["sabund"] = tempOutNames;
		outputTypes["column"] = tempOutNames;
        outputTypes["name"] = tempOutNames;
        outputTypes["file"] = tempOutNames;
        outputTypes["sensspec"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "ClusterSplitCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
//This function checks to make sure the cluster command has no errors and then clusters based on the method chosen.
ClusterSplitCommand::ClusterSplitCommand(string option)  {
	try{
		abort = false; calledHelp = false;   
		format = "";
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter("cluster.split");
		
			//check to make sure all parameters are valid for command
			map<string,string>::iterator it;
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {
					abort = true;
				}
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["list"] = tempOutNames;
			outputTypes["rabund"] = tempOutNames;
			outputTypes["sabund"] = tempOutNames;
			outputTypes["column"] = tempOutNames;
            outputTypes["file"] = tempOutNames;
            outputTypes["sensspec"] = tempOutNames;
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";		}
			
				//if the user changes the input directory command factory will send this info to us in the output parameter 
			inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("phylip");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["phylip"] = inputDir + it->second;		}
				}
				
				it = parameters.find("column");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["column"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
				
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
				
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
                
                it = parameters.find("file");
				//user has given a template file
				if(it != parameters.end()){
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["file"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			file = validParameter.validFile(parameters, "file");
			if (file == "not open") { file = ""; abort = true; }
			else if (file == "not found") { file = ""; }
            else { distfile = file; }
            
            phylipfile = validParameter.validFile(parameters, "phylip");
			if (phylipfile == "not open") { abort = true; }
			else if (phylipfile == "not found") { phylipfile = ""; }
			else {  distfile = phylipfile;  format = "phylip"; 	current->setPhylipFile(phylipfile); }
			
			columnfile = validParameter.validFile(parameters, "column");
			if (columnfile == "not open") { abort = true; }	
			else if (columnfile == "not found") { columnfile = ""; }
			else {  distfile = columnfile; format = "column";	current->setColumnFile(columnfile); }
			
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { abort = true; namefile = "";}	
			else if (namefile == "not found") { namefile = "";  }
			else { current->setNameFile(namefile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { abort = true; countfile = "";}	
			else if (countfile == "not found") { countfile = "";  }
			else { current->setCountFile(countfile); }
			
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not open") { abort = true; }	
			else if (fastafile == "not found") { fastafile = ""; }
			else { distfile = fastafile;  splitmethod = "fasta";  current->setFastaFile(fastafile); }
			
			taxFile = validParameter.validFile(parameters, "taxonomy");
			if (taxFile == "not open") { taxFile = ""; abort = true; }	
			else if (taxFile == "not found") { taxFile = ""; }
			else {
                current->setTaxonomyFile(taxFile);
                if (splitmethod != "fasta")         { splitmethod = "classify";     }
            }
			
			if ((phylipfile == "") && (columnfile == "") && (fastafile == "") && (file == "")) {
				//is there are current file available for either of these?
				//give priority to column, then phylip, then fasta
				columnfile = current->getColumnFile(); 
				if (columnfile != "") {  format = "column"; m->mothurOut("Using " + columnfile + " as input file for the column parameter."); m->mothurOutEndLine(); }
				else { 
					phylipfile = current->getPhylipFile(); 
					if (phylipfile != "") {  format = "phylip"; m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter."); m->mothurOutEndLine(); }
					else { 
						fastafile = current->getFastaFile(); 
						if (fastafile != "") {   m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
						else { 
							m->mothurOut("No valid current files. When executing a cluster.split command you must enter a file, phylip or a column or fastafile."); m->mothurOutEndLine();
							abort = true; 
						}
					}
				}
			}
			else if ((phylipfile != "") && (columnfile != "") && (fastafile != "") && (file != "")) { m->mothurOut("When executing a cluster.split command you must enter ONLY ONE of the following: file, fasta, phylip or column."); m->mothurOutEndLine(); abort = true; }
            
            if ((countfile != "") && (namefile != "")) { m->mothurOut("When executing a cluster.split command you must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
            
			if (columnfile != "") {
				if ((namefile == "") && (countfile == "")) { 
					namefile = current->getNameFile(); 
					if (namefile != "") {  m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
					else { 
						countfile = current->getCountFile();
                        if (countfile != "") {  m->mothurOut("Using " + countfile + " as input file for the count parameter."); m->mothurOutEndLine(); }
                        else { 
                            m->mothurOut("You need to provide a namefile or countfile if you are going to use the column format."); m->mothurOutEndLine(); 
                            abort = true; 
                        }	
					}	
				}
			}
			
			if (fastafile != "") {
				if (taxFile == "") { 
					taxFile = current->getTaxonomyFile(); 
					if (taxFile != "") {  m->mothurOut("Using " + taxFile + " as input file for the taxonomy parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("You need to provide a taxonomy file if you are if you are using a fasta file to generate the split."); m->mothurOutEndLine(); 
						abort = true; 
					}	
				}
				
				if ((namefile == "") && (countfile == "")) { 
					namefile = current->getNameFile(); 
					if (namefile != "") {  m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
					else { 
						countfile = current->getCountFile();
                        if (countfile != "") {  m->mothurOut("Using " + countfile + " as input file for the count parameter."); m->mothurOutEndLine(); }
                        else { 
                            m->mothurOut("You need to provide a namefile or countfile if you are going to use the fasta file to generate the split."); m->mothurOutEndLine(); 
                            abort = true; 
                        }	
					}	
				}
			}
					
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			//get user cutoff and precision or use defaults
			string temp;
			temp = validParameter.valid(parameters, "precision");
			if (temp == "not found") { temp = "100"; }
			//saves precision legnth for formatting below
			length = temp.length();
			util.mothurConvert(temp, precision); 
			
			temp = validParameter.valid(parameters, "large");			if (temp == "not found") { temp = "F"; }
			large = util.isTrue(temp);
            
			temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
			
			temp = validParameter.valid(parameters, "splitmethod");
			if ((splitmethod != "fasta") && (splitmethod != "classify")) {
				if (temp == "not found")  { splitmethod = "distance"; }
				else {  splitmethod = temp; }
			}
			
            temp = validParameter.valid(parameters, "classic");			if (temp == "not found") { temp = "F"; }
			classic = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "runsensspec");			if (temp == "not found") { temp = "T"; }
            runsensSpec = util.isTrue(temp);
            
            //not using file option and don't have fasta method with classic
            if (((splitmethod != "fasta") && classic) && (file == "")) { m->mothurOut("[ERROR]: splitmethod must be fasta to use cluster.classic, or you must use the file option.\n"); abort=true; }
			
			temp = validParameter.valid(parameters, "taxlevel");		if (temp == "not found")  { temp = "3"; }
			util.mothurConvert(temp, taxLevelCutoff);
            
            temp = validParameter.valid(parameters, "iters");		if (temp == "not found")  { temp = "100"; }
            util.mothurConvert(temp, maxIters);
            
            temp = validParameter.valid(parameters, "delta");		if (temp == "not found")  { temp = "0.0001"; }
            util.mothurConvert(temp, stableMetric);
			
            metricName = validParameter.valid(parameters, "metric");		if (metricName == "not found") { metricName = "mcc"; }
            
            if ((metricName == "mcc") || (metricName == "sens") || (metricName == "spec") || (metricName == "tptn") || (metricName == "tp") || (metricName == "tn") || (metricName == "fp") || (metricName == "fn") || (metricName == "f1score") || (metricName == "accuracy") || (metricName == "ppv") || (metricName == "npv") || (metricName == "fdr") || (metricName == "fpfn") ){ }
            else { m->mothurOut("[ERROR]: Not a valid metric.  Valid metrics are mcc, sens, spec, tp, tn, fp, fn, tptn, fpfn, f1score, accuracy, ppv, npv, fdr."); m->mothurOutEndLine(); abort = true; }
            
            initialize = validParameter.valid(parameters, "initialize");		if (initialize == "not found") { initialize = "singleton"; }
            
            if ((initialize == "singleton") || (initialize == "oneotu")){ }
            else { m->mothurOut("[ERROR]: Not a valid initialization.  Valid initializations are singleton and oneotu."); m->mothurOutEndLine(); abort = true; }

			method = validParameter.valid(parameters, "method");		if (method == "not found") { method = "opti";  }
			
            if ((method == "furthest") || (method == "nearest") || (method == "average") || (method == "weighted") || (method == "agc") || (method == "dgc") || (method == "opti")) { }
            else { m->mothurOut("[ERROR]: Not a valid clustering method.  Valid clustering algorithms are furthest, nearest, average, weighted, agc, dgc and opti."); m->mothurOutEndLine(); abort = true; }
            
            if ((method == "agc") || (method == "dgc")) {
                if (fastafile == "") { m->mothurOut("[ERROR]: You must provide a fasta file when using the agc or dgc clustering methods, aborting\n."); abort = true;}
                if (classic) { m->mothurOut("[ERROR]: You cannot use cluster.classic with the agc or dgc clustering methods, aborting\n."); abort = true; }
            }
            
            cutoffNotSet = false;
            temp = validParameter.valid(parameters, "cutoff");
            if (temp == "not found") { cutoffNotSet = true; if ((method == "opti") || (method == "agc") || (method == "dgc")) { temp = "0.03"; }else { temp = "0.15"; } }
            util.mothurConvert(temp, cutoff);
            
			if ((splitmethod == "distance") || (splitmethod == "classify") || (splitmethod == "fasta")) { }
			else { m->mothurOut("[ERROR]: " + splitmethod + " is not a valid splitting method.  Valid splitting algorithms are distance, classify or fasta."); m->mothurOutEndLine(); abort = true; }
			
			if ((splitmethod == "classify") && (taxFile == "")) {  m->mothurOut("[ERROR]: You need to provide a taxonomy file if you are going to use the classify splitmethod."); m->mothurOutEndLine(); abort = true;  }

			temp = validParameter.valid(parameters, "showabund");
			if (temp == "not found") { temp = "T"; }
            showabund = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "cluster");  if (temp == "not found") { temp = "T"; }
            runCluster = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "islist");  if (temp == "not found") { temp = "F"; }
            isList = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "dist");  if (temp == "not found") { temp = "F"; }
            makeDist = util.isTrue(temp);
            if (method == "opti") { makeDist = true; }
            if (((phylipfile != "") || (columnfile != "")) && (method == "opti")) { makeDist = false; }
            
            if (((phylipfile != "") || (columnfile != "")) && makeDist) { m->mothurOut("[ERROR]: You already provided a distance matrix. Mothur will ignore the dist parameter.\n"); makeDist = false; }
            if (classic && makeDist) { m->mothurOut("[ERROR]: You cannot use the dist parameter with the classic parameter. Mothur will ignore the dist parameter.\n"); makeDist = false; }

			timing = validParameter.valid(parameters, "timing");
			if (timing == "not found") { timing = "F"; }
			
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "ClusterSplitCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

int ClusterSplitCommand::execute(){
	try {
	
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        time_t estart;
        vector<string> listFileNames;
        vector< map<string, string> > distName;
        set<string> labels;
        string singletonName = "";
        
        double saveCutoff = cutoff;

        if (file != "") {
            deleteFiles = false; estart = time(NULL);
            singletonName = readFile(distName);

            if (isList) {
                
                //set list file as new current listfile
                string currentName = "";
                itTypes = outputTypes.find("list");
                if (itTypes != outputTypes.end()) {
                    if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setListFile(currentName); }
                }
                
                m->mothurOutEndLine();
                m->mothurOut("Output File Names: "); m->mothurOutEndLine();
                for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
                m->mothurOutEndLine();
                
                return 0;
            }
                    
        }else {
		          
            //****************** file prep work ******************************//
            
                //if user gave a phylip file convert to column file
                if (format == "phylip") {
                    estart = time(NULL);
                    m->mothurOut("Converting to column format..."); m->mothurOutEndLine();
                    
                    ReadCluster* convert = new ReadCluster(distfile, cutoff, outputDir, false);
                    
                    NameAssignment* nameMap = NULL;
                    convert->setFormat("phylip");
                    convert->read(nameMap);
                    
                    if (m->getControl_pressed()) {  delete convert;  return 0;  }
                    
                    distfile = convert->getOutputFile();
                    
                    //if no names file given with phylip file, create it
                    ListVector* listToMakeNameFile =  convert->getListVector();
                    if ((namefile == "") && (countfile == "")) {  //you need to make a namefile for split matrix
                        ofstream out;
                        namefile = phylipfile + ".names";
                        util.openOutputFile(namefile, out);
                        for (int i = 0; i < listToMakeNameFile->getNumBins(); i++) {
                            string bin = listToMakeNameFile->get(i);
                            out << bin << '\t' << bin << endl;
                        }
                        out.close();
                    }
                    delete listToMakeNameFile;
                    delete convert;
                    
                    m->mothurOut("It took " + toString(time(NULL) - estart) + " seconds to convert the distance file."); m->mothurOutEndLine();
                }
                if (m->getControl_pressed()) { return 0; }
                
                estart = time(NULL);
                m->mothurOut("Splitting the file..."); m->mothurOutEndLine();
                current->setMothurCalling(true);
            
                //split matrix into non-overlapping groups
                SplitMatrix* split;
                if (splitmethod == "distance")			{	split = new SplitMatrix(distfile, namefile, countfile, taxFile, cutoff, splitmethod, large);							}
                else if (splitmethod == "classify")		{	split = new SplitMatrix(distfile, namefile, countfile, taxFile, taxLevelCutoff, splitmethod, large);					}
                else if (splitmethod == "fasta")		{
                    if ((method == "agc") || (method == "dgc")) {
                        if (!findVsearch()) { m->mothurOut("[ERROR] cannot find vsearch executable, aborting.\n"); return 0; }
                        split = new SplitMatrix(fastafile, namefile, countfile, taxFile, taxLevelCutoff, cutoff, "vsearch", processors, classic, outputDir, "fasta");
                    }else{
                        split = new SplitMatrix(fastafile, namefile, countfile, taxFile, taxLevelCutoff, cutoff, splitmethod, processors, classic, outputDir, "distance");
                    }
                }
                else { m->mothurOut("Not a valid splitting method.  Valid splitting algorithms are distance, classify or fasta."); m->mothurOutEndLine(); return 0;		}
                split->split();
                if (fastafile != "") {  current->setFastaFile(fastafile);  }

                if (m->getControl_pressed()) { delete split; return 0; }
                
                singletonName = split->getSingletonNames();
                numSingletons = split->getNumSingleton();
                distName = split->getDistanceFiles();  //returns map of distance files -> namefile sorted by distance file size
                delete split;
                current->setMothurCalling(false);
                if (m->getDebug()) { m->mothurOut("[DEBUG]: distName.size() = " + toString(distName.size()) + ".\n"); }
                
                //output a merged distance file
                if (makeDist)		{ createMergedDistanceFile(distName); }
				
                if (m->getControl_pressed()) { return 0; }
                
                m->mothurOut("It took " + toString(time(NULL) - estart) + " seconds to split the distance file."); m->mothurOutEndLine();
                estart = time(NULL);

                if (!runCluster) {
                    string filename = printFile(singletonName, distName);
                    
                    m->mothurOutEndLine();
                    m->mothurOut("Output File Names:\n\n"); m->mothurOut(filename); m->mothurOutEndLine();
                    for (int i = 0; i < distName.size(); i++) {	m->mothurOut(distName[i].begin()->first); m->mothurOutEndLine(); m->mothurOut(distName[i].begin()->second); m->mothurOutEndLine();	}
                    m->mothurOutEndLine();

                    return 0;
                }
                deleteFiles = true;
            }
		//****************** break up files between processes and cluster each file set ******************************//
		
        listFileNames = createProcesses(distName, labels);
        
        if (deleteFiles) {
            //delete the temp files now that we are done
            for (int i = 0; i < distName.size(); i++) {
                string thisNamefile = distName[i].begin()->second;
                string thisDistFile = distName[i].begin()->first;
                util.mothurRemove(thisNamefile);
                util.mothurRemove(thisDistFile);
            }
        }
		
		if (m->getControl_pressed()) { for (int i = 0; i < listFileNames.size(); i++) { util.mothurRemove(listFileNames[i]); } return 0; }
		
		if (!util.isEqual(saveCutoff, cutoff)) { m->mothurOut("\nCutoff was " + toString(saveCutoff) + " changed cutoff to " + toString(cutoff)); m->mothurOutEndLine();  }
		
		m->mothurOut("It took " + toString(time(NULL) - estart) + " seconds to cluster\n");
		
		//****************** merge list file and create rabund and sabund files ******************************//
		estart = time(NULL);
		m->mothurOut("Merging the clustered files...\n");

		ListVector* listSingle;
		map<double, int> labelBins = completeListFile(listFileNames, singletonName, labels, listSingle); //returns map of label to numBins
		
		if (m->getControl_pressed()) { if (listSingle != NULL) { delete listSingle; } for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); } return 0; }
		
		mergeLists(listFileNames, labelBins, listSingle);

		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); } return 0; }
        
        //delete after all are complete incase a crash happens
        if (!deleteFiles) { for (int i = 0; i < distName.size(); i++) {	util.mothurRemove(distName[i].begin()->first); util.mothurRemove(distName[i].begin()->second); 	} }
		
		m->mothurOut("It took " + toString(time(NULL) - estart) + " seconds to merge.\n");
        
        if ((method == "opti") && (runsensSpec)) { runSensSpec();  }
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); } return 0; }
        
		//set list file as new current listfile
		string currentName = "";
		itTypes = outputTypes.find("list");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setListFile(currentName); }
		}
		
		//set rabund file as new current rabundfile
		itTypes = outputTypes.find("rabund");
		if (itTypes != outputTypes.end()) { if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setRabundFile(currentName); } }
		
		//set sabund file as new current sabundfile
		itTypes = outputTypes.find("sabund");
		if (itTypes != outputTypes.end()) { if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setSabundFile(currentName); } }
        
        //set sabund file as new current sabundfile
        itTypes = outputTypes.find("column");
        if (itTypes != outputTypes.end()) { if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setColumnFile(currentName); } }
				
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
map<double, int> ClusterSplitCommand::completeListFile(vector<string> listNames, string singleton, set<string>& userLabels, ListVector*& listSingle){
	try {
		map<double, int> labelBin;
		vector<double> orderFloat;
		int numSingleBins;
        
		//read in singletons
		if (singleton != "none") {
            
            ifstream in;
            util.openInputFile(singleton, in);
				
			string firstCol, secondCol;
			listSingle = new ListVector();
            
            if (countfile != "") { util.getline(in); util.gobble(in); }
            
			while (!in.eof()) {
				in >> firstCol >> secondCol;
                util.getline(in);
				if (countfile == "") { listSingle->push_back(secondCol); }
                else { listSingle->push_back(firstCol); }
                util.gobble(in);
			}
            
			in.close();
			util.mothurRemove(singleton);
            
			numSingleBins = listSingle->getNumBins();
        }else{  listSingle = NULL; numSingleBins = 0;  }
		
        //go through users set and make them floats so we can sort them
        for(set<string>::iterator it = userLabels.begin(); it != userLabels.end(); ++it) {
            double temp = -10.0;
            
            if ((*it != "unique") && (convertTestFloat(*it, temp) ))	{	util.mothurConvert(*it, temp);	}
            else if (*it == "unique")										{	temp = -1.0;		}
            
            if ((temp < cutoff) || util.isEqual(cutoff, temp)) {
                orderFloat.push_back(temp);
                labelBin[temp] = numSingleBins; //initialize numbins
            }
        }
	
		//sort order
		sort(orderFloat.begin(), orderFloat.end());
		userLabels.clear();
        
		//get the list info from each file
		for (int k = 0; k < listNames.size(); k++) {
            
			if (m->getControl_pressed()) {  
				if (listSingle != NULL) { delete listSingle; listSingle = NULL; util.mothurRemove(singleton);  }
				for (int i = 0; i < listNames.size(); i++) {   util.mothurRemove(listNames[i]);  }
				return labelBin;
			}
			
			InputData* input = new InputData(listNames[k], "list", nullVector);
			ListVector* list = input->getListVector();
			string lastLabel = list->getLabel();
            
			string filledInList = listNames[k] + "filledInTemp";
			ofstream outFilled;
			util.openOutputFile(filledInList, outFilled);
            bool printHeaders = true;
            
            
			//for each label needed
			for(int l = 0; l < orderFloat.size(); l++){
                
				string thisLabel;
				if (util.isEqual(orderFloat[l],-1)) { thisLabel = "unique"; }
				else { thisLabel = toString(orderFloat[l],  length-1);  } 
                
				//this file has reached the end
				if (list == NULL) { 
					list = input->getListVector(lastLabel, true); 
				}else{	//do you have the distance, or do you need to fill in
						
					float labelFloat;
					if (list->getLabel() == "unique") {  labelFloat = -1.0;  }
					else { convert(list->getLabel(), labelFloat); }
                    
					//check for missing labels
					if (labelFloat > orderFloat[l]) { //you are missing the label, get the next smallest one
						//if its bigger get last label, otherwise keep it
						delete list;
						list = input->getListVector(lastLabel, true);  //get last list vector to use, you actually want to move back in the file
					}
					lastLabel = list->getLabel();
				}
				
				//print to new file
				list->setLabel(thisLabel);
                list->setPrintedLabels(printHeaders);
                list->print(outFilled, true); printHeaders = false;
		
				//update labelBin
				labelBin[orderFloat[l]] += list->getNumBins();
									
				delete list;
									
				list = input->getListVector();
			}
			
			if (list != NULL) { delete list; }
			delete input;
			
			outFilled.close();
			util.mothurRemove(listNames[k]);
			rename(filledInList.c_str(), listNames[k].c_str());
		}
		
		return labelBin;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "completeListFile");
		exit(1);
	}
}
//**********************************************************************************************************************
int ClusterSplitCommand::mergeLists(vector<string> listNames, map<double, int> userLabels, ListVector* listSingle){
	try {
		if (outputDir == "") { outputDir += util.hasPath(distfile); }
		fileroot = outputDir + util.getRootName(util.getSimpleName(distfile));
		
        map<string, string> variables; 
        variables["[filename]"] = fileroot;
        variables["[clustertag]"] = tag;
        string sabundFileName = getOutputFileName("sabund", variables);
        string rabundFileName = getOutputFileName("rabund", variables);
        //if (countfile != "") { variables["[tag2]"] = "unique_list"; }
        string listFileName = getOutputFileName("list", variables);
        
        map<string, int> counts;
        ofstream outList, outRabund, outSabund;
        if (countfile == "") {
            util.openOutputFile(sabundFileName,	outSabund);
            util.openOutputFile(rabundFileName,	outRabund);
            outputNames.push_back(sabundFileName); outputTypes["sabund"].push_back(sabundFileName);
            outputNames.push_back(rabundFileName); outputTypes["rabund"].push_back(rabundFileName);
            
        }else {
            if (file == "") {
                CountTable ct;
                ct.readTable(countfile, false, false);
                counts = ct.getNameMap();
            }
        }
        
		util.openOutputFile(listFileName,	outList);
        outputNames.push_back(listFileName); outputTypes["list"].push_back(listFileName);
        bool printHeaders = true;

		//for each label needed
		for(map<double, int>::iterator itLabel = userLabels.begin(); itLabel != userLabels.end(); itLabel++) {
			
			string thisLabel;
			if (util.isEqual(itLabel->first,-1)) { thisLabel = "unique"; }
			else { thisLabel = toString(itLabel->first,  length-1);  } 
			
			//outList << thisLabel << '\t' << itLabel->second << '\t';
            
            RAbundVector* rabund = NULL;
            ListVector completeList;
            completeList.setLabel(thisLabel);
            
            if (countfile == "") {
                rabund = new RAbundVector();
                rabund->setLabel(thisLabel);
            }

			//add in singletons
			if (listSingle != NULL) {
				for (int j = 0; j < listSingle->getNumBins(); j++) {
					//outList << listSingle->get(j) << '\t';
                    completeList.push_back(listSingle->get(j));
					if (countfile == "") { rabund->push_back(util.getNumNames(listSingle->get(j))); }
				}
			}
			
			//get the list info from each file
			for (int k = 0; k < listNames.size(); k++) {
	
				if (m->getControl_pressed()) {  if (listSingle != NULL) { delete listSingle;   } for (int i = 0; i < listNames.size(); i++) { util.mothurRemove(listNames[i]);  } if (rabund != NULL) { delete rabund; } return 0; }
				
				InputData* input = new InputData(listNames[k], "list", nullVector);
				ListVector* list = input->getListVector(thisLabel);
				
				//this file has reached the end
				if (list == NULL) { m->mothurOut("Error merging listvectors in file " + listNames[k]); m->mothurOutEndLine();  }	
				else {		
					for (int j = 0; j < list->getNumBins(); j++) {
                        completeList.push_back(list->get(j));
						if (countfile == "") { rabund->push_back(util.getNumNames(list->get(j))); }
					}
					delete list;
				}
				delete input;
			}
			
            if (countfile == "") {
                SAbundVector sabund = rabund->getSAbundVector();
                sabund.print(outSabund);
                rabund->print(outRabund);
            }

            completeList.setPrintedLabels(printHeaders);
            if (countfile == "") { completeList.print(outList);  printHeaders = false; }
            else if ((file == "") && (countfile != "")) { completeList.print(outList, counts);  printHeaders = false; }
            else { completeList.print(outList);  printHeaders = false; }
			
			if (rabund != NULL) { delete rabund; }
		}
		
		outList.close();
        if (countfile == "") {
            outRabund.close();
            outSabund.close();
		}
		if (listSingle != NULL) { delete listSingle;  }
		
		for (int i = 0; i < listNames.size(); i++) {  util.mothurRemove(listNames[i]);  }
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "mergeLists");
		exit(1);
	}
}
/**************************************************************************************************/
struct clusterData {
    MothurOut* m;
    Utils util;
    int count, precision, length, numSingletons, maxIters;
    bool showabund, classic, useName, useCount, deleteFiles, cutoffNotSet;
    double cutoff, stableMetric;
    ofstream outList, outRabund, outSabund;
    string tag, method, outputDir, vsearchLocation, metricName, initialize;
    vector< map<string, string> > distNames;
    set<string> labels;
    vector<string> listFileNames;
    
    clusterData(){}
    clusterData(bool showab, bool cla, bool df, vector< map<string, string> > dN, bool cns, double cu, int prec, int len, string meth, string opd, string vl) {
        showabund = showab;
        distNames = dN;
        cutoff = cu;
        classic = cla;
        method = meth;
        precision = prec;
        length = len;
        outputDir = opd;
        vsearchLocation = vl;
        deleteFiles = df;
        cutoffNotSet = cns;
        m = MothurOut::getInstance();
        count = 0;
        useName = false;
        useCount = false;
        numSingletons = 0;
    }
    void setOptiOptions(string metn, double stabMet, string init, int mxi ) {
        metricName = metn;
        stableMetric = stabMet;
        maxIters = mxi;
        initialize = init;
    }
    void setNamesCount(string nmf, string cnf) {
        useName = false;
        useCount = false;
        if (nmf != "") { useName = true;  }
        if (cnf != "") { useCount = true; }
    }
};
//**********************************************************************************************************************
int createRabund(CountTable*& ct, ListVector*& list, RAbundVector*& rabund, clusterData* params){
    try {
        rabund->setLabel(list->getLabel());
        for(int i = 0; i < list->getNumBins(); i++) {
            if (params->m->getControl_pressed()) { break; }
            vector<string> binNames;
            string bin = list->get(i);
            params->util.splitAtComma(bin, binNames);
            int total = 0;
            for (int j = 0; j < binNames.size(); j++) { total += ct->getNumSeqs(binNames[j]);  }
            rabund->push_back(total);
        }
        return 0;
    }
    catch(exception& e) {
        params->m->errorOut(e, "ClusterSplitCommand", "createRabund");
        exit(1);
    }
    
}
//**********************************************************************************************************************
string clusterClassicFile(string thisDistFile, string thisNamefile, double& smallestCutoff, clusterData* params){
    try {
        string listFileName = "";
        
        ListVector* list = NULL;
        ListVector oldList;
        RAbundVector* rabund = NULL;
        
        params->m->mothurOut("\nReading " + thisDistFile + "\n");
        
        //reads phylip file storing data in 2D vector, also fills list and rabund
        bool sim = false;
        ClusterClassic cluster(params->cutoff, params->method, sim);
        
        NameAssignment* nameMap = NULL;
        CountTable* ct = NULL;
        if(params->useName){
            nameMap = new NameAssignment(thisNamefile);
            nameMap->readMap();
            cluster.readPhylipFile(thisDistFile, nameMap);
        }else if (params->useCount) {
            ct = new CountTable();
            ct->readTable(thisNamefile, false, false);
            cluster.readPhylipFile(thisDistFile, ct);
        }
        params->tag = cluster.getTag();
        
        if (params->m->getControl_pressed()) { if(params->useName){	delete nameMap; }
            else if (params->useCount) { delete ct; } return listFileName; }
        
        list = cluster.getListVector();
        rabund = cluster.getRAbundVector();
        
        string thisOutputDir = params->outputDir;
        if (params->outputDir == "") { thisOutputDir += params->util.hasPath(thisDistFile); }
        string fileroot = thisOutputDir + params->util.getRootName(params->util.getSimpleName(thisDistFile));
        listFileName = fileroot+ params->tag + ".list";
        
        ofstream listFile;
        params->util.openOutputFile(fileroot+ params->tag + ".list",	listFile);
        
        float previousDist = 0.00000;
        float rndPreviousDist = 0.00000;
        bool printHeaders = true;
        oldList = *list;
        
        params->m->mothurOut("\nClustering " + thisDistFile + "\n");
        
        while ((cluster.getSmallDist() < params->cutoff) && (cluster.getNSeqs() > 1)){
            if (params->m->getControl_pressed()) {  delete list; delete rabund; listFile.close();  if(params->useName){	delete nameMap; }
                else if (params->useCount) { delete ct; } return listFileName;  }
            
            cluster.update(params->cutoff);
            
            float dist = cluster.getSmallDist();
            float rndDist = params->util.ceilDist(dist, params->precision);
            
            if(previousDist <= 0.0000 && !params->util.isEqual(dist, previousDist)){
                oldList.setLabel("unique");
                oldList.setPrintedLabels(printHeaders);
                oldList.print(listFile); printHeaders = false;
                if (params->labels.count("unique") == 0) {  params->labels.insert("unique");  }
            }
            else if(!params->util.isEqual(rndDist, rndPreviousDist)){
                oldList.setLabel(toString(rndPreviousDist,  params->length-1));
                oldList.setPrintedLabels(printHeaders);
                oldList.print(listFile); printHeaders = false;
                if (params->labels.count(toString(rndPreviousDist,  params->length-1)) == 0) { params->labels.insert(toString(rndPreviousDist,  params->length-1)); }
            }
            
            
            previousDist = dist;
            rndPreviousDist = rndDist;
            oldList = *list;
        }
        
        if(previousDist <= 0.0000){
            oldList.setLabel("unique");
            oldList.setPrintedLabels(printHeaders);
            oldList.print(listFile); printHeaders = false;
            if (params->labels.count("unique") == 0) { params->labels.insert("unique"); }
        }
        else if(rndPreviousDist<params->cutoff){
            oldList.setLabel(toString(rndPreviousDist,  params->length-1));
            oldList.setPrintedLabels(printHeaders);
            oldList.print(listFile); printHeaders = false;
            if (params->labels.count(toString(rndPreviousDist,  params->length-1)) == 0) { params->labels.insert(toString(rndPreviousDist,  params->length-1)); }
        }
        listFile.close();
        
        delete list; delete rabund;
        if(params->useName)         {	delete nameMap; }
        else if (params->useCount)  {   delete ct;      }
        
        if (params->deleteFiles) {
            params->util.mothurRemove(thisDistFile);
            params->util.mothurRemove(thisNamefile);
        }
        return listFileName;
    }
    catch(exception& e) {
        params->m->errorOut(e, "ClusterSplitCommand", "clusterClassicFile");
        exit(1);
    }
}
//**********************************************************************************************************************
string runOptiCluster(string thisDistFile, string thisNamefile, double& smallestCutoff, clusterData* params){
    try {
        if (params->cutoffNotSet) {  params->m->mothurOut("\nYou did not set a cutoff, using 0.03.\n"); params->cutoff = 0.03;  }
        
        string nameOrCount = "count";
        if (params->useName) { nameOrCount = "name"; }
        
        OptiMatrix matrix(thisDistFile, thisNamefile, nameOrCount, "column", params->cutoff, false);
        
        ClusterMetric* metric = NULL;
        if (params->metricName == "mcc")             { metric = new MCC();              }
        else if (params->metricName == "sens")       { metric = new Sensitivity();      }
        else if (params->metricName == "spec")       { metric = new Specificity();      }
        else if (params->metricName == "tptn")       { metric = new TPTN();             }
        else if (params->metricName == "tp")         { metric = new TP();               }
        else if (params->metricName == "tn")         { metric = new TN();               }
        else if (params->metricName == "fp")         { metric = new FP();               }
        else if (params->metricName == "fn")         { metric = new FN();               }
        else if (params->metricName == "f1score")    { metric = new F1Score();          }
        else if (params->metricName == "accuracy")   { metric = new Accuracy();         }
        else if (params->metricName == "ppv")        { metric = new PPV();              }
        else if (params->metricName == "npv")        { metric = new NPV();              }
        else if (params->metricName == "fdr")        { metric = new FDR();              }
        else if (params->metricName == "fpfn")       { metric = new FPFN();             }
        
        OptiCluster cluster(&matrix, metric, 0);
        params->tag = cluster.getTag();
        
        params->m->mothurOut("\nClustering " + thisDistFile + "\n");
        
        string thisOutputDir = params->outputDir;
        if (params->outputDir == "") { thisOutputDir += params->util.hasPath(thisDistFile); }
        string fileroot = thisOutputDir + params->util.getRootName(params->util.getSimpleName(thisDistFile));
        string listFileName = fileroot+ params->tag + ".list";

        int iters = 0;
        double listVectorMetric = 0; //worst state
        double delta = 1;
        
        cluster.initialize(listVectorMetric, true, params->initialize);
        
        while ((delta > params->stableMetric) && (iters < params->maxIters)) {
            
            if (params->m->getControl_pressed()) { if (params->deleteFiles) { params->util.mothurRemove(thisDistFile);  params->util.mothurRemove(thisNamefile); } return listFileName; }
            double oldMetric = listVectorMetric;
            cluster.update(listVectorMetric);
            
            delta = abs(oldMetric - listVectorMetric);
            iters++;
        }
        
        if (params->m->getControl_pressed()) { delete metric; metric = NULL; return 0; }
        
        ListVector* list = cluster.getList();
        list->setLabel(toString(smallestCutoff));
        //params->cutoff = params->util.ceilDist(params->cutoff, params->precision);
        params->labels.insert(toString(smallestCutoff));
        
        ofstream listFile;
        params->util.openOutputFile(listFileName,	listFile);
        list->print(listFile);
        listFile.close();
        
        if (params->deleteFiles) {
            params->util.mothurRemove(thisDistFile);
            params->util.mothurRemove(thisNamefile);
        }
        
        long long tp, tn, fp, fn;
        params->m->mothurOut("\ntp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");
        vector<double> results = cluster.getStats(tp, tn, fp, fn);
        params->m->mothurOut(toString(tp) + "\t" + toString(tn) + "\t" + toString(fp) + "\t" + toString(fn) + "\t");
        for (int i = 0; i < results.size(); i++) { params->m->mothurOut(toString(results[i]) + "\t");  }
        params->m->mothurOut("\n\n");
        return listFileName;
        
    }
    catch(exception& e) {
        params->m->errorOut(e, "ClusterSplitCommand", "runOptiCluster");
        exit(1);
    }
}
//**********************************************************************************************************************

int vsearchDriver(string inputFile, string ucClusteredFile, string logfile, double cutoff, clusterData* params){
    try {
        
        //vsearch --maxaccepts 16 --usersort --id 0.97 --minseqlength 30 --wordlength 8 --uc $ROOT.clustered.uc --cluster_smallmem $ROOT.sorted.fna --maxrejects 64 --strand both --log $ROOT.clustered.log --sizeorder
        
        
        ucClusteredFile = params->util.getFullPathName(ucClusteredFile);
        inputFile = params->util.getFullPathName(inputFile);
        logfile = params->util.getFullPathName(logfile);
        
        //to allow for spaces in the path
        ucClusteredFile = "\"" + ucClusteredFile + "\"";
        inputFile = "\"" + inputFile + "\"";
        logfile = "\"" + logfile + "\"";
        
        vector<char*> cPara;
        
        string vsearchCommand = params->vsearchLocation;
        vsearchCommand = "\"" + vsearchCommand + "\" ";
        
        vector<char*> vsearchParameters;
        char* vsearchParameter = new char[vsearchCommand.length()+1];  vsearchParameter[0] = '\0'; strncat(vsearchParameter, vsearchCommand.c_str(), vsearchCommand.length());
        vsearchParameters.push_back(vsearchParameter);
        
        //--maxaccepts=16
        char* maxaccepts = new char[16];  maxaccepts[0] = '\0'; strncat(maxaccepts, "--maxaccepts=16", 15);
        vsearchParameters.push_back(maxaccepts);
        
        //--threads=1
        char* threads = new char[12];  threads[0] = '\0'; strncat(threads, "--threads=1", 11);
        vsearchParameters.push_back(threads);
        
        //--usersort
        char* usersort = new char[11];  usersort[0] = '\0'; strncat(usersort, "--usersort", 10);
        vsearchParameters.push_back(usersort);
        
        //--id=0.97
        cutoff = abs(1.0 - cutoff); string cutoffString = toString(cutoff);
        if (cutoffString.length() > 4) {  cutoffString = cutoffString.substr(0, 4);  }
        else if (cutoffString.length() < 4)  {  for (int i = cutoffString.length(); i < 4; i++)  { cutoffString += "0";  } }
        
        cutoffString = "--id=" +  cutoffString;
        char* cutoffParameter = new char[cutoffString.length()+1];  cutoffParameter[0] = '\0'; strncat(cutoffParameter, cutoffString.c_str(), cutoffString.length());
        vsearchParameters.push_back(cutoffParameter);
        
        //--minseqlength=30
        char* minseqlength = new char[18];  minseqlength[0] = '\0'; strncat(minseqlength, "--minseqlength=30", 17);
        vsearchParameters.push_back(minseqlength);
        
        //--wordlength=8
        char* wordlength = new char[15];  wordlength[0] = '\0'; strncat(wordlength, "--wordlength=8", 14);
        vsearchParameters.push_back(wordlength);
        
        //--uc=$ROOT.clustered.uc
        string tempIn = "--uc=" + ucClusteredFile;
        char* uc = new char[tempIn.length()+1];  uc[0] = '\0'; strncat(uc, tempIn.c_str(), tempIn.length());
        vsearchParameters.push_back(uc);
        
        //--cluster_smallmem $ROOT.sorted.fna
        string tempSorted = "--cluster_smallmem=" + inputFile;
        char* cluster_smallmen = new char[tempSorted.length()+1];  cluster_smallmen[0] = '\0'; strncat(cluster_smallmen, tempSorted.c_str(), tempSorted.length());
        vsearchParameters.push_back(cluster_smallmen);
        
        //--maxrejects=64
        char* maxrejects = new char[16];  maxrejects[0] = '\0'; strncat(maxrejects, "--maxrejects=64", 15);
        vsearchParameters.push_back(maxrejects);
        
        //--strand=both
        char* strand = new char[14];  strand[0] = '\0'; strncat(strand, "--strand=both", 13);
        vsearchParameters.push_back(strand);
        
        //--log=$ROOT.clustered.log
        string tempLog = "--log=" + logfile;
        char* log = new char[tempLog.length()+1];  log[0] = '\0'; strncat(log, tempLog.c_str(), tempLog.length());
        vsearchParameters.push_back(log);
        
        if (params->method == "agc") {
            //--sizeorder
            char* sizeorder = new char[12];  sizeorder[0] = '\0'; strncat(sizeorder, "--sizeorder", 11);
            vsearchParameters.push_back(sizeorder);
        }
        
        if (params->m->getDebug()) {  params->m->mothurOut("[DEBUG]: "); for(int i = 0; i < vsearchParameters.size(); i++)  { params->m->mothurOut(toString(vsearchParameters[i]) + "\t"); } params->m->mothurOut("\n");  }
        
        string commandString = "";
        for (int i = 0; i < vsearchParameters.size(); i++) {    commandString += toString(vsearchParameters[i]) + " "; }
        
#if defined NON_WINDOWS
#else
        commandString = "\"" + commandString + "\"";
#endif
        if (params->m->getDebug()) {  params->m->mothurOut("[DEBUG]: vsearch cluster command = " + commandString + ".\n"); }
        
        system(commandString.c_str());
        
        //free memory
        for(int i = 0; i < vsearchParameters.size(); i++)  {  delete vsearchParameters[i];  }
        
        //remove "" from filenames
        ucClusteredFile = ucClusteredFile.substr(1, ucClusteredFile.length()-2);
        inputFile = inputFile.substr(1, inputFile.length()-2);
        logfile = logfile.substr(1, logfile.length()-2);
        
        return 0;
        
    }
    catch(exception& e) {
        params->m->errorOut(e, "ClusterSplitCommand", "vsearchDriver");
        exit(1);
    }
}
//**********************************************************************************************************************
string runVsearchCluster(string thisDistFile, string thisNamefile, double& smallestCutoff, clusterData* params){
    try {
        
        params->m->mothurOut("\nClustering " + thisDistFile + "\n");
        
        string vsearchFastafile = ""; VsearchFileParser* vParse;
        if (params->useName)                    { vParse = new VsearchFileParser(thisDistFile, thisNamefile, "name");       }
        else if (params->useCount)              { vParse = new VsearchFileParser(thisDistFile, thisNamefile, "count");      }
        else                                    { params->m->mothurOut("[ERROR]: Opps, should never get here. ClusterSplitCommand::runVsearchCluster() \n"); params->m->setControl_pressed(true); }
        
        if (params->m->getControl_pressed()) {  return ""; }
        
        vsearchFastafile = vParse->getVsearchFile();
        
        if (params->cutoff > 1.0) {  params->m->mothurOut("You did not set a cutoff, using 0.03.\n"); params->cutoff = 0.03; }
        
        //Run vsearch
        string ucVsearchFile = params->util.getSimpleName(vsearchFastafile) + ".clustered.uc";
        string logfile = params->util.getSimpleName(vsearchFastafile) + ".clustered.log";
        vsearchDriver(vsearchFastafile, ucVsearchFile, logfile, smallestCutoff, params);
        
        if (params->m->getControl_pressed()) { params->util.mothurRemove(ucVsearchFile); params->util.mothurRemove(logfile);  params->util.mothurRemove(vsearchFastafile); return ""; }
        
        string thisOutputDir = params->outputDir;
        if (params->outputDir == "") { thisOutputDir += params->util.hasPath(thisDistFile); }
        params->tag = params->method;
        string listFileName = thisOutputDir + params->util.getRootName(params->util.getSimpleName(thisDistFile)) + params->tag + ".list";
        
        //Convert outputted *.uc file into a list file
        map<string, int> counts;
        ListVector list = vParse->createListFile(ucVsearchFile, vParse->getNumBins(logfile), toString(params->cutoff), counts);
        
        ofstream out;
        params->util.openOutputFile(listFileName,	out);

        list.DataVector::printHeaders(out);
        
        if (params->useCount) { list.print(out, counts); }
        else { list.print(out); } delete vParse;

        //remove temp files
        params->util.mothurRemove(ucVsearchFile); params->util.mothurRemove(logfile);  params->util.mothurRemove(vsearchFastafile);
        
        if (params->deleteFiles) {
            params->util.mothurRemove(thisDistFile);
            params->util.mothurRemove(thisNamefile);
        }
        params->labels.insert(toString(params->cutoff));
        
        return listFileName;
    }
    catch(exception& e) {
        params->m->errorOut(e, "ClusterSplitCommand", "runVsearchCluster");
        exit(1);
    }
}
//**********************************************************************************************************************
string clusterFile(string thisDistFile, string thisNamefile, double& smallestCutoff, clusterData* params){
    try {
        string listFileName = "";
        
        if ((params->method == "agc") || (params->method == "dgc")) {  listFileName = runVsearchCluster(thisDistFile, thisNamefile, smallestCutoff, params);  }
        else if (params->method == "opti")                          {  listFileName = runOptiCluster(thisDistFile, thisNamefile, smallestCutoff, params);     }
        else {
            
            Cluster* cluster = NULL;
            SparseDistanceMatrix* matrix = NULL;
            ListVector* list = NULL;
            ListVector oldList;
            RAbundVector* rabund = NULL;
            
            if (params->m->getControl_pressed()) { return listFileName; }
            
            params->m->mothurOut("\nReading " + thisDistFile + "\n");
            
            ReadMatrix* read = new ReadColumnMatrix(thisDistFile);
            read->setCutoff(params->cutoff);
            
            NameAssignment* nameMap = NULL;
            CountTable* ct = NULL;
            if(params->useName){
                nameMap = new NameAssignment(thisNamefile);
                nameMap->readMap();
                read->read(nameMap);
            }else if (params->useCount) {
                ct = new CountTable();
                ct->readTable(thisNamefile, false, false);
                read->read(ct);
            }else { read->read(nameMap); }
            
            list = read->getListVector();
            matrix = read->getDMatrix();
            
            if(params->useCount) {
                rabund = new RAbundVector();
                createRabund(ct, list, rabund, params); //creates an rabund that includes the counts for the unique list
                delete ct;
            }else { rabund = new RAbundVector(list->getRAbundVector()); }
            
            delete read;  read = NULL;
            if (params->useName) { delete nameMap; nameMap = NULL; }
            
            params->m->mothurOut("\nClustering " + thisDistFile + "\n");
            
            //create cluster
            float adjust = -1.0;
            if (params->method == "furthest")	{	cluster = new CompleteLinkage(rabund, list, matrix, params->cutoff, params->method, adjust); }
            else if(params->method == "nearest"){	cluster = new SingleLinkage(rabund, list, matrix, params->cutoff, params->method, adjust); }
            else if(params->method == "average"){	cluster = new AverageLinkage(rabund, list, matrix, params->cutoff, params->method, adjust);	}
            params->tag = cluster->getTag();
            
            string thisOutputDir = params->outputDir;
            if (params->outputDir == "") { thisOutputDir += params->util.hasPath(thisDistFile); }
            string fileroot = thisOutputDir + params->util.getRootName(params->util.getSimpleName(thisDistFile));
            listFileName = fileroot+ params->tag + ".list";

            ofstream listFile;
            params->util.openOutputFile(listFileName,	listFile);

            float previousDist = 0.00000;
            float rndPreviousDist = 0.00000;
            bool printHeaders = true;
            
            oldList = *list;
            
            double saveCutoff = params->cutoff;
            
            while (matrix->getSmallDist() < params->cutoff && matrix->getNNodes() > 0){
                
                if (params->m->getControl_pressed()) { //clean up
                    delete matrix; delete list;	delete cluster; delete rabund;
                    listFile.close();
                    params->util.mothurRemove(listFileName);
                    return listFileName;
                }
                
                cluster->update(saveCutoff);
                
                float dist = matrix->getSmallDist();
                float rndDist = params->util.ceilDist(dist, params->precision);
                
                if(previousDist <= 0.0000 && !params->util.isEqual(dist, previousDist)){
                    oldList.setLabel("unique");
                    oldList.setPrintedLabels(printHeaders);
                    oldList.print(listFile); printHeaders = false;
                    if (params->labels.count("unique") == 0) {  params->labels.insert("unique");  }
                }
                else if(!params->util.isEqual(rndDist, rndPreviousDist)){
                    oldList.setPrintedLabels(printHeaders);
                    oldList.setLabel(toString(rndPreviousDist,  params->length-1));
                    oldList.setPrintedLabels(printHeaders);
                    oldList.print(listFile); printHeaders = false;
                    if (params->labels.count(toString(rndPreviousDist,  params->length-1)) == 0) { params->labels.insert(toString(rndPreviousDist,  params->length-1)); }
                }
                
                previousDist = dist;
                rndPreviousDist = rndDist;
                oldList = *list;
            }
            
            
            if(previousDist <= 0.0000){
                oldList.setLabel("unique");
                oldList.setPrintedLabels(printHeaders);
                oldList.print(listFile); printHeaders = false;
                if (params->labels.count("unique") == 0) { params->labels.insert("unique"); }
            }
            else if(rndPreviousDist<params->cutoff){
                oldList.setLabel(toString(rndPreviousDist,  params->length-1));
                oldList.setPrintedLabels(printHeaders);
                oldList.print(listFile); printHeaders = false;
                if (params->labels.count(toString(rndPreviousDist,  params->length-1)) == 0) { params->labels.insert(toString(rndPreviousDist,  params->length-1)); }
            }
            
            delete matrix; delete list;	delete cluster; delete rabund;
            matrix = NULL; list = NULL; cluster = NULL; rabund = NULL;
            listFile.close();
            
            if (params->m->getControl_pressed()) { //clean up
                params->util.mothurRemove(listFileName);
                return listFileName;
            }
            
            if (params->deleteFiles) {
                params->util.mothurRemove(thisDistFile);
                params->util.mothurRemove(thisNamefile);
            }
            
            if (!params->util.isEqual(saveCutoff, params->cutoff)) {
                saveCutoff = params->util.ceilDist(saveCutoff, params->precision);
                params->m->mothurOut("Cutoff was " + toString(params->cutoff) + " changed cutoff to " + toString(saveCutoff) + "\n");
            }
            
            if (saveCutoff < smallestCutoff) { smallestCutoff = saveCutoff;  }
        }
        return listFileName;
    }
    catch(exception& e) {
        params->m->errorOut(e, "ClusterSplitCommand", "clusterFile");
        exit(1);
    }
}
//**********************************************************************************************************************
void cluster(clusterData* params){
    try {
        vector<string> listFileNames;
        double smallestCutoff = params->cutoff;
        
        //cluster each distance file
        for (int i = 0; i < params->distNames.size(); i++) {
            
            string thisNamefile = params->distNames[i].begin()->second;
            string thisDistFile = params->distNames[i].begin()->first;
            
            string listFileName = "";
            if (params->classic)    {  listFileName = clusterClassicFile(thisDistFile, thisNamefile, smallestCutoff, params);   }
            else                    {  listFileName = clusterFile(thisDistFile, thisNamefile, smallestCutoff, params);          }
            
            if (params->m->getControl_pressed()) { //clean up
                for (int i = 0; i < listFileNames.size(); i++) {	params->util.mothurRemove(listFileNames[i]); 	}
                params->listFileNames.clear(); break;
            }
            params->listFileNames.push_back(listFileName);
        }
        params->cutoff = smallestCutoff;
    }
    catch(exception& e) {
        params->m->errorOut(e, "ClusterSplitCommand", "cluster");
        exit(1);
    }
    
    
}
//**********************************************************************************************************************
void printData(ListVector* oldList, clusterData* params){
	try {
		string label = oldList->getLabel();
		RAbundVector oldRAbund = oldList->getRAbundVector();
		
		oldRAbund.setLabel(label);
		if (params->showabund) {
			oldRAbund.getSAbundVector().print(cout);
		}
		oldRAbund.print(params->outRabund);
		oldRAbund.getSAbundVector().print(params->outSabund);
	
		oldList->print(params->outList, true);
	}
	catch(exception& e) {
		params->m->errorOut(e, "ClusterSplitCommand", "printData");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string>  ClusterSplitCommand::createProcesses(vector< map<string, string> > distName, set<string>& labels){
	try {
        //sanity check
        if (processors > distName.size()) { processors = distName.size(); }
        deleteFiles = false; //so if we need to recalc the processors the files are still there
        vector<string> listFiles;
        vector < vector < map<string, string> > > dividedNames; //distNames[1] = vector of filenames for process 1...
        dividedNames.resize(processors);
        
        //for each file group figure out which process will complete it
        //want to divide the load intelligently so the big files are spread between processes
        for (int i = 0; i < distName.size(); i++) {
            int processToAssign = (i+1) % processors; 
            if (processToAssign == 0) { processToAssign = processors; }
            
            dividedNames[(processToAssign-1)].push_back(distName[i]);
            if ((processToAssign-1) == 1) { m->mothurOut(distName[i].begin()->first + "\n"); }
        }
        
        //now lets reverse the order of ever other process, so we balance big files running with little ones
        for (int i = 0; i < processors; i++) {
            int remainder = ((i+1) % processors);
            if (remainder) {  reverse(dividedNames[i].begin(), dividedNames[i].end());  }
        }
        
        if (m->getControl_pressed()) { return listFiles; }
        
        
        //create array of worker threads
        vector<thread*> workerThreads;
        vector<clusterData*> data;
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            clusterData* dataBundle = new clusterData(showabund, classic, deleteFiles, dividedNames[i+1], cutoffNotSet, cutoff, precision, length, method, outputDir, vsearchLocation);
            dataBundle->setOptiOptions(metricName, stableMetric, initialize, maxIters);
            dataBundle->setNamesCount(namefile, countfile);
            data.push_back(dataBundle);
            
            workerThreads.push_back(new thread(cluster, dataBundle));
        }
        
        
        clusterData* dataBundle = new clusterData(showabund, classic, deleteFiles, dividedNames[0], cutoffNotSet, cutoff, precision, length, method, outputDir, vsearchLocation);
        dataBundle->setOptiOptions(metricName, stableMetric, initialize, maxIters);
        dataBundle->setNamesCount(namefile, countfile);
        cluster(dataBundle);
        listFiles = dataBundle->listFileNames;
        tag = dataBundle->tag;
        cutoff = dataBundle->cutoff;
        labels = dataBundle->labels;
        
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            
            listFiles.insert(listFiles.end(), data[i]->listFileNames.begin(), data[i]->listFileNames.end());
            labels.insert(data[i]->labels.begin(), data[i]->labels.end());
            if (data[i]->cutoff < cutoff) { cutoff = data[i]->cutoff; }
            
            delete data[i];
            delete workerThreads[i];
        }
        delete dataBundle;
        deleteFiles = true;
        
        return listFiles;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "createProcesses");
		exit(1);
	}
}
//**********************************************************************************************************************

int ClusterSplitCommand::createMergedDistanceFile(vector< map<string, string> > distNames) {
	try{
		string thisOutputDir = outputDir;
		if (outputDir == "") { thisOutputDir = util.hasPath(fastafile); }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastafile));
		string outputFileName = getOutputFileName("column", variables);
		util.mothurRemove(outputFileName);
		
		
		for (int i = 0; i < distNames.size(); i++) {
			if (m->getControl_pressed()) {  return 0; }
			
			string thisDistFile = distNames[i].begin()->first;
			
			util.appendFiles(thisDistFile, outputFileName);
		}	
			
		outputTypes["column"].push_back(outputFileName); outputNames.push_back(outputFileName);
			
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "createMergedDistanceFile");
		exit(1);
	}
}
//**********************************************************************************************************************

int ClusterSplitCommand::runSensSpec() {
    try{
        string listFile = "";
        itTypes = outputTypes.find("list");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { listFile = (itTypes->second)[0];  }
        }
        
        string columnFile = ""; string phylipFile = "";
        if (makeDist) {
            itTypes = outputTypes.find("column");
            if (itTypes != outputTypes.end()) {
                if ((itTypes->second).size() != 0) { columnFile = (itTypes->second)[0];  }
            }
        }else if (columnfile != "") { columnFile = columnfile; }
        else { phylipFile = phylipfile; }
    
        string inputString = "cutoff=" + toString(cutoff) + ", list=" + listFile;
        if (columnFile != "") { inputString += ", column=" + columnFile;  }
        else if (phylipfile != "")   { inputString += ", phylip=" + phylipfile; }
        else { m->mothurOut("[WARNING]: Cannot run sens.spec analysis without a phylip or column file, skipping."); return 0;  }

        if (namefile != "")         {  inputString += ", name=" + namefile; }
        else if (countfile != "")   { inputString += ", count=" + countfile; }
        else { m->mothurOut("[WARNING]: Cannot run sens.spec analysis without a name or count file, skipping."); return 0;  }
        
        m->mothurOut("/******************************************/\n");
        m->mothurOut("Running command: sens.spec(" + inputString + ")\n");
        current->setMothurCalling(true);
        
        Command* sensspecCommand = new SensSpecCommand(inputString);
        sensspecCommand->execute();
        
        map<string, vector<string> > filenames = sensspecCommand->getOutputFiles();
        
        delete sensspecCommand;
         current->setMothurCalling(false);
        
        string outputFileName = filenames["sensspec"][0];
        
        outputTypes["sensspec"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        
        m->mothurOut("/******************************************/\n"); 
        m->mothurOut("Done.\n\n"); m->mothurOutEndLine();
        
        ifstream in;
        util.openInputFile(outputFileName, in);
        
        while(!in.eof()){
            if (m->getControl_pressed()) { break; }
            
            m->mothurOut(util.getline(in)+"\n"); util.gobble(in);
        }
        in.close();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterSplitCommand", "runSensSpec");
        exit(1);
    }
}
//**********************************************************************************************************************
string ClusterSplitCommand::printFile(string singleton, vector< map<string, string> >& distName){
    try {
        ofstream out;
        map<string, string> variables;
        string thisOutputDir = outputDir;
		if (outputDir == "") { thisOutputDir = util.hasPath(distfile); }
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(distfile));
		string outputFileName = getOutputFileName("file", variables);
        util.openOutputFile(outputFileName, out);
        outputTypes["file"].push_back(outputFileName); outputNames.push_back(outputFileName);
        current->setFileFile(outputFileName);
        
        out << singleton << endl;
        if (namefile != "") { out << "name" << endl; }
        else if (countfile != "") { out << "count" << endl; }
        else { out << "unknown" << endl; }
        
        for (int i = 0; i < distName.size(); i++) {    out << distName[i].begin()->first << '\t' << distName[i].begin()->second << endl;	}
        out.close();
        
        return outputFileName;
    }
    catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "printFile");
		exit(1);
	}
    
}
//**********************************************************************************************************************
string ClusterSplitCommand::readFile(vector< map<string, string> >& distName){
    try {
        
        string singleton, thiscolumn, thisname, type;
        
        ifstream in;
        util.openInputFile(file, in);
        
        in >> singleton; util.gobble(in);
        
        string path = util.hasPath(singleton);
        if (path == "") {  singleton = inputDir + singleton;  }
        
        in >> type; util.gobble(in);
        
        if (type == "name") { namefile = "name"; }
        else if (type == "count") { countfile = "count"; }
        else {  m->mothurOut("[ERROR]: unknown file type. Are the files in column 2 of the file name files or count files? Please change unknown to name or count.\n"); m->setControl_pressed(true); }
        
        if (isList) {

            vector<string> listFileNames;
            string thisListFileName = "";
            set<string> listLabels;
            
            while(!in.eof()) {
                if (m->getControl_pressed()) { break; }
                
                in >> thisListFileName; util.gobble(in);
                
                string path = util.hasPath(thisListFileName);
                if (path == "") {  thisListFileName = inputDir + thisListFileName;  }
                
                getLabels(thisListFileName, listLabels);
                listFileNames.push_back(thisListFileName);
            }
            
            ListVector* listSingle;
            map<double, int> labelBins = completeListFile(listFileNames, singleton, listLabels, listSingle);
            
            mergeLists(listFileNames, labelBins, listSingle);
        
        }else {
            
            while(!in.eof()) {
                if (m->getControl_pressed()) { break; }
                
                in >> thiscolumn; util.gobble(in);
                in >> thisname; util.gobble(in);
                
                map<string, string> temp;
                temp[thiscolumn] = thisname;
                distName.push_back(temp);
            }
        }
        
        in.close();
        
        return singleton;
    }
    catch(exception& e) {
		m->errorOut(e, "ClusterSplitCommand", "readFile");
		exit(1);
	}
    
}
//**********************************************************************************************************************
int ClusterSplitCommand::getLabels(string file, set<string>& listLabels){
    try {
        ifstream in;
        util.openInputFile(file, in);

        //read headers
        util.getline(in); util.gobble(in);
        
        string label;
        while(!in.eof()) {
            if (m->getControl_pressed()) { break; }
            
            in >> label; util.getline(in); util.gobble(in);
            
            listLabels.insert(label);
        }
        
        in.close();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterSplitCommand", "getLabels");
        exit(1);
    }
    
}
//**********************************************************************************************************************
bool ClusterSplitCommand::findVsearch(){
    try {
        
        abort = false;
        
        if (cutoffNotSet) {  m->mothurOut("\nYou did not set a cutoff, using 0.03.\n"); cutoff = 0.03; }
        
        //look for vsearch exe

        string path = current->getProgramPath();
      
        string vsearchCommand = path + PATH_SEPARATOR;
        vsearchCommand += "vsearch";  vsearchCommand += EXECUTABLE_EXT;
#if defined NON_WINDOWS
        if (m->getDebug()) {
            m->mothurOut("[DEBUG]: vsearch location using \"which vsearch\" = ");
            Command* newCommand = new SystemCommand("which vsearch"); m->mothurOutEndLine();
            newCommand->execute(); delete newCommand;
            m->mothurOut("[DEBUG]: Mothur's location using \"which mothur\" = ");
            newCommand = new SystemCommand("which mothur"); m->mothurOutEndLine();
            newCommand->execute(); delete newCommand;
        }
#endif
        
        //test to make sure vsearch exists
        ifstream in;
        vsearchCommand = util.getFullPathName(vsearchCommand);
        bool ableToOpen = util.openInputFile(vsearchCommand, in, "no error"); in.close();
        if(!ableToOpen) {
            m->mothurOut(vsearchCommand + " file does not exist. Checking path... \n");
            
            ifstream in2;
            string programName = "vsearch"; programName += EXECUTABLE_EXT;
            string uLocation = util.findProgramPath(programName);
            ableToOpen = util.openInputFile(uLocation+programName, in2, "no error"); in2.close();
            
            if(!ableToOpen) { m->mothurOut("[ERROR]: " + uLocation + " file does not exist. mothur requires the vsearch executable.\n");  m->setControl_pressed(true);  }
            else {  m->mothurOut("Found vsearch in your path, using " + uLocation + "\n");vsearchLocation = uLocation+programName; }
        }else {  vsearchLocation = vsearchCommand; }
        
        vsearchLocation = util.getFullPathName(vsearchLocation);
        
        if (m->getDebug()) { m->mothurOut("[DEBUG]: vsearch location using " + vsearchLocation + "\n"); }
        
        if (!abort) { return true; }
        
        return false;

    }
    catch(exception& e) {
        m->errorOut(e, "ClusterSplitCommand", "findVsearch");
        exit(1);
    }
    
}
//**********************************************************************************************************************
