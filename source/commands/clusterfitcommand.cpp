//
//  clusterfitcommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 1/22/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "clusterfitcommand.hpp"
#include "readphylip.h"
#include "readcolumn.h"
#include "readmatrix.hpp"
#include "sequence.hpp"
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
vector<string> ClusterFitCommand::setParameters(){
    try {
        CommandParameter plist("reflist", "InputTypes", "", "", "", "", "","",false,true,true); parameters.push_back(plist);
        CommandParameter pfasta("fasta", "InputTypes", "", "", "", "", "","list",false,true,true); parameters.push_back(pfasta);
        CommandParameter prepfasta("reffasta", "InputTypes", "", "", "", "", "","",false,true,true); parameters.push_back(prepfasta);
        //CommandParameter pfastamap("fastamap", "InputTypes", "", "", "", "", "","",false,true,true); parameters.push_back(pfastamap);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none","","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount", "none", "","",false,false,true); parameters.push_back(pcount);
        CommandParameter prefname("refname", "InputTypes", "", "", "RefNameCount", "none","","",false,false,true); parameters.push_back(prefname);
        CommandParameter prefcount("refcount", "InputTypes", "", "", "RefNameCount", "none", "","",false,false,true); parameters.push_back(prefcount);
        CommandParameter prefcolumn("refcolumn", "InputTypes", "", "", "PhylipColumnRef", "", "ColumnName","",false,false,true); parameters.push_back(prefcolumn);
        CommandParameter prefphylip("refphylip", "InputTypes", "", "", "PhylipColumnRef", "", "","",false,false,true); parameters.push_back(prefphylip);
        CommandParameter pcolumn("column", "InputTypes", "", "", "PhylipColumn", "", "ColumnName","",false,false,true); parameters.push_back(pcolumn);
        CommandParameter pcutoff("cutoff", "Number", "", "0.03", "", "", "","",false,false,true); parameters.push_back(pcutoff);
        CommandParameter pprecision("precision", "Number", "", "100", "", "", "","",false,false); parameters.push_back(pprecision);
        CommandParameter pmethod("method", "Multiple", "closed-open", "closed", "", "", "","",false,false,true); parameters.push_back(pmethod);
        CommandParameter pinitialize("initialize", "Multiple", "oneotu-singleton", "singleton", "", "", "","",false,false,true); parameters.push_back(pinitialize);
        CommandParameter pmetric("metric", "Multiple", "mcc-sens-spec-tptn-fpfn-tp-tn-fp-fn-f1score-accuracy-ppv-npv-fdr", "mcc", "", "", "","",false,false,true); parameters.push_back(pmetric);
        CommandParameter pmetriccutoff("delta", "Number", "", "0.0001", "", "", "","",false,false,true); parameters.push_back(pmetriccutoff);
        CommandParameter piters("iters", "Number", "", "100", "", "", "","",false,false,true); parameters.push_back(piters);
        CommandParameter pdenovoiters("denovoiters", "Number", "", "100", "", "", "","",false,false,true); parameters.push_back(pdenovoiters);
        CommandParameter pfitpercent("fitpercent", "Number", "", "10", "", "", "","",false,false,true); parameters.push_back(pfitpercent);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        vector<string> myArray;
        for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
        return myArray;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterFitCommand", "setParameters");
        exit(1);
    }
}
//**********************************************************************************************************************
string ClusterFitCommand::getHelpString(){
    try {
        string helpString = "";
        helpString += "The cluster.fit command parameter options are reflist, refcolumn, refphylip, refname, refcount, fasta, name, count, column, method, cutoff, precent, metric, iters, initialize, denovoiters.\n";
        helpString += "The refcolumn parameter allow you to enter your reference data distance file, to reduce processing time. \n";
        helpString += "The refphylip parameter allow you to enter your reference data distance file, to reduce processing time. \n";
        helpString += "The column parameter allow you to enter your data distance file, to reduce processing time. \n";
        helpString += "The fasta parameter allows you to enter your fasta file. \n";
        helpString += "The reffasta parameter allows you to enter your fasta file for your reference dataset. \n";
        helpString += "The reflist parameter allows you to enter your list file for your reference dataset. \n";
        helpString += "The name parameter allows you to enter your name file. \n";
        helpString += "The count parameter allows you to enter your count file.\nA count or name file is required if your distance file is in column format.\n";
        helpString += "The refname parameter allows you to enter your reference name file. \n";
        helpString += "The refcount parameter allows you to enter your reference count file.\nA refcount or refname file is required if your reference distance file is in column format.\n";
        helpString += "The iters parameter allow you to set the maxiters for the opticluster method. \n";
        helpString += "The denovoiters parameter allow you to set the number of randomizations to perform. \n";
        helpString += "The fitpercent parameter allow you to set percentage of reads to be fitted. Default=10. Max=100, min=0.01.\n";
        helpString += "The metric parameter allows to select the metric in the opticluster method. Options are Matthews correlation coefficient (mcc), sensitivity (sens), specificity (spec), true positives + true negatives (tptn), false positives + false negatives (fpfn), true positives (tp), true negative (tn), false positive (fp), false negative (fn), f1score (f1score), accuracy (accuracy), positive predictive value (ppv), negative predictive value (npv), false discovery rate (fdr). Default=mcc.\n";
        helpString += "The initialize parameter allows to select the initial randomization for the opticluster method. Options are singleton, meaning each sequence is randomly assigned to its own OTU, or oneotu meaning all sequences are assigned to one otu. Default=singleton.\n";
        helpString += "The delta parameter allows to set the stable value for the metric in the opticluster method (delta=0.0001). \n";
        helpString += "The method parameter allows you to enter your clustering method. Options are closed and open. Default=closed.\n";
        helpString += "The cluster.fit command should be in the following format: \n";
        helpString += "cluster.fit(list=yourreflist, reffasta=yourReferenceFasta, fasta=yourFastaFile, count=yourCountFile) \n";
        return helpString;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterFitCommand", "getHelpString");
        exit(1);
    }
}
//**********************************************************************************************************************
string ClusterFitCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "list") {  pattern = "[filename],[clustertag],list-[filename],[clustertag],[tag2],list"; }
        else if (type == "sensspec") {  pattern = "[filename],sensspec"; }
        else if (type == "steps") {  pattern = "[filename],[clustertag],steps"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterFitCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ClusterFitCommand::ClusterFitCommand(){
    try {
        abort = true; calledHelp = true;
        setParameters();
        vector<string> tempOutNames;
        outputTypes["list"] = tempOutNames;
        outputTypes["sensspec"] = tempOutNames;
        outputTypes["steps"] = tempOutNames;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterFitCommand", "ClusterFitCommand");
        exit(1);
    }
}
//**********************************************************************************************************************
//This function checks to make sure the cluster command has no errors and then clusters based on the method chosen.
ClusterFitCommand::ClusterFitCommand(string option)  {
    try{
        abort = false; calledHelp = false;
        
        //allow user to run help
        if(option == "help") { help(); abort = true; calledHelp = true; }
        else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        
        else {
            vector<string> myArray = setParameters();
            
            OptionParser parser(option);
            map<string,string> parameters = parser.getParameters();
            map<string,string>::iterator it;
            
            ValidParameters validParameter;
            
            //check to make sure all parameters are valid for command
            for (it = parameters.begin(); it != parameters.end(); it++) {
                if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {
                    abort = true;
                }
            }
            
            //initialize outputTypes
            vector<string> tempOutNames;
            outputTypes["list"] = tempOutNames;
            outputTypes["sensspec"] = tempOutNames;
            outputTypes["steps"] = tempOutNames;
            
            //if the user changes the output directory command factory will send this info to us in the output parameter
            outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";		}
            
            inputDir = validParameter.valid(parameters, "inputdir");
            if (inputDir == "not found"){	inputDir = "";		}
            else {
                string path;
                
                it = parameters.find("refcolumn");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["refcolumn"] = inputDir + it->second;		}
                }
                
                it = parameters.find("refphylip");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["refphylip"] = inputDir + it->second;		}
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
                
                it = parameters.find("count");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["count"] = inputDir + it->second;		}
                }
                
                it = parameters.find("fasta");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
                }
                
                it = parameters.find("reflist");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["reflist"] = inputDir + it->second;		}
                }
                
                it = parameters.find("refname");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["refname"] = inputDir + it->second;		}
                }
                
                it = parameters.find("refcount");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["refcount"] = inputDir + it->second;		}
                }
                
                it = parameters.find("reffasta");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["reffasta"] = inputDir + it->second;		}
                }
            }
            
            selfReference = true;
            refdistfile = "";
            distfile = "";
            
            //check for required parameters
            reffastafile = validParameter.validFile(parameters, "reffasta");
            if (reffastafile == "not open") { abort = true; }
            else if (reffastafile == "not found") { reffastafile = "";  }
            else {  selfReference = false; }
            
            refcolumnfile = validParameter.validFile(parameters, "refcolumn");
            if (refcolumnfile == "not open") { refcolumnfile = ""; abort = true; }
            else if (refcolumnfile == "not found") { refcolumnfile = ""; }
            else {  refdistfile = refcolumnfile; refformat = "column"; current->setColumnFile(refcolumnfile);	selfReference = false; }
            
            refphylipfile = validParameter.validFile(parameters, "refphylip");
            if (refphylipfile == "not open") { refphylipfile = ""; abort = true; }
            else if (refphylipfile == "not found") { refphylipfile = ""; }
            else {  refdistfile = refphylipfile; refformat = "phylip"; current->setPhylipFile(refphylipfile);	selfReference = false; }
            
            reflistfile = validParameter.validFile(parameters, "reflist");
            if (reflistfile == "not open") { abort = true; }
            else if (reflistfile == "not found") { reflistfile = ""; }
            else { selfReference = false; }
            
            refnamefile = validParameter.validFile(parameters, "refname");
            if (refnamefile == "not open") { abort = true; }
            else if (refnamefile == "not found") { refnamefile = ""; }
            else { selfReference = false; }
            
            refcountfile = validParameter.validFile(parameters, "refcount");
            if (refcountfile == "not open") { abort = true;  }
            else if (refcountfile == "not found") { refcountfile = ""; }
            else { selfReference = false; }
            
            
            if (!selfReference) { //if you are providing reference files, lets make sure we have all of them
                if ((refdistfile == "") || (reffastafile == "") || (reflistfile == "")) { m->mothurOut("[ERROR]: When providing a reference file, you must provide a reffasta, refcolumn or refphylip, reflist and refcount or refname, aborting.\n");  abort = true; }
            }
            
            fastafile = validParameter.validFile(parameters, "fasta");
            if (fastafile == "not open") { abort = true; }
            else if (fastafile == "not found") { //if there is a current fasta file, use it
                if (!selfReference) {
                    fastafile = current->getFastaFile();
                    if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter.\n");  }
                    else { 	m->mothurOut("[ERROR]: You have no current fastafile and the fasta parameter is required.\n");  abort = true; }
                }else { fastafile = ""; }
            }else { current->setFastaFile(fastafile); }
            
            namefile = validParameter.validFile(parameters, "name");
            if (namefile == "not open") { abort = true; }
            else if (namefile == "not found") { namefile = ""; }
            else { current->setNameFile(namefile); }
            
            countfile = validParameter.validFile(parameters, "count");
            if (countfile == "not open") { abort = true; countfile = ""; }
            else if (countfile == "not found") { countfile = ""; }
            else { current->setCountFile(countfile); }
            
            columnfile = validParameter.validFile(parameters, "column");
            if (columnfile == "not open") { columnfile = ""; abort = true; }
            else if (columnfile == "not found") { columnfile = ""; }
            else {  distfile = columnfile;  current->setColumnFile(columnfile);	}
            
            method = validParameter.valid(parameters, "method");
            if (method == "not found") {  method = "open";}
            
            if ((method == "closed") || (method == "open")) { }
            else { m->mothurOut("[ERROR]: " + method + " is not a valid cluster fitting method.  Valid options are closed and open.\n"); abort = true; }
            
            if ((countfile != "") && (namefile != "")) { m->mothurOut("When executing a cluster.fit command you must enter ONLY ONE of the following: count or name.\n"); abort = true; }
            
            if (!selfReference) {
                if ((columnfile == "") && (fastafile == "")) {
                    //is there are current file available for either of these?
                    //give priority to column, then phylip
                    columnfile = current->getColumnFile();
                    if (columnfile != "") {  distfile = columnfile;  m->mothurOut("Using " + columnfile + " as input file for the column parameter.\n");  }
                    else {
                        fastafile = current->getFastaFile();
                        if (fastafile != "") {  distfile = fastafile;  m->mothurOut("Using " + fastafile + " as input file for the fasta parameter.\n");  }
                        else {
                            m->mothurOut("No valid current files. You must column or fasta file before you can use the cluster.fit command.\n");
                            abort = true;
                        }
                    }
                }
            }else {
                if (columnfile == "") {
                    //is there are current file available for either of these?
                    columnfile = current->getColumnFile();
                    if (columnfile != "") {  distfile = columnfile;  m->mothurOut("Using " + columnfile + " as input file for the column parameter.\n");  }
                    else {
                        m->mothurOut("No valid current files. You must provide a column file before you can use the cluster.fit command.\n");
                        abort = true;
                    }
                }
            }
            
            if (columnfile != "") {
                if ((namefile == "") && (countfile == "")) {
                    namefile = current->getNameFile();
                    if (namefile != "") {  m->mothurOut("Using " + namefile + " as input file for the name parameter.\n"); }
                    else {
                        countfile = current->getCountFile();
                        if (countfile != "") {  m->mothurOut("Using " + countfile + " as input file for the count parameter.\n"); }
                        else {  m->mothurOut("[ERROR]: You need to provide a namefile or countfile if you are going to use the column format.\n");  abort = true; }
                    }
                }
            }
            
            string temp = validParameter.valid(parameters, "precision");
            if (temp == "not found") { temp = "100"; }
            length = temp.length(); ////saves precision length for formatting below
            util.mothurConvert(temp, precision);
            
            temp = validParameter.valid(parameters, "delta");		if (temp == "not found")  { temp = "0.0001"; }
            util.mothurConvert(temp, stableMetric);
            
            metricName = validParameter.valid(parameters, "metric");		if (metricName == "not found") { metricName = "mcc"; }
            
            if ((metricName == "mcc") || (metricName == "sens") || (metricName == "spec") || (metricName == "tptn") || (metricName == "tp") || (metricName == "tn") || (metricName == "fp") || (metricName == "fn") || (metricName == "f1score") || (metricName == "accuracy") || (metricName == "ppv") || (metricName == "npv") || (metricName == "fdr") || (metricName == "fpfn") ){ }
            else { m->mothurOut("[ERROR]: Not a valid metric.  Valid metrics are mcc, sens, spec, tp, tn, fp, fn, tptn, fpfn, f1score, accuracy, ppv, npv, fdr."); m->mothurOutEndLine(); abort = true; }
            
            initialize = validParameter.valid(parameters, "initialize");		if (initialize == "not found") { initialize = "singleton"; }
            
            if (initialize == "singleton"){ }
            else { m->mothurOut("[ERROR]: Not a valid initialization.  Valid initialization is singleton.\n");  abort = true; }
            
            temp = validParameter.valid(parameters, "iters");		if (temp == "not found")  { temp = "100"; }
            util.mothurConvert(temp, maxIters);
            
            temp = validParameter.valid(parameters, "denovoiters");
            if (temp == "not found")  {
                if (selfReference) { temp = "10"; }
                else { temp = "1";  }
            }
            util.mothurConvert(temp, denovoIters);
            
            temp = validParameter.valid(parameters, "fitpercent");		if (temp == "not found")  { temp = "50.0"; }
            util.mothurConvert(temp, fitPercent);
            
            if ((fitPercent > 100) || (fitPercent < 0.01)) { abort=true; m->mothurOut("[ERROR]: fitpercent must be less than 100, and more than 0.01.\n"); }
            
            adjust=-1.0;
            temp = validParameter.valid(parameters, "cutoff");
            if (temp == "not found") { temp = "0.03"; }
            util.mothurConvert(temp, cutoff);
            
        }
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterFitCommand", "ClusterFitCommand");
        exit(1);
    }
}
//**********************************************************************************************************************
ClusterFitCommand::~ClusterFitCommand(){}
//**********************************************************************************************************************
int ClusterFitCommand::execute(){
    try {
        
        if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        time_t estart = time(NULL);
        
        ClusterMetric* metric = NULL;
        if (metricName == "mcc")             { metric = new MCC();              }
        else if (metricName == "sens")       { metric = new Sensitivity();      }
        else if (metricName == "spec")       { metric = new Specificity();      }
        else if (metricName == "tptn")       { metric = new TPTN();             }
        else if (metricName == "tp")         { metric = new TP();               }
        else if (metricName == "tn")         { metric = new TN();               }
        else if (metricName == "fp")         { metric = new FP();               }
        else if (metricName == "fn")         { metric = new FN();               }
        else if (metricName == "f1score")    { metric = new F1Score();          }
        else if (metricName == "accuracy")   { metric = new Accuracy();         }
        else if (metricName == "ppv")        { metric = new PPV();              }
        else if (metricName == "npv")        { metric = new NPV();              }
        else if (metricName == "fdr")        { metric = new FDR();              }
        else if (metricName == "fpfn")       { metric = new FPFN();             }
    
        map<string, int> counts;
        string dupsFile = countfile; nameOrCount = "count";
        if (namefile != "") { dupsFile = namefile; nameOrCount = "name"; }
        else { CountTable ct; ct.readTable(countfile, false, false); counts = ct.getNameMap();  }
        
        if (outputDir == "") { outputDir += util.hasPath(distfile); }
        fileroot = outputDir + util.getRootName(util.getSimpleName(distfile));
        
        map<string, string> variables;
        variables["[filename]"] = fileroot;
        variables["[clustertag]"] = "optifit_" + metric->getName();
        string outputName = getOutputFileName("steps", variables);
        string listFile = "";
        
        if (selfReference) {
        
            //distfile, distFormat, dupsFile, dupsFormat, cutoff, percentage to be fitseqs - will randomly assign as fit
            OptiData* matrix = new OptiRefMatrix(distfile, "column", dupsFile, nameOrCount, cutoff, fitPercent);
            
            listFile = runDenovoOptiCluster(matrix, metric, counts, outputName);
            
        }else {
            
            createReferenceNameCount(); //creates reference name or count file if needed
            
            calcDists();  //calc distance matrix for fasta file and distances between fasta file and reffasta file
            
            //calc sens.spec values for reference
            InputData input(reflistfile, "list", nullVector);
            ListVector* list = input.getListVector();
            
            //add tag to OTULabels to indicate the reference
            vector<string> refListLabels = list->getLabels();
            for (int i = 0; i < refListLabels.size(); i++) { refListLabels[i] = "Ref_" + refListLabels[i];  }
            list->setLabels(refListLabels);
            
            string refDupsFile = refcountfile;
            if (refNameOrCount == "name") { refDupsFile = refnamefile; }
            
            OptiData* matrix = new OptiRefMatrix(refdistfile, refDupsFile, refNameOrCount, refformat, cutoff, distfile, dupsFile, nameOrCount, "column", comboDistFile, "column");
        
            listFile = runRefOptiCluster(matrix, metric, list, counts, outputName);
            listFiles.push_back(listFile);
        }
    
        if (m->getControl_pressed()) { 	for (int j = 0; j < outputNames.size(); j++) { util.mothurRemove(outputNames[j]); }  return 0; }
        
        //evaluate results
        string bestListFileName = runSensSpec(listFile, columnfile, dupsFile, nameOrCount, metric, listFile);
        
        delete metric;
        outputNames.push_back(outputName); outputTypes["steps"].push_back(outputName);
        outputNames.push_back(bestListFileName); outputTypes["list"].push_back(bestListFileName);

        m->mothurOut("It took " + toString(time(NULL) - estart) + " seconds to fit sequences to reference OTUs.\n");
        
        //set list file as new current listfile
        string currentName = "";
        itTypes = outputTypes.find("list");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setListFile(currentName); }
        }
        
        m->mothurOut("\nOutput File Names: \n");
        for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]+"\n"); 	}
        m->mothurOutEndLine();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterFitCommand", "execute");
        exit(1);
    }
}
//**********************************************************************************************************************
string ClusterFitCommand::runDenovoOptiCluster(OptiData*& matrix, ClusterMetric*& metric, map<string, int>& counts, string outStepFile){
    try {
        m->mothurOut("\nClustering " + distfile + "\n");
        bool printStepsHeader = true;
        
        for (int i = 0; i < denovoIters; i++) {
            
            OptiFitCluster cluster(matrix, metric, 0);
            tag = cluster.getTag();
            
            int iters = 0;
            double listVectorMetric = 0; //worst state
            double delta = 1;
            
            //get "ref" seqs for initialize inputs
            OptiData* refMatrix = matrix->extractRefMatrix();
            
            ListVector* refList = clusterRefs(refMatrix, metric);
            
            delete refMatrix;
            
            vector<vector<string> > otus;
            for (int i = 0; i < refList->getNumBins(); i++) {
                vector<string> binNames;
                string bin = refList->get(i);
                if (bin != "") {
                    util.splitAtComma(bin, binNames);
                    otus.push_back(binNames);
                }
            }
            
            //add tag to OTULabels to indicate the reference
            vector<string> refListLabels = refList->getLabels();
            for (int i = 0; i < refListLabels.size(); i++) { refListLabels[i] = "Ref_" + refListLabels[i];  }
            refList->setLabels(refListLabels);
            
            cluster.initialize(listVectorMetric, true, otus, refList->getLabels(), method, true);
            
            delete refList;
            
            long long numBins = cluster.getNumBins();
            long long tp, tn, fp, fn;
            vector<double> results = cluster.getStats(tp, tn, fp, fn);
            
            long long fittp, fittn, fitfp, fitfn;
            long long numFitBins = cluster.getNumFitBins();
            vector<double> fitresults = cluster.getFitStats(fittp, fittn, fitfp, fitfn);
            
            m->mothurOut("\nFitting " + toString(matrix->getNumFitSeqs()+matrix->getNumFitSingletons()+matrix->getNumFitTrueSingletons()) + " sequences to reference otus.\n");
            
            m->mothurOut("\n\nlist\tstate\titer\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");  
            
            outputSteps(outStepFile, printStepsHeader, tp, tn, fp, fn, results, numBins, fittp, fittn, fitfp, fitfn, fitresults, numFitBins, 0, false, 0);
            
            while ((delta > stableMetric) && (iters < maxIters)) { //
                
                if (m->getControl_pressed()) { break; }
                double oldMetric = listVectorMetric;
                
                cluster.update(listVectorMetric);
                
                delta = abs(oldMetric - listVectorMetric);
                iters++;
                
                results = cluster.getStats(tp, tn, fp, fn);
                numBins = cluster.getNumBins();
                numFitBins = cluster.getNumFitBins();
                fitresults = cluster.getFitStats(fittp, fittn, fitfp, fitfn);
                
                outputSteps(outStepFile, printStepsHeader, tp, tn, fp, fn, results, numBins, fittp, fittn, fitfp, fitfn, fitresults, numFitBins, iters, false, i);
            }
            outputSteps(outStepFile, printStepsHeader, tp, tn, fp, fn, results, numBins, fittp, fittn, fitfp, fitfn, fitresults, numFitBins, iters, true, i);
            m->mothurOutEndLine(); m->mothurOutEndLine();
            
            if (m->getControl_pressed()) {  return 0; }
            
            ofstream listFile;
            tag = "optifit_" + metric->getName() + "_denovo." + toString(i+1);
            string listFileName = fileroot+ tag + ".list";
            util.openOutputFile(listFileName,	listFile);
            
            ListVector* list = cluster.getFittedList(toString(cutoff));
            list->setLabel(toString(cutoff));
            list->setLabels(nullVector);
            
            if(countfile != "") { list->print(listFile, counts); }
            else { list->print(listFile); }
            
            listFile.close();
            listFiles.push_back(listFileName);
            
            delete list;
            
            matrix->randomizeRefs();
        }
        
        tag = "optifit_" + metric->getName() + "_denovo";
        string listFileName = fileroot+ tag + ".list";
        
        return listFileName;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterFitCommand", "runDenovoOptiCluster");
        exit(1);
    }
    
}
/***********************************************************************/
ListVector* ClusterFitCommand::clusterRefs(OptiData*& refsMatrix, ClusterMetric*& metric) {
    try {
        m->mothurOut("\nClustering " + toString(refsMatrix->getNumSeqs()+refsMatrix->getNumSingletons()) + " reference sequences.\n");
        
        ListVector* list = NULL;
        
        OptiCluster cluster(refsMatrix, metric, 0);
        
        int iters = 0;
        double listVectorMetric = 0; //worst state
        double delta = 1;
        
        cluster.initialize(listVectorMetric, true, "singleton");
        
        long long numBins = cluster.getNumBins();
        m->mothurOut("\n\niter\ttime\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");
        
        long long tp, tn, fp, fn;
        vector<double> results = cluster.getStats(tp, tn, fp, fn);
        m->mothurOut("0\t0\t" + toString(cutoff) + "\t" + toString(numBins) + "\t"+ toString(cutoff) + "\t" + toString(tp) + "\t" + toString(tn) + "\t" + toString(fp) + "\t" + toString(fn) + "\t");
        
        for (int i = 0; i < results.size(); i++) { m->mothurOut(toString(results[i]) + "\t");  }
        m->mothurOutEndLine();
        
        while ((delta > 0.0001) && (iters < maxIters)) {
            
            long start = time(NULL);
            
            if (m->getControl_pressed()) { break; }
            double oldMetric = listVectorMetric;
            
            cluster.update(listVectorMetric);
            
            delta = abs(oldMetric - listVectorMetric);
            iters++;
            
            results = cluster.getStats(tp, tn, fp, fn);
            numBins = cluster.getNumBins();
            
            m->mothurOut(toString(iters) + "\t" + toString(time(NULL) - start) + "\t" + toString(cutoff) + "\t" + toString(numBins) + "\t" + toString(cutoff) + "\t"+ toString(tp) + "\t" + toString(tn) + "\t" + toString(fp) + "\t" + toString(fn) + "\t");
            
            for (int i = 0; i < results.size(); i++) { m->mothurOut(toString(results[i]) + "\t");  }
            m->mothurOutEndLine();
            
        }
        m->mothurOutEndLine(); m->mothurOutEndLine();
        
        if (m->getControl_pressed()) { return list; }
        
        list = cluster.getList();
        list->setLabel(toString(cutoff));
        
        return list;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiFitCluster", "clusterRefs");
        exit(1);
    }
}
//**********************************************************************************************************************
string ClusterFitCommand::runRefOptiCluster(OptiData*& matrix, ClusterMetric*& metric, ListVector*& refList, map<string, int>& counts, string outStepFile){
    try {
        OptiFitCluster cluster(matrix, metric, 0);
        tag = cluster.getTag();
        
        m->mothurOut("\nClustering " + distfile + "\n"); 
        
        int iters = 0;
        double listVectorMetric = 0; //worst state
        double delta = 1;
        
        vector<vector<string> > otus;
        for (int i = 0; i < refList->getNumBins(); i++) {
            vector<string> binNames;
            string bin = refList->get(i);
            if (bin != "") {
                util.splitAtComma(bin, binNames);
                otus.push_back(binNames);
            }
        }
        
        cluster.initialize(listVectorMetric, true, otus, refList->getLabels(), method, false);
        
        long long numBins = cluster.getNumBins();
        long long tp, tn, fp, fn;
        vector<double> results = cluster.getStats(tp, tn, fp, fn);
        
        long long fittp, fittn, fitfp, fitfn;
        long long numFitBins = cluster.getNumFitBins();
        vector<double> fitresults = cluster.getFitStats(fittp, fittn, fitfp, fitfn);
        
        bool printStepsHeader = true;
        outputSteps(outStepFile, printStepsHeader, tp, tn, fp, fn, results, numBins, fittp, fittn, fitfp, fitfn, fitresults, numFitBins, 0, true, 0);
        
        while ((delta > stableMetric) && (iters < maxIters)) { //
            
            if (m->getControl_pressed()) { break; }
            double oldMetric = listVectorMetric;
            
            cluster.update(listVectorMetric);
            
            delta = abs(oldMetric - listVectorMetric);
            iters++;
            
            results = cluster.getStats(tp, tn, fp, fn);
            numBins = cluster.getNumBins();
            numFitBins = cluster.getNumFitBins();
            fitresults = cluster.getFitStats(fittp, fittn, fitfp, fitfn);
            
            outputSteps(outStepFile, printStepsHeader, tp, tn, fp, fn, results, numBins, fittp, fittn, fitfp, fitfn, fitresults, numFitBins, iters, true, 0);
        }
        m->mothurOutEndLine(); m->mothurOutEndLine();
        
        if (m->getControl_pressed()) {  return 0; }
        
        ListVector* list = cluster.getFittedList(toString(cutoff));
        list->setLabel(toString(cutoff));
        
        ofstream listFile;
        string listFileName = fileroot+ tag + ".list";
        util.openOutputFile(listFileName,	listFile);
        outputNames.push_back(listFileName); outputTypes["list"].push_back(listFileName);

        if(countfile != "") { list->print(listFile, counts); }
        else { list->print(listFile); }
        listFile.close();
        
        delete list;

        return listFileName;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterFitCommand", "runRefOptiCluster");
        exit(1);
    }
    
}
//**********************************************************************************************************************
string ClusterFitCommand::runSensSpec(string listFileName, string distFName, string dupsfile, string dupsFormat, ClusterMetric*& userMetric, string listFile) {
    try {
        
        ofstream sensSpecFile;
        map<string, string> variables;
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(listFile));
        string sensSpecFileName = getOutputFileName("sensspec",variables);
        util.openOutputFile(sensSpecFileName, sensSpecFile);
        outputNames.push_back(sensSpecFileName); outputTypes["sensspec"].push_back(sensSpecFileName);
        
        sensSpecFile << "iter\tlabel\tcutoff\tnumotus\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";
        m->mothurOut("iter\tlabel\tcutoff\tnumotus\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");
        
        double bestStat = 0; int bestResult = 0;
        
        for (int i = 0; i < listFiles.size(); i++) {
            
            string thislistFileName = listFiles[i];
            string thisDistFile = distFName;
            string thisDupsFile = dupsfile;
            
            if (method != "open") {
                //extract only distances related to the list file
                string inputString = "list=" + thislistFileName;
                m->mothurOut("/******************************************/\n");
                m->mothurOut("Running command: list.seqs(" + inputString + ")\n");
                current->setMothurCalling(true);
                
                Command* listSeqsCommand = new ListSeqsCommand(inputString);
                listSeqsCommand->execute();
                
                map<string, vector<string> > filenames = listSeqsCommand->getOutputFiles();
                
                delete listSeqsCommand;
                current->setMothurCalling(false);
                
                string accnosFileName = filenames["accnos"][0];
                
                inputString = "column=" + thisDistFile + ", accnos=" + accnosFileName;
                m->mothurOut("\n/***** NOTE: Please ignore warnings for get.dists command *****/\n");
                m->mothurOut("Running command: get.dists(" + inputString + ")\n");
                current->setMothurCalling(true);
                
                Command* getDistsCommand = new GetDistsCommand(inputString);
                getDistsCommand->execute();
                
                filenames = getDistsCommand->getOutputFiles();
                
                delete getDistsCommand;
                current->setMothurCalling(false);
                
                thisDistFile = filenames["column"][0];
                
                inputString = "accnos=" + accnosFileName;
                inputString += ", " + dupsFormat + "=" + thisDupsFile;
                
                m->mothurOut("/nRunning command: get.seqs(" + inputString + ")\n");
                current->setMothurCalling(true);
                
                Command* getSeqsCommand = new GetSeqsCommand(inputString);
                getSeqsCommand->execute();
                
                filenames = getSeqsCommand->getOutputFiles();
                
                if (dupsFormat == "name")         {  thisDupsFile = filenames["name"][0];       }
                else if (dupsFormat == "count")   { thisDupsFile = filenames["count"][0];       }
                
                util.mothurRemove(accnosFileName);
                
                delete getSeqsCommand;
                current->setMothurCalling(false);
            }
            
            InputData input(thislistFileName, "list", nullVector);
            ListVector* list = input.getListVector();

            string label = list->getLabel();
            int numBins = list->getNumBins();
        
            OptiMatrix matrix(thisDistFile, thisDupsFile, dupsFormat, "column", cutoff, false);
            SensSpecCalc senscalc(matrix, list);
            long long truePositives, trueNegatives, falsePositives, falseNegatives;
            senscalc.getResults(matrix, truePositives, trueNegatives, falsePositives, falseNegatives);
            
            long long tp =  truePositives;
            long long fp =  falsePositives;
            long long tn =  trueNegatives;
            long long fn =  falseNegatives;
            
            Sensitivity sens;   double sensitivity = sens.getValue(tp, tn, fp, fn);
            Specificity spec;   double specificity = spec.getValue(tp, tn, fp, fn);
            PPV ppv;            double positivePredictiveValue = ppv.getValue(tp, tn, fp, fn);
            NPV npv;            double negativePredictiveValue = npv.getValue(tp, tn, fp, fn);
            FDR fdr;            double falseDiscoveryRate = fdr.getValue(tp, tn, fp, fn);
            Accuracy acc;       double accuracy = acc.getValue(tp, tn, fp, fn);
            MCC mcc;            double matthewsCorrCoef = mcc.getValue(tp, tn, fp, fn);
            F1Score f1;         double f1Score = f1.getValue(tp, tn, fp, fn);
            
            sensSpecFile << i+1 << '\t' << label << '\t' << cutoff << '\t' << numBins << '\t';
            sensSpecFile << truePositives << '\t' << trueNegatives << '\t' << falsePositives << '\t' << falseNegatives << '\t';
            sensSpecFile << setprecision(4);
            sensSpecFile << sensitivity << '\t' << specificity << '\t' << positivePredictiveValue << '\t' << negativePredictiveValue << '\t';
            sensSpecFile << falseDiscoveryRate << '\t' << accuracy << '\t' << matthewsCorrCoef << '\t' << f1Score << endl;
            
            m->mothurOut(toString(i+1) + "\t" + label + "\t" + toString(cutoff) + "\t" + toString(numBins) + "\t"+ toString(truePositives) + "\t" + toString(trueNegatives) + "\t" + toString(falsePositives) + "\t" + toString(falseNegatives) + "\t");
            m->mothurOut(toString(sensitivity) + "\t" + toString(specificity) + "\t" + toString(positivePredictiveValue) + "\t" + toString(negativePredictiveValue) + "\t");
            m->mothurOut(toString(falseDiscoveryRate) + "\t" + toString(accuracy) + "\t" + toString(matthewsCorrCoef) + "\t" + toString(f1Score) + "\n\n");
            
            double userStat = userMetric->getValue(tp, tn, fp, fn);
            if (userStat > bestStat) { bestStat = userStat; bestResult = i; }
        }
        
        sensSpecFile.close();
        
        return listFiles[bestResult];
        
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterFitCommand", "runSensSpec");
        exit(1);
    }
}

//**********************************************************************************************************************
void ClusterFitCommand::outputSteps(string outputName, bool& printHeaders, long long tp, long long tn, long long fp, long long fn, vector<double> results, long long numBins, long long fittp, long long fittn, long long fitfp, long long fitfn, vector<double> fitresults, long long numFitBins, int iter, bool printToFile, int denovoIter) {
    try {

        if (!selfReference) { //writes to file as well
            if (printHeaders) {  m->mothurOut("\n\nstate\titer\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");  }
            
            m->mothurOut("combo\t" + toString(iter) + "\t" + toString(cutoff) + "\t" + toString(numBins) + "\t"+ toString(cutoff) + "\t" + toString(tp) + "\t" + toString(tn) + "\t" + toString(fp) + "\t" + toString(fn) + "\t");
            for (int i = 0; i < results.size(); i++) { m->mothurOut(toString(results[i]) + "\t");  } m->mothurOutEndLine();
            
            m->mothurOut("fit\t" + toString(iter) + "\t" + toString(cutoff) + "\t" + toString(numFitBins) + "\t"+ toString(cutoff) + "\t" + toString(fittp) + "\t" + toString(fittn) + "\t" + toString(fitfp) + "\t" + toString(fitfn) + "\t");
            for (int i = 0; i < fitresults.size(); i++) { m->mothurOut(toString(fitresults[i]) + "\t");  } m->mothurOutEndLine();
            
            ofstream outStep;
            if (printHeaders)   {
                util.openOutputFile(outputName, outStep);
                outStep << "state\titer\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";
                printHeaders = false;
            }else                { util.openOutputFileAppend(outputName, outStep);   }
            
            outStep << "combo\t" + toString(iter) + "\t" + toString(cutoff) + "\t" + toString(numBins) + "\t" + toString(cutoff) + "\t" << tp << '\t' << tn << '\t' << fp << '\t' << fn << '\t';
            for (int i = 0; i < results.size(); i++) {  outStep << results[i] << "\t"; } outStep << endl;
            
            outStep << "fit\t" + toString(iter) + "\t" + toString(cutoff) + "\t" + toString(numFitBins) + "\t" + toString(cutoff) + "\t" << fittp << '\t' << fittn << '\t' << fitfp << '\t' << fitfn << '\t';
            for (int i = 0; i < fitresults.size(); i++) {  outStep << fitresults[i] << "\t"; } outStep << endl;
        }else {
            //print results for each iter???
            if (printToFile) {
                ofstream outStep;
                if (printHeaders)   {
                    util.openOutputFile(outputName, outStep);
                    outStep << "list\t\tstate\titer\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";
                    printHeaders = false;
                }else                { util.openOutputFileAppend(outputName, outStep);   }
                
                outStep << toString(denovoIter+1) + "\tcombo\t" + toString(iter) + "\t" + toString(cutoff) + "\t" + toString(numBins) + "\t" + toString(cutoff) + "\t" << tp << '\t' << tn << '\t' << fp << '\t' << fn << '\t';
                for (int i = 0; i < results.size(); i++) {  outStep << results[i] << "\t"; } outStep << endl;
                
                outStep << toString(denovoIter+1) + "\tfit\t" + toString(iter) + "\t" + toString(cutoff) + "\t" + toString(numFitBins) + "\t" + toString(cutoff) + "\t" << fittp << '\t' << fittn << '\t' << fitfp << '\t' << fitfn << '\t';
                for (int i = 0; i < fitresults.size(); i++) {  outStep << fitresults[i] << "\t"; } outStep << endl;
            }else {
                m->mothurOut(toString(denovoIter+1) + "\t" + "combo\t" + toString(iter) + "\t" + toString(cutoff) + "\t" + toString(numBins) + "\t"+ toString(cutoff) + "\t" + toString(tp) + "\t" + toString(tn) + "\t" + toString(fp) + "\t" + toString(fn) + "\t");
                for (int i = 0; i < results.size(); i++) { m->mothurOut(toString(results[i]) + "\t");  } m->mothurOutEndLine();
                
                m->mothurOut(toString(denovoIter+1) + "\t" +"fit\t" + toString(iter) + "\t" + toString(cutoff) + "\t" + toString(numFitBins) + "\t"+ toString(cutoff) + "\t" + toString(fittp) + "\t" + toString(fittn) + "\t" + toString(fitfp) + "\t" + toString(fitfn) + "\t");
                for (int i = 0; i < fitresults.size(); i++) { m->mothurOut(toString(fitresults[i]) + "\t");  } m->mothurOutEndLine();
            }
        }
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterFitCommand", "outputSteps");
        exit(1);
    }
}
//**********************************************************************************************************************
void ClusterFitCommand::createReferenceNameCount() {
    try {
        
        if (refcountfile != "")     { refNameOrCount = "count";  }
        else if (refnamefile != "") { refNameOrCount = "name"; }
        else { //create count file
            current->setMothurCalling(true);
            string options = "fasta=" + reffastafile + ", format=count";
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Running command: unique.seqs(" + options + ")\n");
            Command* deconvoluteCommand = new DeconvoluteCommand(options);
            
            deconvoluteCommand->execute();
            map<string, vector<string> > filenames = deconvoluteCommand->getOutputFiles();
            refcountfile = filenames["count"][0];
            refNameOrCount = "count";
            
            delete deconvoluteCommand;
            m->mothurOut("/******************************************/\n");
        }
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterFitCommand", "createReferenceNameCount");
        exit(1);
    }
}
//**********************************************************************************************************************
string ClusterFitCommand::calcDists() {
    try {
        if (columnfile == "") { //calc user distances
            string options = "fasta=" + fastafile + ", cutoff=" + toString(cutoff);
            if (outputDir != "")                            { options += ", outputdir=" + outputDir;    }
            current->setMothurCalling(true);
            
            //calc dists for fastafile
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Running command: dist.seqs(" + options + ")\n");
            Command* distCommand = new DistanceCommand(options);
            
            distCommand->execute();
            map<string, vector<string> > filenames = distCommand->getOutputFiles();
            distfile = filenames["column"][0];
            columnfile = distfile;
            
            delete distCommand;
            m->mothurOut("/******************************************/\n");
        }
        
        map<string, vector<string> > filenames;
        int refAlignLength = util.getAlignmentLength(reffastafile);
        int alignLength = util.getAlignmentLength(fastafile);
        
        if (refAlignLength == alignLength) {
            string options = "fitcalc=t, fasta=" + reffastafile + ", oldfasta=" + fastafile + ", cutoff=" + toString(cutoff) + ", column=" + distfile;
            if (outputDir != "")                            { options += ", outputdir=" + outputDir;    }
            
            //dists between reffasta and fastafile
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Running command: dist.seqs(" + options + ")\n");
            Command* distCommand = new DistanceCommand(options);
            
            distCommand->execute();
            filenames = distCommand->getOutputFiles();
            comboDistFile = filenames["column"][0];
            
            delete distCommand;
            m->mothurOut("/******************************************/\n");
            current->setMothurCalling(false);
        }else {
            //filter each file to improve distance calc time
            string options = "fasta=" + reffastafile + ", vertical=t";
            if (outputDir != "")                            { options += ", outputdir=" + outputDir;    }
            
            m->mothurOut("\nRunning vertical filter to improve distance calculation time\n\n");
            
            //filter reffasta
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Running command: filter.seqs(" + options + ")\n");
            Command* filterCommand = new FilterSeqsCommand(options);
            
            filterCommand->execute();
            map<string, vector<string> > filenames = filterCommand->getOutputFiles();
            string filteredRef = filenames["fasta"][0];
            
            delete filterCommand;
            m->mothurOut("/******************************************/\n");
            
            options = "fasta=" + fastafile + ", reference=" + filteredRef;
            if (outputDir != "")                            { options += ", outputdir=" + outputDir;    }
            
            //align fasta to refFasta
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Running command: align.seqs(" + options + ")\n");
            Command* alignCommand = new AlignCommand(options);
            
            alignCommand->execute();
            filenames = alignCommand->getOutputFiles();
            string alignedFasta = filenames["fasta"][0];
            
            delete alignCommand;
            m->mothurOut("/******************************************/\n");
            
            options = "fitcalc=t, fasta=" + filteredRef + ", oldfasta=" + alignedFasta + ", cutoff=" + toString(cutoff) + ", column=" + distfile;
            if (outputDir != "")                            { options += ", outputdir=" + outputDir;    }
            
            //dists between reffasta and fastafile
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Running command: dist.seqs(" + options + ")\n");
            Command* distCommand = new DistanceCommand(options);
            
            distCommand->execute();
            filenames = distCommand->getOutputFiles();
            comboDistFile = filenames["column"][0];
            
            delete distCommand;
            m->mothurOut("/******************************************/\n");
            current->setMothurCalling(false);
        }
        
        return comboDistFile;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterFitCommand", "calcDists");
        exit(1);
    }
}
//**********************************************************************************************************************


