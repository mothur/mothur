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
        CommandParameter plist("list", "InputTypes", "", "", "", "", "","",false,true,true); parameters.push_back(plist);
        CommandParameter pfasta("fasta", "InputTypes", "", "", "", "", "","list",false,true,true); parameters.push_back(pfasta);
        CommandParameter prepfasta("reffasta", "InputTypes", "", "", "", "", "","",false,true,true); parameters.push_back(prepfasta);
        CommandParameter pfastamap("fastamap", "InputTypes", "", "", "", "", "","",false,true,true); parameters.push_back(pfastamap);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none","","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount", "none", "","",false,false,true); parameters.push_back(pcount);
        CommandParameter prefname("refname", "InputTypes", "", "", "RefNameCount", "none","","",false,false,true); parameters.push_back(prefname);
        CommandParameter prefcount("refcount", "InputTypes", "", "", "RefNameCount", "none", "","",false,false,true); parameters.push_back(prefcount);
        CommandParameter prefcolumn("refcolumn", "InputTypes", "", "", "", "", "ColumnName","",false,false,true); parameters.push_back(prefcolumn);
        CommandParameter pcutoff("cutoff", "Number", "", "0.03", "", "", "","",false,false,true); parameters.push_back(pcutoff);
        CommandParameter pprecision("precision", "Number", "", "100", "", "", "","",false,false); parameters.push_back(pprecision);
        CommandParameter pmethod("method", "Multiple", "opti", "opti", "", "", "","",false,false,true); parameters.push_back(pmethod);
        CommandParameter pinitialize("initialize", "Multiple", "oneotu-singleton", "singleton", "", "", "","",false,false,true); parameters.push_back(pinitialize);
        CommandParameter pmetric("metric", "Multiple", "mcc-sens-spec-tptn-fpfn-tp-tn-fp-fn-f1score-accuracy-ppv-npv-fdr", "mcc", "", "", "","",false,false,true); parameters.push_back(pmetric);
        CommandParameter pmetriccutoff("delta", "Number", "", "0.0001", "", "", "","",false,false,true); parameters.push_back(pmetriccutoff);
        CommandParameter piters("iters", "Number", "", "100", "", "", "","",false,false,true); parameters.push_back(piters);
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
        helpString += "The cluster.fit command parameter options are list, refcolumn, refname, refcount, fasta, name, count, method, cutoff, precision, metric, iters, initialize.\n";
        helpString += "The refcolumn parameter allow you to enter your reference data distance file, to reduce processing time. \n";
        helpString += "The fasta parameter allows you to enter your fasta file. \n";
        helpString += "The reffasta parameter allows you to enter your fasta file for your reference dataset. \n";
        helpString += "The list parameter allows you to enter your list file for your reference dataset. \n";
        helpString += "The name parameter allows you to enter your name file. \n";
        helpString += "The count parameter allows you to enter your count file.\nA count or name file is required if your distance file is in column format.\n";
        helpString += "The refname parameter allows you to enter your reference name file. \n";
        helpString += "The refcount parameter allows you to enter your reference count file.\nA refcount or refname file is required if your reference distance file is in column format.\n";
        helpString += "The iters parameter allow you to set the maxiters for the opticluster method. \n";
        helpString += "The metric parameter allows to select the metric in the opticluster method. Options are Matthews correlation coefficient (mcc), sensitivity (sens), specificity (spec), true positives + true negatives (tptn), false positives + false negatives (fpfn), true positives (tp), true negative (tn), false positive (fp), false negative (fn), f1score (f1score), accuracy (accuracy), positive predictive value (ppv), negative predictive value (npv), false discovery rate (fdr). Default=mcc.\n";
        helpString += "The initialize parameter allows to select the initial randomization for the opticluster method. Options are singleton, meaning each sequence is randomly assigned to its own OTU, or oneotu meaning all sequences are assigned to one otu. Default=singleton.\n";
        helpString += "The delta parameter allows to set the stable value for the metric in the opticluster method (delta=0.0001). \n";
        helpString += "The method parameter allows you to enter your clustering mothod. Options are opti. Default=opti.\n";
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
        else if (type == "sensspec") {  pattern = "[filename],[clustertag],sensspec"; }
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
        m->errorOut(e, "ClusterCommand", "ClusterCommand");
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
                
                it = parameters.find("list");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["list"] = inputDir + it->second;		}
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
                
                it = parameters.find("fastamap");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["fastamap"] = inputDir + it->second;		}
                }

            }
            
            //for testing
            fastaMapFile = validParameter.validFile(parameters, "fastamap");
            if (fastaMapFile == "not open") { abort = true; }
            else if (fastaMapFile == "not found") { fastaMapFile = ""; translate = false; }
            else { translate = true; }
            
            //check for required parameters
            reffastafile = validParameter.validFile(parameters, "reffasta");
            if (reffastafile == "not open") { abort = true; }
            else if (reffastafile == "not found") { m->mothurOut("[ERROR]: You must provide a reffasta file before you can use the cluster.fit command\n"); abort = true; reffastafile = ""; }
            else { distfile = "";  format="column";  }
            
            columnfile = validParameter.validFile(parameters, "refcolumn");
            if (columnfile == "not open") { columnfile = ""; abort = true; }
            else if (columnfile == "not found") { columnfile = ""; }
            else {  distfile = columnfile; format = "column"; current->setColumnFile(columnfile);	}
            
            reflistfile = validParameter.validFile(parameters, "list");
            if (reflistfile == "not open") { abort = true; }
            else if (reflistfile == "not found") { m->mothurOut("[ERROR]: You must provide a reference list file before you can use the cluster.fit command\n"); abort = true; reflistfile = ""; }
            
            refnamefile = validParameter.validFile(parameters, "refname");
            if (refnamefile == "not open") { abort = true; }
            else if (refnamefile == "not found") { refnamefile = ""; }
            
            refcountfile = validParameter.validFile(parameters, "refcount");
            if (refcountfile == "not open") { abort = true;  }
            else if (refcountfile == "not found") { refcountfile = ""; }
            
            fastafile = validParameter.validFile(parameters, "fasta");
            if (fastafile == "not open") { abort = true; }
            else if (fastafile == "not found") { //if there is a current fasta file, use it
                fastafile = current->getFastaFile();
                if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter.\n");  }
                else { 	m->mothurOut("[ERROR]: You have no current fastafile and the fasta parameter is required.\n");  abort = true; } }
            else { current->setFastaFile(fastafile); }
            
            namefile = validParameter.validFile(parameters, "name");
            if (namefile == "not open") { abort = true; }
            else if (namefile == "not found") { namefile = ""; }
            else { current->setNameFile(namefile); }
            
            countfile = validParameter.validFile(parameters, "count");
            if (countfile == "not open") { abort = true; countfile = ""; }
            else if (countfile == "not found") { countfile = ""; }
            else { current->setCountFile(countfile); }
            
            method = validParameter.valid(parameters, "method");
            if (method == "not found") {  method = "opti";}
            
            if (method == "opti") { }
            else { m->mothurOut("[ERROR]: " + method + " is not a valid cluster fitting method.  Valid algorithm is opti.\n"); abort = true; }
            
            if (columnfile != "") {
                if ((refnamefile == "") && (refcountfile == "")){
                    m->mothurOut("You need to provide a reference namefile or reference countfile if you are going to use the column format, ignoring distance matrix.\n"); columnfile = "";
                }
            }
            if ((countfile != "") && (namefile != "")) { m->mothurOut("When executing a cluster command you must enter ONLY ONE of the following: count or name.\n"); abort = true; }
            
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
            
            if ((initialize == "singleton") || (initialize == "oneotu")){ }
            else { m->mothurOut("[ERROR]: Not a valid initialization.  Valid initializations are singleton and oneotu."); m->mothurOutEndLine(); abort = true; }
            
            temp = validParameter.valid(parameters, "iters");		if (temp == "not found")  { temp = "100"; }
            util.mothurConvert(temp, maxIters);
            
            adjust=-1.0;
            temp = validParameter.valid(parameters, "cutoff");
            if (temp == "not found") { temp = "0.03"; }
            util.mothurConvert(temp, cutoff);
        }
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterCommand", "ClusterCommand");
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
        
        //calc sens.spec values for reference
        InputData input(reflistfile, "list", nullVector);
        ListVector* list = input.getListVector();
        
        calcReferenceValues(list); //creates distance file if needed

        //calc distance matrix for fasta file and distances between fasta file and reffasta file
        string newDistFile = calcDists(list);
        
        OptiMatrix matrix(newDistFile, combinedNameFile, combinedNameOrCount, "column", cutoff, false); util.mothurRemove(combinedNameFile);
        
        runOptiCluster(matrix, list);
        
        //create list file for new data
        
        if (m->getControl_pressed()) { 	for (int j = 0; j < outputNames.size(); j++) { util.mothurRemove(outputNames[j]); }  return 0; }
        
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
void ClusterFitCommand::calcReferenceValues(ListVector* list) {
    try {
        string nameOrCount = "";
        string thisNamefile = "";
        if (refcountfile != "")     { nameOrCount = "count"; thisNamefile = refcountfile;   }
        else if (refnamefile != "") { nameOrCount = "name"; thisNamefile = refnamefile;     }
        
        distfile = columnfile;
        if (distfile == "") {
            string options = "fasta=" + reffastafile + ", cutoff=" + toString(cutoff);
            if (outputDir != "")        { options += ", outputdir=" + outputDir;        }
            if (thisNamefile == "")     { options += ", output=lt"; format="phylip";    }else { format="column"; }
            
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Running command: dist.seqs(" + options + ")\n");
            Command* distCommand = new DistanceCommand(options);
            
            distCommand->execute();
            map<string, vector<string> > filenames = distCommand->getOutputFiles();
            distfile = filenames[format][0];
            
            delete distCommand;
            current->setMothurCalling(false);
            m->mothurOut("/******************************************/\n");
        }
        
        OptiMatrix matrix(distfile, thisNamefile, nameOrCount, format, cutoff, false); //matrix memory is released before the second matrix is read in
        
        SensSpecCalc senscalc(matrix, list);
        senscalc.getResults(matrix, truePositives, trueNegatives, falsePositives, falseNegatives);
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterFitCommand", "calcReferenceValues");
        exit(1);
    }
}
//**********************************************************************************************************************
string ClusterFitCommand::calcDists(ListVector* list) {
    try {
        string options = "fasta=" + fastafile + ", cutoff=" + toString(cutoff);
        if (outputDir != "")                            { options += ", outputdir=" + outputDir;    }
        format="column";
        
        //calc dists for fastafile
        m->mothurOut("/******************************************/\n");
        m->mothurOut("Running command: dist.seqs(" + options + ")\n");
        Command* distCommand = new DistanceCommand(options);
        
        distCommand->execute();
        map<string, vector<string> > filenames = distCommand->getOutputFiles();
        distfile = filenames[format][0];
        
        delete distCommand;
        current->setMothurCalling(false);
        m->mothurOut("/******************************************/\n");
        
        options = "fitcalc=t, fasta=" + reffastafile + ", oldfasta=" + fastafile + ", cutoff=" + toString(cutoff) + ", column=" + distfile;
        if (outputDir != "")                            { options += ", outputdir=" + outputDir;    }
        
        //dists between reffasta and fastafile
        m->mothurOut("/******************************************/\n");
        m->mothurOut("Running command: dist.seqs(" + options + ")\n");
        distCommand = new DistanceCommand(options);
        
        distCommand->execute();
        filenames = distCommand->getOutputFiles();
        distfile = filenames[format][0];
        
        delete distCommand;
        current->setMothurCalling(false);
        m->mothurOut("/******************************************/\n");

        //create name of count file for combined distance file
        if (countfile != "")     {
            combinedNameOrCount = "count"; combinedNameFile = countfile+".temp";
            CountTable ct; ct.readTable(countfile, false, false);
            
            //add in reffasta seqs
            for (int i = 0; i < list->getNumBins(); i++) {
                vector<string> names;
                string binNames = list->get(i);
                util.splitAtComma(binNames, names);
                for (int j = 0; j < names.size(); j++) { ct.push_back(names[j], 1); }
            }
            
            //print new count table for optiCluster
            ct.printTable(combinedNameFile);
            counts = ct.getNameMap();
        }else if (namefile != "") {
            combinedNameOrCount = "name"; combinedNameFile = namefile+".temp";
            
            string command = "copy ";
#if defined NON_WINDOWS
            command = "cp ";
#endif
            command += namefile + " " + combinedNameFile;
            system(command.c_str());
            
            //add in reffasta seqs
            ofstream out; util.openOutputFileAppend(combinedNameFile, out);
            for (int i = 0; i < list->getNumBins(); i++) {
                vector<string> names;
                string binNames = list->get(i);
                util.splitAtComma(binNames, names);
                for (int j = 0; j < names.size(); j++) { out << names[i] << '\t' << names[i] << endl;  }
            }
            out.close();
        }else {
            //create count file
            CountTable ct; ct.readTable(fastafile, "fasta");
            
            combinedNameOrCount = "count"; combinedNameFile = fastafile+".temp";
            
            //add in reffasta seqs
            for (int i = 0; i < list->getNumBins(); i++) {
                vector<string> names;
                string binNames = list->get(i);
                util.splitAtComma(binNames, names);
                for (int j = 0; j < names.size(); j++) { ct.push_back(names[j], 1); }
            }
            
            //print new count table for optiCluster
            ct.printTable(combinedNameFile);
            counts = ct.getNameMap();
        }
        
        return distfile;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterFitCommand", "calcDists");
        exit(1);
    }
}
//**********************************************************************************************************************

int ClusterFitCommand::runOptiCluster(OptiMatrix& matrix, ListVector* list){
    try {
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
        
        OptiCluster cluster(&matrix, metric, 0);
        tag = cluster.getTag() + ".fit";
        
        m->mothurOut("\nClustering " + distfile + "\n"); 
        
        if (outputDir == "") { outputDir += util.hasPath(distfile); }
        fileroot = outputDir + util.getRootName(util.getSimpleName(distfile));
        
        string listFileName = fileroot+ tag + ".list";
        
        ofstream listFile;
        util.openOutputFile(listFileName,	listFile);
        outputNames.push_back(listFileName); outputTypes["list"].push_back(listFileName);
        
        map<string, string> variables;
        variables["[filename]"] = fileroot;
        variables["[clustertag]"] = tag;
        string outputName = getOutputFileName("steps", variables);
        outputNames.push_back(outputName); outputTypes["steps"].push_back(outputName);
        ofstream outStep;
        util.openOutputFile(outputName, outStep);
        
        int iters = 0;
        double listVectorMetric = 0; //worst state
        double delta = 1;
        
        vector<vector<string> > otus;
        for (int i = 0; i < list->getNumBins(); i++) {
            vector<string> binNames;
            string bin = list->get(i);
            util.splitAtComma(bin, binNames);
            otus.push_back(binNames);
        }
        
        cluster.initialize(listVectorMetric, true, initialize, otus, list->getLabels(), truePositives, trueNegatives, falsePositives, falseNegatives);
        
        long long numBins = cluster.getNumBins();
        m->mothurOut("\n\niter\ttime\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");
        outStep << "iter\ttime\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";
        long long tp, tn, fp, fn;
        vector<double> results = cluster.getStats(tp, tn, fp, fn);
        m->mothurOut("0\t0\t" + toString(cutoff) + "\t" + toString(numBins) + "\t"+ toString(cutoff) + "\t" + toString(tp) + "\t" + toString(tn) + "\t" + toString(fp) + "\t" + toString(fn) + "\t");
        outStep << "0\t0\t" + toString(cutoff) + "\t" + toString(numBins) + "\t" + toString(cutoff) + "\t" << tp << '\t' << tn << '\t' << fp << '\t' << fn << '\t';
        for (int i = 0; i < results.size(); i++) { m->mothurOut(toString(results[i]) + "\t"); outStep << results[i] << "\t"; }
        m->mothurOutEndLine();
        outStep << endl;
        
        while ((delta > stableMetric) && (iters < maxIters)) {
            
            long start = time(NULL);
            
            if (m->getControl_pressed()) { break; }
            double oldMetric = listVectorMetric;
            
            cluster.update(listVectorMetric);
            
            delta = abs(oldMetric - listVectorMetric);
            iters++;
            
            results = cluster.getStats(tp, tn, fp, fn);
            numBins = cluster.getNumBins();
            
            m->mothurOut(toString(iters) + "\t" + toString(time(NULL) - start) + "\t" + toString(cutoff) + "\t" + toString(numBins) + "\t" + toString(cutoff) + "\t"+ toString(tp) + "\t" + toString(tn) + "\t" + toString(fp) + "\t" + toString(fn) + "\t");
            outStep << (toString(iters) + "\t" + toString(time(NULL) - start) + "\t" + toString(cutoff) + "\t" + toString(numBins) + "\t" + toString(cutoff) + "\t") << tp << '\t' << tn << '\t' << fp << '\t' << fn << '\t';
            for (int i = 0; i < results.size(); i++) { m->mothurOut(toString(results[i]) + "\t"); outStep << results[i] << "\t"; }
            m->mothurOutEndLine();
            outStep << endl;
        }
        m->mothurOutEndLine(); m->mothurOutEndLine();
        
        if (m->getControl_pressed()) { delete metric; metric = NULL; return 0; }
        
        set<string> unfittedSeqs;
        ListVector* list = cluster.getList(unfittedSeqs);
        list->setLabel(toString(cutoff));
        if (translate) {
            map<string, string> renameMap;
            //read map file
            ifstream in; util.openInputFile(fastaMapFile, in);
            while (!in.eof()) {
                string name, fixedName;
                in >> name >> fixedName; util.gobble(in);
                renameMap[name] = fixedName;
            }
            in.close();
            
            for (int i = 0; i < list->getNumBins(); i++) {
                string bin = list->get(i);
                if (bin != "") {
                    vector<string> binNames; util.splitAtComma(bin, binNames);
                    string newBin = renameMap[binNames[0]];
                    for (int j = 1; j < binNames.size(); j++) {  newBin += ","+renameMap[binNames[j]];  }
                    list->set(i, newBin);
                }
            }
        }

        if(countfile != "") { list->print(listFile, counts); }
        else { list->print(listFile); }
        listFile.close();
        
        m->mothurOut("\nFitted " + toString(list->getNumSeqs()) + " sequences to existing OTUs. \n");
        delete list;
        
        if (unfittedSeqs.size() != 0) {
            string unfittedListFileName = fileroot+ cluster.getTag() + ".unfitted" + ".list";
            m->mothurOut("Note: Could not fit " + toString(unfittedSeqs.size()) + " sequences. \n");
            ListVector unfittedList; unfittedList.setLabel(toString(cutoff));
            for (set<string>::iterator it = unfittedSeqs.begin(); it != unfittedSeqs.end(); it++) { unfittedList.push_back(*it); }
            ofstream listFile;
            util.openOutputFile(unfittedListFileName,	listFile);
            outputNames.push_back(unfittedListFileName); outputTypes["list"].push_back(unfittedListFileName);
            unfittedList.print(listFile);
            listFile.close();
        }
        
        
        string sensspecFilename = fileroot+ tag + ".sensspec";
        ofstream sensFile;
        util.openOutputFile(sensspecFilename,	sensFile);
        outputNames.push_back(sensspecFilename); outputTypes["sensspec"].push_back(sensspecFilename);
        
        sensFile << "label\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";
        
        results = cluster.getStats(tp, tn, fp, fn);
        
        sensFile << cutoff << '\t' << cutoff << '\t' << tp << '\t' << tn << '\t' << fp << '\t' << fn << '\t';
        for (int i = 0; i < results.size(); i++) {  sensFile << results[i] << '\t'; }
        sensFile << '\n';
        sensFile.close();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterFitCommand", "runOptiCluster");
        exit(1);
    }
    
}
//**********************************************************************************************************************


