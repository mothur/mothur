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
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none","","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount", "none", "","",false,false,true); parameters.push_back(pcount);
        CommandParameter prefname("refname", "InputTypes", "", "", "RefNameCount", "none","","",false,false,true); parameters.push_back(prefname);
        CommandParameter prefcount("refcount", "InputTypes", "", "", "RefNameCount", "none", "","",false,false,true); parameters.push_back(prefcount);
        CommandParameter prefcolumn("refcolumn", "InputTypes", "", "", "PhylipColumnRef", "", "ColumnName","",false,false,true); parameters.push_back(prefcolumn);
        CommandParameter pcolumn("column", "InputTypes", "", "", "PhylipColumn", "", "ColumnName","",false,false,true); parameters.push_back(pcolumn);
        CommandParameter paccnos("accnos", "InputTypes", "", "", "", "", "","",false,false,true); parameters.push_back(paccnos);
        CommandParameter pcutoff("cutoff", "Number", "", "0.03", "", "", "","",false,false,true); parameters.push_back(pcutoff);
        CommandParameter pprecision("precision", "Number", "", "100", "", "", "","",false,false); parameters.push_back(pprecision);
        CommandParameter pmethod("method", "Multiple", "closed-open", "closed", "", "", "","",false,false,true); parameters.push_back(pmethod);
        CommandParameter pcrit("criteria", "Multiple", "fit-combo-both", "both", "", "", "","",false,false,true); parameters.push_back(pcrit);
        CommandParameter prefweight("refweight", "Multiple", "none-abundance-connectivity", "none", "", "", "","",false,false,true); parameters.push_back(prefweight);
        CommandParameter pmetric("metric", "Multiple", "mcc-sens-spec-tptn-fpfn-tp-tn-fp-fn-f1score-accuracy-ppv-npv-fdr", "mcc", "", "", "","",false,false,true); parameters.push_back(pmetric);
        CommandParameter pmetriccutoff("delta", "Number", "", "0.0001", "", "", "","",false,false,true); parameters.push_back(pmetriccutoff);
        CommandParameter piters("iters", "Number", "", "100", "", "", "","",false,false,true); parameters.push_back(piters);
        CommandParameter pdenovoiters("denovoiters", "Number", "", "100", "", "", "","",false,false,true); parameters.push_back(pdenovoiters);
        CommandParameter pfitpercent("fitpercent", "Number", "", "10", "", "", "","",false,false,true); parameters.push_back(pfitpercent);
        CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter prefprint("printref", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(prefprint);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["list"] = tempOutNames;
        outputTypes["sensspec"] = tempOutNames;
        outputTypes["steps"] = tempOutNames;
        
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
        helpString += "The cluster.fit command parameter options are reflist, refcolumn, refname, refcount, fasta, name, count, column, accnos, method, cutoff, precent, metric, iters, initialize, denovoiters.\n";
        helpString += "The refcolumn parameter allow you to enter your reference data distance file, to reduce processing time. \n";
        helpString += "The column parameter allow you to enter your data distance file, to reduce processing time. \n";
        helpString += "The fasta parameter allows you to enter your fasta file. \n";
        helpString += "The reffasta parameter allows you to enter your fasta file for your reference dataset. \n";
        helpString += "The reflist parameter allows you to enter your list file for your reference dataset. \n";
        helpString += "The name parameter allows you to enter your name file. \n";
        helpString += "The count parameter allows you to enter your count file.\nA count or name file is required if your distance file is in column format.\n";
        helpString += "The refname parameter allows you to enter your reference name file. \n";
        helpString += "The refcount parameter allows you to enter your reference count file.\nA refcount or refname file is required if your reference distance file is in column format.\n";
        helpString += "The accnos parameter allows you to assign reference seqeunces by name. This can save time by allowing you to provide a distance matrix containing all the sequence distances rather than a sample matrix and reference matrix and mothur calculating the distances between the sample and reference.\n";
        helpString += "The iters parameter allow you to set the maxiters for the opticluster method. \n";
        helpString += "The denovoiters parameter allow you to set the number of randomizations to perform. \n";
        helpString += "The fitpercent parameter allow you to set percentage of reads to be fitted. Default=50. Max=100, min=0.01.\n";
        helpString += "The criteria parameter allows you to indicate which metric will influence the fitting. Options are fit, combo and both. Default=both. Using fit means a sequence will be fitted to an OTU if the fit makes the metric for the fitted sequences better (only considers metric value generated by fit seqs). Using combo means a sequence will be fitted to an OTU if the fit makes the metric for the fitted and the reference sequences better (considers metric value generated by all reference and fit sequences). Using both means a sequence will be fitted to an OTU if it makes the metric for the fitted sequences better (fit) or the metric for the combo better (combo). \n";
        helpString += "The refweight parameter is used with the denovo method to allows you weight the selection of reference sequences. Options none, abundance and connectivity. Default=none.\n";

        helpString += "The metric parameter allows to select the metric in the opticluster method. Options are Matthews correlation coefficient (mcc), sensitivity (sens), specificity (spec), true positives + true negatives (tptn), false positives + false negatives (fpfn), true positives (tp), true negative (tn), false positive (fp), false negative (fn), f1score (f1score), accuracy (accuracy), positive predictive value (ppv), negative predictive value (npv), false discovery rate (fdr). Default=mcc.\n";
        helpString += "The printref parameter allows to indicate whether you want the reference seqs printed with the fit seqs. For example, if you are trying to see how a new patient's data changes the clustering, you want to set printref=t so the old patient and new patient OTUs are printed together. If you want to see how your data would fit with a reference like silva, setting printref=f would output only your sequences to the list file. By default printref=t for denovo clustering and printref=f when using a reference.\n";
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
//This function checks to make sure the cluster command has no errors and then clusters based on the method chosen.
ClusterFitCommand::ClusterFitCommand(string option)  {
    try{
        //allow user to run help
        if(option == "help") { help(); abort = true; calledHelp = true; }
        else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
        
        else {
            OptionParser parser(option, setParameters());
            map<string,string> parameters = parser.getParameters();
            
            ValidParameters validParameter;
            selfReference = true; createAccnos = false; refdistfile = ""; distfile = "";
            
            //check for required parameters
            reffastafile = validParameter.validFile(parameters, "reffasta");
            if (reffastafile == "not open") { abort = true; }
            else if (reffastafile == "not found") { reffastafile = "";  }
            else {  selfReference = false; }
            
            refdistfile = validParameter.validFile(parameters, "refcolumn");
            if (refdistfile == "not open") { refdistfile = ""; abort = true; }
            else if (refdistfile == "not found") { refdistfile = ""; }
            else {  refformat = "column"; selfReference = false; }
            
            //allow ref list to be entered with denovo and accnos file
            reflistfile = validParameter.validFile(parameters, "reflist");
            if (reflistfile == "not open") { abort = true; }
            else if (reflistfile == "not found") { reflistfile = ""; }
            //else { selfReference = false; }
            
            refnamefile = validParameter.validFile(parameters, "refname");
            if (refnamefile == "not open") { abort = true; }
            else if (refnamefile == "not found") { refnamefile = ""; }
            else { selfReference = false; }
            
            refcountfile = validParameter.validFile(parameters, "refcount");
            if (refcountfile == "not open") { abort = true;  }
            else if (refcountfile == "not found") { refcountfile = ""; }
            else { selfReference = false; }
            
            
            if (!selfReference) { //if you are providing reference files, lets make sure we have all of them
                if ((refdistfile == "") || (reffastafile == "") || (reflistfile == "")) { m->mothurOut("[ERROR]: When providing a reference file, you must provide a reffasta, refcolumn, reflist and refcount or refname, aborting.\n");  abort = true; }
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
            
            accnosfile = validParameter.validFile(parameters, "accnos");
            if (accnosfile == "not open") { accnosfile = ""; abort = true; }
            else if (accnosfile == "not found") { accnosfile = ""; }
            else {   current->setAccnosFile(accnosfile); createAccnos = false;	}
            
            //extract reference names from reflist instead of accnos fil
            if (selfReference) {
                if ((reflistfile != "") && (accnosfile == "")) { createAccnos = true; }
            }
            
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
            else { m->mothurOut("[ERROR]: Not a valid metric.  Valid metrics are mcc, sens, spec, tp, tn, fp, fn, tptn, fpfn, f1score, accuracy, ppv, npv, fdr.\n");  abort = true; }
            
            criteria = validParameter.valid(parameters, "criteria");		if (criteria == "not found") { criteria = "both"; }
            if ((criteria == "fit") || (criteria == "combo") || (criteria == "both")){ }
            else { m->mothurOut("[ERROR]: Not a valid criteria.  Valid criteria are fit, combo and both.\n"); abort = true; }
            
            refWeight = validParameter.valid(parameters, "refweight");		if (refWeight == "not found") { refWeight = "none"; }
            if ((refWeight == "none") || (refWeight == "abundance") || (refWeight == "connectivity")){ }
            else { m->mothurOut("[ERROR]: Not a valid reference weight.  Valid refweight options are none, abundance and connectivity.\n"); abort = true; }
            
            initialize = "singleton";
            
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
            
            temp = validParameter.valid(parameters, "processors");    if (temp == "not found"){    temp = current->getProcessors();    }
            processors = current->setProcessors(temp);
            
            adjust=-1.0;
            temp = validParameter.valid(parameters, "cutoff");
            if (temp == "not found") { temp = "0.03"; }
            util.mothurConvert(temp, cutoff);
            
            temp = validParameter.valid(parameters, "printref");			if (temp == "not found") { if (selfReference) { temp = "t"; }else { temp = "f"; } }
            printref = util.isTrue(temp);
            
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
        
        if (outputdir == "") { outputdir += util.hasPath(distfile); }
        fileroot = outputdir + util.getRootName(util.getSimpleName(distfile));
        
        map<string, string> variables;
        variables["[filename]"] = fileroot;
        variables["[clustertag]"] = "optifit_" + metric->getName();
        string outputName = getOutputFileName("steps", variables);
        string listFile = "";
        string bestListFileName = "";
        
        if (selfReference) { //de novo
            
            if ((accnosfile == "") && (!createAccnos)) { //denovo with mothur randomly assigning references
                
                m->mothurOut("\nRandomly assigning reads from " + distfile + " as reference sequences\n");
                
                //distfile, distFormat, dupsFile, dupsFormat, cutoff, percentage to be fitseqs - will randomly assign as fit
                OptiData* matrix = new OptiRefMatrix(distfile, "column", dupsFile, nameOrCount, cutoff, fitPercent, refWeight);
                
                runDenovoOptiCluster(matrix, metric, counts, outputName);
                
                string sensspecFilename = fileroot+ tag + ".sensspec";
                ofstream sensFile;
                util.openOutputFile(sensspecFilename,    sensFile);
                outputNames.push_back(sensspecFilename); outputTypes["sensspec"].push_back(sensspecFilename);
         
                //evaluate results
                bestListFileName = compareSensSpec(matrix, metric, sensFile);
                
                delete matrix;
                
            }else { //reference with accnos file or reference list file assigning references
                
                set<string> refNames; vector<string> refLabels; vector< vector<string> > otus;
                
                if (accnosfile != "") { //use accnos file to assign references
                    
                    m->mothurOut("\nUsing sequences from " + accnosfile + " as reference sequences\n");
                    
                    refNames = util.readAccnos(accnosfile);
                    
                }else if (createAccnos) { //assign references based on reflist parameter
                    
                    m->mothurOut("\nUsing OTUs from " + reflistfile + " as reference OTUs\n");
                    
                    InputData input(reflistfile, "list", nullVector);
                    set<string> processedLabels, userLabels;
                    string lastLabel = "";
                    
                    ListVector* reflist = util.getNextList(input, true, userLabels, processedLabels, lastLabel);
                    
                    refLabels = reflist->getLabels();
                    for (int i = 0; i < refLabels.size(); i++) { refLabels[i] = "Ref_" + refLabels[i];  }
                    
                    refNames = util.getSetFromList(reflist, otus); delete reflist;
                }
                
                //distfile, distFormat, dupsFile, dupsFormat, cutoff, accnos containing refseq name
                OptiData* matrix = new OptiRefMatrix(distfile, "column", dupsFile, nameOrCount, cutoff, refNames);
                
                //fit seqs
                ListVector* list = runUserRefOptiCluster(matrix, metric, counts, outputName, refLabels, otus);
                
                ofstream listFile; string listFileName = fileroot+ tag + ".list";
                util.openOutputFile(listFileName,    listFile);
                
                if(countfile != "") { list->print(listFile, counts); }
                else { list->print(listFile); }
                listFile.close();
                
                listFiles.push_back(listFileName);
                bestListFileName = listFileName;
                
                delete list; delete matrix;
            }
        }else { //reference with files containing reference seqs
            
            createReferenceNameCount(); //creates reference name or count file if needed
            
            calcDists();  //calc distance matrix for fasta file and distances between fasta file and reffasta file
            
            m->mothurOut("\nUsing OTUs from " + reflistfile + " as reference OTUs\n");
            
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
            
            bestListFileName = listFile;
           
            delete matrix;
        }
        delete metric;
        
        if (m->getControl_pressed()) { 	for (int j = 0; j < outputNames.size(); j++) { util.mothurRemove(outputNames[j]); }  return 0; }
        
        outputNames.push_back(outputName); outputTypes["steps"].push_back(outputName);
        outputNames.push_back(bestListFileName); outputTypes["list"].push_back(bestListFileName);
        
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
string ClusterFitCommand::runDenovoOptiCluster(OptiData*& matrix, ClusterMetric*& metric, map<string, int>& counts, string outStepFile){
    try {
        m->mothurOut("\nClustering " + distfile + "\n");
        bool printStepsHeader = true;
        
        for (int i = 0; i < denovoIters; i++) {
            
            OptiFitCluster cluster(matrix, metric, 0, criteria);
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
            double tp, tn, fp, fn;
            vector<double> results = cluster.getStats(tp, tn, fp, fn);
            
            double fittp, fittn, fitfp, fitfn;
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
            
            ListVector* list = cluster.getFittedList(toString(cutoff), printref);
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
//**********************************************************************************************************************
ListVector* ClusterFitCommand::runUserRefOptiCluster(OptiData*& matrix, ClusterMetric*& metric, map<string, int>& counts, string outStepFile, vector<string> refListLabels, vector<vector<string> > otus){
    try {
        
        bool printStepsHeader = true;

        OptiFitCluster cluster(matrix, metric, 0, criteria);
        tag = cluster.getTag();
        
        int iters = 0;
        double listVectorMetric = 0; //worst state
        double delta = 1;
        
        if (!createAccnos) {
            
            m->mothurOut("\nClustering references from " + distfile + "\n");
            
            //get "ref" seqs for initialize inputs
            OptiData* refMatrix = matrix->extractRefMatrix();
            ListVector* refList = clusterRefs(refMatrix, metric); delete refMatrix;
            
            for (int i = 0; i < refList->getNumBins(); i++) {
                vector<string> binNames;
                string bin = refList->get(i);
                if (bin != "") {
                    util.splitAtComma(bin, binNames);
                    otus.push_back(binNames);
                }
            }
            
            //add tag to OTULabels to indicate the reference
            refListLabels = refList->getLabels();
            for (int i = 0; i < refListLabels.size(); i++) { refListLabels[i] = "Ref_" + refListLabels[i];  }
            refList->setLabels(refListLabels);
            delete refList;
        }
        
        cluster.initialize(listVectorMetric, true, otus, refListLabels, method, false);
        
        long long numBins = cluster.getNumBins();
        double tp, tn, fp, fn;
        vector<double> results = cluster.getStats(tp, tn, fp, fn);
        
        double fittp, fittn, fitfp, fitfn;
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
            
            outputSteps(outStepFile, printStepsHeader, tp, tn, fp, fn, results, numBins, fittp, fittn, fitfp, fitfn, fitresults, numFitBins, iters, false, 0);
        }
        m->mothurOutEndLine(); m->mothurOutEndLine();
        
        if (m->getControl_pressed()) {  return 0; }
        
        ListVector* list = cluster.getFittedList(toString(cutoff), printref);
        list->setLabel(toString(cutoff));
        
        string sensspecFilename = fileroot+ tag + ".sensspec";
        ofstream sensFile;
        util.openOutputFile(sensspecFilename,    sensFile);
        outputNames.push_back(sensspecFilename); outputTypes["sensspec"].push_back(sensspecFilename);
        
        if (method == "closed") {
            sensFile << "label\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";
            
            if (printref) { //combo
                results = cluster.getStats(tp, tn, fp, fn);
                sensFile << cutoff << '\t' << cutoff << '\t' << tp << '\t' << tn << '\t' << fp << '\t' << fn << '\t';
                for (int i = 0; i < results.size(); i++) {  sensFile << results[i] << '\t'; } sensFile << '\n';
            }else { //fit
                fitresults = cluster.getFitStats(fittp, fittn, fitfp, fitfn);
                sensFile << cutoff << '\t' << cutoff << '\t' << fittp << '\t' << fittn << '\t' << fitfp << '\t' << fitfn << '\t';
                for (int i = 0; i < fitresults.size(); i++) {  sensFile << fitresults[i] << "\t"; } sensFile << endl;
            }
        }else {
            runSensSpec(matrix, metric, list, counts, sensFile);
        }
        sensFile.close();
        
        return list;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterFitCommand", "runUserRefOptiCluster");
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
        
        double tp, tn, fp, fn;
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
        OptiFitCluster cluster(matrix, metric, 0, criteria);
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
        
        map<string, int> refCounts;
        if (refcountfile != "") {
            CountTable refct; refct.readTable(refcountfile, false, false);
            refCounts = refct.getNameMap();
        }else if (refnamefile != "") { refCounts = util.readNames(refnamefile); }
        else { //assume unique
            for (int i = 0; i < otus.size(); i++) { for (int j = 0; j < otus[i].size(); j++) { refCounts[otus[i][j]] = 1; } }
        }
        counts.insert(refCounts.begin(), refCounts.end());
        
        cluster.initialize(listVectorMetric, true, otus, refList->getLabels(), method, false);
        
        long long numBins = cluster.getNumBins();
        double tp, tn, fp, fn;
        vector<double> results = cluster.getStats(tp, tn, fp, fn);
        
        double fittp, fittn, fitfp, fitfn;
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
        
        ListVector* list = cluster.getFittedList(toString(cutoff), printref);
        list->setLabel(toString(cutoff));
        
        ofstream listFile;
        string listFileName = fileroot+ tag + ".list";
        util.openOutputFile(listFileName,	listFile);
        
        if(countfile != "") { list->print(listFile, counts); }
        else { list->print(listFile); }
        listFile.close();
        
        string sensspecFilename = fileroot+ tag + ".sensspec";
        ofstream sensFile;
        util.openOutputFile(sensspecFilename,    sensFile);
        outputNames.push_back(sensspecFilename); outputTypes["sensspec"].push_back(sensspecFilename);
        
        if (method == "closed") {
            sensFile << "label\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";
         
            if (printref) { //combo
                results = cluster.getStats(tp, tn, fp, fn);
                sensFile << cutoff << '\t' << cutoff << '\t' << tp << '\t' << tn << '\t' << fp << '\t' << fn << '\t';
                for (int i = 0; i < results.size(); i++) {  sensFile << results[i] << '\t'; } sensFile << '\n';
            }else { //fit
                fitresults = cluster.getFitStats(fittp, fittn, fitfp, fitfn);
                sensFile << cutoff << '\t' << cutoff << '\t' << fittp << '\t' << fittn << '\t' << fitfp << '\t' << fitfn << '\t';
                for (int i = 0; i < fitresults.size(); i++) {  sensFile << fitresults[i] << "\t"; } sensFile << endl;
            }
        }else {
            runSensSpec(matrix, metric, list, counts, sensFile);
        }
        sensFile.close();
        
        delete list;

        return listFileName;
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterFitCommand", "runRefOptiCluster");
        exit(1);
    }
    
}
//**********************************************************************************************************************
string ClusterFitCommand::compareSensSpec(OptiData*& matrix, ClusterMetric*& userMetric, ofstream& sensSpecFile) {
    try {
           
        sensSpecFile << "iter\tlabel\tcutoff\tnumotus\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";
        m->mothurOut("iter\tlabel\tcutoff\tnumotus\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");
        
        double bestStat = 0; int bestResult = 0;
        
        if ((method == "open") && (printref)) {
            
            for (int i = 0; i < listFiles.size(); i++) {
                string thislistFileName = listFiles[i];
                
                InputData input(thislistFileName, "list", nullVector);
                ListVector* list = input.getListVector();
                
                string label = list->getLabel();
                int numBins = list->getNumBins();
                
                SensSpecCalc senscalc(*matrix, list);
                double truePositives, trueNegatives, falsePositives, falseNegatives;
                senscalc.getResults(*matrix, truePositives, trueNegatives, falsePositives, falseNegatives);
                
                double tp =  truePositives;
                double fp =  falsePositives;
                double tn =  trueNegatives;
                double fn =  falseNegatives;
                
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
        }else {
            
            for (int i = 0; i < listFiles.size(); i++) {
               
                InputData input(listFiles[i], "list", nullVector);
                ListVector* list = input.getListVector();
                
                string label = list->getLabel();
                int numBins = list->getNumBins();
                
                //extract seqs in list file from matrix
                set<string> listNames;
                for (int i = 0; i < list->getNumBins(); i++){
                    string bin = list->get(i);
                    if (bin != "") {
                        vector<string> binSeqs; util.splitAtComma(bin, binSeqs);
                        for (int j = 0; j < binSeqs.size(); j++) {
                            listNames.insert(binSeqs[j]);
                        }
                    }
                }
                
                OptiData* fitMatrix = matrix->extractMatrixSubset(listNames);
                SensSpecCalc senscalc(*fitMatrix, list);
                double truePositives, trueNegatives, falsePositives, falseNegatives;
                senscalc.getResults(*fitMatrix, truePositives, trueNegatives, falsePositives, falseNegatives);
                delete fitMatrix;
                
                double tp =  truePositives;
                double fp =  falsePositives;
                double tn =  trueNegatives;
                double fn =  falseNegatives;
                
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
        }
        
        sensSpecFile.close();
        
        return listFiles[bestResult];
        
    }
    catch(exception& e) {
        m->errorOut(e, "ClusterFitCommand", "compareSensSpec");
        exit(1);
    }
}
//**********************************************************************************************************************
void ClusterFitCommand::runSensSpec(OptiData*& matrix, ClusterMetric*& userMetric, ListVector*& list, map<string, int>& counts, ofstream& sensSpecFile) {
    try {
        
        sensSpecFile << "label\tcutoff\tnumotus\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";
        m->mothurOut("label\tcutoff\tnumotus\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");
        
        if (method == "open") {
            double truePositives, trueNegatives, falsePositives, falseNegatives;
            string label = list->getLabel();
            int numBins = list->getNumBins();
            
            if (printref) { //pass whole matrix
                SensSpecCalc senscalc(*matrix, list);
                senscalc.getResults(*matrix, truePositives, trueNegatives, falsePositives, falseNegatives);
            }else { //pass subset matrix
                vector<long long> fSeqs = matrix->getFitSeqs();
                set<long long> fitSeqs = util.mothurConvert(fSeqs);
                OptiData* fitMatrix = matrix->extractMatrixSubset(fitSeqs);
                SensSpecCalc senscalc(*fitMatrix, list);
                senscalc.getResults(*fitMatrix, truePositives, trueNegatives, falsePositives, falseNegatives);
                delete fitMatrix;
            }
            
            double tp =  truePositives;
            double fp =  falsePositives;
            double tn =  trueNegatives;
            double fn =  falseNegatives;
            
            Sensitivity sens;   double sensitivity = sens.getValue(tp, tn, fp, fn);
            Specificity spec;   double specificity = spec.getValue(tp, tn, fp, fn);
            PPV ppv;            double positivePredictiveValue = ppv.getValue(tp, tn, fp, fn);
            NPV npv;            double negativePredictiveValue = npv.getValue(tp, tn, fp, fn);
            FDR fdr;            double falseDiscoveryRate = fdr.getValue(tp, tn, fp, fn);
            Accuracy acc;       double accuracy = acc.getValue(tp, tn, fp, fn);
            MCC mcc;            double matthewsCorrCoef = mcc.getValue(tp, tn, fp, fn);
            F1Score f1;         double f1Score = f1.getValue(tp, tn, fp, fn);
            
            sensSpecFile << label << '\t' << cutoff << '\t' << numBins << '\t';
            sensSpecFile << truePositives << '\t' << trueNegatives << '\t' << falsePositives << '\t' << falseNegatives << '\t';
            sensSpecFile << setprecision(4);
            sensSpecFile << sensitivity << '\t' << specificity << '\t' << positivePredictiveValue << '\t' << negativePredictiveValue << '\t';
            sensSpecFile << falseDiscoveryRate << '\t' << accuracy << '\t' << matthewsCorrCoef << '\t' << f1Score << endl;
            
            m->mothurOut(label + "\t" + toString(cutoff) + "\t" + toString(numBins) + "\t"+ toString(truePositives) + "\t" + toString(trueNegatives) + "\t" + toString(falsePositives) + "\t" + toString(falseNegatives) + "\t");
            m->mothurOut(toString(sensitivity) + "\t" + toString(specificity) + "\t" + toString(positivePredictiveValue) + "\t" + toString(negativePredictiveValue) + "\t");
            m->mothurOut(toString(falseDiscoveryRate) + "\t" + toString(accuracy) + "\t" + toString(matthewsCorrCoef) + "\t" + toString(f1Score) + "\n\n");
        }else {
            m->mothurOut("[ERROR]: should never get here... \n");
        }

    }
    catch(exception& e) {
        m->errorOut(e, "ClusterFitCommand", "runSensSpec");
        exit(1);
    }
}
//**********************************************************************************************************************
void ClusterFitCommand::outputSteps(string outputName, bool& printHeaders, double tp, double tn, double fp, double fn, vector<double> results, long long numBins, double fittp, double fittn, double fitfp, double fitfn, vector<double> fitresults, long long numFitBins, int iter, bool printToFile, int denovoIter) {
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


