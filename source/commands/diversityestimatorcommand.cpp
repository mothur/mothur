//
//  diversityestimatorcommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/4/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "diversityestimatorcommand.hpp"
#include "erarefaction.hpp"
#include "metroig.hpp"
#include "metrolognormal.hpp"
#include "metrologstudent.hpp"
#include "metrosichel.hpp"
#include "igabundance.hpp"
#include "igrarefaction.hpp"
#include "lnabundance.hpp"
#include "lnrarefaction.hpp"
#include "lnshift.hpp"
#include "lsabundance.hpp"
#include "lsrarefaction.hpp"
#include "siabundance.hpp"
#include "sirarefaction.hpp"
#include "sishift.hpp"

//**********************************************************************************************************************
vector<string> EstimatorSingleCommand::setParameters(){
    try {
        CommandParameter plist("list", "InputTypes", "", "", "LRSS", "LRSS", "none","",false,false,true); parameters.push_back(plist);
        CommandParameter prabund("rabund", "InputTypes", "", "", "LRSS", "LRSS", "none","",false,false,true); parameters.push_back(prabund);
        CommandParameter psabund("sabund", "InputTypes", "", "", "LRSS", "LRSS", "none","",false,false,true); parameters.push_back(psabund);
        CommandParameter psample("sample", "InputTypes", "", "", "LRSS", "LRSS", "none","",false,false,true); parameters.push_back(psample);
        CommandParameter pshared("shared", "InputTypes", "", "", "LRSS", "LRSS", "none","",false,false,true); parameters.push_back(pshared);
        CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
        CommandParameter pfreq("freq", "Number", "", "100", "", "", "","",false,false); parameters.push_back(pfreq);
        CommandParameter pcalc("calc", "Multiple", "erarefact-ig-ln-ls-si-igabund-igrarefact-lnrarefact-lnabund-lnshift-lsabund-lsrarefact-siabund-sirarefact-sishift", "ig", "", "", "","",false,false,true); parameters.push_back(pcalc); //lnabund
        CommandParameter palpha("sigmaa", "Number", "", "0.1", "", "", "","",false,false,true); parameters.push_back(palpha);
        CommandParameter pbeta("sigmab", "Number", "", "0.1", "", "", "","",false,false); parameters.push_back(pbeta);
        CommandParameter psigman("sigman", "Number", "", "0.1", "", "", "","",false,false); parameters.push_back(psigman);
        CommandParameter psigmas("sigmas", "Number", "", "100", "", "", "","",false,false); parameters.push_back(psigmas);
        CommandParameter pburn("burn", "Number", "", "2000000", "", "", "","",false,false); parameters.push_back(pburn);
        CommandParameter pcoverage("coverage", "Number", "", "0.8", "", "", "","",false,false); parameters.push_back(pcoverage);
        CommandParameter pfit("fit", "Number", "", "20", "", "", "","",false,false); parameters.push_back(pfit);
        CommandParameter psamplenum("burnsample", "Number", "", "1000", "", "", "","",false,false); parameters.push_back(psamplenum);
        CommandParameter piters("iters", "Number", "", "1000", "", "", "","",false,false); parameters.push_back(piters);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["erarefact"] = tempOutNames;
        outputTypes["igrarefact"] = tempOutNames;
        outputTypes["igabund"] = tempOutNames;
        outputTypes["lnabund"] = tempOutNames;
        outputTypes["lnrarefact"] = tempOutNames;
        outputTypes["lnshift"] = tempOutNames;
        outputTypes["lsabund"] = tempOutNames;
        outputTypes["lsrarefact"] = tempOutNames;
        outputTypes["siabund"] = tempOutNames;
        outputTypes["sirarefact"] = tempOutNames;
        outputTypes["sishift"] = tempOutNames;
        outputTypes["ig"] = tempOutNames;
        outputTypes["ln"] = tempOutNames;
        outputTypes["ls"] = tempOutNames;
        outputTypes["si"] = tempOutNames;
        outputTypes["sample"] = tempOutNames;
        
        vector<string> myArray;
        for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
        return myArray;
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "setParameters");
        exit(1);
    }
}
//**********************************************************************************************************************
string EstimatorSingleCommand::getHelpString(){
    try {
        string helpString = "\n";
        ValidCalculators validCalculator;
        helpString += "The estimator.single command parameters are " + getCommandParameters() + ". You may only choose one calculator at a time.\n";
        helpString += "The estimator.single command should be in the following format: \n";
        helpString += "estimator.single(list=yourListFile, calc=yourEstimators).\n";
        helpString += "Example estimator.single(list=final.opti_mcc.list, calc=erarefaction).\n";
        helpString += "The freq parameter is used indicate when to output your data, by default it is set to 100. But you can set it to a percentage of the number of sequence. For example freq=0.10, means 10%. \n";
        helpString += "The sample file is used to provide mcmc sampling to the calculators.\n";
        helpString += "The default values for freq is 100, and calc is erarefaction.\n";
        helpString += "The sigmaa parameter is used to set the std. dev. of alpha / X / mean prop. distn for MetroIG / MetroLogNormal / MetroLogStudent / MetroSichel, respectively. Default = 0.10. n";
        helpString += "The sigmab parameter is used to set the std. dev. of beta / Y / V prop. distn for MetroIG / MetroLogNormal / MetroLogStudent / MetroSichel, respectively. Default = 0.10. n";
        helpString += "The sigman parameter is used to set the std. dev. of N / Gamma prop. distn for MetroLogStudent / MetroSichel, respectively. Default = 0.10.\n";
        helpString += "The sigmas parameter is used to set the std. dev. of S prop. distn for MetroIG / MetroLogNormal / MetroLogStudent / MetroSichel. Default = 100. n";
        helpString += "The coverage parameter allows you to the desired coverage.  Default=0.8.\n";
        helpString += "The iters parameter allows you to set number of mcmc samples to generate.  The default is 250000.\n";
        helpString += "The burn parameter allows ignore part of the sampling file.  Default = 200000 / 100000 for IGAbundance, LNShift, LSAbundance / IGRarefaction, LNRarefaction, LSRarefaction, SIAbundance, SIRarefaction, SIShift respectively.\n";
        helpString += "The burnsample parameter allows you to set sampling frequency.  The default is 1000 / 100 for IGAbundance, LNShift, LSAbundance / IGRarefaction, LNRarefaction, LSRarefaction, SIAbundance, SIRarefaction, SIShift respectively.\n";
        helpString += "The fit parameter is used to indicate to mothur you want mothur to auto adjust the sampling data parameters. default=10, meaning try fitting 10 times. \n";
        helpString += validCalculator.printCalc("estimator");
        helpString += "Be sure to use the correct sampling estimator with your calculator. IG is used for igabund and igrarefact. LN is used for lnabund, lnshift and lnrarefact. LS is used for lsabund and lsrarefaction. SI is used for siabund, sirarefact and sishift.\n";
        helpString += "The label parameter is used to analyze specific labels in your input.\n";
        
        getCommonQuestions();
        
        return helpString;
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "getHelpString");
        exit(1);
    }
}
//**********************************************************************************************************************
string EstimatorSingleCommand::getCommonQuestions(){
    try {
        vector<string> questions, issues, qanswers, ianswers, howtos, hanswers;
        
        string howto = "How do you create the sampling files?"; howtos.push_back(howto);
        string hanswer = "\tRun a short trial MCMC run of 1000 iterations with guessed std. dev.s for the proposal distributions say about 10% of the parameter values. Adjust the std. dev.s until the acceptance ratios are about 0.5. Then perform a longer run of say 250,000 iterations (mothur's default). Three data files with posterior samples for three different sets of parameter values will be generated.\n"; hanswers.push_back(hanswer);
        
        string commonQuestions = util.getFormattedHelp(questions, qanswers, issues, ianswers, howtos, hanswers);
        
        return commonQuestions;
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "getCommonQuestions");
        exit(1);
    }
}

//**********************************************************************************************************************
string EstimatorSingleCommand::getOutputPattern(string type) {
    try {
        string pattern = "[filename],[distance]," + type;
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
EstimatorSingleCommand::EstimatorSingleCommand(string option) : Command()  {
    try {
        allLines = true;
        
        //allow user to run help
        if(option == "help") { help(); calledHelp = true; abort = true; }
        else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
        
        else {
            OptionParser parser(option, setParameters());
            map<string,string> parameters = parser.getParameters();
             
            ValidParameters validParameter;
            listfile = validParameter.validFile(parameters, "list");
            if (listfile == "not open") { listfile = ""; abort = true; }
            else if (listfile == "not found") { listfile = ""; }
            else {  format = "list"; inputfile = listfile; current->setListFile(listfile); }
            
            sabundfile = validParameter.validFile(parameters, "sabund");
            if (sabundfile == "not open") { sabundfile = ""; abort = true; }
            else if (sabundfile == "not found") { sabundfile = ""; }
            else {  format = "sabund"; inputfile = sabundfile; current->setSabundFile(sabundfile); }
            
            rabundfile = validParameter.validFile(parameters, "rabund");
            if (rabundfile == "not open") { rabundfile = ""; abort = true; }
            else if (rabundfile == "not found") { rabundfile = ""; }
            else {  format = "rabund"; inputfile = rabundfile; current->setRabundFile(rabundfile); }
            
            sharedfile = validParameter.validFile(parameters, "shared");
            if (sharedfile == "not open") { sharedfile = ""; abort = true; }
            else if (sharedfile == "not found") { sharedfile = ""; }
            else {  format = "sharedfile"; inputfile = sharedfile; current->setSharedFile(sharedfile); }
            
            bool hasSample = false;
            samplefile = validParameter.validFile(parameters, "sample");
            if (samplefile == "not open") { samplefile = ""; abort = true; }
            else if (samplefile == "not found") { samplefile = ""; }
            else { hasSample = true; current->setSampleFile(samplefile);  }
            
            
            
            
            
            if ((sharedfile == "") && (listfile == "") && (rabundfile == "") && (sabundfile == "")) {
                //is there are current file available for any of these?
                //give priority to shared, then list, then rabund, then sabund
                //if there is a current shared file, use it
                sharedfile = current->getSharedFile();
                if (sharedfile != "") { inputfile = sharedfile; format = "sharedfile"; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n"); }
                else {
                    listfile = current->getListFile();
                    if (listfile != "") { inputfile = listfile; format = "list"; m->mothurOut("Using " + listfile + " as input file for the list parameter.\n");  }
                    else {
                        rabundfile = current->getRabundFile();
                        if (rabundfile != "") { inputfile = rabundfile; format = "rabund"; m->mothurOut("Using " + rabundfile + " as input file for the rabund parameter.\n");  }
                        else {
                            sabundfile = current->getSabundFile();
                            if (sabundfile != "") { inputfile = sabundfile; format = "sabund"; m->mothurOut("Using " + sabundfile + " as input file for the sabund parameter.\n");  }
                            else {
                                m->mothurOut("[ERROR]: No valid current files. You must provide a list, sabund, rabund or shared file before you can use the estimator.single command.\n");
                                abort = true;
                            }
                        }
                    }
                }
            }
            
            //check for optional parameter and set defaults
            // ...at some point should added some additional type checking...
            label = validParameter.valid(parameters, "label");
            if (label == "not found") { label = ""; }
            else {
                if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
                else { allLines = true;  }
            }
            
            //NOTE: if you add new calc options, don't forget to add them to the parameter initialize in setParameters or the gui won't be able to use them
            ValidCalculators validCalculator;
            calc = validParameter.valid(parameters, "calc"); if (calc == "not found") { calc = "ig";  }
            
            samplingCalcs.insert("ig");
            samplingCalcs.insert("ln");
            samplingCalcs.insert("ls");
            samplingCalcs.insert("si");
            
            rarefactCalcs.push_back("igrarefact"); calcToSamplingCalc["igrarefact"] = "ig";
            rarefactCalcs.push_back("lsrarefact"); calcToSamplingCalc["lsrarefact"] = "ls";
            rarefactCalcs.push_back("lnrarefact"); calcToSamplingCalc["lnrarefact"] = "ln";
            rarefactCalcs.push_back("sirarefact"); calcToSamplingCalc["sirarefact"] = "si";
            
            abundCalcs.push_back("igabund"); calcToSamplingCalc["igabund"] = "ig";
            abundCalcs.push_back("lnabund"); calcToSamplingCalc["lnabund"] = "ln";
            abundCalcs.push_back("lsabund"); calcToSamplingCalc["lsabund"] = "ls";
            abundCalcs.push_back("siabund"); calcToSamplingCalc["siabund"] = "si";
            abundCalcs.push_back("sishift"); calcToSamplingCalc["sishift"] = "si";
            abundCalcs.push_back("lnshift"); calcToSamplingCalc["lnshift"] = "ln";
            abundCalcs.push_back("erarefact");
            
            smallBurn.push_back("erarefact");
            smallBurn.push_back("siabund");
            smallBurn.push_back("sishift");
            smallBurn.insert(smallBurn.end(), rarefactCalcs.begin(), rarefactCalcs.end());

            //remove any typo calcs
            createSampling = false;
            if (validCalculator.isValidCalculator("estimator", calc) ) {
                
                bool ignore = false;
                if (!hasSample) { //if you didn't provide a mcmc sample file, but are trying to run a calc that needs it, then ignore
                    if (samplingCalcs.count(calc) == 0) { ignore = true; }
                    if (calc == "erarefact") { ignore = false; }
                }
                
                if (ignore) { m->mothurOut("\n[WARNING]: " + calc + " requires a mcmc sampling file and you have not provided one. You can produce a sampling file using the ig (metroig), ln (metroln), ls (metrols) or si (metrosichel) calculators. I will create the sampling file for you using the " + calcToSamplingCalc[calc] + " calculator.\n"); createSampling = true; }
            }
           
            if (calc == "") { abort = true; m->mothurOut("[ERROR]: no valid estimators, aborting.\n"); }

            string temp;
            temp = validParameter.valid(parameters, "freq");			if (temp == "not found") { temp = "100"; }
            util.mothurConvert(temp, freq);
            
            temp = validParameter.valid(parameters, "sigmaa");		if (temp == "not found") { temp = "0.1"; }
            util.mothurConvert(temp, sigmaAlpha);
            
            temp = validParameter.valid(parameters, "sigmab");		if (temp == "not found") { temp = "0.1"; }
            util.mothurConvert(temp, sigmaBeta);
            
            temp = validParameter.valid(parameters, "sigman");		if (temp == "not found") { temp = "0.1"; }
            util.mothurConvert(temp, sigmaN);
            
            temp = validParameter.valid(parameters, "sigmas");		if (temp == "not found") { temp = "100.0"; }
            util.mothurConvert(temp, sigmaS);
            
            itersSet = true;
            temp = validParameter.valid(parameters, "iters");		if (temp == "not found") { temp = "250000"; itersSet = false; }
            util.mothurConvert(temp, iters);
            
            temp = validParameter.valid(parameters, "fit");		if (temp == "not found") { temp = "10"; }
            util.mothurConvert(temp, fitIters);
            
            temp = validParameter.valid(parameters, "burn");
            if (temp == "not found") {
                if (util.inUsersGroups(calc, smallBurn)) { temp = "100000"; }
                else {  temp = "2000000";  }
            }
            util.mothurConvert(temp, burn);
            
            temp = validParameter.valid(parameters, "burnsample");
            if (temp == "not found") {
                if (util.inUsersGroups(calc, smallBurn)) { temp = "100"; }
                else {  temp = "1000";  }
            }
            util.mothurConvert(temp, burnSample);
            if (burnSample <= 0)  {  m->mothurOut("[ERROR]: Burn sample must be greater than 0. Aborting.\n"); abort=true; }
            
            temp = validParameter.valid(parameters, "coverage");		if (temp == "not found") { temp = "0.8"; }
            util.mothurConvert(temp, coverage);
            
            if ((util.isEqual(coverage, -1)) && ((calc == "igrarefact") || (calc == "lnrarefact") || (calc == "lsrarefact") || (calc == "sirarefact"))) {
                m->mothurOut("[ERROR]: You must set the coverage parameter to run the igrarefact, lsrarefact, lnrarefact or sirarefact estimator. Aborting.\n"); abort=true;
            }
            
            #ifdef USE_GSL
            #else
            
            m->mothurOut("[ERROR]: You did not build mothur with the GNU Scientific Library which is required before you can use the estimator.single command. Aborting.\n");
            abort = true;
            
            
            #endif
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "EstimatorSingleCommand");
        exit(1);
    }
}
//**********************************************************************************************************************

int EstimatorSingleCommand::execute(){
    try {
        
        if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        if (createSampling) {
            string savedCalc = calc;
            int savedIters = iters;
            
            calc = calcToSamplingCalc[savedCalc];
            if (!itersSet) { iters = 250000; }
            
            if (format != "sharedfile") { processSingleSample(); } //handles multiple label values
            else { processSharedFile(); } //handles multiple label values and multiple samples
            
            vector<string> samplingFiles = outputTypes[calc];
            
            if (samplingFiles.size() != 0) {
                samplefile = samplingFiles[0];
                outputTypes["sample"].push_back(samplefile);
                calc = savedCalc;
                iters = savedIters;
            }else { return 0; }
        }
        
        if (format != "sharedfile") { processSingleSample(); } //handles multiple label values
        else { processSharedFile(); } //handles multiple label values and multiple samples
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	} return 0; }
        
        //set column file as new current columnfile
        itTypes = outputTypes.find("sample");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { string currentName = (itTypes->second)[0]; current->setSampleFile(currentName); }
        }
        
        m->mothurOut("\nOutput File Names: \n"); 
        for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "execute");
        exit(1);
    }
}
//**********************************************************************************************************************
int EstimatorSingleCommand::processSharedFile() {
    try {
        vector<string> Groups;
        InputData input(inputfile, format, Groups);
        set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
        SharedRAbundVectors* shared = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
        Groups = shared->getNamesGroups();
        
        if (outputdir == "") { outputdir += util.hasPath(inputfile); }
        string fileNameRoot = outputdir + util.getRootName(util.getSimpleName(inputfile));
        
        map<string, string> variables;
        variables["[filename]"] = fileNameRoot;
        variables["[distance]"] = shared->getLabel();
        if (util.inUsersGroups(calc, samplingCalcs)) {  variables["[distance]"] = shared->getLabel() + ".0"; }
        string outputFileName = getOutputFileName(calc, variables);
        
        vector<ofstream*> out;
        outputNames.push_back(outputFileName); outputTypes[calc].push_back(outputFileName);
        
        if (util.inUsersGroups(calc, samplingCalcs)) {
            out.resize(3);
            out[0] = new ofstream();
            util.openOutputFile(outputFileName, *out[0]); //format output
            out[0]->setf(ios::fixed, ios::floatfield); out[0]->setf(ios::showpoint);
            
            variables["[distance]"] = shared->getLabel() + ".1";
            string outputFileName1 = getOutputFileName(calc, variables);
            outputNames.push_back(outputFileName1); outputTypes[calc].push_back(outputFileName1);
            out[1] = new ofstream();
            util.openOutputFile(outputFileName1, *out[1]); //format output
            out[1]->setf(ios::fixed, ios::floatfield); out[1]->setf(ios::showpoint);
            
            variables["[distance]"] = shared->getLabel() + ".2";
            string outputFileName2 = getOutputFileName(calc, variables);
            outputNames.push_back(outputFileName2); outputTypes[calc].push_back(outputFileName2);
            out[2] = new ofstream();
            util.openOutputFile(outputFileName2, *out[2]); //format output
            out[2]->setf(ios::fixed, ios::floatfield); out[2]->setf(ios::showpoint);
            
            *out[0] << "#Be sure to use the correct sampling estimator with your calculator. IG is used for igabund and igrarefact. LN is used for lnabund, lnshift and lnrarefact. LS is used for lsabund and lsrarefaction. SI is used for siabund, sirarefact and sishift.\n";
            
            *out[1] << "#Be sure to use the correct sampling estimator with your calculator. IG is used for igabund and igrarefact. LN is used for lnabund, lnshift and lnrarefact. LS is used for lsabund and lsrarefaction. SI is used for siabund, sirarefact and sishift.\n";
            
            *out[2] << "#Be sure to use the correct sampling estimator with your calculator. IG is used for igabund and igrarefact. LN is used for lnabund, lnshift and lnrarefact. LS is used for lsabund and lsrarefaction. SI is used for siabund, sirarefact and sishift.\n";
            
        }else if (util.inUsersGroups(calc, rarefactCalcs)) {
            out.resize(1);
            out[0] = new ofstream();
            util.openOutputFile(outputFileName, *out[0]); //format output
            out[0]->setf(ios::fixed, ios::floatfield); out[0]->setf(ios::showpoint);
            
            *out[0] << "label\tgroup\t" << calc << "_Lower\t" << calc << "_Median\t" << calc << "_Upper\n";
        }else if (util.inUsersGroups(calc, abundCalcs)) {
            out.resize(1);
            out[0] = new ofstream();
            util.openOutputFile(outputFileName, *out[0]); //format output
            out[0]->setf(ios::fixed, ios::floatfield); out[0]->setf(ios::showpoint);
            
            *out[0] << "label\tgroup\tnum\t" << calc << "\n";
        }
        
        while (shared != nullptr) {
            
            if (m->getControl_pressed()) { delete shared; break; }
            
            processShared(shared, out, fileNameRoot); delete shared;
            
            shared = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
        }
        out[0]->close(); delete out[0];
        if (util.inUsersGroups(calc, samplingCalcs)) { out[1]->close(); out[2]->close(); delete out[1]; delete out[2]; }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "processSharedFile");
        exit(1);
    }
}
//**********************************************************************************************************************
int EstimatorSingleCommand::processSingleSample() {
    try {
        InputData input(inputfile, format, nullVector);
        set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
        SAbundVector* sabund = util.getNextSAbund(input, allLines, userLabels, processedLabels, lastLabel);
        
        if (outputdir == "") { outputdir += util.hasPath(inputfile); }
        string fileNameRoot = outputdir + util.getRootName(util.getSimpleName(inputfile));
        
        map<string, string> variables;
        variables["[filename]"] = fileNameRoot;
        variables["[distance]"] = sabund->getLabel();
        if (util.inUsersGroups(calc, samplingCalcs)) {  variables["[distance]"] = sabund->getLabel() + ".0"; }
        string outputFileName = getOutputFileName(calc, variables);
        
        vector<ofstream*> out;
        outputNames.push_back(outputFileName); outputTypes[calc].push_back(outputFileName);
        
        if (util.inUsersGroups(calc, samplingCalcs)) {
            
            out.resize(3);
            out[0] = new ofstream();
            util.openOutputFile(outputFileName, *out[0]); //format output
            out[0]->setf(ios::fixed, ios::floatfield); out[0]->setf(ios::showpoint);
            
            variables["[distance]"] = sabund->getLabel() + ".1";
            string outputFileName1 = getOutputFileName(calc, variables);
            outputNames.push_back(outputFileName1); outputTypes[calc].push_back(outputFileName1);
            out[1] = new ofstream();
            util.openOutputFile(outputFileName1, *out[1]); //format output
            out[1]->setf(ios::fixed, ios::floatfield); out[1]->setf(ios::showpoint);
            
            variables["[distance]"] = sabund->getLabel() + ".2";
            string outputFileName2 = getOutputFileName(calc, variables);
            outputNames.push_back(outputFileName2); outputTypes[calc].push_back(outputFileName2);
            out[2] = new ofstream();
            util.openOutputFile(outputFileName2, *out[2]); //format output
            out[2]->setf(ios::fixed, ios::floatfield); out[2]->setf(ios::showpoint);
            
            *out[0] << "#Be sure to use the correct sampling estimator with your calculator. IG is used for igabund and igrarefact. LN is used for lnabund, lnshift and lnrarefact. LS is used for lsabund and lsrarefaction. SI is used for siabund, sirarefact and sishift.\n";
            
            *out[1] << "#Be sure to use the correct sampling estimator with your calculator. IG is used for igabund and igrarefact. LN is used for lnabund, lnshift and lnrarefact. LS is used for lsabund and lsrarefaction. SI is used for siabund, sirarefact and sishift.\n";
            
            *out[2] << "#Be sure to use the correct sampling estimator with your calculator. IG is used for igabund and igrarefact. LN is used for lnabund, lnshift and lnrarefact. LS is used for lsabund and lsrarefaction. SI is used for siabund, sirarefact and sishift.\n";
        }else if (util.inUsersGroups(calc, rarefactCalcs)) {
            out.resize(1);
            out[0] = new ofstream();
            util.openOutputFile(outputFileName, *out[0]); //format output
            out[0]->setf(ios::fixed, ios::floatfield); out[0]->setf(ios::showpoint);
            
             *out[0] << "label\t" << calc << "_Lower\t" << calc << "_Median\t" << calc << "_Upper\n";
        }else if (util.inUsersGroups(calc, abundCalcs)) {
            out.resize(1);
            out[0] = new ofstream();
            util.openOutputFile(outputFileName, *out[0]); //format output
            out[0]->setf(ios::fixed, ios::floatfield); out[0]->setf(ios::showpoint);
            
            *out[0] << "label\tnum\t" << calc << "\n";
        }
         
        while (sabund != nullptr) {
                   
            if (m->getControl_pressed()) { delete sabund; break; }
                   
            processSingle(sabund, sabund->getLabel(), out, fileNameRoot); delete sabund;
                  
            sabund = util.getNextSAbund(input, allLines, userLabels, processedLabels, lastLabel);
        }
        
        out[0]->close(); delete out[0];
        if (util.inUsersGroups(calc, samplingCalcs)) { out[1]->close(); out[2]->close(); delete out[1]; delete out[2]; }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "processSingleSample");
        exit(1);
    }
}
//**********************************************************************************************************************
int EstimatorSingleCommand::processShared(SharedRAbundVectors*& shared, vector<ofstream*>& out, string fileRoot) {
    try {
        string outputFileName = "";
        string label = shared->getLabel();
        vector<string> groupNames = shared->getNamesGroups();
        vector< vector<double> > abundResults;
        int maxSize = 0;
        
        for (int i = 0; i < groupNames.size(); i++) {
            
            m->mothurOut("\nProcessing group " + groupNames[i] + "\n\n");
            
            string groupName = groupNames[i] + " " + label;
            SAbundVector* sabund = new SAbundVector(shared->getSAbundVector(groupNames[i]));
            
            if (m->getControl_pressed()) { delete sabund; break; }
            
            if (util.inUsersGroups(calc, samplingCalcs)) {
                
                *out[0] << "#" << groupName << endl;
                *out[1] << "#" << groupName << endl;
                *out[2] << "#" << groupName << endl;
                
                vector<string> outputFileNames = runSamplingCalcs(sabund, fileRoot);
                util.appendFiles(outputFileNames[0], *out[0]); util.mothurRemove(outputFileNames[0]);
                util.appendFiles(outputFileNames[1], *out[1]); util.mothurRemove(outputFileNames[1]);
                util.appendFiles(outputFileNames[2], *out[2]); util.mothurRemove(outputFileNames[2]);
                
            }else if (util.inUsersGroups(calc, rarefactCalcs)) {
                
                *out[0] << label << '\t'; if (groupNames[i] != "") { *out[0] << groupNames[i] << '\t'; }
                runRarefactCalcs(sabund->getNumSeqs(), groupName, *out[0]);
            
            }else if (util.inUsersGroups(calc, abundCalcs)) {
                vector<double> results = runAbundCalcs(sabund, groupName);
                
                abundResults.push_back(results);
                
                if (results.size() > maxSize) { maxSize = results.size(); }
            }
            
            delete sabund;
        }
        
        if (abundResults.size() != 0) {//ran an abund calc on several samples, combine results into one file
            //find smallest largest size
            
            for (int i = 0; i < abundResults.size(); i++) {
                if (abundResults[i].size() < maxSize) {
                    abundResults[i].resize(maxSize, -1); //fill blanks with NA
                }
            }
            
            
            for (int i = 0; i < maxSize; i++) {
                
                for (int j = 0; j < groupNames.size(); j++) {
                    
                    *out[0] << label << '\t' << groupNames[j] << '\t' << (i+1);
                    
                    if (abundResults[j][i] == -1) {  *out[0] << "\tNA" << endl;  }
                    else {  *out[0] << '\t' << abundResults[j][i] << endl;  }
                }
            }
            
        }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "processShared");
        exit(1);
    }
}

//**********************************************************************************************************************
int EstimatorSingleCommand::processSingle(SAbundVector*& sabund, string groupName, vector<ofstream*>& out, string fileRoot) {
    try {
        string label = sabund->getLabel();
        vector<string> outputFileNames;
        
        if (util.inUsersGroups(calc, rarefactCalcs)) {
            
            *out[0] << label << '\t';
            runRarefactCalcs(sabund->getNumSeqs(), groupName, *out[0]);
            
        }else if (util.inUsersGroups(calc, samplingCalcs)) {
            
            *out[0] << "#" << groupName << endl;
            *out[1] << "#" << groupName << endl;
            *out[2] << "#" << groupName << endl;
            vector<string> outputFileNames = runSamplingCalcs(sabund, fileRoot);
            util.appendFiles(outputFileNames[0], *out[0]); util.mothurRemove(outputFileNames[0]);
            util.appendFiles(outputFileNames[1], *out[1]); util.mothurRemove(outputFileNames[1]);
            util.appendFiles(outputFileNames[2], *out[2]); util.mothurRemove(outputFileNames[2]);
            
        }else if (util.inUsersGroups(calc, abundCalcs)) {
            vector<double> results = runAbundCalcs(sabund, groupName);
            
            for (int i = 0; i < results.size(); i++) {
                if (m->getControl_pressed()) { break; }
                
                *out[0] << label << '\t' << (i+1) << '\t' << results[i] << endl;
            }
        }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "process");
        exit(1);
    }
}
//**********************************************************************************************************************
int EstimatorSingleCommand::runRarefactCalcs(int numSeqs, string groupName, ofstream& out) {
    try {
        vector<double> results;
        vector<mcmcSample> thisGroupSample;
        
        if ((calc == "igrarefact") || (calc == "lnrarefact"))       { if (samplefile != "") { fillSampling(burn, burnSample); } }
        else if ((calc == "lsrarefact") || (calc == "sirarefact"))  { if (samplefile != "") { fillSampling(burn, burnSample, true); } }
        
        it = sampling.find(groupName);
        if (it != sampling.end()) { thisGroupSample = it->second; }
        else { m->mothurOut("[ERROR]: Unable to find sampling info for group " + groupName + ", quitting. Do you need to adjust the iters, burn or burnsample parameters?\n"); m->setControl_pressed(true); return 0; }
        
        DiversityCalculator* diversityCalc;
        if (calc == "igrarefact")       { diversityCalc = new IGRarefaction(coverage); }
        else if (calc == "lnrarefact")  { diversityCalc = new LNRarefaction(coverage); }
        else if (calc == "lsrarefact")  { diversityCalc = new LSRarefaction(coverage); }
        else if (calc == "sirarefact")  { diversityCalc = new SIRarefaction(coverage); }
       
        results = diversityCalc->getValues(numSeqs, sampling[groupName]);
        delete diversityCalc;
        
        for (int i = 0; i < results.size(); i++) { out << results[i] << '\t'; } out << endl;
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "runRarefactCalcs");
        exit(1);
    }
}
//**********************************************************************************************************************
vector<string> EstimatorSingleCommand::runSamplingCalcs(SAbundVector*& sabund, string fileRoot) {
    try {
        vector<string> resultFiles;
        
        map<string, string> variables;
        variables["[filename]"] = fileRoot;
        variables["[distance]"] = sabund->getLabel();
        variables["[tag]"] = calc;
        string outputFileStub = variables["[filename]"] + variables["[distance]"] + variables["[tag]"];
        
        DiversityCalculator* diversityCalc;
        if (calc == "ig")       { diversityCalc = new MetroIG(fitIters, sigmaAlpha, sigmaBeta, sigmaS, iters, outputFileStub);                 }
        else if (calc == "ln")  { diversityCalc = new MetroLogNormal(fitIters, sigmaAlpha, sigmaBeta, sigmaS, iters, outputFileStub);          }
        else if (calc == "ls")  { diversityCalc = new MetroLogStudent(fitIters, sigmaAlpha, sigmaBeta, sigmaN, sigmaS, iters, outputFileStub); }
        else if (calc == "si")  { diversityCalc = new MetroSichel(fitIters, sigmaAlpha, sigmaBeta, sigmaN, sigmaS, iters, outputFileStub);     }
        
        resultFiles = diversityCalc->getValues(sabund);
        
        if (m->getControl_pressed()) {
           m->mothurOut("\nHow do you create the sampling files?\n\nRun a short trial MCMC run of 1000 iterations with guessed std. dev.s for the proposal distributions say about 10% of the parameter values. Adjust the std. dev.s untill the acceptance ratios are about 0.5. Then perform a longer run of say 250,000 iterations (mothur's default). Three data files with posterior samples for three different sets of parameter values will be generated.\n\n");
        }
        delete diversityCalc;
        
        return resultFiles;
        
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "runSamplingCalcs");
        exit(1);
    }
}
//**********************************************************************************************************************
vector<double> EstimatorSingleCommand::runAbundCalcs(SAbundVector*& sabund, string groupName) {
    try {
        
        int maxRank = sabund->getMaxRank();
        int numSeqs = sabund->getNumSeqs();
        
        vector<double> results;
        vector<mcmcSample> thisGroupSample;
        
        if ((calc == "igabund") || (calc == "lnshift") || (calc == "lnabund"))       { if (samplefile != "") { fillSampling(burn, burnSample); } }
        else if ((calc == "siabund") || (calc == "sishift") || (calc == "lsabund"))  { if (samplefile != "") { fillSampling(burn, burnSample, true); } }
        
        if (calc != "erarefact") {
            it = sampling.find(groupName);
            if (it != sampling.end()) { thisGroupSample = it->second; }
            else { m->mothurOut("[ERROR]: Unable to find sampling info for group " + groupName + ", quitting. Do you need to adjust the iters, burn or burnsample parameters?\n"); m->setControl_pressed(true); return results; }
        }
        
        //convert freq percentage to number
        int increment = 1; if (freq < 1.0) {  increment = numSeqs * freq;  } else { increment = freq;  }
        
        DiversityCalculator* diversityCalc;
        if (calc == "igabund")       { diversityCalc = new IGAbundance();           }
        else if (calc == "lnshift")  { diversityCalc = new LNShift();               }
        else if (calc == "lnabund")  { diversityCalc = new LNAbundance();           }
        else if (calc == "siabund")  { diversityCalc = new SIAbundance();           }
        else if (calc == "sishift")  { diversityCalc = new SIShift();               }
        else if (calc == "lsabund")  { diversityCalc = new LSAbundance();           }
        else if (calc == "erarefact"){ diversityCalc = new ERarefaction(increment); }
        
        if (calc == "erarefact") {
            diversityCalc->getValues(sabund, results);
        }
        else if ((calc == "igabund") || (calc == "siabund") || (calc == "lnabund") || (calc == "lsabund"))       {
            results = diversityCalc->getValues(maxRank, sampling[groupName]);
        }
        else { //sishift, lnshift
            results = diversityCalc->getValues(numSeqs, sampling[groupName]);
        }
        
        delete diversityCalc;
        
        return results;
        
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "runAbundCalcs");
        exit(1);
    }
}
//**********************************************************************************************************************
int EstimatorSingleCommand::fillSampling(int burnValue, int burnSampleValue, bool filldNu) {
    try {
        sampling.clear();
        
        int numPiecesExpected = 5;
        if (filldNu) { numPiecesExpected = 6; }
        
        ifstream in; util.openInputFile(samplefile, in);
        
        util.getline(in); gobble(in); //grab header
        string groupName = "";
        
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            string line = util.getline(in); gobble(in);
            
            if (line != "") {
                
                if (line[0] == '#') {
                    groupName = line.substr(1); //looks like #groupName label ie. #F000D000 0.03
                }else {
            
                    vector<string> pieces; util.splitAtComma(line, pieces);
                    
                    if (pieces.size() == numPiecesExpected) {
                        
                        int sampleSize, ns;
                        util.mothurConvert(pieces[0], sampleSize);
                        
                        if ((sampleSize > burnValue) && (sampleSize % burnSampleValue == 0)) {
                            
                            double alpha = 0, beta = 0, dNu = 0;
                            if (!filldNu) {  util.mothurConvert(pieces[3], ns);  }
                            else {
                                util.mothurConvert(pieces[3], dNu);
                                util.mothurConvert(pieces[4], ns);
                            }
                            
                            util.mothurConvert(pieces[1], alpha);
                            util.mothurConvert(pieces[2], beta);
                            
                            mcmcSample entry(alpha, beta, dNu, ns);
                            sampling[groupName].push_back(entry);
                        }
                        
                    }else {
                        m->mothurOut("\n[WARNING]: Unexpected format in sampling file, ignoring. Expected " + toString(numPiecesExpected) + " values separated by commas, found " + toString(pieces.size()) + ". Expecting something like: '0,7.419861e-01,4.695223e+00,5773,337.552846' or 0,-1.787922,6.348652,4784302.925302,8806,331.214377 for each line.\n\n");
                        sampling.clear(); break;
                    }
                }
            }
        }
        in.close();
        
        return ((int)sampling.size());
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "fillSampling");
        exit(1);
    }
}
//**********************************************************************************************************************

