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
        CommandParameter pcalc("calc", "Multiple", "erarefaction-metroig-metroln-metrols-metrosichel-igabund-igrarefaction-lnrarefaction-lnabund-lnshift", "erarefaction", "", "", "","",true,false,true); parameters.push_back(pcalc); //lnabund
        CommandParameter pabund("abund", "Number", "", "10", "", "", "","",false,false); parameters.push_back(pabund);
        CommandParameter palpha("sigmaa", "Number", "", "0.1", "", "", "","",false,false,true); parameters.push_back(palpha);
        CommandParameter pbeta("sigmab", "Number", "", "0.1", "", "", "","",false,false); parameters.push_back(pbeta);
        CommandParameter psigman("sigman", "Number", "", "0.1", "", "", "","",false,false); parameters.push_back(psigman);
        CommandParameter psigmas("sigmas", "Number", "", "100", "", "", "","",false,false); parameters.push_back(psigmas);
        CommandParameter pburn("burn", "Number", "", "2000000", "", "", "","",false,false); parameters.push_back(pburn);
        CommandParameter pcoverage("coverage", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pcoverage);
        CommandParameter psamplenum("burnsample", "Number", "", "1000", "", "", "","",false,false); parameters.push_back(psamplenum);
        CommandParameter piters("iters", "Number", "", "1000", "", "", "","",false,false); parameters.push_back(piters);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
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
        string helpString = "";
        ValidCalculators validCalculator;
        helpString += "The estimator.single command parameters are " + getCommandParameters() + ".\n";
        helpString += "The estimator.single command should be in the following format: \n";
        helpString += "estimator.single(list=yourListFile, calc=yourEstimators).\n";
        helpString += "Example estimator.single(list=final.opti_mcc.list, calc=erarefaction).\n";
        helpString += "The freq parameter is used indicate when to output your data, by default it is set to 100. But you can set it to a percentage of the number of sequence. For example freq=0.10, means 10%. \n";
        helpString += "The sample file is used to provide mcmc sampling to the following calculators: ...............\n";
        helpString += "The default values for freq is 100, and calc is erarefaction.\n";
        helpString += "The sigmaa parameter is used to set the std. dev. of alpha / X / mean prop. distn for MetroIG / MetroLogNormal / MetroLogStudent / MetroSichel, respectively. Default = 0.10. n";
        helpString += "The sigmab parameter is used to set the std. dev. of beta / Y / V prop. distn for MetroIG / MetroLogNormal / MetroLogStudent / MetroSichel, respectively. Default = 0.10. n";
        helpString += "The sigman parameter is used to set the std. dev. of N / Gamma prop. distn for MetroLogStudent / MetroSichel, respectively. Default = 0.10. n";
        helpString += "The sigmas parameter is used to set the std. dev. of S prop. distn for MetroIG / MetroLogNormal / MetroLogStudent / MetroSichel. Default = 100. n";
        helpString += "The coverage parameter allows you to the desired coverage.  It is required for the igrarefaction calculator.\n";
        helpString += "The iters parameter allows you to set number of mcmc samples to generate.  The default is 1000.\n";
        helpString += "The burn parameter allows ignore part of the sampling file.  Default = 200000 / 100000 for IGAbundance, LNShift / IGRarefaction, LNRarefaction, respectively.\n";
        helpString += "The burnsample parameter allows you to set sampling frequency.  The default is 1000 / 100 for IGAbundance, LNShift / IGRarefaction, LNRarefaction, respectively.\n";
        helpString += validCalculator.printCalc("single");
        helpString += "The label parameter is used to analyze specific labels in your input.\n";
        
        return helpString;
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "getHelpString");
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
EstimatorSingleCommand::EstimatorSingleCommand(){
    try {
        abort = true; calledHelp = true;
        setParameters();
        vector<string> tempOutNames;
        outputTypes["erarefaction"] = tempOutNames;
        outputTypes["igrarefaction"] = tempOutNames;
        outputTypes["igabund"] = tempOutNames;
        outputTypes["lnabund"] = tempOutNames;
        outputTypes["lnrarefaction"] = tempOutNames;
        outputTypes["lnshift"] = tempOutNames;
        outputTypes["metroig"] = tempOutNames;
        outputTypes["metroln"] = tempOutNames;
        outputTypes["metrols"] = tempOutNames;
        outputTypes["metrosichel"] = tempOutNames;
        
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "EstimatorSingleCommand");
        exit(1);
    }
}
//**********************************************************************************************************************
EstimatorSingleCommand::EstimatorSingleCommand(string option)  {
    try {
        abort = false; calledHelp = false;
        allLines = true;
        
        //allow user to run help
        if(option == "help") { help(); calledHelp = true; abort = true; }
        else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        
        else {
            vector<string> myArray = setParameters();
            
            OptionParser parser(option);
            map<string,string> parameters = parser.getParameters();
            map<string,string>::iterator it;
            
            ValidParameters validParameter;
            
            //check to make sure all parameters are valid for command
            for (it = parameters.begin(); it != parameters.end(); it++) {
                if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
            }
            
            //initialize outputTypes
            vector<string> tempOutNames;
            outputTypes["erarefaction"] = tempOutNames;
            outputTypes["igabund"] = tempOutNames;
            outputTypes["lnabund"] = tempOutNames;
            outputTypes["metroig"] = tempOutNames;
            outputTypes["metroln"] = tempOutNames;
            outputTypes["metrols"] = tempOutNames;
            outputTypes["metrosichel"] = tempOutNames;
            outputTypes["igrarefaction"] = tempOutNames;
            outputTypes["lnrarefaction"] = tempOutNames;
            outputTypes["lnshift"] = tempOutNames;
            
            //if the user changes the input directory command factory will send this info to us in the output parameter
            string inputDir = validParameter.valid(parameters, "inputdir");
            if (inputDir == "not found"){	inputDir = "";		}
            else {
                string path;
                it = parameters.find("shared");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["shared"] = inputDir + it->second;		}
                }
                
                it = parameters.find("rabund");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["rabund"] = inputDir + it->second;		}
                }
                
                it = parameters.find("sabund");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["sabund"] = inputDir + it->second;		}
                }
                
                it = parameters.find("list");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["list"] = inputDir + it->second;		}
                }
                
                it = parameters.find("sample");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["sample"] = inputDir + it->second;		}
                }
            }
            
            //check for required parameters
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
            else { hasSample = true; }
            
            
            //if the user changes the output directory command factory will send this info to us in the output parameter
            outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";		}
            
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
            calc = validParameter.valid(parameters, "calc");
            if (calc == "not found") { calc = "erarefaction";  }
            else {
                if (calc == "default")  {  calc = "erarefaction";  }
            }
            util.splitAtDash(calc, Estimators);
            if (util.inUsersGroups("citation", Estimators)) {
                validCalculator.printCitations(Estimators);
                //remove citation from list of calcs
                for (int i = 0; i < Estimators.size(); i++) { if (Estimators[i] == "citation") {  Estimators.erase(Estimators.begin()+i); break; } }
            }
            
            set<string> estimatorsThatRequireSampleFile;
            estimatorsThatRequireSampleFile.insert("igabund");
            estimatorsThatRequireSampleFile.insert("igrarefaction");
            //estimatorsThatRequireSampleFile.insert("lnabund");
            estimatorsThatRequireSampleFile.insert("lnrarefaction");
            
            //remove any typo calcs
            vector<string> validEstimates;
            for (int i=0; i<Estimators.size(); i++) {
                if (validCalculator.isValidCalculator("estimator", Estimators[i]) ) {
                    
                    bool ignore = false;
                    if (!hasSample) { //if you didn't provide a mcmc sample file, but are trying to run a calc that needs it, then ignore
                        if (estimatorsThatRequireSampleFile.count(Estimators[i]) != 0) { ignore = true; }
                    }
                    
                    if (!ignore) { validEstimates.push_back(Estimators[i]); }
                    else { m->mothurOut("[WARNING]: " + Estimators[i] + " requires a mcmc sampling file and you have not provided one, ignoring estimator. You can produce a sampling file using the metroig, metroln, metrols or metrosichel calculator.\n"); }
                }
            }
            Estimators = validEstimates;
            if (Estimators.size() == 0) { abort = true; m->mothurOut("[ERROR]: no valid estimators, aborting.\n"); }
            else {  sort(Estimators.begin(), Estimators.end()); }

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
            
            temp = validParameter.valid(parameters, "iters");		if (temp == "not found") { temp = "1000"; }
            util.mothurConvert(temp, iters);
            
            burnSet = false;
            temp = validParameter.valid(parameters, "burn");		if (temp == "not found") { temp = "2000000"; }else { burnSet = true; }
            util.mothurConvert(temp, burn);
            
            if (burnSet) { //user did not set the parameter
                if ((util.inUsersGroups("lnrarefaction", Estimators)) && (util.inUsersGroups("igrarefaction", Estimators)) && (util.inUsersGroups("igabund", Estimators))) {
                    m->mothurOut("[WARNING]: You set the burn parameter, and the igrarefaction and igabund calulators have different default values. IGAbund burnsample default is 2000000, but IGRarefaction and LNRarection's default is 100000. Are you sure you meant to set them to the same value? If so, ignore this warning.\n");
                }
            }
            
            temp = validParameter.valid(parameters, "coverage");		if (temp == "not found") { temp = "-1"; }
            util.mothurConvert(temp, coverage);
            
            if ((util.isEqual(coverage, -1)) && ((util.inUsersGroups("igrarefaction", Estimators)) || (util.inUsersGroups("lnrarefaction", Estimators)))) {
                m->mothurOut("[ERROR]: You must set the coverage parameter to run the igrarefaction or lnrarefaction estimator. Aborting.\n"); abort=true;
            }
            
            burnSampleSet = false;
            temp = validParameter.valid(parameters, "burnsample");		if (temp == "not found") { temp = "1000"; }else { burnSampleSet = true; }
            util.mothurConvert(temp, burnSample);
            
            if (burnSampleSet) { //user did not set the parameter
                if ((util.inUsersGroups("lnrarefaction", Estimators)) && (util.inUsersGroups("igrarefaction", Estimators)) && (util.inUsersGroups("igabund", Estimators))) {
                    m->mothurOut("[WARNING]: You set the burnsample parameter, and the igrarefaction and igabund calulators have different default values. IGAbund burnsample default is 1000, but IGRarefaction and LNRarection's default is 100. Are you sure you meant to set them to the same value? If so, ignore this warning.\n");
                }
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
        
        vector<string> inputFileNames;
        if (format != "sharedfile") { inputFileNames.push_back(inputfile);  }
        else {  inputFileNames = parseSharedFile(sharedfile);  format = "rabund"; }
        
        for (int p = 0; p < inputFileNames.size(); p++) {
            
            if (m->getControl_pressed()) {  outputTypes.clear(); for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	}    return 0; }
            
            if (outputDir == "") { outputDir += util.hasPath(inputFileNames[p]); }
            string fileNameRoot = outputDir + util.getRootName(util.getSimpleName(inputFileNames[p]));
            
            if (inputFileNames.size() > 1) { m->mothurOut("\nProcessing group " + groups[p] + "\n\n"); }
            
            InputData input(inputFileNames[p], format, nullVector);
            SAbundVector* sabund = input.getSAbundVector();
            string lastLabel = sabund->getLabel();
            
            //if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
            set<string> processedLabels;
            set<string> userLabels = labels;
            
            while((sabund != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
                
                if (m->getControl_pressed()) { delete sabund; break; }
                
                if(allLines == 1 || labels.count(sabund->getLabel()) == 1){
                    
                    m->mothurOut(sabund->getLabel() + "\n"); processedLabels.insert(sabund->getLabel()); userLabels.erase(sabund->getLabel());
                    
                    process(sabund, fileNameRoot);
                }
                //you have a label the user want that is smaller than this label and the last label has not already been processed
                if ((util.anyLabelsToProcess(sabund->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
                    string saveLabel = sabund->getLabel();
                    
                    delete sabund;
                    sabund = (input.getSAbundVector(lastLabel));
                    
                    m->mothurOut(sabund->getLabel() + "\n"); processedLabels.insert(sabund->getLabel()); userLabels.erase(sabund->getLabel());
                    
                    process(sabund, fileNameRoot);
                    
                    //restore real lastlabel to save below
                    sabund->setLabel(saveLabel);
                }
                
                lastLabel = sabund->getLabel();
                
                delete sabund;
                sabund = input.getSAbundVector();
            }
            
            
            if (m->getControl_pressed()) {  for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	} outputTypes.clear(); return 0; }
            
            //output error messages about any remaining user labels
            set<string>::iterator it;
            bool needToRun = false;
            for (it = userLabels.begin(); it != userLabels.end(); it++) {  
                m->mothurOut("Your file does not include the label " + *it); 
                if (processedLabels.count(lastLabel) != 1)  { m->mothurOut(". I will use " + lastLabel + ".\n"); needToRun = true;  }
                else                                        {  m->mothurOut(". Please refer to " + lastLabel + ".\n");              }
            }
            
            //run last label if you need to
            if (needToRun)  {
                if (sabund != NULL) {	delete sabund;	}
                sabund = input.getSAbundVector(lastLabel);
                
                m->mothurOut(sabund->getLabel() + "\n");
                
                process(sabund, fileNameRoot); delete sabund;
            }
        }
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	} return 0; }
        
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
int EstimatorSingleCommand::process(SAbundVector*& sabund, string fileRoot) {
    try {
        
        for (int i = 0; i < Estimators.size(); i++) {
            
            if (m->getControl_pressed()) { break; }
            
            if (Estimators[i] == "erarefaction")         { runErarefaction(sabund, fileRoot);    }
            else if (Estimators[i] == "metroig")         { runMetroIG(sabund, fileRoot);         }
            else if (Estimators[i] == "metroln")         { runMetroLogNormal(sabund, fileRoot);  }
            else if (Estimators[i] == "metrols")         { runMetroLogStudent(sabund, fileRoot); }
            else if (Estimators[i] == "metrosichel")     { runMetroSichel(sabund, fileRoot);     }
            else if (Estimators[i] == "igabund")         { runIGAbund(sabund, fileRoot);         }
            else if (Estimators[i] == "igrarefaction")   { runIGRarefaction(sabund, fileRoot);   }
            else if (Estimators[i] == "lnabund")         { runLNAbund(sabund, fileRoot);         }
            else if (Estimators[i] == "lnrarefaction")   { runLNRarefaction(sabund, fileRoot);   }
            else if (Estimators[i] == "lnshift")         { runLNShift(sabund, fileRoot);         }
        }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "process");
        exit(1);
    }
}
//**********************************************************************************************************************
string EstimatorSingleCommand::runErarefaction(SAbundVector*& sabund, string fileRoot) {
    try {
        map<string, string> variables;
        variables["[filename]"] = fileRoot;
        variables["[distance]"] = sabund->getLabel();
        string outputFileName = getOutputFileName("erarefaction", variables);
        outputNames.push_back(outputFileName); outputTypes["erarefaction"].push_back(outputFileName);
        
        ofstream out; util.openOutputFile(outputFileName, out); //format output
        out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
        
        out << "NumSampled\tERarefaction\n";
        
        ERarefaction erare;
        
        int numSeqs = sabund->getNumSeqs();
        
        //convert freq percentage to number
        int increment = 1; if (freq < 1.0) {  increment = numSeqs * freq;  } else { increment = freq;  }
        
        for (int i = 1; i < numSeqs; i++) {
            if((i % increment) == 0){
                double result = erare.getValues(sabund, i);
                out << i << '\t' << result << endl;
            }
        }
        
        double result = erare.getValues(sabund, numSeqs);
        out << numSeqs << '\t' << result << endl;
        
        out.close();
        
        return outputFileName;
        
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "runErarefaction");
        exit(1);
    }
}
//**********************************************************************************************************************
string EstimatorSingleCommand::runIGRarefaction(SAbundVector*& sabund, string fileRoot) {
    try {
        map<string, string> variables;
        variables["[filename]"] = fileRoot;
        variables["[distance]"] = sabund->getLabel();
        string outputFileName = getOutputFileName("igrarefaction", variables);
        outputNames.push_back(outputFileName); outputTypes["igrarefaction"].push_back(outputFileName);
        
        ofstream out; util.openOutputFile(outputFileName, out); //format output
        out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
        
        out << "IGRarefaction_Lower\tIGRarefaction_Median\tIGRarefaction_Upper\n";
        
        IGRarefaction igRare(coverage);
        
        if (samplefile != "") {
            int burnValue = burn;
            if (!burnSet) { burnValue = 100000; }
            
            int burnSampleValue = burnSample;
            if (!burnSampleSet) { burnSampleValue = 100; }
            
            fillSampling(burnValue, burnSampleValue);
        }
        
        vector<double> results = igRare.getValues(sabund, sampling);
        
        for (int i = 0; i < results.size(); i++) {  out << results[i] << '\t';  }
        out << endl;
        
        out.close();
        
        return outputFileName;
        
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "runIGRarefaction");
        exit(1);
    }
}
//**********************************************************************************************************************
string EstimatorSingleCommand::runLNRarefaction(SAbundVector*& sabund, string fileRoot) {
    try {
        map<string, string> variables;
        variables["[filename]"] = fileRoot;
        variables["[distance]"] = sabund->getLabel();
        string outputFileName = getOutputFileName("lnrarefaction", variables);
        outputNames.push_back(outputFileName); outputTypes["lnrarefaction"].push_back(outputFileName);
        
        ofstream out; util.openOutputFile(outputFileName, out); //format output
        out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
        
        out << "LNRarefaction_Lower\tILNRarefaction_Median\tLNRarefaction_Upper\n";
        
        LNRarefaction lnRare(coverage);
        
        if (samplefile != "") {
            int burnValue = burn;
            if (!burnSet) { burnValue = 100000; }
            
            int burnSampleValue = burnSample;
            if (!burnSampleSet) { burnSampleValue = 100; }
            
            fillSampling(burnValue, burnSampleValue);
        }
        
        vector<double> results = lnRare.getValues(sabund, sampling);
        
        for (int i = 0; i < results.size(); i++) {  out << results[i] << '\t';  }
        out << endl;
        
        out.close();
        
        return outputFileName;
        
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "runLNRarefaction");
        exit(1);
    }
}

//**********************************************************************************************************************
string EstimatorSingleCommand::runMetroIG(SAbundVector*& sabund, string fileRoot) {
    try {
        map<string, string> variables;
        variables["[filename]"] = fileRoot;
        variables["[distance]"] = sabund->getLabel();
        variables["[tag]"] = ".metroig";
        string outputFileStub = variables["[filename]"] + variables["[distance]"] + variables["[tag]"];
        
        MetroIG metroIG(sigmaAlpha, sigmaBeta, sigmaS, iters, outputFileStub);
        
        vector<string> resultFiles = metroIG.getValues(sabund);
        
        for (int i = 0; i < resultFiles.size(); i++) { outputNames.push_back(resultFiles[i]); outputTypes["metroig"].push_back(resultFiles[i]); }
        
        return outputFileStub;
        
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "runMetroIG");
        exit(1);
    }
}
//**********************************************************************************************************************
string EstimatorSingleCommand::runMetroLogNormal(SAbundVector*& sabund, string fileRoot) {
    try {
        map<string, string> variables;
        variables["[filename]"] = fileRoot;
        variables["[distance]"] = sabund->getLabel();
        variables["[tag]"] = ".metroln";
        string outputFileStub = variables["[filename]"] + variables["[distance]"] + variables["[tag]"];

        MetroLogNormal metroLN(sigmaAlpha, sigmaBeta, sigmaS, iters, outputFileStub);
        
        vector<string> resultFiles = metroLN.getValues(sabund);
        
        for (int i = 0; i < resultFiles.size(); i++) { outputNames.push_back(resultFiles[i]); outputTypes["metroln"].push_back(resultFiles[i]); }
        
        return outputFileStub;
        
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "runMetroLogNormal");
        exit(1);
    }
}
//**********************************************************************************************************************
string EstimatorSingleCommand::runMetroLogStudent(SAbundVector*& sabund, string fileRoot) {
    try {
        map<string, string> variables;
        variables["[filename]"] = fileRoot;
        variables["[distance]"] = sabund->getLabel();
        variables["[tag]"] = ".metrols";
        string outputFileStub = variables["[filename]"] + variables["[distance]"] + variables["[tag]"];
        
        MetroLogStudent metroLS(sigmaAlpha, sigmaBeta, sigmaN, sigmaS, iters, outputFileStub);
        
        vector<string> resultFiles = metroLS.getValues(sabund);
        
        for (int i = 0; i < resultFiles.size(); i++) { outputNames.push_back(resultFiles[i]); outputTypes["metrols"].push_back(resultFiles[i]); }
        
        return outputFileStub;
        
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "runMetroLogStudent");
        exit(1);
    }
}
//**********************************************************************************************************************
string EstimatorSingleCommand::runMetroSichel(SAbundVector*& sabund, string fileRoot) {
    try {
        map<string, string> variables;
        variables["[filename]"] = fileRoot;
        variables["[distance]"] = sabund->getLabel();
        variables["[tag]"] = ".metrosichel";
        string outputFileStub = variables["[filename]"] + variables["[distance]"] + variables["[tag]"];
        
        MetroSichel metroSichel(sigmaAlpha, sigmaBeta, sigmaN, sigmaS, iters, outputFileStub);
        
        vector<string> resultFiles = metroSichel.getValues(sabund);
        
        for (int i = 0; i < resultFiles.size(); i++) { outputNames.push_back(resultFiles[i]); outputTypes["metrosichel"].push_back(resultFiles[i]); }
        
        return outputFileStub;
        
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "runMetroSichel");
        exit(1);
    }
}
//**********************************************************************************************************************
int EstimatorSingleCommand::runIGAbund(SAbundVector*& sabund, string fileRoot) {
    try {
        map<string, string> variables;
        variables["[filename]"] = fileRoot;
        variables["[distance]"] = sabund->getLabel();
        string outputFileName = getOutputFileName("igabund", variables);
        outputNames.push_back(outputFileName); outputTypes["igabund"].push_back(outputFileName);
        
        ofstream out; util.openOutputFile(outputFileName, out); //format output
        out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
        
        if (samplefile != "") {
            int burnValue = burn;
            if (!burnSet) { burnValue = 200000; }
            
            int burnSampleValue = burnSample;
            if (!burnSampleSet) { burnSampleValue = 1000; }
            
            fillSampling(burnValue, burnSampleValue);
        }
        
        IGAbundance igAbund;
        
        vector<double> results = igAbund.getValues(sabund, sampling);
        
        for (int i = 0; i < results.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            out << i+1 << '\t' << results[i] << endl;
        }
        out.close();
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "runIGAbund");
        exit(1);
    }
}
//**********************************************************************************************************************
int EstimatorSingleCommand::runLNAbund(SAbundVector*& sabund, string fileRoot) {
    try {
        map<string, string> variables;
        variables["[filename]"] = fileRoot;
        variables["[distance]"] = sabund->getLabel();
        string outputFileName = getOutputFileName("igabund", variables);
        outputNames.push_back(outputFileName); outputTypes["igabund"].push_back(outputFileName);
        
        ofstream out; util.openOutputFile(outputFileName, out); //format output
        out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
        
        if (samplefile != "") {
            int burnValue = burn;
            if (!burnSet) { burnValue = 200000; }
            
            int burnSampleValue = burnSample;
            if (!burnSampleSet) { burnSampleValue = 1000; }
            
            fillSampling(burnValue, burnSampleValue);
        }
        
        IGAbundance igAbund;
        
        vector<double> results = igAbund.getValues(sabund, sampling);
        
        for (int i = 0; i < results.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            out << i+1 << '\t' << results[i] << endl;
        }
        out.close();
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "runLNAbund");
        exit(1);
    }
}
//**********************************************************************************************************************
int EstimatorSingleCommand::runLNShift(SAbundVector*& sabund, string fileRoot) {
    try {
        map<string, string> variables;
        variables["[filename]"] = fileRoot;
        variables["[distance]"] = sabund->getLabel();
        string outputFileName = getOutputFileName("lnshift", variables);
        outputNames.push_back(outputFileName); outputTypes["lnshift"].push_back(outputFileName);
        
        ofstream out; util.openOutputFile(outputFileName, out); //format output
        out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
        
        if (samplefile != "") {
            int burnValue = burn;
            if (!burnSet) { burnValue = 200000; }
            
            int burnSampleValue = burnSample;
            if (!burnSampleSet) { burnSampleValue = 1000; }
            
            fillSampling(burnValue, burnSampleValue);
        }
        
        LNShift lnShift;
        
        vector<double> results = lnShift.getValues(sabund, sampling);
        
        for (int i = 0; i < results.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            out << i+1 << '\t' << results[i] << endl;
        }
        out.close();
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "runLNAbund");
        exit(1);
    }
}
//**********************************************************************************************************************
vector<string> EstimatorSingleCommand::parseSharedFile(string filename) {
    try {
        vector<string> filenames;
        
        map<string, string> files;
        map<string, string>::iterator it3;
        
        InputData input(filename, "sharedfile", groups);
        SharedRAbundVectors* shared = input.getSharedRAbundVectors();
        
        string sharedFileRoot = util.getRootName(filename);
        groups = shared->getNamesGroups();
        
        //clears file before we start to write to it below
        for (int i=0; i<groups.size(); i++) {
            ofstream temp;
            string group = groups[i];
            util.openOutputFile((sharedFileRoot + group + ".rabund"), temp);
            temp.close();
            filenames.push_back((sharedFileRoot + group + ".rabund"));
            files[group] = (sharedFileRoot + group + ".rabund");
        }
        
        while(shared != NULL) {
            
            vector<SharedRAbundVector*> lookup = shared->getSharedRAbundVectors();
            for (int i = 0; i < lookup.size(); i++) {
                ofstream temp;
                string group = groups[i];
                util.openOutputFileAppend(files[group], temp);
                lookup[i]->getRAbundVector().print(temp);
                temp.close();
            }
            
            for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) { delete lookup[i]; } lookup[i] = NULL; }
            shared = input.getSharedRAbundVectors();
        }
        return filenames;
    }
    catch(exception& e) {
        m->errorOut(e, "EstimatorSingleCommand", "parseSharedFile");
        exit(1);
    }
}
//**********************************************************************************************************************
int EstimatorSingleCommand::fillSampling(int burnValue, int burnSampleValue) {
    try {
        sampling.clear();
        
        ifstream in;
        util.openInputFile(samplefile, in);
        
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            string line = util.getline(in); util.gobble(in);
            
            vector<string> pieces; util.splitAtComma(line, pieces);
            
            if (pieces.size() == 5) {
                
                int sampleSize, ns;
                util.mothurConvert(pieces[0], sampleSize);
                
                if ((sampleSize > burnValue) && (sampleSize % burnSampleValue == 0)) {
                
                    util.mothurConvert(pieces[3], ns);
                    
                    double alpha, beta;
                    util.mothurConvert(pieces[1], alpha);
                    util.mothurConvert(pieces[2], beta);
                    
                    mcmcSample entry(alpha, beta, ns);
                    sampling.push_back(entry);
                }
                
            }else {
                m->mothurOut("\n[WARNING]: Unexpected format in sampling file, ignoring. Expecting something like: '0,7.419861e-01,4.695223e+00,5773,337.552846' for each line.\n\n");
                sampling.clear(); break;
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

