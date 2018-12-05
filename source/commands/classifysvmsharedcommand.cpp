//
//  classifysvmsharedcommand.cpp
//  Mothur
//
//  Created by Joshua Lynch on 6/28/2013.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "classifysvmsharedcommand.h"

//**********************************************************************************************************************
vector<string> ClassifySvmSharedCommand::setParameters() {
    try {
        //CommandParameter pprocessors("processors", "Number", "", "1", "", "", "",false,false); parameters.push_back(pprocessors);
        CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none", "summary", false, true, true);
        parameters.push_back(pshared);
        CommandParameter pdesign("design", "InputTypes", "", "", "none", "none", "none", "", false, true, true);
        parameters.push_back(pdesign);

        // RFE or classification?
        // mode should be either 'rfe' or 'classify'
        CommandParameter mode("mode", "String", "", "", "", "", "", "", false, false);
        parameters.push_back(mode);

        // cross validation parameters
        CommandParameter evaluationFoldCountParam("evaluationfolds", "Number", "", "3", "", "", "", "", false, false);
        parameters.push_back(evaluationFoldCountParam);
        CommandParameter trainingFoldCountParam("trainingfolds", "Number", "", "10", "", "", "", "", false, false);
        parameters.push_back(trainingFoldCountParam);

        CommandParameter smoc("smoc", "Number", "", "3", "", "", "", "", false, false);
        parameters.push_back(smoc);

        // Support Vector Machine parameters
        CommandParameter kernelParam("kernel", "String", "", "", "", "", "", "", false, false);
        parameters.push_back(kernelParam);

        // data transformation parameters
        // transform should be 'zeroone' or 'zeromean' ('zeromean' is default)
        CommandParameter transformParam("transform", "String", "", "", "", "", "", "", false, false);
        parameters.push_back(transformParam);

        CommandParameter verbosityParam("verbose", "Number", "", "0", "", "", "", "", false, false);
        parameters.push_back(verbosityParam);


        // want this parameter to behave like the one in classify.rf
        CommandParameter pstdthreshold("stdthreshold", "Number", "", "0.0", "", "", "", "", false, false);
        parameters.push_back(pstdthreshold);
        // pruning params end

        CommandParameter pgroups("groups", "String", "", "", "", "", "", "", false, false);
        parameters.push_back(pgroups);
        CommandParameter plabel("label", "String", "", "", "", "", "", "", false, false);
        parameters.push_back(plabel);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "", "", false, false);
        parameters.push_back(pinputdir);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "", "", false, false);
        parameters.push_back(poutputdir);

        vector<string> myArray;
        for (int i = 0; i < parameters.size(); i++) {
            myArray.push_back(parameters[i].name);
        }
        return myArray;
    }
    catch (exception& e) {
        m->errorOut(e, "ClassifySvmSharedCommand", "setParameters");
        exit(1);
    }
}
//**********************************************************************************************************************
string ClassifySvmSharedCommand::getHelpString() {
    try {
        string helpString = "";
        helpString += "The classifysvm.shared command allows you to ....\n";
        helpString += "The classifysvm.shared command parameters are: shared, design, label, groups.\n";
        helpString += "The label parameter is used to analyze specific labels in your input.\n";
        helpString +=
                "The groups parameter allows you to specify which of the groups in your designfile you would like analyzed.\n";
        helpString += "The classifysvm.shared should be in the following format: \n";
        helpString += "classifysvm.shared(shared=yourSharedFile, design=yourDesignFile)\n";
        return helpString;
    }
    catch (exception& e) {
        m->errorOut(e, "ClassifySvmSharedCommand", "getHelpString");
        exit(1);
    }
}
//**********************************************************************************************************************
string ClassifySvmSharedCommand::getOutputPattern(string type) {
    try {
        string pattern = "";

        if (type == "summary") {
            pattern = "[filename],[distance],summary";
        } //makes file like: amazon.0.03.fasta
        else {
            m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n");
            m->setControl_pressed(true);
        }

        return pattern;
    }
    catch (exception& e) {
        m->errorOut(e, "ClassifySvmSharedCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************

ClassifySvmSharedCommand::ClassifySvmSharedCommand() {
    try {
        abort = true;
        calledHelp = true;
        setParameters();
        vector<string> tempOutNames;
        outputTypes["summary"] = tempOutNames;
    }
    catch (exception& e) {
        m->errorOut(e, "ClassifySvmSharedCommand", "ClassifySvmSharedCommand");
        exit(1);
    }
}

//**********************************************************************************************************************
ClassifySvmSharedCommand::ClassifySvmSharedCommand(string option) {
    try {
        abort = false;
        calledHelp = false;
        allLines = 1;

        //allow user to run help
        if (option == "help") {
            help();
            abort = true;
            calledHelp = true;
        }
        else if (option == "citation") {
            citation();
            abort = true;
            calledHelp = true;
        }
        else {
            //valid parameters for this command
            vector<string> myArray = setParameters();

            OptionParser parser(option);
            map<string, string> parameters = parser.getParameters();

            ValidParameters validParameter;
            map<string, string>::iterator it;
            //check to make sure all parameters are valid for command
            for (it = parameters.begin(); it != parameters.end(); it++) {
                if (!validParameter.isValidParameter(it->first, myArray, it->second)) {
                    abort = true;
                }
            }
            vector<string> tempOutNames;
            outputTypes["summary"] = tempOutNames;

            //if the user changes the input directory command factory will send this info to us in the output parameter
            string inputDir = validParameter.valid(parameters, "inputdir");
            if (inputDir == "not found") {
                inputDir = "";
            }
            else {
                string path;
                it = parameters.find("shared");
                //user has given a shared file
                if (it != parameters.end()) {
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {
                        parameters["shared"] = inputDir + it->second;
                    }
                }

                it = parameters.find("design");
                //user has given a design file
                if (it != parameters.end()) {
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {
                        parameters["design"] = inputDir + it->second;
                    }
                }

            }
            //check for parameters
            //get shared file, it is required
            sharedfile = validParameter.validFile(parameters, "shared");
            if (sharedfile == "not open") {
                sharedfile = "";
                abort = true;
            }
            else if (sharedfile == "not found") {
                //if there is a current shared file, use it
                sharedfile = current->getSharedFile();
                if (sharedfile != "") {
                    m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.");
                    m->mothurOutEndLine();
                }
                else {
                    m->mothurOut("You have no current sharedfile and the shared parameter is required.");
                    m->mothurOutEndLine();
                    abort = true;
                }
            }
            else {
                current->setSharedFile(sharedfile);
            }

            //get design file, it is required
            designfile = validParameter.validFile(parameters, "design");
            if (designfile == "not open") {
                sharedfile = "";
                abort = true;
            }
            else if (designfile == "not found") {
                //if there is a current shared file, use it
                designfile = current->getDesignFile();
                if (designfile != "") {
                    m->mothurOut("Using " + designfile + " as input file for the design parameter.");
                    m->mothurOutEndLine();
                }
                else {
                    m->mothurOut("You have no current designfile and the design parameter is required.");
                    m->mothurOutEndLine();
                    abort = true;
                }
            }
            else {
                current->setDesignFile(designfile);
            }

            //if the user changes the output directory command factory will send this info to us in the output parameter
            outputDir = validParameter.valid(parameters, "outputdir");
            if (outputDir == "not found") {
                outputDir = util.hasPath(sharedfile); //if user entered a file with a path then preserve it
            }

            //Groups must be checked later to make sure they are valid.
            //SharedUtilities has functions of check the validity, just make to so m->setGroups() after the checks.
            //If you are using these with a shared file no need to check the SharedRAbundVector class will call SharedUtilites for you,
            //kinda nice, huh?
            string groups = validParameter.valid(parameters, "groups");
            if (groups == "not found") { groups = ""; }
            else {  util.splitAtDash(groups, Groups); if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } } }

            //Commonly used to process list, rabund, sabund, shared and relabund files.
            //Look at "smart distancing" examples below in the execute function.
            string label = validParameter.valid(parameters, "label");
            if (label == "not found") {
                label = "";
            }
            else {
                if (label != "all") {
                    util.splitAtDash(label, labels);
                    allLines = 0;
                }
                else {
                    allLines = 1;
                }
            }

            string modeOption = validParameter.valid(parameters, "mode");
            if ( modeOption == "not found" || modeOption == "rfe" ) {
                mode = "rfe";
            }
            else if ( modeOption == "classify" ) {
                mode = "classify";
            }
            else {
                m->mothurOut("the mode option " + modeOption + " is not recognized -- must be 'rfe' or 'classify'");
                m->mothurOutEndLine();
                abort = true;
            }

            string ef = validParameter.valid(parameters, "evaluationfolds");
            if ( ef == "not found") {
                evaluationFoldCount = 3;
            }
            else {
                util.mothurConvert(ef, evaluationFoldCount);
            }
            string tf = validParameter.valid(parameters, "trainingfolds");
            if ( tf == "not found") {
                trainingFoldCount = 5;
            }
            else {
                util.mothurConvert(tf, trainingFoldCount);
            }

            string smocOption = validParameter.valid(parameters, "smoc");
            smocList.clear();
            if ( smocOption == "not found" ) {
                //smocOption = "0.001,0.01,0.1,1.0,10.0,100.0,1000.0";
            }
            else {
                vector<string> smocOptionList;
                //split(smocOption, ';', smocOptionList);
                util.splitAtDash(smocOption, smocOptionList);
                for (vector<string>::iterator i = smocOptionList.begin(); i != smocOptionList.end(); i++) {
                    smocList.push_back(atof(i->c_str()));
                }
                
            }

            // kernel specification
            // start with default parameter ranges for all kernels
            kernelParameterRangeMap.clear();
            getDefaultKernelParameterRangeMap(kernelParameterRangeMap);
            // get the kernel option
            string kernelOption = validParameter.valid(parameters, "kernel");
            // if the kernel option is "not found" then use all kernels with default parameter ranges
            // otherwise use only kernels listed in the kernelOption string
            if ( kernelOption == "not found" ) {

            }
            else {
                // if the kernel option has been specified then
                // remove kernel parameters from the kernel parameter map if
                // they are not listed in the kernel option
                // at this point the kernelParameterRangeMap looks like this:
                //    linear_key     : [
                //        smoc_key                : smoc parameter range
                //        linear_constant_key     : linear constant range
                //    ]
                //    rbf_key        : [
                //        smoc_key                : smoc parameter range
                //        rbf_gamma_key           : rbf gamma range
                //    ]
                //    polynomial_key : [
                //        smoc_key                : smoc parameter range
                //        polynomial_degree_key   : polynomial degree range
                //        polynomial_constant_key : polynomial constant range
                //    ]
                vector<string> kernelList;
                vector<string> unspecifiedKernelList;
                //split(kernelOption, '-', kernelList);
                util.splitAtDash(kernelOption, kernelList);
                set<string> kernelSet(kernelList.begin(), kernelList.end());
                // make a list of strings that are keys in the kernel parameter range map
                // but are not in the kernel list
                for (KernelParameterRangeMap::iterator i = kernelParameterRangeMap.begin(); i != kernelParameterRangeMap.end(); i++) {
                    //should be kernelList here
                    string kernelKey = i->first;
                    if ( kernelSet.find(kernelKey) == kernelSet.end() ) {
                        unspecifiedKernelList.push_back(kernelKey);
                    }
                }
                for (vector<string>::iterator i = unspecifiedKernelList.begin(); i != unspecifiedKernelList.end(); i++) {
                    m->mothurOut("removing kernel " + *i ); m->mothurOutEndLine();
                    kernelParameterRangeMap.erase(*i);
                }
            }

            // go through the kernel parameter range map and check for options for each kernel
            for (KernelParameterRangeMap::iterator i = kernelParameterRangeMap.begin(); i != kernelParameterRangeMap.end(); i++) {
                string kernelKey = i->first;
                ParameterRangeMap& kernelParameters = i->second;
                for (ParameterRangeMap::iterator j = kernelParameters.begin(); j != kernelParameters.end(); j++) {
                    string parameterKey = j->first;
                    ParameterRange& kernelParameterRange = j->second;
                    // has an option for this kernel parameter been specified?
                    string kernelParameterKey = kernelKey + parameterKey;
                    //m->mothurOut("looking for option " << kernelParameterKey << endl;
                    string kernelParameterOption = validParameter.valid(parameters, kernelParameterKey);
                    if (kernelParameterOption == "not found") {
                        // we already have default values in the kernel parameter map
                    }
                    else {
                        // replace the default parameters with the specified parameters
                        kernelParameterRange.clear();
                        vector<string> parameterList;
                        //split(kernelParameterOption, ';', parameterList);
                        util.splitAtDash(kernelParameterOption, parameterList);
                        for (vector<string>::iterator k = parameterList.begin(); k != parameterList.end(); k++) {
                            kernelParameterRange.push_back(atof(k->c_str()));
                        }
                    }
                }
            }

            // get the normalization option
            string transformOption = validParameter.valid(parameters, "transform");
            if ( transformOption == "not found" || transformOption == "unitmean") {
                transformName = "unitmean";
            }
            else if ( transformOption == "zeroone" ) {
                transformName = "zeroone";
            }
            else {
                m->mothurOut("the transform option " + transformOption + " is not recognized -- must be 'unitmean' or 'zeroone'");
                m->mothurOutEndLine();
                abort = true;
            }

            // get the verbosity option
            string verbosityOption = validParameter.valid(parameters, "verbose");
            if ( verbosityOption == "not found") {
                verbosity = 0;
            }
            else {
                util.mothurConvert(tf, verbosity);
                if (verbosity < OutputFilter::QUIET || verbosity > OutputFilter::TRACE) {
                    m->mothurOut("verbose set to unsupported value " + verbosityOption  + " -- must be between 0 and 3");
                }
            }

            // get the std threshold option
            string stdthresholdOption = validParameter.valid(parameters, "stdthreshold");
            if ( stdthresholdOption == "not found" ) {
                stdthreshold = -1.0;
            }
            else {
                util.mothurConvert(stdthresholdOption, stdthreshold);
                if ( stdthreshold <= 0.0 ) {
                    m->mothurOut("stdthreshold set to unsupported value " + stdthresholdOption  + " -- must be greater than 0.0");
                }
            }
        }

    }
    catch (exception& e) {
        m->errorOut(e, "ClassifySvmSharedCommand", "ClassifySvmSharedCommand");
        exit(1);
    }
}

//**********************************************************************************************************************
int ClassifySvmSharedCommand::execute() {
  try {

    if (abort) { if (calledHelp) { return 0; }  return 2;   }

    InputData input(sharedfile, "sharedfile", Groups);
    SharedRAbundVectors* lookup = input.getSharedRAbundVectors();
    vector<string> currentLabels = lookup->getOTUNames();
      Groups = lookup->getNamesGroups();
    //read design file
    designMap.read(designfile);

    string lastLabel = lookup->getLabel();
    set<string> processedLabels;
    set<string> userLabels = labels;

      //as long as you are not at the end of the file or done wih the lines you want
    while((lookup != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {

        if (m->getControl_pressed()) { delete lookup;  return 0; }

      if(allLines == 1 || labels.count(lookup->getLabel()) == 1){

        m->mothurOut(lookup->getLabel()+"\n"); 
          vector<SharedRAbundVector*> data = lookup->getSharedRAbundVectors();
        processSharedAndDesignData(data, currentLabels);
          for (int i = 0; i < data.size(); i++) { delete data[i]; } data.clear();

        processedLabels.insert(lookup->getLabel());
        userLabels.erase(lookup->getLabel());
      }

      if ((util.anyLabelsToProcess(lookup->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
        string saveLabel = lookup->getLabel();

        delete lookup;
        lookup = input.getSharedRAbundVectors(lastLabel);
        m->mothurOut(lookup->getLabel()+"\n"); 
          vector<SharedRAbundVector*> data = lookup->getSharedRAbundVectors();
          processSharedAndDesignData(data, currentLabels);
          for (int i = 0; i < data.size(); i++) { delete data[i]; } data.clear();

        processedLabels.insert(lookup->getLabel());
        userLabels.erase(lookup->getLabel());

          //restore real lastlabel to save below
        lookup->setLabels(saveLabel);
      }

      lastLabel = lookup->getLabel();
        //prevent memory leak
      delete lookup;

      if (m->getControl_pressed()) { return 0; }

        //get next line to process
      lookup = input.getSharedRAbundVectors();
    }

    if (m->getControl_pressed()) {  return 0; }

      //output error messages about any remaining user labels
    set<string>::iterator it;
    bool needToRun = false;
    for (it = userLabels.begin(); it != userLabels.end(); it++) {
      m->mothurOut("Your file does not include the label " + *it);
      if (processedLabels.count(lastLabel) != 1) {
        m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
        needToRun = true;
      }else {
        m->mothurOut(". Please refer to " + lastLabel + "."); m->mothurOutEndLine();
      }
    }

      //run last label if you need to
    if (needToRun )  {
      delete lookup;
      lookup = input.getSharedRAbundVectors(lastLabel);

      m->mothurOut(lookup->getLabel()+"\n"); 

        vector<SharedRAbundVector*> data = lookup->getSharedRAbundVectors();
        processSharedAndDesignData(data, currentLabels);
        for (int i = 0; i < data.size(); i++) { delete data[i]; } data.clear();

      delete lookup;

    }

      m->mothurOutEndLine();
      m->mothurOut("Output File Names: "); m->mothurOutEndLine();
      for (int i = 0; i < outputNames.size(); i++) {    m->mothurOut(outputNames[i]); m->mothurOutEndLine();    }
      m->mothurOutEndLine();

    return 0;

  }
  catch(exception& e) {
    m->errorOut(e, "ClassifySharedCommand", "execute");
    exit(1);
  }
}
//**********************************************************************************************************************
// This static function is intended to read all the necessary information from
// a pair of shared and design files needed for SVM classification.  This information
// is used to build a LabeledObservationVector.  Each element of the LabeledObservationVector
// looks like this:
//   LabeledObservationVector[0] = pair("label 0", &vector[10.0, 21.0, 13.0])
// where the vector in the second position of the pair records OTU abundances.
void ClassifySvmSharedCommand::readSharedAndDesignFiles(const string& sharedFilePath, const string& designFilePath, LabeledObservationVector& labeledObservationVector, FeatureVector& featureVector) {
    InputData input(sharedFilePath, "sharedfile", Groups);
    SharedRAbundVectors* lookup = input.getSharedRAbundVectors();
    Groups = lookup->getNamesGroups();

    DesignMap designMap;
    designMap.read(designFilePath);

    while ( lookup != NULL ) {
        vector<SharedRAbundVector*> data = lookup->getSharedRAbundVectors();
        readSharedRAbundVectors(data, designMap, labeledObservationVector, featureVector, lookup->getOTUNames());
        for (int i = 0; i < data.size(); i++) { delete data[i]; } data.clear();
        delete lookup;
        lookup = input.getSharedRAbundVectors();
    }
}

void ClassifySvmSharedCommand::readSharedRAbundVectors(vector<SharedRAbundVector*>& lookup, DesignMap& designMap, LabeledObservationVector& labeledObservationVector, FeatureVector& featureVector, vector<string> currentLabels) {
    for ( int j = 0; j < lookup.size(); j++ ) {
        //i++;
        //vector<individual> data = lookup[j]->getData();
        Observation* observation = new Observation(lookup[j]->getNumBins(), 0.0);
        string sharedGroupName = lookup[j]->getGroup();
        string treatmentName = designMap.get(sharedGroupName);
        
        //labeledObservationVector.push_back(make_pair(treatmentName, observation));
        labeledObservationVector.push_back(LabeledObservation(j, treatmentName, observation));
        
        for (int k = 0; k < lookup[j]->size(); k++) {
            
            observation->at(k) = double(lookup[j]->get(k));
            if ( j == 0) {
                featureVector.push_back(Feature(k, currentLabels[k]));
            }
        }
        
        // let this happen later?
        //delete lookup[j];
    }
}

void printPerformanceSummary(MultiClassSVM* s, ostream& output) {
    output << "multiclass SVM accuracy: " << s->getAccuracy() << endl;

    output << "two-class SVM performance" << endl;
    int labelFieldWidth = 2 + max_element(s->getLabels().begin(), s->getLabels().end())->size();
    int performanceFieldWidth = 10;
    int performancePrecision = 3;
    output << setw(labelFieldWidth) << "class 1"
              << setw(labelFieldWidth) << "class 2"
              << setw(performanceFieldWidth) << "precision"
              << setw(performanceFieldWidth) << "recall"
              << setw(performanceFieldWidth) << "f"
              << setw(performanceFieldWidth) << "accuracy" << endl;
    for ( SvmVector::const_iterator svm = s->getSvmList().begin(); svm != s->getSvmList().end(); svm++ ) {
        SvmPerformanceSummary sps = s->getSvmPerformanceSummary(**svm);
        output << setw(labelFieldWidth) << setprecision(performancePrecision) << sps.getPositiveClassLabel()
                  << setw(labelFieldWidth) << setprecision(performancePrecision) << sps.getNegativeClassLabel()
                  << setw(performanceFieldWidth) << setprecision(performancePrecision) << sps.getPrecision()
                  << setw(performanceFieldWidth) << setprecision(performancePrecision) << sps.getRecall()
                  << setw(performanceFieldWidth) << setprecision(performancePrecision) << sps.getF()
                  << setw(performanceFieldWidth) << setprecision(performancePrecision) << sps.getAccuracy() << endl;
    }

}
//**********************************************************************************************************************

void ClassifySvmSharedCommand::processSharedAndDesignData(vector<SharedRAbundVector*> lookup, vector<string> currentLabels) {
    try {
        OutputFilter outputFilter(verbosity);

        LabeledObservationVector labeledObservationVector;
        FeatureVector featureVector;
        readSharedRAbundVectors(lookup, designMap, labeledObservationVector, featureVector, currentLabels);

        // optionally remove features with low standard deviation
        if ( stdthreshold > 0.0 ) {
            FeatureVector removedFeatureVector = applyStdThreshold(stdthreshold, labeledObservationVector, featureVector);
            if (removedFeatureVector.size() > 0) {
                m->mothurOut(toString(removedFeatureVector.size()) + " OTUs were below the stdthreshold of " + toString(stdthreshold) + " and were removed"); m->mothurOutEndLine();
                if ( outputFilter.debug() ) {
                    m->mothurOut("the following OTUs were below the standard deviation threshold of " + toString(stdthreshold) ); m->mothurOutEndLine();
                    for (FeatureVector::iterator i = removedFeatureVector.begin(); i != removedFeatureVector.end(); i++) {
                        m->mothurOut("  " + toString(i->getFeatureLabel()) ); m->mothurOutEndLine();
                    }
                }
            }
        }

        // apply [0,1] standardization
        if ( transformName == "zeroone") {
            m->mothurOut("transforming data to lie within range [0,1]"); m->mothurOutEndLine();
            transformZeroOne(labeledObservationVector);
        }
        else {
            m->mothurOut("transforming data to have zero mean and unit variance"); m->mothurOutEndLine();
            transformZeroMeanUnitVariance(labeledObservationVector);
        }

        SvmDataset svmDataset(labeledObservationVector, featureVector);

        OneVsOneMultiClassSvmTrainer trainer(svmDataset, evaluationFoldCount, trainingFoldCount, outputFilter);

        if ( mode == "rfe" ) {
            SvmRfe svmRfe;
            ParameterRange& linearKernelConstantRange = kernelParameterRangeMap["linear"]["constant"];
            ParameterRange& linearKernelSmoCRange = kernelParameterRangeMap["linear"]["smoc"];
            RankedFeatureList rankedFeatureList = svmRfe.getOrderedFeatureList(svmDataset, trainer, linearKernelConstantRange, linearKernelSmoCRange);

            map<string, string> variables;
            variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(sharedfile));
            variables["[distance]"] = lookup[0]->getLabel();
            string filename = getOutputFileName("summary", variables);
            outputNames.push_back(filename);
            outputTypes["summary"].push_back(filename);
            m->mothurOutEndLine();

            ofstream outputFile(filename.c_str());

            int n = 0;
            int rfeRoundCount = rankedFeatureList.front().getRank();
            m->mothurOut("ordered features:" ); m->mothurOutEndLine();
            m->mothurOut("index\tOTU\trank"); m->mothurOutEndLine();
            outputFile << setw(5)  << "index"
                       << setw(12) << "OTU"
                       << setw(5)  << "rank"
                       << endl;
            for (RankedFeatureList::iterator i = rankedFeatureList.begin(); i != rankedFeatureList.end(); i++) {
                n++;
                int rank = rfeRoundCount - i->getRank() + 1;
                outputFile << setw(5)  << n
                           << setw(12) << i->getFeature().getFeatureLabel()
                           << setw(5)  << rank
                           ; m->mothurOutEndLine();
                if ( n <= 20 ) {
                    m->mothurOut(toString(n)
                              + toString(i->getFeature().getFeatureLabel())
                              + toString(rank)
                              ); m->mothurOutEndLine();
                }
            }
            outputFile.close();
        }
        else {
            MultiClassSVM* mcsvm = trainer.train(kernelParameterRangeMap);

            map<string, string> variables;
            variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(sharedfile));
            variables["[distance]"] = lookup[0]->getLabel();
            string filename = getOutputFileName("summary", variables);
            outputNames.push_back(filename);
            outputTypes["summary"].push_back(filename);
            m->mothurOutEndLine();

            ofstream outputFile(filename.c_str());

            printPerformanceSummary(mcsvm, cout);
            printPerformanceSummary(mcsvm, outputFile);

            outputFile << "actual  predicted" << endl;
            for ( LabeledObservationVector::const_iterator i = labeledObservationVector.begin(); i != labeledObservationVector.end(); i++ ) {
                Label actualLabel = i->getLabel();
                outputFile << i->getDatasetIndex() << " " << actualLabel << " ";
                try {
                    Label predictedLabel = mcsvm->classify(*(i->getObservation()));
                    outputFile << predictedLabel << endl;
                }
                catch ( MultiClassSvmClassificationTie& e ) {
                    outputFile << "tie" << endl;
                    m->mothurOut("classification tie for observation " + toString(i->datasetIndex) + " with label " + toString(i->first)); m->mothurOutEndLine();
                }

            }
            outputFile.close();
            delete mcsvm;
        }

    }
    catch (exception& e) {
        m->errorOut(e, "ClassifySvmSharedCommand", "processSharedAndDesignData");
        exit(1);
    }
}
//**********************************************************************************************************************

void ClassifySvmSharedCommand::trainSharedAndDesignData(vector<SharedRAbundVector*> lookup, vector<string> currentLabels) {
    try {
        LabeledObservationVector labeledObservationVector;
        FeatureVector featureVector;
        readSharedRAbundVectors(lookup, designMap, labeledObservationVector, featureVector, currentLabels);
        SvmDataset svmDataset(labeledObservationVector, featureVector);
        int evaluationFoldCount = 3;
        int trainFoldCount = 5;
        OutputFilter outputFilter(2);
        OneVsOneMultiClassSvmTrainer t(svmDataset, evaluationFoldCount, trainFoldCount, outputFilter);
        KernelParameterRangeMap kernelParameterRangeMap;
        getDefaultKernelParameterRangeMap(kernelParameterRangeMap);
        t.train(kernelParameterRangeMap);

        m->mothurOut("done training" ); m->mothurOutEndLine();

        map<string, string> variables;
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(sharedfile));
        variables["[distance]"] = lookup[0]->getLabel();
        string filename = getOutputFileName("summary", variables);
        outputNames.push_back(filename);
        outputTypes["summary"].push_back(filename);

        m->mothurOutEndLine();
        m->mothurOut("leaving processSharedAndDesignData" ); m->mothurOutEndLine();
    }
    catch (exception& e) {
        m->errorOut(e, "ClassifySvmSharedCommand", "trainSharedAndDesignData");
        exit(1);
    }
}

//**********************************************************************************************************************
