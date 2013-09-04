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

        // Support Vector Machine parameters
        CommandParameter linearKernel("linear", "Boolean", "", "T", "", "", "", "", false, false);
        parameters.push_back(linearKernel);
        //CommandParameter potupersplit("otupersplit", "Multiple", "log2-squareroot", "log2", "", "", "","",false,false); parameters.push_back(potupersplit);
        //CommandParameter psplitcriteria("splitcriteria", "Multiple", "gainratio-infogain", "gainratio", "", "", "","",false,false); parameters.push_back(psplitcriteria);
        //CommandParameter pnumtrees("numtrees", "Number", "", "100", "", "", "","",false,false); parameters.push_back(pnumtrees);

        // parameters related to pruning
        //CommandParameter pdopruning("prune", "Boolean", "", "T", "", "", "", "", false, false); parameters.push_back(pdopruning);
        //CommandParameter ppruneaggrns("pruneaggressiveness", "Number", "", "0.9", "", "", "", "", false, false); parameters.push_back(ppruneaggrns);
        //CommandParameter pdiscardhetrees("discarderrortrees", "Boolean", "", "T", "", "", "", "", false, false); parameters.push_back(pdiscardhetrees);
        //CommandParameter phetdiscardthreshold("errorthreshold", "Number", "", "0.4", "", "", "", "", false, false); parameters.push_back(phetdiscardthreshold);
        //CommandParameter psdthreshold("stdthreshold", "Number", "", "0.0", "", "", "", "", false, false); parameters.push_back(psdthreshold);
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
            m->control_pressed = true;
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
                if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {
                    abort = true;
                }
            }
            vector<string> tempOutNames;
            outputTypes["summary"] = tempOutNames;

            //if the user changes the input directory command factory will send this info to us in the output parameter
            string inputDir = validParameter.validFile(parameters, "inputdir", false);
            if (inputDir == "not found") {
                inputDir = "";
            }
            else {
                string path;
                it = parameters.find("shared");
                //user has given a shared file
                if (it != parameters.end()) {
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {
                        parameters["shared"] = inputDir + it->second;
                    }
                }

                it = parameters.find("design");
                //user has given a design file
                if (it != parameters.end()) {
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {
                        parameters["design"] = inputDir + it->second;
                    }
                }

            }
            //check for parameters
            //get shared file, it is required
            sharedfile = validParameter.validFile(parameters, "shared", true);
            if (sharedfile == "not open") {
                sharedfile = "";
                abort = true;
            }
            else if (sharedfile == "not found") {
                //if there is a current shared file, use it
                sharedfile = m->getSharedFile();
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
                m->setSharedFile(sharedfile);
            }

            //get design file, it is required
            designfile = validParameter.validFile(parameters, "design", true);
            if (designfile == "not open") {
                sharedfile = "";
                abort = true;
            }
            else if (designfile == "not found") {
                //if there is a current shared file, use it
                designfile = m->getDesignFile();
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
                m->setDesignFile(designfile);
            }

            //if the user changes the output directory command factory will send this info to us in the output parameter
            outputDir = validParameter.validFile(parameters, "outputdir", false);
            if (outputDir == "not found") {
                outputDir = m->hasPath(sharedfile); //if user entered a file with a path then preserve it
            }

            //Groups must be checked later to make sure they are valid.
            //SharedUtilities has functions of check the validity, just make to so m->setGroups() after the checks.
            //If you are using these with a shared file no need to check the SharedRAbundVector class will call SharedUtilites for you,
            //kinda nice, huh?
            string groups = validParameter.validFile(parameters, "groups", false);
            if (groups == "not found") {
                groups = "";
            }
            else {
                m->splitAtDash(groups, Groups);
            }
            m->setGroups(Groups);

            //Commonly used to process list, rabund, sabund, shared and relabund files.
            //Look at "smart distancing" examples below in the execute function.
            string label = validParameter.validFile(parameters, "label", false);
            if (label == "not found") {
                label = "";
            }
            else {
                if (label != "all") {
                    m->splitAtDash(label, labels);
                    allLines = 0;
                }
                else {
                    allLines = 1;
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
// This static function is intended to read all the necessary information from
// a pair of shared and design files needed for SVM classification.  This information
// is used to build a LabeledObservationVector.  Each element of the LabeledObservationVector
// looks like this:
//   LabeledObservationVector[0] = pair("label 0", &vector[10.0, 21.0, 13.0])
// where the vector in the second position of the pair records OTU abundances.
void ClassifySvmSharedCommand::readSharedAndDesignFiles(const std::string& sharedFilePath, const std::string& designFilePath, LabeledObservationVector& labeledObservationVector, FeatureVector& featureVector) {
    InputData input(sharedFilePath, "sharedfile");
    vector<SharedRAbundVector*> lookup = input.getSharedRAbundVectors();

    GroupMap designMap;
    designMap.readDesignMap(designFilePath);

    while ( lookup[0] != NULL ) {
        readSharedRAbundVectors(lookup, designMap, labeledObservationVector, featureVector);
        lookup = input.getSharedRAbundVectors();
    }
    /*
    int i = 0;
    while ( lookup[0] != NULL ) {
        for ( int j = 0; j < lookup.size(); j++ ) {
            i++;
            vector<individual> data = lookup[j]->getData();
            Observation* observation = new Observation(data.size(), 0.0);
            string sharedGroupName = lookup[j]->getGroup();
            string treatmentName = designMap.getGroup(sharedGroupName);
            //std::cout << "shared group name: " << sharedGroupName << " treatment name: " << treatmentName << std::endl;
            labeledObservationVector.push_back(std::make_pair(treatmentName, observation));
            //std::cout << "i=" << i << " j=" << j << " label : " << lookup[j]->getLabel() << " group: " << lookup[j]->getGroup();
            for (int k = 0; k < data.size(); k++) {
                //std::cout << " abundance " << data[k].abundance;
                observation->at(k) = double(data[k].abundance);
            }
            //std::cout << std::endl;
            delete lookup[j];
        }
        //lookup[0] = NULL;
        lookup = input.getSharedRAbundVectors();
    }
    */
}

void ClassifySvmSharedCommand::readSharedRAbundVectors(vector<SharedRAbundVector*>& lookup, GroupMap& designMap, LabeledObservationVector& labeledObservationVector, FeatureVector& featureVector) {
    //std::cout << " lookup.size() = " << lookup.size() << std::endl;
    //int i = 0;
    //while ( lookup[0] != NULL ) {
        for ( int j = 0; j < lookup.size(); j++ ) {
            //i++;
            vector<individual> data = lookup[j]->getData();
            Observation* observation = new Observation(data.size(), 0.0);
            string sharedGroupName = lookup[j]->getGroup();
            string treatmentName = designMap.getGroup(sharedGroupName);
            //std::cout << "shared group name: " << sharedGroupName << " treatment name: " << treatmentName << std::endl;
            //labeledObservationVector.push_back(std::make_pair(treatmentName, observation));
            labeledObservationVector.push_back(LabeledObservation(j, treatmentName, observation));
            //std::cout << " j=" << j << " label : " << lookup[j]->getLabel() << " group: " << lookup[j]->getGroup();
            for (int k = 0; k < data.size(); k++) {
                //std::cout << " abundance " << data[k].abundance;
                observation->at(k) = double(data[k].abundance);
                if ( j == 0) {
                    featureVector.push_back(Feature(k, m->currentBinLabels[k]));
                }
            }
            //std::cout << std::endl;
            delete lookup[j];
        }
        //lookup[0] = NULL;
        //lookup = input.getSharedRAbundVectors();
    //}
}

//**********************************************************************************************************************
int ClassifySvmSharedCommand::execute() {
    try {

        if (abort == true) {
            if (calledHelp) {
                return 0;
            }
            return 2;
        }

        InputData input(sharedfile, "sharedfile");
        vector<SharedRAbundVector*> lookup = input.getSharedRAbundVectors();

        //read design file
        designMap.readDesignMap(designfile);

        string lastLabel = lookup[0]->getLabel();
        set<string> processedLabels;
        set<string> userLabels = labels;

        //as long as you are not at the end of the file or done with the lines you want
        while ((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {

            if (m->control_pressed) {
                for (int i = 0; i < lookup.size(); i++) {
                    delete lookup[i];
                }
                return 0;
            }

            if (allLines == 1 || labels.count(lookup[0]->getLabel()) == 1) {

                m->mothurOut(lookup[0]->getLabel());
                m->mothurOutEndLine();

                processSharedAndDesignData(lookup);

                processedLabels.insert(lookup[0]->getLabel());
                userLabels.erase(lookup[0]->getLabel());
            }

            if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true)
                    && (processedLabels.count(lastLabel) != 1)) {
                string saveLabel = lookup[0]->getLabel();

                for (int i = 0; i < lookup.size(); i++) {
                    delete lookup[i];
                }
                lookup = input.getSharedRAbundVectors(lastLabel);
                m->mothurOut(lookup[0]->getLabel());
                m->mothurOutEndLine();
                processSharedAndDesignData(lookup);

                processedLabels.insert(lookup[0]->getLabel());
                userLabels.erase(lookup[0]->getLabel());

                //restore real lastlabel to save below
                lookup[0]->setLabel(saveLabel);
            }

            lastLabel = lookup[0]->getLabel();
            //prevent memory leak
            for (int i = 0; i < lookup.size(); i++) {
                delete lookup[i];
                lookup[i] = NULL;
            }

            if (m->control_pressed) {
                return 0;
            }

            //get next line to process
            lookup = input.getSharedRAbundVectors();
        }

        if (m->control_pressed) {
            return 0;
        }

        //output error messages about any remaining user labels
        set<string>::iterator it;
        bool needToRun = false;
        for (it = userLabels.begin(); it != userLabels.end(); it++) {
            m->mothurOut("Your file does not include the label " + *it);
            if (processedLabels.count(lastLabel) != 1) {
                m->mothurOut(". I will use " + lastLabel + ".");
                m->mothurOutEndLine();
                needToRun = true;
            }
            else {
                m->mothurOut(". Please refer to " + lastLabel + ".");
                m->mothurOutEndLine();
            }
        }

        //run last label if you need to
        if (needToRun == true) {
            for (int i = 0; i < lookup.size(); i++) {
                if (lookup[i] != NULL) {
                    delete lookup[i];
                }
            }
            lookup = input.getSharedRAbundVectors(lastLabel);

            m->mothurOut(lookup[0]->getLabel());
            m->mothurOutEndLine();

            processSharedAndDesignData(lookup);

            for (int i = 0; i < lookup.size(); i++) {
                delete lookup[i];
            }

        }

        m->mothurOutEndLine();
        m->mothurOut("Output File Names: ");
        m->mothurOutEndLine();
        for (int i = 0; i < outputNames.size(); i++) {
            m->mothurOut(outputNames[i]);
            m->mothurOutEndLine();
        }
        m->mothurOutEndLine();

        return 0;

    }
    catch (exception& e) {
        m->errorOut(e, "ClassifySvmSharedCommand", "execute");
        exit(1);
    }
}
//**********************************************************************************************************************

void ClassifySvmSharedCommand::processSharedAndDesignData(vector<SharedRAbundVector*> lookup) {
    try {
        LabeledObservationVector labeledObservationVector;
        FeatureVector featureVector;
        readSharedRAbundVectors(lookup, designMap, labeledObservationVector, featureVector);
        SvmDataset svmDataset(labeledObservationVector, featureVector);
        int evaluationFoldCount = 3;
        int trainFoldCount = 5;
        OneVsOneMultiClassSvmTrainer trainer(svmDataset, evaluationFoldCount, trainFoldCount, *this);
        //KernelParameterRangeMap kernelParameterRangeMap;
        //getDefaultKernelParameterRangeMap(kernelParameterRangeMap);
        //trainer.train(kernelParameterRangeMap);

        SvmRfe svmRfe;
        FeatureList orderedFeatureList = svmRfe.getOrderedFeatureList(svmDataset, trainer, LinearKernelFunction::defaultConstantRange, SmoTrainer::defaultCRange);

        int n = 0;
        std::cout << "ordered features:" << std::endl;
        for (FeatureList::iterator i = orderedFeatureList.begin(); i != orderedFeatureList.end(); i++) {
            std::cout << i->getFeatureLabel() << std::endl;
            n++;
            if (n > 20) break;
        }
        std::cout << "done training" << std::endl;

        map<string, string> variables;
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(sharedfile));
        variables["[distance]"] = lookup[0]->getLabel();
        string filename = getOutputFileName("summary", variables);
        outputNames.push_back(filename);
        outputTypes["summary"].push_back(filename);

        m->mothurOutEndLine();
        std::cout << "leaving processSharedAndDesignData" << std::endl;
    }
    catch (exception& e) {
        m->errorOut(e, "ClassifySvmSharedCommand", "processSharedAndDesignData");
        exit(1);
    }
}
//**********************************************************************************************************************

void ClassifySvmSharedCommand::trainSharedAndDesignData(vector<SharedRAbundVector*> lookup) {
    try {
        LabeledObservationVector labeledObservationVector;
        FeatureVector featureVector;
        readSharedRAbundVectors(lookup, designMap, labeledObservationVector, featureVector);
        SvmDataset svmDataset(labeledObservationVector, featureVector);
        int evaluationFoldCount = 3;
        int trainFoldCount = 5;
        OneVsOneMultiClassSvmTrainer t(svmDataset, evaluationFoldCount, trainFoldCount, *this);
        KernelParameterRangeMap kernelParameterRangeMap;
        getDefaultKernelParameterRangeMap(kernelParameterRangeMap);
        t.train(kernelParameterRangeMap);

        std::cout << "done training" << std::endl;

        map<string, string> variables;
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(sharedfile));
        variables["[distance]"] = lookup[0]->getLabel();
        string filename = getOutputFileName("summary", variables);
        outputNames.push_back(filename);
        outputTypes["summary"].push_back(filename);

        m->mothurOutEndLine();
        std::cout << "leaving processSharedAndDesignData" << std::endl;
    }
    catch (exception& e) {
        m->errorOut(e, "ClassifySvmSharedCommand", "processSharedAndDesignData");
        exit(1);
    }
}
//**********************************************************************************************************************
