//
//  lefsecommand.cpp
//  Mothur
//
//  Created by SarahsWork on 6/12/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "lefsecommand.h"
#include "linearalgebra.h"


//**********************************************************************************************************************
vector<string> LefseCommand::setParameters(){
	try {
        CommandParameter pdesign("design", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(pdesign);
        CommandParameter pshared("shared", "InputTypes", "", "", "shared-clr", "none", "none","metastats",false,false,true); parameters.push_back(pshared);
        //CommandParameter pclr("clr", "InputTypes", "", "", "shared-clr", "none", "none","metastats",false,false,true); parameters.push_back(pclr);
        CommandParameter pclass("class", "String", "", "", "", "", "","",false,false); parameters.push_back(pclass);
        CommandParameter psubclass("subclass", "String", "", "", "", "", "","",false,false); parameters.push_back(psubclass);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		//CommandParameter pclasses("classes", "String", "", "", "", "", "","",false,false); parameters.push_back(pclasses);
        CommandParameter palpha("aalpha", "Number", "", "0.05", "", "", "","",false,false); parameters.push_back(palpha);
        CommandParameter pwalpha("walpha", "Number", "", "0.05", "", "", "","",false,false); parameters.push_back(pwalpha);
        CommandParameter psets("sets", "String", "", "", "", "", "","",false,false); parameters.push_back(psets);
        CommandParameter plda("lda", "Number", "", "2.0", "", "", "","",false,false); parameters.push_back(plda);
        CommandParameter pwilc("wilc", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pwilc);
        CommandParameter pnormmillion("norm", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pnormmillion);
        CommandParameter piters("iters", "Number", "", "30", "", "", "","",false,false); parameters.push_back(piters);
        //CommandParameter pwilcsamename("wilcsamename", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pwilcsamename);
        CommandParameter pcurv("curv", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pcurv);
        CommandParameter pfiters("fboots", "Number", "", "0.67", "", "", "","",false,false); parameters.push_back(pfiters);
        CommandParameter pstrict("strict", "Multiple", "0-1-2", "0", "", "", "","",false,false); parameters.push_back(pstrict);
        CommandParameter pminc("minc", "Number", "", "10", "", "", "","",false,false); parameters.push_back(pminc);
        CommandParameter pmulticlass_strat("multiclass", "Multiple", "onevone-onevall", "onevall", "", "", "","",false,false); parameters.push_back(pmulticlass_strat);
        //CommandParameter psubject("subject", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(psubject);


        //not used in their current code, but in parameters
        //CommandParameter pnlogs("nlogs", "Number", "", "3", "", "", "","",false,false); parameters.push_back(pnlogs);
        //CommandParameter pranktec("ranktec", "Multiple", "lda-svm", "lda", "", "", "","",false,false); parameters.push_back(pranktec); // svm not implemented in their source yet.
        //CommandParameter psvmnorm("svmnorm", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(psvmnorm); //not used because svm not implemented yet.

        
        //every command must have inputdir and outputdir.  This allows mothur users to redirect input and output files.
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false; allLines = true;
        
        vector<string> tempOutNames;
        outputTypes["summary"] = tempOutNames;
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "LefseCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string LefseCommand::getHelpString(){
	try {
		string helpString = "";
		helpString += "The lefse command allows you to ....\n";
		helpString += "The lefse command parameters are: shared, design, class, subclass, label, walpha, aalpha, lda, wilc, iters, curv, fboots, strict, minc, multiclass and norm.\n";
		helpString += "The class parameter is used to indicate the which category you would like used for the Kruskal Wallis analysis. If none is provided first category is used.\n";
        helpString += "The subclass parameter is used to indicate the .....If none is provided, second category is used, or if only one category subclass is ignored. \n";
        helpString += "The aalpha parameter is used to set the alpha value for the Krukal Wallis Anova test Default=0.05. \n";
        helpString += "The walpha parameter is used to set the alpha value for the Wilcoxon test. Default=0.05. \n";
        helpString += "The lda parameter is used to set the threshold on the absolute value of the logarithmic LDA score. Default=2.0. \n";
        helpString += "The wilc parameter is used to indicate whether to perform the Wilcoxon test. Default=T. \n";
        helpString += "The iters parameter is used to set the number of bootstrap iteration for LDA. Default=30. \n";
        //helpString += "The wilcsamename parameter is used to indicate whether perform the wilcoxon test only among the subclasses with the same name. Default=F. \n";
        helpString += "The sets parameter allows you to specify which of the sets in your designfile you would like to analyze. The set names are separated by dashes. THe default is all sets in the designfile.\n";
        helpString += "The curv parameter is used to set whether perform the wilcoxon testing the Curtis's approach [BETA VERSION] Default=F. \n";
        helpString += "The norm parameter is used to multiply relative abundances by 1000000. Recommended when very low values are present. Default=T. \n";
        helpString += "The fboots parameter is used to set the subsampling fraction value for each bootstrap iteration. Default=0.67. \n";
        helpString += "The strict parameter is used to set the multiple testing correction options. 0 no correction (more strict, default), 1 correction for independent comparisons, 2 correction for independent comparison. Options = 0,1,2. Default=0. \n";
        helpString += "The minc parameter is used to minimum number of samples per subclass for performing wilcoxon test. Default=10. \n";
        helpString += "The multiclass parameter is used to (for multiclass tasks) set whether the test is performed in a one-against-one ( onevone - more strict!) or in a one-against-all setting ( onevall - less strict). Default=onevall. \n";
        //helpString += "The classes parameter is used to indicate the classes you would like to use. Classes should be inputted in the following format: classes=label<value1|value2|value3>-label<value1|value2>. For example to include groups from treatment with value early or late and age= young or old.  class=treatment<Early|Late>-age<young|old>.\n";
        helpString += "The label parameter is used to indicate which distances in the shared file you would like to use. labels are separated by dashes.\n";
		helpString += "The lefse command should be in the following format: lefse(shared=final.an.shared, design=final.design, class=treatment, subclass=age).\n";
        return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "LefseCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string LefseCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "summary") {  pattern = "[filename],[distance],lefse_summary"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "LefseCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
LefseCommand::LefseCommand(string option)  {
	try {

		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
            sharedfile = validParameter.validFile(parameters, "shared");
            if (sharedfile == "not open") { abort = true; }
            else if (sharedfile == "not found") { sharedfile = "";  }
            else { current->setSharedFile(sharedfile); inputfile = sharedfile; format = "sharedfile";  }
            
            if ((sharedfile == "")) {
                //is there are current file available for any of these?
                //give priority to shared, then list, then rabund, then sabund
                //if there is a current shared file, use it
                sharedfile = current->getSharedFile();
                if (sharedfile != "") { inputfile = sharedfile; format = "sharedfile"; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n"); }
                else {
                    m->mothurOut("No valid current files. You must provide a shared file.\n"); abort = true;
                }
            }
            
            //get shared file, it is required
			designfile = validParameter.validFile(parameters, "design");
			if (designfile == "not open") { designfile = ""; abort = true; }
			else if (designfile == "not found") {
				//if there is a current shared file, use it
				designfile = current->getDesignFile();
				if (designfile != "") { m->mothurOut("Using " + designfile + " as input file for the design parameter.\n");  }
				else { 	m->mothurOut("You have no current design file and the design parameter is required.\n");  abort = true; }
			}else { current->setDesignFile(designfile); }
            
            
            if (outputdir == ""){
				outputdir = util.hasPath(inputfile); 
			}
            
            string label = validParameter.valid(parameters, "label");
			if (label == "not found") { label = ""; }
			else {
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				else { allLines = true;  }
			}
            
            mclass = validParameter.valid(parameters, "class");
			if (mclass == "not found") { mclass = ""; }
			
            subclass = validParameter.valid(parameters, "subclass");
			if (subclass == "not found") { subclass = mclass; }
            
            string temp = validParameter.valid(parameters, "aalpha");
			if (temp == "not found") { temp = "0.05"; }
			util.mothurConvert(temp, anovaAlpha);
            
            temp = validParameter.valid(parameters, "walpha");
			if (temp == "not found") { temp = "0.05"; }
			util.mothurConvert(temp, wilcoxonAlpha);
            
            temp = validParameter.valid(parameters, "wilc");
			if (temp == "not found") { temp = "T"; }
			wilc = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "norm");
			if (temp == "not found") { temp = "T"; }
			normMillion = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "lda");
			if (temp == "not found") { temp = "2.0"; }
			util.mothurConvert(temp, ldaThreshold);
            
            temp = validParameter.valid(parameters, "iters");
			if (temp == "not found") { temp = "30"; }
			util.mothurConvert(temp, iters);
            
            temp = validParameter.valid(parameters, "fboots");
			if (temp == "not found") { temp = "0.67"; }
			util.mothurConvert(temp, fBoots);
            
            temp = validParameter.valid(parameters, "curv");
			if (temp == "not found") { temp = "F"; }
			curv = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "strict");
            if (temp == "not found"){	temp = "0";		}
			if ((temp != "0") && (temp != "1") && (temp != "2")) { m->mothurOut("[ERROR]: Invalid strict option: choices are 0, 1 or 2.\n");  abort=true; }
            else {  util.mothurConvert(temp, strict); }
            
            temp = validParameter.valid(parameters, "minc");
			if (temp == "not found") { temp = "10"; }
			util.mothurConvert(temp, minC);
            
            sets = validParameter.valid(parameters, "sets");
            if (sets == "not found") { sets = ""; }
            else {
                util.splitAtDash(sets, Sets);
            }
            
            multiClassStrat = validParameter.valid(parameters, "multiclass");
            if (multiClassStrat == "not found"){	multiClassStrat = "onevall";		}
			if ((multiClassStrat != "onevall") && (multiClassStrat != "onevone")) { m->mothurOut("Invalid multiclass option: choices are onevone or onevall.\n");  abort=true; }
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "LefseCommand", "LefseCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int LefseCommand::execute(){
	try {
        unsigned int holdRandom = m->getRandomSeed();
        m->setRandomSeed(1982);
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        DesignMap designMap(designfile); if (m->getControl_pressed()) { return 0; }
        
        //if user did not select class use first column
        if (mclass == "") {  mclass = designMap.getDefaultClass(); m->mothurOut("\nYou did not provide a class, using " + mclass +".\n\n"); if (subclass == "") { subclass = mclass; } }
        
        vector<string> Groups;
        if (Sets.size() != 0) { //user has picked sets find groups to include from lookup
            designMap.setDefaultClass(mclass);
            Groups = designMap.getNamesGroups(Sets);
        }
        
        InputData input(inputfile, format, Groups);
        set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
        SharedRAbundFloatVectors* lookup = NULL; SharedCLRVectors* clr = NULL;
        
        if (format == "sharedfile") {
            lookup = util.getNextRelabund(input, allLines, userLabels, processedLabels, lastLabel);
            Groups = lookup->getNamesGroups();
        }
        
        while ((lookup != NULL) || (clr != NULL)){
            
            if (m->getControl_pressed()) { if (lookup != NULL) { delete lookup; } if (clr != NULL) { delete clr; }break; }
            
            process(lookup, clr, designMap);
            
            if (format == "sharedfile") { delete lookup; lookup = util.getNextRelabund(input, allLines, userLabels, processedLabels, lastLabel); }
        }
		
        //output files created by command
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
        m->setRandomSeed(holdRandom);
        return 0;
		
    }
	catch(exception& e) {
		m->errorOut(e, "LefseCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int LefseCommand::process(SharedRAbundFloatVectors*& lookup, SharedCLRVectors*& clr, DesignMap& designMap) {
	try {
        vector<string> classes;
        vector<string> subclasses;
        map<string, string> subclass2Class;
        map<string, set<string> > class2SubClasses; //maps class name to vector of its subclasses
        map<string, vector<int> > subClass2GroupIndex; //maps subclass name to vector of indexes in lookup from that subclass. old -> 1,2,3 means groups in location 1,2,3 of lookup are from old.  Saves time below.
        map<string, vector<int> > class2GroupIndex; //maps subclass name to vector of indexes in lookup from that class. old -> 1,2,3 means groups in location 1,2,3 of lookup are from old.  Saves time below.
        if (normMillion) {  normalize(lookup, clr);  }
        vector<string> namesOfGroups;
        if (lookup != NULL) { namesOfGroups = lookup->getNamesGroups();  }
        else                { namesOfGroups = clr->getNamesGroups();     }
        
        for (int j = 0; j < namesOfGroups.size(); j++) {
            string group = namesOfGroups[j];
            string treatment = designMap.get(group, mclass); //get value for this group in this category
            string thisSub = designMap.get(group, subclass);
            map<string, string>::iterator it = subclass2Class.find(thisSub);
            if (it == subclass2Class.end()) {
                subclass2Class[thisSub] = treatment;
                vector<int> temp; temp.push_back(j);
                subClass2GroupIndex[thisSub] = temp;
            }else {
                if (it->second != treatment) {
                    //m->mothurOut("[WARNING]: subclass " + thisSub + " has members in " + it->second + " and " + treatment + ". Subclass members must be from the same class for Wilcoxon. Changing " + thisSub + " to " + treatment + "_" + thisSub + ".\n");
                    thisSub = treatment + "_" + thisSub;
                    subclass2Class[thisSub] = treatment;
                    vector<int> temp; temp.push_back(j);
                    subClass2GroupIndex[thisSub] = temp;
                }else { subClass2GroupIndex[thisSub].push_back(j); }
            }
            
            map<string, set<string> >::iterator itClass = class2SubClasses.find(treatment);
            if (itClass == class2SubClasses.end()) {
                set<string> temp; temp.insert(thisSub);
                class2SubClasses[treatment] = temp;
                vector<int> temp2; temp2.push_back(j);
                class2GroupIndex[treatment] = temp2;
                classes.push_back(treatment);
            }else{
                class2SubClasses[treatment].insert(thisSub);
                class2GroupIndex[treatment].push_back(j);
            }
        }
        //sort classes so order is right
        sort(classes.begin(), classes.end());
        
        vector< vector<double> > means = getMeans(lookup, clr, class2GroupIndex); //[numOTUs][classes] - classes in same order as class2GroupIndex
        
        //run kruskal wallis on each otu
        map<int, double> significantOtuLabels = runKruskalWallis(lookup, clr, designMap);
        
        int numSigBeforeWilcox = significantOtuLabels.size();
        
        if (m->getDebug()) { m->mothurOut("[DEBUG]: completed Kruskal Wallis\n"); } 
        
        //check for subclass
        string wilcoxString = "";
        if ((subclass != "") && wilc) {  significantOtuLabels = runWilcoxon(lookup, clr, designMap, significantOtuLabels, class2SubClasses, subClass2GroupIndex, subclass2Class);  wilcoxString += " ( " + toString(numSigBeforeWilcox) + " ) before internal wilcoxon"; }
        
        int numSigAfterWilcox = significantOtuLabels.size();
        
        if (m->getDebug()) { m->mothurOut("[DEBUG]: completed Wilcoxon\n"); } 
        
        m->mothurOut("\nNumber of significantly discriminative features: " + toString(numSigAfterWilcox) + wilcoxString + ".\n"); 
        
        map<int, double> sigOTUSLDA;
        if (numSigAfterWilcox > 0) {
            sigOTUSLDA = testLDA(lookup, clr, significantOtuLabels, class2GroupIndex, subClass2GroupIndex);
            m->mothurOut("Number of discriminative features with abs LDA score > " + toString(ldaThreshold) + " : " + toString(sigOTUSLDA.size()) + ".\n");
        }
        else { m->mothurOut("No features with significant differences between the classes.\n"); }
        
        if (m->getDebug()) { m->mothurOut("[DEBUG]: completed lda\n"); } 
        
        string label;
        vector<string> otuNames;
        if (lookup != NULL) {   label = lookup->getLabel();  otuNames = lookup->getOTUNames(); }
        else                {   label = clr->getLabel();  otuNames = clr->getOTUNames();       }
        
        printResults(means, significantOtuLabels, sigOTUSLDA, label, classes, otuNames);
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "LefseCommand", "process");
		exit(1);
	}
}
//**********************************************************************************************************************
int LefseCommand::normalize(SharedRAbundFloatVectors*& lookup, SharedCLRVectors*& clr) {
	try {
        vector<double> mul;
        vector<string> namesOfGroups;
        int numBins = 0;
        int numSamples = 0;
        
        if (lookup != NULL) {
            namesOfGroups = lookup->getNamesGroups();
            numBins = lookup->getNumBins();
            numSamples = lookup->size();
        }else {
            namesOfGroups = clr->getNamesGroups();
            numBins = clr->getNumBins();
            numSamples = clr->size();
        }
        
        for (int i = 0; i < numSamples; i++) {
            double sum = 0;
            if (lookup != NULL) {  sum = lookup->getNumSeqs(namesOfGroups[i]);    }
            else                {  sum = clr->getNumSeqs(namesOfGroups[i]);       }
            mul.push_back(1000000.0/sum);
        }
        
        for (int i = 0; i < numSamples; i++) {
            for (int j = 0; j < numBins; j++) {
                if (lookup != NULL) { lookup->set(j, lookup->get(j, namesOfGroups[i])*mul[i], namesOfGroups[i]);    }
                else                { clr->set(j, clr->get(j, namesOfGroups[i])*mul[i], namesOfGroups[i]);       }
            }
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "LefseCommand", "normalize");
		exit(1);
	}
}
//**********************************************************************************************************************
map<int, double> LefseCommand::runKruskalWallis(SharedRAbundFloatVectors*& lookup, SharedCLRVectors*& clr, DesignMap& designMap) {
	try {
        vector<string> namesOfGroups;
        int numBins = 0;
        
        if (lookup != NULL) {
            namesOfGroups = lookup->getNamesGroups();
            numBins = lookup->getNumBins();
        }else {
            namesOfGroups = clr->getNamesGroups();
            numBins = clr->getNumBins();
        }
        map<int, double> significantOtuLabels;
        
        //sanity check to make sure each treatment has a group in the shared file
        set<string> treatments;
        for (int j = 0; j < namesOfGroups.size(); j++) {
            string group = namesOfGroups[j];
            string treatment = designMap.get(group, mclass); //get value for this group in this category
            treatments.insert(treatment);
        }
        if (treatments.size() < 2) { m->mothurOut("[ERROR]: need at least 2 things for classes to compare, quitting.\n"); m->setControl_pressed(true); }
        
        LinearAlgebra linear;
        for (int i = 0; i < numBins; i++) {
            if (m->getControl_pressed()) { break; }
            
            vector<spearmanRank> values;
            vector<float> abunds;
            if (lookup != NULL) { abunds = lookup->getOTU(i);   }
            else                { abunds = clr->getOTU(i);      }
            
            for (int j = 0; j < namesOfGroups.size(); j++) {
                string group = namesOfGroups[j];
                string treatment = designMap.get(group, mclass); //get value for this group in this category
                spearmanRank temp(treatment, abunds[j]);
                values.push_back(temp);
            }
            
            double pValue = 0.0;
            linear.calcKruskalWallis(values, pValue);
             
            if (pValue < anovaAlpha) {  significantOtuLabels[i] = pValue;  }
        }
        
        return significantOtuLabels;
    }
	catch(exception& e) {
		m->errorOut(e, "LefseCommand", "runKruskalWallis");
		exit(1);
	}
}
//**********************************************************************************************************************
//assumes not neccessarily paired
map<int, double> LefseCommand::runWilcoxon(SharedRAbundFloatVectors*& lookup, SharedCLRVectors*& clr, DesignMap& designMap, map<int, double> bins, map<string, set<string> >& class2SubClasses, map<string, vector<int> >& subClass2GroupIndex, map<string, string> subclass2Class) {
    try {
        map<int, double> significantOtuLabels;
        map<int, double>::iterator it;
        //if it exists and meets the following requirements run Wilcoxon
        /*
         1. Subclass members all belong to same main class
         anything else
        */
        int numBins = 0;
        if (lookup != NULL) { numBins = lookup->getNumBins();   }
        else                { numBins = clr->getNumBins();      }
        
        for (int i = 0; i < numBins; i++) {
            if (m->getControl_pressed()) { break; }
            
            it = bins.find(i);
            if (it != bins.end()) { //flagged in Kruskal Wallis
                
                vector<float> abunds;
                if (lookup != NULL) { abunds = lookup->getOTU(i);   }
                else                { abunds = clr->getOTU(i);      }
                
                bool sig = testOTUWilcoxon(class2SubClasses, abunds, subClass2GroupIndex, subclass2Class);
                if (sig) { significantOtuLabels[i] = it->second;  }
                
            }//bins flagged from kw
        }//for bins
        
        return significantOtuLabels;
    }
    catch(exception& e) {
        m->errorOut(e, "LefseCommand", "runWilcoxon");
        exit(1);
    }
}
//**********************************************************************************************************************
//lefse.py - test_rep_wilcoxon_r function
bool LefseCommand::testOTUWilcoxon(map<string, set<string> >& class2SubClasses, vector<float> abunds, map<string, vector<int> >& subClass2GroupIndex, map<string, string> subclass2Class) {
    try {
        int totalOk = 0;
        double alphaMtc = wilcoxonAlpha;
        vector< set<string> > allDiffs;
        LinearAlgebra linear;
        
        //for each subclass comparision
        map<string, set<string> >::iterator itB;
        for(map<string, set<string> >::iterator it=class2SubClasses.begin();it!=class2SubClasses.end();it++){
            itB = it;itB++;
            for(;itB!=class2SubClasses.end();itB++){
                if (m->getControl_pressed()) { return false; }
                bool first = true;
                int dirCmp = 0; // not set?? dir_cmp = "not_set" # 0=notset or none, 1=true, 2=false.
                int curv_sign = 0;
                int ok = 0;
                int count = 0;
                for (set<string>::iterator itClass1 = (it->second).begin(); itClass1 != (it->second).end(); itClass1++) {
                    bool br = false;
                    for (set<string>::iterator itClass2 = (itB->second).begin(); itClass2 != (itB->second).end(); itClass2++) {
                        string subclass1 = *itClass1;
                        string subclass2 = *itClass2;
                        count++;
                        
                        if (m->getDebug()) { m->mothurOut( "[DEBUG comparing " + it->first + "-" + *itClass1 + " to " + itB->first + "-" + *itClass2 + "\n"); }
                        
                        string treatment1 = subclass2Class[subclass1];
                        string treatment2 = subclass2Class[subclass2];
                        int numSubs1 = class2SubClasses[treatment1].size();
                        int numSubs2 = class2SubClasses[treatment2].size();
                        
                        //if mul_cor != 0: alpha_mtc = th*l_subcl1*l_subcl2 if mul_cor == 2 else 1.0-math.pow(1.0-th,l_subcl1*l_subcl2)
                        if (strict != 0) { alphaMtc = wilcoxonAlpha * numSubs1 * numSubs2 ; }
                        if (strict == 2) {}else{ alphaMtc = 1.0-pow((1.0-wilcoxonAlpha),(double)(numSubs1 * numSubs2)); }
                        
                        //fill x and y with this comparisons data
                        vector<double> x; vector<double> y;
                        
                        //fill x and y
                        vector<int> xIndexes = subClass2GroupIndex[subclass1]; //indexes in lookup for this subclass
                        vector<int> yIndexes = subClass2GroupIndex[subclass2]; //indexes in lookup for this subclass
                        for (int k = 0; k < yIndexes.size(); k++) { y.push_back(abunds[yIndexes[k]]);  }
                        for (int k = 0; k < xIndexes.size(); k++) { x.push_back(abunds[xIndexes[k]]);  }
                        
                        // med_comp = False
                        //if len(cl1) < min_c or len(cl2) < min_c:
                        //med_comp = True
                        bool medComp = false; // are there enough samples per subclass
                        if ((xIndexes.size() < minC) || (yIndexes.size() < minC)) { medComp = true; }
                        
                        double sx = util.median(x);
                        double sy = util.median(y);
                       
                        //if cl1[0] == cl2[0] and len(set(cl1)) == 1 and  len(set(cl2)) == 1:
                        //tres, first = False, False
                        double pValue = 0.0;
                        double H = 0.0;
                        bool tres = true; //don't think this is set in the python source.  Not sure how that is handled, but setting it here.
                        if (util.isEqual(x[0],y[0]) && (x.size() == 1) && (y.size() == 1)) { tres = false; first = false; }
                        else if (!medComp) {
                            H = linear.calcWilcoxon(x, y, pValue);
                            if (pValue < (alphaMtc*2.0)) { tres = true; }
                            else { tres = false; }
                        }
                        /*if first:
                         first = False
                         if not curv and ( med_comp or tres ):
                         dir_cmp = sx < sy
                         if sx == sy: br = True
                         elif curv:
                         dir_cmp = None
                         if med_comp or tres:
                         curv_sign += 1
                         dir_cmp = sx < sy
                         else: br = True
                         elif not curv and med_comp:
                         if ((sx < sy) != dir_cmp or sx == sy): br = True
                         elif curv:
                         if tres and dir_cmp == None:
                         curv_sign += 1
                         dir_cmp = sx < sy
                         if tres and dir_cmp != (sx < sy):
                         br = True
                         curv_sign = -1
                         elif not tres or (sx < sy) != dir_cmp or sx == sy: br = True
                         */
                        int sxSy = 2; //false
                        if (sx<sy) {  sxSy = 1; } //true
                        
                        if (first) {
                            first = false;
                            if ((!curv) && (medComp || tres)) {
                                dirCmp = 2; if (sx<sy) { dirCmp = 1; } //dir_cmp = sx < sy
                                if (util.isEqual(sx,sy)) { br = true; }
                            }else if (curv) {
                                dirCmp = 0;
                                if (medComp || tres) {
                                    curv_sign++;
                                    dirCmp = 2; if (sx<sy) { dirCmp = 1; } //dir_cmp = sx < sy
                                }
                            }else { br = true; }
                        }else if (!curv && medComp) {
                            if (sxSy != dirCmp || util.isEqual(sx,sy)) { br = true; }
                        }else if (curv) {
                            if (tres && dirCmp == 0) { curv_sign++; }
                            dirCmp = 2; if (sx<sy) { dirCmp = 1; } //dir_cmp = sx < sy
                            if (tres && dirCmp != sxSy) { //if tres and dir_cmp != (sx < sy):
                                br = true;
                                curv_sign = -1;
                            }
                        }else if (!tres || sxSy != dirCmp || util.isEqual(sx,sy)) { br = true; } //elif not tres or (sx < sy) != dir_cmp or sx == sy: br = True
                        if (br) { break; }
                        ok++;
                    }//for class2 subclasses
                    if (br) { break; }
                }//for class1 subclasses
                bool diff = false;
                if (curv) { diff = false; if (curv_sign > 0) { diff = true; } } //if curv: diff = curv_sign > 0
                else { //else: diff = (ok == len(cl_hie[pair[1]])*len(cl_hie[pair[0]]))
                    diff = false;
                    if (ok == count) { diff = true; }
                }
                if (diff) { totalOk++; }
                if (!diff && (multiClassStrat == "onevone")) { return false; }
                if (diff && (multiClassStrat == "onevall")) { //all_diff.append(pair)
                    set<string> pair; pair.insert(it->first); pair.insert(itB->first);
                    allDiffs.push_back(pair);
                }
            }//classes
        }//classes
        
        if (multiClassStrat == "onevall") {
            int tot_k = class2SubClasses.size();
            for(map<string, set<string> >::iterator it=class2SubClasses.begin();it!=class2SubClasses.end();it++){
                if (m->getControl_pressed()) { return false; }
                int nk = 0;
                //is this class okay in all comparisons
                for (int h = 0; h < allDiffs.size(); h++) {
                    if (allDiffs[h].count(it->first) != 0) {  nk++; }
                }
                if (nk == (tot_k-1)) { return true;  }//if nk == tot_k-1: return True
            }
            return false;
        }
        
        return true;
    }
    catch(exception& e) {
        m->errorOut(e, "LefseCommand", "testOTUWilcoxon");
        exit(1);
    }
}
//**********************************************************************************************************************
//modelled after lefse.py test_lda_r function
map<int, double> LefseCommand::testLDA(SharedRAbundFloatVectors*& lookup, SharedCLRVectors*& clr, map<int, double> bins, map<string, vector<int> >& class2GroupIndex, map<string, vector<int> >& subClass2GroupIndex) {
    try {
        map<int, double> sigOTUS;
        map<int, double>::iterator it;
        LinearAlgebra linear; Utils util;
    
        int numBins = 0;
        int numGroups = 0;
        if (lookup != NULL) {
            numBins = lookup->getNumBins();
            numGroups = lookup->size(); //lfk
        }else                {
            numBins = clr->getNumBins();
            numGroups = clr->size(); //lfk
        }
        vector< vector<double> > adjustedLookup;
        
        for (int i = 0; i < numBins; i++) {
            if (m->getControl_pressed()) { break; }
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: bin = " + toString(i) + "\n."); }
            
            it = bins.find(i);
            if (it != bins.end()) { //flagged in Kruskal Wallis and Wilcoxon(if we ran it)
                
                if (m->getDebug()) { m->mothurOut("[DEBUG]:flagged bin = " + toString(i) + "\n."); }
                
                //fill x with this OTUs abundances
                vector<float> tempx;
                if (lookup != NULL) { tempx = lookup->getOTU(i);   }
                else                { tempx = clr->getOTU(i);      }
                vector<double> x; for (int h = 0; h < tempx.size(); h++) { x.push_back((double)tempx[h]); }
                
                //go through classes
                for (map<string, vector<int> >::iterator it = class2GroupIndex.begin(); it != class2GroupIndex.end(); it++) {
                    
                    if (m->getDebug()) { m->mothurOut("[DEBUG]: class = " + it->first + "\n."); }
                    
                    //max(float(feats['class'].count(c))*0.5,4)
                    //max(numGroups in this class*0.5, 4.0)
                    double necessaryNum = ((double)((it->second).size())*0.5);
                    if (4.0 > necessaryNum) { necessaryNum = 4.0; }
                    
                    set<double> uniques;
                    for (int j = 0; j < (it->second).size(); j++) { uniques.insert(x[(it->second)[j]]); }
                    
                    //if len(set([float(v[1]) for v in ff if v[0] == c])) > max(float(feats['class'].count(c))*0.5,4): continue
                    if ((double)(uniques.size()) > necessaryNum) {  }
                    else {
                        //feats[k][i] = math.fabs(feats[k][i] + lrand.normalvariate(0.0,max(feats[k][i]*0.05,0.01)))
                        for (int j = 0; j < (it->second).size(); j++) { //(it->second) contains indexes of abundance for this class
                            double sigma = max((x[(it->second)[j]]*0.05), 0.01);
                            x[(it->second)[j]] = abs(x[(it->second)[j]] + linear.normalvariate(0.0, sigma));
                        }
                    }
                }
                adjustedLookup.push_back(x); 
            }
        }
                
        //go through classes
        int minCl = MOTHURMAX;
        map<int, string> indexToClass;
        vector<string> classes;
        for (map<string, vector<int> >::iterator it = class2GroupIndex.begin(); it != class2GroupIndex.end(); it++) {
            //class with minimum number of groups
            if ((it->second).size() < minCl) { minCl = (it->second).size(); }
            for (int i = 0; i < (it->second).size(); i++) { indexToClass[(it->second)[i]] = it->first; }
            classes.push_back(it->first);
        }
        
        
        int fractionNumGroups = numGroups * fBoots; //rfk
        minCl = (int)((float)(minCl*fBoots*fBoots*0.05));
        minCl = max(minCl, 1);
 
        if (m->getDebug()) { m->mothurOut("[DEBUG]: about to start iters. FractionGroups = " + toString(fractionNumGroups) + "\n."); }
        
        vector< vector< vector<double> > > results;//[iters][numComparison][numOTUs]
        for (int j = 0; j < iters; j++) {
            if (m->getControl_pressed()) { return sigOTUS; }
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: iter = " + toString(j) + "\n."); }
            
            //find "good" random vector
            vector<int> rand_s;
            int save = 0;
            
            for (int h = 0; h < 1000; h++) { //generate a vector of length fractionNumGroups with range 0 to numGroups-1
                save = h;
                rand_s.clear();
                
                for (int k = 0; k < fractionNumGroups; k++) {  int index = util.getRandomIndex(numGroups-1); rand_s.push_back(index); }
                if (!contastWithinClassesOrFewPerClass(adjustedLookup, rand_s, minCl, class2GroupIndex, indexToClass)) { h+=1000; save += 1000; } //break out of loop
            }
            
            if (m->getControl_pressed()) { return sigOTUS; }
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: after 1000. \n."); }
            
            if (save < 1000) { m->mothurOut("[WARNING]: Skipping iter " + toString(j+1) + " in LDA test. This can be caused by too few groups per class or not enough contrast within the classes. \n"); }
            else {
                //for each pair of classes
                vector< vector<double> > temp = lda(adjustedLookup, rand_s, indexToClass, classes); //[numComparison][numOTUs]
                if (temp.size() != 0) { results.push_back(temp); }
                if (m->getDebug()) { m->mothurOut("[DEBUG]: after lda. \n."); }
            }
        }
        
        if (results.size() == 0) { return sigOTUS; }
        
        if (m->getControl_pressed()) { return sigOTUS; }
        
        //m = max([numpy.mean([means[k][kk][p] for kk in range(boots)]) for p in range(len(pairs))])
        int k = 0;
        for (it = bins.begin(); it != bins.end(); it++) { //[numOTUs] - need to go through bins so we can tie adjustedLookup back to the binNumber. adjustedLookup[0] ->bins entry[0]. 
            vector<double> averageForEachComparison; averageForEachComparison.resize(results[0].size(), 0.0);
            double maxM = 0.0; //max of averages for each comparison
            for (int j = 0; j < results[0].size(); j++) { //numComparisons
                for (int i = 0; i < results.size(); i++) { //iters
                    averageForEachComparison[j]+= results[i][j][k];
                }
                averageForEachComparison[j] /= (double) results.size();
                if (averageForEachComparison[j] > maxM) { maxM = averageForEachComparison[j]; }
            }
            //res[k] = math.copysign(1.0,m)*math.log(1.0+math.fabs(m),10) 
            double multiple = 1.0; if (maxM < 0.0) { multiple = -1.0; }
            double resK = multiple * log10(1.0+abs(maxM));
            if (resK > ldaThreshold) { sigOTUS[it->first] = resK; }
            k++;
        }
        
        return sigOTUS;
    }
    catch(exception& e) {
        m->errorOut(e, "LefseCommand", "testLDA");
        exit(1);
    }
}
//**********************************************************************************************************************
vector< vector<double> > LefseCommand::getMeans(SharedRAbundFloatVectors*& lookup, SharedCLRVectors*& clr, map<string, vector<int> >& class2GroupIndex) {
    try {
        int numBins = 0;
        if (lookup != NULL) {  numBins = lookup->getNumBins();   }
        else                {  numBins = clr->getNumBins();      }
        
        int numClasses = class2GroupIndex.size();
        vector< vector<double> > means; //[numOTUS][classes]
        means.resize(numBins);
        for (int i = 0; i < means.size(); i++) {  means[i].resize(numClasses, 0.0); }
        
        map<int, string> indexToClass;
        int count = 0;
        //shortcut for vectors below
        map<string, int> quickIndex;
        vector<int> classCounts;
        for (map<string, vector<int> >::iterator it = class2GroupIndex.begin(); it != class2GroupIndex.end(); it++) {
            for (int i = 0; i < (it->second).size(); i++) { indexToClass[(it->second)[i]] = it->first; }
            quickIndex[it->first] = count; count++;
            classCounts.push_back((it->second).size());
        }
        
        for (int i = 0; i < numBins; i++) {
            vector<float> abunds;
            if (lookup != NULL) { abunds = lookup->getOTU(i);   }
            else                { abunds = clr->getOTU(i);      }
            
            for (int j = 0; j < abunds.size(); j++) {
                if (m->getControl_pressed()) { return means; }
                means[i][quickIndex[indexToClass[j]]] += abunds[j];
            }
        }
        
        for (int i = 0; i < numBins; i++) {
            for (int j = 0; j < numClasses; j++) { means[i][j] /= (double) classCounts[j];  }
        }
        
        return means;
    }
    catch(exception& e) {
        m->errorOut(e, "LefseCommand", "getMeans");
        exit(1);
    }
}
//**********************************************************************************************************************
vector< vector<double> > LefseCommand::lda(vector< vector<double> >& adjustedLookup, vector<int> rand_s, map<int, string>& indexToClass, vector<string> classes) {
    try {
        //shortcut for vectors below
        map<string, int> quickIndex;
        for (int i = 0; i < classes.size(); i++) { quickIndex[classes[i]] = i; }
        
        vector<string> randClass; //classes for rand sample
        vector<int> counts; counts.resize(classes.size(), 0);
        for (int i = 0; i < rand_s.size(); i++) {
            string thisClass = indexToClass[rand_s[i]];
            randClass.push_back(thisClass);
            counts[quickIndex[thisClass]]++;
        }

        vector< vector<double> > a; //[numOTUs][numSampled]
        for (int i = 0; i < adjustedLookup.size(); i++) {
            vector<double> temp;
            for (int j = 0; j < rand_s.size(); j++) {
                temp.push_back(adjustedLookup[i][rand_s[j]]);
            }
            a.push_back(temp);
        }
        
        LinearAlgebra linear;
        vector< vector<double> > means; bool ignore;
        vector< vector<double> > scaling = linear.lda(a, randClass, means, ignore); //means are returned sorted, quickIndex sorts as well since it uses a map. means[class][otu] =
        if (ignore) { scaling.clear(); return scaling; }
        if (m->getControl_pressed()) { return scaling; }
        
        vector< vector<double> > w; w.resize(a.size()); //w.unit <- w/sqrt(sum(w^2))
        double denom = 0.0;
        for (int i = 0; i < scaling.size(); i++) { w[i].push_back(scaling[i][0]); denom += (w[i][0]*w[i][0]); }
        denom = sqrt(denom);
        for (int i = 0; i < w.size(); i++) {  w[i][0] /= denom;  } //[numOTUs][1] - w.unit
        
        //robjects.r('LD <- xy.matrix%*%w.unit') [numSampled][numOtus] * [numOTUs][1]
        vector< vector<double> > LD = linear.matrix_mult(linear.transpose(a), w);
        
        //find means for each groups LDs
        vector<double> LDMeans; LDMeans.resize(classes.size(), 0.0); //means[0] -> average for [group0].
        for (int i = 0; i < LD.size(); i++) {  LDMeans[quickIndex[randClass[i]]] += LD[i][0]; } 
        for (int i = 0; i < LDMeans.size(); i++) { LDMeans[i] /= (double) counts[i];  }
   
		//calculate for each comparisons i.e. with groups A,B,C = AB, AC, BC = 3;
        vector< vector<double> > results;// [numComparison][numOTUs]
		for (int i = 0; i < LDMeans.size(); i++) {
			for (int l = 0; l < i; l++) {
                
                if (m->getControl_pressed()) { return scaling; }
                //robjects.r('effect.size <- abs(mean(LD[sub_d[,"class"]=="'+p[0]+'"]) - mean(LD[sub_d[,"class"]=="'+p[1]+'"]))')
                double effectSize = abs(LDMeans[i] - LDMeans[l]);
                //scal = robjects.r('wfinal <- w.unit * effect.size')
                vector<double> compResults;
                for (int j = 0; j < w.size(); j++) { //[numOTUs][1]
                    //coeff = [abs(float(v)) if not math.isnan(float(v)) else 0.0 for v in scal]
                    double coeff = abs(w[j][0]*effectSize); if (isnan(coeff) || isinf(coeff)) { coeff = 0.0; }
                    //gm = abs(res[p[0]][j] - res[p[1]][j]) - res is the means for each group for each otu
                    double gm = abs(means[i][j] - means[l][j]);
                    //means[k][i].append((gm+coeff[j])*0.5)
                    compResults.push_back((gm+coeff)*0.5);
                }
                results.push_back(compResults);
            }
		}
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "LefseCommand", "lda");
        exit(1);
    }
}

//**********************************************************************************************************************
//modelled after lefse.py contast_within_classes_or_few_per_class function
bool LefseCommand::contastWithinClassesOrFewPerClass(vector< vector<double> >& lookup, vector<int> rands, int minCl, map<string, vector<int> > class2GroupIndex, map<int, string> indexToClass) {
    try {
        map<string, int> cls;
        int countFound = 0;
        
        for (int i = 0; i < rands.size(); i++) { //fill cls with the classes represented in the random selection
            for (map<string, vector<int> >::iterator it = class2GroupIndex.begin(); it != class2GroupIndex.end(); it++) {
                if (util.inUsersGroups(rands[i], (it->second))) {
                    map<string, int>::iterator itClass = cls.find(it->first);
                    if (itClass != cls.end()) { itClass->second++; }
                    else { cls[it->first] = 1;  }
                    countFound++;
                }
            }
        }
        
        //sanity check
        if (rands.size() != countFound) { m->mothurOut("oops, should never get here, missing something.\n"); }
        
        if (cls.size() < class2GroupIndex.size()) { return true; } //some classes are not present in sampling
        
        for (map<string, int>::iterator itClass = cls.begin(); itClass != cls.end(); itClass++) {
            if (itClass->second < minCl) { return true; } //this sampling has class count below minimum
        }
        
        //for this otu
        int numBins = lookup.size();
        for (int i = 0; i < numBins; i++) {
            if (m->getControl_pressed()) { break; }
                
            //break up random sampling by class
            map<string, set<double> > class2Values; //maps class name -> set of abunds present in random sampling. F003Early -> 0.001, 0.003... 
            for (int j = 0; j < rands.size(); j++) {
                class2Values[indexToClass[rands[j]]].insert(lookup[i][rands[j]]);
                //rands[j] = index of randomly selected group in lookup, randIndex2Class[rands[j]] = class this group belongs to. lookup[rands[j]]->getAbundance(i) = abundance of this group for this OTU.
            }
            //are the unique values less than we want
            //if (len(set(col)) <= min_cl and min_cl > 1) or (min_cl == 1 and len(set(col)) <= 1):
            for (map<string, set<double> >::iterator it = class2Values.begin(); it != class2Values.end(); it++) {
                if (((it->second).size() <= minCl && minCl > 1) || (minCl == 1 && (it->second).size() <= 1)) {  return true; }
            }
        }
        
        return false;
    }
    catch(exception& e) {
        m->errorOut(e, "LefseCommand", "contastWithinClassesOrFewPerClass");
        exit(1);
    }
}
//**********************************************************************************************************************
int LefseCommand::printResults(vector< vector<double> > means, map<int, double> sigKW, map<int, double> sigLDA, string label, vector<string> classes, vector<string> currentLabels) {
    try {
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(inputfile));
        variables["[distance]"] = label;
        string outputFileName = getOutputFileName("summary",variables);
		ofstream out;
		util.openOutputFile(outputFileName, out);
		outputNames.push_back(outputFileName); outputTypes["summary"].push_back(outputFileName);
        
        //output headers
        out << "OTU\tLogMaxMean\tClass\tLDA\tpValue\n";
        
        string temp = "";
        
        for (int i = 0; i < means.size(); i++) { //[numOTUs][classes]
            //find max mean of classes
            double maxMean = -1.0; string maxClass = "none";
            for (int j = 0; j < means[i].size(); j++) {   if (means[i][j] > maxMean) { maxMean = means[i][j]; maxClass = classes[j]; } }
            
            //str(math.log(max(max(v),1.0),10.0))
            double logMaxMean = 1.0;
            if (maxMean > logMaxMean) { logMaxMean = maxMean; }
            logMaxMean = log10(logMaxMean);
            
            out << currentLabels[i] << '\t' << logMaxMean << '\t';
            if (m->getDebug()) { temp = currentLabels[i] + '\t' + toString(logMaxMean) + '\t'; }
            
            map<int, double>::iterator it = sigLDA.find(i);
            if (it != sigLDA.end()) {
                out << maxClass << '\t' << it->second << '\t' << sigKW[i] << endl; //sigLDA is a subset of sigKW so no need to look
                if (m->getDebug()) { temp += maxClass + '\t' + toString(it->second) + '\t' + toString(sigKW[i]) + '\n'; m->mothurOut(temp); temp = ""; }
            }else { out << '-' << endl; }
        }
        
        out.close();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "LefseCommand", "printResults");
        exit(1);
    }
}
//**********************************************************************************************************************
//printToCoutForRTesting(adjustedLookup, rand_s, class2GroupIndex, numBins);
bool LefseCommand::printToCoutForRTesting(vector< vector<double> >& adjustedLookup, vector<int> rand_s, map<string, vector<int> >& class2GroupIndex, map<int, double> bins, map<string, vector<int> >& subClass2GroupIndex, vector<string> groups, vector<string> currentLabels) {
    try {
        cout << "rand_s = ";
        for (int h = 0; h < rand_s.size(); h++) { cout << rand_s[h] << '\t'; } cout << endl;
        
        //print otu data
        int count = 0;
        for (map<int, double>::iterator it = bins.begin(); it != bins.end(); it++) {
            if (m->getControl_pressed()) { break; }
            
            cout << currentLabels[it->first] << " <- c(";
            for (int h = 0; h < rand_s.size()-1; h++) {  cout << (adjustedLookup[count][rand_s[h]]) << ", "; }
            cout << (adjustedLookup[count][rand_s[rand_s.size()-1]]) << ")\n";
            count++;
        }
        
        string tempOutput = "treatments <- c(";
        for (int h = 0; h < rand_s.size(); h++) {
            //find class this index is in
            for (map<string, vector<int> >::iterator it = class2GroupIndex.begin(); it!= class2GroupIndex.end(); it++) {
                if (util.inUsersGroups(rand_s[h], (it->second)) ) {   tempOutput += "\"" +it->first + "\"" + ","; } //"\"" +it->first + "\""
            }
        }
        tempOutput = tempOutput.substr(0, tempOutput.length()-1);
        tempOutput += ")\n";
        cout << tempOutput;
        
        //print data frame
        tempOutput = "dat <- data.frame(";
        for (map<int, double>::iterator it = bins.begin(); it != bins.end(); it++) {
            if (m->getControl_pressed()) { break; }
            
            tempOutput += "\"" + currentLabels[it->first] + "\"=" + currentLabels[it->first] + ",";
        }
        
        tempOutput += " class=treatments";
        tempOutput += ")\n";
        cout << tempOutput;
                
        tempOutput = "z <- suppressWarnings(mylda(as.formula(class ~ ";
        for (map<int, double>::iterator it = bins.begin(); it != bins.end(); it++) {
            if (m->getControl_pressed()) { break; }
            
            tempOutput +=  currentLabels[it->first] + "+";
        }
        tempOutput = tempOutput.substr(0, tempOutput.length()-1); //rip off extra plus sign
        tempOutput += "), data = dat, tol = 1e-10))";
        cout << tempOutput + "\nz\n";
        cout << "w <- z$scaling[,1]\n"; //robjects.r('w <- z$scaling[,1]')
        cout << "w.unit <- w/sqrt(sum(w^2))\n"; //robjects.r('w.unit <- w/sqrt(sum(w^2))')
        cout << "ss <- dat[,-match(\"class\",colnames(dat))]\n"; //robjects.r('ss <- sub_d[,-match("class",colnames(sub_d))]')
        cout << "xy.matrix <- as.matrix(ss)\n"; //robjects.r('xy.matrix <- as.matrix(ss)')
        cout << "LD <- xy.matrix%*%w.unit\n"; //robjects.r('LD <- xy.matrix%*%w.unit')
        cout << "effect.size <- abs(mean(LD[dat[,\"class\"]==\"'+p[0]+'\"]) - mean(LD[dat[,\"class\"]==\"'+p[1]+'\"]))\n"; //robjects.r('effect.size <- abs(mean(LD[sub_d[,"class"]=="'+p[0]+'"]) - mean(LD[sub_d[,"class"]=="'+p[1]+'"]))')
        cout << "wfinal <- w.unit * effect.size\n"; //scal = robjects.r('wfinal <- w.unit * effect.size')
        cout << "mm <- z$means\n"; //rres = robjects.r('mm <- z$means')

        return true;
    }
    catch(exception& e) {
        m->errorOut(e, "LefseCommand", "printToCoutForRTesting");
        exit(1);
    }
}
/******************************************/






