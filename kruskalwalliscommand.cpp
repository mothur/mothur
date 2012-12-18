/* 
 * File:   kruskalwalliscommand.cpp
 * Author: kiverson
 *
 * Created on June 26, 2012, 11:06 AM
 */

#include "kruskalwalliscommand.h"

//**********************************************************************************************************************
vector<string> KruskalWallisCommand::setParameters(){	
	try {
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false,true); parameters.push_back(pgroups);
        CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none","summary",false,true,true); parameters.push_back(pshared);	
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "KruskalWallisCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string KruskalWallisCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The kruskalwallis command parameter options are \n";
        helpString += "Kruskal–Wallis one-way analysis of variance is a non-parametric method for testing whether samples originate from the same distribution.";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "KruskalWallisCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string KruskalWallisCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "summary") {  pattern = "[filename],cooccurence.summary"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "KruskalWallisCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
KruskalWallisCommand::KruskalWallisCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
        vector<string> tempOutNames;
		outputTypes["summary"] = tempOutNames;

	}
	catch(exception& e) {
		m->errorOut(e, "KruskalWallisCommand", "KruskalWallisCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
KruskalWallisCommand::KruskalWallisCommand(string option) {
	try {
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
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
            
            //get shared file
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }	
			else if (sharedfile == "not found") { 
				//if there is a current shared file, use it
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current sharedfile and the shared parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setSharedFile(sharedfile); }
            
            //if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(sharedfile);		}
                    
            groups = validParameter.validFile(parameters, "groups", false);   
            if (groups == "not found") { groups = "";   }
            else { 
            m->splitAtDash(groups, Groups); 
            }   
            m->setGroups(Groups);
				
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
			}
		
            vector<string> tempOutNames;
            outputTypes["summary"] = tempOutNames;


		}

	}
	catch(exception& e) {
		m->errorOut(e, "KruskalWallisCommand", "KruskalWallisCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int KruskalWallisCommand::execute(){
	try {
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        InputData* input = new InputData(sharedfile, "sharedfile");
        vector<SharedRAbundVector*> lookup = input->getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
        
	
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;

        ofstream out;
        map<string,string> variables;
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(sharedfile));
		string outputFileName = getOutputFileName("summary",variables);
        m->openOutputFile(outputFileName, out);
        outputNames.push_back(outputFileName);  outputTypes["summary"].push_back(outputFileName);
        out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
        out << "H\tpvalue\n";
        
        //math goes here
        
        int N = m->getNumGroups();
        double H;
        double tmp = 0.0;
        vector<groupRank> vec;
        vector<string> groups = m->getGroups();
        string group;
        int count;
        double sum;
                
        //merge all groups into a vector
        
        
        
        //rank function here
        assignRank(vec);
        
        //populate counts and ranSums vectors
        for (int i=0;i<N;i++) {
            count = 0;
            sum = 0;
            group = groups[i];
            for(int j;j<vec.size();j++) {
                if (vec[j].group == group) {
                    count++;
                    sum = sum + vec[j].rank;
                }
            }
            counts[i] = count;
            rankSums[i] = sum;
        }
        
        //test statistic
        for (int i=0;i<N;i++) { tmp = tmp + (pow(rankSums[i],2) / counts[i]); }
        
        H = (12 / (N*(N+1))) * tmp - (3*(N+1));
        
        //ss = tmp - pow(accumulate(rankSums.begin(), rankSums.end(), 0), 2);
        
        //H = ss / ( (N * (N + 1))/12 );
        
        //correction for ties?
        
        //p-value calculation
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "KruskalWallisCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
void KruskalWallisCommand::assignRank(vector<groupRank> &vec) {
    try {
        double rank = 1;
        double numRanks, avgRank, j;
        vector<groupRank>::iterator it, oldit;

        sort (vec.begin(), vec.end(), comparevalue);

        it = vec.begin();

        while ( it != vec.end() ) {
            j = rank;
            oldit = it;
            if (!equalvalue(*it, *(it+1))) {
                (*it).rank = rank; 
                rank = rank+1; 
                it++; }
            else {
                while(equalrank(*it, *(it+1))) {
                    j = j + (j+1);
                    rank++;
                    it++;
                }
                numRanks = double (distance(oldit, it));
                avgRank = j / numRanks;
                while(oldit != it) {
                    (*oldit).rank = avgRank;
                    oldit++;
                }
            }

        }
        

    }
    catch(exception& e) {
		m->errorOut(e, "KruskalWallisCommand", "getRank");
		exit(1);
	}
    
}
//**********************************************************************************************************************
void KruskalWallisCommand::assignValue(vector<groupRank> &vec) {
    
}
//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************