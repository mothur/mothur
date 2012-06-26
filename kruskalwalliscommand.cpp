/* 
 * File:   kruskalwalliscommand.cpp
 * Author: kiverson
 *
 * Created on June 26, 2012, 11:06 AM
 */
#include "kruskalwalliscommand.h"

//**********************************************************************************************************************
class groupRank {
public:
    string group;
    double value;
    double rank;
};
//**********************************************************************************************************************
vector<string> KruskalWallisCommand::setParameters(){	
	try {
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
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
string KruskalWallisCommand::getOutputFileNameTag(string type, string inputName=""){	
	try {
        string outputFileName = "";
		map<string, vector<string> >::iterator it;
        
        //is this a type this command creates
        it = outputTypes.find(type);
        if (it == outputTypes.end()) {  m->mothurOut("[ERROR]: this command doesn't create a " + type + " output file.\n"); }
        else {
            if (type == "summary") {  outputFileName =  "cooccurence.summary"; }
            else { m->mothurOut("[ERROR]: No definition for type " + type + " output file tag.\n"); m->control_pressed = true;  }
        }
        return outputFileName;
	}
	catch(exception& e) {
		m->errorOut(e, "KruskalWallisCommand", "getOutputFileNameTag");
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
        
        //math goes here
        
        int N;
        double ss, H;
        double tmp = 0.0;
                
        //merge all groups into a vector
        //rank function here
        
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
multimap<double,double> KruskalWallisCommand::getRank(vector<groupRank> vec) {
    try {
        multimap<double,double> rankMap;
        double rank = 1;
        double previous;
        double tie = 0.0;
        int tiecount = 0;

        sort (vec.begin(), vec.end());

        for (int i=0;i<vec.size();i++) {
            if (vec[i] != previous) { rankMap[rank] = vec[i]; }
            else {tie = tie + rank; tiecount++;}
            rank++;
            previous = vec[i];
        }
    }
    catch(exception& e) {
		m->errorOut(e, "KruskalWallisCommand", "getRank");
		exit(1);
	}
    
}
//**********************************************************************************************************************

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************