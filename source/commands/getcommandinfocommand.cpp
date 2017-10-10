/*
 *  getcommandinfo.cpp
 *  Mothur
 *
 *  Created by westcott on 4/6/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "getcommandinfocommand.h"

//**********************************************************************************************************************
vector<string> GetCommandInfoCommand::setParameters(){	
	try {
		CommandParameter poutput("output", "String", "", "", "", "", "","",false,false); parameters.push_back(poutput);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetCommandInfoCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetCommandInfoCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "This command is used by the gui to get the information about current commands available in mothur.\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "GetCommandInfoCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************

GetCommandInfoCommand::GetCommandInfoCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			output = validParameter.validFile(parameters, "output", false);			
			if (output == "not found") {  output = ""; m->mothurOut("You must provide an output filename."); m->mothurOutEndLine(); abort=true; } 
			
		}
	}
	catch(exception& e) {
		m->errorOut(e, "GetCommandInfoCommand", "GetCommandInfoCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetCommandInfoCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		commandFactory = CommandFactory::getInstance();
		
		ofstream out;
		m->openOutputFile(output+".temp", out);
		
		int numNonHidden = 0;
		
		out << "mothurLocation=" << m->getProgramPath() << endl;
		out << "mothurVersion=" << m->getVersion() << endl;
		
		map<string, string> commands = commandFactory->getListCommands();
		map<string, string>::iterator it;
		
		//loop through each command outputting info
		for (it = commands.begin(); it != commands.end(); it++) {
			
			if (m->getControl_pressed()) { m->mothurOut("[ERROR]: did not complete making the file."); m->mothurOutEndLine(); out.close(); m->mothurRemove((output+".temp")); }
			
			Command* thisCommand = commandFactory->getCommand(it->first);
			
			//don't add hidden commands
			if (thisCommand->getCommandCategory() != "Hidden") {
				numNonHidden++;
				
				//general info
				out << "commandName=" << thisCommand->getCommandName() << endl;
				//cout << thisCommand->getCommandName() << " current citation = " << thisCommand->getCitation() << endl;
				out << "commandCategory=" << thisCommand->getCommandCategory() << endl;
				
				//remove /n from help string since gui reads line by line
				string myhelpString = thisCommand->getHelpString();
				string newHelpString = "";
				for (int i = 0; i < myhelpString.length(); i++) { 
					if (myhelpString[i] != '\n') { newHelpString += myhelpString[i]; }
				}
				out << "help=" << newHelpString << endl;
				
				//remove /n from citation string since gui reads line by line
				string mycitationString = thisCommand->getCitation();
				string newCitationString = "";
				for (int i = 0; i < mycitationString.length(); i++) { 
					if (mycitationString[i] != '\n') { newCitationString += mycitationString[i]; }
				}
				out << "citation=" << newCitationString << endl;
				
				out << "description=" << thisCommand->getDescription() << endl;
				
				//outputTypes - makes something like outputTypes=fasta-name-qfile
				map<string, vector<string> > thisOutputTypes = thisCommand->getOutputFiles();
				map<string, vector<string> >::iterator itTypes;
				
				if (thisOutputTypes.size() == 0) { out << "outputTypesNames=0" << endl; }
				else {
					//string types = "";
					//for (itTypes = thisOutputTypes.begin(); itTypes != thisOutputTypes.end(); itTypes++) {	types += itTypes->first + "-";	}
					//rip off last -
					//types = types.substr(0, types.length()-1);
					out << "outputTypesNames=" << thisOutputTypes.size() << endl;
                    
                    for (itTypes = thisOutputTypes.begin(); itTypes != thisOutputTypes.end(); itTypes++) {
                        out << itTypes->first << "=" << thisCommand->getOutputPattern(itTypes->first) << endl;
                    }
				}
				
				vector<string> booleans; vector<string> numbers; vector<string> multiples; vector<string> Strings;
				vector<string> inputGroupNames; map<string, string> inputTypes;
				
				getInfo(thisCommand->getParameters(), booleans, numbers, multiples, Strings, inputGroupNames, inputTypes);
				
				//output booleans
				out << "Boolean=" << booleans.size() << endl;
				for (int i = 0; i < booleans.size(); i++) { out << booleans[i] << endl; }
				
				//output mulitples
				out << "Multiple=" << multiples.size() << endl;
				for (int i = 0; i < multiples.size(); i++) { out << multiples[i] << endl; }
				
				//output numbers
				out << "Numbers=" << numbers.size() << endl;
				for (int i = 0; i < numbers.size(); i++) { out << numbers[i] << endl; }
				
				//output strings
				out << "String=" << Strings.size() << endl;
				for (int i = 0; i < Strings.size(); i++) { out << Strings[i] << endl; }
				
				//output groups
				out << "inputGroupNames=" << inputGroupNames.size() << endl;
				for (int i = 0; i < inputGroupNames.size(); i++) { out << inputGroupNames[i] << endl; }
				
				//output input types
				if (inputTypes.size() == 0) { out << "inputTypes=" << endl; }
				else {
					string types = "";
					for (map<string, string>::iterator it2 = inputTypes.begin(); it2 != inputTypes.end(); it2++) {	types += it2->first + "-";	}
					//rip off last -
					types = types.substr(0, types.length()-1);
					out << "inputTypes=" << types << endl;
					
					for (map<string, string>::iterator it2 = inputTypes.begin(); it2 != inputTypes.end(); it2++) {	
						out << it2->first << "=" << it2->second << endl;
					}
				}
				
			}
		}
		
		out.close();
		
		ofstream out2;
		m->openOutputFile(output, out2);
		out2 << numNonHidden << endl; 
		out2.close();
		
		m->appendFiles(output+".temp", output);
		m->mothurRemove((output+".temp"));
	
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOut(output); m->mothurOutEndLine();	
		m->mothurOutEndLine();
		
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "GetCommandInfoCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetCommandInfoCommand::getInfo(vector<CommandParameter> para, vector<string>& booleans, vector<string>& numbers, vector<string>& multiples, vector<string>& strings, vector<string>& inputGroupNames, map<string, string>& inputTypes){
	try {
		
		map<string, set<string> > groups;
		map<string,  set<string> >::iterator itGroups;
		
		for (int i = 0; i < para.size(); i++) {
			if ((para[i].name == "inputdir") || (para[i].name == "outputdir")) {} //ignore
			else {
                string important = "|F";
                if (para[i].important || para[i].required) { important = "|T"; }
                
                string outputType = "|none";
                if (para[i].outputTypes != "") { outputType = "|" + para[i].outputTypes; }
                
				if (para[i].type == "Boolean") {
					string temp = para[i].name + "=" + para[i].optionsDefault + important + outputType;
					booleans.push_back(temp);
				}else if (para[i].type == "Multiple") {
					string multAllowed = "F";
					if (para[i].multipleSelectionAllowed) { multAllowed = "T"; }
					string temp = para[i].name + "=" + para[i].options + "|" + para[i].optionsDefault + "|" + multAllowed + important + outputType;
					multiples.push_back(temp);
				}else if (para[i].type == "Number") {
					string temp = para[i].name + "=" + para[i].optionsDefault + important + outputType;
					numbers.push_back(temp);
				}else if (para[i].type == "String") {
					string temp = para[i].name + "=" + para[i].optionsDefault + important + outputType;
					strings.push_back(temp);
				}else if (para[i].type == "InputTypes") {
					string required = "F";
					if (para[i].required) { required = "T"; }
					string temp = required + important + "|" + para[i].chooseOnlyOneGroup + "|" + para[i].chooseAtLeastOneGroup + "|" + para[i].linkedGroup + outputType;
					inputTypes[para[i].name] = temp;
					
					//add choose only one groups
                    vector<string> tempGroups;
                    m->splitAtDash(para[i].chooseOnlyOneGroup, tempGroups);
                    for (int l = 0; l < tempGroups.size(); l++) {
                        groups[tempGroups[l]].insert(para[i].name);
                    }
					tempGroups.clear();
                    
					//add at least one group names
                    m->splitAtDash(para[i].chooseAtLeastOneGroup, tempGroups);
                    for (int l = 0; l < tempGroups.size(); l++) {
                        groups[tempGroups[l]].insert(para[i].name);
                    }
					tempGroups.clear();
					
					
					//add at linked group names
                    m->splitAtDash(para[i].linkedGroup, tempGroups);
                    for (int l = 0; l < tempGroups.size(); l++) {
                        groups[tempGroups[l]].insert(para[i].name);
                    }
					tempGroups.clear();
						  
				}else { m->mothurOut("[ERROR]: " + para[i].type + " is an unknown parameter type, please correct."); m->mothurOutEndLine(); }
			}
		}
		
		for (itGroups = groups.begin(); itGroups != groups.end(); itGroups++) {
			if (itGroups->first != "none") {
				set<string> tempNames = itGroups->second;
				string temp = itGroups->first + "=";
				for (set<string>::iterator itNames = tempNames.begin(); itNames != tempNames.end(); itNames++) {
					temp += *itNames + "-";
				}
				//rip off last -
				temp = temp.substr(0, temp.length()-1);
				inputGroupNames.push_back(temp);
			}
		}
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "GetCommandInfoCommand", "getInfo");
		exit(1);
	}
}
//**********************************************************************************************************************/
