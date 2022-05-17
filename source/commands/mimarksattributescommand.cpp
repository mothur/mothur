//
//  mimarksattributescommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/17/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#include "mimarksattributescommand.h"

//**********************************************************************************************************************
vector<string> MimarksAttributesCommand::setParameters(){
    try {
        CommandParameter pxml("xml", "InputTypes", "", "", "none", "none", "none","summary",false,false,true); parameters.push_back(pxml);
        CommandParameter psets("package", "String", "", "", "", "", "","",false,false); parameters.push_back(psets);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false;
        
        vector<string> myArray;
        for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
        return myArray;
    }
    catch(exception& e) {
        m->errorOut(e, "MimarksAttributesCommand", "setParameters");
        exit(1);
    }
}
//**********************************************************************************************************************
string MimarksAttributesCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "source") {  pattern = "[filename],source"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "MimarksAttributesCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
string MimarksAttributesCommand::getHelpString(){
    try {
        string helpString = "";
        helpString += "Reads bioSample Attributes xml and generates source for get.mimarkspackage command. Only parameter required is xml.\n";
        helpString += "The package parameter allows you to set the package you want. Default MIMARKS.survey.\n";
        return helpString;
    }
    catch(exception& e) {
        m->errorOut(e, "MimarksAttributesCommand", "getHelpString");
        exit(1);
    }
}
//**********************************************************************************************************************
MimarksAttributesCommand::MimarksAttributesCommand(string option) : Command()  {
    try {

        if(option == "help") { help(); abort = true; calledHelp = true; }
        else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
        
        else {
            OptionParser parser(option, setParameters());
            map<string,string> parameters = parser.getParameters();
            
            ValidParameters validParameter;
            xmlFile = validParameter.validFile(parameters, "xml");
            if (xmlFile == "not open") { abort = true; }
            else if (xmlFile == "not found") {  xmlFile = ""; abort=true; m->mothurOut("You must provide an xml file. It is required.\n");  }
            
            selectedPackage = validParameter.valid(parameters, "package");
            if (selectedPackage == "not found") { selectedPackage = "MIMARKS.survey."; }
            
            if (outputdir == ""){    outputdir = util.hasPath(xmlFile);		}
        }
    }
    catch(exception& e) {
        m->errorOut(e, "MimarksAttributesCommand", "MimarksAttributesCommand");
        exit(1);
    }
}
//**********************************************************************************************************************

int MimarksAttributesCommand::execute(){
    try {
        
        if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        ifstream in; util.openInputFile(xmlFile, in);
        string header = util.getline(in); gobble(in);
        
        if (header != "<BioSampleAttributes>") { m->mothurOut("[ERROR]: " + header + " is not a bioSample attribute file.\n"); m->setControl_pressed(true); }
        
        map<string, Group> categories;
        map<string, Group>::iterator it;
        
        while (!in.eof()) {
            if (m->getControl_pressed()) { in.close(); return 0; }
            
            Attribute attribute = readAttribute(in);
            
            if (attribute.name != "") {
                if (m->getDebug()) {
                    m->mothurOut("[DEBUG]: name=" + attribute.name + " harmonizedName=" + attribute.harmonizedName + " format=" + attribute.format + " description=" + attribute.description + " package=" + attribute.getPackagesString() + "\n");
                }
                
                if (attribute.format == "") { attribute.format = "{none}"; }
                if (attribute.description == "") { attribute.description = "none"; }
                
                
                for (int i = 0; i < attribute.packages.size(); i++) {
                    for (int j = 0; j < attribute.packages[i].groupName.length(); j++) {
                        if (attribute.packages[i].groupName[j] == '-') { attribute.packages[i].groupName[j] = '_'; }
                    }
                    
                    it = categories.find(attribute.packages[i].groupName);
                    if (it != categories.end()) { //we already have this category, ie air, soil...
                        if (attribute.packages[i].name == (it->second).packageName) { //add attribute to category
                            (it->second).values[attribute.harmonizedName].required = attribute.packages[i].required;
                            (it->second).values[attribute.harmonizedName].format = attribute.format;
                            string newDescription = "";
                            for (int j = 0; j < attribute.description.length(); j++) {
                                if (attribute.description[j] == '"') { newDescription += "\\\""; }
                                else { newDescription += attribute.description[j]; }
                            }
                            (it->second).values[attribute.harmonizedName].description = newDescription;
                        }
                    }else {
                        if ((attribute.packages[i].groupName == "\"Built\"") || (attribute.packages[i].groupName == "\"Nucleic Acid Sequence Source\"")) {}
                        else {
                            Group thisGroup(attribute.packages[i].name);
                            thisGroup.values[attribute.harmonizedName].required = attribute.packages[i].required;
                            thisGroup.values[attribute.harmonizedName].format = attribute.format;
                            string newDescription = "";
                            for (int j = 0; j < attribute.description.length(); j++) {
                                if (attribute.description[j] == '"') { newDescription += "\\\""; }
                                else { newDescription += attribute.description[j]; }
                            }
                            thisGroup.values[attribute.harmonizedName].description = newDescription;
                            categories[attribute.packages[i].groupName] = thisGroup;
                        }
                    }
                }
            }
        }
        in.close();
        
        string requiredByALL = "*sample_name\t*description\t*sample_title\t*seq_methods\t*organism";
        string rFormatALL = "#{text}\t{text}\t{text}\t{text}\t{controlled vacabulary}";
        string rDescriptionALL = "#{sample name}\t{description of sample}\t{sample title}\t{description of library_construction_protocol}\t{http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=408169&lvl=3&keep=1&srchmode=1&unlock}";
        string environment = "\"Environment\"";
        it = categories.find(environment);
        if (it != categories.end()) {
            map<string, Value>::iterator itValue = (it->second).values.begin();
            if (itValue->second.required) {
                requiredByALL += "\t*" + itValue->first;
                rFormatALL += "\t{" + (itValue->second.format) + "}";
                 rDescriptionALL += "\t{" + (itValue->second.description) + "}";
            }
            itValue++;
            
            for (; itValue != (it->second).values.end(); itValue++) {
                if (itValue->second.required) {
                    requiredByALL += "\t*" + itValue->first;
                    rFormatALL += "\t{" + (itValue->second.format)  + "}";
                    rDescriptionALL += "\t{" + (itValue->second.description) + "}";
                }
            }
        }
        
        ofstream out;
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(xmlFile));
        string outputFileName = getOutputFileName("source",variables);
        outputNames.push_back(outputFileName); outputTypes["source"].push_back(outputFileName);
        util.openOutputFile(outputFileName, out);
        
        //create outputs
        string requiredValues = requiredByALL; string nonRequiredValues = "";
        string rFormat = rFormatALL; string nonRFormat = "";
        string rDescription = rDescriptionALL; string nonRDescription = "";
        it = categories.begin();
        map<string, Value>::iterator itValue = (it->second).values.begin();
        if (itValue->second.required) {
            requiredValues += "\t*" + itValue->first;
            rFormat += "\t{" + (itValue->second.format) + "}";
            rDescription += "\t{" + (itValue->second.description) + "}";
        }else {
            nonRequiredValues += itValue->first;
            nonRFormat += "{" + itValue->second.format + "}";
            nonRDescription += "{" + (itValue->second.description) + "}";
        }
        itValue++;
        for (; itValue != (it->second).values.end(); itValue++) {
            if (itValue->second.required) {
                requiredValues += "\t*" + itValue->first;
                rFormat += "\t{" + (itValue->second.format) + "}";
                rDescription += "\t{" + (itValue->second.description) + "}";
            }else {
                nonRequiredValues += "\t" + itValue->first;
                nonRFormat += "\t{" + itValue->second.format + "}";
                nonRDescription += "\t{" + (itValue->second.description) + "}";
            }
        }
        
        out << "if (package == " + it->first + ") {\n";
        out << "\tout << \"#" + it->second.packageName + "\" << endl;\n";
        out << "\t if (requiredonly) {\n";
        out << "\t\tout << \"" + rDescription + "\" << endl;\n";
        out << "\t\tout << \"" + rFormat + "\" << endl;\n";
        out << "\t\tout << \"" + requiredValues + "\" << endl;\n";
        out << "\t}else {\n";
        out << "\t\tout << \"" + rDescription + '\t' + nonRDescription + "\" << endl;\n";
        out << "\t\tout << \"" + rFormat + '\t' + nonRFormat + "\" << endl;\n";
        out << "\t\tout << \"" + requiredValues + '\t' + nonRequiredValues + "\" << endl;\n";
        out << "\t}\n";
        out << "}";
        
        it++;
        for (; it != categories.end(); it++) {
            if ((it->first  == "\"Environment\"")) {}
            else {
                //create outputs
                string requiredValues = requiredByALL; string nonRequiredValues = "";
                string rFormat = rFormatALL; string nonRFormat = "";
                string rDescription = rDescriptionALL; string nonRDescription = "";
                map<string, Value>::iterator itValue = (it->second).values.begin();
                if (itValue->second.required) {
                    requiredValues += "\t*" + itValue->first;
                    rFormat += "\t{" + (itValue->second.format)+ "}";
                    rDescription += "\t{" + (itValue->second.description) + "}";
                }else {
                    nonRequiredValues += itValue->first;
                    nonRFormat += "{" + itValue->second.format+ "}";
                    nonRDescription += "{" + (itValue->second.description) + "}";
                }
                itValue++;
                for (; itValue != (it->second).values.end(); itValue++) {
                    if (itValue->second.required) {
                        requiredValues += "\t*" + itValue->first;
                        rFormat += "\t{" + (itValue->second.format)+ "}";
                        rDescription += "\t{" + (itValue->second.description) + "}";
                    }else {
                        nonRequiredValues += "\t" + itValue->first;
                        nonRFormat +=  "\t{" + itValue->second.format+ "}";
                        nonRDescription += "\t{" + (itValue->second.description) + "}";
                    }
                }
                
                out << "else if (package == " + it->first + ") {\n";
                out << "\tout << \"#" + it->second.packageName + "\" << endl;\n";
                out << "\t if (requiredonly) {\n";
                out << "\t\tout << \"" + rDescription + "\" << endl;\n";
                out << "\t\tout << \"" + rFormat + "\" << endl;\n";
                out << "\t\tout << \"" + requiredValues + "\" << endl;\n";
                out << "\t}else {\n";
                out << "\t\tout << \"" + rDescription + '\t' + nonRDescription + "\" << endl;\n";
                out << "\t\tout << \"" + rFormat + '\t' + nonRFormat + "\" << endl;\n";
                out << "\t\tout << \"" + requiredValues + '\t' + nonRequiredValues + "\" << endl;\n";
                out << "\t}\n";
                out << "}";
            }
        }
        
        out << endl << endl;
        it = categories.begin();
        out << "if ((package == " << it->first << ") ";
        it++;
        for (; it != categories.end(); it++) {
            out << "|| (package == " << it->first << ") ";
        }
        out << ") {}\n\n";
        
        out << "vector<string> requiredFieldsForPackage;\n";
        vector<string> rAll;
        util.splitAtChar(requiredByALL, rAll, '\t');
        for (int i = 0; i < rAll.size(); i++) {
            out << "requiredFieldsForPackage.push_back(\"" + rAll[i].substr(1) + "\");\n";
        }
        out << "\n\n";
        
        
        for (it = categories.begin(); it != categories.end(); it++) {
            out << "if (packageType == \"" << it->second.packageName << "\") {";
            for (map<string, Value>::iterator itValue = (it->second).values.begin(); itValue != (it->second).values.end(); itValue++) {
                if (itValue->second.required) {
                    out << "\trequiredFieldsForPackage.push_back(\"" +  itValue->first + "\");";
                }
            }
            out << "}\n";
        }
        out.close();
        
        m->mothurOutEndLine();
        m->mothurOut("Output File Names: \n"); 
        for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}	
        m->mothurOutEndLine();
        
        return 0;		
    }
    
    catch(exception& e) {
        m->errorOut(e, "MimarksAttributesCommand", "execute");
        exit(1);
    }
}
//**********************************************************************************************************************

Attribute MimarksAttributesCommand::readAttribute(ifstream& in){
    try {
        
        
        //read <Attribute>
        string header = util.getline(in); gobble(in);
        
        if (header == "</BioSampleAttributes>") { Attribute temp; return temp; }
        
        if (header != "<Attribute>") {  m->mothurOut("[ERROR]: " + header + ", expected '<Attribute>' in file.\n"); m->setControl_pressed(true);  }
        
        //read name
        //<Name>wastewater type</Name>
        gobble(in);
        string name = util.getline(in); gobble(in);
        trimTags(name);
        
        //read hamonized name
        //<HarmonizedName>wastewater_type</HarmonizedName>
        gobble(in);
        string hname = util.getline(in); gobble(in);
        trimTags(hname);
        
        //read description
        //<Description>
        //the origin of wastewater such as human waste, rainfall, storm drains, etc.
        //</Description>
        string description = "";
        unsigned long long spot = in.tellg();
        gobble(in);
        char c = in.get(); c = in.get();
        if (c == 'D') { //description
            description += "<D";
            while (!in.eof()) {
                gobble(in);
                string thisLine = util.getline(in); gobble(in);
                description += thisLine;
                if (thisLine.find("</Description>") != string::npos)  { break; }
            }
            trimTags(description);
        }else { //package
            in.seekg(spot);
        }
        
        //read format
        //<Format>{text}</Format>
        spot = in.tellg();
        gobble(in);
        c = in.get(); c = in.get();
        string format = "";
        if (c == 'F') { //format
            format += "<F" + util.getline(in); gobble(in);
            if (format.find("</Format>") == string::npos) { //format is not on oneline
                while (!in.eof()) {
                    gobble(in);
                    string thisLine = util.getline(in); gobble(in);
                    format += thisLine;
                    if (thisLine.find("</Format>") != string::npos) { break; }
                }
            }
            trimTags(format);
        }else { //package
            in.seekg(spot);
        }
        
        Attribute attribute(hname, description, name, format);
        
        //read Synonym - may be none
        //<Synonym>ref biomaterial</Synonym>
        bool FirstTime = true;
        while (!in.eof()) {
            unsigned long long thisspot = in.tellg();
            gobble(in);
            char c = in.get(); c = in.get();
            if (c == 'S') { //synonym
                FirstTime = false;
                util.getline(in); gobble(in);
            }else { //package
                if (FirstTime) { in.seekg(spot); }
                else { in.seekg(thisspot); }
                break;
            }
        }

        
        //read packages - may be none
        //<Package use="optional" group_name="Air">MIGS.ba.air.4.0</Package>
        while (!in.eof()) {
            string package = util.getline(in); gobble(in);
            if (package == "</Attribute>") { break; }
            else {
                Package thisPackage = parsePackage(package);
                if (thisPackage.groupName != "ignore") { attribute.packages.push_back(thisPackage); }
            }
        }
        
        return attribute;
    }
    
    catch(exception& e) {
        m->errorOut(e, "MimarksAttributesCommand", "execute");
        exit(1);
    }
}
//**********************************************************************************************************************

Package MimarksAttributesCommand::parsePackage(string package){
    try {
        string openingTag = trimTags(package);
        Package thispackage; thispackage.name = package;
        
        //only care about packages from our selection
        if (thispackage.name.find(selectedPackage) == string::npos) {  thispackage.groupName = "ignore"; return thispackage; }
        
        int pos = openingTag.find("use");
        if (pos != string::npos) {
            //read required or not
            string use = openingTag.substr(openingTag.find_first_of("\""), 11);
            if (use == "\"mandatory\"") { thispackage.required = true; }
        }else {
            m->mothurOut("[ERROR]: parsing error - " + openingTag + ". Expeacted something like <Package use=\"optional\" group_name=\"Air\"> in file.\n"); m->setControl_pressed(true);  return thispackage;
        }
        
        //selectedPackage = MIMARKS.survey.
        pos = package.find(selectedPackage);
        if (pos != string::npos) {
            //read groupname
            string group = package.substr(pos+15);
            group = group.substr(0, (group.find_first_of(".")));
            thispackage.groupName = "\"" + group + "\"";
        }else {
             thispackage.groupName = "ignore";
        }
        
        return thispackage;
    }
    
    catch(exception& e) {
        m->errorOut(e, "MimarksAttributesCommand", "parsePackage");
        exit(1);
    }
}
//**********************************************************************************************************************

string MimarksAttributesCommand::trimTags(string& value){
    try {
        string forwardTag = "";
        string thisValue = "";
        int openCarrot = 0;
        int closedCarrot = 0;
        
        
        for (int i = 0; i < value.length(); i++) {
            if (m->getControl_pressed()) { return forwardTag; }
            
            if (value[i] == '<')         { openCarrot++;     }
            else if (value[i] == '>')    { closedCarrot++;   }
            
            //you are reading front tag
            if ((openCarrot == 1) && (closedCarrot == 0)) { forwardTag += value[i];  }
            
            if (openCarrot == closedCarrot) { //reading value
                if (value[i] != '>') { thisValue += value[i]; }
            }
            
            if (openCarrot > 1) { break; }
        }
    
        value = thisValue;
        return (forwardTag + '>');
    }
    
    catch(exception& e) {
        m->errorOut(e, "MimarksAttributesCommand", "trimTags");
        exit(1);
    }
}
//**********************************************************************************************************************


