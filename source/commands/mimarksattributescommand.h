//
//  mimarksattributescommand.h
//  Mothur
//
//  Created by Sarah Westcott on 3/17/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__mimarksattributescommand__
#define __Mothur__mimarksattributescommand__

#include "command.hpp"

struct Package {
    bool required;
    string groupName;
    string name;
    Package() { required=false; groupName=""; name=""; }
    Package(bool r, string g, string n) : required(r), groupName(g), name(n) {}
    ~Package() {}
    
    string getPackageString() {
        string r = "mandatory";
        if (!required) { r = "optional";  }
        string packageString = name + '\t' + groupName + '\t' + r;
        return packageString;
    }
};

struct Value {
    bool required;
    string format, description;
    
    Value() { format=""; description=""; required=false; }
    Value(bool r, string d, string f) : format(f), description(d), required(r) {}
    ~Value() {}

};

struct Group {
    string packageName;
    map<string, Value> values;
    
    Group() { packageName = ""; }
    Group(string p) :  packageName(p) {}
    ~Group() {}
};

struct Attribute {
    string name, harmonizedName, description, format;
    vector<Package> packages;
    
    string getPackagesString() {
        string packagesString = "";
        for (int i = 0; i < packages.size(); i++) {
            packagesString += packages[i].getPackageString() + "\n";
        }
        return packagesString;
    }
    
    Attribute() { format=""; description=""; harmonizedName=""; name=""; }
    Attribute(string hn, string d, string n, string f) : format(f), harmonizedName(hn), name(n), description(d) {}
    ~Attribute() {}
};

/**************************************************************************************************/

class MimarksAttributesCommand : public Command {
public:
    MimarksAttributesCommand(string);
   ~MimarksAttributesCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "mimarks.attributes";			}
    string getCommandCategory()		{ return "Hidden";		}
    
    string getOutputPattern(string);
    
    string getHelpString();
    string getCitation() { return "http://www.mothur.org/wiki/Mimarks.attributes"; }
    string getDescription()		{ return "Reads bioSample Attributes xml and generates source for get.mimarkspackage command."; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:

    Attribute readAttribute(ifstream& in);
    Package parsePackage(string package);
    string trimTags(string& value);
    
    bool abort;
    string  xmlFile, selectedPackage;
    vector<string> outputNames;
};

/**************************************************************************************************/



#endif /* defined(__Mothur__mimarksattributescommand__) */
