//
//  biomsimple.cpp
//  Mothur
//
//  Created by Sarah Westcott on 10/26/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "biomsimple.hpp"

/**************************************************************************************************/
BiomSimple::BiomSimple() : Biom(){
    try {
       
        version = "Biological Observation Matrix 0.9.1";
    }
    catch(exception& e) {
        m->errorOut(e, "BiomSimple", "BiomSimple");
        exit(1);
    }
}

/**************************************************************************************************/
BiomSimple::BiomSimple(string fname, string l) : Biom("Biological Observation Matrix 0.9.1"){
    try {
        label = l;
        read(fname);
    }
    catch(exception& e) {
        m->errorOut(e, "BiomSimple", "BiomSimple");
        exit(1);
    }
}
/**************************************************************************************************/
void BiomSimple::read(string fname){
    try {
       
        /*{
         "id":"/Users/SarahsWork/Desktop/release/temp.job2.shared-unique",
         "format": "Biological Observation Matrix 0.9.1",
         "format_url": "http://biom-format.org",
         "type": "OTU table",
         "generated_by": "mothur1.44.0",
         "date": "Tue Apr 17 13:12:07 2020",
         
         rows represent OTUS
         columns represent samples
         
         */
        
        ifstream in; util.openInputFile(fname, in);
        
        matrixFormat = ""; matrixElementType = "";
        vector<string> otuNames;  vector<string> groupNames;
        map<string, string> fileLines;
        //vector<string> names;
        int numOTUs, numCols;
        bool hasTaxonomy;
       
        numOTUs = 0; numCols = 0; maxLevel = 0;
        int shapeNumRows = 0; int shapeNumCols = 0;
        
        int countOpenBrace = 0; int countClosedBrace = 0;
        int closeParen = 0; int openParen = -1; //account for opening brace
        bool ignoreCommas = false; bool atComma = false;
        
        string line = "";
        bool printHeaders = true;
        
        while (!in.eof()) { //split file by tags, so each "line" will have something like "id":"/Users/SarahsWork/Desktop/release/final.tx.1.subsample.1.pick.shared-1"
            if (m->getControl_pressed()) { break; }
            
            char c = in.get(); util.gobble(in);
            
            if (c == '[')               { countOpenBrace++;     }
            else if (c == ']')          { countClosedBrace++;   }
            else if (c == '{')          { openParen++;          }
            else if (c == '}')          { closeParen++;         }
            else if ((!ignoreCommas) && (c == ','))          { atComma = true;       }
            
            if ((countOpenBrace != countClosedBrace) && (countOpenBrace != countClosedBrace)) { ignoreCommas = true;  }
            else if ((countOpenBrace == countClosedBrace) && (countOpenBrace == countClosedBrace)) { ignoreCommas = false;  }
            if (atComma && !ignoreCommas) {
                if (fileLines.size() == 0) { //clip first {
                    line = line.substr(1);
                }
                string tag = getTag(line);
                fileLines[tag] = line;
                
                line = "";
                atComma = false;
                ignoreCommas = false;
                
            }else {  line += c;  }
            
        }
        if (line != "") {
            line = line.substr(0, line.length()-1);
            string tag = getTag(line);
            fileLines[tag] = line;
        }
        in.close();
        
        //check for required fields
        map<string, string>::iterator it;
        it = fileLines.find("type");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a type provided.\n"); }
        else {
            string thisLine = it->second;
            tableType = getTag(thisLine);
        }
        
        if (m->getControl_pressed()) { return; }
        
        it = fileLines.find("matrix_type");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a matrix_type provided.\n"); }
        else {
            string thisLine = it->second;
            matrixFormat = getTag(thisLine);
            if ((matrixFormat != "sparse") && (matrixFormat != "dense")) { m->mothurOut("[ERROR]: " + matrixFormat + " is not a valid biom matrix_type for mothur. Types allowed are sparse and dense.\n"); m->setControl_pressed(true); }
        }
        
        if (m->getControl_pressed()) { return; }
        
        it = fileLines.find("matrix_element_type");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a matrix_element_type provided.\n"); }
        else {
            string thisLine = it->second;
            matrixElementType = getTag(thisLine);
            if ((matrixElementType != "int") && (matrixElementType != "float")) { m->mothurOut("[ERROR]: " + matrixElementType + " is not a valid biom matrix_element_type for mothur. Types allowed are int and float.\n"); m->setControl_pressed(true); }
        }
        
        if (m->getControl_pressed()) { return; }
        
        map<string, string> otuTaxonomies;
        it = fileLines.find("rows");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a rows provided.\n"); }
        else {
            maxLevel = 0;
            string thisLine = it->second;
            
            bool hasTaxonomy = false;
            vector< vector<string> > results = extractTaxonomyData(thisLine, numOTUs, hasTaxonomy);
            
            if ((tableType == "Taxon table") || (tableType == "Taxontable")) {
                vector<string> taxonomies = results[0];
                
                //create OTU names
                string snumBins = toString(numOTUs);
                for (int i = 0; i < numOTUs; i++) {
                    
                    //if there is a bin label use it otherwise make one
                    string binLabel = "OTU";
                    string sbinNumber = toString(i+1);
                    if (sbinNumber.length() < snumBins.length()) {
                        int diff = snumBins.length() - sbinNumber.length();
                        for (int h = 0; h < diff; h++) { binLabel += "0"; }
                    }
                    binLabel += sbinNumber;
                    
                    otuNames.push_back(binLabel);
                    otuTaxonomies[otuNames[i]] = taxonomies[i];
                }
               
            }else{
                otuNames = results[0];
                if (hasTaxonomy) {
                    for (int i = 0; i < otuNames.size(); i++) { otuTaxonomies[otuNames[i]] = results[1][i]; }
                }
            }
        }
        
        if (m->getControl_pressed()) {  return; }
        
        it = fileLines.find("columns");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a columns provided.\n"); }
        else {
            string thisLine = it->second;
            
            //read sample names
            maxLevel = 0;
            bool hasTaxonomy = false;
            vector< vector<string> > results = extractTaxonomyData(thisLine, numCols, hasTaxonomy);
            groupNames = results[0];
            if (hasTaxonomy) {
                for (int i = 0; i < results[1].size(); i++) {
                    if (m->getControl_pressed()) { break; }
                    
                    string completeTax = util.addUnclassifieds(results[1][i], maxLevel, false);
                    groupTaxonomies[results[0][i]] = completeTax;
                }
            }
        }
        
        if (m->getControl_pressed()) {  return; }
        
        it = fileLines.find("shape");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a shape provided.\n"); }
        else {
            string thisLine = it->second;
            getDims(thisLine, shapeNumRows, shapeNumCols);
            
            //check shape
            if (shapeNumCols != numCols) { m->mothurOut("[ERROR]: shape indicates " + toString(shapeNumCols) + " columns, but I only read " + toString(numCols) + " columns.\n"); m->setControl_pressed(true); }
            
            if (shapeNumRows != numOTUs) { m->mothurOut("[ERROR]: shape indicates " + toString(shapeNumRows) + " rows, but I only read " + toString(numOTUs) + " rows.\n"); m->setControl_pressed(true); }
        }
        
        if (m->getControl_pressed()) {  return; }
        
        it = fileLines.find("data");
        if (it == fileLines.end()) { m->mothurOut("[ERROR]: you file does not have a data provided.\n"); }
        else {
            string thisLine = it->second;
            
            if (shared != NULL) { delete shared; }
           
            shared = extractOTUData(thisLine, groupNames, numOTUs);
            shared->setOTUNames(otuNames);
            m->mothurOut("\n"+shared->getLabel()+"\n");
            
            if (otuTaxonomies.size() != 0) {
                //sanity check
                if ((shared->getNumBins() == otuTaxonomies.size()) && (shared->getNumBins() == numOTUs)) {
                    
                    for (int i = 0; i < shared->getNumBins(); i++) {
                        if (m->getControl_pressed()) { break; }
                        
                        string thisOTUsTax = otuTaxonomies[otuNames[i]];
                        string newTax = util.addUnclassifieds(thisOTUsTax, maxLevel, false);
                        Taxonomy thisOTUsTaxonomy(otuNames[i], newTax, shared->getOTUTotal(i));
                        consTax.push_back(thisOTUsTaxonomy);
                    }
                }
            }
        }
    
    }
    catch(exception& e) {
        m->errorOut(e, "BiomSimple", "read");
        exit(1);
    }
}
//**********************************************************************************************************************
//designed for things like "type": "OTU table", returns type
string BiomSimple::getTag(string& line) {
    try {
        bool inQuotes = false;
        string tag = "";
        char c = '\"';
        
        for (int i = 0; i < line.length(); i++) {
            
            //you want to ignore any ; until you reach the next '
            if ((line[i] == c) && (!inQuotes)) {  inQuotes = true;  }
            else if ((line[i] == c) && (inQuotes)) {
                inQuotes= false;
                line = line.substr(i+1);
                return tag;
            }
            
            if (inQuotes) {  if (line[i] != c) { tag += line[i]; }  }
        }
        
        return tag;
    }
    catch(exception& e) {
        m->errorOut(e, "BiomSimple", "getTag");
        exit(1);
    }
}
//**********************************************************************************************************************
//readRows
vector< vector<string> > BiomSimple::extractTaxonomyData(string line, int& numOTUs, bool& hasTaxonomy) {
    try {
        /*"rows":[
         {"id":"Otu01", "metadata":{"taxonomy":["Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Porphyromonadaceae", "unclassified"], "bootstrap":[100, 100, 100, 100, 100, 100]}},
         {"id":"Otu02", "metadata":{"taxonomy":["Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Rikenellaceae", "Alistipes"], "bootstrap":[100, 100, 100, 100, 100, 100]}},
         ...
         
         "rows":[{"id": "k__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae", "metadata": null},
         {"id": "k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae", "metadata": null}
         ....
         
         make look like above
         
         
         ],*/
        
        vector< vector<string> > results; results.resize(2);
        int countOpenBrace = 0; int countClosedBrace = 0; int openParen = 0; int closeParen = 0;
        string nextRow = "";
        bool end = false; bool allBlank = true;
        
        for (int i = 0; i < line.length(); i++) {
            
            if (m->getControl_pressed()) { return results; }
            
            if (line[i] == '[')         { countOpenBrace++;     }
            else if (line[i] == ']')    { countClosedBrace++;   }
            else if (line[i] == '{')    { openParen++;          }
            else if (line[i] == '}')    { closeParen++;         }
            else if (openParen != 0)    { nextRow += line[i];   }  //you are reading the row info
            
            //you have reached the end of the rows info
            if ((countOpenBrace == countClosedBrace) && (countClosedBrace != 0)) { end = true; break; }
            if ((openParen == closeParen) && (closeParen != 0)) { //process row
                numOTUs++;
                
                vector<string> result = getNamesAndTaxonomies(nextRow);
                if (result.size() != 0) { results[0].push_back(result[0]); results[1].push_back(result[1]); if (result[1] != "") { allBlank = false; } }
                
                nextRow = ""; openParen = 0; closeParen = 0;
            }
        }
        
        if (allBlank) { hasTaxonomy = false; }
        else { hasTaxonomy = true; }
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "BiomSimple", "extractTaxonomyData");
        exit(1);
    }
}
//**********************************************************************************************************************
//items[0] = id, items[1] = taxonomy, if items[2] then thats the taxonomy bootstrap values
vector<string> BiomSimple::getNamesAndTaxonomies(string line) {
    try {
        /*"rows":[
         {"id":"Otu01", "metadata":{"taxonomy":["Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Porphyromonadaceae", "unclassified"], "bootstrap":[100, 100, 100, 100, 100, 100]}},
         {"id":"Otu02", "metadata":{"taxonomy":["Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Rikenellaceae", "Alistipes"], "bootstrap":[100, 100, 100, 100, 100, 100]}},
         ...
         
         "rows":[{"id": "k__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae", "metadata": null},
         {"id": "k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae", "metadata": null}
         ....
         
         make look like above
         
         
         ],*/
        
        vector<string> results;
        if (line == "") { return results; }
        
        int pos = line.find_first_of(',');
        if (pos == string::npos) { //some kind of error?? we expect at least metadata : null, just grab name
            results.push_back(getName(line)); results.push_back("");
        }else {
            string value;
            util.splitAtComma(value, line);  //value hold name portion ("id":"Otu01") line holds rest
            results.push_back(getName(value));
            
            string taxonomy = ""; string bootstrap = "";
            int pos = line.find("taxonomy");
            if (pos != string::npos) { //no taxonomy info given
                int pos2 = line.find("bootstrap");
                if (pos2 != string::npos) { //no taxonomy info given
                    taxonomy = line.substr(pos, (pos2-pos));
                    taxonomy = taxonomy.substr(0, taxonomy.find_last_of(','));
                    bootstrap = line.substr(pos2);
                }else {
                    taxonomy = line.substr(pos);
                }
            }
            
            results.push_back(getTaxonomy(taxonomy, bootstrap));
        }
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "BiomSimple", "getNamesAndTaxonomies");
        exit(1);
    }
}
//**********************************************************************************************************************
string BiomSimple::getName(string line) {
    try {
        vector<string> nameItems;
        util.splitAtChar(line, nameItems, ':'); //split part we want containing the ids
        string name = nameItems[1];
        
        //remove "" if needed
        int pos = name.find("\"");
        if (pos != string::npos) {
            string newName = "";
            for (int k = 0; k < name.length(); k++) {
                if (name[k] != '\"') { newName += name[k]; }
            }
            name = newName;
        }
        
        return name;
    }
    catch(exception& e) {
        m->errorOut(e, "BiomSimple", "getName");
        exit(1);
    }
}
//**********************************************************************************************************************
//"taxonomy":"Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Porphyromonadaceae", "unclassified",
//"bootstrap":100, 100, 100, 100, 100, 100
string BiomSimple::getTaxonomy(string taxonomy, string bootstrap) {
    try {
        vector<string> results;
        
        if (taxonomy != "") {
            vector<string> taxItems;
            util.splitAtChar(taxonomy, taxItems, ':'); //split part we want containing the ids
            string taxons = taxItems[1];
            
            string taxon;
            while((taxons.find_first_of(',') != -1)) {
                if (m->getControl_pressed()) {break;}
                util.splitAtComma(taxon, taxons);
                results.push_back(taxon);
            }
            if (!util.stringBlank(taxons)) { results.push_back(taxons); }
        }
        
        if (bootstrap != "") {
            vector<string> bootItems;
            util.splitAtChar(bootstrap, bootItems, ':'); //split part we want containing the ids
            string bootValues = bootItems[1];
            
            string bootValue;
            int i = 0;
            while((bootValues.find_first_of(',') != -1)) {
                if (m->getControl_pressed()) {break;}
                util.splitAtComma(bootValue, bootValues);
                results[i]+="("+bootValue+")";
                i++;
            }
            if (!util.stringBlank(bootValues)) { results[i]+="("+bootValues+")"; }
        }
        
        string result = "";
        for (int i = 0; i < results.size(); i++) {
            if (m->getControl_pressed()) {result = ""; break;}
            result += results[i] + ";";
        }
        
        if (results.size() > maxLevel) { maxLevel = results.size(); }
       
        return result;
    }
    catch(exception& e) {
        m->errorOut(e, "BiomSimple", "getTaxonomy");
        exit(1);
    }
}
//**********************************************************************************************************************
void BiomSimple::getDims(string line, int& shapeNumRows, int& shapeNumCols) {
    try {
        //get shape
        bool inBar = false;
        string num = "";
        
        for (int i = 0; i < line.length(); i++) {
            
            //you want to ignore any ; until you reach the next '
            if ((line[i] == '[') && (!inBar)) {  inBar = true; i++;  if (!(i < line.length())) { break; } }
            else if ((line[i] == ']') && (inBar)) {
                inBar= false;
                util.mothurConvert(num, shapeNumCols);
                break;
            }
            
            if (inBar) {
                if (line[i] == ',') {
                    util.mothurConvert(num, shapeNumRows);
                    num = "";
                }else { if (!isspace(line[i])) { num += line[i]; }  }
            }
        }
    }
    catch(exception& e) {
        m->errorOut(e, "BiomSimple", "getDims");
        exit(1);
    }
}
//**********************************************************************************************************************
//readData
SharedRAbundVectors* BiomSimple::extractOTUData(string line, vector<string>& groupNames, int numOTUs) {
    try {
        SharedRAbundVectors* lookup = new SharedRAbundVectors();
        
        //creates new sharedRAbunds
        for (int i = 0; i < groupNames.size(); i++) {
            SharedRAbundVector* temp = new SharedRAbundVector(numOTUs); //sets all abunds to 0
            temp->setLabel(label);
            temp->setGroup(groupNames[i]);
            lookup->push_back(temp);
        }
        
        bool dataStart = false;
        bool inBrackets = false;
        string num = "";
        vector<int> nums;
        int otuCount = 0;
        for (int i = 0; i < line.length(); i++) {
            
            if (m->getControl_pressed()) { return lookup; }
            
            //look for opening [ to indicate data is starting
            if ((line[i] == '[') && (!dataStart)) { dataStart = true; i++;  if (!(i < line.length())) { break; } }
            else if ((line[i] == ']') && dataStart && (!inBrackets)) { break; } //we are done reading data
            
            if (dataStart) {
                if ((line[i] == '[') && (!inBrackets)) { inBrackets = true; i++;  if (!(i < line.length())) { break; } }
                else if ((line[i] == ']') && (inBrackets)) {
                    inBrackets = false;
                    int temp;
                    float temp2;
                    if (matrixElementType == "float") { util.mothurConvert(num, temp2); temp = (int)temp2; }
                    else { util.mothurConvert(num, temp); }
                    nums.push_back(temp);
                    num = "";
                    
                    //save info to vectors
                    if (matrixFormat == "dense") {
                        
                        //sanity check
                        if (nums.size() != lookup->size()) { m->mothurOut("[ERROR]: trouble parsing OTU data.  OTU " + toString(otuCount) + " causing errors.\n"); m->setControl_pressed(true); }
                        
                        //set abundances for this otu
                        //nums contains [abundSample0, abundSample1, abundSample2, ...] for current OTU
                        for (int j = 0; j < groupNames.size(); j++) { lookup->set(otuCount, nums[j], groupNames[j]); }
                        
                        otuCount++;
                    }else {
                        //sanity check
                        if (nums.size() != 3) { m->mothurOut("[ERROR]: trouble parsing OTU data.\n"); m->setControl_pressed(true); }
                        
                        //nums contains [otuNum, sampleNum, abundance]
                        lookup->set(nums[0], nums[2], groupNames[nums[1]]);
                    }
                    nums.clear();
                }
                
                if (inBrackets) {
                    if (line[i] == ',') {
                        int temp;
                        util.mothurConvert(num, temp);
                        nums.push_back(temp);
                        num = "";
                    }else { if (!isspace(line[i])) { num += line[i]; }  }
                }
            }
        }
        
        return lookup;
    }
    catch(exception& e) {
        m->errorOut(e, "BiomSimple", "readData");
        exit(1);
    }
}
//**********************************************************************************************************************
void BiomSimple::printHeading(ofstream& out, string mothurVersionString, string sharedfile) {
    try {
        //string mothurString = "mothur" + toString(current->getVersion());
        
        time_t rawtime; struct tm * timeinfo;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        string dateString = asctime (timeinfo);
        int pos = dateString.find('\n');
        if (pos != string::npos) { dateString = dateString.substr(0, pos);}
        string spaces = "      ";
        
        //standard
        out << "{\n" + spaces + "\"id\":\"" + util.getSimpleName(sharedfile) + "-" + label + "\",\n" + spaces + "\"format\": \"" + version + "\",\n" + spaces + "\"format_url\": \"http://biom-format.org\",\n";
        out << spaces + "\"type\": \"OTU table\",\n" + spaces + "\"generated_by\": \"" << mothurVersionString << "\",\n" + spaces + "\"date\": \"" << dateString << "\",\n";
        
   
    }
    catch(exception& e) {
        m->errorOut(e, "BiomSimple", "printHeadings");
        exit(1);
    }
}
//**********************************************************************************************************************
void BiomSimple::print(ofstream& out) {
    try {
       /* vector<string> metadata = getMetaData(lookup);
        int numBins = lookup->getNumBins();
        
        if (m->getControl_pressed()) {  out.close(); return 0; }
        
        //get row info
        /*"rows":[
                {"id":"GG_OTU_1", "metadata":null},
                {"id":"GG_OTU_2", "metadata":null},
                {"id":"GG_OTU_3", "metadata":null},
                {"id":"GG_OTU_4", "metadata":null},
                {"id":"GG_OTU_5", "metadata":null}
                ],*/
       /* out << spaces + "\"rows\":[\n";
        string rowFront = spaces + spaces + "{\"id\":\"";
        string rowBack = "\", \"metadata\":";
        vector<string> currentLabels = lookup->getOTUNames();
        for (int i = 0; i < numBins-1; i++) {
            if (m->getControl_pressed()) {  out.close(); return 0; }
            out << rowFront << currentLabels[i] << rowBack << metadata[i] << "},\n";
        }
        out << rowFront << currentLabels[(numBins-1)] << rowBack << metadata[(numBins-1)] << "}\n" + spaces + "],\n";
       
        //get column info
        /*"columns": [
                    {"id":"Sample1", "metadata":null},
                    {"id":"Sample2", "metadata":null},
                    {"id":"Sample3", "metadata":null},
                    {"id":"Sample4", "metadata":null},
                    {"id":"Sample5", "metadata":null},
                    {"id":"Sample6", "metadata":null}
                    ],*/
        /*
        string colBack = "\", \"metadata\":";
        out << spaces + "\"columns\":[\n";
        vector<string> namesOfGroups = lookup->getNamesGroups();
        for (int i = 0; i < namesOfGroups.size()-1; i++) {
            if (m->getControl_pressed()) {  out.close(); return 0; }
            out << rowFront << namesOfGroups[i] << colBack << sampleMetadata[i] << "},\n";
        }
        out << rowFront << namesOfGroups[(namesOfGroups.size()-1)] << colBack << sampleMetadata[lookup->size()-1] << "}\n" + spaces + "],\n";
        
        out << spaces + "\"matrix_type\": \"" << format << "\",\n" + spaces + "\"matrix_element_type\": \"int\",\n";
        out <<  spaces + "\"shape\": [" << numBins << "," << lookup->size() << "],\n";
        out << spaces + "\"data\":  [";
        
        vector<string> dataRows;
        if (format == "sparse") {
            /*"data":[[0,2,1],
             [1,0,5],
             [1,1,1],
             [1,3,2],
             [1,4,3],
             [1,5,1],
             [2,2,1],
             [2,3,4],
             [2,4,2],
             [3,0,2],
             [3,1,1],
             [3,2,1],
             [3,5,1],
             [4,1,1],
             [4,2,1]
             ]*/
        /*    string output = "";
            for (int i = 0; i < lookup->getNumBins(); i++) {
                
                if (m->getControl_pressed()) { out.close(); return 0; }
                vector<int> binAbunds = lookup->getOTU(i);
                
                for (int j = 0; j < binAbunds.size(); j++) {
                    int abund = binAbunds[j];
                    string binInfo = "[" + toString(i) + "," + toString(j) + "," + toString(abund) + "]";
                    //only print non zero values
                    if (abund != 0) { dataRows.push_back(binInfo); }
                }
            }
        }else {
            
            /* "matrix_type": "dense",
             "matrix_element_type": "int",
             "shape": [5,6],
             "data":  [[0,0,1,0,0,0],
             [5,1,0,2,3,1],
             [0,0,1,4,2,0],
             [2,1,1,0,0,1],
             [0,1,1,0,0,0]]*/
            /*
            for (int i = 0; i < lookup->getNumBins(); i++) {
                
                if (m->getControl_pressed()) { out.close(); return 0; }
                
                string binInfo = "[";
                vector<int> binAbund = lookup->getOTU(i);
                for (int j = 0; j < binAbund.size()-1; j++) {  binInfo += toString(binAbund[j]) + ","; }
                binInfo += toString(binAbund[binAbund.size()-1]) + "]";
                dataRows.push_back(binInfo);
            }
        }
        
        for (int i = 0; i < dataRows.size()-1; i++) {
            out << dataRows[i] << ",\n" + spaces  + spaces;
        }
        out << dataRows[dataRows.size()-1] << "]\n";
        
        out << "}\n";*/
   
    }
    catch(exception& e) {
        m->errorOut(e, "BiomSimple", "print");
        exit(1);
    }
}
/**********************************************************************************************************************
vector<string> BiomSimple::getMetaData(SharedRAbundVectors*& lookup){
    try {
        vector<string> metadata;
        
        if (consTax.size() == 0) { for (int i = 0; i < lookup->getNumBins(); i++) {  metadata.push_back("null");  } }
        else {
            
            //should the labels be Otu001 or PhyloType001
            string firstBin = lookup->getOTUNames()[0];
            string binTag = "Otu";
            if ((firstBin.find("Otu")) == string::npos) { binTag = "PhyloType";  }
            
            //convert list file bin labels to shared file bin labels
            //parse tax strings
            //save in map
            map<string, string> labelTaxMap;
            string snumBins = toString(otuLabels.size());
            for (int i = 0; i < otuLabels.size(); i++) {
                
                if (m->getControl_pressed()) { return metadata; }
                
                //if there is a bin label use it otherwise make one
                if (util.isContainingOnlyDigits(otuLabels[i])) {
                    string binLabel = binTag;
                    string sbinNumber = otuLabels[i];
                    if (sbinNumber.length() < snumBins.length()) {
                        int diff = snumBins.length() - sbinNumber.length();
                        for (int h = 0; h < diff; h++) { binLabel += "0"; }
                    }
                    binLabel += sbinNumber;
                    binLabel = util.getSimpleLabel(binLabel);
                    labelTaxMap[binLabel] = taxs[i];
                }else {
                    map<string, string>::iterator it = labelTaxMap.find(util.getSimpleLabel(otuLabels[i]));
                    if (it == labelTaxMap.end()) {
                        labelTaxMap[util.getSimpleLabel(otuLabels[i])] = taxs[i];
                    }else {
                        m->mothurOut("[ERROR]: Cannot add OTULabel " +  otuLabels[i] + " because it's simple label " + util.getSimpleLabel(otuLabels[i]) + " has already been added and will result in downstream errors. Have you mixed mothur labels and non mothur labels? To make the files work well together and backwards compatible mothur treats 1, OTU01, OTU001, OTU0001 all the same. We do this by removing any non numeric characters and leading zeros. For eaxample: Otu000018 and OtuMY18 both map to 18.\n"); m->setControl_pressed(true);
                    }
                }
            }
            
            //sanity check for file issues - do you have the same number of bins in the shared and constaxonomy file
            if (lookup->getNumBins() != labelTaxMap.size()) {
                m->mothurOut("[ERROR]: Your constaxonomy file contains " + toString(labelTaxMap.size()) + " otus and your shared file contain " + toString(lookup->getNumBins()) + " otus, cannot continue.\n"); m->setControl_pressed(true); return metadata;
            }
            
            //merges OTUs classified to same gg otuid, sets otulabels to gg otuids, averages confidence scores of merged otus.  overwritting of otulabels is fine because constaxonomy only allows for one label to be processed.  If this assumption changes, could cause bug.
            if (picrust) {  getGreenGenesOTUIDs(lookup, labelTaxMap);  }
            
            //{"taxonomy":["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriaceae", "g__Escherichia", "s__"]}
            
            //traverse the binLabels forming the metadata strings and saving them
            //make sure to sanity check
            map<string, string>::iterator it;
            vector<string> currentLabels = lookup->getOTUNames();
            for (int i = 0; i < lookup->getNumBins(); i++) {
                
                if (m->getControl_pressed()) { return metadata; }
                
                it = labelTaxMap.find(util.getSimpleLabel(currentLabels[i]));
                
                if (it == labelTaxMap.end()) { m->mothurOut("[ERROR]: can't find taxonomy information for " + currentLabels[i] + ".\n"); m->setControl_pressed(true); }
                else {
                    vector<string> bootstrapValues;
                    string data = "{\"taxonomy\":[";
            
                    vector<string> scores;
                    vector<string> taxonomies = parseTax(it->second, scores);
                    
                    for (int j = 0; j < taxonomies.size()-1; j ++) { data += "\"" + taxonomies[j] + "\", "; }
                    data += "\"" + taxonomies[taxonomies.size()-1] + "\"]";
                    
                    //add bootstrap values if available
                    if (scores[0] != "null") {
                        data += ", \"bootstrap\":[";
                        
                        for (int j = 0; j < scores.size()-1; j ++) { data += scores[j] + ", "; }
                        data += scores[scores.size()-1] + "]";

                    }
                    data += "}";
                    
                    metadata.push_back(data);
                }
            }
        }
        
        return metadata;
        
    }
    catch(exception& e) {
        m->errorOut(e, "BiomSimple", "getMetadata");
        exit(1);
    }

}
/**************************************************************************************************/
