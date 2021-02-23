//
//  constaxonomy.cpp
//  Mothur
//
//  Created by Sarah Westcott on 1/13/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "taxonomy.hpp"

/***********************************************************************/
Taxonomy::Taxonomy(){
    m = MothurOut::getInstance();
    containsConfidence = false;
}
/***********************************************************************/
Taxonomy::Taxonomy(string otuname, string consensusTax, int num) {
    try {
        m = MothurOut::getInstance();
        containsConfidence = false;
        
        name = otuname;
        numReps = num;
        taxonomy = parseTax(consensusTax);
        
    }
    catch(exception& e) {
        m->errorOut(e, "Taxonomy", "Taxonomy");
        exit(1);
    }
}
/***********************************************************************/
Taxonomy::Taxonomy(string otuname, string consensusTax) {
    try {
        m = MothurOut::getInstance();
        containsConfidence = false;
        
        name = otuname;
        numReps = 1;
        taxonomy = parseTax(consensusTax);
        
    }
    catch(exception& e) {
        m->errorOut(e, "Taxonomy", "Taxonomy");
        exit(1);
    }
}
/***********************************************************************/
Taxonomy::Taxonomy(ifstream& in) {
    try {
        m = MothurOut::getInstance();
        containsConfidence = false;
        
        string otu = ""; string consensusTax = "unknown";
        int size = 0;

        in >> otu; util.gobble(in);
        in >> size; util.gobble(in);
        consensusTax = util.getline(in); util.gobble(in);
        
        name = otu;
        numReps = size;
        taxonomy = parseTax(consensusTax);
        
    }
    catch(exception& e) {
        m->errorOut(e, "Taxonomy", "Taxonomy");
        exit(1);
    }
}
/***********************************************************************/
void Taxonomy::setTaxons(string consensusTax){
    try {
        
        taxonomy = parseTax(consensusTax);
        
    }catch(exception& e) {
            m->errorOut(e, "Taxonomy", "setTaxons");
            exit(1);
    }
}
/***********************************************************************/
string Taxonomy::getInlineConsTaxonomy(){
    try {
        string otuConsensus = "";
        
        otuConsensus += name + '\t' + toString(numReps) + '\t' + getConsTaxString(true);
        
        return otuConsensus;
        
    }catch(exception& e) {
            m->errorOut(e, "Taxonomy", "getInlineConsTaxonomy");
            exit(1);
    }
}
/***********************************************************************/
vector<string> Taxonomy::getSimpleTaxons(bool includeConfidence) { //pass in true to include confidences
    try {
        
        if (!containsConfidence) { includeConfidence = false; }
        vector<string> items;
        
        for (int i = 0; i < taxonomy.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            string conTax = taxonomy[i].name;
            
            if (includeConfidence) { conTax += "(" + toString(taxonomy[i].confidence) + ")"; }
            
            items.push_back(conTax);
        }
        
        return items;
        
    }catch(exception& e) {
            m->errorOut(e, "Taxonomy", "getSimpleTaxons");
            exit(1);
    }
}
/***********************************************************************/
string Taxonomy::getConsTaxString(bool includeConfidence) { //pass in true to include confidences
    try {
        
        string conTax = "";
        if (!containsConfidence) { includeConfidence = false; }
        
        for (int i = 0; i < taxonomy.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            conTax += taxonomy[i].name;
            
            if (includeConfidence) { conTax += "(" + toString(taxonomy[i].confidence) + ")"; }
            
            conTax += ";";
        }
        
        return conTax;
        
    }catch(exception& e) {
            m->errorOut(e, "Taxonomy", "getConsTaxString");
            exit(1);
    }
}
/***********************************************************************/
vector<Taxon> Taxonomy::parseTax(string tax){
    try {
        string taxon = "";
        vector<Taxon> consTaxs;
        
        for(int i=0;i<tax.length();i++){
            
            if (m->getControl_pressed()) { break; }
            
            if(tax[i] == ';'){
                
                string newtaxon = taxon; float confidence = 0;
                containsConfidence = util.hasConfidenceScore(newtaxon, confidence);
                
                Taxon thisTax(newtaxon, confidence);
                consTaxs.push_back(thisTax);
                
                taxon = "";
            }
            else{ taxon += tax[i]; }
        }
        
        return consTaxs;
       
    }catch(exception& e) {
            m->errorOut(e, "Taxonomy", "parseTax");
            exit(1);
    }
}
/***********************************************************************/
void Taxonomy::printConsTax(ostream& out){
    try {
        
        out << getInlineConsTaxonomy() << endl;
        
    }catch(exception& e) {
            m->errorOut(e, "Taxonomy", "printConsTax");
            exit(1);
    }
}
/***********************************************************************/
void Taxonomy::printConsTax(OutputWriter* out){
    try {
        
        out->write(getInlineConsTaxonomy()+"\n");
        
    }catch(exception& e) {
            m->errorOut(e, "Taxonomy", "printConsTax");
            exit(1);
    }
}
/***********************************************************************/
void Taxonomy::printConsTaxNoConfidence(ostream& out){
    try {
        
        string otuConsensus = name + '\t' + toString(numReps) + '\t' + getConsTaxString(false);
        
        out << otuConsensus << endl;
        
    }catch(exception& e) {
            m->errorOut(e, "Taxonomy", "printConsTaxNoConfidence");
            exit(1);
    }
}
/***********************************************************************/
