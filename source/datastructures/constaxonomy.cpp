//
//  constaxonomy.cpp
//  Mothur
//
//  Created by Sarah Westcott on 1/13/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "constaxonomy.hpp"

/***********************************************************************/
ConsTaxonommy::ConsTaxonommy(){
    m = MothurOut::getInstance();
}
/***********************************************************************/
ConsTaxonommy::ConsTaxonommy(string otuname, string consensusTax, int num) {
    try {
        m = MothurOut::getInstance();
        
        name = otuname;
        numReps = num;
        taxonomy = parseTax(consensusTax);
        
    }
    catch(exception& e) {
        m->errorOut(e, "ConsTaxonommy", "ConsTaxonommy");
        exit(1);
    }
}
/***********************************************************************/
ConsTaxonommy::ConsTaxonommy(ifstream& in) {
    try {
        m = MothurOut::getInstance();
        
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
        m->errorOut(e, "ConsTaxonommy", "ConsTaxonommy");
        exit(1);
    }
}
/***********************************************************************/
void ConsTaxonommy::setTaxons(string consensusTax){
    try {
        
        taxonomy = parseTax(consensusTax);
        
    }catch(exception& e) {
            m->errorOut(e, "ConsTaxonommy", "setTaxons");
            exit(1);
    }
}
/***********************************************************************/
string ConsTaxonommy::getInlineConsTaxonomy(){
    try {
        string otuConsensus = "";
        
        otuConsensus += name + '\t' + toString(numReps) + '\t' + getConsTaxString(true);
        
        return otuConsensus;
        
    }catch(exception& e) {
            m->errorOut(e, "ConsTaxonommy", "getInlineConsTaxonomy");
            exit(1);
    }
}
/***********************************************************************/
string ConsTaxonommy::getConsTaxString(bool includeConfidence) { //pass in true to include confidences
    try {
        
        string conTax = "";
        
        for (int i = 0; i < taxonomy.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            conTax += taxonomy[i].name;
            
            if (includeConfidence) { conTax += "(" + toString(taxonomy[i].confidence) + ")"; }
            
            conTax += ";";
        }
        
        return conTax;
        
    }catch(exception& e) {
            m->errorOut(e, "ConsTaxonommy", "getConsTaxString");
            exit(1);
    }
}
/***********************************************************************/
vector<Taxon> ConsTaxonommy::parseTax(string tax){
    try {
        string taxon = "";
        vector<Taxon> consTaxs;
        
        for(int i=0;i<tax.length();i++){
            
            if (m->getControl_pressed()) { break; }
            
            if(tax[i] == ';'){
                
                string newtaxon, confidence;
                
                int openParen = taxon.find_last_of('(');
                int closeParen = taxon.find_last_of(')');
                
                if ((openParen != string::npos) && (closeParen != string::npos)) {
                    string confidenceScore = taxon.substr(openParen+1, (closeParen-(openParen+1)));
                    if (util.isNumeric1(confidenceScore)) {  //its a confidence
                        newtaxon = taxon.substr(0, openParen); //rip off confidence
                        confidence = taxon.substr((openParen+1), (closeParen-openParen-1));
                    }else { //its part of the taxon
                        newtaxon = taxon;
                        confidence = "0";
                    }
                }else{
                    newtaxon = taxon;
                    confidence = "-1";
                }
                
                float con = 0; convert(confidence, con);
                
                Taxon thisTax(taxon, con);
                consTaxs.push_back(thisTax);
                
                taxon = "";
            }
            else{ taxon += tax[i]; }
        }
        
        return consTaxs;
       
    }catch(exception& e) {
            m->errorOut(e, "ConsTaxonommy", "parseTax");
            exit(1);
    }
}
/***********************************************************************/
void ConsTaxonommy::printConsTax(ostream& out){
    try {
        
        out << getInlineConsTaxonomy() << endl;
        
    }catch(exception& e) {
            m->errorOut(e, "ConsTaxonommy", "printConsTax");
            exit(1);
    }
}
/***********************************************************************/
void ConsTaxonommy::printConsTax(OutputWriter* out){
    try {
        
        out->write(getInlineConsTaxonomy()+"\n");
        
    }catch(exception& e) {
            m->errorOut(e, "ConsTaxonommy", "printConsTax");
            exit(1);
    }
}
/***********************************************************************/
void ConsTaxonommy::printConsTaxNoConfidence(ostream& out){
    try {
        
        string otuConsensus = name + '\t' + toString(numReps) + '\t' + getConsTaxString(false);
        
        out << otuConsensus << endl;
        
    }catch(exception& e) {
            m->errorOut(e, "ConsTaxonommy", "printConsTaxNoConfidence");
            exit(1);
    }
}
/***********************************************************************/
