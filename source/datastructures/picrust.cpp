//
//  picrust.cpp
//  Mothur
//
//  Created by Sarah Westcott on 11/16/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "picrust.hpp"

/**************************************************************************************************/

Picrust::Picrust(string ref, string otumapfile){
    try {
        m = MothurOut::getInstance();
        phyloTree = nullptr;
        
        read(ref, otumapfile);
        
    }
    catch(exception& e) {
        m->errorOut(e, "Picrust", "Picrust");
        exit(1);
    }
}
/**************************************************************************************************/

Picrust::Picrust(){
    try {
        m = MothurOut::getInstance();
        phyloTree = nullptr;
        
    }
    catch(exception& e) {
        m->errorOut(e, "Picrust", "Picrust");
        exit(1);
    }
}
/**************************************************************************************************/

Picrust::~Picrust(){
    try {
        if (phyloTree != nullptr) { delete phyloTree; }
        
    }
    catch(exception& e) {
        m->errorOut(e, "Picrust", "Picrust");
        exit(1);
    }
}
/**************************************************************************************************/

void Picrust::read(string ref, string otumapfile){
    try {
        
        //read reftaxonomy
        phyloTree = new PhyloTree(ref);
        
        //read otu map file
        readGGOtuMap(otumapfile); //maps reference ID -> OTU ID
        
    }
    catch(exception& e) {
        m->errorOut(e, "Picrust", "read");
        exit(1);
    }
}
//**********************************************************************************************************************
void Picrust::setGGOTUIDs(map<string, string>& labelTaxMap, SharedRAbundFloatVectors*& lookup){
    try {
        
        map<string, vector<string> > ggOTUIDs;
        
        //loop through otu taxonomies
        for (map<string, string>::iterator it = labelTaxMap.begin(); it != labelTaxMap.end(); it++) { //maps label -> consensus taxonomy
            if (m->getControl_pressed()) { break; }
            
            string OTUTaxonomy = it->second;
            
            //remove confidences
            util.removeConfidences(OTUTaxonomy);
            
            //remove unclassifieds to match template
            int thisPos = OTUTaxonomy.find("unclassified;"); //"Porphyromonadaceae"_unclassified;
            if (thisPos != string::npos) {
                OTUTaxonomy = OTUTaxonomy.substr(0, thisPos);
                thisPos = OTUTaxonomy.find_last_of(";"); //remove rest of parent taxon
                if (thisPos != string::npos) {
                    OTUTaxonomy = OTUTaxonomy.substr(0, thisPos+1);
                }
            }
            
            //get list of reference ids that map to this taxonomy
            vector<string> referenceIds = phyloTree->getSeqs(OTUTaxonomy);
            
            if (m->getControl_pressed()) { break; }
            
            //look for each one in otu map to find match
            string otuID = "not found";
            string referenceString = "";
            for (int i = 0; i < referenceIds.size(); i++) {
                referenceString += referenceIds[i] + " ";
                map<string, string>::iterator itMap = otuMap.find(referenceIds[i]);
                if (itMap != otuMap.end()) { //found it
                    otuID = itMap->second;
                    i += referenceIds.size(); //stop looking
                }
            }
            
            //if found, add otu to ggOTUID list
            if (otuID != "not found") {
                map<string, vector<string> >::iterator itGG = ggOTUIDs.find(otuID);
                if (itGG == ggOTUIDs.end()) {
                    vector<string> temp; temp.push_back(it->first); //save mothur OTU label
                    ggOTUIDs[otuID] = temp;
                }else { ggOTUIDs[otuID].push_back(it->first); } //add mothur OTU label to list
            }else {  m->mothurOut("[ERROR]: could not find OTUId for " + it->second + ". Its reference sequences are " + referenceString + ".\n"); m->setControl_pressed(true); }
            
        }
        
        vector<SharedRAbundFloatVector*> newLookup;
        vector<string> namesOfGroups = lookup->getNamesGroups();
        for (int i = 0; i < namesOfGroups.size(); i++) {
            SharedRAbundFloatVector* temp = new SharedRAbundFloatVector();
            temp->setLabel(lookup->getLabel());
            temp->setGroup(namesOfGroups[i]);
            newLookup.push_back(temp);
        }
        
        map<string, int> labelIndex;
        vector<string> currentLabels = lookup->getOTUNames();
        for (int i = 0; i < currentLabels.size(); i++) {  labelIndex[util.getSimpleLabel(currentLabels[i])] = i; }
        
        vector<string> newBinLabels;
        map<string, string> newLabelTaxMap;
        //loop through ggOTUID list combining mothur otus and adjusting labels
        //ggOTUIDs = 16097 -> <OTU01, OTU10, OTU22>
        
        for (map<string, vector<string> >::iterator itMap = ggOTUIDs.begin(); itMap != ggOTUIDs.end(); itMap++) {
            if (m->getControl_pressed()) { for (int j = 0; j < newLookup.size(); j++) {  delete newLookup[j];  } return; }
            
            //set new gg otu id to taxonomy. OTU01 -> k__Bacteria becomes 16097 -> k__Bacteria
            //find taxonomy of this otu
            map<string, string>::iterator it = labelTaxMap.find(util.getSimpleLabel(itMap->second[0]));
            vector<string> scores;
            vector<string> taxonomies = util.parseTax(it->second, scores);
            
            //merge/set OTU abundances
            vector<float> abunds; abunds.resize(lookup->size(), 0.0);
            string mergeString = "";
            vector<float> boots; boots.resize(scores.size(), 0.0);
            bool scoresnullptr = false;
            for (int j = 0; j < itMap->second.size(); j++) { //<OTU01, OTU10, OTU22>
                
                if (scores[0] != "null") {
                    //merge bootstrap scores
                    vector<string> scores;
                    vector<string> taxonomies = util.parseTax(it->second, scores);
                    for (int i = 0; i < boots.size(); i++) {
                        if (scores[i] == "null") { scoresnullptr = true; break; }
                        else {
                            float tempScore; util.mothurConvert(scores[i], tempScore);
                            boots[i] += tempScore;
                        }
                    }
                }else { scoresnullptr = true; }
                
                //merge abunds
                mergeString += (itMap->second)[j] + " ";
                for (int i = 0; i < lookup->size(); i++) { abunds[i] += lookup->get(labelIndex[util.getSimpleLabel((itMap->second)[j])], namesOfGroups[i]); }
            }
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: merging " + mergeString + " for ggOTUid = " + itMap->first + ".\n");  }
            
            //average scores
            //add merged otu to new lookup
            string newTaxString = "";
            if (!scoresnullptr) {
                for (int j = 0; j < boots.size(); j++) { boots[j] /= (float) itMap->second.size(); }
                
                //assemble new taxomoy
                for (int j = 0; j < boots.size(); j++) {
                    newTaxString += taxonomies[j] + "(" + toString(boots[j]) + ");";
                }
            }else {
                //assemble new taxomoy
                for (int j = 0; j < taxonomies.size(); j++) {
                    newTaxString += taxonomies[j] + ";";
                }
            }
            
            //set new gg otu id to taxonomy. OTU01 -> k__Bacteria becomes 16097 -> k__Bacteria
            //find taxonomy of this otu
            newLabelTaxMap[itMap->first] = newTaxString;
            
            //add merged otu to new lookup
            for (int j = 0; j < abunds.size(); j++) { newLookup[j]->push_back(abunds[j]); }
            
            //saved otu label
            newBinLabels.push_back(itMap->first);
        }
        
        lookup->clear();
        for (int i = 0; i < newLookup.size(); i++) { lookup->push_back(newLookup[i]);  }
        lookup->eliminateZeroOTUS();
        
        lookup->setOTUNames(newBinLabels);
        labelTaxMap = newLabelTaxMap;
        
        return;
    }
    catch(exception& e) {
        m->errorOut(e, "Picrust", "setGGOTUIDs");
        exit(1);
    }
    
}
//**********************************************************************************************************************
void Picrust::setGGOTUIDs(map<string, string>& labelTaxMap, SharedRAbundVectors*& lookup){
    try {
        
        map<string, vector<string> > ggOTUIDs;
        
        //loop through otu taxonomies
        for (map<string, string>::iterator it = labelTaxMap.begin(); it != labelTaxMap.end(); it++) { //maps label -> consensus taxonomy
            if (m->getControl_pressed()) { break; }
            
            string OTUTaxonomy = it->second;
            
            //remove confidences
            util.removeConfidences(OTUTaxonomy);
            
            //remove unclassifieds to match template
            int thisPos = OTUTaxonomy.find("unclassified;"); //"Porphyromonadaceae"_unclassified;
            if (thisPos != string::npos) {
                OTUTaxonomy = OTUTaxonomy.substr(0, thisPos);
                thisPos = OTUTaxonomy.find_last_of(";"); //remove rest of parent taxon
                if (thisPos != string::npos) {
                    OTUTaxonomy = OTUTaxonomy.substr(0, thisPos+1);
                }
            }
            
            //get list of reference ids that map to this taxonomy
            vector<string> referenceIds = phyloTree->getSeqs(OTUTaxonomy);
            
            if (m->getControl_pressed()) { break; }
            
            //look for each one in otu map to find match
            string otuID = "not found";
            string referenceString = "";
            for (int i = 0; i < referenceIds.size(); i++) {
                referenceString += referenceIds[i] + " ";
                map<string, string>::iterator itMap = otuMap.find(referenceIds[i]);
                if (itMap != otuMap.end()) { //found it
                    otuID = itMap->second;
                    i += referenceIds.size(); //stop looking
                }
            }
            
            //if found, add otu to ggOTUID list
            if (otuID != "not found") {
                map<string, vector<string> >::iterator itGG = ggOTUIDs.find(otuID);
                if (itGG == ggOTUIDs.end()) {
                    vector<string> temp; temp.push_back(it->first); //save mothur OTU label
                    ggOTUIDs[otuID] = temp;
                }else { ggOTUIDs[otuID].push_back(it->first); } //add mothur OTU label to list
            }else {  m->mothurOut("[ERROR]: could not find OTUId for " + it->second + ". Its reference sequences are " + referenceString + ".\n"); m->setControl_pressed(true); }
            
        }
        
        vector<SharedRAbundVector*> newLookup;
        vector<string> namesOfGroups = lookup->getNamesGroups();
        for (int i = 0; i < namesOfGroups.size(); i++) {
            SharedRAbundVector* temp = new SharedRAbundVector();
            temp->setLabel(lookup->getLabel());
            temp->setGroup(namesOfGroups[i]);
            newLookup.push_back(temp);
        }
        
        map<string, int> labelIndex;
        vector<string> currentLabels = lookup->getOTUNames();
        for (int i = 0; i < currentLabels.size(); i++) {  labelIndex[util.getSimpleLabel(currentLabels[i])] = i; }
        
        vector<string> newBinLabels;
        map<string, string> newLabelTaxMap;
        //loop through ggOTUID list combining mothur otus and adjusting labels
        //ggOTUIDs = 16097 -> <OTU01, OTU10, OTU22>
        
        for (map<string, vector<string> >::iterator itMap = ggOTUIDs.begin(); itMap != ggOTUIDs.end(); itMap++) {
            if (m->getControl_pressed()) { for (int j = 0; j < newLookup.size(); j++) {  delete newLookup[j];  } return; }
            
            //set new gg otu id to taxonomy. OTU01 -> k__Bacteria becomes 16097 -> k__Bacteria
            //find taxonomy of this otu
            map<string, string>::iterator it = labelTaxMap.find(util.getSimpleLabel(itMap->second[0]));
            vector<string> scores;
            vector<string> taxonomies = util.parseTax(it->second, scores);
            
            //merge/set OTU abundances
            vector<float> abunds; abunds.resize(lookup->size(), 0.0);
            string mergeString = "";
            vector<float> boots; boots.resize(scores.size(), 0.0);
            bool scoresnullptr = false;
            for (int j = 0; j < itMap->second.size(); j++) { //<OTU01, OTU10, OTU22>
                
                if (scores[0] != "null") {
                    //merge bootstrap scores
                    vector<string> scores;
                    vector<string> taxonomies = util.parseTax(it->second, scores);
                    for (int i = 0; i < boots.size(); i++) {
                        if (scores[i] == "null") { scoresnullptr = true; break; }
                        else {
                            float tempScore; util.mothurConvert(scores[i], tempScore);
                            boots[i] += tempScore;
                        }
                    }
                }else { scoresnullptr = true; }
                
                //merge abunds
                mergeString += (itMap->second)[j] + " ";
                for (int i = 0; i < lookup->size(); i++) { abunds[i] += lookup->get(labelIndex[util.getSimpleLabel((itMap->second)[j])], namesOfGroups[i]); }
            }
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: merging " + mergeString + " for ggOTUid = " + itMap->first + ".\n");  }
            
            //average scores
            //add merged otu to new lookup
            string newTaxString = "";
            if (!scoresnullptr) {
                for (int j = 0; j < boots.size(); j++) { boots[j] /= (float) itMap->second.size(); }
                
                //assemble new taxomoy
                for (int j = 0; j < boots.size(); j++) {
                    newTaxString += taxonomies[j] + "(" + toString(boots[j]) + ");";
                }
            }else {
                //assemble new taxomoy
                for (int j = 0; j < taxonomies.size(); j++) {
                    newTaxString += taxonomies[j] + ";";
                }
            }
            
            //set new gg otu id to taxonomy. OTU01 -> k__Bacteria becomes 16097 -> k__Bacteria
            //find taxonomy of this otu
            newLabelTaxMap[itMap->first] = newTaxString;
            
            //add merged otu to new lookup
            for (int j = 0; j < abunds.size(); j++) { newLookup[j]->push_back(abunds[j]); }
            
            //saved otu label
            newBinLabels.push_back(itMap->first);
        }
        
        lookup->clear();
        for (int i = 0; i < newLookup.size(); i++) { lookup->push_back(newLookup[i]);  }
        lookup->eliminateZeroOTUS();
        
        lookup->setOTUNames(newBinLabels);
        labelTaxMap = newLabelTaxMap;
        
        return;
    }
    catch(exception& e) {
        m->errorOut(e, "Picrust", "setGGOTUIDs");
        exit(1);
    }
    
}
//**********************************************************************************************************************
void Picrust::readGGOtuMap(string otumapfile){
    try {
        
        ifstream in; util.openInputFile(otumapfile, in);
        
        //map referenceIDs -> otuIDs
        //lines look like:
        //16097    671376    616121    533566    683683    4332909    4434717    772666    611808    695209
        while(!in.eof()) {
            if (m->getControl_pressed()) { break; }
            
            string line = util.getline(in); gobble(in);
            vector<string> pieces = util.splitWhiteSpace(line);
            
            if (pieces.size() != 0) {
                string otuID = pieces[1];
                for (int i = 1; i < pieces.size(); i++) {  otuMap[pieces[i]] = otuID; }
            }
        }
        in.close();
        
    }
    catch(exception& e) {
        m->errorOut(e, "Picrust", "readGGOtuMap");
        exit(1);
    }
    
}

/**************************************************************************************************/
