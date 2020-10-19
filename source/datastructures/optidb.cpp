//
//  optidb.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/26/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "optidb.hpp"

/**************************************************************************************************/

OptiDB::OptiDB(string referenceFileName, string v) : Database() {
    alignedLength = 0;
    baseMap['A'] = 0;
    baseMap['T'] = 1;
    baseMap['G'] = 2;
    baseMap['C'] = 3;
    baseMap['-'] = 4;
    baseMap['N'] = 5;
    numBases = baseMap.size(); //A,T,G,C,-,N
    
    version = v;
    optiDBName = referenceFileName.substr(0,referenceFileName.find_last_of(".")+1) + "optidb";

}
/**************************************************************************************************/
vector< vector<int> > OptiDB::get(int i, char& allSame)  {
    try {
        vector< vector<int> > thisDistribution;
        allSame = 'x';
    
        vector<char> thisColumn;
        if (alignedLength < i) {
            m->mothurOut("[ERROR]: The reference alignment length is " + toString(alignedLength) + ", but you are requesting column " + toString(i) + ", please correct.\n"); m->setControl_pressed(true);
        }else {
            thisColumn = reference.otuData[i];
        }
        
        if (m->getControl_pressed()) { } //error in get column
        else {
            if (thisColumn.size() == 1) { //all sequences are the same in this column
                allSame = thisColumn[0];
            }else {
                thisDistribution.resize(numBases);
                for (int i = 0; i < thisColumn.size(); i++) {
                    thisDistribution[baseMap[thisColumn[i]]].push_back(i);
                }
            }
        }
        
        return thisDistribution;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiDB", "get");
        exit(1);
    }
}
/**************************************************************************************************/
//adds otu with seq as only reference included
void OptiDB::addSequence(Sequence seq)  {
    try {
        lengths.insert(seq.getAligned().length());
        refs.push_back(seq);
    }
    catch(exception& e) {
        m->errorOut(e, "OptiDB", "addSequence");
        exit(1);
    }
}
/**************************************************************************************************/
void OptiDB::convertSequences()  {
    try {
        
        if (refs.size() == 0) { return; } //sanity check
        else {

            vector<string> seqs;
            for (int i = 0; i < refs.size(); i++) {
                lengths.insert(refs[i].getAligned().length());
                
                //convert '.' gaps to '-'
                string aligned = refs[i].getAligned();
                for (int i = 0; i < aligned.length(); i++) { if (aligned[i] == '.') { aligned[i] = '-'; } }
                
                seqs.push_back(aligned);
                numSeqs++;
            }
            
            refs.clear();
            reference.readSeqs(seqs);
            
            if (reference.numSeqs == 0) { m->mothurOut("[ERROR]: mothur expects the reference for opti_classifier to be aligned, please correct.\n"); m->setControl_pressed(true); alignedLength = 0; }
        }
    }
    catch(exception& e) {
        m->errorOut(e, "OptiDB", "convertSequences");
        exit(1);
    }
}
/**************************************************************************************************/
void OptiDB::generateDB()  {
    try {
        //check to make sure actually aligned
        if (lengths.size() == 1) {  alignedLength = *lengths.begin(); longest = alignedLength-1;  } //database stores longest for aligner (longest = longest+1) so remove one.
        else {
            m->mothurOut("[ERROR]: mothur expects the reference for opti_classifier to be aligned, please correct.\n"); m->setControl_pressed(true); return;
        }
        
        //creates reference
        convertSequences();
         
        //finds columns in alignment with little noise that are included in query filter
        calcIndicatorColumns();
        
        //create shortcut file for reading next time
        ofstream out; util.openOutputFile(optiDBName, out);
        
        //output version
        out << "#" << version << endl;
        
        //output reference aligned length
        out << alignedLength << endl;
        
        //output number of seqs in reference
        out << reference.numSeqs << endl;
        
        out << util.getStringFromVector(indicatorColumns, ",") << endl;
        
        for (int i = 0; i < reference.otuData.size(); i++) { //for each alignment location
            out << i << '\t' << reference.otuData[i].size() << '\t';
            
            for (int j = 0; j < reference.otuData[i].size(); j++) { //for each reference, if all bases are the same in this location, size = 1; saves space
                out << reference.otuData[i][j];
            }
            out << endl;
        }
        out.close();
    }
    catch(exception& e) {
        m->errorOut(e, "OptiDB", "generateDB");
        exit(1);
    }
}
/**************************************************************************************************/

void OptiDB::readDB(ifstream& optiDBFile){
    try {
        optiDBFile.seekg(0);
        
        //read version
        string line = util.getline(optiDBFile); util.gobble(optiDBFile);
        
        //read alignedLength
        optiDBFile >> alignedLength; util.gobble(optiDBFile);
        
        longest = alignedLength-1;
        
        //read numSeqs in reference
        int numSeqs = 0;
        optiDBFile >> numSeqs; util.gobble(optiDBFile);
        
        line = util.getline(optiDBFile); util.gobble(optiDBFile);
        vector<string> iCols; util.splitAtComma(line, iCols);
        for (int i = 0; i < iCols.size(); i++) {
            int temp; util.mothurConvert(iCols[i], temp);
            //this column is significant
            indicatorColumns.push_back(temp);
        }
        
        vector<vector<char> > refDistrib; refDistrib.resize(alignedLength);
        int location, size; string bases;
        for (int i = 0; i < alignedLength; i++) { //for each alignment location
            
            optiDBFile >> location >> size >> bases; util.gobble(optiDBFile);
            
            char base;
            for (int j = 0; j < size; j++) { //for each reference, if all bases are the same in this location, size = 1; saves space
                
                refDistrib[location].push_back(bases[j]);
            }
        }
        optiDBFile.close();
        
        reference.readSeqs(refDistrib, numSeqs);
    }
    catch(exception& e) {
        m->errorOut(e, "OptiDB", "readDB");
        exit(1);
    }
}
/**************************************************************************************************/
//remove columns from list to process if they are not present in the query filter
map<int, int> OptiDB::filterIndicatorColumns(string filter, vector<int>& filteredICols){
    try {
        map<int, int> colsMap;
        
        filteredICols.clear();
        
        if (filter == "") { filter.resize(alignedLength, '1');  }
        
        //sanity check
        if (filter.length() != alignedLength) {  m->mothurOut("[ERROR]: Your filter indicates your alignment length is " + toString(filter.length()) + ", but your reference files indicate an alignment length of " + toString(alignedLength) + ". Cannot continue.\n");  m->setControl_pressed(true); return colsMap; }
        
        //process filter information
        map<int, int> colsPresentInQueryFiles; map<int, int>::iterator it;
        int filterCount = 0;
        for (int i = 0; i < alignedLength; i++) {
            if (filter[i] == '1') { //cols to keep
                colsPresentInQueryFiles[i] = filterCount;
                filterCount++;
            }
        }
        
        set<int> indicatorColsInTemplate;
        for (int i = 0; i < indicatorColumns.size(); i++) {
            indicatorColsInTemplate.insert(indicatorColumns[i]);
        }
        
        for (int i = 0; i < alignedLength; i++) {
            
            colsMap[i] = -1.0; //ignore col
            
            if (indicatorColsInTemplate.count(i) != 0) { //this is a template indicator column
                
                it = colsPresentInQueryFiles.find(i);
                if (it != colsPresentInQueryFiles.end()) { //this indicator column is present in the filtered query
                    filteredICols.push_back(it->second);
                    colsMap[i] = it->second; //use col
                }
            }
        }
        
        return colsMap;
        
    }
     catch(exception& e) {
         m->errorOut(e, "OptiDB", "filterIndicatorColumns");
         exit(1);
     }
 }
/**************************************************************************************************/
//an indicator column must have at least 50% of the bases the same
void OptiDB::calcIndicatorColumns(){
    try {
        
        for (int i = 0; i < reference.otuData.size(); i++) { //for each alignment location
            
            vector<char> thisColumn = reference.otuData[i];
            
            if (thisColumn.size() == 1) { } //all sequences are the same in this column, ignore
            else {
                vector<double> counts;
                counts.resize(numBases, 0.0);
                
                //find occurances of each base
                for (int j = 0; j < thisColumn.size(); j++) { counts[baseMap[thisColumn[j]]]++; }
                
                //find percentages
                for (int k = 0; k < counts.size(); k++) {
                    
                    counts[k] /= numSeqs;
                    
                    if ((counts[k] > 0.50) && (counts[k] < 0.95)) {
                        indicatorColumns.push_back(i);
                        break;
                    }
                }
            }
        }
    }
    catch(exception& e) {
        m->errorOut(e, "OptiDB", "calcIndicatorColumns");
        exit(1);
    }
}
/**************************************************************************************************/

