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
        
        //create shortcut file for reading next time
        ofstream out; util.openOutputFile(optiDBName, out);
        
        //output version
        out << "#" << version << endl;
        
        //output reference aligned length
        out << alignedLength << endl;
        
        //output number of seqs in reference
        out << reference.numSeqs << endl;
        
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
        
        //read numSeqs in reference
        int numSeqs = 0;
        optiDBFile >> numSeqs; util.gobble(optiDBFile);
        
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
//an indicator column must have at least 50% of the bases the same
vector<int> OptiDB::getIndicatorColumns(){
    try {
        
        vector<int> indicatorCols;
        
        for (int i = 0; i < reference.otuData.size(); i++) { //for each alignment location
            
            vector<char> thisColumn = reference.otuData[i];
            
            if (thisColumn.size() == 1) { } //all sequences are the same in this column, ignore
            else {
                vector<double> counts;
                counts.resize(numBases, 0.0);
                
                //find occurances of each base
                for (int j = 0; j < thisColumn.size(); j++) { counts[baseMap[thisColumn[j]]]++; }
                
                //find percentages
                for (int j = 0; j < counts.size(); j++) {
                    
                    counts[j] /= numSeqs;
                    
                    if ((counts[j] >= 0.50) && (counts[j] < 0.95)) {
                        indicatorCols.push_back(i);
                        j+=counts.size(); //break
                    }
                }
            }
        }
        return indicatorCols;
        
    }
    catch(exception& e) {
        m->errorOut(e, "OptiDB", "getIndicatorColumns");
        exit(1);
    }
}
/**************************************************************************************************/

