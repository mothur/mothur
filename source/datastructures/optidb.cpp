//
//  optidb.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/26/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "optidb.hpp"

/**************************************************************************************************/
OptiDB::OptiDB() : Database() {
    alignedLength = 0;
    baseMap['A'] = 0;
    baseMap['T'] = 1;
    baseMap['G'] = 2;
    baseMap['C'] = 3;
    baseMap['-'] = 4;
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
vector<char> OptiDB::getColumn(int i)  {
    try {
        vector<char> thisColumn;
        
        if (alignedLength < i) {
            m->mothurOut("[ERROR]: The reference alignment length is " + toString(alignedLength) + ", but you are requesting column " + toString(i) + ", please correct.\n"); m->setControl_pressed(true);
        }else {
            thisColumn = reference.otuData[i];
        }
        
        return thisColumn;
    }
    catch(exception& e) {
        m->errorOut(e, "OptiDB", "getColumn");
        exit(1);
    }
}
/**************************************************************************************************/
vector< vector<int> > OptiDB::get(int i, char& allSame)  {
    try {
        vector< vector<int> > thisDistribution;
        allSame = 'x';
    
        vector<char> thisColumn = getColumn(i);
        
        if (m->getControl_pressed()) { } //error in get column
        else {
            if (thisColumn.size() == 1) { //all sequences are the same in this column
                allSame = thisColumn[0];
            }else {
                thisDistribution.resize(5);
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
void OptiDB::convertSequences()  {
    try {
        
        if (refs.size() == 0) { return; } //sanity check
        else {
            vector<string> seqs;
            for (int i = 0; i < refs.size(); i++) {
                lengths.insert(refs[i].getAligned().length());
                seqs.push_back(refs[i].getAligned());
                numSeqs++;
            }
            
            refs.clear();
            reference.readSeqs(seqs);
            
            if (reference.numSeqs == 0) { m->mothurOut("[ERROR]: mothur expects the reference for opti_classifier to be aligned, please correct.\n"); aligned = false; m->setControl_pressed(true); alignedLength = 0; }
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
            m->mothurOut("[ERROR]: mothur expects the reference for opti_classifier to be aligned, please correct.\n"); aligned = false; m->setControl_pressed(true);
        }
        
        convertSequences();
        //
    }
    catch(exception& e) {
        m->errorOut(e, "OptiDB", "generateDB");
        exit(1);
    }
}
/**************************************************************************************************/

void OptiDB::readDB(ifstream& optiDBFile){
    try {
        optiDBFile.seekg(0);                                    //    start at the beginning of the file
        
        //read version
        string line = util.getline(optiDBFile); util.gobble(optiDBFile);
        
        string seqName;
        int seqNumber;

        
        /*
         
         /***** TODO read shortcut file *******
         
         
        for(int i=0;i<maxKmer;i++){
            int numValues = 0;
            optiDBFile >> seqName >> numValues;
            
            for(int j=0;j<numValues;j++){                        //    for each kmer number get the...
                optiDBFile >> seqNumber;                        //        1. number of sequences with the kmer number
                kmerLocations[i].push_back(seqNumber);            //        2. sequence indices
            }
        }
         
         */
        optiDBFile.close();
        
    }
    catch(exception& e) {
        m->errorOut(e, "OptiDB", "readDB");
        exit(1);
    }
}
/**************************************************************************************************/

