//
//  optidb.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/26/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "optidb.hpp"

/**************************************************************************************************/
OptiDB::OptiDB(double c) : Database(), cutoff(c), aligned(true), alignedLength(0) {}
/**************************************************************************************************/
//adds otu with seq as only reference included
void OptiDB::addSequence(Sequence seq)  {
    try {
        
        vector<Sequence> seqs; seqs.push_back(seq);
        addSequences(seqs);
    }
    catch(exception& e) {
        m->errorOut(e, "OptiDB", "addSequence");
        exit(1);
    }
}
/**************************************************************************************************/
//adds otu with seqs in it
void OptiDB::addSequences(vector<Sequence> seqs)  {
    try {
        
        if (seqs.size() == 0) { return; } //sanity check
        else {
            if (alignedLength == 0) { //first otu, lets set it
                alignedLength = seqs[0].getAligned().length();
            }
        }
        
        vector<vector<char> > otuData;
        
        if (seqs.size() == 1) {
            string aligned = seqs[0].getAligned();
            
            if (alignedLength != aligned.length()) { aligned = false; }
            
            for (int i = 0; i < aligned.length(); i++) {
                vector<char> thisSpot;
                thisSpot.push_back(aligned[i]);
                
                otuData.push_back(thisSpot);
            }
        }else {
            otuData.resize(alignedLength);
            
            for (int i = 0; i < seqs.size(); i++) {
                string aligned = seqs[i].getAligned();
                
                if (alignedLength != aligned.length()) { aligned = false; break; }
                
                for (int j = 0; j < alignedLength; j++) { otuData[j].push_back(aligned[j]); }
            }
            
            //check for identical columns
            for (int i = 0; i < otuData.size(); i++) { //for each alignment column
                bool identical = true;
                char thisChar = otuData[i][0]; //set it first seq in otu
                
                for (int j = 1; j < otuData[i].size(); j++) { //for each seq in the otu
                    if (otuData[i][j] != thisChar) { identical = false; j += otuData[i].size(); }
                }
                
                //if all seqs are identical in this column, reduce to 1 to save space
                if (identical) { otuData[i].clear(); otuData[i].push_back(thisChar); }
            }
        }
        
        if (!aligned) { m->mothurOut("[ERROR]: mothur expects the reference for opti_classifier to be aligned, please correct.\n"); m->setControl_pressed(true); }
        
        classifierOTU thisOtu(otuData);
        reference.push_back(thisOtu);
    }
    catch(exception& e) {
        m->errorOut(e, "OptiDB", "addSequences");
        exit(1);
    }
}
/**************************************************************************************************/
void OptiDB::generateDB()  {
    try {
        
        
    }
    catch(exception& e) {
        m->errorOut(e, "OptiDB", "generateDB");
        exit(1);
    }
}
/**************************************************************************************************/

