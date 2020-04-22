//
//  optidb.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/26/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "optidb.hpp"
#include "onegapdist.h"
#include "ignoregaps.h"
#include "eachgapdist.h"
#include "eachgapignore.h"
#include "onegapdist.h"
#include "onegapignore.h"

/**************************************************************************************************/
OptiDB::OptiDB(double c, string calcMethod, bool countends) : Database(), cutoff(c), aligned(true), alignedLength(0) {
    
    if (countends) {
        if (calcMethod == "nogaps")         {    distCalculator = new ignoreGaps(cutoff);       }
        else if (calcMethod == "eachgap")   {    distCalculator = new eachGapDist(cutoff);      }
        else if (calcMethod == "onegap")    {    distCalculator = new oneGapDist(cutoff);       }
    }else {
        if (calcMethod == "nogaps")         {    distCalculator = new ignoreGaps(cutoff);                }
        else if (calcMethod == "eachgap")   {    distCalculator = new eachGapIgnoreTermGapDist(cutoff);  }
        else if (calcMethod == "onegap")    {    distCalculator = new oneGapIgnoreTermGapDist(cutoff);   }
    }
}
/**************************************************************************************************/
//adds otu with seq as only reference included
void OptiDB::addSequence(Sequence seq)  {
    try {
        classifierOTU thisOtu(seq.getAligned());
        reference.push_back(thisOtu);
        
        numSeqs++; //this is the number of otus
    }
    catch(exception& e) {
        m->errorOut(e, "OptiDB", "addSequence");
        exit(1);
    }
}
/**************************************************************************************************/
//adds otu with seqs in it
void OptiDB::addSequences(vector<Sequence> refs)  {
    try {
        
        if (refs.size() == 0) { return; } //sanity check
        else {
            vector<string> seqs;
            for (int i = 0; i < refs.size(); i++) { seqs.push_back(refs[i].getAligned()); }
            
            classifierOTU thisOtu(seqs);
            reference.push_back(thisOtu);
            
            numSeqs++; //this is the number of otus
            
            if (thisOtu.numSeqs == 0) { m->mothurOut("[ERROR]: mothur expects the reference for opti_classifier to be aligned, please correct.\n"); aligned = false; m->setControl_pressed(true); }
        }
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
OptiData* OptiDB::findClosestSequences(Sequence*, int n) const  {
    try {
        
    }
    catch(exception& e) {
        m->errorOut(e, "OptiDB", "findClosestSequences");
        exit(1);
    }
}
/**************************************************************************************************/

