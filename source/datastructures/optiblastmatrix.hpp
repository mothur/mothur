//
//  optiblastmatrix.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/10/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef optiblastmatrix_hpp
#define optiblastmatrix_hpp

#include "optidata.hpp"


class OptiBlastMatrix : public OptiData {
    
public:
    
    OptiBlastMatrix(string, string, string, bool, double c, int l, float p, bool min);
    ~OptiBlastMatrix(){ }
    
    vector< set<int> > getBlastOverlap() { return blastOverlap; }
    string getOverlapName(int); //name from nameMap index
    
protected:
    vector<string> overlapNameMap;
    vector< set<int> > blastOverlap;  //empty unless reading a blast file.
    
    string distFile, namefile, countfile, format;
    bool sim, minWanted;
    float penalty;
    int length;
    
    int readBlast();
    int readBlastNames(map<string, int>& nameAssignment);
};


#endif /* optiblastmatrix_hpp */
