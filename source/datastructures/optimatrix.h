//
//  optimatrix.h
//  Mothur
//
//  Created by Sarah Westcott on 4/20/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__optimatrix__
#define __Mothur__optimatrix__

#include "mothurout.h"
#include "listvector.hpp"
#include "sparsedistancematrix.h"

class OptiMatrix {
    
#ifdef UNIT_TEST
    friend class TestOptiMatrix;
    friend class FakeOptiMatrix;
#endif

    
public:
    
    OptiMatrix() { m = MothurOut::getInstance(); }
    OptiMatrix(string, string, double, bool); //distfile, distformat, cutoff, sim
    OptiMatrix(string, string, string, string, double, bool); //distfile, name or count, format, distformat, cutoff, sim
    //OptiMatrix(SparseDistanceMatrix*, string, string); //sparse matrix, name or count, format
    ~OptiMatrix(){ }
    
    int readFile(string, string, string, string, double, bool); //distfile, name or count, format, distformat, cutoff, sim
    set<int> getCloseSeqs(int i) { return closeness[i]; }
    bool isClose(int, int);
    int getNumClose(int index) { return closeness[index].size(); }
    int getNumSeqs() { return closeness.size(); }
    int getNumSingletons() { return singletons.size(); }
    
    string getName(int); //name from nameMap index
    ListVector* getListSingle();
    long int print(ostream&);
    
protected:
    
    MothurOut* m;
    vector< set<int> > closeness;  //closeness[0] contains indexes of seqs "close" to seq 0.
    vector< set<int> > blastOverlap;  //empty unless reading a blast file.
    vector<string> singletons;
    vector<string> nameMap;
    vector<string> overlapNameMap;
    
    string distFile, namefile, countfile, format, distFormat;
    double cutoff;
    bool sim, minWanted;
    float penalty;
    int length;

    int readPhylip();
    int readColumn();
    int readBlast();
    int readBlastNames(map<string, int>& nameAssignment);
    
};


#endif /* defined(__Mothur__optimatrix__) */
