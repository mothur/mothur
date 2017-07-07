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

class OptiMatrix {
    
#ifdef UNIT_TEST
    friend class TestOptiMatrix;
    friend class FakeOptiMatrix;
#endif

    
public:
    OptiMatrix() { m = MothurOut::getInstance(); }
    OptiMatrix(string, string, double, bool); //distfile, distformat, cutoff, sim
    OptiMatrix(string, string, string, string, double, bool); //distfile, name or count, format, distformat, cutoff, sim
    ~OptiMatrix(){ }
    
    set<int> getCloseSeqs(int i) { return closeness[i]; }
    //int get(int i, int j) { return closeness[i][j]; }
    bool isClose(int, int);
    int getNumClose(int index) { return closeness[index].size(); }
    int getNumSeqs() { return closeness.size(); }
    int getNumSingletons() { return singletons.size(); }
    //map<int, string> getNameMap() { return nameMap; }
    string getName(int); //name from nameMap index
    ListVector* getListSingle();
    long int print(ostream&);
    int readFile(string, string, string, string, double, bool); //distfile, name or count, format, distformat, cutoff, sim
    
protected:
    
    vector< set<int> > closeness;  //closeness[0] contains indexes of seqs "close" to seq 0.
    vector<string> singletons;
    vector<string> nameMap;
    
    string distFile, namefile, countfile, format, distFormat;
    double cutoff;
    bool sim;
    
    MothurOut* m;
    
    string findDistFormat(string distFile);
    int readPhylip();
    int readColumn();
    map<string, int> readNames(string namefile, vector<string>&);
    
};


#endif /* defined(__Mothur__optimatrix__) */
