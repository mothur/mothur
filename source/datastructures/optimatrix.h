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
    
public:
    OptiMatrix(string, double, bool); //distfile, cutoff, sim
    OptiMatrix(string, string, string, double, bool); //distfile, name or count, format, cutoff, sim
    ~OptiMatrix(){ }
    
    vector<int> getCloseSeqs(int i) { return closeness[i]; }
    int get(int i, int j) { return closeness[i][j]; }
    bool isClose(int, int);
    int getNumClose(int index) { return closeness[index].size(); }
    int getNumSeqs() { return closeness.size(); }
    map<int, string> getNameMap() { return nameMap; }
    ListVector* getList();
    
private:
    
    vector< vector<int> > closeness;  //closeness[0] contains indexes of seqs "close" to seq 0.
    vector<int> singletons;
    map<int, string> nameMap;
    
    string distFile, namefile, countfile, format, distFormat;
    double cutoff;
    bool sim;
    
    MothurOut* m;
    ListVector* list;
    ListVector* singlelist;
    
    string findDistFormat(string distFile);
    int readPhylip();
    int readColumn();
};


#endif /* defined(__Mothur__optimatrix__) */
