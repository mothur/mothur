//
//  optidata.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/10/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef optidata_hpp
#define optidata_hpp

#include "mothurout.h"
#include "listvector.hpp"
#include "sparsedistancematrix.h"
#include "counttable.h"

class OptiData {
    
public:
    
    OptiData(double c)  { m = MothurOut::getInstance(); cutoff = c; }
    virtual ~OptiData(){}
    
    set<int> getCloseSeqs(int i);// { return closeness[i]; }
    bool isClose(int, int);
    int getNumClose(int);
    map<string, int> getNameIndexMap();
    string getName(int); //name from nameMap index
    
    long long getNumSeqs() { return closeness.size(); }
    long long getNumSingletons() { return singletons.size(); }
    virtual long long getNumDists(); //number of distances under cutoff
    ListVector* getListSingle();
    
    //for mgcluster - reading blast files
    virtual vector< set<int> > getBlastOverlap() { vector< set<int> > blank; return blank; }
    virtual string getOverlapName(int) { return ""; } //name from nameMap index
    
    virtual vector<int> getTranslatedBins(vector<vector<string> >&, vector< vector<int> >&) { vector<int> temp; return temp;  }
    virtual OptiData* extractUnFitted(set<int>&) { OptiData* temp = NULL; return temp;  }
    virtual long long getNumFitSingletons() { return 0; } //user singletons
    virtual long long getNumRefSingletons() { return 0; } //reference singletons
    
    virtual long long getNumFitDists() { return 0; } //user distances under cutoff
    virtual long long getNumRefDists() { return 0; } //ref distances under cutoff
    
    //virtual ListVector* getRefListSingle();
    virtual ListVector* getFitListSingle() { ListVector* list = NULL; return list; }
    
    virtual vector<int> getRefSeqs() { vector<int> temp; return temp;  }
    virtual vector<int> getFitSeqs() { vector<int> temp; return temp;  }
    virtual long long getNumUniqueFitSeqs() { return 0; }
    virtual int getNumFitClose(int) { return 0;  }
    virtual int getNumRefClose(int) { return 0;  }
    virtual set<int> getCloseFitSeqs(int i) { set<int> temp; return temp;  }
    virtual set<int> getCloseRefSeqs(int i) { set<int> temp; return temp;  }
    virtual bool isCloseFit(int j, int i, bool&) { return false; }
    
protected:
    
    Utils util;
    MothurOut* m;
    vector< set<int> > closeness;  //closeness[0] contains indexes of seqs "close" to seq 0.
    vector<string> singletons; //name of seqs with NO distances in matrix, if name file is given then it contains 2nd column of namefile
    vector<string> nameMap;  //name of seqs with distances in matrix, if name file is given then it contains 2nd column of namefile
    double cutoff;
    
    long int print(ostream&);
    
};


#endif /* optidata_hpp */
