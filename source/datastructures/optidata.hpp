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
    
#ifdef UNIT_TEST
    friend class TestOptiMatrix;
    friend class FakeOptiMatrix;
#endif
    
public:
    
    OptiData(double c)  { m = MothurOut::getInstance(); cutoff = c; }
    virtual ~OptiData(){}
    
    set<long long> getCloseSeqs(long long i);// { return closeness[i]; }
    bool isClose(long long, long long);
    long long getNumClose(long long);
    map<string, long long> getNameIndexMap();
    string getName(long long); //name from nameMap index
    set<string> getNames(set<long long>); //name from nameMap index
    
    long long getNumSeqs() { return closeness.size(); }
    long long getNumSingletons() { return singletons.size(); }
    virtual long long getNumDists(); //number of distances under cutoff
    ListVector* getListSingle();
    
    //for mgcluster - reading blast files
    virtual vector< set<long long> > getBlastOverlap() { vector< set<long long> > blank; return blank; }
    virtual string getOverlapName(long long) { return ""; } //name from nameMap index
    
    virtual void randomizeRefs() {}
    virtual vector<string> getRefSingletonNames() { vector<string> temp; return temp;  }
    virtual vector<long long> getTranslatedBins(vector<vector<string> >&, vector< vector<long long> >&) { vector<long long> temp; return temp;  }
    virtual OptiData* extractRefMatrix() { OptiData* temp = NULL; return temp;  }
    virtual OptiData* extractMatrixSubset(set<long long>&) { OptiData* temp = NULL; return temp;  }
    virtual OptiData* extractMatrixSubset(set<string>&) { OptiData* temp = NULL; return temp;  }
    virtual long long getNumFitSingletons() { return 0; } //user singletons
    
    virtual long long getNumFitDists() { return 0; } //user distances under cutoff
    virtual long long getNumRefDists() { return 0; } //ref distances under cutoff
    
    virtual ListVector* getFitListSingle() { ListVector* list = NULL; return list; }
    virtual long long getNumFitTrueSingletons() { return 0; }
    
    virtual vector<long long> getRefSeqs() { vector<long long> temp; return temp;  }
    virtual vector<long long> getFitSeqs() { vector<long long> temp; return temp;  }
    virtual long long getNumFitSeqs() { return 0; }
    virtual long long getNumFitClose(long long) { return 0;  }
    virtual long long getNumRefClose(long long) { return 0;  }
    virtual set<long long> getCloseFitSeqs(long long i) { set<long long> temp; return temp;  }
    virtual set<long long> getCloseRefSeqs(long long i) { set<long long> temp; return temp;  }
    virtual bool isCloseFit(long long j, long long i, bool&) { return false; }
    virtual long long print(ostream&);
    
protected:
    
    Utils util;
    MothurOut* m;
    vector< set<long long> > closeness;  //closeness[0] contains indexes of seqs "close" to seq 0.
    vector<string> singletons; //name of seqs with NO distances in matrix, if name file is given then it contains 2nd column of namefile
    vector<string> nameMap;  //name of seqs with distances in matrix, if name file is given then it contains 2nd column of namefile
    double cutoff;
    
    set<long long> getIndexes(set<string> seqs);
    
};


#endif /* optidata_hpp */
