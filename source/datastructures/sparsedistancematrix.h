#ifndef Mothur_sparsedistancematrix_h
#define Mothur_sparsedistancematrix_h

//
//  sparsedistancematrix.h
//  Mothur
//
//  Created by Sarah Westcott on 7/16/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "mothur.h"
#include "mothurout.h"


class ListVector;

/***********************************************************************/

class SparseDistanceMatrix {
	
public:
	SparseDistanceMatrix();
	~SparseDistanceMatrix(){ clear(); }
	int getNNodes();
	ull getSmallestCell(ull& index);		//Return the cell with the smallest distance
	float getSmallDist();
	
	int rmCell(ull, ull);
    int updateCellCompliment(ull, ull);
    void resize(ull n) { seqVec.resize(n);  }
    void clear();
	void addCell(ull, PDistCell);
    int addCellSorted(ull, PDistCell);
    vector<vector<PDistCell> > seqVec;
    
    
private:
	PDistCell smallCell;				//The cell with the smallest distance
	int numNodes;
    
    bool sorted;
    int sortSeqVec();
    int sortSeqVec(int);
	float smallDist, aboveCutoff;
    
	MothurOut* m;

};

/***********************************************************************/



#endif
