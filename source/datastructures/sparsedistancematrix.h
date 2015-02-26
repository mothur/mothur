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


/* For each distance in a sparse matrix we have a row, column and distance.  
 The PDistCell consists of the column and distance.
 We know the row by the distances row in the seqVec matrix.  
 SeqVec is square and each row is sorted so the column values are ascending to save time in the search for the smallest distance. */

/***********************************************************************/
struct PDistCellMin{
	ull row;
    ull col;
	//PDistCell* cell;
	PDistCellMin(ull r, ull c) :  col(c), row(r) {}
};
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
