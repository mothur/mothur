#ifndef SPLITMATRIX_H
#define SPLITMATRIX_H
/*
 *  splitmatrix.h
 *  Mothur
 *
 *  Created by westcott on 5/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "mothur.h"
#include "mothurout.h"

/******************************************************/

class SplitMatrix  {
	
	public:

		SplitMatrix(string, string, float); //column formatted distance file, namesfile, cutoff
		~SplitMatrix();
		int split();
		vector< map<string, string> > getDistanceFiles();  //returns map of distance files -> namefile sorted by distance file size
		string getSingletonNames() { return singleton; } //returns namesfile containing singletons
	
	private:
		string distFile, namefile, singleton;
		vector< map< string, string> > dists;
		float cutoff;
		
		MothurOut* m;
};

/******************************************************/

#endif

