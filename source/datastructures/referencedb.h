#ifndef MYREFERENCEDB_H
#define MYREFERENCEDB_H

/*
 *  referencedb.h
 *  Mothur
 *
 *  Created by westcott on 6/29/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */


#include "mothur.h"
#include "sequence.hpp"

/***********************************************/

class ReferenceDB {
	
	public:
	
		static ReferenceDB* getInstance();
		void clearMemory();
	
		bool save;
		vector<Sequence> referenceSeqs;
		vector< vector<float> > wordGenusProb;
		vector<diffPair> WordPairDiffArr;
	
		string getSavedReference()			{ return referencefile;		}
		void setSavedReference(string p)	{ referencefile = p;		}
		string getSavedTaxonomy()			{ return taxonomyfile;		}
		void setSavedTaxonomy(string p)		{ taxonomyfile = p;			}
	
	private:
	
		static ReferenceDB* myInstance;
		ReferenceDB() { referencefile = ""; taxonomyfile = ""; save = false; }
		ReferenceDB(const ReferenceDB&){}// Disable copy constructor
		void operator=(const ReferenceDB&){} // Disable assignment operator
		~ReferenceDB(){ myInstance = 0; }
	
		string referencefile, taxonomyfile;	
};
/***********************************************/

#endif

