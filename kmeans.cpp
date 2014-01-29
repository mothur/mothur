//
//  kmeans.cpp
//  Mothur
//
//  Created by SarahsWork on 12/4/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "kmeans.h"

/**************************************************************************************************/

KMeans::KMeans(vector<vector<int> > cm, int p) : CommunityTypeFinder() {
    try {
        countMatrix = cm;
        numSamples = (int)countMatrix.size();
        numOTUs = (int)countMatrix[0].size();
        numPartitions = p;
        
        findkMeans();
    }
	catch(exception& e) {
		m->errorOut(e, "KMeans", "KMeans");
		exit(1);
	}
}
/**************************************************************************************************/




