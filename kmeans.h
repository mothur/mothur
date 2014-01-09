//
//  kmeans.h
//  Mothur
//
//  Created by SarahsWork on 12/4/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef Mothur_kmeans_h
#define Mothur_kmeans_h

#include "communitytype.h"

/**************************************************************************************************/

class KMeans : public CommunityTypeFinder {
    
public:
    KMeans(vector<vector<int> >, int);
    vector<double> calcSilhouettes(vector< vector< double> >);
    double calcCHIndex(vector< vector< double> >);
    
private:

    int findSecondClosest(vector<int>&, vector<vector<double> >&, map<int, int>);
    double calcScore(int sample, int partition, vector<vector<double> >&, map<int, int>);

};

/**************************************************************************************************/

#endif
