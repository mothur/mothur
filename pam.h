//
//  pam.h
//  Mothur
//
//  Created by SarahsWork on 12/10/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef Mothur_pam_h
#define Mothur_pam_h

#include "communitytype.h"

/**************************************************************************************************/

class Pam : public CommunityTypeFinder {
    
public:
    Pam(vector<vector<int> >, vector<vector<double> >, int);
    vector<double> calcSilhouettes(vector< vector< double> >);
    double calcCHIndex(vector< vector< double> >);
    
private:
    set<int> medoids;
    map<int, int> medoid2Partition;
    double largestDist;
    vector<vector<double> > dists;
    vector<vector< double> > Dp; // [numSamples][2] - It contains Dp and Ep. Dp is in [numSamples][0] and Ep is in [numSamples][1]. Dp is the distance between p and the closest sample in S and Ep is the distance between p and the second closest object in S. Both are used in the build and swap phases.
    
    int buildPhase();
    int swapPhase();
    int updateDp();
    
    
    
/**************************************************************************************************/
};


#endif
