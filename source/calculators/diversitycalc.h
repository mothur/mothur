//
//  diversitycalc.h
//  Mothur
//
//  Created by Sarah Westcott on 5/23/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef diversitycalc_h
#define diversitycalc_h

#include "mothurout.h"
#include "sabundvector.hpp"
#include "utils.hpp"




/***********************************************************************/
struct acceptRatioPos  {
    double acceptRatio;
    int pos;
    
    acceptRatioPos() { pos = 0; acceptRatio = 1.0; }
    acceptRatioPos(double ac, int po) : acceptRatio(ac), pos(po) {}
    ~acceptRatioPos() {}
};

/***********************************************************************/

inline bool operator< (const acceptRatioPos& lhs, const acceptRatioPos& rhs){ return rhs.acceptRatio > lhs.acceptRatio; }
inline bool operator> (const acceptRatioPos& lhs, const acceptRatioPos& rhs){ return rhs.acceptRatio < lhs.acceptRatio; }
inline bool operator<=(const acceptRatioPos& lhs, const acceptRatioPos& rhs){ return !(lhs.acceptRatio > rhs.acceptRatio); }
inline bool operator>=(const acceptRatioPos& lhs, const acceptRatioPos& rhs){ return !(lhs.acceptRatio < rhs.acceptRatio); }


/***********************************************************************/


class DiversityCalculator {
    
public:
    DiversityCalculator(bool rs){ m = MothurOut::getInstance();  requiresSamples = rs; }
    virtual ~DiversityCalculator(){};
    
    virtual string getTag() = 0;
    virtual bool requiresSample() { return requiresSamples; }
    virtual vector<double> getValues(int, vector<mcmcSample>&)  { return results;   }
    virtual vector<string> getValues(SAbundVector* rank)        { return outputs;   }
    virtual int getValues(SAbundVector* rank, vector<double>& ) { return 0;         }
    
protected:
    Utils util;
    MothurOut* m;
    
    bool requiresSamples;
    vector<double> results;
    vector<string> outputs;
    
    
};
/***********************************************************************/


#endif /* diversitycalc_h */
