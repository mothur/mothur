//
//  shannonrange.h
//  Mothur
//
//  Created by SarahsWork on 1/3/14.
//  Copyright (c) 2014 Schloss Lab. All rights reserved.
//

/*
 1] Haegeman, B., Hamelin, J., Moriarty, J., Neal, P., Dushoff, J., & Weitz, J. S. (2013). Robust estimation of microbial diversity in theory and in practice. The ISME journal, 7(6), 1092–1101.
 [2] Hill, M. O. (1973). Diversity and evenness: A unifying notation and its consequences. Ecology, 54(2), 427–432.
 [3] Orlitsky, A., Santhanam, N. P., & Zhang, J. (2003). Always Good Turing: Asymptoti- cally optimal probability estimation. Science, 302(5644), 427–431.
 [4] Roesch, L. F., Fulthorpe, R. R., Riva, A., Casella, G., Hadwin, A. K., Kent, A. D., et al. (2007). Pyrosequencing enumerates and contrasts soil microbial diversity. The ISME Journal, 1(4), 283–290.
 */

#ifndef Mothur_shannonrange_h
#define Mothur_shannonrange_h

#include "calculator.h"

/***********************************************************************/

class RangeShannon : public Calculator  {
	
public:
	RangeShannon(int a) : alpha(a), Calculator("rangeshannon", 3, false) {};
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(vector<RAbundVector*>) {return data;};
	string getCitation() { return "Haegeman, B., Hamelin, J., Moriarty, J., Neal, P., Dushoff, J., & Weitz, J. S. (2013). Robust estimation of microbial diversity in theory and in practice. The ISME journal, 7(6), 1092–1101., http://www.mothur.org/wiki/rangeshannon"; }
private:
    int alpha;
};

/***********************************************************************/



#endif
