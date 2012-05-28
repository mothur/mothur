//
//  macros.h
//  rrf-fs-prototype
//
//  Created by Abu Zaher Faridee on 5/28/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef rrf_fs_prototype_macros_h
#define rrf_fs_prototype_macros_h

#define DEBUGMSG_LOCATION (cout << "DEBUGMSG " << __PRETTY_FUNCTION__ << "\nDEBUGMSG " << __FILE__ <<  "#"  << __LINE__ << endl)
#define DEBUGMSG_VAR(X) (cout << "DEBUGMSG " << __PRETTY_FUNCTION__ << "\nDEBUGMSG " <<#X << " -> " << X << endl << endl)

#endif
