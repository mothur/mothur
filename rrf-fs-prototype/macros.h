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
#define DEBUGMSG_VAR(X) (cout << "DEBUGMSG " << __PRETTY_FUNCTION__ << "\nDEBUGMSG " << #X << " -> " << X << endl << endl)
#define NAME_VALUE_PAIR(VAR, OSTREAM) (OSTREAM << #VAR << " : " << VAR)

ostream& operator <<(ostream& os, vector<int>& integers){
  os << "[ ";
  for (unsigned i = 0; i < integers.size(); i++) {
    os << integers[i] << " ";
  }
  os << "]";
  return os;
}

ostream& operator <<(ostream& os, vector<bool>& booleans){
  os << "[ ";
  for (unsigned i = 0; i < booleans.size(); i++) {
    os << booleans[i] << " ";
  }
  os << "]";
  return os;
}


#endif
