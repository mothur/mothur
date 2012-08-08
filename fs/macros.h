//
//  macros.h
//  rrf-fs-prototype
//
//  Created by Abu Zaher Faridee on 5/28/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef rrf_fs_prototype_macros_h
#define rrf_fs_prototype_macros_h

#define DEBUG_LEVEL_1
#define DEBUG_LEVEL_2
//#define DEBUG_LEVEL_3
//#define DEBUG_LEVEL_4

#define DEBUGMSG_LOCATION (cout << "DEBUGMSG " << __PRETTY_FUNCTION__ << "\nDEBUGMSG " << __FILE__ <<  "#"  << __LINE__ << endl)
#define DEBUGMSG_FUNC (cout << __PRETTY_FUNCTION__ << endl \
        << "--------------------------------------------------------------------------------" << endl)
#define DEBUGMSG_VAR(X) (cout << "DEBUGMSG " << __PRETTY_FUNCTION__ << "\nDEBUGMSG " << #X << " -> " << X << endl << endl)
#define PRINT_VAR(X) (cout << #X << " -> " << X << endl)
#define PRINT_MSG(Y, X) (cout <<  Y << " " << #X << " -> " << X << endl << endl)
#define DEBUGMSG(X) (cout << X << endl)
#define NAME_VALUE_PAIR(VAR, OSTREAM) (OSTREAM << #VAR << " : " << VAR)

#define FEATURE_DISCARD_SD_THRESHOLD  0

using namespace std;

class OptimumFeatureSubsetSelector{
public:
  OptimumFeatureSubsetSelector(string selectionType = "log2"): selectionType(selectionType){
  }
  
  int getOptimumFeatureSubsetSize(int numFeatures){
#ifdef DEBUG_MODE
    DEBUGMSG_LOCATION;
#endif
    if (selectionType == "log2"){ return (int)ceil(log2(numFeatures)); }
    else if (selectionType == "squareRoot"){ return (int)ceil(sqrt(numFeatures)); } 
    return -1;
  }
private:
  string selectionType;
};

// function for calculating standard deviation
double getStandardDeviation(vector<int> featureVector){
  int sum = accumulate(featureVector.begin(), featureVector.end(), 0);
    //    int zeroCount = count(featureVectors[i].begin(), featureVectors[i].end(), 0);
  double mean = (double) sum / (double) featureVector.size();
  vector<double> differenceFromMean(featureVector.size());
  transform(featureVector.begin(), featureVector.end(), differenceFromMean.begin(), bind2nd(minus<double>(), mean));
  double squaredSum = inner_product(differenceFromMean.begin(), differenceFromMean.end(), differenceFromMean.begin(), 0.0);
  double standardDeviation = sqrt(squaredSum / (double)featureVector.size());    
  return standardDeviation;
}


/* overrding "cout <<" for vector of integers */
ostream& operator <<(ostream& os, vector<int>& integers){
  os << "[ ";
  for (int i = 0; i < integers.size(); i++) {
    os << integers[i] << " ";
  }
  os << "]";
  return os;
}

/* overrding "cout <<" 2d matrix of itergers whuch uses vectors */
ostream& operator <<(ostream& os, vector< vector<int> > matrix){
//  os << "[ " << endl;
  os << "[ ";
  for (int i = 0; i < matrix.size(); i++) {
//    os << "\t" << i << " : " << matrix[i] << endl;
    os << i << " : " << matrix[i] << ", ";
  }
  os << "]";
  return os;
}

/* overrding "cout <<" for vector of booleans */
ostream& operator <<(ostream& os, vector<bool>& booleans){
  os << "[ ";
  for (int i = 0; i < booleans.size(); i++) {
    os << booleans[i] << " ";
  }
  os << "]";
  return os;
}

/* overrding "cout <<" for vector of double */
ostream& operator <<(ostream& os, vector<double>& doubles){
  os << "[ ";
  for (int i = 0; i < doubles.size(); i++) {
    os << doubles[i] << " ";
  }
  os << "]";
  return os;
}

/* overrding "cout <<" for map of int and vector<int> */
ostream& operator <<(ostream& os, map<int, vector<int> >& keyValuePairs){
//  os << "[ " << endl;
  os << "[ ";
  for (map<int, vector<int> >::iterator it = keyValuePairs.begin(); it != keyValuePairs.end(); it++) {
//    os << "\t" << it->first << " => " << it->second << endl;
    os << it->first << " => " << it->second << ", ";
  }
  os << "]";
  return os;
}

/* overrding "cout <<" for map of int and int */  
ostream& operator <<(ostream& os, map<int, int>& keyValuePairs){
  os << "[ " << endl;
  for (map<int, int>::iterator it = keyValuePairs.begin(); it != keyValuePairs.end(); it++) {
    os << "\t" << it->first << " => " << it->second << endl;
  }
  os << "]";
  return os;
}
  
  
#endif
