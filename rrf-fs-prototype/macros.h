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

using namespace std;

class OptimumFeatureSubsetSelector{
public:
  OptimumFeatureSubsetSelector(string selectionType = "log2"): selectionType(selectionType){
  }
  
  int getOptimumFeatureSubsetSize(unsigned numFeatures){
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
  unsigned sum = accumulate(featureVector.begin(), featureVector.end(), 0);
    //    unsigned zeroCount = count(featureVectors[i].begin(), featureVectors[i].end(), 0);
  double mean = (double) sum / (double) featureVector.size();
  vector<double> differenceFromMean(featureVector.size());
  transform(featureVector.begin(), featureVector.end(), differenceFromMean.begin(), bind2nd(minus<double>(), mean));
  double squaredSum = inner_product(differenceFromMean.begin(), differenceFromMean.end(), differenceFromMean.begin(), 0.0);
  double standardDeviation = sqrt(squaredSum / featureVector.size());    
  return standardDeviation;
}


/* overrding "cout <<" for vector of integers */
ostream& operator <<(ostream& os, vector<int>& integers){
  os << "[ ";
  for (unsigned i = 0; i < integers.size(); i++) {
    os << integers[i] << " ";
  }
  os << "]";
  return os;
}

/* overrding "cout <<" 2d matrix of itergers whuch uses vectors */
ostream& operator <<(ostream& os, vector< vector<int> > matrix){
  os << "[ ";
  for (unsigned i = 0; i < matrix.size(); i++) {
    os << "ROW " << i << ":" << matrix[i] << endl;
  }
  os << "]";
  return os;
}

/* overrding "cout <<" for vector of booleans */
ostream& operator <<(ostream& os, vector<bool>& booleans){
  os << "[ ";
  for (unsigned i = 0; i < booleans.size(); i++) {
    os << booleans[i] << " ";
  }
  os << "]";
  return os;
}

/* overrding "cout <<" for vector of double */
ostream& operator <<(ostream& os, vector<double>& doubles){
  os << "[ ";
  for (unsigned i = 0; i < doubles.size(); i++) {
    os << doubles[i] << " ";
  }
  os << "]";
  return os;
}

#endif
