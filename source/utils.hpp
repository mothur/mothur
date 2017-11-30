//
//  utils.hpp
//  Mothur
//
//  Created by Sarah Westcott on 11/13/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef utils_hpp
#define utils_hpp

#include "mothurout.h"

class OrderVector;
class SharedOrderVector;
class RAbundVector;
class SharedRAbundVector;

class Utils {
    
public:
    
    Utils(); 
    ~Utils() {}
    
    //random operations
    int getRandomIndex(int); //highest
    int getRandomNumber();
    double getRandomDouble0to1();
    void mothurRandomShuffle(vector<int>&);
    void mothurRandomShuffle(vector< vector<double> >&);
    void mothurRandomShuffle(vector<string>&);
    void mothurRandomShuffle(vector<item>&);
    void mothurRandomShuffle(vector<PCell*>&);
    void mothurRandomShuffle(vector<PDistCellMin>&);
    void mothurRandomShuffle(OrderVector&);
    void mothurRandomShuffle(SharedOrderVector&);
    void mothurRandomShuffle(vector<SharedRAbundVector*>&);
    
    //checks
    bool isTrue(string);
    bool isContainingOnlyDigits(string);
    bool containsAlphas(string);
    bool isAllAlphas(string);
    bool isAllAlphaNumerics(string);
    bool isNumeric1(string);
    bool isNumeric1(char);
    bool isInteger(string);
    bool allSpaces(string);
    bool isLabelEquivalent(string, string);
    unsigned long long getRAMUsed();
    unsigned long long getTotalRAM();
    
    //file operations
    bool anyLabelsToProcess(string, set<string>&, string);
    bool appendBinaryFiles(string, string);
    int appendFiles(string, string);
    int appendFilesWithoutHeaders(string, string);
    vector<bool> allGZFiles(vector<string>&);
    bool appendSFFFiles(string, string);
    bool checkLocations(string&, vector<string>);  //filename, inputDir, outputDir. checks for file in ./, inputdir, outputdir, default and mothur's exe location.  Returns false if cant be found. If found completes name with location
    bool dirCheck(string&); //completes path, appends appropriate / or \, makes sure dir is writable.
    bool dirCheck(string&, string); //completes path, appends appropriate / or \, makes sure dir is writable. - no error
    vector<unsigned long long> divideFile(string, int&); //divides splitting unevenness by sequence
    vector<unsigned long long> divideFile(string filename, int& proc, char delimChar);
    vector<unsigned long long> divideFilePerLine(string, int&); //divides splitting unevenness at line breaks
    int divideFile(string, int&, vector<string>&);
    bool fileExists(string name);
    string findProgramPath(string programName);
    string getExtension(string);
    string getFullPathName(string);
    string getPathName(string);
    string getRootName(string);
    string getSimpleName(string);
    int getTimeStamp(string filename);
    string hasPath(string);
    bool isBlank(string);
    
    vector<bool> isGZ(string); //checks existence and format - will fail for either or both.
    bool mkDir(string&); //completes path, appends appropriate / or \. //returns true it exits or if we can make it
    bool mothurRemove(string);
    bool openInputFile(string, ifstream&, string); //no error given
    bool openOutputFile(string, ofstream&);
    bool openOutputFileAppend(string, ofstream&);
    bool openOutputFileBinary(string, ofstream&);
    bool openOutputFileBinaryAppend(string, ofstream&);
    bool openInputFile(string, ifstream&);
    bool openInputFileBinary(string, ifstream&);
    bool openInputFileBinary(string, ifstream&, string);
#ifdef USE_BOOST
    bool openInputFileBinary(string, ifstream&, boost::iostreams::filtering_istream&);
    bool openInputFileBinary(string, ifstream&, boost::iostreams::filtering_istream&, string);
#endif
    
    int printVsearchFile(vector<seqPriorityNode>&, string, string, string); //sorts and prints by abundance adding /ab=xxx/
    int renameFile(string, string); //oldname, newname
    vector<unsigned long long> setFilePosEachLine(string, long long&);
    vector<unsigned long long> setFilePosEachLine(string, unsigned long long&);
    vector<unsigned long long> setFilePosFasta(string, long long&);
    vector<unsigned long long> setFilePosFasta(string, long long&, char);
    string sortFile(string, string);
    
    //file reads
    bool checkReleaseVersion(ifstream&, string);
    string getline(ifstream&);
    string getline(istringstream&);
    void getNumSeqs(ifstream&, int&);
    int getNumSeqs(ifstream&);
    void gobble(istream&);
    void gobble(istringstream&);
    set<string> readAccnos(string);
    int readAccnos(string, vector<string>&);
    int readAccnos(string, vector<string>&, string);
    vector<consTax> readConsTax(string);
    int readConsTax(string, map<int, consTax2>&);
    map<string, int> readNames(string);
    map<string, int> readNames(string, unsigned long int&);
    int readNames(string, map<string, string>&, map<string, int>&);
    int readNames(string, map<string, string>&);
    int readNames(string, map<string, string>&, bool);
    int readNames(string, map<string, string>&, int);
    int readNames(string, map<string, vector<string> >&);
    int readNames(string, vector<seqPriorityNode>&, map<string, string>&);
    int readTax(string, map<string, string>&, bool);
    void zapGremlins(istream&);
    void zapGremlins(istringstream&);
    

    //math
    int average(vector<int>);
    vector<vector<double> > binomial(int);
    float ceilDist(float, int);
    int factorial(int num);
    unsigned int fromBase36(string);
    vector<double> getAverages(vector< vector<double> >&);
    double getAverage(vector<double>);
    vector< vector<seqDist> > getAverages(vector< vector< vector<seqDist> > >&, string);
    vector< vector<seqDist> > getAverages(vector< vector< vector<seqDist> > >&);
    double getStandardDeviation(vector<int>&);
    vector<double> getStandardDeviation(vector< vector<double> >&);
    vector<double> getStandardDeviation(vector< vector<double> >&, vector<double>&);
    vector< vector<seqDist> > getStandardDeviation(vector< vector< vector<seqDist> > >&);
    vector< vector<seqDist> > getStandardDeviation(vector< vector< vector<seqDist> > >&, vector< vector<seqDist> >&);
    double median(vector<double>);
    int median(vector<int>);
    float roundDist(float, int);
    int sum(vector<int> v) { return (accumulate(v.begin(), v.end(), 0)); }
    double sum(vector<double> v ) { return (accumulate(v.begin(), v.end(), 0.0)); }
    int max(vector<int> v) {  int max = 0; vector<int>::iterator it = max_element(v.begin(), v.end()); max = *it;  return max; }
    float max(vector<float> v) {  float max = 0; vector<float>::iterator it = max_element(v.begin(), v.end()); max = *it;  return max; }
    double max(vector<double> v) {  double max = 0; vector<double>::iterator it = max_element(v.begin(), v.end()); max = *it;  return max; }
    
    //type conversion
    bool mothurConvert(char, int&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
    bool mothurConvert(string, int&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
    bool mothurConvert(string, intDist&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
    bool mothurConvert(string, float&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
    bool mothurConvert(string, double&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
    bool mothurConvert(char, string&);

    
    
    //string manipulation
    string addUnclassifieds(string tax, int maxlevel, bool probs);
    bool checkGroupName(string name);
    int checkName(string&);
    void getCombos(vector<string>& groupComb, vector<string> userGroups, int& numComp);
    int getNumNames(string);
    int getNumChar(string, char);
    string getSimpleLabel(string);
    string getStringFromVector(vector<string>&, string); //creates string like "v[0], v[1], ... v[n]" where ', ' is string.
    string getStringFromVector(vector<int>&, string); //creates string like "v[0], v[1], ... v[n]" where ', ' is string.
    string getStringFromVector(vector<double>&, string); //creates string like "v[0], v[1], ... v[n]" where ', ' is string.
    bool inUsersGroups(vector<string>, vector<string>); //returns true if any of the strings in first vector are in second vector
    bool inUsersGroups(vector<int>, vector< vector<int> >);
    bool inUsersGroups(string, vector<string>);
    bool inUsersGroups(int, vector<int>);
    bool isSubset(vector<string>, vector<string>); //bigSet, subset
    string makeList(vector<string>&);
    map<string, vector<string> > parseClasses(string);
    int removeBlanks(vector<string>&);
    float removeConfidences(string&);
    string removeQuotes(string);
    bool stringBlank (string);
    void splitAtComma(string&, string&);
    void splitAtComma(string&, vector<string>&);
    void splitAtDash(string&, set<int>&);
    void splitAtDash(string&, set<string>&);
    void splitAtDash(string&, vector<string>&);
    void splitAtChar(string&, vector<string>&, char);
    void splitAtChar(string&, string&, char);
    void splitAtEquals(string&, string&);
    vector<string> splitWhiteSpaceWithQuotes(string);
    int splitWhiteSpaceWithQuotes(string, vector<string>&);
    vector<string> splitWhiteSpace(string& rest, char[], int);
    vector<string> splitWhiteSpace(string);
    int splitWhiteSpace(string, vector<float>&, int);
    int getOTUNames(vector<string>&, int, string);
    string getTag(string); //filename
    string findEdianness();

    
    
private:
    MothurOut* m;
    bool modifyNames;
    mt19937_64 mersenne_twister_engine;
    
};

#endif /* utils_hpp */
