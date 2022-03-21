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
#include "utf8.h"

class OrderVector;
class SharedOrderVector;
class RAbundVector;
class SAbundVector;
class SharedRAbundVector;
class SharedRAbundVectors;
class SharedCLRVectors;
class Tree;
class PhyloTree;
class Taxonomy;
class InputData;
class ListVector;
class SharedRAbundFloatVectors;

class Utils {
    
public:
    
    Utils(); 
    ~Utils() {}
    
    //random operations
    int getRandomIndex(int); //highest
    long long getRandomIndex(long long); //highest
    int getRandomNumber();
    float randomUniform();
    float randomExp();
    float randomNorm();
    float randomGamma(float);
    vector<float> randomDirichlet(vector<float> alphas);
    double getRandomDouble0to1();
    void mothurRandomShuffle(vector<int>&);
    void mothurRandomShuffle(vector<long long>&);
    void mothurRandomShuffle(vector< vector<double> >&);
    void mothurRandomShuffle(vector<string>&);
    void mothurRandomShuffle(vector<intPair>&);
    void mothurRandomShuffle(vector<PCell*>&);
    void mothurRandomShuffle(vector<PDistCellMin>&);
    void mothurRandomShuffle(vector<colDist>&);
    void mothurRandomShuffle(OrderVector&);
    void mothurRandomShuffle(SharedOrderVector&);
    void mothurRandomShuffle(vector<SharedRAbundVector*>&);
    void mothurRandomShuffle(Tree* t, vector<string> g);
    void mothurRandomShuffle(vector<weightedSeq>&);
    
    //checks
    bool isTrue(string);
    bool isContainingOnlyDigits(string);
    bool containsAlphas(string);
    bool isASCII(string);
    bool isAllAlphas(string);
    bool isAllAlphaNumerics(string);
    bool isNumeric1(string);
    bool isNumeric1(char);
    bool isInteger(string);
    bool allSpaces(string);
    bool isLabelEquivalent(string, string);
    double getRAMUsed();
    double getTotalRAM();
    void getCurrentDate(string& thisYear, string& thisMonth, string& thisDay);
    
    //file operations
    bool mothurInitialPrep(vector<string>& defaultPath, vector<string>& tools, string& mothurVersion, string& releaseDate, string& OS);
    bool anyLabelsToProcess(string, set<string>&, string);
    bool appendBinaryFiles(string, string);
    int appendFiles(string first, string second); //first is appending to the end of second. 
    void appendFiles(string, ofstream&);
    int appendFilesFront(string, string);
    int appendFilesWithoutHeaders(string, string);
    vector<bool> allGZFiles(vector<string>&);
    bool appendSFFFiles(string, string);
    bool findTool(string& toolName, string&, vector<string>&, vector< vector<string> > locations);
    bool checkSpecificLocations(string&, vector<string>, string silent);
    bool checkLocations(string&, vector< vector<string> >, string silent);
    bool checkLocations(string&, vector< vector<string> >);  //filename, locations to check.  Returns false if cant be found. If found completes name with location
    bool checkLocationsGZ(string&, vector< vector<string> >);
    bool dirCheckWritable(string&); //completes path, appends appropriate / or \, makes sure dir is writable.
    bool dirCheckExists(string&);
    bool dirCheckExists(string&, bool); //completes path, appends appropriate / or \, makes sure dir is present.
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
    int getAlignmentLength(string);
    
    vector<bool> isGZ(string); //checks existence and format - will fail for either or both.
    bool isHDF5(string);
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
    bool openOutputFileBinary(string fileName, ofstream& file, ostream*& out, boost::iostreams::filtering_streambuf<boost::iostreams::output>& outBoost);
    string getline(boost::iostreams::filtering_istream& fileHandle);
#endif
    
    int printVsearchFile(vector<seqPriorityNode>&, string, string, string); //sorts and prints by abundance adding /ab=xxx/
    int renameFile(string, string); //oldname, newname
    int copyFile(string, string); //oldname, newname
    vector<double> setFilePosEachLine(string, long long&);
    vector<double> setFilePosEachLine(string, unsigned long long&);
    vector<double> setFilePosFasta(string, long long&);
    vector<double> setFilePosFasta(string, long long&, char);
    vector<double> divideFile(string, int&); //divides splitting unevenness by sequence
    vector<double> divideFile(string filename, int& proc, char delimChar);
    vector<double> divideFilePerLine(string, int&); //divides splitting unevenness at line breaks
    int divideFile(string, int&, vector<string>&);
    string sortFile(string, string);
    
    //file reads
    bool checkReleaseVersion(string, string);
    bool isVsearchVersionValid(string, string); //check to make sure minimum version requirements are made
    string getline(ifstream&);
    string getline(istringstream&);
    void getline(ifstream&, vector<string>&);
    void getNumSeqs(ifstream&, int&);
    int getNumSeqs(ifstream&);
    void gobble(istream&);
    void gobble(istringstream&);
    vector<string> parseTreeFile(string filename); //returns treenames
    set<string> readAccnos(string);
    int readAccnos(string, vector<string>&);
    int readAccnos(string, vector<string>&, string);
    void printAccnos(string, set<string>&);
    void printAccnos(string, vector<string>&);
    vector<consTax> readConsTax(string);
    void readConsTax(string, vector<Taxonomy>&);
    vector<Taxonomy> readConsTax(string inputfile, PhyloTree& tree); //fills tree
    int readConsTax(string, map<int, consTax2>&);
    void readNames(string, map<string, long long>&);
    map<string, int> readNames(string);
    map<string, int> readNames(string, unsigned long int&);
    int readNames(string, map<string, string>&, map<string, int>&);
    int readNames(string, map<string, string>&);
    int readNames(string, map<string, string>&, set<string>&);
    int readNames(string, map<string, string>&, bool);
    int readNames(string, map<string, string>&, int);
    int readNames(string, map<string, vector<string> >&);
    int readNames(string, vector<seqPriorityNode>&, map<string, string>&);
    int scanNames(string namefile); //return totalnum seqs
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
    float sum(vector<float> v ) { return (accumulate(v.begin(), v.end(), 0.0)); }
    int max(vector<int> v) {  int max = 0; vector<int>::iterator it = max_element(v.begin(), v.end()); max = *it;  return max; }
    float max(vector<float> v) {  float max = 0; vector<float>::iterator it = max_element(v.begin(), v.end()); max = *it;  return max; }
    double max(vector<double> v) {  double max = 0; vector<double>::iterator it = max_element(v.begin(), v.end()); max = *it;  return max; }
    long long max(vector<long long> v) {  long long max = 0; vector<long long>::iterator it = max_element(v.begin(), v.end()); max = *it;  return max; }
    bool isPositiveNumeric(string);
    bool isEqual(float, float); //handles approximate equal
    bool isEqual(double, double); //handles approximate equal
    double geometricMean(vector<float>&, double);
    
    //type conversion
    bool mothurConvert(char, int&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
    bool mothurConvert(string, int&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
    bool mothurConvert(string, intDist&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
    bool mothurConvert(string, float&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
    bool mothurConvert(string, double&); //use for converting user inputs. Sets commandInputsConvertError to true if error occurs. Engines check this.
    bool mothurConvert(char, string&);
    set<string> mothurConvert(vector<string>&);
    vector<string> mothurConvert(set<string>&);
    set<long long> mothurConvert(vector<long long>&);
    vector<long long > mothurConvert(set<long long>&);
    char* mothurConvert(string); //convert string to char*

    
    //string manipulation
    string addUnclassifieds(string tax, int maxlevel, bool probs);
    string trimTax(string tax, int trimLevel);
    bool checkGroupName(string& name);
    bool checkGroupNames(vector<string>& name);
    void checkName(string&);
    void getCombos(vector<string>& groupComb, vector<string> userGroups, int& numComp);
    int getNumNames(string);
    string getSimpleLabel(string);
    string toUpper(string);
    string toLower(string);
    
    string getStringFromVector(vector<string>&, string); //creates string like "v[0], v[1], ... v[n]" where ', ' is string.
    string getStringFromVector(vector<int>&, string); //creates string like "v[0], v[1], ... v[n]" where ', ' is string.
    string getStringFromVector(vector<double>&, string); //creates string like "v[0], v[1], ... v[n]" where ', ' is string.
    string getStringFromSet(set<int>&, string); //creates string like "v[0], v[1], ... v[n]" where ', ' is string.
    string getStringFromSet(set<string>&, string); //creates string like "v[0], v[1], ... v[n]" where ', ' is string.
    set<string> getSetFromList(ListVector*&, vector< vector<string> >&); 
    string getFormattedHelp(vector<string> question, vector<string> aquestion, vector<string> issue, vector<string> aissue, vector<string> howto,vector<string> ahowto);
    string trimStringEnd(string, int); //string, number of chars to remove from end.
    
    bool inUsersGroups(vector<string>, vector<string>); //returns true if any of the strings in first vector are in second vector
    bool inUsersGroups(vector<int>, vector< vector<int> >);
    bool inUsersGroups(string, vector<string>);
    bool inUsersGroups(string, set<string>);
    bool inUsersGroups(int, vector<int>);
    bool isSubset(vector<string>, vector<string>); //bigSet, subset
    
    string makeList(vector<string>&);
    map<string, vector<string> > parseClasses(string);
    int removeBlanks(vector<string>&);
    float removeConfidences(string&);
    string removeQuotes(string);
    void removeQuotes(vector<Taxon>& tax);
    vector<string> parseTax(string tax, vector<string>& scores);
    bool stringBlank (string);
    
    
    //file reading
    SharedRAbundVectors* getNextShared(InputData&, bool, set<string>&, set<string>&, string&, string optionalOutput = "");//input, allLines, userLabels, processedLabels, lastLabel
    SharedRAbundFloatVectors* getNextRelabund(InputData&, bool, set<string>&, set<string>&, string&);//input, allLines, userLabels, processedLabels, lastLabel
    SharedCLRVectors* getNextCLR(InputData&, bool, set<string>&, set<string>&, string&);//input, allLines, userLabels, processedLabels, lastLabel
    SharedOrderVector* getNextSharedOrder(InputData&, bool, set<string>&, set<string>&, string&);//input, allLines, userLabels, processedLabels, lastLabel
    ListVector* getNextList(InputData&, bool, set<string>&, set<string>&, string&);//input, allLines, userLabels, processedLabels, lastLabel
    RAbundVector* getNextRAbund(InputData&, bool, set<string>&, set<string>&, string&);//input, allLines, userLabels, processedLabels, lastLabel
    SAbundVector* getNextSAbund(InputData&, bool, set<string>&, set<string>&, string&);//input, allLines, userLabels, processedLabels, lastLabel
    OrderVector* getNextOrder(InputData&, bool, set<string>&, set<string>&, string&);//input, allLines, userLabels, processedLabels, lastLabel

    
    void splitAtComma(string&, string&);
    void splitAtComma(string&, vector<string>&);
    void splitAtComma(string&, vector<int>&);
    void splitAtDash(string&, set<int>&);
    void splitAtDash(string&, set<string>&);
    void splitAtDash(string&, vector<string>&);
    void splitAtChar(string&, set<string>&, char);
    void splitAtChar(string&, vector<string>&, char);
    void splitAtChar(string&, string&, char);
    void splitAtEquals(string&, string&);
    vector<string> splitWhiteSpaceWithQuotes(string);
    vector<string> splitWhiteSpace(string& rest, char[], int);
    vector<string> splitWhiteSpace(string);
    int splitWhiteSpace(string, vector<float>&, int);
    string trimWhiteSpace(string input);
    
    int getOTUNames(vector<string>&, int, string);
    string getTag(string); //filename
    string findEdianness();
    string removeNs(string seq);
    string reverseOligo(string);
    bool hasConfidenceScore(string&, float&); //taxon, confidence score. Returns taxon with confidence removed and confidence score.  If no confidence score, then confidence=0
    bool hasConfidenceScore(vector<Taxon> taxons); //returns true if taxons have positve confidences
    vector<Taxon> getTaxons(string, bool&);
    bool findTaxon(vector<Taxon> tax, vector<Taxon> searchTax);
    bool searchTax(vector<Taxon>, vector<bool> taxonsHasConfidence, vector< vector<Taxon> > searchTaxons);

private:
    MothurOut* m;
    bool modifyNames, modifyGroups;
    string homePath, currentWorkingDirectory;
    mt19937_64 mersenne_twister_engine;
    vector<string> paths; //paths stored in environment varaibale PATH
    
    vector<string> readTreeString(ifstream& filehandle);
};

#endif /* utils_hpp */


