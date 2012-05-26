//
//  sharedanddesignfilereader.h
//  rrf-fs-prototype
//
//  Created by Abu Zaher Faridee on 5/23/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef rrf_fs_prototype_sharedanddesignfilereader_h
#define rrf_fs_prototype_sharedanddesignfilereader_h

using namespace std;

class SharedAndDesignFileReader{
  
public:
  SharedAndDesignFileReader(string filePath){
    this->filePath = filePath;
    
  }
  void readFileContent(){
    inputFileStream.open(filePath.c_str());
    while (!inputFileStream.eof()) {
      getline(inputFileStream, lineBuffer);
      stringStream << lineBuffer;
      while (!stringStream.eof()) {
        string temp;
        getline(stringStream, temp, '\t');
        lineStrings.push_back(temp);
        temp.clear();
      }
      stringStream.clear();
      if (lineStrings.size() > 1){  // there was a succesful split, in case of 
                                    // the end of the file, there will be no 
                                    // split, we need to discard that
        fileContent.push_back(lineStrings);
      }
      lineStrings.clear();
    }
    inputFileStream.close();
  }
  
  void printFileContent(){
    for (unsigned i = 0; i < fileContent.size(); ++i) {
      vector<string> temp = fileContent[i];
      for (unsigned j = 0; j < temp.size(); ++j) {
        cout << temp[j] << "\t";
      }
      cout << endl << endl;
    }
  }
  
  vector< vector<string> > getFileContent(){ return fileContent; }
  
private:
  string filePath;
  ifstream inputFileStream;
  stringstream stringStream;
  string lineBuffer;
  vector<string> lineStrings; 
  vector< vector<string> > fileContent;
};


#endif
