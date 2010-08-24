#ifndef MOTHUROUT_H
#define MOTHUROUT_H

/*
 *  mothurOut.h
 *  Mothur
 *
 *  Created by westcott on 2/25/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"

/***********************************************/

class MothurOut {
	
	public:
		static MothurOut* getInstance();
		void setFileName(string);
		
		void mothurOut(string);
		void mothurOutEndLine();
		void mothurOutJustToLog(string);
		void errorOut(exception&, string, string);
		void closeLog();
		string getDefaultPath() { return defaultPath; }
		void setDefaultPath(string);
		
		string getReleaseDate() { return releaseDate; }
		void setReleaseDate(string r) { releaseDate = r; }
		string getVersion() { return version; }
		void setVersion(string r) { version = r; }
		
		//functions from mothur.h
		//file operations
		vector<unsigned long int> divideFile(string, int&);
		vector<unsigned long int> setFilePosEachLine(string, int&);
		vector<unsigned long int> setFilePosFasta(string, int&);
		string sortFile(string, string);
		void appendFiles(string, string);
		int renameFile(string, string); //oldname, newname
		string getFullPathName(string);
		string hasPath(string);
		string getExtension(string);
		string getPathName(string);
		string getSimpleName(string);
		string getRootName(string);
		bool isBlank(string);
		int openOutputFile(string, ofstream&);
		int openOutputFileAppend(string, ofstream&);
		int openInputFile(string, ifstream&);
		int openInputFile(string, ifstream&, string); //no error given 
		string getline(ifstream&);
		string getline(istringstream&);
		void gobble(istream&);
		void gobble(istringstream&);
		
		//searchs and checks
		bool checkReleaseVersion(ifstream&, string);
		bool anyLabelsToProcess(string, set<string>&, string);
		bool inUsersGroups(vector<string>, vector<string>);
		bool inUsersGroups(string, vector<string>);
		void getNumSeqs(ifstream&, int&);
		int getNumSeqs(ifstream&);
		int getNumNames(string);
		bool isTrue(string);
	
		
		//string manipulation
		void splitAtEquals(string&, string&);
		void splitAtComma(string&, string&);	
		void splitAtComma(string&, vector<string>&);
		void splitAtDash(string&, set<int>&);
		void splitAtDash(string&, set<string>&);
		void splitAtDash(string&, vector<string>&);
		void splitAtChar(string&, vector<string>&, char);
		
		//math operation
		int factorial(int num);
		vector<vector<double> > binomial(int);
		float ceilDist(float, int);
		float roundDist(float, int);

		int control_pressed;
		bool executing;
		

	private:
		static MothurOut* _uniqueInstance;
		MothurOut( const MothurOut& ); // Disable copy constructor
		void operator=( const MothurOut& ); // Disable assignment operator
		MothurOut() { control_pressed = false; defaultPath=""; };
		~MothurOut();

		string logFileName;
		string defaultPath;
		string releaseDate, version;
		
		ofstream out;
		
		int mem_usage(double&, double&);

};
/***********************************************/

#endif

