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

