/*
 *  m->mothurOut.cpp
 *  Mothur
 *
 *  Created by westcott on 2/25/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothurout.h"

/******************************************************/
MothurOut* MothurOut::getInstance() {
	if( _uniqueInstance == 0) {
		_uniqueInstance = new MothurOut();
	}
	return _uniqueInstance;
}
/*********************************************************************************************/
void MothurOut::setFileName(string filename)  {
	try {
		logFileName = filename;
		
		#ifdef USE_MPI
			int pid;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should output to screen
		#endif
		
		openOutputFile(filename, out);
		
		#ifdef USE_MPI
			}
		#endif
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "setFileName");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::closeLog()  {
	try {
		
		#ifdef USE_MPI
			int pid;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should output to screen
		#endif
		
		out.close();
		
		#ifdef USE_MPI
			}
		#endif
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "closeLog");
		exit(1);
	}
}

/*********************************************************************************************/
MothurOut::~MothurOut() {
	try {
		_uniqueInstance = 0;
		
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOut");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOut(string output) {
	try {
		
		#ifdef USE_MPI
			int pid;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should output to screen
		#endif
		
		cout << output;
		out << output;
		
		#ifdef USE_MPI
			}
		#endif
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOut");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOutEndLine() {
	try {
		#ifdef USE_MPI
			int pid;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should output to screen
		#endif
		
		cout << endl;
		out << endl;
		
		#ifdef USE_MPI
			}
		#endif
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOutEndLine");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::mothurOutJustToLog(string output) {
	try {
		#ifdef USE_MPI
			int pid;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
			if (pid == 0) { //only one process should output to screen
		#endif
		
		out << output;
		
		#ifdef USE_MPI
			}
		#endif
	}
	catch(exception& e) {
		errorOut(e, "MothurOut", "MothurOutJustToLog");
		exit(1);
	}
}
/*********************************************************************************************/
void MothurOut::errorOut(exception& e, string object, string function) {
	double vm, rss;
	mem_usage(vm, rss);
	
	mothurOut("Error: ");
	mothurOut(toString(e.what()));
	mothurOut(" has occurred in the " + object + " class function " + function + ". Please contact Pat Schloss at mothur.bugs@gmail.com, and be sure to include the mothur.logFile with your inquiry.");
	mothurOutEndLine();
}
/*********************************************************************************************/
//The following was originally from http://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-run-time-in-c 
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0
int MothurOut::mem_usage(double& vm_usage, double& resident_set) {
  #if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
  
	   vm_usage     = 0.0;
	   resident_set = 0.0;

	   // 'file' stat seems to give the most reliable results
	   //
	   ifstream stat_stream("/proc/self/stat",ios_base::in);

	   // dummy vars for leading entries in stat that we don't care about
	   //
	   string pid, comm, state, ppid, pgrp, session, tty_nr;
	   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	   string utime, stime, cutime, cstime, priority, nice;
	   string O, itrealvalue, starttime;

	   // the two fields we want
	   //
	   unsigned long vsize;
	   long rss;

	   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
				   >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
				   >> utime >> stime >> cutime >> cstime >> priority >> nice
				   >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

	   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	   vm_usage     = vsize / 1024.0;
	   resident_set = rss * page_size_kb;
	   
	   mothurOut("Memory Usage: vm = " + toString(vm_usage) + " rss = " + toString(resident_set) + "\n");
		return 0;

	#else
/*		//windows memory usage
		// Get the list of process identifiers.
		DWORD aProcesses[1024], cbNeeded, cProcesses;
		
		if ( !EnumProcesses( aProcesses, sizeof(aProcesses), &cbNeeded ) ){ return 1; }

		// Calculate how many process identifiers were returned.
		cProcesses = cbNeeded / sizeof(DWORD);

		// Print the memory usage for each process
		for (int i = 0; i < cProcesses; i++ ) {
			DWORD processID = aProcesses[i];
			
			PROCESS_MEMORY_COUNTERS pmc;

			HANDLE hProcess = OpenProcess((PROCESS_QUERY_INFORMATION | PROCESS_VM_READ), FALSE, processID);

			// Print the process identifier.
			printf( "\nProcess ID: %u\n", processID);
			
			if (NULL != hProcess) {

				if ( GetProcessMemoryInfo( hProcess, &pmc, sizeof(pmc)) ) {
					printf( "\tPageFaultCount: 0x%08X\n", pmc.PageFaultCount );
					printf( "\tPeakWorkingSetSize: 0x%08X\n", pmc.PeakWorkingSetSize );
					printf( "\tWorkingSetSize: 0x%08X\n", pmc.WorkingSetSize );
					printf( "\tQuotaPeakPagedPoolUsage: 0x%08X\n", pmc.QuotaPeakPagedPoolUsage );
					printf( "\tQuotaPagedPoolUsage: 0x%08X\n", pmc.QuotaPagedPoolUsage );
					printf( "\tQuotaPeakNonPagedPoolUsage: 0x%08X\n", pmc.QuotaPeakNonPagedPoolUsage );
					printf( "\tQuotaNonPagedPoolUsage: 0x%08X\n", pmc.QuotaNonPagedPoolUsage );
					printf( "\tPagefileUsage: 0x%08X\n", pmc.PagefileUsage ); 
					printf( "\tPeakPagefileUsage: 0x%08X\n", pmc.PeakPagefileUsage );
				}
				CloseHandle(hProcess);
			}
		}
*/
			return 0;

	#endif
}
/*********************************************************************************************/





