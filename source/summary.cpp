//
//  summary.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/27/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "summary.hpp"


//**********************************************************************************************************************
Summary::Summary(string n) { //name or count file to include in counts
    try {
        m = MothurOut::getInstance(); total = 0; numUniques = 0; hasNameOrCount = true;
        
        if (m->isCountFile(n)) {
            CountTable ct;
            ct.readTable(n, false, false);
            nameMap = ct.getNameMap();
            type = "count";
        }else { nameMap = m->readNames(n); type = "name"; }
        
        nameCountNumUniques = nameMap.size();
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "Summary");
        exit(1);
    }
}
//**********************************************************************************************************************
int Summary::addSeq(Sequence seq) {
    try {
        long long num = 1;
        if (hasNameOrCount) {
            //make sure this sequence is in the namefile, else error
            itFindName = nameMap.find(seq.getName());
            
            if (itFindName == nameMap.end()) { m->mothurOut("[ERROR]: '" + seq.getName() + "' is not in your name or count file, please correct."); m->mothurOutEndLine(); m->control_pressed = true; }
            else { num = it->second; }
        }
        
        int thisStartPosition = seq.getStartPos();
        it = startPosition.find(thisStartPosition);
        if (it == startPosition.end()) { startPosition[thisStartPosition] = num; } //first finding of this start position, set count.
        else { it->second += num; } //add counts
        
        int thisEndPosition = seq.getEndPos();
        it = endPosition.find(thisEndPosition);
        if (it == endPosition.end()) { endPosition[thisEndPosition] = num; } //first finding of this end position, set count.
        else { it->second += num; } //add counts
        
        int thisSeqLength = seq.getNumBases();
        it = seqLength.find(thisSeqLength);
        if (it == seqLength.end()) { seqLength[thisSeqLength] = num; } //first finding of this length, set count.
        else { it->second += num; } //add counts
        
        int thisAmbig = seq.getAmbigBases();
        it = ambigBases.find(thisAmbig);
        if (it == ambigBases.end()) { ambigBases[thisAmbig] = num; } //first finding of this ambig, set count.
        else { it->second += num; } //add counts
        
        int thisHomoP = seq.getLongHomoPolymer();
        it = longHomoPolymer.find(thisHomoP);
        if (it == longHomoPolymer.end()) { longHomoPolymer[thisHomoP] = num; } //first finding of this homop, set count.
        else { it->second += num; } //add counts

        numUniques++;
        total += num;
        
        return num;
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "addSeq");
        exit(1);
    }
}
//**********************************************************************************************************************
vector<long long> Summary::getDefaults() {
    try {
        vector<long long> locations;
        
        long long ptile0_25	= 1+(long long)(total * 0.025); //number of sequences at 2.5%
        long long ptile25		= 1+(long long)(total * 0.250); //number of sequences at 25%
        long long ptile50		= 1+(long long)(total * 0.500);
        long long ptile75		= 1+(long long)(total * 0.750);
        long long ptile97_5	= 1+(long long)(total * 0.975);
        long long ptile100	= (long long)(total);
        
        locations.push_back(0); locations.push_back(ptile0_25); locations.push_back(ptile25); locations.push_back(ptile50);
        locations.push_back(ptile75); locations.push_back(ptile97_5); locations.push_back(ptile100);
        
        return locations;
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "getDefaults");
        exit(1);
    }
}
//**********************************************************************************************************************
vector<long long> Summary::getStart() {
    try {
        vector<long long> defaults = getDefaults();
        vector<long long> starts; starts.resize(7,0);
        long long meanStartPosition; meanStartPosition = 0;
        long long totalSoFar = 0;
        int lastValue = 0;
        
        //minimum
        if ((startPosition.begin())->first == -1) { starts[0] = 0; }
        else {starts[0] = (startPosition.begin())->first; }
        starts[1] = starts[0]; starts[2] = starts[0]; starts[3] = starts[0]; starts[4] = starts[0]; starts[5] = starts[0];
        
        for (map<int, long long>::iterator it = startPosition.begin(); it != startPosition.end(); it++) {
            int value = it->first; if (value == -1) { value = 0; }
            meanStartPosition += (value*it->second);
            totalSoFar += it->second;
            if (((totalSoFar <= defaults[1]) && (totalSoFar > 1)) || ((lastValue < defaults[1]) && (totalSoFar > defaults[1]))){  starts[1] = value;   } //save value
            if (((totalSoFar <= defaults[2]) && (totalSoFar > defaults[1])) ||  ((lastValue < defaults[2]) && (totalSoFar > defaults[2]))) { starts[2] = value;  } //save value
            if (((totalSoFar <= defaults[3]) && (totalSoFar > defaults[2])) ||  ((lastValue < defaults[3]) && (totalSoFar > defaults[3]))) {  starts[3] = value; } //save value
            if (((totalSoFar <= defaults[4]) && (totalSoFar > defaults[3])) ||  ((lastValue < defaults[4]) && (totalSoFar > defaults[4]))) {  starts[4] = value; } //save value
            if (((totalSoFar <= defaults[5]) && (totalSoFar > defaults[4])) ||  ((lastValue < defaults[5]) && (totalSoFar > defaults[5]))) {  starts[5] = value;  } //save value
            if ((totalSoFar <= defaults[6]) && (totalSoFar > defaults[5])) {  starts[6] = value; } //save value
            lastValue = totalSoFar;
        }
        starts[6] = (startPosition.rbegin())->first;
        starts.push_back(meanStartPosition);
        
        return starts;
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "getStart");
        exit(1);
    }
}
//**********************************************************************************************************************
long long Summary::getStart(double value) {
    try {
        long long percentage = 1+(long long)(total * value * 0.01);
        long long start = 0;
        long long totalSoFar = 0;
        long long lastValue = 0;
        
        //minimum
        if ((startPosition.begin())->first == -1) { start = 0; }
        else {start = (startPosition.begin())->first; }
    
        for (map<int, long long>::iterator it = startPosition.begin(); it != startPosition.end(); it++) {
            int value = it->first; if (value == -1) { value = 0; }
            totalSoFar += it->second;
            if (((totalSoFar <= percentage) && (totalSoFar > 1)) || ((lastValue < percentage) && (totalSoFar > percentage))){  start = value;   } //save value
            lastValue = totalSoFar;
        }
        
        return start;
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "getStart");
        exit(1);
    }
}
//**********************************************************************************************************************
vector<long long> Summary::getEnd() {
    try {
        vector<long long> defaults = getDefaults();
        vector<long long> ends; ends.resize(7,0);
        long long meanEndPosition; meanEndPosition = 0;
        long long lastValue = 0;
        long long totalSoFar = 0;
        
        if ((endPosition.begin())->first == -1) { ends[0] = 0; }
        else {ends[0] = (endPosition.begin())->first; }
        ends[1] = ends[0]; ends[2] = ends[0]; ends[3] = ends[0]; ends[4] = ends[0]; ends[5] = ends[0];
        
        for (map<int, long long>::iterator it = endPosition.begin(); it != endPosition.end(); it++) {
            int value = it->first; if (value == -1) { value = 0; }
            meanEndPosition += (value*it->second);
            totalSoFar += it->second;
            
            if (((totalSoFar <= defaults[1]) && (totalSoFar > 1)) || ((lastValue < defaults[1]) && (totalSoFar > defaults[1]))){  ends[1] = value;   } //save value
            if (((totalSoFar <= defaults[2]) && (totalSoFar > defaults[1])) ||  ((lastValue < defaults[2]) && (totalSoFar > defaults[2]))) { ends[2] = value;  } //save value
            if (((totalSoFar <= defaults[3]) && (totalSoFar > defaults[2])) ||  ((lastValue < defaults[3]) && (totalSoFar > defaults[3]))) {  ends[3] = value; } //save value
            if (((totalSoFar <= defaults[4]) && (totalSoFar > defaults[3])) ||  ((lastValue < defaults[4]) && (totalSoFar > defaults[4]))) {  ends[4] = value; } //save value
            if (((totalSoFar <= defaults[5]) && (totalSoFar > defaults[4])) ||  ((lastValue < defaults[5]) && (totalSoFar > defaults[5]))) {  ends[5] = value;  } //save value
            if ((totalSoFar <= defaults[6]) && (totalSoFar > defaults[5])) {  ends[6] = value; } //save value
            lastValue = totalSoFar;
        }
        ends[6] = (endPosition.rbegin())->first;
        ends.push_back(meanEndPosition);
        
        return ends;
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "getEnd");
        exit(1);
    }
}
//**********************************************************************************************************************
long long Summary::getEnd(double value) {
    try {
        long long percentage = 1+(long long)(total * value * 0.01);
        long long end = 0;
        long long totalSoFar = 0;
        long long lastValue = 0;
        
        //minimum
        if ((endPosition.begin())->first == -1) { end = 0; }
        else {end = (endPosition.begin())->first; }
        
        for (map<int, long long>::iterator it = endPosition.begin(); it != endPosition.end(); it++) {
            int value = it->first; if (value == -1) { value = 0; }
            totalSoFar += it->second;
            if (((totalSoFar <= percentage) && (totalSoFar > 1)) || ((lastValue < percentage) && (totalSoFar > percentage))){  end = value;   } //save value
            lastValue = totalSoFar;
        }
        
        return end;
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "getEnd");
        exit(1);
    }
}
//**********************************************************************************************************************
int Summary::summarize(string fastafile) {
    try {
        vector<unsigned long long> positions;
        vector<linePair> lines;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        int processors = 1;
        m->mothurConvert(m->getProcessors(), processors);
        positions = m->divideFile(fastafile, processors);
        for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(linePair(positions[i], positions[(i+1)]));	}
#else
        positions = m->setFilePosFasta(fastafile, numSeqs);
        if (numSeqs < processors) { processors = numSeqs; }
        
        //figure out how many sequences you have to process
        int numSeqsPerProcessor = numSeqs / processors;
        for (int i = 0; i < processors; i++) {
            int startIndex =  i * numSeqsPerProcessor;
            if(i == (processors - 1)){	numSeqsPerProcessor = numSeqs - i * numSeqsPerProcessor; 	}
            lines.push_back(new linePair(positions[startIndex], numSeqsPerProcessor));
        }
#endif
        
        int process = 1;
        int num = 0;
        vector<int> processIDS;
        bool recalc = false;
        
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        
        //loop through and create all the processes you want
        while (process != processors) {
            pid_t pid = fork();
            
            if (pid > 0) {
                processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                process++;
            }else if (pid == 0){
                num = driverSummarize(fastafile, lines[process]);
                
                //pass numSeqs to parent
                ofstream out;
                string tempFile = m->mothurGetpid(process) + ".num.temp";
                m->openOutputFile(tempFile, out);
                
                out << num << endl;
                out << startPosition.size() << endl;
                for (map<int,  long long>::iterator it = startPosition.begin(); it != startPosition.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                out << endPosition.size() << endl;
                for (map<int,  long long>::iterator it = endPosition.begin(); it != endPosition.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                out << seqLength.size() << endl;
                for (map<int,  long long>::iterator it = seqLength.begin(); it != seqLength.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                out << ambigBases.size() << endl;
                for (map<int,  long long>::iterator it = ambigBases.begin(); it != ambigBases.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                out << longHomoPolymer.size() << endl;
                for (map<int,  long long>::iterator it = longHomoPolymer.begin(); it != longHomoPolymer.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                out.close();
                
                exit(0);
            }else {
                m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(process) + "\n"); processors = process;
                for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                //wait to die
                for (int i=0;i<processIDS.size();i++) {
                    int temp = processIDS[i];
                    wait(&temp);
                }
                m->control_pressed = false;
                for (int i=0;i<processIDS.size();i++) {
                    m->mothurRemove((toString(processIDS[i]) + ".num.temp"));
                }
                recalc = true;
                break;
            }
        }
        
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->control_pressed = false;  for (int i=0;i<processIDS.size();i++) {m->mothurRemove(fastafile + (toString(processIDS[i]) + ".num.temp"));}processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            
            //redo file divide
            lines.clear();
            vector<unsigned long long> positions = m->divideFile(fastafile, processors);
            for (int i = 0; i < (positions.size()-1); i++) {  lines.push_back(linePair(positions[i], positions[(i+1)]));  }
            
            startPosition.clear();
            endPosition.clear();
            seqLength.clear();
            ambigBases.clear();
            longHomoPolymer.clear();
            
            num = 0;
            processIDS.resize(0);
            process = 1;
            
            //loop through and create all the processes you want
            while (process != processors) {
                pid_t pid = fork();
                
                if (pid > 0) {
                    processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                    process++;
                }else if (pid == 0){
                    num = driverSummarize(fastafile, lines[process]);
                    
                    //pass numSeqs to parent
                    ofstream out;
                    string tempFile = m->mothurGetpid(process) + ".num.temp";
                    m->openOutputFile(tempFile, out);
                    
                    out << num << endl;
                    out << total << endl;
                    out << startPosition.size() << endl;
                    for (map<int,  long long>::iterator it = startPosition.begin(); it != startPosition.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                    out << endPosition.size() << endl;
                    for (map<int,  long long>::iterator it = endPosition.begin(); it != endPosition.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                    out << seqLength.size() << endl;
                    for (map<int,  long long>::iterator it = seqLength.begin(); it != seqLength.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                    out << ambigBases.size() << endl;
                    for (map<int,  long long>::iterator it = ambigBases.begin(); it != ambigBases.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                    out << longHomoPolymer.size() << endl;
                    for (map<int,  long long>::iterator it = longHomoPolymer.begin(); it != longHomoPolymer.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                    out.close();
                    
                    exit(0);
                }else {
                    m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine();
                    for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                    exit(0);
                }
            }
        }
        
        
        //do your part
        num = driverSummarize(fastafile, lines[0]);
        
        //force parent to wait until all the processes are done
        for (int i=0;i<processIDS.size();i++) {
            int temp = processIDS[i];
            wait(&temp);
        }
        
        //parent reads in and combine Filter info
        for (int i = 0; i < processIDS.size(); i++) {
            string tempFilename = toString(processIDS[i]) + ".num.temp";
            ifstream in;
            m->openInputFile(tempFilename, in);
            
            long long  tempNum;
            in >> tempNum; m->gobble(in); num += tempNum;
            in >> tempNum; m->gobble(in); total += tempNum;
            in >> tempNum; m->gobble(in);
            for (int k = 0; k < tempNum; k++)			{
                long long first, second;
                in >> first; m->gobble(in); in >> second; m->gobble(in);
                map<int,  long long>::iterator it = startPosition.find(first);
                if (it == startPosition.end()) { startPosition[first] = second; } //first finding of this start position, set count.
                else { it->second += second; } //add counts
            }
            m->gobble(in);
            in >> tempNum; m->gobble(in);
            for (int k = 0; k < tempNum; k++)			{
                long long first, second;
                in >> first; m->gobble(in); in >> second; m->gobble(in);
                map<int,  long long>::iterator it = endPosition.find(first);
                if (it == endPosition.end()) { endPosition[first] = second; } //first finding of this end position, set count.
                else { it->second += second; } //add counts
            }
            m->gobble(in);
            in >> tempNum; m->gobble(in);
            for (int k = 0; k < tempNum; k++)			{
                long long first, second;
                in >> first; m->gobble(in); in >> second; m->gobble(in);
                map<int,  long long>::iterator it = seqLength.find(first);
                if (it == seqLength.end()) { seqLength[first] = second; } //first finding of this end position, set count.
                else { it->second += second; } //add counts
            }
            m->gobble(in);
            in >> tempNum; m->gobble(in);
            for (int k = 0; k < tempNum; k++)			{
                long long first, second;
                in >> first; m->gobble(in); in >> second; m->gobble(in);
                map<int,  long long>::iterator it = ambigBases.find(first);
                if (it == ambigBases.end()) { ambigBases[first] = second; } //first finding of this end position, set count.
                else { it->second += second; } //add counts
            }
            m->gobble(in);
            in >> tempNum; m->gobble(in);
            for (int k = 0; k < tempNum; k++)			{
                long long first, second;
                in >> first; m->gobble(in); in >> second; m->gobble(in);
                map<int,  long long>::iterator it = longHomoPolymer.find(first);
                if (it == longHomoPolymer.end()) { longHomoPolymer[first] = second; } //first finding of this end position, set count.
                else { it->second += second; } //add counts
            }
            m->gobble(in);
            
            in.close();
            m->mothurRemove(tempFilename);
        }
        
#else
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        //Windows version shared memory, so be careful when passing variables through the seqSumData struct.
        //Above fork() will clone, so memory is separate, but that's not the case with windows,
        //Taking advantage of shared memory to allow both threads to add info to vectors.
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        
        vector<seqSumData*> pDataArray;
        DWORD   dwThreadIdArray[processors-1];
        HANDLE  hThreadArray[processors-1];
        
        //Create processor worker threads.
        for( int i=0; i<processors-1; i++ ){
            
            // Allocate memory for thread data.
            seqSumData* tempSum = new seqSumData(fastafile, m, lines[i]->start, lines[i]->end, hasNameOrCount, nameMap);
            pDataArray.push_back(tempSum);
            
            //MySeqSumThreadFunction is in header. It must be global or static to work with the threads.
            //default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
            hThreadArray[i] = CreateThread(NULL, 0, MySeqSumThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);
        }
        
        //do your part
        num = driverSummarize(fastafile, lines[processors-1]);
        processIDS.push_back(processors-1);
        
        //Wait until all threads have terminated.
        WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
        
        //Close all thread handles and free memory allocations.
        for(int i=0; i < pDataArray.size(); i++){
            num += pDataArray[i]->count;
            total += pDataArray[i]->total;
            if (pDataArray[i]->count != pDataArray[i]->end) {
                m->mothurOut("[ERROR]: process " + toString(i) + " only processed " + toString(pDataArray[i]->count) + " of " + toString(pDataArray[i]->end) + " sequences assigned to it, quitting. \n"); m->control_pressed = true;
            }
            for (map<int, long long>::iterator it = pDataArray[i]->startPosition.begin(); it != pDataArray[i]->startPosition.end(); it++)		{
                map<int, long long>::iterator itMain = startPosition.find(it->first);
                if (itMain == startPosition.end()) { //newValue
                    startPosition[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = pDataArray[i]->endPosition.begin(); it != pDataArray[i]->endPosition.end(); it++)		{
                map<int, long long>::iterator itMain = endPosition.find(it->first);
                if (itMain == endPosition.end()) { //newValue
                    endPosition[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = pDataArray[i]->seqLength.begin(); it != pDataArray[i]->seqLength.end(); it++)		{
                map<int, long long>::iterator itMain = seqLength.find(it->first);
                if (itMain == seqLength.end()) { //newValue
                    seqLength[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = pDataArray[i]->ambigBases.begin(); it != pDataArray[i]->ambigBases.end(); it++)		{
                map<int, long long>::iterator itMain = ambigBases.find(it->first);
                if (itMain == ambigBases.end()) { //newValue
                    ambigBases[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = pDataArray[i]->longHomoPolymer.begin(); it != pDataArray[i]->longHomoPolymer.end(); it++)		{
                map<int, long long>::iterator itMain = longHomoPolymer.find(it->first);
                if (itMain == longHomoPolymer.end()) { //newValue
                    longHomoPolymer[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            CloseHandle(hThreadArray[i]);
            delete pDataArray[i];
        }
#endif
        
        if (hasNameOrCount) {
            if (nameCountNumUniques != num) { // do fasta and name/count files match
                m->mothurOut("[ERROR]: Your " + type + " file contains " + toString(nameCountNumUniques) + " unique sequences, but your fasta file contains " + toString(num) + ". File mismatch detected, quitting command.\n"); m->control_pressed = true;
            }
        }
        
        return num; //number of uniques
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "summarize");
        exit(1);
    }
}
//**********************************************************************************************************************
int Summary::driverSummarize(string fastafile, linePair lines) {
    try {
        ifstream in;
        m->openInputFile(fastafile, in);
        
        in.seekg(lines.start);
        
        //print header if you are process 0
        if (lines.start == 0) { m->zapGremlins(in); m->gobble(in); }
        
        bool done = false;
        long long count = 0;
        
        while (!done) {
            
            if (m->control_pressed) { in.close(); return 1; }
            
            Sequence current(in); m->gobble(in);

            if (current.getName() != "") {
                
                if (m->debug) { m->mothurOut("[DEBUG]: " + current.getName() + '\t' + toString(current.getNumBases()) + "\n");  }
    
                addSeq(current); count++;
            }
            
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
            unsigned long long pos = in.tellg();
            if ((pos == -1) || (pos >= lines.end)) { break; }
#else
            if (in.eof()) { break; }
#endif
        }
        
        in.close();
        
        return count;
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "driverSummarize");
        exit(1);
    }
}
//**********************************************************************************************************************


