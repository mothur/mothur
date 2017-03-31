//
//  summary.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/27/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "summary.hpp"


//**********************************************************************************************************************
void Summary::processNameCount(string n) { //name or count file to include in counts
    try {
        nameMap.clear(); nameCountNumUniques = 0;
        if (n != "") {
            hasNameOrCount = true;
            if (m->isCountFile(n)) {
                CountTable ct;
                ct.readTable(n, false, false);
                nameMap = ct.getNameMap();
                type = "count";
            }else { nameMap = m->readNames(n); type = "name"; }
        }
        nameCountNumUniques = nameMap.size();
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "Summary");
        exit(1);
    }
}

//**********************************************************************************************************************
string Summary::addSeq(Sequence seq) {
    try {
        long long num = 1;
        
        if (hasNameOrCount) {
            //make sure this sequence is in the namefile, else error
            itFindName = nameMap.find(seq.getName());
            
            if (itFindName == nameMap.end()) { m->mothurOut("[ERROR]: '" + seq.getName() + "' is not in your name or count file, please correct."); m->mothurOutEndLine(); m->control_pressed = true; }
            else { num = itFindName->second; }
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
        
        string output = "";
        output += seq.getName() + '\t';
        output += toString(thisStartPosition) + '\t' + toString(thisEndPosition) + '\t';
        output += toString(thisSeqLength) + '\t' + toString(thisAmbig) + '\t';
        output += toString(thisHomoP) + '\t' + toString(num);
        
        return output;
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "addSeq - FastaFile");
        exit(1);
    }
}
//**********************************************************************************************************************
//string seqInfo = addSeq(name, start, end, length, ambigs, polymer, numReps);
string Summary::addSeq(string name, int start, int end, int length, int ambigs, int polymer, long long numReps) {
    try {
        
        if ((numReps == 1) && hasNameOrCount) {
            //make sure this sequence is in the namefile, else error
            itFindName = nameMap.find(name);
            
            if (itFindName == nameMap.end()) { m->mothurOut("[ERROR]: '" + name + "' is not in your name or count file, please correct."); m->mothurOutEndLine(); m->control_pressed = true; }
            else { numReps = itFindName->second; }
        }
        
        it = startPosition.find(start);
        if (it == startPosition.end()) { startPosition[start] = numReps; } //first finding of this start position, set count.
        else { it->second += numReps; } //add counts
        
        it = endPosition.find(end);
        if (it == endPosition.end()) { endPosition[end] = numReps; } //first finding of this end position, set count.
        else { it->second += numReps; } //add counts
        
        it = seqLength.find(length);
        if (it == seqLength.end()) { seqLength[length] = numReps; } //first finding of this length, set count.
        else { it->second += numReps; } //add counts
        
        it = ambigBases.find(ambigs);
        if (it == ambigBases.end()) { ambigBases[ambigs] = numReps; } //first finding of this ambig, set count.
        else { it->second += numReps; } //add counts
        
        it = longHomoPolymer.find(polymer);
        if (it == longHomoPolymer.end()) { longHomoPolymer[polymer] = numReps; } //first finding of this homop, set count.
        else { it->second += numReps; } //add counts
        
        numUniques++;
        total += numReps;
        
        string output = "";
        output += name + '\t';
        output += toString(start) + '\t' + toString(end) + '\t';
        output += toString(length) + '\t' + toString(ambigs) + '\t';
        output += toString(polymer) + '\t' + toString(numReps);
        
        return output;
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "addSeq - FastaSummaryFile");
        exit(1);
    }
}
//**********************************************************************************************************************
//string seqInfo = addSeq(name, length, OLength, thisOStart, thisOEnd, numMisMatches, numns);
string Summary::addSeq(string name, int length, int olength, int ostart, int oend, int mismatches, int numns) {
    try {
        long long numReps = 1;
        if (hasNameOrCount) {
            //make sure this sequence is in the namefile, else error
            itFindName = nameMap.find(name);
            
            if (itFindName == nameMap.end()) { m->mothurOut("[ERROR]: '" + name + "' is not in your name or count file, please correct."); m->mothurOutEndLine(); m->control_pressed = true; }
            else { numReps = itFindName->second; }
        }
        
        it = ostartPosition.find(ostart);
        if (it == ostartPosition.end()) { ostartPosition[ostart] = numReps; } //first finding of this start position, set count.
        else { it->second += numReps; } //add counts
        
        it = oendPosition.find(oend);
        if (it == oendPosition.end()) { oendPosition[oend] = numReps; } //first finding of this end position, set count.
        else { it->second += numReps; } //add counts
        
        it = oseqLength.find(olength);
        if (it == oseqLength.end()) { oseqLength[olength] = numReps; } //first finding of this length, set count.
        else { it->second += numReps; } //add counts
        
        it = seqLength.find(length);
        if (it == seqLength.end()) { seqLength[length] = numReps; } //first finding of this length, set count.
        else { it->second += numReps; } //add counts

        it = misMatches.find(mismatches);
        if (it == misMatches.end()) { misMatches[mismatches] = numReps; } //first finding of this ambig, set count.
        else { it->second += numReps; } //add counts
        
        it = numNs.find(numns);
        if (it == numNs.end()) { numNs[numns] = numReps; } //first finding of this homop, set count.
        else { it->second += numReps; } //add counts
        
        numUniques++;
        total += numReps;
        
        string output = "";
        output += name + '\t';
        output += toString(ostart) + '\t' + toString(oend) + '\t';
        output += toString(length) + '\t' + toString(olength) + '\t'+ toString(mismatches) + '\t';
        output += toString(numns) + '\t' + toString(numReps);
        
        return output;
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "addSeq - ContigsSummaryFile");
        exit(1);
    }
}
//**********************************************************************************************************************
//string seqInfo = addSeq(name, sim, score, insert);
string Summary::addSeq(string name, int length, float SimBtwnQueryTemplate, float SearchScore, int LongestInserts) {
    try {
        long long numReps = 1;
        if (hasNameOrCount) {
            //make sure this sequence is in the namefile, else error
            itFindName = nameMap.find(name);
            
            if (itFindName == nameMap.end()) { m->mothurOut("[ERROR]: '" + name + "' is not in your name or count file, please correct."); m->mothurOutEndLine(); m->control_pressed = true; }
            else { numReps = itFindName->second; }
        }
        
        map<float, long long>:: iterator itFloat = sims.find(SimBtwnQueryTemplate);
        if (itFloat == sims.end()) { sims[SimBtwnQueryTemplate] = numReps; } //first finding of this similarity score, set count.
        else { itFloat->second += numReps; } //add counts
        
        itFloat = scores.find(SearchScore);
        if (itFloat == scores.end()) { scores[SearchScore] = numReps; } //first finding of this end position, set count.
        else { itFloat->second += numReps; } //add counts
        
        it = inserts.find(LongestInserts);
        if (it == inserts.end()) { inserts[LongestInserts] = numReps; } //first finding of this length, set count.
        else { it->second += numReps; } //add counts
        
        it = seqLength.find(length);
        if (it == seqLength.end()) { seqLength[length] = numReps; } //first finding of this length, set count.
        else { it->second += numReps; } //add counts
        
        numUniques++;
        total += numReps;
        
        string output = "";
        output += name + '\t';
        output += toString(SimBtwnQueryTemplate) + '\t' + toString(SearchScore) + '\t';
        output += toString(LongestInserts) + '\t' + toString(numReps);
        
        return output;
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "addSeq - ContigsSummaryFile");
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
        
        locations.push_back(1); locations.push_back(ptile0_25); locations.push_back(ptile25); locations.push_back(ptile50);
        locations.push_back(ptile75); locations.push_back(ptile97_5); locations.push_back(ptile100);
        
        return locations;
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "getDefaults");
        exit(1);
    }
}
//**********************************************************************************************************************
vector<long long> Summary::getValues(map<int, long long>& positions) {
    try {
        vector<long long> defaults = getDefaults();
        vector<long long> results; results.resize(7,0);
        long long meanPosition; meanPosition = 0;
        long long totalSoFar = 0;
        int lastValue = 0;
        
        //minimum
        if ((positions.begin())->first == -1) { results[0] = 0; }
        else {results[0] = (positions.begin())->first; }
        results[1] = results[0]; results[2] = results[0]; results[3] = results[0]; results[4] = results[0]; results[5] = results[0];
        
        for (map<int, long long>::iterator it = positions.begin(); it != positions.end(); it++) {
            int value = it->first; if (value == -1) { value = 0; }
            meanPosition += (value*it->second);
            totalSoFar += it->second;
            if (((totalSoFar <= defaults[1]) && (totalSoFar > 1)) || ((lastValue < defaults[1]) && (totalSoFar > defaults[1]))){  results[1] = value;   } //save value
            if (((totalSoFar <= defaults[2]) && (totalSoFar > defaults[1])) ||  ((lastValue < defaults[2]) && (totalSoFar > defaults[2]))) { results[2] = value;  } //save value
            if (((totalSoFar <= defaults[3]) && (totalSoFar > defaults[2])) ||  ((lastValue < defaults[3]) && (totalSoFar > defaults[3]))) {  results[3] = value; } //save value
            if (((totalSoFar <= defaults[4]) && (totalSoFar > defaults[3])) ||  ((lastValue < defaults[4]) && (totalSoFar > defaults[4]))) {  results[4] = value; } //save value
            if (((totalSoFar <= defaults[5]) && (totalSoFar > defaults[4])) ||  ((lastValue < defaults[5]) && (totalSoFar > defaults[5]))) {  results[5] = value;  } //save value
            if ((totalSoFar <= defaults[6]) && (totalSoFar > defaults[5])) {  results[6] = value; } //save value
            lastValue = totalSoFar;
        }
        results[6] = (positions.rbegin())->first;
        
        double meansPosition = meanPosition / (double) total;
        results.push_back(meansPosition);
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "getValues");
        exit(1);
    }
}
//**********************************************************************************************************************
long long Summary::getValue(map<int, long long>& positions, double value) {
    try {
        long long percentage = 1+(long long)(total * value * 0.01);
        long long result = 0;
        long long totalSoFar = 0;
        long long lastValue = 0;
        
        //minimum
        if ((positions.begin())->first == -1) { result = 0; }
        else {result = (positions.begin())->first; }
    
        for (map<int, long long>::iterator it = positions.begin(); it != positions.end(); it++) {
            int value = it->first; if (value == -1) { value = 0; }
            totalSoFar += it->second;
            if (((totalSoFar <= percentage) && (totalSoFar > 1)) || ((lastValue < percentage) && (totalSoFar > percentage))){  result = value;   } //save value
            lastValue = totalSoFar;
        }
        
        return result;
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "getValue");
        exit(1);
    }
}
//**********************************************************************************************************************
vector<long long> Summary::getValues(map<float, long long>& positions) {
    try {
        vector<long long> defaults = getDefaults();
        vector<long long> results; results.resize(7,0);
        long long meanPosition; meanPosition = 0;
        long long totalSoFar = 0;
        int lastValue = 0;
        
        //minimum
        if ((positions.begin())->first == -1) { results[0] = 0; }
        else {results[0] = (positions.begin())->first; }
        results[1] = results[0]; results[2] = results[0]; results[3] = results[0]; results[4] = results[0]; results[5] = results[0];
        
        for (map<float, long long>::iterator it = positions.begin(); it != positions.end(); it++) {
            int value = it->first; if (value == -1) { value = 0; }
            meanPosition += (value*it->second);
            totalSoFar += it->second;
            if (((totalSoFar <= defaults[1]) && (totalSoFar > 1)) || ((lastValue < defaults[1]) && (totalSoFar > defaults[1]))){  results[1] = value;   } //save value
            if (((totalSoFar <= defaults[2]) && (totalSoFar > defaults[1])) ||  ((lastValue < defaults[2]) && (totalSoFar > defaults[2]))) { results[2] = value;  } //save value
            if (((totalSoFar <= defaults[3]) && (totalSoFar > defaults[2])) ||  ((lastValue < defaults[3]) && (totalSoFar > defaults[3]))) {  results[3] = value; } //save value
            if (((totalSoFar <= defaults[4]) && (totalSoFar > defaults[3])) ||  ((lastValue < defaults[4]) && (totalSoFar > defaults[4]))) {  results[4] = value; } //save value
            if (((totalSoFar <= defaults[5]) && (totalSoFar > defaults[4])) ||  ((lastValue < defaults[5]) && (totalSoFar > defaults[5]))) {  results[5] = value;  } //save value
            if ((totalSoFar <= defaults[6]) && (totalSoFar > defaults[5])) {  results[6] = value; } //save value
            lastValue = totalSoFar;
        }
        results[6] = (positions.rbegin())->first;
        
        double meansPosition = meanPosition / (double) total;
        results.push_back(meansPosition);
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "getValues");
        exit(1);
    }
}
//**********************************************************************************************************************
long long Summary::getValue(map<float, long long>& positions, double value) {
    try {
        long long percentage = 1+(long long)(total * value * 0.01);
        long long result = 0;
        long long totalSoFar = 0;
        long long lastValue = 0;
        
        //minimum
        if ((positions.begin())->first == -1) { result = 0; }
        else {result = (positions.begin())->first; }
        
        for (map<float, long long>::iterator it = positions.begin(); it != positions.end(); it++) {
            int value = it->first; if (value == -1) { value = 0; }
            totalSoFar += it->second;
            if (((totalSoFar <= percentage) && (totalSoFar > 1)) || ((lastValue < percentage) && (totalSoFar > percentage))){  result = value;   } //save value
            lastValue = totalSoFar;
        }
        
        return result;
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "getValue");
        exit(1);
    }
}
//**********************************************************************************************************************
long long Summary::summarizeFasta(string fastafile, string n, string output) {
    try {
        //fill namemap
        processNameCount(n);
        return (summarizeFasta(fastafile, output));
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "getHomop");
        exit(1);
    }
}
//**********************************************************************************************************************
long long Summary::summarizeFasta(string fastafile, string output) {
    try {
        
        vector<unsigned long long> positions;
        vector<linePair> lines;
        string p = m->getProcessors();
        int processors = 1; m->mothurConvert(p, processors);
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
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
            lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
        }
#endif
        
        int process = 1;
        long long num = 0;
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
                num = driverSummarize(fastafile, (output + m->mothurGetpid(process) + ".temp"), lines[process]);
                
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
                    string outputName = output + m->mothurGetpid(process) + ".temp";
                    if (output == "") {  outputName = "";  }
                    num = driverSummarize(fastafile, outputName, lines[process]);
                    
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
        num = driverSummarize(fastafile, output, lines[0]);
        
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
            string extension = "";
            if (i != 0) { extension = toString(i) + ".temp"; processIDS.push_back(i); }
            string outputName = output + extension;
            if (output == "") {  outputName = "";  }
            seqSumData* tempSum = new seqSumData(fastafile, outputName, m, lines[i].start, lines[i].end, hasNameOrCount, nameMap);
            pDataArray.push_back(tempSum);
            
            //MySeqSumThreadFunction is in header. It must be global or static to work with the threads.
            //default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
            hThreadArray[i] = CreateThread(NULL, 0, MySeqSumThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);
        }
        
        //do your part
        string extension = toString(processors-1) + ".temp";
        string outputName = output + extension;
        if (output == "") {  outputName = "";  }
        num = driverSummarize(fastafile, outputName, lines[processors-1]);
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
        //append files
        for(int i=0;i<processIDS.size();i++){
            if (output != "") {
                m->appendFiles((output + toString(processIDS[i]) + ".temp"), output);
                m->mothurRemove((output + toString(processIDS[i]) + ".temp"));
            }
        }
        
        if (hasNameOrCount) {
            if (nameCountNumUniques != num) { // do fasta and name/count files match
                m->mothurOut("[ERROR]: Your " + type + " file contains " + toString(nameCountNumUniques) + " unique sequences, but your fasta file contains " + toString(num) + ". File mismatch detected, quitting command.\n"); m->control_pressed = true;
            }
        }
        
        numUniques = num;
        
        return num; //number of uniques
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "summarizeFasta");
        exit(1);
    }
}
//**********************************************************************************************************************
int Summary::driverSummarize(string fastafile, string output, linePair lines) {
    try {
        ofstream out;
        if (output != "") { m->openOutputFile(output, out); }
        
        ifstream in;
        m->openInputFile(fastafile, in);
        
        in.seekg(lines.start);
        
        //print header if you are process 0
        if (lines.start == 0) {
            m->zapGremlins(in); m->gobble(in);
            //print header if you are process 0
            out << "seqname\tstart\tend\tnbases\tambigs\tpolymer\tnumSeqs" << endl;
        }
        
        bool done = false;
        long long count = 0;
        
        while (!done) {
            
            if (m->control_pressed) { in.close(); return 1; }
            
            Sequence current(in); m->gobble(in);

            if (current.getName() != "") {
                
                if (m->debug) { m->mothurOut("[DEBUG]: " + current.getName() + "\t" + toString(current.getStartPos()) + "\t" + toString(current.getEndPos()) + "\t" + toString(current.getNumBases()) + "\n"); }
    
                string seqInfo = addSeq(current); count++;
                if (output != "") { out << seqInfo << endl; }
            }
            
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
            unsigned long long pos = in.tellg();
            if ((pos == -1) || (pos >= lines.end)) { break; }
#else
            if (in.eof()) { break; }
#endif
        }
        
        if (output != "") { out.close(); }
        in.close();
        
        return count;
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "driverSummarize");
        exit(1);
    }
}
//**********************************************************************************************************************
long long Summary::summarizeFastaSummary(string summaryfile, string n) {
    try {
        //fill namemap
        processNameCount(n);
        return (summarizeFastaSummary(summaryfile));
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "summarizeFastaSummary");
        exit(1);
    }
}
/**********************************************************************************************************************/
long long Summary::summarizeFastaSummary(string summaryfile) {
    try {
        
        vector<unsigned long long> positions;
        vector<linePair> lines;
        string p = m->getProcessors();
        int processors = 1; m->mothurConvert(p, processors);
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        positions = m->divideFile(summaryfile, processors);
        for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(linePair(positions[i], positions[(i+1)]));	}
#else
        positions = m->setFilePosFasta(summaryfile, numSeqs);
        if (numSeqs < processors) { processors = numSeqs; }
        
        //figure out how many sequences you have to process
        int numSeqsPerProcessor = numSeqs / processors;
        for (int i = 0; i < processors; i++) {
            int startIndex =  i * numSeqsPerProcessor;
            if(i == (processors - 1)){	numSeqsPerProcessor = numSeqs - i * numSeqsPerProcessor; 	}
            lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
        }
#endif
        
        int process = 1;
        long long num = 0;
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
                num = driverFastaSummarySummarize(summaryfile, lines[process]);
                
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
            vector<unsigned long long> positions = m->divideFile(summaryfile, processors);
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
                    
                    num = driverFastaSummarySummarize(summaryfile, lines[process]);
                    
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
        num = driverFastaSummarySummarize(summaryfile, lines[0]);
        
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
            seqSumData* tempSum = new seqSumData(summaryfile, m, lines[i]->start, lines[i]->end, hasNameOrCount, nameMap);
            pDataArray.push_back(tempSum);
            
            //MySeqSumThreadFunction is in header. It must be global or static to work with the threads.
            //default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
            hThreadArray[i] = CreateThread(NULL, 0, MySeqFastaSumThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);
        }
        
        //do your part
        num = driverFastaSummarySummarize(summaryfile, lines[processors-1]);
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
        
        numUniques = num;
        
        return num; //number of uniques
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "summarizeFastaSummary");
        exit(1);
    }
}
//**********************************************************************************************************************
int Summary::driverFastaSummarySummarize(string fastafile, linePair lines) {
    try {
        ifstream in;
        m->openInputFile(fastafile, in);
        
        in.seekg(lines.start);
        
        //print header if you are process 0
        if (lines.start == 0) { m->zapGremlins(in); m->getline(in); m->gobble(in); }
        
        bool done = false;
        long long count = 0;
        string name;
        int start, end, length, ambigs, polymer;
        long long numReps;
        
        while (!done) {
            
            if (m->control_pressed) { in.close(); return 1; }
            
            //seqname	start	end	nbases	ambigs	polymer	numSeqs
            in >> name >> start >> end >> length >> ambigs >> polymer >> numReps; m->gobble(in);
            
            if (m->debug) { m->mothurOut("[DEBUG]: " + name + "\t" + toString(start) + "\t" + toString(end) + "\t" + toString(length) + "\n"); }
            
            if (name != "") {
                string seqInfo = addSeq(name, start, end, length, ambigs, polymer, numReps); count++;
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
long long Summary::summarizeContigsSummary(string summaryfile, string n) {
    try {
        //fill namemap
        processNameCount(n);
        return (summarizeContigsSummary(summaryfile));
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "summarizeFastaSummary");
        exit(1);
    }
}
/**********************************************************************************************************************/
long long Summary::summarizeContigsSummary(string summaryfile) {
    try {
        
        vector<unsigned long long> positions;
        vector<linePair> lines;
        string p = m->getProcessors();
        int processors = 1; m->mothurConvert(p, processors);
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        positions = m->divideFile(summaryfile, processors);
        for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(linePair(positions[i], positions[(i+1)]));	}
#else
        positions = m->setFilePosFasta(summaryfile, numSeqs);
        if (numSeqs < processors) { processors = numSeqs; }
        
        //figure out how many sequences you have to process
        int numSeqsPerProcessor = numSeqs / processors;
        for (int i = 0; i < processors; i++) {
            int startIndex =  i * numSeqsPerProcessor;
            if(i == (processors - 1)){	numSeqsPerProcessor = numSeqs - i * numSeqsPerProcessor; 	}
            lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
        }
#endif
        
        int process = 1;
        long long num = 0;
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
                num = driverContigsSummarySummarize(summaryfile, lines[process]);
                
                //pass numSeqs to parent
                ofstream out;
                string tempFile = m->mothurGetpid(process) + ".num.temp";
                m->openOutputFile(tempFile, out);
                
                out << num << endl;
                out << total << endl;
                out << ostartPosition.size() << endl;
                for (map<int,  long long>::iterator it = ostartPosition.begin(); it != ostartPosition.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                out << oendPosition.size() << endl;
                for (map<int,  long long>::iterator it = oendPosition.begin(); it != oendPosition.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                out << oseqLength.size() << endl;
                for (map<int,  long long>::iterator it = oseqLength.begin(); it != oseqLength.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                out << seqLength.size() << endl;
                for (map<int,  long long>::iterator it = seqLength.begin(); it != seqLength.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                out << misMatches.size() << endl;
                for (map<int,  long long>::iterator it = misMatches.begin(); it != misMatches.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                out << numNs.size() << endl;
                for (map<int,  long long>::iterator it = numNs.begin(); it != numNs.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
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
            vector<unsigned long long> positions = m->divideFile(summaryfile, processors);
            for (int i = 0; i < (positions.size()-1); i++) {  lines.push_back(linePair(positions[i], positions[(i+1)]));  }
            
            ostartPosition.clear();
            oendPosition.clear();
            oseqLength.clear();
            seqLength.clear();
            misMatches.clear();
            numNs.clear();
            
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
                    
                    num = driverContigsSummarySummarize(summaryfile, lines[process]);
                    
                    //pass numSeqs to parent
                    ofstream out;
                    string tempFile = m->mothurGetpid(process) + ".num.temp";
                    m->openOutputFile(tempFile, out);
                    
                    out << num << endl;
                    out << total << endl;
                    out << ostartPosition.size() << endl;
                    for (map<int,  long long>::iterator it = ostartPosition.begin(); it != ostartPosition.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                    out << oendPosition.size() << endl;
                    for (map<int,  long long>::iterator it = oendPosition.begin(); it != oendPosition.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                    out << oseqLength.size() << endl;
                    for (map<int,  long long>::iterator it = oseqLength.begin(); it != oseqLength.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                    out << seqLength.size() << endl;
                    for (map<int,  long long>::iterator it = seqLength.begin(); it != seqLength.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                    out << misMatches.size() << endl;
                    for (map<int,  long long>::iterator it = misMatches.begin(); it != misMatches.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                    out << numNs.size() << endl;
                    for (map<int,  long long>::iterator it = numNs.begin(); it != numNs.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
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
        num = driverContigsSummarySummarize(summaryfile, lines[0]);
        
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
                map<int,  long long>::iterator it = ostartPosition.find(first);
                if (it == ostartPosition.end()) { ostartPosition[first] = second; } //first finding of this start position, set count.
                else { it->second += second; } //add counts
            }
            m->gobble(in);
            in >> tempNum; m->gobble(in);
            for (int k = 0; k < tempNum; k++)			{
                long long first, second;
                in >> first; m->gobble(in); in >> second; m->gobble(in);
                map<int,  long long>::iterator it = oendPosition.find(first);
                if (it == oendPosition.end()) { oendPosition[first] = second; } //first finding of this end position, set count.
                else { it->second += second; } //add counts
            }
            m->gobble(in);
            in >> tempNum; m->gobble(in);
            for (int k = 0; k < tempNum; k++)			{
                long long first, second;
                in >> first; m->gobble(in); in >> second; m->gobble(in);
                map<int,  long long>::iterator it = oseqLength.find(first);
                if (it == oseqLength.end()) { oseqLength[first] = second; } //first finding of this end position, set count.
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
                map<int,  long long>::iterator it = misMatches.find(first);
                if (it == misMatches.end()) { misMatches[first] = second; } //first finding of this end position, set count.
                else { it->second += second; } //add counts
            }
            m->gobble(in);
            in >> tempNum; m->gobble(in);
            for (int k = 0; k < tempNum; k++)			{
                long long first, second;
                in >> first; m->gobble(in); in >> second; m->gobble(in);
                map<int,  long long>::iterator it = numNs.find(first);
                if (it == numNs.end()) { numNs[first] = second; } //first finding of this end position, set count.
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
            seqSumData* tempSum = new seqSumData(summaryfile, m, lines[i]->start, lines[i]->end, hasNameOrCount, nameMap);
            pDataArray.push_back(tempSum);
            
            //MySeqSumThreadFunction is in header. It must be global or static to work with the threads.
            //default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
            hThreadArray[i] = CreateThread(NULL, 0, MySeqFastaSumThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);
        }
        
        //do your part
        num = driverContigsSummarySummarize(summaryfile, lines[processors-1]);
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
            for (map<int, long long>::iterator it = pDataArray[i]->ostartPosition.begin(); it != pDataArray[i]->ostartPosition.end(); it++)		{
                map<int, long long>::iterator itMain = ostartPosition.find(it->first);
                if (itMain == ostartPosition.end()) { //newValue
                    ostartPosition[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = pDataArray[i]->oendPosition.begin(); it != pDataArray[i]->oendPosition.end(); it++)		{
                map<int, long long>::iterator itMain = oendPosition.find(it->first);
                if (itMain == oendPosition.end()) { //newValue
                    oendPosition[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = pDataArray[i]->oseqLength.begin(); it != pDataArray[i]->oseqLength.end(); it++)		{
                map<int, long long>::iterator itMain = oseqLength.find(it->first);
                if (itMain == oseqLength.end()) { //newValue
                    oseqLength[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = pDataArray[i]->seqLength.begin(); it != pDataArray[i]->seqLength.end(); it++)		{
                map<int, long long>::iterator itMain = seqLength.find(it->first);
                if (itMain == seqLength.end()) { //newValue
                    seqLength[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = pDataArray[i]->misMatches.begin(); it != pDataArray[i]->misMatches.end(); it++)		{
                map<int, long long>::iterator itMain = misMatches.find(it->first);
                if (itMain == misMatches.end()) { //newValue
                    misMatches[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = pDataArray[i]->numNs.begin(); it != pDataArray[i]->numNs.end(); it++)		{
                map<int, long long>::iterator itMain = numNs.find(it->first);
                if (itMain == numNs.end()) { //newValue
                    numNs[it->first] = it->second;
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
        
        numUniques = num;
        
        return num; //number of uniques
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "summarizeFastaSummary");
        exit(1);
    }
}
//**********************************************************************************************************************
long long Summary::summarizeAlignSummary(string summaryfile, string n) {
    try {
        //fill namemap
        processNameCount(n);
        return (summarizeAlignSummary(summaryfile));
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "summarizeFastaSummary");
        exit(1);
    }
}
//**********************************************************************************************************************
int Summary::driverContigsSummarySummarize(string fastafile, linePair lines) {
    try {
        ifstream in;
        m->openInputFile(fastafile, in);
        
        in.seekg(lines.start);
        
        //print header if you are process 0
        if (lines.start == 0) { m->zapGremlins(in); m->getline(in); m->gobble(in); }
        
        bool done = false;
        long long count = 0;
        string name;
        int length, OLength, thisOStart, thisOEnd, numMisMatches, numns; //Name	Length	Overlap_Length	Overlap_Start	Overlap_End	MisMatches	Num_Ns
        
        while (!done) {
            
            if (m->control_pressed) { in.close(); return 1; }
            
            //seqname	start	end	nbases	ambigs	polymer	numSeqs
            in >> name >> length >> OLength >> thisOStart >> thisOEnd >> numMisMatches >> numns; m->gobble(in);
            
            if (m->debug) { m->mothurOut("[DEBUG]: " + name + "\t" + toString(thisOStart) + "\t" + toString(thisOEnd) + "\t" + toString(length) + "\n"); }
            
            if (name != "") {
                string seqInfo = addSeq(name, length, OLength, thisOStart, thisOEnd, numMisMatches, numns); count++;
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
        m->errorOut(e, "Summary", "driverContigsSummarySummarize");
        exit(1);
    }
}
/**********************************************************************************************************************/
long long Summary::summarizeAlignSummary(string summaryfile) {
    try {
        
        vector<unsigned long long> positions;
        vector<linePair> lines;
        string p = m->getProcessors();
        int processors = 1; m->mothurConvert(p, processors);
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        positions = m->divideFile(summaryfile, processors);
        for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(linePair(positions[i], positions[(i+1)]));	}
#else
        positions = m->setFilePosFasta(summaryfile, numSeqs);
        if (numSeqs < processors) { processors = numSeqs; }
        
        //figure out how many sequences you have to process
        int numSeqsPerProcessor = numSeqs / processors;
        for (int i = 0; i < processors; i++) {
            int startIndex =  i * numSeqsPerProcessor;
            if(i == (processors - 1)){	numSeqsPerProcessor = numSeqs - i * numSeqsPerProcessor; 	}
            lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
        }
#endif
        
        int process = 1;
        long long num = 0;
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
                num = driverAlignSummarySummarize(summaryfile, lines[process]);
                
                //pass numSeqs to parent
                ofstream out;
                string tempFile = m->mothurGetpid(process) + ".num.temp";
                m->openOutputFile(tempFile, out);
                
                out << num << endl;
                out << total << endl;
                out << sims.size() << endl;
                for (map<float,  long long>::iterator it = sims.begin(); it != sims.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                out << scores.size() << endl;
                for (map<float,  long long>::iterator it = scores.begin(); it != scores.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                out << inserts.size() << endl;
                for (map<int,  long long>::iterator it = inserts.begin(); it != inserts.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                out << seqLength.size() << endl;
                for (map<int,  long long>::iterator it = seqLength.begin(); it != seqLength.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
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
            vector<unsigned long long> positions = m->divideFile(summaryfile, processors);
            for (int i = 0; i < (positions.size()-1); i++) {  lines.push_back(linePair(positions[i], positions[(i+1)]));  }
            
            sims.clear();
            scores.clear();
            inserts.clear();
            
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
                    
                    num = driverContigsSummarySummarize(summaryfile, lines[process]);
                    
                    //pass numSeqs to parent
                    ofstream out;
                    string tempFile = m->mothurGetpid(process) + ".num.temp";
                    m->openOutputFile(tempFile, out);
                    
                    out << num << endl;
                    out << total << endl;
                    out << sims.size() << endl;
                    for (map<float,  long long>::iterator it = sims.begin(); it != sims.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                    out << scores.size() << endl;
                    for (map<float,  long long>::iterator it = scores.begin(); it != scores.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                    out << inserts.size() << endl;
                    for (map<int,  long long>::iterator it = inserts.begin(); it != inserts.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                    out.close();
                    out << seqLength.size() << endl;
                    for (map<int,  long long>::iterator it = seqLength.begin(); it != seqLength.end(); it++)		{		out << it->first << '\t' << it->second << endl; }
                    
                    exit(0);
                }else {
                    m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine();
                    for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                    exit(0);
                }
            }
        }
        
        //do your part
        num = driverAlignSummarySummarize(summaryfile, lines[0]);
        
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
                float first; long long second;
                in >> first; m->gobble(in); in >> second; m->gobble(in);
                map<float,  long long>::iterator it = sims.find(first);
                if (it == sims.end()) { sims[first] = second; } //first finding of this start position, set count.
                else { it->second += second; } //add counts
            }
            m->gobble(in);
            in >> tempNum; m->gobble(in);
            for (int k = 0; k < tempNum; k++)			{
                float first; long long second;
                in >> first; m->gobble(in); in >> second; m->gobble(in);
                map<float,  long long>::iterator it = scores.find(first);
                if (it == scores.end()) { scores[first] = second; } //first finding of this end position, set count.
                else { it->second += second; } //add counts
            }
            m->gobble(in);
            in >> tempNum; m->gobble(in);
            for (int k = 0; k < tempNum; k++)			{
                long long first, second;
                in >> first; m->gobble(in); in >> second; m->gobble(in);
                map<int,  long long>::iterator it = inserts.find(first);
                if (it == inserts.end()) { inserts[first] = second; } //first finding of this end position, set count.
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
            seqSumData* tempSum = new seqSumData(summaryfile, m, lines[i]->start, lines[i]->end, hasNameOrCount, nameMap);
            pDataArray.push_back(tempSum);
            
            //MySeqSumThreadFunction is in header. It must be global or static to work with the threads.
            //default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
            hThreadArray[i] = CreateThread(NULL, 0, MySeqFastaSumThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);
        }
        
        //do your part
        num = driverAlignSummarySummarize(summaryfile, lines[processors-1]);
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
            for (map<float, long long>::iterator it = pDataArray[i]->sims.begin(); it != pDataArray[i]->sims.end(); it++)		{
                map<int, long long>::iterator itMain = sims.find(it->first);
                if (itMain == sims.end()) { //newValue
                    sims[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<float, long long>::iterator it = pDataArray[i]->scores.begin(); it != pDataArray[i]->scores.end(); it++)		{
                map<int, long long>::iterator itMain = scores.find(it->first);
                if (itMain == scores.end()) { //newValue
                    scores[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = pDataArray[i]->inserts.begin(); it != pDataArray[i]->inserts.end(); it++)		{
                map<int, long long>::iterator itMain = inserts.find(it->first);
                if (itMain == inserts.end()) { //newValue
                    inserts[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = pDataArray[i]->seqLength.begin(); it != pDataArray[i]->seqLength.end(); it++)		{
                map<int, long long>::iterator itMain = seqLength.find(it->first);
                if (itMain == seqLength.end()) { //newValue
                    seqLength[it->first] = it->second;
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
        
        numUniques = num;
        
        return num; //number of uniques
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "summarizeAlignSummary");
        exit(1);
    }
}
//**********************************************************************************************************************
int Summary::driverAlignSummarySummarize(string fastafile, linePair lines) {
    try {
        ifstream in;
        m->openInputFile(fastafile, in);
        
        in.seekg(lines.start);
        
        //print header if you are process 0
        if (lines.start == 0) { m->zapGremlins(in); m->getline(in); m->gobble(in); }
        
        bool done = false;
        long long count = 0;
        string name, TemplateName, SearchMethod, AlignmentMethod;
        int length, TemplateLength,	 QueryStart,	QueryEnd,	TemplateStart,	TemplateEnd,	PairwiseAlignmentLength,	GapsInQuery,	GapsInTemplate,	LongestInsert;
        float SearchScore, SimBtwnQueryTemplate;  //QueryName	QueryLength	TemplateName	TemplateLength	SearchMethod	SearchScore	AlignmentMethod	QueryStart	QueryEnd	TemplateStart	TemplateEnd	PairwiseAlignmentLength	GapsInQuery	GapsInTemplate	LongestInsert	SimBtwnQuery&Template
        //checking for minScore, maxInsert, minSim
        
        while (!done) {
            
            if (m->control_pressed) { in.close(); return 1; }
            
            in >> name >> length >> TemplateName >> TemplateLength >> SearchMethod >> SearchScore >> AlignmentMethod >> QueryStart >> QueryEnd >> TemplateStart >> TemplateEnd >> PairwiseAlignmentLength >> GapsInQuery >> GapsInTemplate >> LongestInsert >> SimBtwnQueryTemplate; m->gobble(in);

            
            if (m->debug) { m->mothurOut("[DEBUG]: " + name + "\t" + toString(TemplateName) + "\t" + toString(SearchScore) + "\t" + toString(length) + "\n"); }
            
            if (name != "") {
                string seqInfo = addSeq(name, length, SimBtwnQueryTemplate, SearchScore, LongestInsert); count++;
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
        m->errorOut(e, "Summary", "driverAlignSummarySummarize");
        exit(1);
    }
}
//**********************************************************************************************************************







