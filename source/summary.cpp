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
            if (isCountFile(n)) {
                CountTable ct;
                ct.readTable(n, false, false);
                nameMap = ct.getNameMap();
                type = "count";
            }else { Utils util; nameMap = util.readNames(n); type = "name"; }
        }
        nameCountNumUniques = nameMap.size();
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "Summary");
        exit(1);
    }
}
//**********************************************************************************************************************
bool Summary::isCountFile(string inputfile){
    try {
        CountTable ct;
        bool isCount = ct.isCountTable(inputfile);
        return isCount;
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "isCountFile");
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
long long Summary::getValue(map<int, long long>& spots, double value) {
    try {
        long long percentage = 1+(long long)(total * value * 0.01);
        long long result = 0;
        long long totalSoFar = 0;
        long long lastValue = 0;

        //minimum
        if ((spots.begin())->first == -1) { result = 0; }
        else {result = (spots.begin())->first; }

        for (it = spots.begin(); it != spots.end(); it++) {
            long long value = it->first; if (value == -1) { value = 0; }
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
        if (util.isEqual((positions.begin())->first, -1)) { results[0] = 0; }
        else {results[0] = (positions.begin())->first; }
        results[1] = results[0]; results[2] = results[0]; results[3] = results[0]; results[4] = results[0]; results[5] = results[0];

        for (map<float, long long>::iterator it = positions.begin(); it != positions.end(); it++) {
            long long value = it->first; if (value == -1) { value = 0; }
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
        if (util.isEqual((positions.begin())->first, -1)) { result = 0; }
        else { result = (positions.begin())->first; }

        for (map<float, long long>::iterator it = positions.begin(); it != positions.end(); it++) {
            long long value = it->first; if (value == -1) { value = 0; }
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

//**************************************************************************************************

int Summary::getMaxAbundance(){

	int max = 0;

	for(map<string,int>::iterator it=nameMap.begin();it!=nameMap.end();it++){
		if(it->second > max){
			max = it->second;
		}
	}

	return max;

}

//**************************************************************************************************

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
void driverSummarize(seqSumData* params) { //(string fastafile, string output, linePair lines) {
    try {
        ofstream out;
        if (params->summaryFile != "") { params->util.openOutputFile(params->summaryFile, out); }

        ifstream in;
        params->util.openInputFile(params->filename, in);

        in.seekg(params->start);

        //print header if you are process 0
        if (params->start == 0) {
            params->util.zapGremlins(in); params->util.gobble(in);
            //print header if you are process 0
            if (params->summaryFile != "") { out << "seqname\tstart\tend\tnbases\tambigs\tpolymer\tnumSeqs" << endl; }
        }

        bool done = false;
        params->count = 0;

        while (!done) {

            if (params->m->getControl_pressed()) {  break; }

            Sequence seq(in); params->util.gobble(in);

            if (seq.getName() != "") {

                if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: " + seq.getName() + "\t" + toString(seq.getStartPos()) + "\t" + toString(seq.getEndPos()) + "\t" + toString(seq.getNumBases()) + "\n"); }

                //string seqInfo = addSeq(current);
                params->count++;

                long long num = 1;

                if (params->hasNameMap) {
                    //make sure this sequence is in the namefile, else error
                    map<string, int>::iterator itFindName = params->nameMap.find(seq.getName());

                    if (itFindName == params->nameMap.end()) { params->m->mothurOut("[ERROR]: '" + seq.getName() + "' is not in your name or count file, please correct."); params->m->mothurOutEndLine(); params->m->setControl_pressed(true); }
                    else { num = itFindName->second; }
                }

                int thisStartPosition = seq.getStartPos();
                map<int, long long>::iterator it = params->startPosition.find(thisStartPosition);
                if (it == params->startPosition.end()) { params->startPosition[thisStartPosition] = num; } //first finding of this start position, set count.
                else { it->second += num; } //add counts

                int thisEndPosition = seq.getEndPos();
                it = params->endPosition.find(thisEndPosition);
                if (it == params->endPosition.end()) { params->endPosition[thisEndPosition] = num; } //first finding of this end position, set count.
                else { it->second += num; } //add counts

                int thisSeqLength = seq.getNumBases();
                it = params->seqLength.find(thisSeqLength);
                if (it == params->seqLength.end()) { params->seqLength[thisSeqLength] = num; } //first finding of this length, set count.
                else { it->second += num; } //add counts

                int thisAmbig = seq.getAmbigBases();
                it = params->ambigBases.find(thisAmbig);
                if (it == params->ambigBases.end()) { params->ambigBases[thisAmbig] = num; } //first finding of this ambig, set count.
                else { it->second += num; } //add counts

                int thisHomoP = seq.getLongHomoPolymer();
                it = params->longHomoPolymer.find(thisHomoP);
                if (it == params->longHomoPolymer.end()) { params->longHomoPolymer[thisHomoP] = num; } //first finding of this homop, set count.
                else { it->second += num; } //add counts

                int numns = seq.getNumNs();
                it = params->numNs.find(numns);
                if (it == params->numNs.end()) { params->numNs[numns] = num; } //first finding of this homop, set count.
                else { it->second += num; } //add counts

                params->total += num;

                string seqInfo = "";
                seqInfo += seq.getName() + '\t';
                seqInfo += toString(thisStartPosition) + '\t' + toString(thisEndPosition) + '\t';
                seqInfo += toString(thisSeqLength) + '\t' + toString(thisAmbig) + '\t';
                seqInfo += toString(thisHomoP) + '\t' + toString(num);

                if (params->summaryFile != "") { out << seqInfo << endl; }
            }

#if defined NON_WINDOWS
            unsigned long long pos = in.tellg();
            if ((pos == -1) || (pos >= params->end)) { break; }
#else
            if (params->count == params->end) { break; }
#endif
        }

        if (params->summaryFile != "") { out.close(); }
        in.close();

    }
    catch(exception& e) {
        params->m->errorOut(e, "Summary", "driverSummarize");
        exit(1);
    }
}
//**********************************************************************************************************************
long long Summary::summarizeFasta(string fastafile, string output) {
    try {
        long long num = 0;
        vector<linePair> lines;
        vector<double> positions;
#if defined NON_WINDOWS
        positions = util.divideFile(fastafile, processors);
        for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(linePair(positions[i], positions[(i+1)]));	 }
#else
        positions = util.setFilePosFasta(fastafile, num);
        if (num < processors) { processors = num; }

        //figure out how many sequences you have to process
        int numSeqsPerProcessor = num / processors;
        for (int i = 0; i < processors; i++) {
            int startIndex =  i * numSeqsPerProcessor;
            if(i == (processors - 1)){	numSeqsPerProcessor = num - i * numSeqsPerProcessor; 	}
            lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
        }
#endif

        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<seqSumData*> data;

        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            string extension = "";
            extension = toString(i) + ".temp";
            string outputName = output + extension;
            if (output == "") {  outputName = "";  }

            seqSumData* dataBundle = new seqSumData(fastafile, outputName, lines[i+1].start, lines[i+1].end, hasNameOrCount, nameMap);
            data.push_back(dataBundle);

            workerThreads.push_back(new std::thread(driverSummarize, dataBundle));
        }

        seqSumData* dataBundle = new seqSumData(fastafile, output, lines[0].start, lines[0].end, hasNameOrCount, nameMap);

        driverSummarize(dataBundle);
        num = dataBundle->count;
        total = dataBundle->total;
        startPosition = dataBundle->startPosition;
        endPosition = dataBundle->endPosition;
        seqLength = dataBundle->seqLength;
        ambigBases = dataBundle->ambigBases;
        longHomoPolymer = dataBundle->longHomoPolymer;
        numNs = dataBundle->numNs;
        delete dataBundle;

        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            num += data[i]->count;
            total += data[i]->total;

            for (map<int, long long>::iterator it = data[i]->startPosition.begin(); it != data[i]->startPosition.end(); it++)		{
                map<int, long long>::iterator itMain = startPosition.find(it->first);
                if (itMain == startPosition.end()) { //newValue
                    startPosition[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = data[i]->endPosition.begin(); it != data[i]->endPosition.end(); it++)		{
                map<int, long long>::iterator itMain = endPosition.find(it->first);
                if (itMain == endPosition.end()) { //newValue
                    endPosition[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = data[i]->seqLength.begin(); it != data[i]->seqLength.end(); it++)		{
                map<int, long long>::iterator itMain = seqLength.find(it->first);
                if (itMain == seqLength.end()) { //newValue
                    seqLength[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = data[i]->ambigBases.begin(); it != data[i]->ambigBases.end(); it++)		{
                map<int, long long>::iterator itMain = ambigBases.find(it->first);
                if (itMain == ambigBases.end()) { //newValue
                    ambigBases[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = data[i]->longHomoPolymer.begin(); it != data[i]->longHomoPolymer.end(); it++)		{
                map<int, long long>::iterator itMain = longHomoPolymer.find(it->first);
                if (itMain == longHomoPolymer.end()) { //newValue
                    longHomoPolymer[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = data[i]->numNs.begin(); it != data[i]->numNs.end(); it++)		{
                map<int, long long>::iterator itMain = numNs.find(it->first);
                if (itMain == numNs.end()) { //newValue
                    numNs[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }

            delete data[i];
            delete workerThreads[i];
        }

        //append files
        for (int i = 0; i < processors-1; i++) {
            string extension = "";
            extension = toString(i) + ".temp";
            string outputName = output + extension;
            if (output == "") {  outputName = "";  }

            if (outputName != "") {
                util.appendFiles((output + toString(i) + ".temp"), output);
                util.mothurRemove((output + toString(i) + ".temp"));
            }
        }

        if (hasNameOrCount) {
            if (nameCountNumUniques != num) { // do fasta and name/count files match
                m->mothurOut("[ERROR]: Your " + type + " file contains " + toString(nameCountNumUniques) + " unique sequences, but your fasta file contains " + toString(num) + ". File mismatch detected, quitting command.\n"); m->setControl_pressed(true);
            }
        }

        numUniques = num;

        return num;

    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "summarizeFasta");
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
//**********************************************************************************************************************
void driverFastaSummarySummarize(seqSumData* params) {
    try {
        ifstream in;
        params->util.openInputFile(params->filename, in);

        in.seekg(params->start);

        //print header if you are process 0
        if (params->start == 0) { params->util.zapGremlins(in); params->util.getline(in); params->util.gobble(in); params->count++; }

        bool done = false;
        string name;
        int start, end, length, ambigs, polymer;
        long long numReps;

        while (!done) {

            if (params->m->getControl_pressed()) {  break; }

            //seqname	start	end	nbases	ambigs	polymer	numSeqs
            in >> name >> start >> end >> length >> ambigs >> polymer >> numReps; params->util.gobble(in);

            if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: " + name + "\t" + toString(start) + "\t" + toString(end) + "\t" + toString(length) + "\n"); }

            if (name != "") {
                if ((numReps == 1) && params->hasNameMap) {
                    //make sure this sequence is in the namefile, else error
                    map<string, int>::iterator itFindName = params->nameMap.find(name);

                    if (itFindName == params->nameMap.end()) { params->m->mothurOut("[ERROR]: '" + name + "' is not in your name or count file, please correct."); params->m->mothurOutEndLine(); params->m->setControl_pressed(true); }
                    else { numReps = itFindName->second; }
                }

                map<int, long long>::iterator it = params->startPosition.find(start);
                if (it == params->startPosition.end()) { params->startPosition[start] = numReps; } //first finding of this start position, set count.
                else { it->second += numReps; } //add counts

                it = params->endPosition.find(end);
                if (it == params->endPosition.end()) { params->endPosition[end] = numReps; } //first finding of this end position, set count.
                else { it->second += numReps; } //add counts

                it = params->seqLength.find(length);
                if (it == params->seqLength.end()) { params->seqLength[length] = numReps; } //first finding of this length, set count.
                else { it->second += numReps; } //add counts

                it = params->ambigBases.find(ambigs);
                if (it == params->ambigBases.end()) { params->ambigBases[ambigs] = numReps; } //first finding of this ambig, set count.
                else { it->second += numReps; } //add counts

                it = params->longHomoPolymer.find(polymer);
                if (it == params->longHomoPolymer.end()) { params->longHomoPolymer[polymer] = numReps; } //first finding of this homop, set count.
                else { it->second += numReps; } //add counts

                params->count++;
                params->total += numReps;
            }

#if defined NON_WINDOWS
            unsigned long long pos = in.tellg();
            if ((pos == -1) || (pos >= params->end)) { break; }
#else
            if (params->end == params->count) { break; }
#endif
        }

        in.close();

    }
    catch(exception& e) {
        params-> m->errorOut(e, "Summary", "driverFastaSummarySummarize");
        exit(1);
    }
}
/**********************************************************************************************************************/
long long Summary::summarizeFastaSummary(string summaryfile) {
    try {
        long long num = 0;
        vector<double> positions;
        vector<linePair> lines;
#if defined NON_WINDOWS
        positions = util.divideFilePerLine(summaryfile, processors);
        for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(linePair(positions[i], positions[(i+1)]));	}
#else
        positions = util.setFilePosEachLine(summaryfile, num);
        if (num < processors) { processors = num; }

        //figure out how many sequences you have to process
        int numSeqsPerProcessor = num / processors;
        for (int i = 0; i < processors; i++) {
            int startIndex =  i * numSeqsPerProcessor;
            if(i == (processors - 1)){	numSeqsPerProcessor = num - i * numSeqsPerProcessor; 	}
            lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
        }
#endif

        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<seqSumData*> data;

        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            seqSumData* dataBundle = new seqSumData(summaryfile, lines[i+1].start, lines[i+1].end, hasNameOrCount, nameMap);
            data.push_back(dataBundle);

            workerThreads.push_back(new std::thread(driverFastaSummarySummarize, dataBundle));
        }

        seqSumData* dataBundle = new seqSumData(summaryfile, lines[0].start, lines[0].end, hasNameOrCount, nameMap);

        driverFastaSummarySummarize(dataBundle);
        num = dataBundle->count-1; //header line
        total = dataBundle->total;
        startPosition = dataBundle->startPosition;
        endPosition = dataBundle->endPosition;
        seqLength = dataBundle->seqLength;
        ambigBases = dataBundle->ambigBases;
        longHomoPolymer = dataBundle->longHomoPolymer;
        delete dataBundle;

        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            num += data[i]->count;
            total += data[i]->total;

            for (map<int, long long>::iterator it = data[i]->startPosition.begin(); it != data[i]->startPosition.end(); it++)		{
                map<int, long long>::iterator itMain = startPosition.find(it->first);
                if (itMain == startPosition.end()) { //newValue
                    startPosition[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = data[i]->endPosition.begin(); it != data[i]->endPosition.end(); it++)		{
                map<int, long long>::iterator itMain = endPosition.find(it->first);
                if (itMain == endPosition.end()) { //newValue
                    endPosition[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = data[i]->seqLength.begin(); it != data[i]->seqLength.end(); it++)		{
                map<int, long long>::iterator itMain = seqLength.find(it->first);
                if (itMain == seqLength.end()) { //newValue
                    seqLength[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = data[i]->ambigBases.begin(); it != data[i]->ambigBases.end(); it++)		{
                map<int, long long>::iterator itMain = ambigBases.find(it->first);
                if (itMain == ambigBases.end()) { //newValue
                    ambigBases[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = data[i]->longHomoPolymer.begin(); it != data[i]->longHomoPolymer.end(); it++)		{
                map<int, long long>::iterator itMain = longHomoPolymer.find(it->first);
                if (itMain == longHomoPolymer.end()) { //newValue
                    longHomoPolymer[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }

            delete data[i];
            delete workerThreads[i];
        }


        if (hasNameOrCount) {
            if (nameCountNumUniques != num) { // do fasta and name/count files match
                m->mothurOut("[ERROR]: Your " + type + " file contains " + toString(nameCountNumUniques) + " unique sequences, but your fasta file contains " + toString(num) + ". File mismatch detected, quitting command.\n"); m->setControl_pressed(true);
            }
        }

        numUniques = num;

        return num;

    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "summarizeFastaSummary");
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
//**********************************************************************************************************************
void driverContigsSummarySummarize(seqSumData* params) {
    try {
        ifstream in;
        params->util.openInputFile(params->filename, in);

        in.seekg(params->start);

        //print header if you are process 0
        if (params->start == 0) { params->util.zapGremlins(in); params->util.getline(in); params->util.gobble(in); params->count++; }

        bool done = false;
        string name;
        int length, OLength, thisOStart, thisOEnd, numMisMatches, numns; //Name	Length	Overlap_Length	Overlap_Start	Overlap_End	MisMatches	Num_Ns
        double expectedErrors;

        while (!done) {

            if (params->m->getControl_pressed()) { break; }

            //seqname	start	end	nbases	ambigs	polymer	numSeqs
            in >> name >> length >> OLength >> thisOStart >> thisOEnd >> numMisMatches >> numns >> expectedErrors; params->util.gobble(in);

            if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: " + name + "\t" + toString(thisOStart) + "\t" + toString(thisOEnd) + "\t" + toString(length) + "\n"); }

            if (name != "") {
                long long numReps = 1;
                if (params->hasNameMap) {
                    //make sure this sequence is in the namefile, else error
                    map<string, int>::iterator itFindName = params->nameMap.find(name);

                    if (itFindName == params->nameMap.end()) { params->m->mothurOut("[ERROR]: '" + name + "' is not in your name or count file, please correct."); params->m->mothurOutEndLine(); params->m->setControl_pressed(true); }
                    else { numReps = itFindName->second; }
                }

                map<int, long long>::iterator it = params->ostartPosition.find(thisOStart);
                if (it == params->ostartPosition.end()) { params->ostartPosition[thisOStart] = numReps; } //first finding of this start position, set count.
                else { it->second += numReps; } //add counts

                it = params->oendPosition.find(thisOEnd);
                if (it == params->oendPosition.end()) { params->oendPosition[thisOEnd] = numReps; } //first finding of this end position, set count.
                else { it->second += numReps; } //add counts

                it = params->oseqLength.find(OLength);
                if (it == params->oseqLength.end()) { params->oseqLength[OLength] = numReps; } //first finding of this length, set count.
                else { it->second += numReps; } //add counts

                it = params->seqLength.find(length);
                if (it == params->seqLength.end()) { params->seqLength[length] = numReps; } //first finding of this length, set count.
                else { it->second += numReps; } //add counts

                it = params->misMatches.find(numMisMatches);
                if (it == params->misMatches.end()) { params->misMatches[numMisMatches] = numReps; } //first finding of this ambig, set count.
                else { it->second += numReps; } //add counts

                it = params->numNs.find(numns);
                if (it == params->numNs.end()) { params->numNs[numns] = numReps; } //first finding of this homop, set count.
                else { it->second += numReps; } //add counts

                params->count++;
                params->total += numReps;
            }

#if defined NON_WINDOWS
            unsigned long long pos = in.tellg();
            if ((pos == -1) || (pos >= params->end)) { break; }
#else
            if (params->end == params->count) { break; }
#endif
        }

        in.close();
    }
    catch(exception& e) {
        params->m->errorOut(e, "Summary", "driverContigsSummarySummarize");
        exit(1);
    }
}
/**********************************************************************************************************************/
long long Summary::summarizeContigsSummary(string summaryfile) {
    try {
        long long num = 0;
        vector<double> positions;
        vector<linePair> lines;
#if defined NON_WINDOWS
        positions = util.divideFilePerLine(summaryfile, processors);
        for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(linePair(positions[i], positions[(i+1)]));	}
#else
        positions = util.setFilePosEachLine(summaryfile, num);
        if (num < processors) { processors = num; }

        //figure out how many sequences you have to process
        int numSeqsPerProcessor = num / processors;
        for (int i = 0; i < processors; i++) {
            int startIndex =  i * numSeqsPerProcessor;
            if(i == (processors - 1)){	numSeqsPerProcessor = num - i * numSeqsPerProcessor; 	}
            lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
        }
#endif

        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<seqSumData*> data;

        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            if (m->getDebug()) { m->mothurOut("[DEBUG]: creating thread " + toString(i+1) + "\n"); }
            seqSumData* dataBundle = new seqSumData(summaryfile, lines[i+1].start, lines[i+1].end, hasNameOrCount, nameMap);
            data.push_back(dataBundle);

            workerThreads.push_back(new std::thread(driverContigsSummarySummarize, dataBundle));
        }

        seqSumData* dataBundle = new seqSumData(summaryfile, lines[0].start, lines[0].end, hasNameOrCount, nameMap);

        driverContigsSummarySummarize(dataBundle);
        num = dataBundle->count-1; //header line
        total = dataBundle->total;
        ostartPosition = dataBundle->ostartPosition;
        oendPosition = dataBundle->oendPosition;
        seqLength = dataBundle->seqLength;
        oseqLength = dataBundle->oseqLength;
        misMatches = dataBundle->misMatches;
        numNs = dataBundle->numNs;
        delete dataBundle;


        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            num += data[i]->count;
            total += data[i]->total;


            for (map<int, long long>::iterator it = data[i]->ostartPosition.begin(); it != data[i]->ostartPosition.end(); it++)		{
                map<int, long long>::iterator itMain = ostartPosition.find(it->first);
                if (itMain == ostartPosition.end()) { //newValue
                    ostartPosition[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = data[i]->oendPosition.begin(); it != data[i]->oendPosition.end(); it++)		{
                map<int, long long>::iterator itMain = oendPosition.find(it->first);
                if (itMain == oendPosition.end()) { //newValue
                    oendPosition[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = data[i]->oseqLength.begin(); it != data[i]->oseqLength.end(); it++)		{
                map<int, long long>::iterator itMain = oseqLength.find(it->first);
                if (itMain == oseqLength.end()) { //newValue
                    oseqLength[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = data[i]->seqLength.begin(); it != data[i]->seqLength.end(); it++)		{
                map<int, long long>::iterator itMain = seqLength.find(it->first);
                if (itMain == seqLength.end()) { //newValue
                    seqLength[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = data[i]->misMatches.begin(); it != data[i]->misMatches.end(); it++)		{
                map<int, long long>::iterator itMain = misMatches.find(it->first);
                if (itMain == misMatches.end()) { //newValue
                    misMatches[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = data[i]->numNs.begin(); it != data[i]->numNs.end(); it++)		{
                map<int, long long>::iterator itMain = numNs.find(it->first);
                if (itMain == numNs.end()) { //newValue
                    numNs[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }

            delete data[i];
            delete workerThreads[i];
        }


        if (hasNameOrCount) {
            if (nameCountNumUniques != num) { // do fasta and name/count files match
                m->mothurOut("[ERROR]: Your " + type + " file contains " + toString(nameCountNumUniques) + " unique sequences, but your fasta file contains " + toString(num) + ". File mismatch detected, quitting command.\n"); m->setControl_pressed(true);
            }
        }

        numUniques = num;

        return num;
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
void driverAlignSummarySummarize(seqSumData* params) {
    try {
        ifstream in;
        params->util.openInputFile(params->filename, in);

        in.seekg(params->start);

        //print header if you are process 0
        if (params->start == 0) { params->util.zapGremlins(in); params->util.getline(in); params->util.gobble(in); params->count++; }

        bool done = false;
        string name, TemplateName, SearchMethod, AlignmentMethod;
        int length, TemplateLength,	 QueryStart,	QueryEnd,	TemplateStart,	TemplateEnd,	PairwiseAlignmentLength,	GapsInQuery,	GapsInTemplate,	LongestInsert;
        float SearchScore, SimBtwnQueryTemplate;

        while (!done) {

            if (params->m->getControl_pressed()) {  break; }

            in >> name >> length >> TemplateName >> TemplateLength >> SearchMethod >> SearchScore >> AlignmentMethod >> QueryStart >> QueryEnd >> TemplateStart >> TemplateEnd >> PairwiseAlignmentLength >> GapsInQuery >> GapsInTemplate >> LongestInsert >> SimBtwnQueryTemplate; params->util.gobble(in);


            if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: " + name + "\t" + toString(TemplateName) + "\t" + toString(SearchScore) + "\t" + toString(length) + "\n"); }

            if (name != "") {
                //string seqInfo = addSeq(name, length, SimBtwnQueryTemplate, SearchScore, LongestInsert);
                long long numReps = 1;
                if (params->hasNameMap) {
                    //make sure this sequence is in the namefile, else error
                     map<string, int>::iterator itFindName = params->nameMap.find(name);

                    if (itFindName == params->nameMap.end()) { params->m->mothurOut("[ERROR]: '" + name + "' is not in your name or count file, please correct."); params->m->mothurOutEndLine(); params->m->setControl_pressed(true); }
                    else { numReps = itFindName->second; }
                }

                map<float, long long>:: iterator itFloat = params->sims.find(SimBtwnQueryTemplate);
                if (itFloat == params->sims.end()) { params->sims[SimBtwnQueryTemplate] = numReps; } //first finding of this similarity score, set count.
                else { itFloat->second += numReps; } //add counts

                itFloat = params->scores.find(SearchScore);
                if (itFloat == params->scores.end()) { params->scores[SearchScore] = numReps; } //first finding of this end position, set count.
                else { itFloat->second += numReps; } //add counts

                map<int, long long>::iterator it = params->inserts.find(LongestInsert);
                if (it == params->inserts.end()) { params->inserts[LongestInsert] = numReps; } //first finding of this length, set count.
                else { it->second += numReps; } //add counts

                it = params->seqLength.find(length);
                if (it == params->seqLength.end()) { params->seqLength[length] = numReps; } //first finding of this length, set count.
                else { it->second += numReps; } //add counts

                params->count++;
                params->total += numReps;

            }

#if defined NON_WINDOWS
            unsigned long long pos = in.tellg();
            if ((pos == -1) || (pos >= params->end)) { break; }
#else
            if (params->end == params->count) { break; }
#endif
        }

        in.close();
    }
    catch(exception& e) {
        params->m->errorOut(e, "Summary", "driverAlignSummarySummarize");
        exit(1);
    }
}
/**********************************************************************************************************************/
long long Summary::summarizeAlignSummary(string summaryfile) {
    try {
        long long num = 0;
        vector<double> positions;
        vector<linePair> lines;
#if defined NON_WINDOWS
        positions = util.divideFilePerLine(summaryfile, processors);
        for (int i = 0; i < (positions.size()-1); i++) {	lines.push_back(linePair(positions[i], positions[(i+1)]));	}
#else
        positions = util.setFilePosEachLine(summaryfile, num);
        if (num < processors) { processors = num; }

        //figure out how many sequences you have to process
        int numSeqsPerProcessor = num / processors;
        for (int i = 0; i < processors; i++) {
            int startIndex =  i * numSeqsPerProcessor;
            if(i == (processors - 1)){	numSeqsPerProcessor = num - i * numSeqsPerProcessor; 	}
            lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
        }
#endif
        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<seqSumData*> data;

        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {

            seqSumData* dataBundle = new seqSumData(summaryfile, lines[i+1].start, lines[i+1].end, hasNameOrCount, nameMap);
            data.push_back(dataBundle);

            workerThreads.push_back(new std::thread(driverAlignSummarySummarize, dataBundle));
        }

        seqSumData* dataBundle = new seqSumData(summaryfile, lines[0].start, lines[0].end, hasNameOrCount, nameMap);

        driverAlignSummarySummarize(dataBundle);
        num = dataBundle->count-1; //header line
        total = dataBundle->total;
        sims = dataBundle->sims;
        scores = dataBundle->scores;
        inserts = dataBundle->inserts;
        seqLength = dataBundle->seqLength;
        delete dataBundle;

        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            num += data[i]->count;
            total += data[i]->total;

            for (map<float, long long>::iterator it = data[i]->sims.begin(); it != data[i]->sims.end(); it++)		{
                map<float, long long>::iterator itMain = sims.find(it->first);
                if (itMain == sims.end()) { //newValue
                    sims[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<float, long long>::iterator it = data[i]->scores.begin(); it != data[i]->scores.end(); it++)		{
                map<float, long long>::iterator itMain = scores.find(it->first);
                if (itMain == scores.end()) { //newValue
                    scores[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = data[i]->inserts.begin(); it != data[i]->inserts.end(); it++)		{
                map<int, long long>::iterator itMain = inserts.find(it->first);
                if (itMain == inserts.end()) { //newValue
                    inserts[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }
            for (map<int, long long>::iterator it = data[i]->seqLength.begin(); it != data[i]->seqLength.end(); it++)		{
                map<int, long long>::iterator itMain = seqLength.find(it->first);
                if (itMain == seqLength.end()) { //newValue
                    seqLength[it->first] = it->second;
                }else { itMain->second += it->second; } //merge counts
            }

            delete data[i];
            delete workerThreads[i];
        }


        if (hasNameOrCount) {
            if (nameCountNumUniques != num) { // do fasta and name/count files match
                m->mothurOut("[ERROR]: Your " + type + " file contains " + toString(nameCountNumUniques) + " unique sequences, but your fasta file contains " + toString(num) + ". File mismatch detected, quitting command.\n"); m->setControl_pressed(true);
            }
        }

        numUniques = num;

        return num;
    }
    catch(exception& e) {
        m->errorOut(e, "Summary", "summarizeAlignSummary");
        exit(1);
    }
}
//**********************************************************************************************************************
