/*
 *  qualityscores.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 7/12/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "qualityscores.h"

/**************************************************************************************************/

QualityScores::QualityScores(){
	try {
		m = MothurOut::getInstance();
		seqName = "";
		seqLength = -1;

	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "QualityScores");
		exit(1);
	}
}
/**************************************************************************************************/

QualityScores::QualityScores(string n, vector<int> s){
	try {
		m = MothurOut::getInstance();
		setName(n);
        setScores(s);
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "QualityScores");
		exit(1);
	}
}
/**************************************************************************************************/

QualityScores::QualityScores(ifstream& qFile){
	try {

		m = MothurOut::getInstance();

		int score;
        seqName = getSequenceName(qFile); util.gobble(qFile); getCommentString(qFile);

        if (m->getDebug()) { m->mothurOut("[DEBUG]: name = '" + seqName + "'\n.");  }

		if (!m->getControl_pressed()) {
            string qScoreString = util.getline(qFile); util.gobble(qFile);

            if (m->getDebug()) { m->mothurOut("[DEBUG]: scores = '" + qScoreString + "'\n.");  }

            while(qFile.peek() != '>' && qFile.peek() != EOF){
                if (m->getControl_pressed()) { break; }
                string temp = util.getline(qFile); util.gobble(qFile);
                qScoreString +=  ' ' + temp;
            }

            istringstream qScoreStringStream(qScoreString);
            int count = 0;
            while(!qScoreStringStream.eof()){
                if (m->getControl_pressed()) { break; }
                string temp;
                qScoreStringStream >> temp;  util.gobble(qScoreStringStream);

                  //check temp to make sure its a number
                if (!util.isContainingOnlyDigits(temp)) { m->mothurOut("[ERROR]: In sequence " + seqName + "'s quality scores, expected a number and got " + temp + ", setting score to 0."); m->mothurOutEndLine(); temp = "0"; }
                convert(temp, score);

                qScores.push_back(score);
                count++;
            }
        }

		seqLength = qScores.size();

	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "QualityScores");
		exit(1);
	}

}
/**************************************************************************************************/
#ifdef USE_BOOST
QualityScores::QualityScores(boost::iostreams::filtering_istream& qFile){
    try {

        m = MothurOut::getInstance();

        int score;
        seqName = getSequenceName(qFile); util.gobble(qFile); getCommentString(qFile);

        if (m->getDebug()) { m->mothurOut("[DEBUG]: name = '" + seqName + "'\n.");  }

        if (!m->getControl_pressed()) {
            string qScoreString = util.getline(qFile); util.gobble(qFile);

            if (m->getDebug()) { m->mothurOut("[DEBUG]: scores = '" + qScoreString + "'\n.");  }

            while(qFile.peek() != '>' && qFile.peek() != EOF){
                if (m->getControl_pressed()) { break; }
                string temp = util.getline(qFile); util.gobble(qFile);
                
                qScoreString +=  ' ' + temp;
            }
           
            istringstream qScoreStringStream(qScoreString);
            int count = 0;
            while(!qScoreStringStream.eof()){
                if (m->getControl_pressed()) { break; }
                string temp;
                qScoreStringStream >> temp;  util.gobble(qScoreStringStream);

                //check temp to make sure its a number
                if (!util.isContainingOnlyDigits(temp)) { m->mothurOut("[ERROR]: In sequence " + seqName + "'s quality scores, expected a number and got " + temp + ", setting score to 0."); m->mothurOutEndLine(); temp = "0"; }
                convert(temp, score);

                
                qScores.push_back(score);
                count++;
            }
        }

        seqLength = qScores.size();

    }
    catch(exception& e) {
        m->errorOut(e, "QualityScores", "QualityScores");
        exit(1);
    }

}
#endif
/**************************************************************************************************/

int QualityScores::read(ifstream& qFile){
    try {
        int score;
        seqName = getSequenceName(qFile); util.gobble(qFile); getCommentString(qFile);

        if (m->getDebug()) { m->mothurOut("[DEBUG]: name = '" + seqName + "'\n.");  }

        if (!m->getControl_pressed()) {
            string qScoreString = util.getline(qFile); util.gobble(qFile);

            if (m->getDebug()) { m->mothurOut("[DEBUG]: scores = '" + qScoreString + "'\n.");  }

            while(qFile.peek() != '>' && qFile.peek() != EOF){
                if (m->getControl_pressed()) { break; }
                string temp = util.getline(qFile); util.gobble(qFile);
                
                qScoreString +=  ' ' + temp;
            }
            
            istringstream qScoreStringStream(qScoreString);
            int count = 0;
            while(!qScoreStringStream.eof()){
                if (m->getControl_pressed()) { break; }
                string temp;
                qScoreStringStream >> temp;  util.gobble(qScoreStringStream);

                

                //check temp to make sure its a number
                if (!util.isContainingOnlyDigits(temp)) { m->mothurOut("[ERROR]: In sequence " + seqName + "'s quality scores, expected a number and got " + temp + ", setting score to 0."); m->mothurOutEndLine(); temp = "0"; }
                convert(temp, score);

                qScores.push_back(score);
                count++;
            }
        }

        seqLength = qScores.size();

        return seqLength;

    }
    catch(exception& e) {
        m->errorOut(e, "QualityScores", "read");
        exit(1);
    }
}
//********************************************************************************************************************
string QualityScores::getSequenceName(ifstream& qFile) {
	try {
		string name = "";

        qFile >> name;

		if (name.length() != 0) {

			name = name.substr(1);

            util.checkName(name);

        }else{ m->mothurOut("Error in reading your qfile, at position " + toString(qFile.tellg()) + ". Blank name."); m->mothurOutEndLine(); m->setControl_pressed(true);  }

		return name;
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "getSequenceName");
		exit(1);
	}
}
//********************************************************************************************************************
#ifdef USE_BOOST
string QualityScores::getSequenceName(boost::iostreams::filtering_istream& qFile) {
    try {
        string name = "";

        qFile >> name; string temp;

        if (name.length() != 0) {

            name = name.substr(1);

            util.checkName(name);

        }else{ m->mothurOut("Error in reading your qfile, at position " + toString(qFile.tellg()) + ". Blank name."); m->mothurOutEndLine(); m->setControl_pressed(true);  }

        return name;
    }
    catch(exception& e) {
        m->errorOut(e, "QualityScores", "getSequenceName");
        exit(1);
    }
}
#endif
//********************************************************************************************************************
//comment can contain '>' so we need to account for that
string QualityScores::getCommentString(ifstream& fastaFile) {
    try {
        char letter;
        string temp = "";

        while(fastaFile){
            letter=fastaFile.get();
            if((letter == '\r') || (letter == '\n') || letter == -1){
                util.gobble(fastaFile);  //in case its a \r\n situation
                break;
            }else {
                temp += letter;
            }
        }

        return temp;
    }
    catch(exception& e) {
        m->errorOut(e, "QualityScores", "getCommentString");
        exit(1);
    }
}
//********************************************************************************************************************
#ifdef USE_BOOST
//comment can contain '>' so we need to account for that
string QualityScores::getCommentString(boost::iostreams::filtering_istream& fastaFile) {
    try {
        char letter;
        string temp = "";

        while(fastaFile){
            letter=fastaFile.get();
            if((letter == '\r') || (letter == '\n') || letter == -1){
                util.gobble(fastaFile);  //in case its a \r\n situation
                break;
            }else {
                temp += letter;
            }
        }

        return temp;
    }
    catch(exception& e) {
        m->errorOut(e, "QualityScores", "getCommentString");
        exit(1);
    }
}
#endif
//********************************************************************************************************************
void QualityScores::setName(string name) {
	try {

        util.checkName(name);
        seqName = name;
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "setName");
		exit(1);
	}
}
/**************************************************************************************************/

string QualityScores::getName(){

	try {
		return seqName;
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "getName");
		exit(1);
	}
}
/**************************************************************************************************/
void QualityScores::printQScores(OutputWriter* qFile){
    try {
        double expected_errors = calculateExpectedErrors();

        string outputQualString = ">";  outputQualString += seqName + '\t' + toString(expected_errors) + '\n';

        for(int i=0;i<seqLength;i++){ outputQualString += qScores[i] + ' '; } outputQualString += '\n';

        qFile->write(outputQualString);
    }
    catch(exception& e) {
        m->errorOut(e, "QualityScores", "printQScores");
        exit(1);
    }
}
/**************************************************************************************************/

void QualityScores::printQScores(ofstream& qFile){
	try {

	    double expected_errors = calculateExpectedErrors();

			qFile << '>' << seqName << '\t' << expected_errors << endl;

			for(int i=0;i<seqLength;i++){ qFile << qScores[i] << ' '; } qFile << endl;
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "printQScores");
		exit(1);
	}
}
/**************************************************************************************************/

void QualityScores::printQScores(ostream& qFile){
    try {

        double expected_errors = calculateExpectedErrors();

        qFile << '>' << seqName << '\t' << expected_errors << endl;

        for(int i=0;i<seqLength;i++){
            qFile << qScores[i] << ' ';
        }
        qFile << endl;

    }
    catch(exception& e) {
        m->errorOut(e, "QualityScores", "printQScores");
        exit(1);
    }
}


/**************************************************************************************************/

void QualityScores::trimQScores(int start, int end){
	try {
		vector<int> hold;
        
		if(end == -1){
			hold = vector<int>(qScores.begin()+start, qScores.end());
			qScores = hold;
		}
		if(start == -1){
			if(qScores.size() > end){
				hold = vector<int>(qScores.begin(), qScores.begin()+end);
				qScores = hold;
			}
		}

		seqLength = qScores.size();
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "trimQScores");
		exit(1);
	}
}

/**************************************************************************************************/

void QualityScores::flipQScores(){
	try {

		vector<int> temp = qScores;
		for(int i=0;i<seqLength;i++){
			qScores[seqLength - i - 1] = temp[i];
		}

	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "flipQScores");
		exit(1);
	}
}

/**************************************************************************************************/

bool QualityScores::stripQualThreshold(Sequence& sequence, double qThreshold){
	try {
		string rawSequence = sequence.getUnaligned();
		int seqLength = sequence.getNumBases();

		if(seqName != sequence.getName()){
			m->mothurOut("sequence name mismatch btwn fasta: " + sequence.getName() + " and qual file: " + seqName);
			m->mothurOutEndLine();	m->setControl_pressed(true);
		}

		int end;
		for(int i=0;i<seqLength;i++){
			end = i;
			if(qScores[i] < qThreshold){
				break;
			}
		}

		//every score passed
		if (end == (seqLength-1)) { end = seqLength; }

		sequence.setUnaligned(rawSequence.substr(0,end));
		trimQScores(-1, end);

		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "flipQScores");
		exit(1);
	}

}

/**************************************************************************************************/

bool QualityScores::stripQualRollingAverage(Sequence& sequence, double qThreshold, bool logTransform){
	try {
		string rawSequence = sequence.getUnaligned();
		int seqLength = sequence.getNumBases();

		if(seqName != sequence.getName()){
			m->mothurOut("sequence name mismatch btwn fasta: " + sequence.getName() + " and qual file: " + seqName);
			m->mothurOutEndLine();
		}

		int end = -1;
		double rollingSum = 0.0000;
        double value = 0.0;

		for(int i=0;i<seqLength;i++){

            if (logTransform)   {
                rollingSum += (double)pow(10.0, qScores[i]);
                value = log10(rollingSum / (double)(i+1));

            } //Sum 10^Q
            else                {
                rollingSum += (double)qScores[i];
                value = rollingSum / (double)(i+1);
            }


			if(value < qThreshold){
				end = i;
				break;
			}
		}

		if(end == -1){	end = seqLength;	}


		sequence.setUnaligned(rawSequence.substr(0,end));
		trimQScores(-1, end);


		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "flipQScores");
		exit(1);
	}

}

/**************************************************************************************************/

bool QualityScores::stripQualWindowAverage(Sequence& sequence, int stepSize, int windowSize, double qThreshold, bool logTransform){
	try {
		string rawSequence = sequence.getUnaligned();
		int seqLength = sequence.getNumBases();

		if(seqName != sequence.getName()){
			m->mothurOut("sequence name mismatch between fasta: " + sequence.getName() + " and qual file: " + seqName);
			m->mothurOutEndLine();
		}

		int end = windowSize;
		int start = 0;

		if(seqLength < windowSize) {	return 0;	}

		while((start+windowSize) < seqLength){
			double windowSum = 0.0000;

			for(int i=start;i<end;i++){
                if (logTransform)   {  windowSum += pow(10.0, qScores[i]);  }
                else                {  windowSum += qScores[i];             }
			}
			double windowAverage = 0.0;
            if (logTransform)   { windowAverage = log10(windowSum / (double)(end-start)); }
            else                { windowAverage = windowSum / (double)(end-start);      }

			if(windowAverage < qThreshold){
				end = end - stepSize;
				break;
			}

			start += stepSize;
			end = start + windowSize;

			if(end >= seqLength){	end = seqLength;	}

		}

		if(end == -1){	end = seqLength;	}

		//failed first window
		if (end < windowSize) { return 0; }

		sequence.setUnaligned(rawSequence.substr(0,end));
		trimQScores(-1, end);

		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "stripQualWindowAverage");
		exit(1);
	}

}

/**************************************************************************************************/

double QualityScores::calculateExpectedErrors(void){

	double expected_errors = 0.0000;

	for(int i=0;i<seqLength;i++){
        expected_errors += pow(10.0, -qScores[i]/10.0);
	}

	return expected_errors;
}

/**************************************************************************************************/

double QualityScores::calculateAverage(bool logTransform){

	double aveQScore = 0.0000;

	for(int i=0;i<seqLength;i++){
        if (logTransform)   {  aveQScore += pow(10.0, qScores[i]);  }
        else                {  aveQScore += qScores[i];             }
	}

    if (logTransform)   {  aveQScore = log10(aveQScore /(double) seqLength);    }
    else                {  aveQScore /= (double) seqLength;                     }

	return aveQScore;
}

/**************************************************************************************************/

bool QualityScores::cullQualAverage(Sequence& sequence, double qAverage, bool logTransform){
	try {
		string rawSequence = sequence.getUnaligned();
		bool success = 0;	//guilty until proven innocent

		if(seqName != sequence.getName())	{
			m->mothurOut("sequence name mismatch btwn fasta: " + sequence.getName() + " and qual file: " + seqName);
			m->mothurOutEndLine();
		}

		double aveQScore = calculateAverage(logTransform);

        if (m->getDebug()) { m->mothurOut("[DEBUG]: " + sequence.getName() + " average = " + toString(aveQScore) + "\n"); }

		if(aveQScore >= qAverage)	{	success = 1;	}
		else						{	success = 0;	}

		return success;
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "cullQualAverage");
		exit(1);
	}
}

/**************************************************************************************************/

void QualityScores::updateQScoreErrorMap(map<char, vector<int> >& qualErrorMap, string errorSeq, int start, int stop, int weight){
	try {

		int seqLength = errorSeq.size();

		int qIndex = start - 1;

		for(int i=0;i<seqLength;i++){

			if(errorSeq[i] == 'm')		{	qualErrorMap['m'][qScores[qIndex]] += weight;	}
			else if(errorSeq[i] == 's')	{	qualErrorMap['s'][qScores[qIndex]] += weight;	}
			else if(errorSeq[i] == 'i')	{	qualErrorMap['i'][qScores[qIndex]] += weight;	}
			else if(errorSeq[i] == 'a')	{	qualErrorMap['a'][qScores[qIndex]] += weight;	}
			else if(errorSeq[i] == 'd')	{	/*	there are no qScores for deletions	*/		}

			if(errorSeq[i] != 'd')		{	qIndex++;	}

			if(qIndex > stop){	break;	}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "updateQScoreErrorMap");
		exit(1);
	}
}

/**************************************************************************************************/

void QualityScores::updateForwardMap(vector<vector<int> >& forwardMap, int start, int stop, int weight){
	try {

		int index = 0;
		for(int i=start-1;i<stop;i++){
			forwardMap[index++][qScores[i]] += weight;
		}

	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "updateForwardMap");
		exit(1);
	}
}

/**************************************************************************************************/

void QualityScores::updateReverseMap(vector<vector<int> >& reverseMap, int start, int stop, int weight){
	try {

		int index = 0;
		for(int i=stop-1;i>=start-1;i--){
			reverseMap[index++][qScores[i]] += weight;
		}

	}
	catch(exception& e) {
		m->errorOut(e, "QualityScores", "updateReverseMap");
		exit(1);
	}
}

/**************************************************************************************************/
