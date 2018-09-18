/*
 *  preclustercommand.cpp
 *  Mothur
 *
 *  Created by westcott on 12/21/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "preclustercommand.h"
#include "deconvolutecommand.h"
#include "summary.hpp"

//**********************************************************************************************************************
vector<string> PreClusterCommand::setParameters(){
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta-name",false,true,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","name",false,false,true); parameters.push_back(pname);
    CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","count",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
		CommandParameter pdiffs("diffs", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pdiffs);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
    CommandParameter palign("align", "Multiple", "needleman-gotoh-blast-noalign", "needleman", "", "", "","",false,false); parameters.push_back(palign);
    CommandParameter pmatch("match", "Number", "", "1.0", "", "", "","",false,false); parameters.push_back(pmatch);
    CommandParameter pmismatch("mismatch", "Number", "", "-1.0", "", "", "","",false,false); parameters.push_back(pmismatch);
    CommandParameter pgapopen("gapopen", "Number", "", "-2.0", "", "", "","",false,false); parameters.push_back(pgapopen);
    CommandParameter pgapextend("gapextend", "Number", "", "-1.0", "", "", "","",false,false); parameters.push_back(pgapextend);
    CommandParameter palpha("alpha", "Number", "", "2.0", "", "", "","",false,false); parameters.push_back(palpha);
    CommandParameter pdelta("delta", "Number", "", "2.0", "", "", "","",false,false); parameters.push_back(pdelta);
    CommandParameter pmethod("method", "Multiple", "simple-unoise-tree-deblur", "simple", "", "", "","",false,false); parameters.push_back(pmethod);

    CommandParameter perror_rate("error_rate", "Number", "", "0.005", "", "", "","",false,false); parameters.push_back(perror_rate);
    CommandParameter pindel_prob("indel_prob", "Number", "", "0.01", "", "", "","",false,false); parameters.push_back(pindel_prob);
    CommandParameter pmax_indels("max_indels", "Number", "", "3", "", "", "","",false,false); parameters.push_back(pmax_indels);
    CommandParameter perror_dist("error_dist", "String", "", "1-0.06-0.02-0.02-0.01-0.005-0.005-0.005-0.001-0.001-0.001-0.0005", "", "", "","",false,false); parameters.push_back(perror_dist);

		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
    CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);

		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string PreClusterCommand::getHelpString(){
	try {
		string helpString = "";
		helpString += "The pre.cluster command groups sequences that are within a given number of base mismatches.\n";
		helpString += "The pre.cluster command outputs a new fasta and name file.\n";
		helpString += "The pre.cluster command parameters are fasta, name, group, count, method, processors and diffs. The fasta parameter is required. \n";
		helpString += "The name parameter allows you to give a list of seqs that are identical. This file is 2 columns, first column is name or representative sequence, second column is a list of its identical sequences separated by commas.\n";
		helpString += "The group parameter allows you to provide a group file so you can cluster by group. \n";
    helpString += "The count parameter allows you to provide a count file so you can cluster by group. \n";
		helpString += "The diffs parameter allows you to specify maximum number of mismatched bases allowed between sequences in a grouping. The default is 1.\n";
    helpString += "The method parameter allows you to specify the algorithm to use to complete the preclusterign step. Possible methods include simple, tree, unoise, and deblur.  Default=simple.\n";
    helpString += "The align parameter allows you to specify the alignment align_method to use.  Your options are: gotoh, needleman, blast and noalign. The default is needleman.\n";
    helpString += "The match parameter allows you to specify the bonus for having the same base. The default is 1.0.\n";
    helpString += "The mistmatch parameter allows you to specify the penalty for having different bases.  The default is -1.0.\n";
    helpString += "The gapopen parameter allows you to specify the penalty for opening a gap in an alignment. The default is -2.0.\n";
    helpString += "The gapextend parameter allows you to specify the penalty for extending a gap in an alignment.  The default is -1.0.\n";
    helpString += "The alpha parameter allows you to specify the alpha value for the beta formula, which is used in the unoise algorithm. The default is 2.0.\n";
    helpString += "The delta parameter allows you to specify the delta value, which describes the amount of amplification between rounds of PCR. It is used in the tree algorithm. The default is 2.0.\n";
		helpString += "The error_rate parameter is used with the deblur algorithm and is the expected mean error rate, as a fraction, of the data going into this command.\n";
		helpString += "The indel_prob parameter is used with the deblur algorithm and is the expected fraction of sequences that have an insertion or deletion, of the data going into this command.\n";
		helpString += "The max_indels parameter is used with the deblur algorithm and is the maximum number of insertions or deletions you expect to be in the data going into this command.\n";
		helpString += "The error_dist parameter is used with the deblur algorithm and is the fraction of sequences you expect to have 0, 1, 2, 3, etc. errors. Should start with 1 and be separated by hyphens (e.g. 1-0.06-0.02-0.02-0.01-0.005-0.005-0.005-0.001-0.001-0.001-0.0005). Alternatively, you can use error_dist=binomial and the command will determine the distribution for you\n";

		helpString += "The pre.cluster command should be in the following format: \n";
		helpString += "pre.cluster(fasta=yourFastaFile, names=yourNamesFile, diffs=yourMaxDiffs) \n";
		helpString += "Example pre.cluster(fasta=amazon.fasta, diffs=2).\n";
		;
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string PreClusterCommand::getOutputPattern(string type) {
    try {
        string pattern = "";

        if (type == "fasta") {  pattern = "[filename],precluster,[extension]"; }
        else if (type == "name") {  pattern = "[filename],precluster.names"; }
        else if (type == "count") {  pattern = "[filename],precluster.count_table"; }
        else if (type == "map") {  pattern =  "[filename],precluster.map"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }

        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "PreClusterCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
PreClusterCommand::PreClusterCommand(){
	try {
		abort = true; calledHelp = true;
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
    outputTypes["count"] = tempOutNames;
		outputTypes["map"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "PreClusterCommand");
		exit(1);
	}
}
//**************************************************************************************************

PreClusterCommand::PreClusterCommand(string option) {
	try {
		abort = false; calledHelp = false;

		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}

		else {
			vector<string> myArray = setParameters();

			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();

			ValidParameters validParameter;
			map<string, string>::iterator it;

			//check to make sure all parameters are valid for command
			for (map<string, string>::iterator it2 = parameters.begin(); it2 != parameters.end(); it2++) {
				if (validParameter.isValidParameter(it2->first, myArray, it2->second) != true) {  abort = true;  }
			}

			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["name"] = tempOutNames;
			outputTypes["map"] = tempOutNames;
            outputTypes["count"] = tempOutNames;

			//if the user changes the input directory command factory will send this info to us in the output parameter
			string inputDir = validParameter.valid(parameters, "inputdir");
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}

				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}

				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}

                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
			}

			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not found") {
				fastafile = current->getFastaFile();
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (fastafile == "not open") { abort = true; }
			else { current->setFastaFile(fastafile); }

			//if the user changes the output directory command factory will send this info to us in the output parameter
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){
				outputDir = "";
				outputDir += util.hasPath(fastafile); //if user entered a file with a path then preserve it
			}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not found") { namefile =  "";  }
			else if (namefile == "not open") { namefile = ""; abort = true; }
			else {  current->setNameFile(namefile); }

			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not found") { groupfile =  "";  bygroup = false; }
			else if (groupfile == "not open") { abort = true; groupfile =  ""; }
			else {   current->setGroupFile(groupfile); bygroup = true;  }

      countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not found") { countfile =  "";   }
			else if (countfile == "not open") { abort = true; countfile =  ""; }
			else {
        current->setCountFile(countfile);
        ct.readTable(countfile, true, false);
        if (ct.hasGroupInfo()) { bygroup = true; }
        else { bygroup = false;  }
      }

      if ((namefile != "") && (countfile != "")) {
        m->mothurOut("[ERROR]: you may only use one of the following: name or count.");
 				m->mothurOutEndLine();
				abort = true;
      }

      if ((groupfile != "") && (countfile != "")) {
        m->mothurOut("[ERROR]: you may only use one of the following: group or count.");
 				m->mothurOutEndLine();
				abort=true;
      }


			string temp	= validParameter.valid(parameters, "diffs"); if(temp == "not found"){	temp = "1"; }
			util.mothurConvert(temp, diffs);

			temp = validParameter.valid(parameters, "processors"); if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);

      temp = validParameter.valid(parameters, "method"); if(temp == "not found"){  temp = "simple"; }
			pc_method = temp;

      temp = validParameter.valid(parameters, "match"); if(temp == "not found"){	temp = "1.0";			}
      util.mothurConvert(temp, match);

      temp = validParameter.valid(parameters, "mismatch"); if(temp == "not found"){	temp = "-1.0";			}
      util.mothurConvert(temp, misMatch);
      if (misMatch > 0) { m->mothurOut("[ERROR]: mismatch must be negative.\n"); abort=true; }

      temp = validParameter.valid(parameters, "gapopen"); if(temp == "not found"){	temp = "-2.0";			}
      util.mothurConvert(temp, gapOpen);
      if (gapOpen > 0) { m->mothurOut("[ERROR]: gapopen must be negative.\n"); abort=true; }

      temp = validParameter.valid(parameters, "gapextend");	if (temp == "not found"){	temp = "-1.0";			}
      util.mothurConvert(temp, gapExtend);
      if (gapExtend > 0) { m->mothurOut("[ERROR]: gapextend must be negative.\n"); abort=true; }

      temp = validParameter.valid(parameters, "alpha");	if (temp == "not found"){	temp = "2.0";			}
      util.mothurConvert(temp, alpha);
      if (alpha < 0) { m->mothurOut("[ERROR]: alpha must be positive.\n"); abort=true; }

      temp = validParameter.valid(parameters, "delta");	if (temp == "not found"){	temp = "2.0";			}
      util.mothurConvert(temp, delta);
      if (delta < 0) { m->mothurOut("[ERROR]: delta must be positive.\n"); abort=true; }

      temp = validParameter.valid(parameters, "error_dist");
			if (temp == "not found"){
				error_dist = {1, 0.06, 0.02, 0.02, 0.01, 0.005, 0.005, 0.005, 0.001, 0.001, 0.001, 0.0005};
			} else if(temp == "binomial"){
				error_dist = {-100};
			} else {
				string probability;
				istringstream probabilityStream(temp);
				while (getline(probabilityStream, probability, '-')){
					error_dist.push_back(stof(probability));
				}
			}

				// vector<double>

      temp = validParameter.valid(parameters, "error_rate");	if (temp == "not found"){	temp = "0.005";			}
      util.mothurConvert(temp, error_rate);
      if (error_rate < 0) { m->mothurOut("[ERROR]: error_rate must be positive.\n"); abort=true; }
      else if (error_rate > 1) { m->mothurOut("[ERROR]: error_rate is a fraction and should be less than 1.\n"); abort=true; }

      temp = validParameter.valid(parameters, "indel_prob");	if (temp == "not found"){	temp = "0.01";			}
      util.mothurConvert(temp, indel_prob);
      if (indel_prob < 0) { m->mothurOut("[ERROR]: indel_prob must be positive.\n"); abort=true; }
      else if (indel_prob > 1) { m->mothurOut("[ERROR]: indel_prob is a fraction and should be less than 1.\n"); abort=true; }

      temp = validParameter.valid(parameters, "max_indels");	if (temp == "not found"){	temp = "3";			}
      util.mothurConvert(temp, max_indels);
      if (indel_prob < 0) { m->mothurOut("[ERROR]: max_indels must be positive.\n"); abort=true; }



      align = validParameter.valid(parameters, "align");		if (align == "not found"){	align = "needleman";	}

      align_method = "unaligned";

      if (countfile == "") {
          if (namefile == "") {
              vector<string> files; files.push_back(fastafile);
              if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
          }
      }


			if(pc_method == "tree"){

				diffs = 1;

			} else {

				Summary summary_data(processors);
				if(countfile != ""){
					summary_data.summarizeFasta(fastafile, countfile, "");
				} else if(namefile != ""){
					summary_data.summarizeFasta(fastafile, namefile, "");
				}
				int median_length = summary_data.getLength()[3];
				int max_abund = summary_data.getMaxAbundance();


				if(pc_method == "unoise"){

					diffs = int((log2(max_abund)-1) / alpha + 1);

				} else if(pc_method == "deblur"){

					if(error_dist[0] == -100){ //construct the binomial distribution

						error_dist.clear();

						for(int i=0;i<100;i++){

							double choose = 1;
							double k = 1;

							for(int j=1;j<=i;j++){
								choose *= (median_length - j + 1);
								k *= j;
							}
							choose = choose/k;

							double a = pow(error_rate, i);
							double b = pow(1 - error_rate, median_length - i);

							if(a != 0 && b != 0){
								error_dist.push_back(choose * a * b);
								if(error_dist[i] < 0.1 / max_abund){ break;	}
							} else {
								break;
							}
						}
						error_dist[0] = 1;
					}

					double mod_factor = pow((1-error_rate), median_length);

					for(int i=0;i<error_dist.size();i++){
						error_dist[i] /= (double)mod_factor;
					}
					diffs = error_dist.size();
				}
			}
		}

	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "PreClusterCommand");
		exit(1);
	}
}
/************************************************************/
struct seqPNode {
    int numIdentical;
    Sequence seq;
    string filteredSeq;
    string names;
    bool active;
    int diffs;
    seqPNode() {}
    seqPNode(int n, Sequence s, string nm) : numIdentical(n), seq(s), names(nm), active(1) { diffs = 0; filteredSeq = "";}
    ~seqPNode() {}
};

/************************************************************/

inline bool comparePriorityTopDown(seqPNode* first, seqPNode* second) {
  if (first->numIdentical > second->numIdentical){
 		return true;
	}
  return false;
}

//**************************************************************************************************

struct preClusterData {
  string fastafile, namefile, groupfile, countfile, pc_method, align_method, align, newMName;
  OutputWriter* newFName;
  OutputWriter* newNName;
  MothurOut* m;
  int start, end, count, diffs, length;
  vector<string> groups;
  bool hasCount, hasName;
  float match, misMatch, gapOpen, gapExtend, alpha, delta, error_rate, indel_prob, max_indels;
  Utils util;
	vector<float> error_dist;
  vector<string> outputNames;
  map<string, vector<string> > outputTypes;
  vector<seqPNode*> alignSeqs; //maps the number of identical seqs to a sequence
  Alignment* alignment;

				// double error_rate = 0.005;error_rate
				// double indel_prob = 0.01;indel_prob
				// double max_indels = 3;max_indels

  ~preClusterData() { if (alignment != NULL) { delete alignment; } }

  preClusterData(){}

  preClusterData(string f, string n, string g, string c, string pcm, string am, OutputWriter* nff,  OutputWriter* nnf, string nmf, vector<string> gr) {
    fastafile = f;
    groupfile = g;
		pc_method = pcm;
		align_method = am;
    newFName = nff;
    newNName = nnf;
    newMName = nmf;
    groups = gr;
    hasName = false;
    namefile = n; if (namefile != "") { hasName = true; }
    hasCount = false;
    countfile = c; if (countfile != "") { hasCount = true; }
    count=0;
    m = MothurOut::getInstance();
  }

  void setVariables(int st, int en, int d, string pcm, string am, string al, float ma, float misma, float gpOp, float gpEx, float a, float del, float me, float ip, float mi, vector<float> ed) {
    start = st;
    end = en;
    diffs = d;
		pc_method = pcm;
    align_method = am;
    align = al;
    match = ma;
    misMatch = misma;
    gapExtend = gpEx;
    gapOpen = gpOp;
		alpha = a;
		delta = del;
		error_rate = me;
		indel_prob = ip;
		max_indels = mi;
		error_dist = ed;
    length = 0;

    if (align_method == "unaligned") {
      if(align == "gotoh")	{	alignment = new GotohOverlap(gapOpen, gapExtend, match, misMatch, 1000);	}
      else if(align == "needleman")	{	alignment = new NeedlemanOverlap(gapOpen, match, misMatch, 1000);			}
      else if(align == "blast")		{	alignment = new BlastAlignment(gapOpen, gapExtend, match, misMatch);		}
      else if(align == "noalign")		{	alignment = new NoAlign();													}
      else {
          m->mothurOut(align + " is not a valid alignment option. I will run the command using needleman.");
          m->mothurOutEndLine();
          alignment = new NeedlemanOverlap(gapOpen, match, misMatch, 1000);
      }
    } else { alignment = NULL; }

  }

};

/**************************************************************************************************/

int calcMisMatches(string seq1, string seq2, preClusterData* params){
  try {
    int numBad = 0;

    if (params->align_method == "unaligned") {
      //align to eachother
      Sequence seqI("seq1", seq1);
      Sequence seqJ("seq2", seq2);

      //align seq2 to seq1 - less abundant to more abundant
      params->alignment->align(seqJ.getUnaligned(), seqI.getUnaligned());
      seq2 = params->alignment->getSeqAAln();
      seq1 = params->alignment->getSeqBAln();

      //chop gap ends
      int startPos = 0;
      int endPos = seq2.length()-1;
      for (int i = 0; i < seq2.length(); i++)		{  if (isalpha(seq2[i])) { startPos = i; break; } }
      for (int i = seq2.length()-1; i >= 0; i--){  if (isalpha(seq2[i])) { endPos = i; break;	} 	}

      //count number of diffs
      for (int i = startPos; i <= endPos; i++) {
        if (seq2[i] != seq1[i]) { numBad++; }
        if (numBad > params->diffs) { return params->length;  } //too far to cluster
      }

    } else {
      //count diffs
      for (int i = 0; i < seq1.length(); i++) {
        //do they match
        if (seq1[i] != seq2[i]) { numBad++; }
        if (numBad > params->diffs) { return params->length;  } //too far to cluster
      }
    }
    return numBad;
  }
  catch(exception& e) {
      params->m->errorOut(e, "PreClusterCommand", "calcMisMatches");
      exit(1);
  }
}

/**************************************************************************************************/

vector<int> calcMisMatchesIndels(string seq1, string seq2, preClusterData* params){
  try {


    int numSubstitutions = 0;
		int numInDels = 0;

    if (params->align_method == "unaligned") {
      //align to eachother
      Sequence seqI("seq1", seq1);
      Sequence seqJ("seq2", seq2);

      //align seq2 to seq1 - less abundant to more abundant
      params->alignment->align(seqJ.getUnaligned(), seqI.getUnaligned());
      seq2 = params->alignment->getSeqAAln();
      seq1 = params->alignment->getSeqBAln();

      //chop gap ends
      int startPos = 0;
      int endPos = seq2.length()-1;
      for (int i = 0; i < seq2.length(); i++)		{  if (isalpha(seq2[i])) { startPos = i; break; } }
      for (int i = seq2.length()-1; i >= 0; i--){  if (isalpha(seq2[i])) { endPos = i; break;	} 	}

      //count number of diffs
      for (int i = startPos; i <= endPos; i++) {

        if (seq2[i] != seq1[i] && (isalpha(seq1[i]) && isalpha(seq2[i]))) { numSubstitutions++; }
				else if(seq2[i] != seq1[i] && (isalpha(seq1[i]) || isalpha(seq2[i]))){ numInDels++;	}

        if (numSubstitutions > params->diffs) { numSubstitutions = params->length; break; }
      }

    } else {
      //count diffs
      for (int i = 0; i < seq1.length(); i++) {
        //do they match
        if (seq2[i] != seq1[i] && (isalpha(seq1[i]) && isalpha(seq2[i]))) { numSubstitutions++; }
				else if(seq2[i] != seq1[i] && (isalpha(seq1[i]) || isalpha(seq2[i]))){ numInDels++;	}

        if (numSubstitutions > params->diffs) { numSubstitutions = params->length; break; }
      }
    }
    return vector<int>{numSubstitutions, numInDels};
  }
  catch(exception& e) {
      params->m->errorOut(e, "PreClusterCommand", "calcMisMatchesIndels");
      exit(1);
  }
}

/**************************************************************************************************/

int process(string group, string newMapFile, preClusterData* params){
  try {
    ofstream out;

		if(params->pc_method != "deblur"){
	    params->util.openOutputFile(newMapFile, out);
			out << "ideal_seq\terror_seq\tabundance\tdiffs\tsequence" << endl;
		}

    int count = 0;
    long long numSeqs = params->alignSeqs.size();

    vector<int> originalCount(numSeqs);

    for (int i = 0; i < numSeqs; i++) {
			originalCount[i] = params->alignSeqs[i]->numIdentical;
		}


		if(params->pc_method == "simple"){
	    for (int i = 0; i < numSeqs; i++) {

	      if (params->alignSeqs[i]->active) {  //this sequence has not been merged yet

	        string chunk = params->alignSeqs[i]->seq.getName() + "\t" + params->alignSeqs[i]->seq.getName() + "\t" + toString(originalCount[i]) + "\t" + toString(0) + "\t" + params->alignSeqs[i]->seq.getAligned() + "\n";

	        //try to merge it with all smaller seqs
	        for (int j = i+1; j < numSeqs; j++) {

	          if (params->m->getControl_pressed()) { out.close(); return 0; }

	          if (params->alignSeqs[j]->active && (originalCount[j] < originalCount[i])) {  //this sequence has not been merged yet
	            //are you within "diff" bases
	            int mismatch = params->length;
	            if (params->align_method == "unaligned") {
								mismatch = calcMisMatches(params->alignSeqs[i]->seq.getAligned(),
	 																				params->alignSeqs[j]->seq.getAligned(), params);
							} else {
								mismatch = calcMisMatches(params->alignSeqs[i]->filteredSeq,
																					params->alignSeqs[j]->filteredSeq, params);
							}

	            if (mismatch <= params->diffs) {
	              //merge
	              params->alignSeqs[i]->names += ',' + params->alignSeqs[j]->names;
	              params->alignSeqs[i]->numIdentical += params->alignSeqs[j]->numIdentical;

	              chunk += params->alignSeqs[i]->seq.getName() + "\t" + params->alignSeqs[j]->seq.getName() + "\t" + toString(originalCount[j]) + "\t" + toString(mismatch) + "\t" + params->alignSeqs[j]->seq.getAligned() + "\n";

	              params->alignSeqs[j]->active = 0;
	              params->alignSeqs[j]->numIdentical = 0;
	              params->alignSeqs[j]->diffs = mismatch;
	              count++;
	            }
	          }//end if j active
	        }//end for loop j

		      //remove from active list
		      params->alignSeqs[i]->active = 0;

		      out << chunk;

		     }//end if active i
	      if(i % 100 == 0)	{
					params->m->mothurOutJustToScreen(group + toString(i) + "\t" + toString(numSeqs - count) + "\t" + toString(count)+"\n");
				}
	    }

			if(numSeqs % 100 != 0)	{ params->m->mothurOut(group + toString(numSeqs) + "\t" + toString(numSeqs - count) + "\t" + toString(count) + "\n"); 	}

 		} else if(params->pc_method == "unoise") {

			vector<double> beta(params->diffs, 0);
			for(int i=0;i<beta.size();i++){
				beta[i] = pow(0.5, params->alpha * i + 1.0);
			}


	    for (int i = 0; i < numSeqs; i++) {

	      if (params->alignSeqs[i]->active) {  //this sequence has not been merged yet

	        string chunk = params->alignSeqs[i]->seq.getName() + "\t" + params->alignSeqs[i]->seq.getName() + "\t" + toString(originalCount[i]) + "\t" + toString(0) + "\t" + params->alignSeqs[i]->seq.getAligned() + "\n";

	        //try to merge it with all smaller seqs
	        for (int j = i+1; j < numSeqs; j++) {

	          if (params->m->getControl_pressed()) { out.close(); return 0; }

	          if (params->alignSeqs[j]->active && originalCount[j] < originalCount[i]) {  //this sequence has not been merged yet

	            int mismatch = params->length;
							double skew = (double)originalCount[j]/(double)originalCount[i];

	            if (params->align_method == "unaligned") {
								mismatch = calcMisMatches(params->alignSeqs[i]->seq.getAligned(),
	 																				params->alignSeqs[j]->seq.getAligned(), params);
							} else {
								mismatch = calcMisMatches(params->alignSeqs[i]->filteredSeq,
																					params->alignSeqs[j]->filteredSeq, params);
							}

	            if (mismatch <= params->diffs && skew <= beta[mismatch]) { //merge
	              params->alignSeqs[i]->names += ',' + params->alignSeqs[j]->names;
	              params->alignSeqs[i]->numIdentical += params->alignSeqs[j]->numIdentical;

	              chunk += params->alignSeqs[i]->seq.getName() + "\t" + params->alignSeqs[j]->seq.getName() + "\t" + toString(originalCount[j]) + "\t" + toString(mismatch) + "\t" + params->alignSeqs[j]->seq.getAligned() + "\n";

	              params->alignSeqs[j]->active = 0;
	              params->alignSeqs[j]->numIdentical = 0;
	              params->alignSeqs[j]->diffs = mismatch;
	              count++;
	            }
	          }//end if j active
	        }//end for loop j

		      //remove from active list
		      params->alignSeqs[i]->active = 0;

		      out << chunk;

		    }//end if active i
	      if(i % 100 == 0)	{ params->m->mothurOutJustToScreen(group + toString(i) + "\t" + toString(numSeqs - count) + "\t" + toString(count)+"\n"); 	}

	    }

	  	if(numSeqs % 100 != 0)	{ params->m->mothurOut(group + toString(numSeqs) + "\t" + 	toString(numSeqs - count) + "\t" + toString(count) + "\n"); 	}

		} else if(params->pc_method == "tree") {
			params->diffs = 1;

			params->m->mothurOutJustToScreen("Determining which sequences can be merged...");
			cout.flush();

			vector<vector<bool> > mergable(numSeqs);
			vector<vector<int> > mismatches(numSeqs);

			mergable[0].resize(numSeqs, FALSE);
			mismatches[0].resize(numSeqs, 1000);

	    for (int i = 1; i < numSeqs; i++) {
				mergable[i].resize(numSeqs, FALSE);
				mismatches[i].resize(numSeqs, 1000);

        for (int j = 0; j < i; j++) {

          if (params->m->getControl_pressed()) { out.close(); return 0; }

          int mismatch = params->length;

					if(originalCount[j] > originalCount[i] * params->delta){
	          if (params->align_method == "unaligned") {
							mismatches[i][j] = calcMisMatches(params->alignSeqs[i]->seq.getAligned(),
																					params->alignSeqs[j]->seq.getAligned(), params);
						} else {
							mismatches[i][j] = calcMisMatches(params->alignSeqs[i]->filteredSeq,
																				params->alignSeqs[j]->filteredSeq, params);
						}

						mergable[i][j] = (mismatches[i][j] == 1);
					}
        }
      }
			params->m->mothurOutJustToScreen(" done\n");

			vector<int> cluster(numSeqs, -1);

			params->m->mothurOutJustToScreen("Clusterng sequences...");
			cout.flush();


			for(int i=0;i<numSeqs;i++){

        if (params->m->getControl_pressed()) { out.close(); return 0; }

				vector<int> indices_to_merge;

				for(int j=0;j<numSeqs;j++){
					if(mergable[j][i] == TRUE){	indices_to_merge.push_back(j); }
				}

				while(indices_to_merge.size() != 0){
					for(int j=0;j<numSeqs;j++){
						bool to_merge = mergable[j][i];

						if(!to_merge){
							for(int k=0;k<indices_to_merge.size();k++){

								if(mergable[j][indices_to_merge[k]]){
									to_merge = TRUE;
									break;
								}

							}
						}
						mergable[i][j] = to_merge;
						mergable[j][i] = mergable[i][j];

						if (params->m->getControl_pressed()) { out.close(); return 0; }
					}

					cluster[i] = i;

					for(int k=0;k<indices_to_merge.size();k++){
						cluster[indices_to_merge[k]] = i;

						for(int l=0;l<numSeqs;l++){
							mergable[indices_to_merge[k]][l] = FALSE;
							mergable[l][indices_to_merge[k]] = FALSE;
						}

		        if (params->m->getControl_pressed()) { out.close(); return 0; }
					}
					mergable[i][i] = FALSE;

					indices_to_merge.clear();

					for(int j=0;j<numSeqs;j++){
						if(mergable[j][i] == TRUE){	indices_to_merge.push_back(j); }
					}

					if (params->m->getControl_pressed()) { out.close(); return 0; }
				}
			}

			params->m->mothurOutJustToScreen("done.\n");

			vector<string> chunk(numSeqs, "");

			for(int i=0;i<numSeqs;i++){

				if(i == cluster[i] || cluster[i] == -1){

	        chunk[i] = params->alignSeqs[i]->seq.getName() + "\t" + params->alignSeqs[i]->seq.getName() + "\t" + toString(originalCount[i]) + "\t" + toString(0) + "\t" + params->alignSeqs[i]->seq.getAligned() + "\n";

				} else {

	        params->alignSeqs[cluster[i]]->names += ',' + params->alignSeqs[i]->names;
	        params->alignSeqs[cluster[i]]->numIdentical += params->alignSeqs[i]->numIdentical;

        	chunk[cluster[i]] += params->alignSeqs[cluster[i]]->seq.getName() + "\t" + params->alignSeqs[i]->seq.getName() + "\t" + toString(originalCount[i]) + "\t" + toString(mismatches[i][cluster[i]]) + "\t" + params->alignSeqs[i]->seq.getAligned() + "\n";

        	params->alignSeqs[i]->active = 0;
        	params->alignSeqs[i]->numIdentical = 0;

        	count++;

				}
			}

			for(int i=0;i<numSeqs;i++){
				if(chunk[i] != ""){
					out << chunk[i];
				}
			}

	  } else if(params->pc_method == "deblur") {

			vector<double> weights(numSeqs, 0);
			for(int i=0;i<numSeqs;i++){
				weights[i] = (double)params->alignSeqs[i]->numIdentical;
			}

			for(int i=0;i<numSeqs;i++){

        if (params->m->getControl_pressed()) { out.close(); return 0; }

				if(i % 100 == 0){ cout << i << endl; }

				if(weights[i] <= 0){	continue;	}

				int max_h_dist = params->error_dist.size();
				vector<double> expected_bad_reads(max_h_dist, 0);
 				for(int j=0;j<max_h_dist;j++){
 					expected_bad_reads[j] = params->error_dist[j] * weights[i];
				}

				if(expected_bad_reads[1] < 0.1){	continue;	}

				for(int j=0;j<numSeqs;j++){

        	if (params->m->getControl_pressed()) { out.close(); return 0; }

					if(i == j) { continue; }
					if(weights[j] <= 0){	continue;	}

					vector<int> nSubsInDels(2, 0);
          if (params->align_method == "unaligned") {
						nSubsInDels = calcMisMatchesIndels(params->alignSeqs[i]->seq.getAligned(),
																				params->alignSeqs[j]->seq.getAligned(), params);
					} else {
						nSubsInDels = calcMisMatchesIndels(params->alignSeqs[i]->filteredSeq,
																				params->alignSeqs[j]->filteredSeq, params);
					}

					if(nSubsInDels[0] >= max_h_dist) { continue; }

					double correction = expected_bad_reads[nSubsInDels[0]];

					if(nSubsInDels[1] > params->max_indels){
						correction = 0;
					} else if(nSubsInDels[1] > 0){
						correction = correction * params->indel_prob;
					}

					weights[j] -= correction;
				}
			}

			for(int i=0;i<numSeqs;i++){
				params->alignSeqs[i]->numIdentical = round(weights[i]);

				if(weights[i] <= 0){
					params->alignSeqs[i]->numIdentical = 0;
					params->alignSeqs[i]->active = 0;
					count++;
				}
			}

		} else {

			cout << "fail!\n";

		}

		if(params->pc_method != "deblur"){
	  	out.close();
		}

		return count;
  }
  catch(exception& e) {
      params->m->errorOut(e, "PreClusterCommand", "process");
      exit(1);
  }
}

/**************************************************************************************************/

void filterSeqs(vector<seqPNode*>& alignSeqs, int length, MothurOut* m){
    try {
        string filterString = "";
        Filters F;

        F.setLength(length);
        F.initialize();
        F.setFilter(string(length, '1'));

        for (int i = 0; i < alignSeqs.size(); i++) { F.getFreqs(alignSeqs[i]->seq); }

        F.setNumSeqs(alignSeqs.size());
        F.doVerticalAllBases();
        filterString = F.getFilter();

        //run filter
        for (int i = 0; i < alignSeqs.size(); i++) {
            if (m->getControl_pressed()) { break; }
            alignSeqs[i]->filteredSeq = "";
            string align = alignSeqs[i]->seq.getAligned();
            for(int j=0;j<length;j++){ if(filterString[j] == '1'){ alignSeqs[i]->filteredSeq += align[j]; } }
        }
    }
    catch(exception& e) {
        m->errorOut(e, "PreClusterCommand", "filterSeqs");
        exit(1);
    }
}

/**************************************************************************************************/

vector<seqPNode*> readFASTA(CountTable ct, preClusterData* params, long long& num){
    try {
        map<string, string> nameMap;
        map<string, string>::iterator it;
        if (params->hasName) { params->util.readNames(params->namefile, nameMap); }

        ifstream inFasta;
        params->util.openInputFile(params->fastafile, inFasta);
        set<int> lengths;
        vector<seqPNode*> alignSeqs;

        while (!inFasta.eof()) {

            if (params->m->getControl_pressed()) { inFasta.close(); break; }

            Sequence seq(inFasta);  params->util.gobble(inFasta);

            if (seq.getName() != "") {  //can get "" if commented line is at end of fasta file

                if (params->hasName) {
                    it = nameMap.find(seq.getName());

                    if (it == nameMap.end()) { params->m->mothurOut("[ERROR]: " + seq.getName() + " is not in your names file, please correct.\n"); exit(1); }
                    else{
                        string second = it->second;
                        int numReps = params->util.getNumNames(second);
                        seqPNode* tempNode = new seqPNode(numReps, seq, second);
                        alignSeqs.push_back(tempNode);
                        lengths.insert(seq.getAligned().length());
                    }
                }else { //no names file, you are identical to yourself
                    int numRep = 1;
                    if (params->hasCount) { numRep = ct.getNumSeqs(seq.getName()); }
                    seqPNode* tempNode = new seqPNode(numRep, seq, seq.getName());
                    alignSeqs.push_back(tempNode);
                    lengths.insert(seq.getAligned().length());
                }
            }
        }
        inFasta.close();

        params->length = *(lengths.begin());

        if (lengths.size() > 1) { params->align_method = "unaligned"; }
        else if (lengths.size() == 1) {  params->align_method = "aligned"; filterSeqs(alignSeqs, params->length, params->m); }

        //sort seqs by number of identical seqs
        sort(alignSeqs.begin(), alignSeqs.end(), comparePriorityTopDown);
        num = alignSeqs.size();

        return alignSeqs;
    }
    catch(exception& e) {
        params->m->errorOut(e, "PreClusterCommand", "readFASTA");
        exit(1);
    }
}

/**************************************************************************************************/

void print(string newfasta, string newname, preClusterData* params){
    try {
        ofstream outFasta;
        ofstream outNames;

        params->util.openOutputFile(newfasta, outFasta);
        params->util.openOutputFile(newname, outNames);

        if (params->countfile != "")  { outNames << "Representative_Sequence\ttotal\n";  }

        if (params->countfile != "") {
            for (int i = 0; i < params->alignSeqs.size(); i++) {
                if (params->alignSeqs[i]->numIdentical != 0) {
                    params->alignSeqs[i]->seq.printSequence(outFasta);
                    outNames << params->alignSeqs[i]->seq.getName() << '\t' << params->alignSeqs[i]->numIdentical << endl;
                }
            }
        }else {
            for (int i = 0; i < params->alignSeqs.size(); i++) {
                if (params->alignSeqs[i]->numIdentical != 0) {
                    params->alignSeqs[i]->seq.printSequence(outFasta);
                    outNames << params->alignSeqs[i]->seq.getName() << '\t' << params->alignSeqs[i]->names << endl;
                }
            }
        }
        outFasta.close();
        outNames.close();

    }
    catch(exception& e) {
        params->m->errorOut(e, "PreClusterCommand", "print");
        exit(1);
    }
}

/**************************************************************************************************/
int PreClusterCommand::execute(){
	try {

		if (abort) { if (calledHelp) { return 0; }  return 2;	}

		long start = time(NULL);

		string fileroot = outputDir + util.getRootName(util.getSimpleName(fastafile));
		map<string, string> variables;
		variables["[filename]"] = fileroot;

		string newNamesFile = getOutputFileName("name",variables);
		string newCountFile = getOutputFileName("count",variables);
		string newMapFile = getOutputFileName("map",variables); //add group name if by group

		variables["[extension]"] = util.getExtension(fastafile);
		string newFastaFile = getOutputFileName("fasta", variables);
		outputNames.push_back(newFastaFile); outputTypes["fasta"].push_back(newFastaFile);

		if (countfile == "") {
			outputNames.push_back(newNamesFile); outputTypes["name"].push_back(newNamesFile);
		}	else {
			outputNames.push_back(newCountFile); outputTypes["count"].push_back(newCountFile);
		}

		if (bygroup) {
			//clear out old files
			ofstream outFasta; util.openOutputFile(newFastaFile, outFasta); outFasta.close();
			ofstream outNames; util.openOutputFile(newNamesFile, outNames);  outNames.close();
			newMapFile = fileroot + "precluster.";

			createProcessesGroups(newFastaFile, newNamesFile, newMapFile);

			if (countfile != "") {
				mergeGroupCounts(newCountFile, newNamesFile, newFastaFile);
      } else {
        //run unique.seqs for deconvolute results
        string inputString = "fasta=" + newFastaFile;
        if (namefile != "") { inputString += ", name=" + newNamesFile; }

        m->mothurOut("\n/******************************************/\n");
        m->mothurOut("Running command: unique.seqs(" + inputString + ")\n");
        current->setMothurCalling(true);

        Command* uniqueCommand = new DeconvoluteCommand(inputString);
        uniqueCommand->execute();

        map<string, vector<string> > filenames = uniqueCommand->getOutputFiles();

        delete uniqueCommand;
        current->setMothurCalling(false);
        m->mothurOut("/******************************************/"); m->mothurOutEndLine();

        util.renameFile(filenames["fasta"][0], newFastaFile);
        util.renameFile(filenames["name"][0], newNamesFile);
			}

      if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	}	 return 0; }

			m->mothurOut("It took " + toString(time(NULL) - start) + " secs to run pre.cluster.\n");

		}else {
            if (processors != 1) { m->mothurOut("When using running without group information mothur can only use 1 processor, continuing."); m->mothurOutEndLine(); processors = 1; }

            preClusterData* params = new preClusterData(fastafile, namefile, groupfile, countfile, pc_method, align_method, NULL, NULL, newMapFile, nullVector);
            params->setVariables(0,0, diffs, pc_method, align_method, align, match, misMatch, gapOpen, gapExtend, alpha, delta, error_rate, indel_prob, max_indels, error_dist);

            //reads fasta file and return number of seqs
            long long numSeqs = 0; params->alignSeqs = readFASTA(ct, params, numSeqs); //fills alignSeqs and makes all seqs active
            length = params->length;

			if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	}  return 0; }

			if (numSeqs == 0) { m->mothurOut("Error reading fasta file...please correct.\n");  return 0;  }
			if (diffs > length) { m->mothurOut("Error: diffs is set to " + toString(diffs) + " which is greater than your sequence length of " + toString(length) + ".\n");   return 0;  }

			int count = process("", newMapFile, params);
			if(params->pc_method != "deblur"){
				outputNames.push_back(newMapFile); outputTypes["map"].push_back(newMapFile);
			}

			if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	}  return 0; }

			m->mothurOut("Total number of sequences before precluster was " + toString(params->alignSeqs.size()) + ".\n");
			m->mothurOut("pre.cluster removed " + toString(count) + " sequences.\n\n");
			if (countfile != "") { newNamesFile = newCountFile; }
            print(newFastaFile, newNamesFile, params);
            for (int i = 0; i < params->alignSeqs.size(); i++) {  delete params->alignSeqs[i]; } params->alignSeqs.clear();
			m->mothurOut("It took " + toString(time(NULL) - start) + " secs to cluster " + toString(numSeqs) + " sequences.\n");
		}

		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	}  return 0; }

		m->mothurOut("\nOutput File Names: \n");
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

		//set fasta file as new current fastafile
		string currentName = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
		}

		itTypes = outputTypes.find("name");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setNameFile(currentName); }
		}

        itTypes = outputTypes.find("count");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
		}

		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "execute");
		exit(1);
	}
}

/**************************************************************************************************/

vector<seqPNode*> loadSeqs(map<string, string>& thisName, vector<Sequence>& thisSeqs, map<string, int>& thisCount, string group, long long& num, bool hasName, bool hasCount, int& length, string& align_method){
    MothurOut* m = MothurOut::getInstance();
    try {
        set<int> lengths;
        bool error = false; num = 0;
        Utils util;
        vector<seqPNode*> alignSeqs;

        for (int i = 0; i < thisSeqs.size(); i++) {

            if (m->getControl_pressed()) { return alignSeqs; }

            if (hasName) {
                map<string, string>::iterator it = thisName.find(thisSeqs[i].getName());

                //should never be true since parser checks for this
                if (it == thisName.end()) { m->mothurOut("[ERROR]: " + thisSeqs[i].getName() + " is not in your names file, please correct.\n"); error = true; }
                else{
                    //get number of reps
                    int numReps = util.getNumNames(it->second);
                    seqPNode* tempNode = new seqPNode(numReps, thisSeqs[i], it->second);
                    alignSeqs.push_back(tempNode);
                    lengths.insert(thisSeqs[i].getAligned().length());
                }
            }else { //no names file, you are identical to yourself
                int numRep = 1;
                if (hasCount) {
                    map<string, int>::iterator it2 = thisCount.find(thisSeqs[i].getName());

                    //should never be true since parser checks for this
                    if (it2 == thisCount.end()) { m->mothurOut("[ERROR]: " + thisSeqs[i].getName() + " is not in your count file, please correct.\n"); error = true; }
                    else { numRep = it2->second;  }
                }
                seqPNode* tempNode = new seqPNode(numRep, thisSeqs[i], thisSeqs[i].getName());
                alignSeqs.push_back(tempNode);
                lengths.insert(thisSeqs[i].getAligned().length());
            }
        }

        length = *(lengths.begin());

        if (lengths.size() > 1) { align_method = "unaligned"; }
        else if (lengths.size() == 1) {  align_method = "aligned"; filterSeqs(alignSeqs, length, m); }

        //sanity check
        if (error) { m->setControl_pressed(true); }

        thisSeqs.clear();

        //sort seqs by number of identical seqs
        sort(alignSeqs.begin(), alignSeqs.end(), comparePriorityTopDown);

        num = alignSeqs.size();

        return alignSeqs;
    }
    catch(exception& e) {
      m->errorOut(e, "PreClusterCommand", "loadSeqs");
      exit(1);
    }
}

/**************************************************************************************************/

void printData(string group, preClusterData* params){
    try {
        if ((params->hasCount) && (group == ""))  {
					params->newNName->write("Representative_Sequence\ttotal\n");
				}

        if (params->hasCount) {
            if (group != "") {

                for (int i = 0; i < params->alignSeqs.size(); i++) {
                    if (params->alignSeqs[i]->numIdentical != 0) {
                        params->alignSeqs[i]->seq.printSequence(params->newFName);

												if(params->pc_method == "deblur"){
                        	params->newNName->write(group + '\t' + params->alignSeqs[i]->seq.getName() + '\t' + toString(params->alignSeqs[i]->numIdentical) + '\n');
												} else {
                        	params->newNName->write(group + '\t' + params->alignSeqs[i]->seq.getName() + '\t' + params->alignSeqs[i]->names + '\n');
												}
                    }
                }
            }
            else {

                for (int i = 0; i < params->alignSeqs.size(); i++) {
                    if (params->alignSeqs[i]->numIdentical != 0) {
                        params->alignSeqs[i]->seq.printSequence(params->newFName);
                        params->newNName->write(params->alignSeqs[i]->seq.getName()  + '\t' + toString(params->alignSeqs[i]->numIdentical) + '\n');
                    }
                }
            }
        }else {
            for (int i = 0; i < params->alignSeqs.size(); i++) {
                if (params->alignSeqs[i]->numIdentical != 0) {
                    params->alignSeqs[i]->seq.printSequence(params->newFName);
                    params->newNName->write(params->alignSeqs[i]->seq.getName() + '\t' + params->alignSeqs[i]->names + '\n');
                }
            }
        }

    }
    catch(exception& e) {
        params->m->errorOut(e, "PreClusterCommand", "printData");
        exit(1);
    }
}
/**************************************************************************************************/

long long driverGroups(preClusterData* params){
	try {
    vector<string> subsetGroups;
    for (int i = params->start; i < params->end; i++) {  subsetGroups.push_back(params->groups[i]);  }

    //parse fasta and name file by group
    SequenceCountParser* cparser = NULL;
    SequenceParser* parser = NULL;
    if (params->hasCount) {
        cparser = new SequenceCountParser(params->countfile, params->fastafile, subsetGroups);
    }else {
        if (params->hasName) { parser = new SequenceParser(params->groupfile, params->fastafile, params->namefile, subsetGroups);	}
        else				{ parser = new SequenceParser(params->groupfile, params->fastafile, subsetGroups);                      }
    }

		long long numSeqs = 0;

		//precluster each group
		for (int i = params->start; i < params->end; i++) {
			if (params->m->getControl_pressed()) { if (params->hasCount) { delete cparser; }else { delete parser; } return numSeqs; }

      params->m->mothurOut("\nProcessing group " + params->groups[i] + ":\n");

			time_t start = time(NULL);
			map<string, string> thisNameMap;
      vector<Sequence> thisSeqs;

			if (params->groupfile != "")        {  thisSeqs = parser->getSeqs(params->groups[i]);       }
      else if (params->hasCount)          { thisSeqs = cparser->getSeqs(params->groups[i]);       }

      if (params->hasName)                {  thisNameMap = parser->getNameMap(params->groups[i]); }

      map<string, int> thisCount;
      if (params->hasCount) { thisCount = cparser->getCountTable(params->groups[i]);  }

      long long num = 0;
			params->alignSeqs = loadSeqs(thisNameMap, thisSeqs, thisCount, params->groups[i], num, params->hasName, params->hasCount, params->length, params->align_method);
      numSeqs += num;

			if (params->m->getControl_pressed()) {   return 0; }

      if (params->align_method == "aligned") {
				if (params->diffs > params->length) {
					params->m->mothurOut("[ERROR]: diffs is greater than your sequence length.\n");
					params->m->setControl_pressed(true);
					return 0;
				}
			}

      string extension = params->groups[i]+".map";
			long long count = process(params->groups[i]+"\t", params->newMName+extension, params);

			if(params->pc_method != "deblur"){
				params->outputNames.push_back(params->newMName+extension);
				params->outputTypes["map"].push_back(params->newMName+extension);
			}

			if (params->m->getControl_pressed()) {  return 0; }

			params->m->mothurOut("Total number of sequences before pre.cluster was " + toString(params->alignSeqs.size()) + ".\n");
			params->m->mothurOut("pre.cluster removed " + toString(count) + " sequences.\n\n");
			printData(params->groups[i], params);
      for (int i = 0; i < params->alignSeqs.size(); i++) {  delete params->alignSeqs[i]; } params->alignSeqs.clear();

			params->m->mothurOut("It took " + toString(time(NULL) - start) + " secs to cluster " + toString(num) + " sequences.\n");
		}

    if (params->hasCount) { delete cparser; }else { delete parser; }

		return numSeqs;
	}
	catch(exception& e) {
		params->m->errorOut(e, "PreClusterCommand", "driverGroups");
		exit(1);
	}
}

/**************************************************************************************************/

int PreClusterCommand::mergeGroupCounts(string newcount, string newname, string newfasta){
	try {
		ifstream inNames;
    util.openInputFile(newname, inNames);

		if(pc_method != "deblur"){
	    string group, first, second;
	    set<string> uniqueNames;

	    while (!inNames.eof()) {

	      if (m->getControl_pressed()) { break; }

	      inNames >> group; util.gobble(inNames);
	      inNames >> first; util.gobble(inNames);
	      inNames >> second; util.gobble(inNames);

	      vector<string> names;
	      util.splitAtComma(second, names);

	      uniqueNames.insert(first);

	      int total = ct.getGroupCount(first, group);
	      for (int i = 1; i < names.size(); i++) {
	          total += ct.getGroupCount(names[i], group);
	          ct.setAbund(names[i], group, 0);
	      }
	      ct.setAbund(first, group, total);
	    }
	    inNames.close();

		} else { //for deblur

	    string group, unique_sequence;
			int count;

			ct.clearTable();

	    while (!inNames.eof()) {

	      if (m->getControl_pressed()) { break; }

	      inNames >> group; util.gobble(inNames);
	      inNames >> unique_sequence; util.gobble(inNames);
	      inNames >> count; util.gobble(inNames);

				ct.setAbund(unique_sequence, group, count);

			}
		}

    vector<string> namesOfSeqs = ct.getNamesOfSeqs();
    for (int i = 0; i < namesOfSeqs.size(); i++) {
        if (ct.getNumSeqs(namesOfSeqs[i]) == 0) {
            ct.remove(namesOfSeqs[i]);
        }
    }

    ct.printTable(newcount);
    // util.mothurRemove(newname);

    if (bygroup) { //if by group, must remove the duplicate seqs that are named the same
      ifstream in;
      util.openInputFile(newfasta, in);

      ofstream out;
      util.openOutputFile(newfasta+"temp", out);

      int count = 0;
      set<string> already;
      while(!in.eof()) {
        if (m->getControl_pressed()) { break; }

        Sequence seq(in); util.gobble(in);

        if (seq.getName() != "") {
          count++;
          if (already.count(seq.getName()) == 0) {
            seq.printSequence(out);
            already.insert(seq.getName());
          }
        }
      }
      in.close();
      out.close();
      util.mothurRemove(newfasta);
      util.renameFile(newfasta+"temp", newfasta);
    }

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "mergeGroupCounts");
		exit(1);
	}
}

/**************************************************************************************************/

void PreClusterCommand::createProcessesGroups(string newFName, string newNName, string newMFile) {
  try {
    //parse fasta and name file by group
    vector<string> groups;
    if (countfile != "") { CountTable ct; ct.testGroups(countfile, groups); }
    else { GroupMap gp; gp.readMap(groupfile); groups = gp.getNamesOfGroups(); }

    //sanity check
    if (groups.size() < processors) { processors = groups.size(); m->mothurOut("Reducing processors to " + toString(groups.size()) + ".\n"); }

    //divide the groups between the processors
    vector<linePair> lines;
    int remainingPairs = groups.size();
    int startIndex = 0;
    for (int remainingProcessors = processors; remainingProcessors > 0; remainingProcessors--) {
      int numPairs = remainingPairs; //case for last processor
      if (remainingProcessors != 1) { numPairs = ceil(remainingPairs / remainingProcessors); }
      lines.push_back(linePair(startIndex, (startIndex+numPairs))); //startIndex, endIndex
      startIndex = startIndex + numPairs;
      remainingPairs = remainingPairs - numPairs;
    }

    //create array of worker threads
    vector<thread*> workerThreads;
    vector<preClusterData*> data;

    auto synchronizedFastaFile = std::make_shared<SynchronizedOutputFile>(newFName);
    auto synchronizedNameFile = std::make_shared<SynchronizedOutputFile>(newNName);

    //Lauch worker threads
    for (int i = 0; i < processors-1; i++) {
      OutputWriter* threadFastaWriter = new OutputWriter(synchronizedFastaFile);
      OutputWriter* threadNameWriter = new OutputWriter(synchronizedNameFile);

      preClusterData* dataBundle = new preClusterData(fastafile, namefile, groupfile, countfile, pc_method, align_method, threadFastaWriter, threadNameWriter, newMFile, groups);
      dataBundle->setVariables(lines[i+1].start, lines[i+1].end, diffs, pc_method, align_method, align, match, misMatch, gapOpen, gapExtend, alpha, delta, error_rate, indel_prob, max_indels, error_dist);
      data.push_back(dataBundle);

      workerThreads.push_back(new thread(driverGroups, dataBundle));
    }

    OutputWriter* threadFastaWriter = new OutputWriter(synchronizedFastaFile);
    OutputWriter* threadNameWriter = new OutputWriter(synchronizedNameFile);

    preClusterData* dataBundle = new preClusterData(fastafile, namefile, groupfile, countfile, pc_method, align_method, threadFastaWriter, threadNameWriter, newMFile, groups);
    dataBundle->setVariables(lines[0].start, lines[0].end, diffs, pc_method, align_method, align, match, misMatch, gapOpen, gapExtend, alpha, delta, error_rate, indel_prob, max_indels, error_dist);

    driverGroups(dataBundle);

    outputNames.insert(outputNames.end(), dataBundle->outputNames.begin(), dataBundle->outputNames.end());
    outputTypes.insert(dataBundle->outputTypes.begin(), dataBundle->outputTypes.end());

    for (int i = 0; i < processors-1; i++) {
      workerThreads[i]->join();

      delete data[i]->newFName;
      delete data[i]->newNName;

      outputNames.insert(outputNames.end(), data[i]->outputNames.begin(), data[i]->outputNames.end());
      outputTypes.insert(data[i]->outputTypes.begin(), data[i]->outputTypes.end());

      delete data[i];
      delete workerThreads[i];
    }
    delete threadFastaWriter;
    delete threadNameWriter;
    delete dataBundle;

  }
  catch(exception& e) {
      m->errorOut(e, "PreClusterCommand", "createProcessesGroups");
      exit(1);
  }
}

/**************************************************************************************************/
