/*
 *  preclustercommand.cpp
 *  Mothur
 *
 *  Created by westcott on 12/21/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "splitgroupscommand.h"
#include "preclustercommand.h"
#include "uniqueseqscommand.h"
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
        CommandParameter palign("align", "Multiple", "needleman-gotoh-noalign", "needleman", "", "", "","",false,false); parameters.push_back(palign);
        CommandParameter pmatch("match", "Number", "", "1.0", "", "", "","",false,false); parameters.push_back(pmatch);
        CommandParameter pmismatch("mismatch", "Number", "", "-1.0", "", "", "","",false,false); parameters.push_back(pmismatch);
        CommandParameter pgapopen("gapopen", "Number", "", "-2.0", "", "", "","",false,false); parameters.push_back(pgapopen);
        CommandParameter pgapextend("gapextend", "Number", "", "-1.0", "", "", "","",false,false); parameters.push_back(pgapextend);
        CommandParameter palpha("alpha", "Number", "", "2.0", "", "", "","",false,false); parameters.push_back(palpha);
        CommandParameter pdelta("delta", "Number", "", "2.0", "", "", "","",false,false); parameters.push_back(pdelta);
        CommandParameter pmethod("method", "Multiple", "simple-unoise-tree-deblur", "simple", "", "", "","",false,false); parameters.push_back(pmethod);
        CommandParameter pclump("clump", "Multiple", "lessthan-lessthanequal", "lessthan", "", "", "","",false,false); parameters.push_back(pclump);
        CommandParameter perror_rate("error_rate", "Number", "", "0.005", "", "", "","",false,false); parameters.push_back(perror_rate);
        CommandParameter pindel_prob("indel_prob", "Number", "", "0.01", "", "", "","",false,false); parameters.push_back(pindel_prob);
        CommandParameter pmax_indels("max_indels", "Number", "", "3", "", "", "","",false,false); parameters.push_back(pmax_indels);
        CommandParameter perror_dist("error_dist", "String", "", "1-0.06-0.02-0.02-0.01-0.005-0.005-0.005-0.001-0.001-0.001-0.0005", "", "", "","",false,false); parameters.push_back(perror_dist);
        
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
        CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["fasta"] = tempOutNames;
        outputTypes["map"] = tempOutNames;
        outputTypes["count"] = tempOutNames;

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
		helpString += "The pre.cluster command outputs a new fasta and count file.\n";
		helpString += "The pre.cluster command parameters are fasta, count, method, processors and diffs. The fasta parameter is required. \n";
		helpString += "The name parameter allows you to give a list of seqs that are identical. This file is 2 columns, first column is name or representative sequence, second column is a list of its identical sequences separated by commas.\n";
		helpString += "The group parameter allows you to provide a group file so you can cluster by group. \n";
    helpString += "The count parameter allows you to provide a count file so you can cluster by group. \n";
		helpString += "The diffs parameter allows you to specify maximum number of mismatched bases allowed between sequences in a grouping. The default is 1.\n";
    helpString += "The method parameter allows you to specify the algorithm to use to complete the preclusterign step. Possible methods include simple, tree, unoise, and deblur.  Default=simple.\n";
    helpString += "The clump parameter allows you to specify which reads can be combined. Possible options include lessthan and lessthanequal. lessthan -> merge reads with less abundance. lessthanequal -> merge reads with less than or equal abundance Default=lessthan.\n";
    helpString += "The align parameter allows you to specify the alignment align_method to use.  Your options are: gotoh, needleman and noalign. The default is needleman.\n";
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
        //else if (type == "name") {  pattern = "[filename],precluster.names"; }
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
//**************************************************************************************************
PreClusterCommand::PreClusterCommand(string option) : Command() {
	try {

		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }

		else {
			OptionParser parser(option, setParameters());
			map<string, string> parameters = parser.getParameters();

			ValidParameters validParameter;
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not found") {
				fastafile = current->getFastaFile();
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter.\n"); }
				else { 	m->mothurOut("[ERROR]: You have no current fastafile and the fasta parameter is required.\n");  abort = true; }
			}
			else if (fastafile == "not open") { abort = true; }
			else { current->setFastaFile(fastafile); }

			 if (outputdir == ""){ outputdir += util.hasPath(fastafile); }

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
                CountTable ct;
                if (ct.testGroups(countfile)) { bygroup = true; } //check for groups without reading
                else { bygroup = false;  }
            }
            
            if ((namefile != "") && (countfile != "")) { m->mothurOut("[ERROR]: you may only use one of the following: name or count.\n"); abort = true; }
            
            if ((groupfile != "") && (countfile != "")) { m->mothurOut("[ERROR]: you may only use one of the following: group or count.\n"); abort=true; }
            
            
            string temp	= validParameter.valid(parameters, "diffs"); if(temp == "not found"){	temp = "1"; }
            util.mothurConvert(temp, diffs);
            
            temp = validParameter.valid(parameters, "processors"); if (temp == "not found"){	temp = current->getProcessors();	}
            processors = current->setProcessors(temp);
            
            temp = validParameter.valid(parameters, "method"); if(temp == "not found"){  temp = "simple"; }
            pc_method = temp;
            
            if ((pc_method == "simple") || (pc_method == "tree") || (pc_method == "unoise") || (pc_method == "deblur")) { }
            else { m->mothurOut("[ERROR]: Not a valid precluster method.  Valid preclustering algorithms include simple, tree, unoise, and deblur. Using simple.\n");  pc_method= "simple"; }
            
            temp = validParameter.valid(parameters, "clump"); if(temp == "not found"){  temp = "lessthan"; }
            clump = temp;
            
            if ((clump == "lessthan") || (clump == "lessthanequal")) { }
            else { m->mothurOut("[ERROR]: Not a valid clump method.  Valid clumping options are lessthan and lessthanequal. Using lessthan.\n");  clump = "lessthan"; }
            
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
            
            if (!abort) {
                if ((namefile != "") || (groupfile != "")) { //convert to count
                    
                    string rootFileName = namefile;
                    if (rootFileName == "") { rootFileName = groupfile; }
                    
                    if (outputdir == "") { outputdir = util.hasPath(rootFileName); }
                    string outputFileName = outputdir + util.getRootName(util.getSimpleName(rootFileName)) + "count_table";
                    
                    CountTable ct; ct.createTable(namefile, groupfile, nullVector); ct.printCompressedTable(outputFileName);
                    outputNames.push_back(outputFileName);
                    
                    current->setCountFile(outputFileName);
                    countfile = outputFileName;
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

inline bool comparePriorityAbundance(seqPNode* first, seqPNode* second) {
  if (first->numIdentical > second->numIdentical){
 		return true;
	}
  return false;
}
//**************************************************************************************************

struct preClusterData {
  string fastafile, countfile, pc_method, align_method, align, newMName, clump;
  OutputWriter* newNName;
  MothurOut* m;
  int start, end, count, diffs, length, numGroups;
  vector<string> groups;
  bool hasCount;
  float match, misMatch, gapOpen, gapExtend, alpha, delta, error_rate, indel_prob, max_indels;
  Utils util;
	vector<float> error_dist;
  vector<string> outputNames;
  map<string, vector<string> > outputTypes;
  vector<seqPNode*> alignSeqs; //maps the number of identical seqs to a sequence. filled and freed by functions
  Alignment* alignment;
    map<string, vector<string> > parsedFiles;

				// double error_rate = 0.005;error_rate
				// double indel_prob = 0.01;indel_prob
				// double max_indels = 3;max_indels

  ~preClusterData() { if (alignment != NULL) { delete alignment; } }

  preClusterData(){}

  preClusterData(map<string, vector<string> > g2f, string f, string c, string pcm, string am,  string cl, OutputWriter* nnf, string nmf, vector<string> gr) {
    fastafile = f;
    pc_method = pcm;
    align_method = am;
    clump = cl;
    newNName = nnf;
    newMName = nmf;
    groups = gr;
    hasCount = false;
    countfile = c; if (countfile != "") { hasCount = true; }
    count=0;
    m = MothurOut::getInstance();
      parsedFiles = g2f;
  }

  void setVariables(int d, string pcm, string am, string al, float ma, float misma, float gpOp, float gpEx, float a, float del, float me, float ip, float mi, vector<float> ed) {
    
      numGroups = groups.size();
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

void mergeSeqs(seqPNode* representative, seqPNode* duplicate, string& chunk, int mismatch, int originalCount, preClusterData* params){
    try {
        //merge
        if (params->hasCount) { representative->clusteredIndexes.clear(); duplicate->clusteredIndexes.clear(); } //we use numIdentical to build the count table, don't need the names.
        else {
            representative->clusteredIndexes.insert(representative->clusteredIndexes.end(), duplicate->clusteredIndexes.begin(), duplicate->clusteredIndexes.end());
            duplicate->clusteredIndexes.clear();
        }
        representative->numIdentical += duplicate->numIdentical;
        
        chunk += representative->name + "\t" + duplicate->name + "\t" + toString(originalCount) + "\t" + toString(mismatch) + "\t" + duplicate->sequence + "\n";
        duplicate->numIdentical = 0;
        duplicate->diffs = mismatch;
    }
    catch(exception& e) {
        params->m->errorOut(e, "PreClusterCommand", "mergeSeqs");
        exit(1);
    }
}
/**************************************************************************************************/

int process(string group, string newMapFile, preClusterData* params){
    try {
        ofstream out;
        
        if(params->pc_method != "deblur"){
            params->util.openOutputFile(newMapFile, out);
            if (params->align_method == "unaligned") {
                out << "ideal_seq\terror_seq\tabundance\tdiffs\tsequence" << endl;
            }else {
                out << "ideal_seq\terror_seq\tabundance\tdiffs\tfiltered_sequence" << endl;
            }
        }
        
        int count = 0;
        long long numSeqs = params->alignSeqs.size();
        
        vector<int> originalCount(numSeqs);
        
        for (int i = 0; i < numSeqs; i++) { originalCount[i] = params->alignSeqs[i]->numIdentical; }
        
        bool lessThan = true;
        if (params->clump == "lessthanequal") { lessThan = false; }
        
        if(params->pc_method == "simple"){
            for (int i = 0; i < numSeqs; i++) {
                
                if (params->alignSeqs[i]->numIdentical != 0) {  //this sequence has not been merged yet
                    
                    string chunk = params->alignSeqs[i]->name + "\t" + params->alignSeqs[i]->name + "\t" + toString(originalCount[i]) + "\t" + toString(0) + "\t" + params->alignSeqs[i]->sequence + "\n";
                    
                    //try to merge it with all smaller seqs
                    for (int j = i+1; j < numSeqs; j++) {
                        
                        if (params->m->getControl_pressed()) { out.close(); return 0; }
                        
                        bool ableToMerge = false;
                        if (lessThan) { //default
                            if (originalCount[j] < originalCount[i]) { ableToMerge = true; }
                        }else { //less than equal to
                            if (originalCount[j] <= originalCount[i]) { ableToMerge = true; }
                        }
                        if ((params->alignSeqs[j]->numIdentical != 0) && (ableToMerge)) {  //this sequence has not been merged yet //
                            //are you within "diff" bases
                            int mismatch = calcMisMatches(params->alignSeqs[i]->sequence, params->alignSeqs[j]->sequence, params);
                            
                            if (mismatch <= params->diffs) { mergeSeqs(params->alignSeqs[i], params->alignSeqs[j], chunk, mismatch, originalCount[j], params); count++; }
                        }
                    }
                    out << chunk;
                }
                if(i % 100 == 0)	{ params->m->mothurOutJustToScreen(group + toString(i) + "\t" + toString(numSeqs - count) + "\t" + toString(count)+"\n"); }
            }
            
            if(numSeqs % 100 != 0)	{ params->m->mothurOut(group + toString(numSeqs) + "\t" + toString(numSeqs - count) + "\t" + toString(count) + "\n"); 	}
            
        } else if(params->pc_method == "unoise") {
            
            vector<double> beta(params->diffs+1, 0);
            for(int i=0;i<beta.size();i++){
                beta[i] = pow(0.5, params->alpha * i + 1.0);
            }
            
            for (int i = 0; i < numSeqs; i++) {
                
                if (params->alignSeqs[i]->numIdentical != 0) {  //this sequence has not been merged yet
                    
                    string chunk = params->alignSeqs[i]->name + "\t" + params->alignSeqs[i]->name + "\t" + toString(originalCount[i]) + "\t" + toString(0) + "\t" + params->alignSeqs[i]->sequence + "\n";
                    
                    //try to merge it with all smaller seqs
                    for (int j = i+1; j < numSeqs; j++) {
                        
                        if (params->m->getControl_pressed()) { out.close(); return 0; }
                        
                        if (params->alignSeqs[j]->numIdentical != 0 && (originalCount[j] < originalCount[i])) {  //this sequence has not been merged yet
                            
                            double skew = (double)originalCount[j]/(double)originalCount[i];
                            
                            int mismatch = calcMisMatches(params->alignSeqs[i]->sequence, params->alignSeqs[j]->sequence, params);
                        
                            if (mismatch <= params->diffs) {
                                if (skew <= beta[mismatch]) {
                                    mergeSeqs(params->alignSeqs[i], params->alignSeqs[j], chunk, mismatch, originalCount[j], params); count++;
                                }
                            }
                        }
                    }
                    out << chunk;
                }
                if(i % 100 == 0)	{ params->m->mothurOutJustToScreen(group + toString(i) + "\t" + toString(numSeqs - count) + "\t" + toString(count)+"\n"); 	}
            }
            
            if(numSeqs % 100 != 0)	{ params->m->mothurOut(group + toString(numSeqs) + "\t" + 	toString(numSeqs - count) + "\t" + toString(count) + "\n"); 	}
            
        } else if(params->pc_method == "tree") {
            
            vector<int> cluster(numSeqs, -1);
            for(int i=0;i<cluster.size();i++){
                cluster[i] = i;
            }
            
            for (int i=0;i<numSeqs-1;i++) {
                
                if (params->m->getControl_pressed()) { out.close(); return 0; }
                
                for (int j=i+1;j<numSeqs;j++) {
                    
                    if (params->m->getControl_pressed()) { out.close(); return 0; }
                    
                    if(originalCount[i] > originalCount[j] * params->delta){
                        int mismatches = calcMisMatches(params->alignSeqs[i]->sequence, params->alignSeqs[j]->sequence, params);
                        
                        if(mismatches == 1){ cluster[j] = cluster[i]; }
                    }
                }
                
                if(cluster[i] == i) count++;
                
                if(i % 100 == 0)	{ params->m->mothurOutJustToScreen(group + toString(i) + "\t" + toString(numSeqs - count) + "\t" + toString(count)+"\n"); 	}
            }
            params->m->mothurOutJustToScreen(group + toString(numSeqs) + "\t" + toString(numSeqs - count) + "\t" + toString(count)+"\n");
            
            
            count = 0;
            vector<string> chunk(numSeqs, "");
            
            for(int i=0;i<numSeqs;i++){
                
                if(i == cluster[i] || cluster[i] == -1){
                    
                    chunk[i] = params->alignSeqs[i]->name + "\t" + params->alignSeqs[i]->name + "\t" + toString(originalCount[i]) + "\t" + toString(0) + "\t" + params->alignSeqs[i]->sequence + "\n";
                    
                } else {
                    
                    params->alignSeqs[cluster[i]]->clusteredIndexes.insert(params->alignSeqs[i]->clusteredIndexes.begin(), params->alignSeqs[i]->clusteredIndexes.end(), params->alignSeqs[cluster[i]]->clusteredIndexes.end());  params->alignSeqs[i]->clusteredIndexes.clear();
                    params->alignSeqs[cluster[i]]->numIdentical += params->alignSeqs[i]->numIdentical;
                    
                    params->diffs = params->length;
                    int mismatches = calcMisMatches(params->alignSeqs[i]->sequence, params->alignSeqs[cluster[i]]->sequence, params);
                    
                    chunk[cluster[i]] += params->alignSeqs[cluster[i]]->name + "\t" + params->alignSeqs[i]->name + "\t" + toString(originalCount[i]) + "\t" + toString(mismatches) + "\t" + params->alignSeqs[i]->sequence + "\n";
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
                    
                    nSubsInDels = calcMisMatchesIndels(params->alignSeqs[i]->sequence, params->alignSeqs[j]->sequence, params);
                    
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
                    count++;
                }
            }
            
        } else { cout << "fail!\n"; }
        
        if(params->pc_method != "deblur"){ out.close(); }
        
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

        for (int i = 0; i < alignSeqs.size(); i++) { F.getFreqs(alignSeqs[i]->sequence); }

        F.setNumSeqs(alignSeqs.size());
        F.doVerticalAllBases();
        filterString = F.getFilter();

        //run filter
        for (int i = 0; i < alignSeqs.size(); i++) {
            if (m->getControl_pressed()) { break; }
            string filteredSeq = "";
            string align = alignSeqs[i]->sequence;
            for(int j=0;j<length;j++){ if(filterString[j] == '1'){ filteredSeq += align[j]; } }
            alignSeqs[i]->sequence = filteredSeq;
        }
    }
    catch(exception& e) {
        m->errorOut(e, "PreClusterCommand", "filterSeqs");
        exit(1);
    }
}

/**************************************************************************************************/
//seqPNode(string na, string seq, int n, string nm) : numIdentical(n), name(na), sequence(seq), clusteredNames(nm) { diffs = 0; active = true; }
vector<seqPNode*> readFASTA(preClusterData* params, long long& num){
    try {
        CountTable ct;
        if (params->hasCount) {  ct.readTable(params->countfile, false, true); } //don't read groups to save space

        ifstream inFasta;
        params->util.openInputFile(params->fastafile, inFasta);
        set<int> lengths;
        vector<seqPNode*> alignSeqs;
        
        while (!inFasta.eof()) {

            if (params->m->getControl_pressed()) { inFasta.close(); break; }

            Sequence seq(inFasta);  params->util.gobble(inFasta);

            if (seq.getName() != "") {  //can get "" if commented line is at end of fasta file

                //no names file, you are identical to yourself
                int numReps = 1;
                if (params->hasCount) { numReps = ct.getNumSeqs(seq.getName()); }
                vector<int> clusteredIndexes; clusteredIndexes.push_back(alignSeqs.size());
                seqPNode* tempNode = new seqPNode(seq.getName(), seq.getAligned(), numReps, clusteredIndexes);
                alignSeqs.push_back(tempNode);
                lengths.insert(seq.getAligned().length());
                
            }
        }
        inFasta.close();

        params->length = *(lengths.begin());

        if (lengths.size() > 1) { params->align_method = "unaligned"; }
        else if (lengths.size() == 1) {  params->align_method = "aligned"; filterSeqs(alignSeqs, params->length, params->m); }

        //sort seqs by number of identical seqs
        sort(alignSeqs.begin(), alignSeqs.end(), comparePriorityAbundance);
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
        ofstream outAccnos;
        string outAccnosFileName = newfasta + ".temp";
        params->util.openOutputFile(outAccnosFileName, outAccnos);
        
        CountTable ct;
        if (params->countfile != "") {
            for (int i = 0; i < params->alignSeqs.size(); i++) {
                if (params->alignSeqs[i]->numIdentical != 0) {
                    outAccnos << params->alignSeqs[i]->name << endl;
                    ct.push_back(params->alignSeqs[i]->name, params->alignSeqs[i]->numIdentical);
                }
            }
        }
        outAccnos.close();
        
        if (params->countfile != "")  { ct.printTable(newname); }

        //use unique.seqs to create new name and fastafile
        string inputString = "fasta=" + params->fastafile + ", accnos=" + outAccnosFileName;
        params->m->mothurOut("/******************************************/\n");
        params->m->mothurOut("Running command: get.seqs(" + inputString + ")\n");
        
        Command* getCommand = new GetSeqsCommand(inputString);
        getCommand->execute();
        
        map<string, vector<string> > filenames = getCommand->getOutputFiles();
        
        delete getCommand;
        
        params->util.renameFile(filenames["fasta"][0], newfasta);
        
        params->m->mothurOut("/******************************************/\nDone.\n");
        
        params->util.mothurRemove(outAccnosFileName);

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
                
                if(util.isEqual(error_dist[0], -100)){ //construct the binomial distribution
                    
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
                        
                        if(!util.isEqual(a, 0) && !util.isEqual(b, 0)){
                            error_dist.push_back(choose * a * b);
                            if(error_dist[i] < 0.1 / max_abund){ break;    }
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
        
        long start = time(NULL);
        
        string numProcessors = current->getProcessors();
        
        string fileroot = outputdir + util.getRootName(util.getSimpleName(fastafile));
        map<string, string> variables;
        variables["[filename]"] = fileroot;
        
        string newCountFile = getOutputFileName("count",variables);
        string newMapFile = getOutputFileName("map",variables); //add group name if by group
        
        variables["[extension]"] = util.getExtension(fastafile);
        string newFastaFile = getOutputFileName("fasta", variables);
        outputNames.push_back(newFastaFile); outputTypes["fasta"].push_back(newFastaFile);
        
        if (countfile != "")    { outputNames.push_back(newCountFile); outputTypes["count"].push_back(newCountFile);    }
        
        if (bygroup) {
            //clear out old files
            ofstream outFasta; util.openOutputFile(newFastaFile, outFasta); outFasta.close();
            ofstream outCount; util.openOutputFile(newCountFile, outCount);  outCount.close();
            
            newMapFile = fileroot + "precluster.";
            string convolutedNamesFile = newCountFile + ".temp";
            
            vector<string> groups;
            map<string, vector<string> > group2Files;
            
            current->setMothurCalling(true);
            SequenceCountParser cparser(countfile, fastafile, nullVector);
            current->setMothurCalling(false);
            //cout << " groups = "<< cparser.getNamesOfGroups().size() << endl;
            groups = cparser.getNamesOfGroups();
            group2Files = cparser.getFiles();
            
            createProcessesGroups(group2Files, groups, convolutedNamesFile, newMapFile);
            
            string accnosFile;
            if (countfile != "") {  accnosFile = mergeGroupCounts(newCountFile, convolutedNamesFile); }
            
            printFasta(newFastaFile, accnosFile);
            util.mothurRemove(convolutedNamesFile);
            
            if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	}	 return 0; }
            
            m->mothurOut("It took " + toString(time(NULL) - start) + " secs to run pre.cluster.\n");
            
        }else {
            if (processors != 1) { m->mothurOut("When using running without group information mothur can only use 1 processor, continuing.\n");  processors = 1; }
            
            vector<string> groups;
            map<string, vector<string> > group2Files;
            preClusterData* params = new preClusterData(group2Files, fastafile, countfile, pc_method, align_method, clump, NULL, newMapFile, nullVector);
            params->setVariables(diffs, pc_method, align_method, align, match, misMatch, gapOpen, gapExtend, alpha, delta, error_rate, indel_prob, max_indels, error_dist);
            
            //reads fasta file and return number of seqs
            long long numSeqs = 0; params->alignSeqs = readFASTA(params, numSeqs); //fills alignSeqs and makes all seqs active
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
            
            print(newFastaFile, newCountFile, params);
            for (int i = 0; i < params->alignSeqs.size(); i++) {  delete params->alignSeqs[i]; } params->alignSeqs.clear();
            m->mothurOut("It took " + toString(time(NULL) - start) + " secs to cluster " + toString(numSeqs) + " sequences.\n");
        }
        
        current->setProcessors(numProcessors);
        
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

void PreClusterCommand::printFasta(string newFastaFileName, string accnosFile){
    try {
        string inputString = "fasta=" + fastafile + ", accnos=" + accnosFile;
        
        m->mothurOut("\n/******************************************/\n");
        m->mothurOut("Running command: get.seqs(" + inputString + ")\n");
        current->setMothurCalling(true);
        
        Command* getCommand = new GetSeqsCommand(inputString);
        getCommand->execute();
        
        map<string, vector<string> > filenames = getCommand->getOutputFiles();
        
        delete getCommand;
        current->setMothurCalling(false);
        m->mothurOut("/******************************************/\n");
        
        util.renameFile(filenames["fasta"][0], newFastaFileName);
    }
    catch(exception& e) {
        m->errorOut(e, "PreClusterCommand", "printFasta");
        exit(1);
    }
}
/**************************************************************************************************/

void printData(string group, preClusterData* params, map<string, string>& optionalNameMap){
    try {
        if ((params->hasCount) && (group == ""))  {
            params->newNName->write("Representative_Sequence\ttotal\n");
        }
        
        if (params->hasCount) {
            if (group != "") {
                
                for (int i = 0; i < params->alignSeqs.size(); i++) {
                    if (params->alignSeqs[i]->numIdentical != 0) {
                        params->newNName->write(group + '\t' + params->alignSeqs[i]->name + '\t' + toString(params->alignSeqs[i]->numIdentical) + '\n');
                    }
                }
            }
            else {
                for (int i = 0; i < params->alignSeqs.size(); i++) {
                    if (params->alignSeqs[i]->numIdentical != 0) {
                        params->newNName->write(params->alignSeqs[i]->name  + '\t' + toString(params->alignSeqs[i]->numIdentical) + '\n');
                    }
                }
            }
        }else {
            for (int i = 0; i < params->alignSeqs.size(); i++) {
                if (params->alignSeqs[i]->numIdentical != 0) {
                    string clusteredNames = "";
                    for (int j = 0; j < params->alignSeqs[i]->clusteredIndexes.size(); j++) {
                        int indexOfSeq = params->alignSeqs[i]->clusteredIndexes[j];
            
                        string repName = params->alignSeqs[indexOfSeq]->name;
                        string dupNames = "";
                        map<string, string>::iterator itDupName = optionalNameMap.find(repName);
                        if (itDupName != optionalNameMap.end()) {
                            dupNames = itDupName->second;
                            clusteredNames += dupNames + ",";
                        }
                    }
                    
                    clusteredNames = clusteredNames.substr(0,clusteredNames.length()-1); //remove last comma
                    
                    params->newNName->write(params->alignSeqs[i]->name + '\t' + clusteredNames + '\n');
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

bool fillWeighted(preClusterData* params, string fastafileName, string groupOrCountFile){
    try {
        set<int> lengths;
        
        map<string, int> counts;
        if (params->hasCount) {
            CountTable ct; ct.readTable(groupOrCountFile, false, true); //don't read groups because it only contains one group
            counts = ct.getNameMap();
            params->util.mothurRemove(groupOrCountFile);
        }
        
        ifstream in;
        params->util.openInputFile(fastafileName, in);
        
        while (!in.eof()) {
            
            if (params->m->getControl_pressed()) { break; }
            
            Sequence seq(in); params->util.gobble(in);
            
            if (seq.getName() != "") {
                
                map<string, int>::iterator it = counts.find(seq.getName());
                if (it != counts.end()) {
                    vector<int> clusteredIndexes; //don't use indexes for precluster with count file
                    if (!params->hasCount) {  clusteredIndexes.push_back(params->alignSeqs.size());  }
                    
                    seqPNode* tempNode = new seqPNode(seq.getName(), seq.getAligned(), it->second, clusteredIndexes);
                    params->alignSeqs.push_back(tempNode);
                    lengths.insert(seq.getAligned().length());
                }
            }
        }
        in.close();
        
        params->length = *(lengths.begin());
        
        if (lengths.size() > 1) { params->align_method = "unaligned";  return false; } //unaligned
        else if (lengths.size() == 1) {  params->align_method = "aligned"; filterSeqs(params->alignSeqs, params->length, params->m);  return true; } //aligned
        
        return true;

    }
    catch(exception& e) {
        params->m->errorOut(e, "PreClusterCommand", "fillWeighted");
        exit(1);
    }
}
/**************************************************************************************************/

long long driverGroups(preClusterData* params){
    try {
        long long numSeqs = 0;
    
        //precluster each group
        for (map<string, vector<string> >::iterator it = params->parsedFiles.begin(); it != params->parsedFiles.end(); it++) {
            if (params->m->getControl_pressed()) { return numSeqs; }
            
            string thisGroup = it->first;
            params->m->mothurOut("\nProcessing group " + thisGroup + ":\n");
            
            time_t start = time(NULL);
            bool aligned = false;
            
            string thisGroupsFasta = it->second[0];
            
            map<string, string> thisGroupsNameMap;
            if (params->hasCount)          {
                string thisGroupsCount = it->second[1];
                
                aligned = fillWeighted(params, thisGroupsFasta, thisGroupsCount);
                
                params->util.mothurRemove(thisGroupsCount);
            }
            params->util.mothurRemove(thisGroupsFasta);
            
            //sort seqs by number of identical seqs
            sort(params->alignSeqs.begin(), params->alignSeqs.end(), comparePriorityAbundance);
            
            long long num = params->alignSeqs.size();
            numSeqs += num;
            
            if (params->m->getControl_pressed()) {   return 0; }
            
            if (params->align_method == "aligned") {
                if (params->diffs > params->length) {
                    params->m->mothurOut("[ERROR]: diffs is greater than your sequence length.\n");
                    params->m->setControl_pressed(true);
                    return 0;
                }
            }
            
            string extension = thisGroup+".map";
            long long count = process(thisGroup+"\t", params->newMName+extension, params);
            
            if(params->pc_method != "deblur"){
                params->outputNames.push_back(params->newMName+extension);
                params->outputTypes["map"].push_back(params->newMName+extension);
            }
            
            if (params->m->getControl_pressed()) {  return 0; }
            
            params->m->mothurOut("Total number of sequences before pre.cluster was " + toString(num) + ".\n");
            params->m->mothurOut("pre.cluster removed " + toString(count) + " sequences.\n\n");
                
            printData(thisGroup, params, thisGroupsNameMap);
            
            for (int i = 0; i < params->alignSeqs.size(); i++) {  delete params->alignSeqs[i]; } params->alignSeqs.clear();
            
            params->m->mothurOut("It took " + toString(time(NULL) - start) + " secs to cluster " + toString(num) + " sequences.\n");
        }
        
        return numSeqs;
    }
    catch(exception& e) {
        params->m->errorOut(e, "PreClusterCommand", "driverGroups");
        exit(1);
    }
}

/**************************************************************************************************/
//only called with count table including groups
string PreClusterCommand::mergeGroupCounts(string newcount, string newname){
	try {
        m->mothurOut("\nDeconvoluting count table results...\n");
        
        CountTable ct; vector<string> groups;
        ct.testGroups(countfile, groups);
        ct.readTable(countfile, false, true); //read table no groups
        for (int i = 0; i < groups.size(); i++)  { ct.addGroup(groups[i]); } //add groups
        ct.zeroOutTable();
        
        ifstream inNames;
        util.openInputFile(newname, inNames);
        
        /*
         newname looks like:
         
        groupName seqName seqCountForGroup
         
         FDF6  seq1  35
         
         seq1 has an abundance of 35 in group FDF6
         
         */
        
        time_t start = time(NULL);
        long long count = 0;
        
        string group, unique_sequence;
        int numDups;
        
        //build table
        vector<string> namesOfSeqs;
        while (!inNames.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            inNames >> group; util.gobble(inNames);
            inNames >> unique_sequence; util.gobble(inNames);
            inNames >> numDups; util.gobble(inNames);
            
            ct.setAbund(unique_sequence, group, numDups); count++;
            namesOfSeqs.push_back(unique_sequence);
            
            //report progress
            if((count) % 1000 == 0){	m->mothurOutJustToScreen(toString(count) + "\n"); 		}
        }
        //report progress
        if((count) % 1000 != 0){	m->mothurOutJustToScreen(toString(count) + "\n"); 		}
        
        inNames.close();
        
        m->mothurOut("It took " + toString(time(NULL) - start) + " secs to merge " + toString(count) + " sequences group data.");
        start = time(NULL);
        
        ct.printTable(newcount);
        util.mothurRemove(newname);
        
        ofstream outAccnos; util.openOutputFile(newname, outAccnos);
        for (int i = 0; i < namesOfSeqs.size(); i++) { outAccnos << namesOfSeqs[i] << endl; }
        outAccnos.close();

		return newname;
	}
	catch(exception& e) {
		m->errorOut(e, "PreClusterCommand", "mergeGroupCounts");
		exit(1);
	}
}
/**************************************************************************************************/

void PreClusterCommand::createProcessesGroups(map<string, vector<string> >& parsedFiles, vector<string> groups, string newNName, string newMFile) {
  try {
    
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
    vector<std::thread*> workerThreads;
    vector<preClusterData*> data;
      
    auto synchronizedNameFile = std::make_shared<SynchronizedOutputFile>(newNName);

    //Lauch worker threads
    for (int i = 0; i < processors-1; i++) {
      OutputWriter* threadNameWriter = new OutputWriter(synchronizedNameFile);
        
        vector<string> thisGroups;
        map<string, vector<string> > thisGroupsParsedFiles;
        for (int j = lines[i+1].start; j < lines[i+1].end; j++) {
            
            map<string, vector<string> >::iterator it = parsedFiles.find(groups[j]);
            if (it != parsedFiles.end()) {
                thisGroupsParsedFiles[groups[j]] = (it->second);
                thisGroups.push_back(groups[j]);
            }
            else { m->mothurOut("[ERROR]: missing files for group " + groups[j] + ", skipping\n"); }
        }
        
      preClusterData* dataBundle = new preClusterData(thisGroupsParsedFiles, fastafile, countfile, pc_method, align_method, clump, threadNameWriter, newMFile, thisGroups);
      dataBundle->setVariables(diffs, pc_method, align_method, align, match, misMatch, gapOpen, gapExtend, alpha, delta, error_rate, indel_prob, max_indels, error_dist);
      data.push_back(dataBundle);

      workerThreads.push_back(new std::thread(driverGroups, dataBundle));
    }
    OutputWriter* threadNameWriter = new OutputWriter(synchronizedNameFile);
      
      vector<string> thisGroups;
      map<string, vector<string> > thisGroupsParsedFiles;
      for (int j = lines[0].start; j < lines[0].end; j++) {
          
          map<string, vector<string> >::iterator it = parsedFiles.find(groups[j]);
          if (it != parsedFiles.end()) {
              thisGroupsParsedFiles[groups[j]] = (it->second);
              thisGroups.push_back(groups[j]);
          }
          else { m->mothurOut("[ERROR]: missing files for group " + groups[j] + ", skipping\n"); }
      }
    preClusterData* dataBundle = new preClusterData(thisGroupsParsedFiles, fastafile, countfile, pc_method, align_method, clump, threadNameWriter, newMFile, thisGroups);
    dataBundle->setVariables(diffs, pc_method, align_method, align, match, misMatch, gapOpen, gapExtend, alpha, delta, error_rate, indel_prob, max_indels, error_dist);

    driverGroups(dataBundle);

    outputNames.insert(outputNames.end(), dataBundle->outputNames.begin(), dataBundle->outputNames.end());
    for (itTypes = dataBundle->outputTypes.begin(); itTypes != dataBundle->outputTypes.end(); itTypes++) {
        outputTypes[itTypes->first].insert(outputTypes[itTypes->first].end(), itTypes->second.begin(), itTypes->second.end());
    }

    for (int i = 0; i < processors-1; i++) {
      workerThreads[i]->join();

      delete data[i]->newNName;

      outputNames.insert(outputNames.end(), data[i]->outputNames.begin(), data[i]->outputNames.end());
      for (itTypes = data[i]->outputTypes.begin(); itTypes != data[i]->outputTypes.end(); itTypes++) {
          outputTypes[itTypes->first].insert(outputTypes[itTypes->first].end(), itTypes->second.begin(), itTypes->second.end());
      }

      delete data[i];
      delete workerThreads[i];
    }
    delete threadNameWriter;
    delete dataBundle;

  }
  catch(exception& e) {
      m->errorOut(e, "PreClusterCommand", "createProcessesGroups");
      exit(1);
  }
}

/**************************************************************************************************/
