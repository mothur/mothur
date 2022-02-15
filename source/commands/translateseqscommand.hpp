//
//  translateseqscommand.hpp
//  Mothur
//
//  Created by Sarah Westcott on 11/8/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#ifndef translateseqscommand_hpp
#define translateseqscommand_hpp

#include "command.hpp"
#include "sequence.hpp"
#include "protein.hpp"

/*
 This command would take...

  * DNA sequences and translate it to an amino acid sequence

 By default it would use the first frame
 The user could also specify the frame (1, 2, 3, -1, -2, -3) or possibly use all 6 frames
 Another option would be stop=T/F. If T, then if the translation hits a stop codon, it stops before that codon. If F, it returns the full translation with a * as the stop codon
 Output as *.aa#.fasta where # is the frame
 
  * Amino acid sequences and translate it to a DNA sequence

 Because of degeneracies there will be non-ATGC IUPAC codes in the output sequence
 Output as *.dna.fasta
 
  * Unaligned DNA and unaligned/aligned Amino acid sequences

 Back translate the amino acid sequence to the DNA sequence so that the DNA is aligned. This should result in the DNA bases being clustered in groups of 3 corresponding to each amino acid codon
 Hopefully the DNA sequence and the amino acid sequence will be in the same frame
 Output alignment as *.dna.align
 
 */

/**************************************************************************************************/

class TranslateSeqsCommand : public Command {
    
public:
    TranslateSeqsCommand(string);
    ~TranslateSeqsCommand(){}
    
    vector<string> setParameters();
    string getCommandName()            { return "tranlate.seqs";              }
    string getCommandCategory()        { return "Sequence Processing";        }
    
    string getHelpString();
    string getCommonQuestions();
    string getOutputPattern(string);
    string getCitation() { return "http://www.mothur.org/wiki/translate.seqs"; }
    string getDescription()        { return "tranlate dna to amino acids or amino acids to dna"; }
    
    int execute();
    void help() { m->mothurOut(getHelpString()); }
    
private:
    bool abort, stop, aminoAligned, dnaAligned;
    string fastafile, aminofile;
    int processors;
    vector<string> outputNames;
    vector<int> frames;
    vector<linePair> lines;
    vector<linePair> aLines;
    
    bool setLines(); //returns true if error free
    void translateDNAtoAmino();
    void alignDNAAmino();
    double createProcessesTranslateDNAtoAminoAcids(string, vector<linePair>, int);
    double createProcessesAlign(string);
};

//**********************************************************************************************************************
struct translateSeqsStruct {
    OutputWriter* outputWriter;
    string inputFilename;
    bool stop;
    int frame;
    double numSeqs;
    
    linePair filePos;
    MothurOut* m; Utils util;
    
    translateSeqsStruct (linePair fP, OutputWriter* oFName, string fname, bool st, int f) {
        
        //passed in
        filePos.start = fP.start;
        filePos.end = fP.end;
        outputWriter = oFName;
        inputFilename = fname;
        frame = f;
        stop = st;
        
        //initialized
        numSeqs = 0;
        m = MothurOut::getInstance();
    }
    ~translateSeqsStruct() {}
};
//**********************************************************************************************************************
struct alignAminoStruct {
    OutputWriter* outputWriter;
    string fastaFilename, aminoFilename;
    bool stop, aminoAligned, dnaAligned;
    double numSeqs;
        
    linePair fastaPos;
    linePair aminoPos;
    MothurOut* m; Utils util;
        
    alignAminoStruct (linePair fP, linePair aP, OutputWriter* oFName, string fname, string aname, bool st, bool da, bool aa) {
            
        //passed in
        fastaPos.start = fP.start;
        fastaPos.end = fP.end;
        aminoPos.start = aP.start;
        aminoPos.end = aP.end;
        outputWriter = oFName;
        fastaFilename = fname;
        aminoFilename = aname;
        dnaAligned = da;
        aminoAligned = aa;
        stop = st;
            
        //initialized
        numSeqs = 0;
        m = MothurOut::getInstance();
    }
    ~alignAminoStruct() {}
};
//**********************************************************************************************************************


#endif /* translateseqscommand_hpp */
