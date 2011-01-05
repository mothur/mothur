#ifndef PCOACOMMAND_H
#define PCOACOMMAND_H

/*
 *  pcoacommand.h
 *  Mothur
 *
 *  Created by westcott on 1/4/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */
 
#include "command.hpp"


/*****************************************************************/
class PCOACommand : public Command {
	
public:
	PCOACommand(string);	
	PCOACommand();
	~PCOACommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();	
	void help();
	
private:

	bool abort, metric;
	string phylipfile, columnfile, namefile, format, filename, fbase, outputDir;
	float cutoff, precision;
	vector<string> outputNames;
	map<string, vector<string> > outputTypes;
	
	void get_comment(istream&, char, char);
	int read_phylip(istream&, int, vector<string>&, vector<vector<double> >&);
	void read(string, vector<string>&, vector<vector<double> >&);
	double pythag(double, double);
	void matrix_mult(vector<vector<double> >, vector<vector<double> >, vector<vector<double> >&);
	void recenter(double, vector<vector<double> >, vector<vector<double> >&);
	void tred2(vector<vector<double> >&, vector<double>&, vector<double>&);
	void qtli(vector<double>&, vector<double>&, vector<vector<double> >&);
	void output(string, vector<string>, vector<vector<double> >&, vector<double>);
	vector< vector<double> > calculateEuclidianDistance(vector<vector<double> >&, int);
	double calcPearson(vector<vector<double> >&, vector<vector<double> >&);
	
};
	
/*****************************************************************/
	
#endif

