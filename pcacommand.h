#ifndef PCACOMMAND_H
#define PCACOMMAND_H

/*
 *  pcacommand.h
 *  Mothur
 *
 *  Created by westcott on 1/4/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */
 
#include "command.hpp"


/*****************************************************************/
class PCACommand : public Command {
	
public:
	PCACommand(string);	
	~PCACommand();
	int execute();	
	void help();
	
private:

	bool abort;
	string phylipfile, columnfile, namefile, format, filename, fbase;
	float cutoff, precision;
	
	void get_comment(istream&, char, char);
	void read_mega(istream&, vector<string>&, vector<vector<double> >&);
	void read_phylip(istream&, int, vector<string>&, vector<vector<double> >&);
	void read(string, vector<string>&, vector<vector<double> >&);
	double pythag(double, double);
	void matrix_mult(vector<vector<double> >, vector<vector<double> >, vector<vector<double> >&);
	void recenter(double, vector<vector<double> >, vector<vector<double> >&);
	void tred2(vector<vector<double> >&, vector<double>&, vector<double>&);
	void qtli(vector<double>&, vector<double>&, vector<vector<double> >&);
	void output(string, vector<string>, vector<vector<double> >, vector<double>);
	void print_matrix(vector<vector<double> >);
	
};
	
/*****************************************************************/
	
#endif

