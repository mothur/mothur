#ifndef PROGRESS_H
#define PROGRESS_H

#include <string>
using namespace std;

class Progress {
	
public:
	Progress(string, int);
	void update(int);
	void finish();
	
private:
	int nTicks;
	int finalPos;	
};

#endif
