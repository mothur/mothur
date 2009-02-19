using namespace std;

#include "nameassignment.hpp"

//**********************************************************************************************************************

NameAssignment::NameAssignment(string nameMapFile){
	
	openInputFile(nameMapFile, fileHandle);
	
}

//**********************************************************************************************************************

void NameAssignment::readMap(int colA, int colB){
	try{
		string firstCol, secondCol, skip;
	//	int index = 0;
	
		int skipNCols = colB-colA-1;
	
		map<string, string> data;
	
		while(fileHandle){
			fileHandle >> firstCol;				//read from first column
		
			for(int i=0;i<skipNCols;i++){		//allows for anticipated file format
				fileHandle >> skip;
			}
		
			fileHandle >> secondCol;			//read from second column
		
			data[firstCol] = secondCol;			//store data in map
		
			gobble(fileHandle);
		}
		fileHandle.close();
	
		int rowIndex = 0;
		map<string, string>::iterator it = data.begin();
		for(it;it!=data.end();it++){
			list.push_back(it->second);		//adds data's value to list
			(*this)[it->first] = rowIndex;
			rowIndex++;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the NameAssignment class Function readMap. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the NameAssignment class function readMap. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

//**********************************************************************************************************************

ListVector NameAssignment::getListVector(void){

	return list;
	
}

//**********************************************************************************************************************

void NameAssignment::print(void){
	try {
		map<string,int>::iterator it = (*this).begin();
		for(it;it!=(*this).end();it++){
			cout << it->first << '\t' << it->second << endl;  //prints out keys and values of the map this.
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the NameAssignment class Function print. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the NameAssignment class function print. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

//**********************************************************************************************************************

int NameAssignment::get(string key){
	
	return	(*this)[key];	

}

//**********************************************************************************************************************

