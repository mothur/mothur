

#include "nameassignment.hpp"

//**********************************************************************************************************************

NameAssignment::NameAssignment(string nameMapFile){
	m = MothurOut::getInstance();
	openInputFile(nameMapFile, fileHandle);
	
}

//**********************************************************************************************************************

void NameAssignment::readMap(){
	try{
		string firstCol, secondCol, skip;
	//	int index = 0;
	
		
		map<string, int>::iterator itData;
		int rowIndex = 0;
		
		while(fileHandle){
			fileHandle >> firstCol;				//read from first column
			fileHandle >> secondCol;			//read from second column
						
			itData = (*this).find(firstCol);
			if (itData == (*this).end()) {
			
				(*this)[firstCol] = rowIndex++;
				list.push_back(secondCol);		//adds data's value to list
				reverse[rowIndex] = firstCol;
				
			}else{	m->mothurOut(firstCol + " is already in namesfile. I will use first definition."); m->mothurOutEndLine();  }
			
			gobble(fileHandle);
		}
		fileHandle.close();
	
	}
	catch(exception& e) {
		m->errorOut(e, "NameAssignment", "readMap");
		exit(1);
	}
}
//**********************************************************************************************************************
void NameAssignment::push_back(string name) {
	try{
	
		int num = (*this).size();
		(*this)[name] = num;
		reverse[num] = name;
		
		list.push_back(name);
	}
	catch(exception& e) {
		m->errorOut(e, "NameAssignment", "push_back");
		exit(1);
	}
}

//**********************************************************************************************************************

ListVector NameAssignment::getListVector(void){

	return list;
	
}

//**********************************************************************************************************************

void NameAssignment::print(ostream& out){
	try {
		map<string,int>::iterator it;
//cout << (*this).size() << endl;
		for(it = (*this).begin(); it!=(*this).end(); it++){
			out << it->first << '\t' <<  it->second << endl;  //prints out keys and values of the map this.
		}
	}
	catch(exception& e) {
		m->errorOut(e, "NameAssignment", "print");
		exit(1);
	}
}

//**********************************************************************************************************************

int NameAssignment::get(string key){
	
	return	(*this)[key];	

}
//**********************************************************************************************************************

string NameAssignment::get(int key){
	
	return	reverse[key];	

}
//**********************************************************************************************************************

