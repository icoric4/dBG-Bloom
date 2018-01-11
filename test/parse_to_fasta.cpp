#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;



int main (int argc, char *argv[]) {
	if (argc != 2) {
		cout << "Please enter a name of file to parse!" << endl;
		exit(1);
	}
	int i = 0;
	string s = argv[1], d;
	ifstream f(s);
	ofstream of(s.substr(0,s.find(".")) + ".fasta");
	string line;
	
	if (f.is_open()) {
		while ( getline (f,line) )	{
			of << ">r" + to_string(++i) + "\n";
			of << line << "\n";
		}
    	f.close();
    	of.close();
  	} else {
  		cout << "File can not be opened!" << endl;
  	}
}