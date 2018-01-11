#include <iostream>
#include <string>
#include <fstream>
using namespace std;



int main(int argc, char *argv[]) {

	if (argc != 2) {
		cout << "Please enter a name of file with data count!" << endl;
		exit(1);
	}
	ifstream ifs(argv[1]);
	string s = argv[1];
	ofstream of(s.substr(0,s.find("."))+"S.txt");

	if (ifs.is_open())	{
		char c;
		string seq;
		int count;
		while ( ifs >> c >> count >> seq)	{
			if (count >= 1) {
				of << seq << "\n";
			}
		}
    	ifs.close();
    	of.close();
  	} else {
  		cout << "File can not be opened!" << endl;
  	}
}