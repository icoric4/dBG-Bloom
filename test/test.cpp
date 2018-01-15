#include "../src/dbg_bloom.hpp"
#include <ctime>

using namespace std;

/*
 * Program that test the implementatio of dbg class.
 * It takes 2 arguments from command line: name of file
 * with solid k-mers and name of file from which the 
 * k-mers was created
 */
int main (int argc, char *argv[]) {

	dbg d(argv[1]);
	d.compute_cFP();
	d.traverse_graph();
	ifstream ifs("contigs.txt");
	ifstream ifs2(argv[2]);
	string dna;
	getline(ifs2, dna); // read >r1
	getline(ifs2, dna);
	string s;
	bool read = false;
	while (ifs >> s) {
		read = true;
		if (dna.find(s)==string::npos && dna.find(d.reverse_complement(s))==string::npos) {
			bool passed = false;
			for (int i = 1; i != 2; ++i) { 	// maybe contig came to the end of original dna and unmarked kmer existed, 
											// such that it managed to add a nucleotide to its end, this can happen more times, 
											// but eventually it stops (after 5 times we hope) 
				string ss = s.substr(0,s.size()-i);
				if (dna.find(ss)!=string::npos || dna.find(d.reverse_complement(ss))!=string::npos) {
					passed = true;
					break;
				}
				ss = s.substr(i);
				if (dna.find(ss)!=string::npos || dna.find(d.reverse_complement(ss))!=string::npos) {
					passed = true;
					break;
				}

			}
			if (!passed) {
				cout << "FAILED" << endl;
				exit(0);
			}
		}
	}

	if (read) {
		cout << "PASS" << endl;
	} else {
		cout << "contigs.txt empty" << endl; 
	}

	return 0;
}