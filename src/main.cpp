#include "dbg_bloom.hpp" 
#include "SimpleOpt.h"
#include <iomanip>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

using namespace std;

std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
            result += buffer.data();
    }
    return result;
}


void help () {
	cout << "Specify -in option." << endl;
	cout << "[dbg options]" << endl;

	cout << setw (10) << "-in : " << "Path to the input file. File must be in fasta format." << endl;
	cout << setw (10) << "-out : " << "Name of the output file. (default: contigs.txt)" << endl;
	cout << setw (10) << "-k : " << "k-mer size. (default: 31)" << endl;
	cout << setw (10) << "--help : " << "Option usage." << endl;

}

int main(int argc, char *argv[]) {

	clock_t start = clock();
	CSimpleOpt::SOption g_rgOptions[] = {
	{ 0, "-in", SO_REQ_SEP },     // "-a"
	{ 1, "-k", SO_REQ_SEP },     // "-b"
	{ 2, "-out", SO_REQ_SEP },  // "-f ARG"
	{ 3, "--help", SO_NONE }, // "--help"
	{ 4, "-depth-bound", SO_REQ_SEP},
	{ 5, "-width-bound", SO_REQ_SEP},
	{ 6, "-abundance-min", SO_REQ_SEP},
	SO_END_OF_OPTIONS                // END
	};
	CSimpleOpt args(argc, argv, g_rgOptions);
	int k = 31;
	string out = "contigs.txt";
	string in = "";
	int depth = 500;
	int width = 20;
	int abundance_min = 3;
	while (args.Next()) {
		if (args.LastError() == SO_SUCCESS) {
			// handle option: use OptionId(), OptionText() and OptionArg()
			switch (args.OptionId()) {
				case 0 : in = args.OptionArg();
					break;
				case 1 : k = atoi(args.OptionArg());
					break;
				case 2 : out = args.OptionArg();
					break;
				case 3 : help(); exit(0);
					break;
				case 4 : depth = atoi(args.OptionArg());
					break;
				case 5 : width = atoi(args.OptionArg());
					break;	
				case 6 : abundance_min = atoi(args.OptionArg());
					break;	
			}
		} else {
			cout << "ERROR" << endl << "For options usage, please see help using --help option." << endl;
			exit(1); 
		}
	}
	if (in == "") {
		cout << "ERROR" << endl << "For options usage, please see help using --help option." << endl;
		exit(1); 
	}

	string cmd = "./jellyfish count -m " + to_string(k) +" -s 100M -t 10 -C -L " 
									+to_string(abundance_min)+" " + in +" -o tmp.jf";
	exec(cmd.c_str());
	cmd = "./jellyfish dump tmp.jf > tmp.fa";
	exec(cmd.c_str());
	cmd = "rm tmp.jf";
	exec(cmd.c_str());
	cmd = "./../test/create_s tmp.fa";
	exec(cmd.c_str());
	cmd = "rm tmp.fa";
	exec(cmd.c_str());
	
	dbg d("tmpS.txt");
	d.compute_cFP();
  	cmd = "rm tmpS.txt";
	exec(cmd.c_str());

	d.traverse_graph(out, depth, width);
	clock_t end = clock();
	double time = double(end - start) / CLOCKS_PER_SEC;
	cout << "Statistics: " << endl;
	cout << setw (30) << "Time : " << time << " s." << endl;
	d.test_size();

}








