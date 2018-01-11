#ifndef DBG_BLOOM_DBG_HPP
#define DBG_BLOOM_DBG_HPP

#include <iostream>
#include <string>
#include <fstream>
#include <set>
#include "bloom_filter.hpp"
#include <map>
#include <list>

using namespace std;


/* 
 * Class that builds space efficient and exact de Brujin graph
 * based on Bloom filter
 */
class dbg {
public: 

	dbg(string name);
	void compute_cFP(unsigned int M=1e9);
	void traverse_graph(string name="contigs.txt");
	string reverse_complement (string s);

private:
	string input_file_name;
	set<string> S;
	bloom_filter filter;
	set<string> P;
	set<string> cFP;
	vector<string> nodes;
	set<string> marked;
	int k_size;
	set<string> contigs;
	void create_BF();
	void compute_P();
	set<string> branch(string s, bool left, bool right);
	string canonical (string s1);
	string complement(string c);
};

#endif //DBG_HPP

