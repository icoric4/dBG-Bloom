#ifndef DBG_BLOOM_DBG_HPP
#define DBG_BLOOM_DBG_HPP

#include <iostream>
#include <string>
#include <fstream>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <list>
#include <ctime>
#include <sstream>
#include "KmerBloomFilter.hpp"



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
	void test_size();

private:
	string input_file_name;
	vector<string> S;
	//bloom_filter filter;
	KmerBloomFilter filter;
	unordered_set<string> cFP;
	unordered_set<string> marked;
	int k_size;
	vector<string> contigs;

	void create_BF();
	vector<string> compute_P();
	vector<string> branch(string s, bool right);
	string canonical (string s1);
	string complement(string c);
};

#endif //DBG_HPP

