#ifndef DBG_HPP
#define DBG_HPP
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
#include <vector>
#include <math.h>
#include <iomanip>
#include <cstdio>



using namespace std;

/*
 * Class that builds space efficient and exact de Brujin graph
 * based on Bloom filter
 */
class dbg {
public: 

	dbg(string name);
	void compute_cFP(unsigned int M=1e9);
	void traverse_graph(string name="contigs.txt", int depth_bound=500, int width_bound=20);
	string reverse_complement (string s);
	void test_size();

private:
	string input_file_name;
	vector<string> S;
	KmerBloomFilter filter;
	unordered_set<string> cFP;
	unordered_set<string> marked;
	int k_size;
	vector<string> contigs;
	unsigned bloom_size;

	void create_BF();
	vector<string> compute_P();
	vector<char> branch(string s, bool right);
	string canonical (string s1);
	char complement(char c);
    vector<string> construct_contigs(vector<vector<vector<char> > > mem, int depth, bool right);
};

#endif //DBG_HPP

