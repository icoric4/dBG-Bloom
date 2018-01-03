#include <iostream>
#include <string>
#include <fstream>
#include <set>
#include "bloom_filter.hpp"
#include <map>
#include <list>

using namespace std;

class dbg {
public: 

	string input_file_name;
	set<string> S;
	bloom_filter filter;
	set<string> P;
	set<string> cFP;
	vector<string> nodes;
	map<string,bool> marked;
	int k_size;
	set<string> contigs;

	dbg(string name) {
		ifstream ifs(name);
		string s;
		while (ifs >> s) {
			S.insert(s);
		}
		vector<bool> bla(S.size(), false);
		k_size = (*S.begin()).size();
		create_BF();
	}

	void create_BF() {
		bloom_parameters parameters;
		parameters.projected_element_count = S.size();
		parameters.false_positive_probability = 0.0001;
		parameters.random_seed = 0xA5A5A5A5;
		parameters.compute_optimal_parameters();
		filter = bloom_filter(parameters);
		for (set<string>::iterator it = S.begin(); it != S.end(); ++it)	{
			filter.insert(*it);
			filter.insert(reverse_complement(*it));
		}
	}

	void compute_P() {
		string h[] = {"A","C","G","T"};
		set<string> nucleotides(h, h+4);
		for (set<string>::iterator it = S.begin(); it != S.end(); ++it)	{
			string prefix = (*it).substr(0, (*it).size()-1);
			string sufix = (*it).substr(1);
			for (set<string>::iterator nucl = nucleotides.begin(); nucl != nucleotides.end(); ++nucl) {

				if (filter.contains(sufix + *nucl)) {
					P.insert(canonical(sufix + *nucl));
				}

				if (filter.contains(*nucl + prefix)) {
					P.insert(canonical(*nucl + prefix));
				}
			}
		}
		 cout << S.size() << " " << P.size() << endl;
	}

	void compute_cFP(unsigned int M) {
		compute_P();
		set<string> D(P);
		set<string>::iterator it = S.begin();
		while(it != S.end()) {
			set<string> Pi;
			set<string> Dn;

			while (Pi.size() < M && it != S.end()) {
				Pi.insert(*it);
				it++;
			}

			for (set<string>::iterator m = D.begin(); m != D.end(); ++m) {
				if (Pi.count(canonical(*m)) == 0) {
					Dn.insert(canonical(*m));
				} //else izbrisi
			}

			D = set<string>(Dn);
		}


		cFP = set<string>(D);
		cout << cFP.size() << endl;
		vector<string> bla(S.begin(), S.end());
		nodes = bla;

	}

	void print1() {
		int i = 0;
		for (set<string>::iterator it = P.begin(); it != P.end(); ++it)	{
			i++;
			if (filter.contains(*it)) {
				cout << i << ". "  << "BF contains: " << *it << endl;
			}
		}
	}

	void print2() {
		int i = 0;
		for (set<string>::iterator it = S.begin(); it != S.end(); ++it)	{
			i++;
			if (filter.contains(*it)) {
				cout << i << ". "  << "BF contains: " << *it << endl;
			}
		}
	}

	void traverse_graph() {

		int index = 0;

		while (index < nodes.size()) {
			set<string> front;
			front.insert(nodes[index]);
			int depth = 0;
			int reason;

			do {
				// branch to next front
				// cout << depth << ": FRONT" << " "; 
				// for (set<string>::iterator it = front.begin(); it != front.end(); ++it)
				// 	cout << *it << " ";
				// cout << endl;
				set<string> new_front;
				for (set<string>::iterator it = front.begin(); it != front.end(); ++it) {
					if (!marked.count((*it).substr((*it).size()-k_size))) {
						continue;
					}
					set<string> ret = branch(*it);
					new_front.insert(ret.begin(), ret.end());
					marked.insert(pair<string,bool>((*it).substr((*it).size()-k_size),true));
					marked.insert(pair<string,bool>((*it).substr(0,k_size),true));
					marked.insert(pair<string,bool>(reverse_complement((*it).substr((*it).size()-k_size)),true));
					marked.insert(pair<string,bool>(reverse_complement((*it).substr(0,k_size)),true));
				}
				if(new_front.size() == 0) {
					reason = 2;
					break;
				}
				depth++;
				front = new_front;
				// cout << depth <<": NEW_FRONT" << " "; 
				// for (set<string>::iterator it = front.begin(); it != front.end(); ++it)
				// 	cout << *it << " ";
				// cout << endl;
				// if (depth == 40) 
				// 	exit(1);
				int width = new_front.size();

				if (depth > 500 || width > 20) {
					reason = 0;
					break;
				}
				if (new_front.size() == 1 && depth > 2*k_size+1) {
					
					reason = 1;
					break;
				}
			} while (1);

			for (set<string>::iterator it = front.begin(); it != front.end(); ++it) {
				if ((*it).size()>2*k_size+1) {
					contigs.insert(*it);
				}
			}

			index++;
			for (int i = index; i < nodes.size(); ++i) {
				if(!marked[nodes[i]]) {
					index = i;
					break;
				}
			} 
		}
	}

	set<string> branch(string s) {

		char bla[4] = {'A', 'C', 'G', 'T'};
		list<char> l(bla, bla + sizeof(bla) / sizeof(char));
		set<string> ret;
		// cout << "1: " << s;
		for (list<char>::iterator c = l.begin(); c != l.end(); ++c) {
			if (filter.contains(*c+s.substr(0,k_size-1)) && 
				!cFP.count(canonical(*c+s.substr(0,k_size-1)))&&
				!marked.count(*c+s.substr(0,k_size-1)) ) {
				ret.insert(*c+s);
			}
			if (filter.contains(s.substr(s.size()-k_size+1)+*c) && 
				!cFP.count(canonical(s.substr(s.size()-k_size+1)+*c)) &&
				!marked.count(s.substr(s.size()-k_size+1)+*c)) {
				ret.insert(s+*c);
			}	
		}
		// cout << " 2: ";
		// for (auto r: ret) {
		// 	cout << r << " ";
		// }
		// cout << endl;
		return ret;
	}

	string canonical (string s1) {
		string s2 = reverse_complement(s1);
		return s1 < s2 ? s1 : s2;
	}
	string reverse_complement (string s) {
		string tmp = "";
		for(string::iterator it = s.begin(); it != s.end(); ++it) {
			tmp = complement(string(1,*it)) + tmp;
		}
		return tmp;
	}

	string complement(string c) {
		if (c == "C")
			return "G";
		if (c == "T") 
			return "A";
		if (c == "A") 
			return "T";
		if (c == "G")
			return "C";
		return "";
	}
};


