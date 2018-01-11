#include "dbg_bloom.hpp"

/* 
 * Constructor of dbg class.
 * Argument name represents the name of the file which contains solid k-mers.
 * Set S which holds all solid k-mers (true positives) is olso created.
 */
dbg::dbg(string name) {
		ifstream ifs(name);
		string s;
		while (ifs >> s) {
			S.insert(s);
		}
		vector<bool> bla(S.size(), false);
		k_size = (*S.begin()).size();
		create_BF();
	}

/* 
 * Method for creation of bloom filter and adding all solid k-mers in the filter.
 */
void dbg::create_BF() {
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

/* 
 * Method for calculation of set P which represents the set of all extensions of
 * set S for which the Bloom filter answers yes.
 */
void dbg::compute_P() {
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

/* 
 * Method for computation of set cFP which is
 * formally defined as cFP = P\S
 */
void dbg::compute_cFP(unsigned int M) {
	compute_P();
	set<string> D(P);
	set<string>::iterator it = S.begin();
	//int i = 0;
	while(it != S.end()) {
		set<string> Pi;
		set<string> Dn;
		//cout << i << "/" << S.size() << endl;
		while (Pi.size() < M && it != S.end()) {
			Pi.insert(*it);
			it++;
		//	++i;
		}

		for (set<string>::iterator m = D.begin(); m != D.end(); ++m) {
			if (Pi.count(canonical(*m)) == 0) {
				Dn.insert(canonical(*m));
			} //else delete
		}

		D = set<string>(Dn);
	}

	cFP = set<string>(D);
	cout << cFP.size() << endl;
	vector<string> tmp(S.begin(), S.end());
	nodes = tmp;
}

/* 
 * Method for graph traversal. It starts from a solid k-mer 
 * and, using Bloom filter and cFP structure, looks for its neighbours
 * using bounded-breadth, bounded-depth BFS
 */
void dbg::traverse_graph() {

	int index = 0;

	while (index < nodes.size()) {
		set<string> front;
		front.insert(nodes[index]);
		int depth = 0;
		int reason;

		do {
			set<string> new_front;
			for (set<string>::iterator it = front.begin(); it != front.end(); ++it) {
				bool l = true, r = true;

				if (marked.count((*it).substr((*it).size()-k_size))) {
					r = false;
				}
				if (marked.count((*it).substr(0,k_size))) {
					l = false;
				}
				if (!l && !r){
					continue;
				}
				set<string> ret = branch(*it, l, r);
				new_front.insert(ret.begin(), ret.end());
			}
			for (set<string>::iterator it = front.begin(); it != front.end(); ++it) {
				bool l = true, r = true;

				if (marked.count((*it).substr((*it).size()-k_size))) {
					r = false;
				}
				if (marked.count((*it).substr(0,k_size))) {
					l = false;
				}
				if (!l && !r){
					continue;
				}

				if (r) {
					marked.insert((*it).substr((*it).size()-k_size));
					marked.insert(reverse_complement((*it).substr((*it).size()-k_size)));
				}
				if (l) {
					marked.insert((*it).substr(0,k_size));
					marked.insert(reverse_complement((*it).substr(0,k_size)));
				}
			}

			if(new_front.empty()) {
				reason = -1;
				break;
			}
			depth++;
			front = new_front;

			int width = new_front.size();

			if (depth > 500 || width > 20) {
				reason = 0;
				break;
			}
			if (new_front.size() == 1) {
				bool r = true, l = true;
				if (marked.count((*new_front.begin()).substr((*new_front.begin()).size()-k_size))) {
					r = false;
				}
				if (marked.count((*new_front.begin()).substr(0,k_size))) {
					l = false;
				}
				int bla;
				set<string> ret = branch(*(new_front.begin()), l, r);
				if (ret.size() >= 1) {
					continue;
				}
				reason = 1;
				break;
			}
		} while (1);

		if (reason != -1) {
			for (set<string>::iterator it = front.begin(); it != front.end(); ++it) {
				if ((*it).size() > 2*k_size + 1) {
					contigs.insert(*it);
				}
			}
		}

		index++;
		for (int i = index; i < nodes.size(); ++i) {
			if(!marked.count(nodes[i])) {
				index = i;
				break;
			}
		}
	}
}

/* 
 * Method that takes a contig s and tries to expand it 
 * to the left and right if corresponding flags: left and right
 * are true
 */
set<string> dbg::branch(string s, bool left, bool right) {

	char bla[4] = {'A', 'C', 'G', 'T'};
	list<char> l(bla, bla + sizeof(bla) / sizeof(char));
	set<string> ret;
	set<string> tmp;

	if (left && right) {
		for (list<char>::iterator c = l.begin(); c != l.end(); ++c) {
			string h = *c + s.substr(0, k_size - 1);
			if (left && filter.contains(h) && !cFP.count(canonical(h)) && !marked.count(h)) {
				tmp.insert(*c + s);
			}
		}
		for (list<char>::iterator c = l.begin(); c != l.end(); ++c) {
			for (set<string>::iterator it = tmp.begin(); it != tmp.end(); ++it) {
				string h = (*it).substr((*it).size() - k_size + 1) + *c;
				if (right && filter.contains(h) && !cFP.count(canonical(h)) && !marked.count(h)) {
					ret.insert((*it) + *c);
				}
			}
		}
	} else if (left) {
		for (list<char>::iterator c = l.begin(); c != l.end(); ++c) {
			string h = *c + s.substr(0, k_size - 1);
			if (left && filter.contains(h) && !cFP.count(canonical(h)) && !marked.count(h)) {
				ret.insert(*c + s);
			}
		}
	} else if (right) {
		for (list<char>::iterator c = l.begin(); c != l.end(); ++c) {
			string h = s.substr(s.size() - k_size + 1) + *c;
			if (right && filter.contains(h) && !cFP.count(canonical(h)) && !marked.count(h)) {
				ret.insert(s + *c);
			}
		}
	}

	return ret;
}

/* 
 * Method for calcuating caonical k-mer.
 * K-mer is canonical if it is less than
 * its reverse complement in means of ASCII.
 */
string dbg::canonical (string s1) {
	string s2 = reverse_complement(s1);
	return s1 < s2 ? s1 : s2;
}

/* 
 * Method for calculation reverse complement of a k-mer.
 */
string dbg::reverse_complement (string s) {
	string tmp = "";
	for(string::iterator it = s.begin(); it != s.end(); ++it) {
		tmp = complement(string(1,*it)) + tmp;
	}
	return tmp;
}

/* 
 * Method for calculation of a complement of a nucleotid.
 */
string dbg::complement(string c) {
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

int main (int argc, char *argv[]) {

	dbg d(argv[1]);
	d.compute_cFP(1000000000);
	d.traverse_graph();
	//cout << d.contigs.size() << endl;
}